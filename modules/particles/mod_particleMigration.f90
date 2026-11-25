module mod_particleMigration
  !===========================================================================
  ! Particle migration between MPI ranks for the z-slab decomposition
  ! (mod_zslab). After the push, every particle whose cell plane
  ! iz = int(z/hz)+1 has left this rank's slab is removed from its
  ! ParticleSet and sent to another rank, where it is appended.
  !
  ! Two modes:
  !   - neighbour mode (global_mode = .false., every step after the push):
  !     particles go to the left/right neighbour in z (a ring when z is
  !     periodic). A particle that has to travel further than one slab is
  !     forwarded again in a following round; rounds repeat until no rank
  !     holds a foreign particle (one small MPI_Allreduce per round).
  !   - global mode (global_mode = .true., after a restart or an injection,
  !     where particles may sit anywhere): one MPI_Alltoallv straight to the
  !     owner rank.
  !
  ! Payload per particle (9 doubles): x, y, z, vx, vy, vz, w, flag_cex,
  ! species. Particles flagged dead are never sent: they are removed (with
  ! their energy/momentum correction) by the next push on the rank that
  ! killed them, and collisions only act on that rank's own particles.
  ! Received particles are dealt out evenly over the OMP threads' sets.
  ! Cell lists (cell_id/cell_count/cell_start) are left stale, as after a
  ! push: they are rebuilt by the next sort, before any collision.
  !===========================================================================
  use iso_fortran_env, only: int8, int32, int64, real64
  use mpi
  use mod_particles,   only: ParticleSet
  use mod_zslab,       only: zs_comm, zs_rank, zs_size, zs_left, zs_right, &
                             zs_cz_lo, zs_cz_hi, zslab_iz_of_z, &
                             zslab_owner_of_iz, zslab_next_hop
  implicit none
  private

  public :: migrate_particles

  integer, parameter :: NPAY = 9

  type :: LeaverBuf
    integer(int32) :: n = 0
    real(real64), allocatable :: d(:,:)       ! (NPAY, cap)
    integer,      allocatable :: dest(:)      ! destination rank
  end type LeaverBuf

  type(LeaverBuf), allocatable, save :: tb(:)
  real(real64),    allocatable, save :: sendbuf(:), recvbuf(:)

  ! Diagnostics (cumulative, this rank): particles sent, rounds run.
  integer(int64), public, save :: mig_sent_total   = 0_int64
  integer(int64), public, save :: mig_rounds_total = 0_int64
  integer(int64), public, save :: mig_calls_total  = 0_int64

contains

  subroutine migrate_particles(part, global_mode)
    type(ParticleSet), intent(inout) :: part(:,:)
    logical,           intent(in)    :: global_mode

    integer(int64) :: n_foreign, n_foreign_glob
    integer :: ierr, iround

    if (zs_size <= 1) return
    mig_calls_total = mig_calls_total + 1_int64

    iround = 0
    do
      iround = iround + 1
      call migrate_round(part, global_mode, n_foreign)
      mig_rounds_total = mig_rounds_total + 1_int64
      if (global_mode) exit
      call MPI_Allreduce(n_foreign, n_foreign_glob, 1, MPI_INTEGER8, MPI_SUM, zs_comm, ierr)
      if (n_foreign_glob == 0_int64) exit
      if (iround > 4*zs_size) error stop 'migrate_particles: particles keep moving between ranks'
    end do

    if (check_enabled()) call check_all_local(part)
  end subroutine migrate_particles


  logical function check_enabled()
    ! Debug aid: set LEPIC_CHECK_MIGRATION (any value) to verify after every
    ! migration that each live particle lies in this rank's slab.
    integer, save :: state = -1
    integer :: length, stat
    if (state < 0) then
      call get_environment_variable('LEPIC_CHECK_MIGRATION', length=length, status=stat)
      state = merge(1, 0, stat == 0)
    end if
    check_enabled = (state == 1)
  end function check_enabled


  subroutine check_all_local(part)
    type(ParticleSet), intent(in) :: part(:,:)
    integer :: s, t
    integer(int32) :: i, iz
    do t = 1, size(part,2)
      do s = 1, size(part,1)
        if (.not. allocated(part(s,t)%pv)) cycle
        do i = 1, part(s,t)%n
          if (part(s,t)%flag_dead(i) /= 0_int8) cycle
          iz = zslab_iz_of_z(part(s,t)%pv(3,i))
          if (iz < zs_cz_lo .or. iz > zs_cz_hi) then
            write(*,*) 'migrate_particles check: rank', zs_rank, ' species', s, ' thread', t, &
                       ' particle', i, ' iz', iz, ' slab', zs_cz_lo, zs_cz_hi
            error stop 'migrate_particles: particle left outside its slab'
          end if
        end do
      end do
    end do
  end subroutine check_all_local


  !=========================================================================
  ! One exchange. n_foreign_out: received particles that still do not
  ! belong to this rank (neighbour mode, they are forwarded next round).
  !=========================================================================
  subroutine migrate_round(part, global_mode, n_foreign_out)
    type(ParticleSet), intent(inout) :: part(:,:)
    logical,           intent(in)    :: global_mode
    integer(int64),    intent(out)   :: n_foreign_out

    integer :: ntype, nproc, P, ierr, r, s, t, k
    integer :: nslot
    integer, allocatable :: slot_rank(:)
    integer(int64), allocatable :: cnt(:,:)       ! (nslot, nproc)
    integer(int64), allocatable :: off(:,:)       ! (nslot, nproc)
    integer(int64) :: ntot_send, ntot_recv, pos
    integer, allocatable :: scount(:), rcount(:), sdispl(:), rdispl(:)
    integer :: nleft_s, nright_s, nleft_r, nright_r
    integer :: i, j

    ntype = size(part,1)
    nproc = size(part,2)
    P     = zs_size
    n_foreign_out = 0_int64

    if (global_mode) then
      nslot = P
      allocate(slot_rank(P))
      do r = 0, P-1
        slot_rank(r+1) = r
      end do
    else
      nslot = 2
      allocate(slot_rank(2))
      slot_rank(1) = zs_left
      slot_rank(2) = zs_right
    end if

    if (.not. allocated(tb)) then
      allocate(tb(nproc))
    else if (size(tb) < nproc) then
      deallocate(tb)
      allocate(tb(nproc))
    end if

    allocate(cnt(nslot, nproc), off(nslot, nproc))
    cnt = 0_int64

    ! ---- 1. compact every set, collecting the leavers per thread --------
    !$omp parallel do private(t,s,i,j,k) schedule(static)
    do t = 1, nproc
      tb(t)%n = 0
      do s = 1, ntype
        call extract_leavers(part(s,t), s, tb(t), global_mode)
      end do
      do k = 1, tb(t)%n
        cnt(dest_slot(tb(t)%dest(k)), t) = cnt(dest_slot(tb(t)%dest(k)), t) + 1_int64
      end do
    end do
    !$omp end parallel do

    ! ---- 2. send buffer, grouped by destination slot ----------------------
    pos = 0_int64
    do k = 1, nslot
      do t = 1, nproc
        off(k,t) = pos
        pos = pos + cnt(k,t)
      end do
    end do
    ntot_send = pos
    call ensure_buf(sendbuf, NPAY*ntot_send)

    !$omp parallel do private(t,i,k) schedule(static)
    do t = 1, nproc
      block
        integer(int64) :: cur(nslot), q
        cur = off(:,t)
        do i = 1, tb(t)%n
          k = dest_slot(tb(t)%dest(i))
          q = cur(k)
          sendbuf(NPAY*q+1:NPAY*q+NPAY) = tb(t)%d(:,i)
          cur(k) = q + 1_int64
        end do
      end block
    end do
    !$omp end parallel do

    ! ---- 3. exchange --------------------------------------------------------
    if (global_mode) then
      allocate(scount(0:P-1), rcount(0:P-1), sdispl(0:P-1), rdispl(0:P-1))
      do r = 0, P-1
        if (sum(cnt(r+1,:))*NPAY > int(huge(1),int64)) &
          error stop 'migrate_particles: too many particles for one MPI message'
        scount(r) = int(sum(cnt(r+1,:))) * NPAY
      end do
      call MPI_Alltoall(scount, 1, MPI_INTEGER, rcount, 1, MPI_INTEGER, zs_comm, ierr)
      sdispl(0) = 0; rdispl(0) = 0
      do r = 1, P-1
        sdispl(r) = sdispl(r-1) + scount(r-1)
        rdispl(r) = rdispl(r-1) + rcount(r-1)
      end do
      ntot_recv = (int(rdispl(P-1),int64) + int(rcount(P-1),int64)) / NPAY
      call ensure_buf(recvbuf, NPAY*ntot_recv)
      call MPI_Alltoallv(sendbuf, scount, sdispl, MPI_DOUBLE_PRECISION, &
                         recvbuf, rcount, rdispl, MPI_DOUBLE_PRECISION, zs_comm, ierr)
      deallocate(scount, rcount, sdispl, rdispl)
    else
      nleft_s  = int(sum(cnt(1,:))) * NPAY
      nright_s = int(sum(cnt(2,:))) * NPAY
      ! counts: send to right / receive from left, then the reverse
      call MPI_Sendrecv(nright_s, 1, MPI_INTEGER, zs_right, 801, &
                        nleft_r,  1, MPI_INTEGER, zs_left,  801, zs_comm, MPI_STATUS_IGNORE, ierr)
      call MPI_Sendrecv(nleft_s,  1, MPI_INTEGER, zs_left,  802, &
                        nright_r, 1, MPI_INTEGER, zs_right, 802, zs_comm, MPI_STATUS_IGNORE, ierr)
      if (zs_left  == MPI_PROC_NULL) nleft_r  = 0
      if (zs_right == MPI_PROC_NULL) nright_r = 0
      ntot_recv = (int(nleft_r,int64) + int(nright_r,int64)) / NPAY
      call ensure_buf(recvbuf, NPAY*ntot_recv)
      ! data: slot 1 (left) occupies sendbuf(1:nleft_s), slot 2 follows
      call MPI_Sendrecv(sendbuf(nleft_s+1), nright_s, MPI_DOUBLE_PRECISION, zs_right, 803, &
                        recvbuf(1),         nleft_r,  MPI_DOUBLE_PRECISION, zs_left,  803, &
                        zs_comm, MPI_STATUS_IGNORE, ierr)
      call MPI_Sendrecv(sendbuf(1),         nleft_s,  MPI_DOUBLE_PRECISION, zs_left,  804, &
                        recvbuf(nleft_r+1), nright_r, MPI_DOUBLE_PRECISION, zs_right, 804, &
                        zs_comm, MPI_STATUS_IGNORE, ierr)
    end if

    mig_sent_total = mig_sent_total + ntot_send

    ! ---- 4. append the received particles ----------------------------------
    if (ntot_recv > 0_int64) call append_received(part, ntot_recv, n_foreign_out)

    deallocate(cnt, off, slot_rank)

  contains

    pure integer function dest_slot(dest) result(k)
      integer, intent(in) :: dest
      if (global_mode) then
        k = dest + 1
      else
        ! Ring of two ranks: left == right; leavers were tagged with the
        ! slot itself (see extract_leavers), encoded as -1 / -2.
        k = merge(1, 2, dest == -1)
      end if
    end function dest_slot

  end subroutine migrate_round


  !=========================================================================
  ! Remove from one set the live particles outside this rank's slab and
  ! append them to the thread's leaver buffer. Neighbour mode tags each
  ! leaver with -1 (to the left neighbour) or -2 (to the right one) so the
  ! two directions stay distinct even when both neighbours are the same
  ! rank (2-rank ring); global mode tags it with its owner rank.
  !=========================================================================
  subroutine extract_leavers(p, species, buf, global_mode)
    type(ParticleSet), intent(inout) :: p
    integer,           intent(in)    :: species
    type(LeaverBuf),   intent(inout) :: buf
    logical,           intent(in)    :: global_mode

    integer(int32) :: i, j, k, n_old
    logical :: has_dead, removed_any

    if (.not. allocated(p%pv)) return
    if (p%n <= 0_int32) return
    has_dead = allocated(p%flag_dead)

    ! Leavers are a small fraction: remove each by moving the current last
    ! particle into its slot (O(leavers) writes instead of a stable
    ! compaction that would rewrite every particle after the first one).
    ! Particle order changes slightly; the next sort restores it.
    n_old = p%n
    j = n_old

    if (p%mig_valid .and. .not. global_mode) then
      ! The mover already flagged the survivors that left the slab
      ! (p%mig_idx, increasing); only particles appended after that push
      ! (indices > mig_n_checked, e.g. secondary electrons) are unchecked.
      ! Visiting candidates in DECREASING index order makes the swap-with-
      ! last removal safe: the particle moved into a freed slot always
      ! comes from a higher index that has already been visited.
      do i = n_old, p%mig_n_checked + 1_int32, -1_int32
        removed_any = try_remove(i)
      end do
      do k = p%n_mig, 1_int32, -1_int32
        removed_any = try_remove(p%mig_idx(k))
      end do
    else
      i = 1_int32
      do while (i <= j)
        if (try_remove(i)) cycle   ! slot i now holds another particle
        i = i + 1_int32
      end do
    end if
    p%mig_valid = .false.
    p%n_mig     = 0_int32

    if (has_dead) p%flag_dead(j+1:n_old) = 0_int8
    if (allocated(p%flag_cex)) p%flag_cex(j+1:n_old) = 0_int32
    p%n = j

  contains

    logical function try_remove(ii) result(removed)
      ! Send particle ii away if it is live and outside the slab; its
      ! slot is then refilled with the current last particle j.
      integer(int32), intent(in) :: ii
      integer(int32) :: iz
      integer :: owner, hop

      removed = .false.
      if (ii > j) return
      iz = zslab_iz_of_z(p%pv(3,ii))
      if (iz >= zs_cz_lo .and. iz <= zs_cz_hi) return
      if (has_dead) then
        if (p%flag_dead(ii) /= 0_int8) return
      end if

      owner = zslab_owner_of_iz(iz)
      if (global_mode) then
        hop = owner
      else if (zs_left == zs_right) then
        ! 2-rank ring: both neighbours are the same rank; use the side the
        ! particle left through so the two directions stay paired.
        hop = merge(-2, -1, iz > zs_cz_hi)
      else
        hop = merge(-1, -2, zslab_next_hop(owner) == zs_left)
      end if
      call push_leaver(ii, hop)
      if (ii /= j) call move_last(j, ii)
      j = j - 1_int32
      removed = .true.
    end function try_remove

    subroutine move_last(src, dst)
      integer(int32), intent(in) :: src, dst
      p%pv(:,dst) = p%pv(:,src)
      if (allocated(p%w))         p%w(dst)         = p%w(src)
      if (allocated(p%sp))        p%sp(dst)        = p%sp(src)
      if (allocated(p%cell_id))   p%cell_id(dst)   = p%cell_id(src)
      if (has_dead)               p%flag_dead(dst) = p%flag_dead(src)
      if (allocated(p%flag_cex))  p%flag_cex(dst)  = p%flag_cex(src)
    end subroutine move_last

    subroutine push_leaver(ii, d)
      integer(int32), intent(in) :: ii
      integer,        intent(in) :: d
      real(real64), allocatable :: dnew(:,:)
      integer,      allocatable :: tnew(:)
      integer :: cap

      if (.not. allocated(buf%d)) then
        allocate(buf%d(NPAY, 1024), buf%dest(1024))
      else if (buf%n >= size(buf%dest)) then
        cap = 2*size(buf%dest)
        allocate(dnew(NPAY, cap), tnew(cap))
        dnew(:,1:buf%n) = buf%d(:,1:buf%n)
        tnew(1:buf%n)   = buf%dest(1:buf%n)
        call move_alloc(dnew, buf%d)
        call move_alloc(tnew, buf%dest)
      end if
      buf%n = buf%n + 1
      buf%d(1:6,buf%n) = p%pv(1:6,ii)
      buf%d(7,buf%n)   = 1.0_real64
      if (allocated(p%w)) buf%d(7,buf%n) = p%w(ii)
      buf%d(8,buf%n)   = 0.0_real64
      if (allocated(p%flag_cex)) buf%d(8,buf%n) = real(p%flag_cex(ii), real64)
      buf%d(9,buf%n)   = real(species, real64)
      buf%dest(buf%n)  = d
    end subroutine push_leaver

  end subroutine extract_leavers


  !=========================================================================
  ! Deal the received particles out evenly over the threads of each
  ! species and copy them in (parallel over received particles).
  !=========================================================================
  subroutine append_received(part, nrecv, n_foreign)
    type(ParticleSet), intent(inout) :: part(:,:)
    integer(int64),    intent(in)    :: nrecv
    integer(int64),    intent(out)   :: n_foreign

    integer :: ntype, nproc, s, t
    integer(int64) :: j, q
    integer(int64), allocatable :: rc(:), rank_in_s(:)
    integer(int64), allocatable :: blk0(:,:)     ! (ntype, nproc+1) block starts
    integer(int32), allocatable :: n_old(:,:)
    integer(int32) :: slot
    integer :: iz

    ntype = size(part,1)
    nproc = size(part,2)

    allocate(rc(ntype), rank_in_s(nrecv), blk0(ntype, nproc+1), n_old(ntype, nproc))
    rc = 0_int64
    do j = 1_int64, nrecv
      s = nint(recvbuf(NPAY*(j-1)+9))
      if (s < 1 .or. s > ntype) error stop 'migrate_particles: bad species in received particle'
      rank_in_s(j) = rc(s)
      rc(s) = rc(s) + 1_int64
    end do

    do s = 1, ntype
      do t = 1, nproc+1
        blk0(s,t) = (rc(s) * int(t-1,int64)) / int(nproc,int64)
      end do
    end do

    !$omp parallel do private(t,s) schedule(static)
    do t = 1, nproc
      do s = 1, ntype
        n_old(s,t) = 0_int32
        if (allocated(part(s,t)%pv)) n_old(s,t) = part(s,t)%n
        if (blk0(s,t+1) > blk0(s,t)) then
          call part(s,t)%ensure_capacity(n_old(s,t) + int(blk0(s,t+1) - blk0(s,t), int32))
        end if
      end do
    end do
    !$omp end parallel do

    n_foreign = 0_int64
    !$omp parallel do private(j,s,q,t,slot,iz) reduction(+:n_foreign) schedule(static)
    do j = 1_int64, nrecv
      s = nint(recvbuf(NPAY*(j-1)+9))
      q = rank_in_s(j)
      ! thread t such that blk0(s,t) <= q < blk0(s,t+1)
      t = int((q * int(nproc,int64)) / max(1_int64, rc(s))) + 1
      do while (t > 1)
        if (q >= blk0(s,t)) exit
        t = t - 1
      end do
      do while (t < nproc)
        if (q < blk0(s,t+1)) exit
        t = t + 1
      end do
      slot = n_old(s,t) + int(q - blk0(s,t), int32) + 1_int32
      part(s,t)%pv(1:6,slot)    = recvbuf(NPAY*(j-1)+1:NPAY*(j-1)+6)
      part(s,t)%w(slot)         = recvbuf(NPAY*(j-1)+7)
      part(s,t)%sp(slot)        = s
      part(s,t)%flag_dead(slot) = 0_int8
      part(s,t)%flag_cex(slot)  = nint(recvbuf(NPAY*(j-1)+8), int32)
      part(s,t)%cell_id(slot)   = 0_int32
      iz = zslab_iz_of_z(recvbuf(NPAY*(j-1)+3))
      if (iz < zs_cz_lo .or. iz > zs_cz_hi) n_foreign = n_foreign + 1_int64
    end do
    !$omp end parallel do

    do t = 1, nproc
      do s = 1, ntype
        part(s,t)%n = n_old(s,t) + int(blk0(s,t+1) - blk0(s,t), int32)
      end do
    end do

    deallocate(rc, rank_in_s, blk0, n_old)
  end subroutine append_received


  subroutine ensure_buf(b, nwanted)
    real(real64), allocatable, intent(inout) :: b(:)
    integer(int64), intent(in) :: nwanted
    integer(int64) :: cap, need
    ! NPAY of slack: MPI calls are given b(count+1) as the start of an
    ! empty second message, which must still be a valid element.
    need = nwanted + NPAY
    if (allocated(b)) then
      if (size(b,kind=int64) >= need) return
      cap = max(need, 2_int64*size(b,kind=int64))
      deallocate(b)
    else
      cap = max(need, 65536_int64)
    end if
    allocate(b(cap))
  end subroutine ensure_buf

end module mod_particleMigration
