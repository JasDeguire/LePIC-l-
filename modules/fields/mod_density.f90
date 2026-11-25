module mod_density
  !!
  !! Density reduction / charge-density builder, z-slab domain decomposition.
  !!
  !! Each MPI rank only holds the particles of its z-slab (mod_zslab), and
  !! np_thread only covers that slab's deposit window of node planes
  !! zs_nlo:zs_nhi (= owned cell planes k0+1..k0+m, plus plane k0+m+1 that
  !! their upper corners also reach). Building the density is then:
  !!
  !!   1. sum np_thread over the OMP threads into a rank window buffer that
  !!      spans planes k0 .. k0+m+1 (one halo plane on each side);
  !!   2. periodic stitching in y (in-plane, as legacy);
  !!   3. add the partial top plane k0+m+1 into the right neighbour's first
  !!      plane (open chain: the last rank keeps plane nz+1, a node it owns);
  !!   4. periodic stitching in z between plane 1 (rank 0) and plane nz+1
  !!      (last rank), same rule as legacy;
  !!   5. halo: plane k0 from the left neighbour, plane k0+m+1 from the
  !!      right one, so the window holds final values on k0..k0+m+1 - the
  !!      exact range of the Poisson right-hand side this rank solves for.
  !!
  !! Every operation of the legacy convention (np_thread -> periodic
  !! stitching -> sum) is linear, so applying it to the rank sums instead
  !! of per thread gives the same density. Neighbour exchange of one plane
  !! replaces the former MPI_Allreduce of rho over the whole grid.
  !!
  !! fld%np keeps its global shape: after the per-step reduction it holds
  !! only this rank's owned planes (zeros elsewhere), and
  !! sync_species_density MPI-sums it into the global density for the
  !! readers that need it (collisions, outputs, RF/skin-depth updates).
  !!
  use iso_fortran_env, only: real64, int32, int8
  use mpi
  use mod_constants, only: qe, eps0
  use mod_zslab,     only: zs_comm, zs_rank, zs_size, zs_nlo, zs_nhi, &
                           zs_cz_lo, zs_cz_hi, zs_k0, zs_nz, &
                           zs_left_open, zs_right_open

  implicit none
  private

  public :: reduce_species_density
  public :: reduce_density_and_rho
  public :: sync_species_density
  public :: clear_np_planes
  public :: build_rho_from_np
  public :: density_max_per_species
  public :: average_species_density

  ! Rank window buffer: planes zs_k0 .. zs_nhi (= k0 .. k0+m+1), ntype.
  real(real64), allocatable, save :: win(:,:,:,:)
  ! One-plane (all species) send/receive buffers.
  real(real64), allocatable, save :: sbuf(:,:,:), rbuf(:,:,:)
  ! Whether any node of planes 1 / nz+1 is periodic (bcnd == 0): decided
  ! once from the (replicated, static) bcnd array.
  integer, save :: zstitch_state = -1   ! -1 unknown, 0 no, 1 yes

contains

  subroutine reduce_species_density(n, bcnd, np_thread, ntype, nproc, mpi_comm, np_red)
    ! np_thread -> assembled density -> np_red, made globally valid
    ! (MPI-summed) right away. Used at initial load / restart.
    integer(int32), intent(in)    :: n(3)
    integer,        intent(in)    :: ntype, nproc, mpi_comm
    integer,        intent(in)    :: bcnd(0:n(1)+2,0:n(2)+2,0:n(3)+2)
    real(real64),   intent(inout) :: np_thread(0:n(1)+2,0:n(2)+2,zs_nlo:zs_nhi,ntype,nproc)
    real(real64),   intent(out)   :: np_red(0:n(1)+2,0:n(2)+2,0:n(3)+2,ntype)

    np_red = 0.0_real64
    call assemble_density(n, bcnd, np_thread, ntype, nproc, np_red, np_red_dirty=.false.)
    call sync_species_density(n, ntype, mpi_comm, np_red)
  end subroutine reduce_species_density


  subroutine reduce_density_and_rho(n, bcnd, np_thread, ntype, nproc, mpi_comm, charge, &
                                    np_red, rho, live, np_red_dirty)
    ! Per-step path (legacy calc_rho). On exit:
    !   rho   : final charge density on planes k0..k0+m+1 (this rank's
    !           Poisson range); other planes are not touched.
    !   np_red: this rank's owned planes only (zeros elsewhere) - the
    !           caller must run sync_species_density before reading it.
    ! np_red_dirty: np_red currently holds other ranks' planes (it was
    !   synced since the last reduction), which must be zeroed first.
    integer(int32), intent(in)    :: n(3)
    integer,        intent(in)    :: ntype, nproc, mpi_comm
    integer,        intent(in)    :: bcnd(0:n(1)+2,0:n(2)+2,0:n(3)+2)
    real(real64),   intent(inout) :: np_thread(0:n(1)+2,0:n(2)+2,zs_nlo:zs_nhi,ntype,nproc)
    real(real64),   intent(in)    :: charge(ntype)
    real(real64),   intent(inout) :: np_red(0:n(1)+2,0:n(2)+2,0:n(3)+2,ntype)
    real(real64),   intent(inout) :: rho(0:n(1)+2,0:n(2)+2,0:n(3)+2)
    integer(int8),  intent(in),   optional :: live(zs_nlo:zs_nhi,ntype,nproc)
    logical,        intent(in)    :: np_red_dirty

    integer :: ix, iy, iz, ptype, nx1
    real(real64) :: r

    if (.false.) print *, mpi_comm

    call assemble_density(n, bcnd, np_thread, ntype, nproc, np_red, np_red_dirty, live)

    ! rho on the whole window k0..k0+m+1 (x/y ghosts zero, as legacy).
    nx1 = n(1) + 1
    !$omp parallel do collapse(2) private(iz,iy,ix,ptype,r) schedule(static)
    do iz = zs_k0, zs_nhi
      do iy = 0, n(2)+2
        if (iy < 1 .or. iy > n(2)+1 .or. iz < 1 .or. iz > n(3)+1) then
          rho(:,iy,iz) = 0.0_real64
          cycle
        end if
        rho(0,iy,iz)     = 0.0_real64
        rho(nx1+1,iy,iz) = 0.0_real64
        do ix = 1, nx1
          r = 0.0_real64
          do ptype = 1, ntype
            r = r - charge(ptype) * win(ix,iy,iz,ptype)
          end do
          rho(ix,iy,iz) = r
        end do
      end do
    end do
    !$omp end parallel do
  end subroutine reduce_density_and_rho


  subroutine sync_species_density(n, ntype, mpi_comm, np_red)
    ! MPI-sum a rank-local np_red in place (no-op on a single rank). Each
    ! rank contributes only its owned planes, so the sum is the global
    ! density.
    integer(int32), intent(in)    :: n(3)
    integer,        intent(in)    :: ntype, mpi_comm
    real(real64),   intent(inout) :: np_red(0:n(1)+2,0:n(2)+2,0:n(3)+2,ntype)

    integer :: ierr, mpi_size

    call MPI_Comm_size(mpi_comm, mpi_size, ierr)
    if (mpi_size > 1) then
      call MPI_Allreduce(MPI_IN_PLACE, np_red, &
          (n(1)+3)*(n(2)+3)*(n(3)+3)*ntype, MPI_DOUBLE_PRECISION, MPI_SUM, mpi_comm, ierr)
    end if
  end subroutine sync_species_density


  subroutine clear_np_planes(np_local, n, live)
    ! Zero the planes of one thread's np_thread slice flagged in live (and
    ! reset those flags), leaving the whole slice zero.
    integer(int32), intent(in)    :: n(3)
    real(real64),   intent(inout) :: np_local(0:n(1)+2,0:n(2)+2,zs_nlo:zs_nhi)
    integer(int8),  intent(inout) :: live(zs_nlo:zs_nhi)

    integer :: iz

    do iz = zs_nlo, zs_nhi
      if (live(iz) /= 0_int8) then
        np_local(:,:,iz) = 0.0_real64
        live(iz) = 0_int8
      end if
    end do
  end subroutine clear_np_planes


  !=========================================================================
  ! Core of the reduction - see the module header for the five steps.
  !=========================================================================
  subroutine assemble_density(n, bcnd, np_thread, ntype, nproc, np_red, np_red_dirty, live)
    integer(int32), intent(in)    :: n(3)
    integer,        intent(in)    :: ntype, nproc
    integer,        intent(in)    :: bcnd(0:n(1)+2,0:n(2)+2,0:n(3)+2)
    real(real64),   intent(in)    :: np_thread(0:n(1)+2,0:n(2)+2,zs_nlo:zs_nhi,ntype,nproc)
    real(real64),   intent(inout) :: np_red(0:n(1)+2,0:n(2)+2,0:n(3)+2,ntype)
    logical,        intent(in)    :: np_red_dirty
    integer(int8),  intent(in), optional :: live(zs_nlo:zs_nhi,ntype,nproc)

    integer :: ix, iy, iz, ptype, iproc, nx1, ny1, nz1
    integer :: npl, ierr, k_own_hi
    logical :: ranged

    nx1 = n(1) + 1
    ny1 = n(2) + 1
    nz1 = n(3) + 1
    ranged = present(live)
    npl = (n(1)+3)*(n(2)+3)*ntype

    call ensure_buffers(n, ntype)

    ! ---- 1. thread sum into the window (planes zs_nlo..zs_nhi) ---------
    !$omp parallel do collapse(2) private(iz,iy,ix,ptype,iproc) schedule(static)
    do iz = zs_k0, zs_nhi
      do iy = 0, n(2)+2
        do ptype = 1, ntype
          win(:,iy,iz,ptype) = 0.0_real64
        end do
        if (iz < zs_nlo .or. iy < 1 .or. iy > ny1) cycle
        do ptype = 1, ntype
          do iproc = 1, nproc
            if (ranged) then
              if (live(iz,ptype,iproc) == 0_int8) cycle
            end if
            do ix = 1, nx1
              win(ix,iy,iz,ptype) = win(ix,iy,iz,ptype) + np_thread(ix,iy,iz,ptype,iproc)
            end do
          end do
        end do
      end do
    end do
    !$omp end parallel do

    ! ---- 2. periodic stitching in y (in-plane) --------------------------
    !$omp parallel do collapse(2) private(iz,ix,ptype) schedule(static)
    do iz = max(1, zs_nlo), min(nz1, zs_nhi)
      do ix = 1, nx1
        if (bcnd(ix,1,iz) == 0) then
          do ptype = 1, ntype
            win(ix,1,iz,ptype) = 0.5_real64 * (win(ix,1,iz,ptype) + win(ix,ny1,iz,ptype))
          end do
        end if
        if (bcnd(ix,ny1,iz) == 0) then
          do ptype = 1, ntype
            win(ix,ny1,iz,ptype) = win(ix,1,iz,ptype)
          end do
        end if
      end do
    end do
    !$omp end parallel do

    ! ---- 3. add the partial top plane into the right neighbour ----------
    if (zs_size > 1) then
      call pack_plane(zs_nhi, sbuf)
      call MPI_Sendrecv(sbuf, npl, MPI_DOUBLE_PRECISION, zs_right_open, 701, &
                        rbuf, npl, MPI_DOUBLE_PRECISION, zs_left_open,  701, &
                        zs_comm, MPI_STATUS_IGNORE, ierr)
      if (zs_left_open /= MPI_PROC_NULL) call add_plane(zs_nlo, rbuf)
    end if

    ! ---- 4. periodic stitching in z: plane 1 <-> plane nz+1 -------------
    if (zstitch_state < 0) call detect_zstitch(n, bcnd)
    if (zstitch_state == 1) call stitch_z(n, bcnd, npl)

    ! ---- 5. halo planes k0 (from left) and k0+m+1 (from right) ----------
    if (zs_size > 1) then
      call pack_plane(zs_cz_lo, sbuf)
      call MPI_Sendrecv(sbuf, npl, MPI_DOUBLE_PRECISION, zs_left_open,  702, &
                        rbuf, npl, MPI_DOUBLE_PRECISION, zs_right_open, 702, &
                        zs_comm, MPI_STATUS_IGNORE, ierr)
      if (zs_right_open /= MPI_PROC_NULL) call unpack_plane(zs_nhi, rbuf)

      call pack_plane(zs_cz_hi, sbuf)
      call MPI_Sendrecv(sbuf, npl, MPI_DOUBLE_PRECISION, zs_right_open, 703, &
                        rbuf, npl, MPI_DOUBLE_PRECISION, zs_left_open,  703, &
                        zs_comm, MPI_STATUS_IGNORE, ierr)
      if (zs_left_open /= MPI_PROC_NULL) call unpack_plane(zs_k0, rbuf)
    end if
    ! Legacy convention: planes 0 and nz+2 of the density are always zero.
    if (zs_k0 == 0) win(:,:,0,:) = 0.0_real64

    ! ---- owned planes -> np_red ------------------------------------------
    ! Owned: cell planes k0+1..k0+m, plus node plane nz+1 on the last rank.
    k_own_hi = zs_cz_hi
    if (zs_rank == zs_size-1) k_own_hi = nz1

    if (np_red_dirty) then
      !$omp parallel do private(iz,ptype) schedule(static)
      do iz = 0, n(3)+2
        if (iz >= zs_cz_lo .and. iz <= k_own_hi) cycle
        do ptype = 1, ntype
          np_red(:,:,iz,ptype) = 0.0_real64
        end do
      end do
      !$omp end parallel do
    end if

    !$omp parallel do private(iz,ptype) schedule(static)
    do iz = zs_cz_lo, k_own_hi
      do ptype = 1, ntype
        np_red(:,:,iz,ptype) = win(:,:,iz,ptype)
      end do
    end do
    !$omp end parallel do

  contains

    subroutine pack_plane(iz, buf)
      integer,      intent(in)  :: iz
      real(real64), intent(out) :: buf(0:n(1)+2,0:n(2)+2,ntype)
      integer :: p
      do p = 1, ntype
        buf(:,:,p) = win(:,:,iz,p)
      end do
    end subroutine pack_plane

    subroutine unpack_plane(iz, buf)
      integer,      intent(in) :: iz
      real(real64), intent(in) :: buf(0:n(1)+2,0:n(2)+2,ntype)
      integer :: p
      do p = 1, ntype
        win(:,:,iz,p) = buf(:,:,p)
      end do
    end subroutine unpack_plane

    subroutine add_plane(iz, buf)
      integer,      intent(in) :: iz
      real(real64), intent(in) :: buf(0:n(1)+2,0:n(2)+2,ntype)
      integer :: p
      do p = 1, ntype
        win(:,:,iz,p) = win(:,:,iz,p) + buf(:,:,p)
      end do
    end subroutine add_plane

  end subroutine assemble_density


  subroutine stitch_z(n, bcnd, npl)
    ! Legacy z-stitching, on the assembled density:
    !   where bcnd(ix,iy,1)    == 0 : np(1)    = (np(1) + np(nz+1)) / 2
    !   where bcnd(ix,iy,nz+1) == 0 : np(nz+1) = np(1)   (updated value)
    ! Plane 1 lives on rank 0, plane nz+1 on the last rank.
    integer(int32), intent(in) :: n(3)
    integer,        intent(in) :: bcnd(0:n(1)+2,0:n(2)+2,0:n(3)+2)
    integer,        intent(in) :: npl

    integer :: ix, iy, p, ierr, last, nx1, ny1, nz1, ntype

    nx1 = n(1)+1; ny1 = n(2)+1; nz1 = n(3)+1
    ntype = size(win,4)
    last = zs_size - 1

    if (zs_size == 1) then
      do p = 1, ntype
        do iy = 1, ny1
          do ix = 1, nx1
            if (bcnd(ix,iy,1) == 0) win(ix,iy,1,p) = 0.5_real64*(win(ix,iy,1,p) + win(ix,iy,nz1,p))
          end do
        end do
        do iy = 1, ny1
          do ix = 1, nx1
            if (bcnd(ix,iy,nz1) == 0) win(ix,iy,nz1,p) = win(ix,iy,1,p)
          end do
        end do
      end do
      return
    end if

    if (zs_rank == last) then
      do p = 1, ntype
        sbuf(:,:,p) = win(:,:,nz1,p)
      end do
      call MPI_Send(sbuf, npl, MPI_DOUBLE_PRECISION, 0, 704, zs_comm, ierr)
      call MPI_Recv(rbuf, npl, MPI_DOUBLE_PRECISION, 0, 705, zs_comm, MPI_STATUS_IGNORE, ierr)
      do p = 1, ntype
        do iy = 1, ny1
          do ix = 1, nx1
            if (bcnd(ix,iy,nz1) == 0) win(ix,iy,nz1,p) = rbuf(ix,iy,p)
          end do
        end do
      end do
    else if (zs_rank == 0) then
      call MPI_Recv(rbuf, npl, MPI_DOUBLE_PRECISION, last, 704, zs_comm, MPI_STATUS_IGNORE, ierr)
      do p = 1, ntype
        do iy = 1, ny1
          do ix = 1, nx1
            if (bcnd(ix,iy,1) == 0) win(ix,iy,1,p) = 0.5_real64*(win(ix,iy,1,p) + rbuf(ix,iy,p))
          end do
        end do
        sbuf(:,:,p) = win(:,:,1,p)
      end do
      call MPI_Send(sbuf, npl, MPI_DOUBLE_PRECISION, last, 705, zs_comm, ierr)
    end if
  end subroutine stitch_z


  subroutine detect_zstitch(n, bcnd)
    integer(int32), intent(in) :: n(3)
    integer,        intent(in) :: bcnd(0:n(1)+2,0:n(2)+2,0:n(3)+2)

    zstitch_state = 0
    if (any(bcnd(1:n(1)+1,1:n(2)+1,1) == 0) .or. &
        any(bcnd(1:n(1)+1,1:n(2)+1,n(3)+1) == 0)) zstitch_state = 1
  end subroutine detect_zstitch


  subroutine ensure_buffers(n, ntype)
    integer(int32), intent(in) :: n(3)
    integer,        intent(in) :: ntype

    if (allocated(win)) then
      if (size(win,4) == ntype .and. lbound(win,3) == zs_k0 .and. &
          ubound(win,3) == zs_nhi .and. ubound(win,1) == n(1)+2) return
      deallocate(win, sbuf, rbuf)
    end if
    allocate(win(0:n(1)+2, 0:n(2)+2, zs_k0:zs_nhi, ntype))
    allocate(sbuf(0:n(1)+2, 0:n(2)+2, ntype), rbuf(0:n(1)+2, 0:n(2)+2, ntype))
    win = 0.0_real64
  end subroutine ensure_buffers


  subroutine build_rho_from_np(n, np_red, charge, ntype, rho, bcnd, flag_pbc)
    ! rho from a globally valid (synced) np_red, on the whole grid. Used
    ! once at start-up (build_initial_fields).
    integer(int32), intent(in)  :: n(3)
    integer,        intent(in)  :: ntype
    real(real64),   intent(in)  :: np_red(0:n(1)+2,0:n(2)+2,0:n(3)+2,ntype)
    real(real64),   intent(in)  :: charge(ntype)
    real(real64),   intent(out) :: rho(0:n(1)+2,0:n(2)+2,0:n(3)+2)
    integer(int32), intent(in)  :: bcnd(0:n(1)+2,0:n(2)+2,0:n(3)+2)
    integer(int32), intent(in)  :: flag_pbc

    integer :: ix, iy, iz, ptype
    real(real64) :: r

    rho = 0.0_real64

    !$omp parallel do collapse(3) private(ptype,r) schedule(static) default(shared)
    do iz = 1, n(3)+1
      do iy = 1, n(2)+1
        do ix = 1, n(1)+1
          r = 0.0_real64
          do ptype = 1, ntype
            r = r - charge(ptype) * np_red(ix,iy,iz,ptype)
          end do
          rho(ix,iy,iz) = r
        end do
      end do
    end do
    !$omp end parallel do

  end subroutine build_rho_from_np


  subroutine density_max_per_species(n, bcnd, np_red, ntype, np_mx)
    integer(int32), intent(in)  :: n(3)
    integer,        intent(in)  :: ntype
    integer,        intent(in)  :: bcnd(0:n(1)+2,0:n(2)+2,0:n(3)+2)
    real(real64),   intent(in)  :: np_red(0:n(1)+2,0:n(2)+2,0:n(3)+2,ntype)
    real(real64),   intent(out) :: np_mx(ntype)

    integer :: ix, iy, iz, ptype
    real(real64) :: local_max

    np_mx = 0.0_real64

    do ptype = 1, ntype
      local_max = 0.0_real64

      !$omp parallel do collapse(3) reduction(max:local_max) schedule(static) default(shared)
      do iz = 1, n(3)+1
        do iy = 1, n(2)+1
          do ix = 1, n(1)+1
            local_max = max(local_max, np_red(ix,iy,iz,ptype))
          end do
        end do
      end do
      !$omp end parallel do

      np_mx(ptype) = local_max
    end do

  end subroutine density_max_per_species


  ! Spatial average of np_red(:,:,:,ptype) over the physical node grid.
  ! np_red is populated on nodes 1..n(dim)+1 in each dimension; when a
  ! dimension is periodic, the density stitching makes node n(dim)+1
  ! a duplicate of node 1, so that node is excluded here to avoid double
  ! counting the seam plane.
  real(real64) function average_species_density(n, np_red, ntype, ptype, flag_pbc, flag_pbcz)
    integer(int32), intent(in) :: n(3)
    integer,        intent(in) :: ntype, ptype
    real(real64),   intent(in) :: np_red(0:n(1)+2,0:n(2)+2,0:n(3)+2,ntype)
    integer(int32), intent(in) :: flag_pbc, flag_pbcz

    integer :: ny_uniq, nz_uniq

    ny_uniq = merge(n(2), n(2)+1, flag_pbc  == 1_int32)
    nz_uniq = merge(n(3), n(3)+1, flag_pbcz == 1_int32)

    average_species_density = sum(np_red(1:n(1)+1, 1:ny_uniq, 1:nz_uniq, ptype)) / &
        real((n(1)+1) * ny_uniq * nz_uniq, real64)
  end function average_species_density

end module mod_density
