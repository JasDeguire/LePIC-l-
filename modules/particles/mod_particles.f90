module mod_particles
  use iso_fortran_env, only: real64, int32, int8
  use mod_zslab,       only: zs_c_lo, zs_c_hi
  implicit none
  private
  public :: ParticleSet

  type :: ParticleSet
    integer(int32) :: n    = 0
    integer(int32) :: nmax = 0

    ! Particle kinematic data, packed into one flat (6,nmax) array instead
    ! of 6 separate allocatable arrays (x,y,z,vx,vy,vz individually) - one
    ! array descriptor for all 6 components instead of 6, and one particle's
    ! full state is contiguous (matches legacy's vxp(6,nmax,...) layout and
    ! its own "1==x, 2==y, 3==z, 4==vx, 5==vy, 6==vz" indexing convention,
    ! Src/part_expmover.f90). Compiler vectorization-report comparison
    ! (mod_particleMover.f90's move_and_bc_boris vs legacy's part_mover)
    ! showed the per-field allocatable layout as the one concrete structural
    ! difference between the two movers' hot per-particle loop once the
    ! B-field negligible-skip and vectorization angles were both ruled out.
    real(real64), allocatable :: pv(:,:)  ! (6, nmax): 1=x 2=y 3=z 4=vx 5=vy 6=vz
    real(real64), allocatable :: w(:)
    integer(int32), allocatable :: sp(:)

    ! Sorting / binning metadata
    integer(int32), allocatable :: cell_id(:)     ! size nmax
    integer(int32), allocatable :: cell_count(:)  ! size ncells
    integer(int32), allocatable :: cell_start(:)  ! size ncells
    integer(int32) :: ncells = 0
    ! cell_count/cell_start are zero outside cells cell_lo..cell_hi (empty
    ! when cell_lo > cell_hi). Lets the sort and the collision setup touch
    ! only the cell range a set actually occupies - after the per-species
    ! redistribution that is one z-slab - instead of all ncells.
    integer(int32) :: cell_lo = 1, cell_hi = 0

    ! Collision-related per-particle flags
    integer(int8),  allocatable :: flag_dead(:)   ! size nmax
    integer(int32), allocatable :: flag_cex(:)    ! size nmax

    ! z-slab migration (mod_particleMigration): the movers record here the
    ! indices of the survivors that left this rank's z-slab, in increasing
    ! order, and set mig_valid; mig_n_checked is part%n right after that
    ! push (particles appended later, e.g. secondary electrons, are not in
    ! the list). The migration then only visits these particles instead of
    ! scanning the whole set, and clears mig_valid.
    integer(int32), allocatable :: mig_idx(:)
    integer(int32) :: n_mig = 0
    integer(int32) :: mig_n_checked = 0
    logical        :: mig_valid = .false.

  contains
    procedure :: allocate_pset
    procedure :: ensure_capacity
    procedure :: ensure_cell_storage
    procedure :: push_mig
    procedure :: from_vxp
    procedure :: clear
    procedure :: destroy
  end type ParticleSet


contains

  subroutine allocate_pset(self, nmax_in)
    class(ParticleSet), intent(inout) :: self
    integer(int32),     intent(in)    :: nmax_in

    if (allocated(self%pv)) call self%destroy()

    self%n    = 0_int32
    self%nmax = max(0_int32, nmax_in)
    self%ncells = 0_int32

    if (self%nmax <= 0_int32) return

    allocate(self%pv(6, self%nmax))
    allocate(self%w(self%nmax))
    allocate(self%sp(self%nmax))

    allocate(self%cell_id(self%nmax))
    allocate(self%flag_dead(self%nmax))
    allocate(self%flag_cex(self%nmax))

    self%pv = 0.0_real64
    self%w  = 0.0_real64
    self%sp = 0_int32

    self%cell_id   = 0_int32
    self%flag_dead = 0_int8
    self%flag_cex  = 0_int32
  end subroutine allocate_pset

  subroutine ensure_capacity(self, needed)
    class(ParticleSet), intent(inout) :: self
    integer(int32),     intent(in)    :: needed

    integer(int32) :: new_nmax, old_nmax, ncopy
    real(real64), allocatable :: pv_new(:,:)
    real(real64), allocatable :: w_new(:)
    integer(int32), allocatable :: sp_new(:), cell_id_new(:), flag_cex_new(:)
    integer(int8),  allocatable :: flag_dead_new(:)

    if (needed <= self%nmax) return

    old_nmax = self%nmax
    new_nmax = max(needed, max(1_int32, 2_int32*old_nmax))
    ncopy    = self%n

    allocate(pv_new(6, new_nmax))
    allocate(w_new(new_nmax))
    allocate(sp_new(new_nmax))
    allocate(cell_id_new(new_nmax))
    allocate(flag_dead_new(new_nmax))
    allocate(flag_cex_new(new_nmax))

    pv_new = 0.0_real64
    w_new  = 0.0_real64
    sp_new = 0_int32
    cell_id_new   = 0_int32
    flag_dead_new = 0_int8
    flag_cex_new  = 0_int32

    if (old_nmax > 0_int32) then
      if (allocated(self%pv)) pv_new(:,1:ncopy) = self%pv(:,1:ncopy)
      if (allocated(self%w))  w_new(1:ncopy) = self%w(1:ncopy)
      if (allocated(self%sp)) sp_new(1:ncopy) = self%sp(1:ncopy)

      if (allocated(self%cell_id))   cell_id_new(1:ncopy)   = self%cell_id(1:ncopy)
      if (allocated(self%flag_dead)) flag_dead_new(1:ncopy) = self%flag_dead(1:ncopy)
      if (allocated(self%flag_cex))  flag_cex_new(1:ncopy)  = self%flag_cex(1:ncopy)

      if (allocated(self%pv))        deallocate(self%pv)
      if (allocated(self%w))         deallocate(self%w)
      if (allocated(self%sp))        deallocate(self%sp)
      if (allocated(self%cell_id))   deallocate(self%cell_id)
      if (allocated(self%flag_dead)) deallocate(self%flag_dead)
      if (allocated(self%flag_cex))  deallocate(self%flag_cex)
    end if

    call move_alloc(pv_new, self%pv)
    call move_alloc(w_new, self%w)
    call move_alloc(sp_new, self%sp)
    call move_alloc(cell_id_new, self%cell_id)
    call move_alloc(flag_dead_new, self%flag_dead)
    call move_alloc(flag_cex_new, self%flag_cex)

    self%nmax = new_nmax
  end subroutine ensure_capacity


  subroutine ensure_cell_storage(self)
    ! cell_count/cell_start cover only this MPI rank's z-slab of cells
    ! (mod_zslab), allocated with their GLOBAL cell ids as bounds
    ! (zs_c_lo:zs_c_hi) so callers keep indexing by global cell id.
    class(ParticleSet), intent(inout) :: self
    integer(int32) :: ncells_in

    ncells_in = zs_c_hi - zs_c_lo + 1_int32
    if (self%ncells == ncells_in .and. allocated(self%cell_count)) then
      if (lbound(self%cell_count,1) == zs_c_lo) return
    end if

    if (allocated(self%cell_count)) deallocate(self%cell_count)
    if (allocated(self%cell_start)) deallocate(self%cell_start)

    self%ncells = max(0_int32, ncells_in)

    if (self%ncells <= 0_int32) return

    allocate(self%cell_count(zs_c_lo:zs_c_hi))
    allocate(self%cell_start(zs_c_lo:zs_c_hi))

    self%cell_count = 0_int32
    self%cell_start = 0_int32
    self%cell_lo = zs_c_lo
    self%cell_hi = zs_c_lo - 1_int32
  end subroutine ensure_cell_storage


  subroutine from_vxp(self, vxp_in, npar, species_id)
    class(ParticleSet), intent(inout) :: self
    real(real64),       intent(in)    :: vxp_in(:,:)
    integer(int32),     intent(in)    :: npar
    integer(int32),     intent(in)    :: species_id

    integer(int32) :: i

    if (size(vxp_in,1) /= 6) then
      error stop 'from_vxp: first dimension of vxp_in must be 6'
    end if

    if (size(vxp_in,2) < npar) then
      error stop 'from_vxp: second dimension of vxp_in is smaller than npar'
    end if

    call self%ensure_capacity(npar)

    self%n = npar

    do i = 1, npar
      self%pv(1,i) = vxp_in(1,i)
      self%pv(2,i) = vxp_in(2,i)
      self%pv(3,i) = vxp_in(3,i)
      self%pv(4,i) = vxp_in(4,i)
      self%pv(5,i) = vxp_in(5,i)
      self%pv(6,i) = vxp_in(6,i)
      self%w(i)  = 1.0_real64
      self%sp(i) = species_id
    end do

    self%cell_id(1:npar)   = 0_int32
    self%flag_dead(1:npar) = 0_int8
    self%flag_cex(1:npar)  = 0_int32
  end subroutine from_vxp

  subroutine clear(self)
    class(ParticleSet), intent(inout) :: self

    self%n = 0_int32

    if (allocated(self%cell_id))   self%cell_id   = 0_int32
    if (allocated(self%flag_dead)) self%flag_dead = 0_int8
    if (allocated(self%flag_cex))  self%flag_cex  = 0_int32

    if (allocated(self%cell_count)) self%cell_count = 0_int32
    if (allocated(self%cell_start)) self%cell_start = 0_int32
    self%cell_lo = 1_int32
    self%cell_hi = 0_int32
  end subroutine clear


  subroutine push_mig(self, i)
    ! Append index i to the migration candidate list (rare: only called
    ! for particles leaving the rank's slab).
    class(ParticleSet), intent(inout) :: self
    integer(int32),     intent(in)    :: i
    integer(int32), allocatable :: tmp(:)

    if (.not. allocated(self%mig_idx)) then
      allocate(self%mig_idx(1024))
    else if (self%n_mig >= size(self%mig_idx)) then
      allocate(tmp(2*size(self%mig_idx)))
      tmp(1:self%n_mig) = self%mig_idx(1:self%n_mig)
      call move_alloc(tmp, self%mig_idx)
    end if
    self%n_mig = self%n_mig + 1_int32
    self%mig_idx(self%n_mig) = i
  end subroutine push_mig


  subroutine destroy(self)
    class(ParticleSet), intent(inout) :: self

    if (allocated(self%pv))        deallocate(self%pv)
    if (allocated(self%w))         deallocate(self%w)
    if (allocated(self%sp))        deallocate(self%sp)

    if (allocated(self%cell_id))   deallocate(self%cell_id)
    if (allocated(self%cell_count)) deallocate(self%cell_count)
    if (allocated(self%cell_start)) deallocate(self%cell_start)

    if (allocated(self%flag_dead)) deallocate(self%flag_dead)
    if (allocated(self%flag_cex))  deallocate(self%flag_cex)
    if (allocated(self%mig_idx))   deallocate(self%mig_idx)
    self%n_mig = 0
    self%mig_valid = .false.

    self%n      = 0_int32
    self%nmax   = 0_int32
    self%ncells = 0_int32
  end subroutine destroy

end module mod_particles
