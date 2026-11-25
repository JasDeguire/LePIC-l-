module mod_zslab
  !===========================================================================
  ! Particle domain decomposition along z, shared by the whole code.
  !
  ! Each MPI rank owns the SAME z-slab as the Poisson decomposition
  ! (mod_PoissonDecomposition): with m = nz/nranks and k0 = rank*m,
  !
  !   owned cell planes   iz = k0+1 .. k0+m          (zs_cz_lo .. zs_cz_hi)
  !   owned cell ids      ic = nx*ny*k0+1 .. nx*ny*(k0+m)  (zs_c_lo .. zs_c_hi)
  !   deposit node window iz = k0+1 .. k0+m+1        (zs_nlo .. zs_nhi)
  !
  ! A rank holds exactly the particles whose cell plane iz = int(z/hz)+1
  ! (clamped to 1..nz) lies in its owned range, so all per-thread
  ! grid-sized work arrays (np_thread, cell_count/cell_start, collision
  ! cell tables) only cover that slab: memory per rank is divided by the
  ! number of ranks instead of being a full copy of the grid per thread.
  !
  ! Arrays local to the slab are allocated with GLOBAL index bounds
  ! (e.g. cell_count(zs_c_lo:zs_c_hi), np_thread(:,:,zs_nlo:zs_nhi,...)),
  ! so every routine keeps indexing with global cell/plane numbers.
  !
  ! Fields (phi, E, bcnd, kq, B) stay replicated on every rank (phase 1).
  !===========================================================================
  use iso_fortran_env, only: int32, real64
  use mpi
  implicit none
  private

  public :: zslab_init, zslab_owner_of_iz, zslab_iz_of_z, zslab_next_hop

  logical,        public, protected :: zs_ready = .false.
  integer,        public, protected :: zs_comm  = MPI_COMM_NULL
  integer,        public, protected :: zs_rank  = 0
  integer,        public, protected :: zs_size  = 1
  ! Neighbours in the z chain. zs_left/zs_right are MPI_PROC_NULL at the
  ! ends of a non-periodic chain; with z-periodic particles (flag_pbc==1,
  ! positions wrap in z) the chain is closed into a ring for migration.
  integer,        public, protected :: zs_left  = MPI_PROC_NULL
  integer,        public, protected :: zs_right = MPI_PROC_NULL
  ! Non-periodic neighbours (always open chain): used for the density
  ! halo exchange, where the z-periodic stitching is done separately.
  integer,        public, protected :: zs_left_open  = MPI_PROC_NULL
  integer,        public, protected :: zs_right_open = MPI_PROC_NULL
  logical,        public, protected :: zs_periodic = .false.

  integer(int32), public, protected :: zs_nx = 0, zs_ny = 0, zs_nz = 0
  integer(int32), public, protected :: zs_m  = 0, zs_k0 = 0
  integer(int32), public, protected :: zs_cz_lo = 1, zs_cz_hi = 0
  integer(int32), public, protected :: zs_nlo = 1, zs_nhi = 0
  integer(int32), public, protected :: zs_c_lo = 1, zs_c_hi = 0
  real(real64),   public, protected :: zs_hz = 1.0_real64

contains

  subroutine zslab_init(comm, n, h, k0, m, periodic_z)
    integer,        intent(in) :: comm
    integer(int32), intent(in) :: n(3)
    real(real64),   intent(in) :: h(3)
    integer(int32), intent(in) :: k0, m      ! from PoissonDecomp (same split)
    logical,        intent(in) :: periodic_z
    integer :: ierr

    zs_comm = comm
    call MPI_Comm_rank(comm, zs_rank, ierr)
    call MPI_Comm_size(comm, zs_size, ierr)

    zs_nx = n(1); zs_ny = n(2); zs_nz = n(3)
    zs_hz = h(3)
    zs_k0 = k0
    zs_m  = m
    zs_periodic = periodic_z

    zs_cz_lo = k0 + 1_int32
    zs_cz_hi = k0 + m
    zs_nlo   = zs_cz_lo
    zs_nhi   = zs_cz_hi + 1_int32
    zs_c_lo  = n(1)*n(2)*k0 + 1_int32
    zs_c_hi  = n(1)*n(2)*(k0 + m)

    zs_left_open  = merge(zs_rank-1, MPI_PROC_NULL, zs_rank > 0)
    zs_right_open = merge(zs_rank+1, MPI_PROC_NULL, zs_rank < zs_size-1)
    zs_left  = zs_left_open
    zs_right = zs_right_open
    if (periodic_z .and. zs_size > 1) then
      zs_left  = modulo(zs_rank-1, zs_size)
      zs_right = modulo(zs_rank+1, zs_size)
    end if

    zs_ready = .true.
  end subroutine zslab_init


  pure integer(int32) function zslab_iz_of_z(z) result(iz)
    ! Cell plane of a z position, same convention (and clamp) as
    ! cell_index_from_position in mod_particle_sorting.
    real(real64), intent(in) :: z
    iz = int(z / zs_hz, int32) + 1_int32
    iz = max(1_int32, min(zs_nz, iz))
  end function zslab_iz_of_z


  pure integer function zslab_owner_of_iz(iz) result(r)
    integer(int32), intent(in) :: iz
    r = int((iz - 1_int32) / zs_m)
    r = max(0, min(zs_size-1, r))
  end function zslab_owner_of_iz


  pure integer function zslab_next_hop(owner) result(dest)
    ! Neighbour to forward a particle to so that it gets closer to its
    ! owner rank (shortest way round the ring when z is periodic).
    integer, intent(in) :: owner
    integer :: d

    if (owner == zs_rank) then
      dest = zs_rank
      return
    end if
    if (zs_periodic) then
      d = modulo(owner - zs_rank, zs_size)
      if (2*d <= zs_size) then
        dest = zs_right
      else
        dest = zs_left
      end if
    else
      dest = merge(zs_right, zs_left, owner > zs_rank)
    end if
  end function zslab_next_hop

end module mod_zslab
