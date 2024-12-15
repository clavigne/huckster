! Parameters that are used in calculations
module params
  implicit none
  character(len=4), parameter :: calc_huckel = "eht", calc_promolecule = "pro", calc_skip = "skip"

  type Huckel
    double precision :: K = 1.75
  end type

  type CritEngine
    ! Parameters for Newton-Raphson search
    ! defaults from Rodriguez, Journal of Computational Chemistry 34, 681686 (2013)

    integer :: max_neval = 150 ! maximum number of evaluations for any critical point search
    double precision :: nr_mixing = 0.30d0 ! Newton-Raphson mixing parameter
    double precision :: tol_root = 1d-5 ! acceptable magnitude of |rho'| for a critical point
    double precision :: tol_rank = 1d-5 ! magnitude of an eigenvector of rho'' to be considered "zero" for rank calculations
    double precision :: tol_density_abandon = 1d-3
    double precision :: tol_grad_abandon = 1d-3
    double precision :: atom_density_abandon = 0.40d0 ! give up search if "atom density" is below this number
    double precision :: rtrust = 0.27d0 ! give up search if another critical point is closer than this distance

    ! Driver parameter
    double precision :: search_bond_distance = 6.0 ! search for a CP between every pair of atoms closer than this distance
    ! double precision :: search_ring_distance = 10.0 for ring search, not implemented
    ! double precision :: search_cage_distance = 14.0 for cage search, not implemented
    double precision :: search_grid_spacing = 0.5 ! search for a cp on a grid with this spacing

    ! Parameters for bond searches
    integer ::          path_max_segments = 10 ! maximum number of segments in a path search
    double precision :: path_segment_time = 100.0 ! effective "time" for path propagation
    double precision :: path_dopri5_rtol = 1d-6 ! relative tolerance for dopri5
    double precision :: path_dopri5_atol = 1d-4 ! absolute tolerance for dopri5
    double precision :: path_nudge_factor = 1d-3 ! how far to nudge the initial point from a CP before starting path propagation
  end type CritEngine

  type Parameters
    ! Model
    character(len=4) :: calc_type = calc_huckel
    integer :: charge = 0

    ! Route
    logical :: find_critical_points = .false.
    logical :: connect_graph = .false.

    ! Parametrizations
    type(Huckel) :: huckel_model
    type(CritEngine) :: crit_engine
  end type

contains
  subroutine params_export(ifile, p)
    integer, intent(in) :: ifile
    type(Parameters), intent(in) :: p

    namelist /PARAMS/ p

    write(ifile, PARAMS)
  end subroutine

  subroutine params_import(ifile, p)
    integer, intent(in) :: ifile
    type(Parameters), intent(inout) :: p

    namelist /PARAMS/ p

    read(ifile, PARAMS)
  end subroutine  

end module
