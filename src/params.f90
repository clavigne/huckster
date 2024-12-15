! Parameters that are used in calculations
module params
  use log
  implicit none
  integer, parameter :: io_read = 1, io_write = 2, io_update = 3
  character(len=4), parameter :: calc_huckel = "eht", calc_promolecule = "pro", calc_skip = "skip"

  type Model
    double precision :: Huckel_K = 1.75
  end type Model

  ! All defaults are from Rodriguez, Journal of Computational Chemistry 34, 681686 (2013)

  type CritNR
    ! Parameters for Newton-Raphson search
    integer :: max_neval = 150 ! maximum number of evaluations for any critical point search
    double precision :: mixing = 0.30d0 ! Newton-Raphson mixing parameter
    double precision :: tol_root = 1d-5 ! acceptable magnitude of |rho'| for a critical point
    double precision :: tol_rank = 1d-5 ! magnitude of an eigenvector of rho'' to be considered "zero" for rank calculations

    ! when the density is lower than tol_density_abandon and any of the gradients are higher than tol_grad_abandon, we abandon a search
    double precision :: tol_density_abandon = 1d-3
    double precision :: tol_grad_abandon = 1d-3

    double precision :: min_atom_density = 0.40d0 ! give up search if "atom density" is below this number
    double precision :: rtrust = 0.27d0 ! give up search if another critical point is closer than this distance
  end type CritNR

  type CritSearch
    logical :: search_atoms = .true. ! do a search at every atomic center
    logical :: search_bonds = .false. ! do a search between every pairs of atoms
    logical :: search_grid = .true. ! do a search over a grid around the molecule
    double precision :: max_bond = 6.0 ! search for a CP between every pair of atoms closer than this distance
    ! double precision :: max_ring = 10.0 for ring search, not implemented
    ! double precision :: max_cage = 14.0 for cage search, not implemented
    double precision :: grid = 0.5 ! search for a cp on a grid with this spacing
  end type CritSearch

  type CritPath
    integer ::          max_segments = 10 ! maximum number of segments in a path search
    double precision :: time = 100.0 ! effective "time" for path propagation
    double precision :: dopri5_rtol = 1d-6 ! relative tolerance for dopri5
    double precision :: dopri5_atol = 1d-4 ! absolute tolerance for dopri5
    double precision :: nudge = 1d-3 ! how far to nudge the initial point from a CP before starting path propagation
  end type CritPath

  type Route
    ! Computation model
    character(len=4) :: calc_type = calc_huckel
    integer :: charge = 0
    logical :: find_crits = .false.
    logical :: path_crits = .false.
  end type

  type Parameters
    type(Route) :: route
    type(Model) :: model
    type(CritNR) :: nr
    type(CritSearch) :: search
    type(CritPath) :: path
  end type

contains
  subroutine params_export(ifile, p)
    integer, intent(in) :: ifile
    type(Parameters), intent(in) :: p

    namelist /PARAMS/ p

    write (ifile, PARAMS)
  end subroutine

  subroutine params_print_extended_options()
    type(Parameters) :: defaults
    namelist /PARAMS/ defaults
    call params_io(6, defaults, io_write)
  end subroutine

  subroutine params_import(ifile, p)
    integer, intent(in) :: ifile
    type(Parameters), intent(inout) :: p

    namelist /PARAMS/ p

    read (ifile, PARAMS)
  end subroutine

  subroutine params_io(unit, p, io, iunit)
    character(len=*), intent(in), optional :: iunit
    integer, intent(in) :: unit
    type(Parameters), intent(inout) :: p
    integer, intent(in) :: io

    integer :: info

    ! route
    character(len=4) :: calc_type
    integer :: charge
    logical :: find_crits, path_crits
    namelist /ROUTE/ calc_type, charge, find_crits, path_crits

    ! model
    double precision :: Huckel_K
    namelist /MODEL/ Huckel_K

    !nr
    integer :: max_neval
    double precision :: mixing, tol_root, tol_rank, tol_density_abandon, tol_grad_abandon
    double precision :: min_atom_density, rtrust
    namelist /critnr/ max_neval, mixing, tol_root, tol_rank, tol_density_abandon, tol_grad_abandon, min_atom_density, rtrust

    !search
    logical :: search_atoms, search_bonds, search_grid
    double precision :: max_bond, grid
    namelist /critsearch/ search_atoms, search_bonds, search_grid, max_bond, grid

    !pathing
    integer :: max_segments
    double precision :: time, dopri5_rtol, dopri5_atol, nudge
    namelist /critpath/ max_segments, time, dopri5_rtol, dopri5_atol, nudge

    if (present(iunit) .and. (io .ne. io_read)) then
      call log_err("params_io", "illegal use of iunit with writes")
      error stop -1
    end if

    calc_type = p%route%calc_type
    charge = p%route%charge
    find_crits = p%route%find_crits
    path_crits = p%route%path_crits

    Huckel_K = p%model%Huckel_K

    max_neval = p%nr%max_neval
    mixing = p%nr%mixing
    tol_root = p%nr%tol_root
    tol_rank = p%nr%tol_rank
    tol_density_abandon = p%nr%tol_density_abandon
    tol_grad_abandon = p%nr%tol_grad_abandon
    min_atom_density = p%nr%min_atom_density
    rtrust = p%nr%rtrust

    search_atoms = p%search%search_atoms
    search_bonds = p%search%search_bonds
    search_grid = p%search%search_grid
    max_bond = p%search%max_bond
    grid = p%search%grid

    max_segments = p%path%max_segments
    time = p%path%time
    dopri5_rtol = p%path%dopri5_rtol
    dopri5_atol = p%path%dopri5_atol
    nudge = p%path%nudge

    if (io .eq. io_read .or. io .eq. io_update) then
      if (present(iunit))  then
        read (iunit, nml=route, iostat=info)
      else 
        read (unit, nml=route, iostat=info)
        rewind(unit)
      end if
      call check_read("ROUTE")

      if (present(iunit))  then
        read (iunit, nml=model, iostat=info)
      else 
        read (unit, nml=model, iostat=info)
        rewind(unit)
      end if      
      call check_read("MODEL")

      if (present(iunit))  then
        read (iunit, nml=critnr, iostat=info)
      else 
        read (unit, nml=critnr, iostat=info)
        rewind(unit)
      end if            
      call check_read("CRITNR")

      if (present(iunit))  then
        read (iunit, nml=critsearch, iostat=info)
      else 
        read (unit, nml=critsearch, iostat=info)
        rewind(unit)
      end if            
      call check_read("CRITSEARCH")

      if (present(iunit))  then
        read (iunit, nml=critpath, iostat=info)
      else 
        read (unit, nml=critpath, iostat=info)
        rewind(unit)
      end if            
      call check_read("CRITPATH")

      p%route%calc_type = calc_type
      p%route%charge = charge
      p%route%find_crits = find_crits
      p%route%path_crits = path_crits

      p%model%Huckel_K = Huckel_K

      p%nr%max_neval = max_neval
      p%nr%mixing = mixing
      p%nr%tol_root = tol_root
      p%nr%tol_rank = tol_rank
      p%nr%tol_density_abandon = tol_density_abandon
      p%nr%tol_grad_abandon = tol_grad_abandon
      p%nr%min_atom_density = min_atom_density
      p%nr%rtrust = rtrust

      p%search%search_atoms = search_atoms
      p%search%search_bonds = search_bonds
      p%search%search_grid = search_grid
      p%search%max_bond = max_bond
      p%search%grid = grid

      p%path%max_segments = max_segments
      p%path%time = time
      p%path%dopri5_rtol = dopri5_rtol
      p%path%dopri5_atol = dopri5_atol
      p%path%nudge = nudge
    end if

    if (io .eq. io_write .or. io .eq. io_update) then
      write (unit, nml=route, iostat=info)
      if (info .ne. 0) then
        call log_err('params_io', "Error while writing ROUTE")
        error stop - 1
      end if

      write (unit, nml=model, iostat=info)
      if (info .ne. 0) then
        call log_err('params_io', "Error while writing MODEL")
        error stop - 1
      end if

      write (unit, nml=critnr, iostat=info)
      if (info .ne. 0) then
        call log_err('params_io', "Error while writing CRITNR")
        error stop - 1
      end if

      write (unit, nml=critsearch, iostat=info)
      if (info .ne. 0) then
        call log_err('params_io', "Error while writing CRITSEARCH")
        error stop - 1
      end if

      write (unit, nml=critpath, iostat=info)
      if (info .ne. 0) then
        call log_err('params_io', "Error while writing CRITPATH")
        error stop - 1
      end if
    end if

  contains
    subroutine check_read(which)
      character(len=*), intent(in) :: which
      ! Note: we can't log here because we don't know yet if -q was passed or no

      if ((info .eq. 0) .or. (info .eq. -1)) then
        ! -1 == end of file = nml wasnt found, at least on gfortran
      else
        call log_err('params_io', "Error while reading "// which)
        error stop - 1
      end if

    end subroutine

  end subroutine

end module
