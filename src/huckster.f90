program huckster
  use integrals
  use constants
  use hamiltonians
  use ed
  use crits
  use log
  use params
  use models
  implicit none

  ! Command line arguments
  character(len=4096)           :: arg

  ! Parameter sets
  type(Parameters) :: ps
  character(len=4096) :: path
  character(len=:), allocatable :: output_file
  character(len=:), allocatable :: input
  character(len=:), allocatable :: output
  integer :: info, i, j

  ! Result from model
  type(UnrolledMOs) :: umos

  ! outputs
  integer, parameter :: iwfn = 10, icrit = 11, iout = 12

  ! CP output
  type(CritPoint), allocatable :: CPs(:)
  integer :: ncp
  integer, allocatable :: adjacency(:, :)

  ! Load command line arguments
  i = 1
  do while (i .le. command_argument_count())
    call get_command_argument(i, arg)
    select case (arg)

    case ('--')
      ! prepping to receive final argument
      exit

    case ('-h', '--help')
      call print_help()
      stop 0

    case ('-X', "--extended")
      if (command_argument_count() .eq. 1) then
        call params_print_extended_options()
        stop 0
      else
        i = i + 1
        if (i .ge. command_argument_count()) then
          call log_err('huckster', '-X needs an argument')
          error stop - 1
        end if

        call get_command_argument(i, arg)
        call params_io(2, ps, io_read, arg)
      end if

    case ('-x', "--extended-file")
      if (command_argument_count() .eq. 1) then
        call params_print_extended_options()
        stop 0
      else
        i = i + 1
        if (i .ge. command_argument_count()) then
          call log_err('huckster', '-x needs an argument')
          error stop - 1
        end if

        call get_command_argument(i, arg)
        open (unit=2, file=arg, action='READ', iostat=info)
        if (info .ne. 0) then
          call log_err('huckster', 'could not open extended options file '//trim(arg))
          error stop - 1
        end if
        call params_io(2, ps, io_read)
        close (2)
      end if

    case ('-q', '--quiet')
      call set_verbosity(-1)

    case ('-v', '--verbose')
      call set_verbosity(1)

    case ('-vv', '--very-verbose')
      call set_verbosity(2)

    case ('-a', '--aim')
      ps%route%find_crits = .true.

    case ('-g', '--graph')
      ps%route%find_crits = .true.        ! we need CPs for the graph
      ps%route%path_crits = .true.

    case ('-t', '--type')
      call get_command_argument(i + 1, arg)

      select case (arg)
      case ('pro', 'promolecule', 'sad')
        ps%route%calc_type = calc_promolecule

      case ('huckel', 'eht')
        ps%route%calc_type = calc_huckel

      case ('skip')
        ps%route%calc_type = calc_skip

      case default
        call log_err('huckster', trim(arg)//' is not a valid calculation type')
        error stop - 1
      end select
      i = i + 1

    case ('-c', '--charge')
      call get_command_argument(i + 1, arg)
      read (arg, *, iostat=info) ps%route%charge
      if (info .ne. 0) then
        call log_err('huckster', trim(arg)//' is an invalid charge.')
        error stop - 1
      else
        i = i + 1
      end if

    case default
      call log_err('huckster', 'unrecognized command-line option '//trim(arg))
    end select
    i = i + 1
  end do

  ! print out the program banner
  call log_banner

  if (very_verbose()) then
    write (*, *) ""
    call params_io(6, ps, io_write)
    write (*, *) ""
  end if

  if ((i .eq. command_argument_count()) .or. (command_argument_count() .eq. 0)) then
    call log_err('huckster', 'no geometry file passed')
    error stop - 1
  else
    i = i + 1
    call get_command_argument(i, path)
    input = trim(path)

    if (i < command_argument_count()) then
      call get_command_argument(i + 1, path)
      output_file = trim(path)
    else
      output_file = input
    end if
  end if

  ! remove extension
  associate (idot => index(output_file, '.', back=.true.))
    if (idot > 0) then
      output = output_file(1:idot - 1)
    else
      output = output_file(1:len_trim(output_file))
    end if
  end associate

  if (ps%route%calc_type .ne. calc_skip) then
    call models_from_params(input, ps, umos)

    call log_program_step('Saving wavefunction')
    open (iwfn, file=output//'.wfn', action='write', iostat=info)

    if (info .ne. 0) then
      call log_err('huckster', 'Could not open output wfn file')
      error stop - 1
    end if

    call log_program_substep('writing to file: '//output//'.wfn')
    call integrals_write_to_wfn(iwfn, umos)
    close (iwfn)

    call log_program_step_end
  else
    call log_program_step('Reading wavefunction from file')
    open (iwfn, file=output//'.wfn', action='read', iostat=info)

    if (info .ne. 0) then
      call log_err('huckster', 'Could not open wfn file')
      error stop - 1
    end if

    call integrals_read_wfn(iwfn, umos)
    close (iwfn)

    call log_program_step_end
  end if

  if (ps%route%find_crits) then
    call log_program_step('Constructing electron density calculator')
    call ed_initialize(umos)
    call crits_initialize(umos, ps)

    call log_program_step('Searching for critical points')

    if (ps%search%search_atoms) then
      call log_program_substep("at atom positions")
      call crits_do_atoms
    end if

    if (ps%search%search_bonds) then
      call log_program_substep("between pairs of atoms")
      call crits_do_bonds
    end if

    if (ps%search%search_grid) then
      call log_program_substep("over a grid")
      call crits_do_grid
    end if

    call crits_print

    call log_program_substep("Calculating electronic densities at found CP")
    call crits_get_CPs(CPs, ncp)
    call log_program_substep("Saving found critical points")
    if (logging()) write (*, *) '         writing CPs to csv file'

    open (unit=iout, file=output//'.csv')
    write (iout, *) 'index,atom,rank,curv,x,y,z,rho,ellip_x,ellip_y,ellip_z'
    do i = 1, ncp
      write (iout, '(i4,a)', advance='no') i, ','
      if (CPs(i)%atom_id > 0) then
        write (iout, '(a,a)', advance='no') trim(umos%atoms(CPs(i)%atom_id)), ','
      else
        write (iout, '(a)', advance='no') ','
      end if
      write (iout, '(i2,a)', advance='no') CPs(i)%rank, ','
      write (iout, '(i2,a)', advance='no') CPs(i)%curvature, ','
      write (iout, '(E23.15,a)', advance='no') CPs(i)%position(1), ','
      write (iout, '(E23.15,a)', advance='no') CPs(i)%position(2), ','
      write (iout, '(E23.15,a)', advance='no') CPs(i)%position(3), ','
      write (iout, '(E23.15,a)', advance='no') CPs(i)%electron_density, ','
      write (iout, '(E23.15,a)', advance='no') CPs(i)%ellipticity(1), ','
      write (iout, '(E23.15,a)', advance='no') CPs(i)%ellipticity(2), ','
      write (iout, '(E23.15)', advance='no') CPs(i)%ellipticity(3)
      write (iout, *) ''
    end do
    close (unit=iout)

    call log_program_step_end
  end if

  if (ps%route%path_crits) then
    call log_program_step('Connecting critical points into graph')
    neval = 0
    call crits_perceive_graph(adjacency)
    call log_program_step_end
    if (logging()) write (*, *) '         writing adjacency matrix to .mat file'

    open (unit=iout, file=output//'.mat')
    do i = 1, ncp
      do j = 1, ncp
        write (iout, '(i2)', advance='no') adjacency(j, i)
      end do
      write (iout, *) ''
    end do
    close (unit=iout)
    write (*, *) ''
  end if

  if (logging()) write (*, *) '~~~~~END~~~~END~~~~END~~~~END~~~~END~~~~END~~~~END~~~~~~'
  stop 0
contains
  subroutine print_help
    write (*, *) 'USAGE: huckster [options] -- GEOMETRY [OUTPUT]'
    write (*, *) ''
    write (*, *) 'Huckster is a standalone program to perform extended Huckel theory calculations,'
    write (*, *) 'specifically to use in QTAIM calculations. It takes a single mandatory argument,'
    write (*, *) 'a molecular geometry in the XMOL format in file GEOMETRY and outputs a AIM-format'
    write (*, *) 'wavefunction file to OUTPUT (defaults to GEOMETRY with a .wfn extension).'

    write (*, *) 'Standard command-line options'
    write (*, *) '   -c, --charge CHRG          Set the global charge of the computed molecule.'
    write (*, *) '   -t, --type TYPE            Define computation type. Valid choices are:'
    write (*, *) '                               eht  => extended Huckel theory'
    write (*, *) '                               pro  => promolecule electron density'
    write (*, *) '                               skip => read wfn file at GEOMETRY.wfn'
    write (*, *) '   -a, --aim                  Do a QTAIM search for critical points.'
    write (*, *) '   -g, --graph                Build molecular graph from QTAIM analysis.'
    write (*, *) ''
    write (*, *) '   -x, --extended-file OPTS   Read extended options from file OPTS.'
    write (*, *) '   -X, --extended NML         Read extended options as an argument.'
    write (*, *) ''
    write (*, *) '   -h, --help                 Print these instructions.'
    write (*, *) '   -v, --verbose              Print more information.'
    write (*, *) '   -vv, --very-verbose        Print more information.'
    write (*, *) '   -q, --quiet                Silence output'
    write (*, *) ''
    write (*, *) 'Extended options'
    write (*, *) '    Additional options can be passed using Fortran "namelist" format using the'
    write (*, *) '    -x command-line option. A list of those options and their default values is'
    write (*, *) '    printed to STDOUT by running `huckster -x`. Any of these options can be modified'
    write (*, *) '    either in the OPTS file passed to the -x option.'
    write (*, *) ''
    write (*, *) '    Exended options can also be passed directly as argumens using -X,'
    write (*, *) ''
    write (*, *) "       huckster -X '&critsearch search_bonds=t, search_grid=f /' -a -- h2o.xyz"

  end subroutine print_help

end program huckster
