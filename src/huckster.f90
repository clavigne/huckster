program huckster
  use integrals
  use constants
  use hamiltonians
  use ed
  use crits
  use log
  use params
  implicit none

  ! homo lumo
  character(len=4)::homolumo(3) = (/'    ', 'HOMO', 'LUMO'/)
  integer :: hl_index
  integer :: HOMO
  integer :: LUMO

  ! Command line arguments
  character(len=4096)           :: arg

  ! Parameter sets
  type(Parameters) :: ps
  character(len=4096) :: path
  character(len=:), allocatable :: output_file
  character(len=:), allocatable :: input
  character(len=:), allocatable :: output
  integer :: info

  ! Integrals
  type(ElectronicSystem) :: electrons
  type(UnrolledMOs) :: umos

  ! Work parameters
  integer :: i, ii, j, naos, remaining_electrons

  ! Hamiltonian / overlap
  double precision, allocatable :: S_prim(:, :)
  double precision, allocatable :: H(:, :)
  double precision, allocatable :: S(:, :)
  double precision, allocatable :: MO(:, :)
  double precision, allocatable :: MO_AO(:, :)
  double precision, allocatable :: E_MO(:)
  double precision, allocatable :: occ(:)

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
    ! ---------------------------------------------------------------------------------
    ! Initialize integral module
    call integrals_initialize

    ! ---------------------------------------------------------------------------------
    ! Read input geometry
    call log_program_step('Reading geometry file')
    if (logging()) write (*, *) '     file: ', input

    open (unit=2, file=input, action='READ', iostat=info)
    if (info .ne. 0) then
      call log_err('huckster', 'could not open input geometry file: '//input)
      error stop - 1
    end if
    call integrals_init_from_file(2, electrons, ps%route%charge)
    close (unit=2)
    call log_program_step_end

    ! ---------------------------------------------------------------------------------
    ! Generate the basis
    call log_program_step('Loading primitive basis')
    call integrals_build_basis(electrons)
    call log_program_step_end

    call log_program_step('Loading AO coeffs and energies')
    call integrals_build_atomic_orbitals(electrons)
    naos = electrons%naos
    call log_program_step_end

    if (ps%route%calc_type .eq. calc_huckel) then
      ! ---------------------------------------------------------------------------------
      ! Setup the Hamiltonian
      call log_program_step('Building Extended Huckel Hamiltonian')

      call log_program_substep('generating overlaps in the gaussian basis')
      call integrals_overlaps(electrons, S_prim)

      call log_program_substep('projecting S_orb to opt. AOs')
      call integrals_symm_bas2ao(electrons, S_prim, S)

      call log_program_substep('generating Hamiltonian elements')
      allocate (H(naos, naos))
      do i = 1, naos
        H(i, i) = electrons%E_AO(i)
        do j = i + 1, naos
          ! multiply off diag by Ei + Ej / 2
          H(i, j) = S(i, j)*(electrons%E_AO(i) + electrons%E_AO(j))*0.5*ps%model%Huckel_K
        end do
      end do
      call log_program_step_end

      call log_program_step('Eigendecomposition of H')
      if (logging()) write (*, '(a,i5,a,i5)') '      H has dimension', naos, ' x ', naos

      call log_program_substep('diagonalizing')
      call hamiltonians_diag_on_basis(H, S, E_MO, MO_AO)

      call log_program_substep('projecting to the primitive basis')
      call integrals_MO_AO_transform(electrons, MO_AO, MO)

      call log_program_substep('projecting to the primitive basis')
      call integrals_MO_AO_transform(electrons, MO_AO, MO)
      call log_program_step_end

      if (verbose()) then
        write (*, *) '     +------------------------'
        write (*, *) '     MO |  Energy  |  occ'
      end if

      ! todo move that somewhere else
      ! setup occupation
      allocate (occ(naos))
      remaining_electrons = electrons%nelec

      hl_index = 1
      do ii = 1, naos
        if (remaining_electrons > 0) then
          occ(ii) = dble(min(2, remaining_electrons))
          remaining_electrons = remaining_electrons - min(2, remaining_electrons)
          if (remaining_electrons .eq. 0) then
            hl_index = 2
            HOMO = ii
          end if
        else
          select case (hl_index)
          case (3)
            hl_index = 1
          case (2)
            hl_index = 3
            LUMO = ii
          case default
          end select
          occ(ii) = 0d0
        end if

        if (verbose()) then
          write (*, '(a, i4,a,f9.4,a,f3.1, 2a)') '    ', ii, ' | ', E_MO(ii), &
            '|  ', occ(ii), '  ', homolumo(hl_index)
        end if
      end do
      if (logging()) then
        write (*, *) '     +------------------------'
        write (*, '(a,f9.4,a)') '      Δ HOMO-LUMO = ', E_MO(LUMO) - E_MO(HOMO), ' Hartrees'
        write (*, '(a,f9.4,a)') '                    ', (E_MO(LUMO) - E_MO(HOMO))*conv_ev, ' eV'
        write (*, *) ''
      end if

    elseif (ps%route%calc_type .eq. calc_promolecule) then
      call log_program_step('Generating promolecule electron density')
      if (ps%route%charge .ne. 0) then
        call log_err('huckster', 'promolecule density is not compatible with charged species.')
        error stop 1
      end if
      occ = electrons%pop_ao
      MO = electrons%C_AO
      allocate (E_MO(naos))
      E_MO(:) = 0.0
      call log_program_step_end
    end if

    call log_program_step('Saving wavefunction')
    open (iwfn, file=output//'.wfn', action='write', iostat=info)

    if (info .ne. 0) then
      call log_err('huckster', 'Could not open output wfn file')
      error stop - 1
    end if

    ! Put the MOs into the fully unrolled, uncontracted primitive form that goes
    ! into AIM file and related routines.
    call log_program_substep('unrolling MOs to primitive functions')
    call integrals_unroll(electrons, MO, E_MO, occ, umos)
    deallocate (E_MO)
    deallocate (occ)
    deallocate (MO)

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
