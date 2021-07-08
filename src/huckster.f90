program huckster
  use integrals
  use constants
  use hamiltonians
  use density
  use crits
  use log
  implicit none

  ! homo lumo
  character(len=4)::homolumo(3)=(/ '    ', 'HOMO', 'LUMO'/)
  integer :: hl_index
  integer :: HOMO
  integer :: LUMO

  ! Command line arguments
  character(len=150)           :: arg
  integer ::  charge
  character(len=150)::input_file, output_file
  character(len=:), allocatable :: output_name
  integer :: info

  integer, parameter :: calc_huckel=1, calc_promolecule=2
  integer :: calc_type
  logical :: do_crit, do_graph

  ! Integrals
  type(ElectronicSystem) :: electrons
  type(UnrolledMOs) :: umos

  ! Work parameters
  character(len=2) :: atomchar
  integer :: i, ii, j, naos, remaining_electrons

  ! Hamiltonian / overlap
  double precision, allocatable :: S_prim(:,:)
  double precision, allocatable :: H(:,:)
  double precision, allocatable :: S(:,:)
  double precision, allocatable :: MO(:,:)
  double precision, allocatable :: MO_AO(:,:)
  double precision, allocatable :: E_MO(:)
  double precision, allocatable :: occ(:)

  ! wfn output
  integer,parameter :: iwfn = 10
  integer,parameter :: icrit = 11

  ! CP output
  type(CritPoint), allocatable :: CPs(:)
  integer :: ncp
  integer, allocatable :: adjacency(:,:)
  integer :: ifile

  ! Load command line arguments
  charge = 0
  i = 1
  calc_type = calc_huckel
  do_crit = .false.
  do_graph = .false.
  verbosity = 0
  do while (i .le. command_argument_count())
     call get_command_argument(i, arg)
     select case (arg)

     case ('--')
        ! prepping to receive final argument
        exit

     case ('-h', '--help')
        call print_help()
        stop 0

     case ('-q', '--quiet')
        verbosity = -1

     case ('-v', '--verbose')
        verbosity = 1

     case ('-vv', '--very-verbose')
        verbosity = 2

     case ('-a', '--aim')
        do_crit = .true.

     case ('-g', '--graph')
        do_crit = .true.        ! we need CPs for the graph
        do_graph = .true.

     case ('-t', '--type')
        call get_command_argument(i+1,arg)

        select case(arg)
        case('pro', 'promolecule', 'sad')
           calc_type = calc_promolecule

        case('huckel','eht')
           calc_type = calc_huckel

        case default
           call log_err('huckster', trim(arg) // ' is not a valid calculation type')
           error stop -1
        end select
        i = i + 1

     case ('-c', '--charge')
        call get_command_argument(i+1,arg)
        read(arg, *, iostat=info) charge
        if (info .ne. 0) then
           call log_err('huckster', trim(arg) // ' is an invalid charge.')
           error stop -1
        else
           i = i + 1
        end if

     case default
        call log_err('huckster', 'unrecognized command-line option ' // trim(arg))
     end select
     i = i + 1
  end do


  ! print out the program banner
  call log_banner

  if ((i .eq. command_argument_count()) .or. (command_argument_count() .eq. 0)) then
     call log_err('huckster', 'no geometry file passed')
     error stop -1
  else
     i = i + 1
     call get_command_argument(i, input_file)

     if (i < command_argument_count()) then
        call get_command_argument(i+1, output_file)
     else
        output_file = input_file
     end if
  end if

  ! remove extension
  associate (idot => index(output_file, '.', back=.true.))
    if (idot > 0) then
       output_name = output_file(1:idot-1)
    else
       output_name = output_file(1:len_trim(output_file))
    end if
  end associate

  ! ---------------------------------------------------------------------------------
  ! Initialize integral module
  call integrals_initialize

  ! ---------------------------------------------------------------------------------
  ! Read input geometry
  call log_program_step('Reading geometry file')
  write(*,*) '     file: ',input_file

  open(unit=2, file=input_file, action='READ', iostat=info)
  if (info .ne. 0) then
     call log_err('huckster', 'could not open input geometry file: ' // input_file)
     error stop -1
  end if
  call integrals_init_from_file(2, electrons, charge)
  close(unit=2)
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

  if (calc_type .eq. calc_huckel) then
     ! ---------------------------------------------------------------------------------
     ! Setup the Hamiltonian
     call log_program_step('Building Extended Huckel Hamiltonian')

     call log_program_substep('generating overlaps in the gaussian basis')
     call integrals_overlaps(electrons, S_prim)

     call log_program_substep('projecting S_orb to opt. AOs')
     call integrals_symm_bas2ao(electrons, S_prim, S)

     call log_program_substep('generating Hamiltonian elements')
     allocate(H(naos,naos))
     do i = 1, naos
        H(i,i) = electrons%E_AO(i)
        do j =i+1, naos
           ! multiply off diag by Ei + Ej / 2
           H(i,j) = S(i,j) * (electrons%E_AO(i) + electrons%E_AO(j)) * 0.5 * K_PARAMETER
        end do
     end do
     call log_program_step_end

     call log_program_step('Eigendecomposition of H')
     if (verbosity .ge. 0) write(*,'(a,i5,a,i5)') '      H has dimension', naos, ' x ', naos

     call log_program_substep('diagonalizing')
     call hamiltonians_diag_on_basis(H, S, E_MO, MO_AO)

     call log_program_substep('projecting to the primitive basis')
     call integrals_MO_AO_transform(electrons, MO_AO, MO)

     call log_program_substep('projecting to the primitive basis')
     call integrals_MO_AO_transform(electrons, MO_AO, MO)
     call log_program_step_end

     if (verbosity > 0) then 
        write(*,*) '     +------------------------'
        write(*,*) '     MO |  Energy  |  occ'
     end if

     ! todo move that somewhere else
     ! setup occupation
     allocate(occ(naos))
     remaining_electrons = electrons%nelec

     hl_index = 1
     do ii=1,naos
        if (remaining_electrons > 0) then
           occ(ii) = dble(min(2, remaining_electrons))
           remaining_electrons = remaining_electrons - min(2, remaining_electrons)
           if (remaining_electrons .eq. 0) then
              hl_index = 2
              HOMO = ii
           end if
        else
           select case (hl_index)
           case(3)
              hl_index = 1
           case (2)
              hl_index = 3
              LUMO = ii
           case default
           end select
           occ(ii) = 0d0
        end if

        if (verbosity > 0) then 
           write(*, '(a, i4,a,f9.4,a,f3.1, 2a)')      '    ',ii, ' | ', E_MO(ii),&
                '|  ', occ(ii),'  ', homolumo(hl_index)
        end if
     end do
     if (verbosity .ge. 0) then 
        write(*,*) '     +------------------------'
        write(*,'(a,f9.4,a)') '      Δ HOMO-LUMO = ', E_MO(LUMO) - E_MO(HOMO), ' Hartrees'
        write(*,'(a,f9.4,a)') '                    ', (E_MO(LUMO) - E_MO(HOMO)) * conv_ev, ' eV'
        write(*,*) ''
     end if


  elseif (calc_type .eq. calc_promolecule) then
     call log_program_step('Generating promolecule electron density')
     if (charge .ne. 0) then
        call log_err('huckster', 'promolecule density is not compatible with charged species.')
        error stop 1
     end if
     occ = electrons%pop_ao
     MO = electrons%C_AO
     allocate(E_MO(naos))
     E_MO(:) = 0.0
     call log_program_step_end
  end if

  call log_program_step('Saving wavefunction')
  open(iwfn, file=output_name // '.wfn', action='write', iostat=info)

  if (info .ne. 0) then
     call log_err('huckster', 'Could not open output wfn file')
     error stop -1
  end if

  ! Put the MOs into the fully unrolled, uncontracted primitive form that goes
  ! into AIM file and related routines.
  call log_program_substep('unrolling MOs to primitive functions')
  call integrals_unroll(electrons, MO, umos)

  call log_program_substep('writing to file: ' // output_name // '.wfn')
  call integrals_write_to_wfn(iwfn, umos, E_MO, occ)
  close(iwfn)

  call log_program_step_end

  if (do_crit) then
     call log_program_step('Finding critical points')
     ! initialize density evaluator
     call density_initialize(umos, occ)
     call crits_initialize(umos)

     call crits_do_grid
     call crits_print
     ! todo: if some are missing, maybe try do_bonds etc.
     call crits_get_CPs(CPs, ncp)
     if (verbosity .ge. 0) write(*,*) '         writing CPs to csv file'

     open(unit=ifile, file= output_name // '.csv')
     write(ifile, *) 'index,atom,rank,curv,x,y,z,rho,ellip_x,ellip_y,ellip_z'
     do i=1, ncp
        write(ifile,'(i4,a)', advance='no') i, ','
        if (CPs(i)%atom_id > 0) then
           write(ifile,'(a,a)', advance='no') umos%atoms(CPs(i)%atom_id), ','
        else
           write(ifile,'(a)', advance='no') ','
        end if
        write(ifile,'(i2,a)', advance='no') CPs(i)%rank, ','
        write(ifile,'(i2,a)', advance='no') CPs(i)%curvature, ','
        write(ifile,'(E23.15,a)', advance='no') CPs(i)%position(1), ','
        write(ifile,'(E23.15,a)', advance='no') CPs(i)%position(2), ','
        write(ifile,'(E23.15,a)', advance='no') CPs(i)%position(3), ','
        write(ifile,'(E23.15,a)', advance='no') CPs(i)%electron_density, ','
        write(ifile,'(E23.15,a)', advance='no') CPs(i)%ellipticity(1), ','
        write(ifile,'(E23.15,a)', advance='no') CPs(i)%ellipticity(2), ','
        write(ifile,'(E23.15)', advance='no') CPs(i)%ellipticity(3)
        write(ifile,*) ''
     end do
     close(unit=ifile)

     call log_program_step_end
  end if

  if (do_graph) then
     call log_program_step('Connecting critical points into graph')
     neval = 0
     call crits_perceive_graph(adjacency)
     call log_program_step_end
     if (verbosity .ge. 0) write(*,*) '         writing adjacency matrix to .mat file'

     open(unit=ifile, file=output_name // '.mat')
     do i=1, ncp
        do j=1,ncp
           write(ifile,'(i2)',advance='no') adjacency(j, i)
        end do
        write(ifile,*) ''
     end do
     close(unit=ifile)
     write(*,*) ''
  end if



  if (verbosity .ge. 0) write(*,*) '~~~~~END~~~~END~~~~END~~~~END~~~~END~~~~END~~~~END~~~~~~'

  stop 0
contains
  subroutine print_help
    write(*,*) 'USAGE: huckster [options] -- GEOMETRY [OUTPUT]'
    write(*,*) ''
    write(*,*) 'Huckster is a standalone program to perform extended Huckel theory calculations,'
    write(*,*) 'specifically to use in QTAIM calculations. It takes a single mandatory argument,'
    write(*,*) 'a molecular geometry in the XMOL format in file GEOMETRY and outputs a AIM-format'
    write(*,*) 'wavefunction file to OUTPUT (defaults to GEOMETRY with a .wfn extension).'

    write(*,*) 'Command-line options'
    write(*,*) '   -c, --charge CHRG          Set the global charge of the computed molecule.'
    write(*,*) '   -t, --type TYPE            Define computation type (TYPE = {eht, pro}.)'
    write(*,*) '   -a, --aim                  Do a QTAIM search for critical points.'
    write(*,*) '   -g, --graph                Build molecular graph from QTAIM analysis.'
    write(*,*) ''

    write(*,*) '   -h, --help                 Print these instructions.'
    write(*,*) '   -v, --verbose              Print more information.'
    write(*,*) '   -vv, --very-verbose        Print more information.'
    write(*,*) '   -q, --quiet                Silence output'

  end subroutine print_help


end program huckster
