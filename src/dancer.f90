program dancer
  use integrals
  use log
  implicit none

  ! Command line arguments
  character(len=150)           :: arg
  integer ::  charge
  character(len=150)::input_file, output_file
  character(len=:), allocatable :: output_name
  integer :: info

  ! wfn output
  integer,parameter :: iwfn = 10
  integer,parameter :: icrit = 11

  ! Electronic density
  type(UnrolledMOs) :: MOs

  ! workspace variables
  integer :: i

  ! Load command line arguments
  i = 1
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

     case default
        call log_err('dancer', 'unrecognized command-line option ' // trim(arg))
     end select

     i = i + 1
  end do


  ! print out the program banner (TODO)
!!$  call log_banner_dancer

  if ((i .eq. command_argument_count()) .or. (command_argument_count() .eq. 0)) then
     call log_err('dancer', 'no geometry file passed')
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
  ! Read input wfn
  call log_program_step('Input')
  call log_program_substep('reading from file: ' // input_file)

  open(unit=2, file=input_file, action='READ', iostat=info)
  if (info .ne. 0) then
     call log_err('dancer', 'could not open input geometry file: ' // input_file)
     error stop -1
  end if

  call integrals_read_wfn(2, MOs)
  close(unit=2)
  call log_program_step_end

!!$
!!$  ! ---------------------------------------------------------------------------------
!!$  ! Generate the basis
!!$  call log_program_step('Loading primitive basis')
!!$  call integrals_build_basis(electrons)
!!$  call log_program_step_end
!!$
!!$  call log_program_step('Loading AO coeffs and energies')
!!$  call integrals_build_atomic_orbitals(electrons)
!!$
!!$  naos = electrons%naos
!!$  call log_program_step_end
!!$
!!$  if (calc_type .eq. calc_huckel) then
!!$     ! ---------------------------------------------------------------------------------
!!$     ! Setup the Hamiltonian
!!$     call log_program_step('Building Extended Huckel Hamiltonian')
!!$
!!$     call log_program_substep('generating overlaps in the gaussian basis')
!!$     call integrals_overlaps(electrons, S_prim)
!!$
!!$     call log_program_substep('projecting S_orb to opt. AOs')
!!$     call integrals_symm_bas2ao(electrons, S_prim, S)
!!$
!!$     call log_program_substep('generating Hamiltonian elements')
!!$     allocate(H(naos,naos))
!!$     do i = 1, naos
!!$        H(i,i) = electrons%E_AO(i)
!!$        do j =i+1, naos
!!$           ! multiply off diag by Ei + Ej / 2
!!$           H(i,j) = S(i,j) * (electrons%E_AO(i) + electrons%E_AO(j)) * 0.5 * K_PARAMETER
!!$        end do
!!$     end do
!!$     call log_program_step_end
!!$
!!$     call log_program_step('Eigendecomposition of H')
!!$     if (verbosity .ge. 0) write(*,'(a,i5,a,i5)') '      H has dimension', naos, ' x ', naos
!!$
!!$     call log_program_substep('diagonalizing')
!!$     call hamiltonians_diag_on_basis(H, S, E_MO, MO_AO)
!!$
!!$     call log_program_substep('projecting to the primitive basis')
!!$     call integrals_MO_AO_transform(electrons, MO_AO, MO)
!!$
!!$     call log_program_step_end
!!$
!!$     if (verbosity > 0) then 
!!$        write(*,*) '     +------------------------'
!!$        write(*,*) '     MO |  Energy  |  occ'
!!$     end if
!!$
!!$     ! todo move that somewhere else
!!$     ! setup occupation
!!$     allocate(occ(naos))
!!$     remaining_electrons = electrons%nelec
!!$
!!$     hl_index = 1
!!$     do ii=1,naos
!!$        if (remaining_electrons > 0) then
!!$           if ((remaining_electrons.le.electrons%nvalence).or.do_core_electrons) then
!!$              occ(ii) = dble(min(2, remaining_electrons))
!!$           end if
!!$           remaining_electrons = remaining_electrons - min(2, remaining_electrons)
!!$           if (remaining_electrons .eq. 0) then
!!$              hl_index = 2
!!$              HOMO = ii
!!$           end if
!!$        else
!!$           select case (hl_index)
!!$           case(3)
!!$              hl_index = 1
!!$           case (2)
!!$              hl_index = 3
!!$              LUMO = ii
!!$           case default
!!$           end select
!!$           occ(ii) = 0d0
!!$        end if
!!$
!!$        if (verbosity > 0) then 
!!$           write(*, '(a, i4,a,f9.4,a,f3.1, 2a)')      '    ',ii, ' | ', E_MO(ii),&
!!$                '|  ', occ(ii),'  ', homolumo(hl_index)
!!$        end if
!!$     end do
!!$     if (verbosity .ge. 0) then 
!!$        write(*,*) '     +------------------------'
!!$        write(*,'(a,f9.4,a)') '      Δ HOMO-LUMO = ', E_MO(LUMO) - E_MO(HOMO), ' Hartrees'
!!$        write(*,'(a,f9.4,a)') '                    ', (E_MO(LUMO) - E_MO(HOMO)) * conv_ev, ' eV'
!!$        write(*,*) ''
!!$     end if
!!$
!!$
!!$  elseif (calc_type .eq. calc_promolecule) then
!!$     call log_program_step('Generating promolecule electron density')
!!$     if (charge .ne. 0) then
!!$        call log_err('huckster', 'promolecule density is not compatible with charged species.')
!!$        error stop 1
!!$     end if
!!$     if (.not.do_core_electrons) then
!!$        call log_err('huckster', 'promolecule density is not compatible valence-only calc.')
!!$        error stop 1
!!$     end if
!!$     occ = electrons%pop_ao
!!$     MO = electrons%C_AO
!!$     allocate(E_MO(naos))
!!$     E_MO(:) = 0.0
!!$     call log_program_step_end
!!$  end if
!!$
!!$  call log_program_step('Saving wavefunction')
!!$  open(iwfn, file=output_name // '.wfn', action='write', iostat=info)
!!$
!!$  if (info .ne. 0) then
!!$     call log_err('huckster', 'Could not open output wfn file')
!!$     error stop -1
!!$  end if
!!$
!!$  ! Put the MOs into the fully unrolled, uncontracted primitive form that goes
!!$  ! into AIM file and related routines.
!!$  call log_program_substep('unrolling MOs to primitive functions')
!!$  call integrals_unroll(electrons, MO, umos)
!!$
!!$  call log_program_substep('writing to file: ' // output_name // '.wfn')
!!$  call integrals_write_to_wfn(iwfn, umos, E_MO, occ)
!!$  close(iwfn)
!!$
!!$  if (verbosity .ge. 0) write(*,*) '~~~~~END~~~~END~~~~END~~~~END~~~~END~~~~END~~~~END~~~~~~'
!!$
!!$  stop 0
contains
  subroutine print_help
    write(*,'(a)') 'USAGE: dancer [options] -- WFN [OUTPUT]'
    write(*,'(a)') ''
    write(*,'(a)') 'Dancer is a program to produce electron density cubes from wfn files.'
  end subroutine print_help


end program dancer
