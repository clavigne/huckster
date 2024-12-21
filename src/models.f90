module models
  use log
  use integrals
  use params
  use constants
  use hamiltonians
  implicit none
contains
  subroutine models_from_params(input, ps, result)
    character(len=*), intent(in) :: input
    type(Parameters), intent(in) :: ps
    type(UnrolledMOs), intent(out) :: result
    type(ElectronicSystem) :: electrons
    integer :: info

    ! ---------------------------------------------------------------------------------
    ! Initialize integral module
    call integrals_initialize

    ! ---------------------------------------------------------------------------------
    ! Read input geometry
    call log_program_step('Reading geometry file')
    if (logging()) write (*, *) '     file: ', input

    open (unit=2, file=input, action='READ', iostat=info)
    if (info .ne. 0) then
      call log_err('models_from_params', 'could not open input geometry file: '//input)
      error stop - 1
    end if
    call integrals_init_from_file(2, electrons, ps%route%charge)
    close (unit=2)
    call log_program_step_end

    if (ps%route%calc_type .eq. calc_huckel) then
      call log_program_step('Extended Huckel Hamiltonian')
      call models_huckel(ps, electrons, result)
      call log_program_step_end

    else if (ps%route%calc_type .eq. calc_promolecule) then
      call log_program_step('Promolecule model')
      if (ps%route%charge .ne. 0) then
        call log_err('models_from_params', 'promolecule density is not compatible with charged species.')
        error stop 1
      end if

      call models_promolecule(ps, electrons, result)
      call log_program_step_end
    end if

  end subroutine

  subroutine models_setup_huckel_basis(electrons)
    type(ElectronicSystem), intent(inout) :: electrons

    call log_program_substep('setting up Huckel basis')
    call integrals_build_basis(electrons)
    call log_program_substep('Loading AO coeffs and energies')
    call integrals_build_atomic_orbitals(electrons)
  end subroutine

  ! This is just the huckel model, but without any diagonalization
  subroutine models_promolecule(ps, electrons, result)
    type(Parameters), intent(in) :: ps
    type(ElectronicSystem), intent(inout) :: electrons
    type(UnrolledMOs), intent(out) :: result
    double precision, allocatable :: MO(:, :)
    double precision, allocatable :: E_MO(:)
    double precision, allocatable :: occ(:)
    integer :: naos
    call models_setup_huckel_basis(electrons)
    naos = electrons%naos
    occ = electrons%pop_ao
    MO = electrons%C_AO
    allocate (E_MO(naos))
    E_MO = electrons%E_AO

    ! Put the MOs into the fully unrolled, uncontracted primitive form that goes
    ! into AIM file and related routines.
    call log_program_substep('unrolling MOs to primitive functions')
    call integrals_unroll(electrons, MO, E_MO, occ, result)
  end subroutine

  subroutine models_huckel(ps, electrons, result)
    type(Parameters), intent(in) :: ps
    type(ElectronicSystem), intent(inout) :: electrons
    type(UnrolledMOs), intent(out) :: result
    ! Hamiltonian / overlap
    double precision, allocatable :: S_prim(:, :)
    double precision, allocatable :: H(:, :)
    double precision, allocatable :: S(:, :)
    double precision, allocatable :: MO(:, :)
    double precision, allocatable :: MO_AO(:, :)
    double precision, allocatable :: E_MO(:)
    double precision, allocatable :: occ(:)
    ! Work parameters
    integer :: i, ii, j, naos, remaining_electrons

    ! homo lumo
    character(len=4)::homolumo(3) = (/'    ', 'HOMO', 'LUMO'/)
    integer :: hl_index
    integer :: HOMO
    integer :: LUMO
    call models_setup_huckel_basis(electrons)
    naos = electrons%naos

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

    call log_program_substep('diagonalizing the Hamiltonian')
    if (logging()) write (*, '(a,i5,a,i5)') '      H has dimension', naos, ' x ', naos
    call hamiltonians_diag_on_basis(H, S, E_MO, MO_AO)

    call log_program_substep('projecting to the primitive basis')
    call integrals_MO_AO_transform(electrons, MO_AO, MO)

    if (verbose()) then
      write (*, *) '     +------------------------'
      write (*, *) '     MO |  Energy  |  occ'
    end if

    ! todo move that somewhere else
    ! setup occupation
    allocate (occ(naos))
    remaining_electrons = electrons%nelec

    HOMO = 0
    LUMO = 0

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
    if (logging() .and. ((HOMO .ne. 0) .and. (LUMO .ne. 0))) then
      write (*, *) '     +------------------------'
      write (*, '(a,f9.4,a)') '      Δ HOMO-LUMO = ', E_MO(LUMO) - E_MO(HOMO), ' Hartrees'
      write (*, '(a,f9.4,a)') '                    ', (E_MO(LUMO) - E_MO(HOMO))*conv_ev, ' eV'
      write (*, *) ''
    end if

    ! Put the MOs into the fully unrolled, uncontracted primitive form that goes
    ! into AIM file and related routines.
    call log_program_substep('unrolling MOs to primitive functions')
    call integrals_unroll(electrons, MO, E_MO, occ, result)
  end subroutine
end module
