module integrals
  use log
  private

  ! Holder for integral data
  type ElectronicSystem
     integer,allocatable :: atm(:,:)
     integer,allocatable :: bas(:,:)
     integer :: natm, nelec, charge
     integer :: nbas, nprim, norb
     integer :: ncore, nvalence
     double precision, allocatable :: C_AO(:,:), E_AO(:)
     double precision, allocatable :: pop_AO(:)
     integer :: naos
  end type ElectronicSystem

  type UnrolledMOs
     character(len=2), allocatable :: atoms(:)
     integer, allocatable :: Z(:)
     double precision, allocatable :: xyz(:,:)
     double precision, allocatable :: coeffs(:,:)
     double precision, allocatable :: exps(:)
     integer, allocatable :: icnt(:), ityp(:)

     integer :: natm, nprim, nmos
  end type UnrolledMOs

  integer :: current_n_atoms = 0

  ! libcint parameters
  integer,parameter :: ATM_SLOTS  = 6
  integer,parameter :: BAS_SLOTS  = 8
  integer,parameter :: AOS_SLOTS  = 5

  integer,parameter :: CHARGE_OF  = 1
  integer,parameter :: PTR_COORD  = 2
  integer,parameter :: NUC_MOD_OF = 3
  integer,parameter :: PTR_ZETA   = 4

  integer,parameter :: ATOM_OF    = 1
  integer,parameter :: ANG_OF     = 2
  integer,parameter :: NPRIM_OF   = 3
  integer,parameter :: NCTR_OF    = 4
  integer,parameter :: KAPPA_OF   = 5
  integer,parameter :: PTR_EXP    = 6
  integer,parameter :: PTR_COEFF  = 7

  ! AOS parameters
  integer,parameter :: PTR_ATOM = 1
  integer,parameter :: NAOS_OF = 2
  integer,parameter :: PTR_AOS_E = 3
  integer,parameter :: PTR_AOS_C = 4
  integer,parameter :: PTR_AOS_OCC = 5

  ! some libcint functions
  integer,external :: CINTcgto_spheric
  integer,external :: CINTcgto_cart
  double precision,external :: CINTgto_norm

  ! data from the python script
  include 'data-header.f90'

  ! spherical - cartesian conversions
  double precision   :: s2cd(6,5),s2cf(10,7),s2cg(15,9),s2ch(21,11)

  public :: integrals_initialize, integrals_init_from_file, integrals_build_basis
  public :: integrals_build_atomic_orbitals, integrals_overlaps
  public :: integrals_symm_ao2bas, integrals_symm_bas2ao
  public :: integrals_MO_AO_transform
  public :: integrals_write_to_wfn, integrals_unroll
  public :: ElectronicSystem, UnrolledMOs

contains

  subroutine integrals_initialize
    implicit none
    include 'data-content.f90'
    call set_s2c_coef
  end subroutine integrals_initialize

  subroutine integrals_init_from_file(unit, system, charge)
    use constants
    implicit none
    ! IO
    integer,intent(in)  :: unit
    integer,intent(in)  :: charge
    type(ElectronicSystem), intent(out) :: system

    ! Workspace
    character(len=2)  :: atomchar
    integer, allocatable :: atm(:,:)
    integer :: off, i, j
    integer :: info


    read(2,*, iostat=info) system%natm
    if (info .ne. 0) then
       call log_err('integrals_init_from_file', 'could not read atom number from file')
       error stop -1
    end if

    if (system%natm + current_n_atoms > MAX_NATOMS) then
       ! TODO : logger
       write(*,*) 'integrals_init_from_file:'
       print '(a,i8,/)', 'program has been compiled with MAX_NATOMS = ', MAX_NATOMS
       print '(a,i8,a,/)', 'your integral system now contains ', system%natm, ' atoms'
       print '(a/)', 'to run such a large system, recompile with a modified MAX_NATOMS parameter.'
       error stop -1
    end if

    read(2,*)                    !skip comment

    system%nelec = -charge
    system%charge = charge

    ! Set the pointer for the geometry
    off = ENV_COORD_PTR + 3*current_n_atoms

    allocate(atm(ATM_SLOTS, system%natm))
    do i = 1, system%natm
       atm(PTR_COORD, i) = off
       read(2,*, iostat=info) atomchar, ENV(off + 1 : off + 3)
       if (info .ne. 0) then
          call log_err('integrals_init_from_file', 'could not read atom or coordinate')
          error stop -1
       end if

       ENV(off+1) = ENV(off+1) / conv_bohr
       ENV(off+2) = ENV(off+2) / conv_bohr
       ENV(off+3) = ENV(off+3) / conv_bohr

       off = off + 3

       do j = 1, 110
          if (atomchar .eq. atomname(j)) then
             atm(CHARGE_OF, i) = j
             system%nelec = system%nelec + j
             system%ncore = system%ncore + num_core_electrons(j)
             exit
          end if
       end do

       if (j .eq. 110) then
          call log_err('integrals_init_from_file', 'atomic name is invalid =' // atomchar)
          error stop -1
       end if
    end do

    system%atm = atm
    system%nvalence = system%nelec - system%ncore

    ! Increment the overall atom counter
    current_n_atoms = current_n_atoms + system%natm

    ! Write out
    if (verbosity .ge. 0) then
       write(*,'(a,i4)') '          Charge:         ', system%charge
       write(*,'(a,i4)') '          N. of atoms:    ', system%natm
       write(*,'(a,i4)') '          N. of electrons:', system%nelec
       write(*,'(a,i4)') '                   (core):', system%ncore
       write(*,'(a,i4)') '                (valence):', system%nvalence
    end if

    if (verbosity > 0) then
       write(*,'(a,i4)') '          N. atoms so far in evaluator:', current_n_atoms
    end if
  end subroutine integrals_init_from_file

  subroutine integrals_build_basis(system)
    use constants, only:atomname
    implicit none
    type(ElectronicSystem), intent(inout) :: system

    ! workspace
    integer :: atom_Z, iatom, bas_ptr
    integer :: nbas=0, nprim=0, norb=0,  di
    integer :: ishell, iprims, iaos
    integer, allocatable :: bas(:,:)

    ! Count all the basis functions
    if (verbosity .ge. 0) then 
       write(*,*) '                 nshells   norbs    nprims'
    end if
    do iatom=1,system%natm
       atom_Z = system%atm(CHARGE_OF,iatom)
       bas_ptr = paos(PTR_ATOM, atom_Z)
       if (bas_ptr < 1) then
          print '(2a)', 'integrals_build_basis: no basis function for element ',atomname(atom_Z)
          error stop -1
       end if
       ishell = 0                       ! nshells
       iprims = 0                       ! primitives
       iaos = 0                       ! aos
       basis_loop:do while (.true.)
          if (pbas(ATOM_OF, bas_ptr) .eq. atom_Z) then
             ishell = ishell + 1
             di = CINTcgto_spheric(bas_ptr-1, pbas)
             iaos = iaos + di
             iprims = iprims + pbas(NPRIM_OF, bas_ptr) * di

             bas_ptr = bas_ptr + 1
          else
             exit basis_loop
          end if
       end do basis_loop

       nbas = nbas + ishell
       norb = norb + iaos
       nprim = nprim + iprims


       if (verbosity > 0) then
          write(*,'(a,i3,2a,i3,a,i3,a,i3,a,i3)') &
               '      ', iatom,atomname(atom_Z),&
               '       ', ishell,&
               '       ', iaos, &
               '       ', iprims
       end if
    end do


    ! allocate shells
    allocate(bas(BAS_SLOTS, nbas))

    ! build basis and atoms 
    ishell = 1
    do iatom=1, system%natm
       atom_Z = system%atm(CHARGE_OF, iatom)
       bas_ptr = paos(PTR_ATOM, atom_Z)
       basis_loop2:do while (.true.)
          if (pbas(ATOM_OF, bas_ptr) .eq. atom_Z) then
             bas(:,ishell) = pbas(:,bas_ptr)
             bas(ATOM_OF,ishell) = iatom-1 ! zero-based indexing

             if (bas(NCTR_OF,ishell) .ne. 1) then
                call log_err('integrals_build_basis', 'uncontracted shells not supported yet!')
                error stop -1
             end if

             bas_ptr = bas_ptr + 1
             ishell = ishell + 1
          else
             exit basis_loop2
          end if
       end do basis_loop2
    end do

    if (verbosity .ge. 0) then
       write(*,'(a,i4,a,i4,a,i4)') '       Tot.      ', &
            nbas, '      ', norb,  '      ', nprim
    end if

    ! Write results to system
    system%bas = bas
    system%nbas = nbas
    system%nprim = nprim
    system%norb = norb
  end subroutine integrals_build_basis

  subroutine integrals_build_atomic_orbitals(system)
    ! Load atomic orbitals for system.
    use constants, only: atomname
    implicit none
    type(ElectronicSystem), intent(inout) :: system

    ! workspace
    integer :: i, atom_Z, n, ii, jj, kk, k, off, j
    double precision :: nel

    ! ouptuts
    double precision, allocatable :: C_AO(:,:), E_AO(:), pop_AO(:)
    integer :: naos=0

    
    do i=1,system%natm
       naos = naos + paos(NAOS_OF, system%atm(CHARGE_OF, i))
    end do

    if (verbosity > 0) then
       write(*,*) '          atom    ao energies | pop '
    end if
   allocate(C_AO(naos, system%norb))
   allocate(E_AO(naos))
   allocate(pop_AO(naos))

   C_AO = 0.0
   pop_AO = 0.0
   k = 1
   do i=1,system%natm
      atom_Z = system%atm(CHARGE_OF, i)
      n = paos(NAOS_OF,atom_Z)
      ii = paos(PTR_AOS_C, atom_Z) + 1
      jj = paos(PTR_AOS_E, atom_Z) + 1
      kk = paos(PTR_AOS_OCC, atom_Z) + 1

      if (verbosity > 0) write(*,'(a,i3,a)')  '          ',i,atomname(atom_Z)
      ! where the atom group is
      off = k                   ! row offset

      do j=0,n-1
         pop_AO(k) = ENV(kk+j)
         C_AO(off:off+n-1,k) = ENV(ii:ii+n -1)
         E_AO(k) = ENV(jj+j)
         ii = ii + n
         k = k + 1
         if (verbosity >0) write(*, '(a, f8.2, f8.2)')      '                     ', ENV(jj+j),pop_AO(k)
      end do
   end do

   system%C_AO = C_AO
   system%E_AO = E_AO
   system%pop_AO = pop_AO
   system%naos = naos
   if (verbosity >0) write(*, '(a, f8.2)')      '                   total pop: ',sum(pop_AO)
 end subroutine integrals_build_atomic_orbitals

  subroutine integrals_overlaps(system, S)
    implicit none
    type(ElectronicSystem), intent(in) :: system
    double precision, allocatable, intent(out) :: S(:,:)

    ! workspace
    double precision, allocatable :: buf1e(:,:), val
    integer :: ii, jj, k, i,j, di, dj, i2,j2
    integer :: shls(4)

    if(allocated(S)) deallocate(S)
    allocate(S(1:system%norb, 1:system%norb))

    ii = 1
    k = 0                           !todo only L or U part of matrix?
    do i=0, system%nbas-1           ! 0 based index
       shls(1) = i; di = CINTcgto_spheric(i, system%bas)

       jj = 1
       do j=0, system%nbas-1    ! 0 based index
          shls(2) = j; dj = CINTcgto_spheric(j, system%bas)

          call cint1e_ovlp_sph(S(ii:ii+di-1,jj:jj+dj-1), shls, &
               system%atm, system%natm, system%bas, system%nbas, ENV)

!!$          allocate(buf1e(di,dj))
!!$          call cint1e_ovlp_sph(buf1e(1:di,1:dj), shls, &
!!$               system%atm, system%natm, system%bas, system%nbas, ENV)
!!$
!!$          do i2=1,di
!!$             do j2=1,dj
!!$                S(ii+i2-1,jj+j2-1) = buf1e(i2,j2)
!!$             end do
!!$          end do
!!$
!!$          deallocate(buf1e)
          k = k + di*dj
          jj = jj + dj
       end do

       ii = ii + di
    end do

    ! Compute norm of basis
    val = 0.0
    do i=1,system%norb
       val = val + S(i,i)
    end do

    write(*,*) '     computed overlaps:',k
    write(*,'(a, f6.4)') '      basis trace = ', val/dble(system%norb)
    ! todo check if basis trace == norbs
  end subroutine integrals_overlaps

  subroutine integrals_symm_bas2ao(system, O, O_proj)
    ! Project a (real and symmetric) operator in the Gaussian basis to the AO basis.
    implicit none
    type(ElectronicSystem), intent(in) :: system
    double precision, intent(in) :: O(:,:)
    double precision,allocatable, intent(out) :: O_proj(:,:)

    ! workspace
    integer :: naos, norb
    double precision, allocatable :: tmp(:,:) 

    naos = system%naos
    norb = system%norb
    allocate(O_proj(naos,naos))
    allocate(tmp(naos,norb))

    ! Compute: tmp = C_AO(naos,norb) * O(norb,norb)
    call dsymm('R','U',naos, norb, 1.0d0, O, norb, system%C_AO, naos, 0.0d0, tmp, naos)

    ! Compute: O_proj = tmp * C_AO^T
    call dgemm('n', 't', naos, naos, norb, 1.0d0, tmp, naos, system%C_AO, naos, 0.0d0, O_proj, naos)
  end subroutine integrals_symm_bas2ao


  subroutine integrals_symm_ao2bas(system, O, O_proj)
    ! Project a (real and symmetric) operator in the AO basis to the Gaussian basis.
    implicit none
    type(ElectronicSystem), intent(in) :: system
    double precision, intent(in) :: O(:,:)
    double precision,allocatable, intent(out) :: O_proj(:,:)

    ! workspace
    integer :: naos, norb
    double precision, allocatable :: tmp(:,:) 

    naos = system%naos
    norb = system%norb
    allocate(O_proj(naos,naos))
    allocate(tmp(naos,norb))

    ! Compute: tmp(naos, norb) = O(naos,naos) * C_AO(naos,norb) 
    call dsymm('L','U', naos, norb, 1.0d0, O, norb, system%C_AO, naos, 0.0d0, tmp, naos)

    ! Compute: O_proj(norb, norb) = C(naos,norb)^T tmp(norb, naos)
    call dgemm('t', 'n', norb, norb, naos, 1.0d0, system%C_AO, naos, tmp, norb, 0.0d0, O_proj, norb)
  end subroutine integrals_symm_ao2bas

  subroutine integrals_MO_AO_transform(system, MO_AOs, MO_bas)
    ! Project MOs in the AO basis to the original basis
    implicit none
    type(ElectronicSystem), intent(in) :: system
    double precision, intent(in) :: MO_AOs(:,:)
    double precision,allocatable, intent(out) :: MO_bas(:,:)

    ! workspace
    integer :: naos, norb

    naos = system%naos
    norb = system%norb

   ! Compute: MO = C_AO(naos, norb).T * eigenvectors(naos, naos) 
   ! Note that the eigenvectors are returned in the old Hamiltonian matrix.
    if (allocated(MO_bas)) deallocate(MO_bas)
    allocate(MO_bas(norb,naos))
    call dgemm('t', 'n', norb, naos, naos, 1.0d0, system%C_AO, naos, MO_AOs, naos, 0.0d0, MO_bas, norb)
  end subroutine integrals_MO_AO_transform

  subroutine integrals_unroll(system, C_MO, unrolled_MOs)
    use constants, only: atomname
    implicit none
    type(ElectronicSystem), intent(in) :: system
    double precision, intent(in) :: C_MO(1:,1:)
    type(UnrolledMOs), intent(out) :: unrolled_MOs

    ! workspace
    double precision, allocatable ::  coeffs(:), coeffs_cart(:)
    integer :: nbas, natm, naos, i, j, k, ii, jj, imo, n
    integer :: ncart, nsph

    naos = system%naos
    natm = system%natm
    nbas = system%nbas

    ! compute how many cartesian basis we will have
    ncart = 0
    do i=1, nbas
       ncart = ncart + CINTcgto_cart(i-1, system%bas) * system%bas(NPRIM_OF, i)
    end do

    unrolled_MOs%nmos = naos
    unrolled_MOs%natm = natm
    unrolled_MOs%nprim = ncart

    allocate(unrolled_MOs%atoms(natm))
    allocate(unrolled_MOs%Z(natm))
    allocate(unrolled_MOs%xyz(natm, 3))

    do i=1,natm
       unrolled_MOs%atoms(i) = atomname(system%atm(CHARGE_OF,i))
       unrolled_MOs%Z(i) = system%atm(CHARGE_OF,i)
       unrolled_MOs%xyz(i,:) = env(system%atm(PTR_COORD,i)+1 : system%atm(PTR_COORD,i)+3)
    end do

    allocate(unrolled_MOs%icnt(ncart))
    allocate(unrolled_MOs%ityp(ncart))
    allocate(unrolled_MOs%exps(ncart))

    ii = 1                          ! Prim index
    jj = 1                          ! AO index
    do i=1,nbas                                  !shells
       do j=1, system%bas(NPRIM_OF, i)             !primitives
          do k=1,CINTcgto_cart(i-1, system%bas)       !number of cartesian functions
             unrolled_MOs%icnt(ii) = system%bas(ATOM_OF,i) + 1
             unrolled_MOs%exps(ii) = env(system%bas(PTR_EXP,i) + j)
             ! this follows the .wfn / AIM  standard (or should anyway!)
             unrolled_MOs%ityp(ii) = angular_start(system%bas(ANG_OF, i)) + k

             ! update primitive counter
             ii = ii + 1
          end do

          ! update AO index
          jj = jj + 1
       end do
    end do


    allocate(coeffs(1:11))
    allocate(coeffs_cart(1:21))
    allocate(unrolled_MOs%coeffs(ncart, naos))

    do imo=1,naos
       ii = 1
       jj = 1
       do i=1,nbas                                  !shells
          nsph = CINTcgto_spheric(i-1, system%bas)
          do j=1, system%bas(NPRIM_OF, i)             !primitives
             ! ii is the spheric primitive index, jj is the cartesian primitive index

             ! coefficient of primitives * norm * MO coefficients
             coeffs(1:nsph) = cintgto_norm(system%bas(ANG_OF,i), & 
                  env(system%bas(PTR_EXP,i) + j) ) &
                  * env(system%bas(PTR_COEFF,i) + j) * C_MO(ii : ii + nsph - 1, imo)

             ! convert to cartesian representation
             ncart = CINTcgto_cart(i-1, system%bas)
             call sph2car(system%bas(ANG_OF,i), coeffs(1:nsph), coeffs_cart(1:ncart))

             unrolled_MOs%coeffs(jj : jj + ncart - 1, imo) = coeffs_cart(1:ncart)
             jj = jj + ncart
          end do

          ii = ii + nsph
       end do
    end do

  contains 
    pure integer function angular_start(i)
      ! Return the starting index for the basis type in a wfn file of a basis
      ! function with angular momentum i
      integer, intent(in) :: i

      select case(i)
      case(0)
         angular_start=0
      case(1)
         angular_start=1
      case(2)
         angular_start=4
      case(3)
         angular_start=10
      case(4)
         angular_start=20
      case(5)
         angular_start=35
      case default
         angular_start = -1
      end select
    end function angular_start
  end subroutine integrals_unroll


  subroutine integrals_write_to_wfn(iwfn, mos, E_MO, occ)
    use constants, only: atomname
    implicit none
    integer, intent(in) :: iwfn
    double precision, intent(in) :: E_MO(1:)
    double precision, intent(in) :: occ(1:)
    type(UnrolledMOs), intent(in) :: mos

    ! workspace
    integer :: i, j, k

    ! This next bit is adapted from molden2aim
    ! thanks open source !
    write(iwfn,*) ' Extended Huckel Theory orbitals'
    write(iwfn,"('GAUSSIAN',8x,i7,' MOL ORBITALS',i7,' PRIMITIVES',i9,' NUCLEI')") &
         mos%nmos, mos%nprim, mos%natm

    do i=1,mos%natm
       write(iwfn, "(2x,A3,i3,4x,'(CENTRE',i3,') ',3f12.8,'  CHARGE =',f5.1)") &
            mos%atoms(i), i, i, &
            mos%xyz(i,1:3), dble(mos%Z(i))
    end do

    write(iwfn,"('CENTRE ASSIGNMENTS  ',20i3)")(mos%icnt(i), i=1, mos%nprim)
    write(iwfn,"('TYPE ASSIGNMENTS    ',20i3)")(mos%ityp(i), i=1, mos%nprim)
    write(iwfn,"('EXPONENTS ',5d14.7)")(mos%exps(i), i=1, mos%nprim)

    do i=1,mos%nmos
       write(iwfn,"('MO',i5,5x,'MO 0.0',8x,'OCC NO =',f13.7,'  ORB. ENERGY =',f12.6)") &
            i, occ(i), E_MO(i)

       write(iwfn,"(5d16.8)")(mos%coeffs(j,i), j=1,mos%nprim)
    end do

    write(iwfn,"('END DATA')")
    write(iwfn,"(' THE  HF ENERGY =',f20.12,' THE VIRIAL(-V/T)=',f13.8)")0.d0,2.d0
    close(iwfn)
  end subroutine integrals_write_to_wfn



  ! from molden2aim -> convert spherical basis to cartesian
  Subroutine sph2car(lq,fi,fo)
    implicit real(kind=8) (a-h,o-z)
    integer :: lq, i, j
    dimension         :: fi(*),fo(*)

    select case(lq)
    case(0)    ! S
       fo(1)=fi(1)
    case(1)    ! P
       fo(1:3)=fi(1:3)
    case(2)    ! D
       fo(1:6)=0.0d0
       do j=1,5
          do i=1,6
             fo(i)=fo(i)+s2cd(i,j)*fi(j)
          end do
       end do
    case(3)    ! F
       fo(1:10)=0.0d0
       do j=1,7
          do i=1,10
             fo(i)=fo(i)+s2cf(i,j)*fi(j)
          end do
       end do
    case(4)    ! G
       fo(1:15)=0.0d0
       do j=1,9
          do i=1,15
             fo(i)=fo(i)+s2cg(i,j)*fi(j)
          end do
       end do
    case(5)    ! H
       fo(1:21)=0.0d0
       do j=1,11
          do i=1,21
             fo(i)=fo(i)+s2ch(i,j)*fi(j)
          end do
       end do
    case default
       write(*,*) 'sph2car: angular momentum invalid', lq
       error stop -1
    end select

    return
  end Subroutine sph2car

Subroutine set_s2c_coef()
 implicit none

   s2cd = 0.0D0
   s2cd( 1, 1) = -0.5000000000000000D+00
   s2cd( 2, 1) = -0.5000000000000000D+00
   s2cd( 3, 1) =  1.0D0
   s2cd( 5, 2) =  1.0D0
   s2cd( 6, 3) =  1.0D0
   s2cd( 1, 4) =  0.8660254037844387D+00
   s2cd( 2, 4) = -0.8660254037844387D+00
   s2cd( 4, 5) =  1.0D0

   s2cf = 0.0D0
   s2cf( 3, 1) =  1.0D0
   s2cf( 6, 1) = -0.6708203932499368D+00
   s2cf( 9, 1) = -0.6708203932499368D+00
   s2cf( 1, 2) = -0.6123724356957946D+00
   s2cf( 4, 2) = -0.2738612787525830D+00
   s2cf( 7, 2) =  0.1095445115010332D+01
   s2cf( 2, 3) = -0.6123724356957946D+00
   s2cf( 5, 3) = -0.2738612787525830D+00
   s2cf( 8, 3) =  0.1095445115010332D+01
   s2cf( 6, 4) =  0.8660254037844386D+00
   s2cf( 9, 4) = -0.8660254037844386D+00
   s2cf(10, 5) =  1.0D0
   s2cf( 1, 6) =  0.7905694150420949D+00
   s2cf( 4, 6) = -0.1060660171779821D+01
   s2cf( 2, 7) = -0.7905694150420949D+00
   s2cf( 5, 7) =  0.1060660171779821D+01

   s2cg = 0.0D0
   s2cg( 1, 1) =  0.3750000000000000D+00
   s2cg( 2, 1) =  0.3750000000000000D+00
   s2cg( 3, 1) =  1.0D0
   s2cg(10, 1) =  0.2195775164134200D+00
   s2cg(11, 1) = -0.8783100656536799D+00
   s2cg(12, 1) = -0.8783100656536799D+00
   s2cg( 5, 2) = -0.8964214570007953D+00
   s2cg( 8, 2) =  0.1195228609334394D+01
   s2cg(14, 2) = -0.4008918628686366D+00
   s2cg( 7, 3) = -0.8964214570007953D+00
   s2cg( 9, 3) =  0.1195228609334394D+01
   s2cg(13, 3) = -0.4008918628686366D+00
   s2cg( 1, 4) = -0.5590169943749475D+00
   s2cg( 2, 4) =  0.5590169943749475D+00
   s2cg(11, 4) =  0.9819805060619659D+00
   s2cg(12, 4) = -0.9819805060619659D+00
   s2cg( 4, 5) = -0.4225771273642583D+00
   s2cg( 6, 5) = -0.4225771273642583D+00
   s2cg(15, 5) =  0.1133893419027682D+01
   s2cg( 5, 6) =  0.7905694150420950D+00
   s2cg(14, 6) = -0.1060660171779821D+01
   s2cg( 7, 7) = -0.7905694150420950D+00
   s2cg(13, 7) =  0.1060660171779821D+01
   s2cg( 1, 8) =  0.7395099728874520D+00
   s2cg( 2, 8) =  0.7395099728874520D+00
   s2cg(10, 8) = -0.1299038105676658D+01
   s2cg( 4, 9) =  0.1118033988749895D+01
   s2cg( 6, 9) = -0.1118033988749895D+01

   s2ch = 0.0D0
   s2ch( 1, 1) =  1.0D0
   s2ch( 3, 1) = -0.1091089451179962D+01
   s2ch( 5, 1) =  0.625D+00
   s2ch(12, 1) = -0.1091089451179962D+01
   s2ch(14, 1) =  0.3659625273557000D+00
   s2ch(19, 1) =  0.625D+00
   s2ch( 7, 2) =  0.1290994448735806D+01
   s2ch( 9, 2) = -0.5669467095138407D+00
   s2ch(11, 2) =  0.1613743060919757D+00
   s2ch(16, 2) = -0.1267731382092775D+01
   s2ch(18, 2) =  0.2112885636821291D+00
   s2ch(21, 2) =  0.4841229182759271D+00
   s2ch( 2, 3) =  0.1290994448735806D+01
   s2ch( 4, 3) = -0.1267731382092775D+01
   s2ch( 6, 3) =  0.4841229182759271D+00
   s2ch(13, 3) = -0.5669467095138407D+00
   s2ch(15, 3) =  0.2112885636821291D+00
   s2ch(20, 3) =  0.1613743060919757D+00
   s2ch( 3, 4) = -0.1118033988749895D+01
   s2ch( 5, 4) =  0.8539125638299665D+00
   s2ch(12, 4) =  0.1118033988749895D+01
   s2ch(19, 4) = -0.8539125638299665D+00
   s2ch( 8, 5) =  0.1290994448735806D+01
   s2ch(10, 5) = -0.6454972243679028D+00
   s2ch(17, 5) = -0.6454972243679028D+00
   s2ch( 9, 6) = -0.1224744871391589D+01
   s2ch(11, 6) =  0.5229125165837972D+00
   s2ch(16, 6) =  0.9128709291752769D+00
   s2ch(18, 6) =  0.2282177322938192D+00
   s2ch(21, 6) = -0.5229125165837972D+00
   s2ch( 4, 7) = -0.9128709291752769D+00
   s2ch( 6, 7) =  0.5229125165837972D+00
   s2ch(13, 7) =  0.1224744871391589D+01
   s2ch(15, 7) = -0.2282177322938192D+00
   s2ch(20, 7) = -0.5229125165837972D+00
   s2ch( 5, 8) =  0.7395099728874520D+00
   s2ch(14, 8) = -0.1299038105676658D+01
   s2ch(19, 8) =  0.7395099728874520D+00
   s2ch(10, 9) = -0.1118033988749895D+01
   s2ch(17, 9) =  0.1118033988749895D+01
   s2ch(11,10) =  0.1169267933366857D+01
   s2ch(18,10) = -0.1530931089239487D+01
   s2ch(21,10) =  0.7015607600201140D+00
   s2ch( 6,11) =  0.7015607600201140D+00
   s2ch(15,11) = -0.1530931089239487D+01
   s2ch(20,11) =  0.1169267933366857D+01
end Subroutine set_s2c_coef


end module integrals
