module density
  use integrals, only: UnrolledMOs
  private
  ! we will only include basis for which exp(-alpha R^2) > this cutoff
  double precision, parameter :: INTEGRAL_CUTOFF = 1d-6

  ! center positions for each orbitals
  double precision, allocatable :: x(:), y(:), z(:)

  ! x^n etc. powers for each orbital
  integer, allocatable :: ix(:), iy(:), iz(:)

  ! orbital exponents
  double precision, allocatable :: alpha(:)

  ! densmat (norb x nmos)
  double precision, allocatable :: densmat(:, :)

  integer :: norb               ! number of orbitals

  ! Scratch variables used in density evaluations
  integer :: npow
  double precision, allocatable :: alphaR2(:), exp_alphaR2(:)
  double precision, allocatable :: dx(:,:), dy(:,:), dz(:,:)
  double precision, allocatable :: ao_C(:), dao_C(:,:), ddao_C(:,:,:)

  ! temp workspaces used to pack integrals after thresholding. here to avoid
  ! extraneous allocations.
  double precision, allocatable :: packed_densmat(:, :)
  integer, allocatable :: packed_ix(:), packed_iy(:), packed_iz(:)
  double precision, allocatable :: packed_alpha(:)
  integer, allocatable :: indices(:)
  double precision, parameter :: lcutoff = log(INTEGRAL_CUTOFF)

  integer :: neval

  public :: density_initialize, density_eval, density_eval_test
  public :: density_write_cube
  public :: neval

contains
  subroutine density_initialize(orbitals)
    type(UnrolledMOs), intent(in) :: orbitals

    ! Temporary arrays reassigned to module
    double precision, allocatable :: my_x(:), my_y(:), my_z(:)
    double precision, allocatable :: my_alpha(:)
    integer, allocatable :: my_ix(:), my_iy(:), my_iz(:)

    double precision, allocatable :: my_dx(:,:), my_dy(:,:), my_dz(:,:)
    double precision, allocatable :: my_ao_C(:), my_dao_C(:,:), my_ddao_C(:,:,:)
    double precision, allocatable :: my_alphaR2(:), my_exp_alphaR2(:)
    double precision, allocatable :: dmat2(:,:)
    integer, allocatable :: inds(:)

    ! workspace
    integer :: natm, iatm, i, j
    double precision, allocatable :: my_C(:,:)

    norb = orbitals%nprim
    natm = orbitals%natm
    neval = 0

    ! build density matrix
    allocate(my_C(norb,norb))
    my_C = 0d0
    do i=1, orbitals%nmos
       do j=1, norb
          my_C(j,:) = my_C(j,:) + orbitals%coeffs(j,i) * orbitals%coeffs(:,i) &
               * orbitals%occ(i)
       end do
    end do

    ! Fun trick -> take all off diag elements of C to 2x so that we can do only
    ! U part when calculating rho. cuts half of the compute time.
    do i=1,norb
       do j=1, norb
          if (i .ne. j) then
             my_C(i,j) = my_C(i,j) * 2d0
          end if
       end do
    end do

    allocate(my_x(norb), my_y(norb), my_z(norb))
    allocate(my_ix(norb), my_iy(norb), my_iz(norb))
    allocate(my_alpha(norb))
    do i=1, norb
       iatm = orbitals%icnt(i)

       ! coordinates
       my_x(i) = orbitals%xyz(iatm,1)
       my_y(i) = orbitals%xyz(iatm,2)
       my_z(i) = orbitals%xyz(iatm,3)

       ! exponent
       my_alpha(i) = orbitals%exps(i)

       ! powers
       call assign_powers(orbitals%ityp(i), my_ix(i), my_iy(i), my_iz(i))

    end do

    ! Allocate to outer module
    x = my_x
    y = my_y
    z = my_z
    ix = my_ix
    iy = my_iy
    iz = my_iz
    alpha = my_alpha
    densmat = my_C

    ! packing arrays
    allocate(dmat2(norb,norb))
    dmat2 = 0
    packed_densmat = dmat2
    allocate(inds(norb))
    inds = 0
    indices = inds              !move to indices
    packed_ix = inds
    packed_iy = inds
    packed_iz = inds
    packed_alpha = my_alpha


    ! assign scratch variables and allocate scratch arrays
    npow = max(maxval(ix), maxval(iy), maxval(iz))
    npow = max(2, npow)

    allocate(my_alphaR2(norb), my_exp_alphaR2(norb))
    allocate(my_ao_C(norb), my_dao_C(norb,3), my_ddao_C(norb, 3, 3))
    allocate(my_dx(norb, -2:npow), my_dy(norb, -2:npow), my_dz(norb, -2:npow))

    my_dx(:,0) = 1.0
    my_dy(:,0) = 1.0
    my_dz(:,0) = 1.0

    alphaR2 = my_alphaR2
    exp_alphaR2 = my_exp_alphaR2
    ao_C = my_ao_C
    dao_C = my_dao_C
    ddao_C = my_ddao_C
    dx = my_dx
    dy = my_dy
    dz = my_dz

  end subroutine density_initialize

  subroutine assign_powers(itype, l,m,n)
    implicit none
    integer, intent(in) :: itype
    integer, intent(out) :: l,m,n

    ! From molden2aim
    integer  :: PATDAT(3,56)
    data PATDAT/  &
    !              1          2          3          4          5              6          7          8          9         10
            0, 0, 0,   1, 0, 0,   0, 1, 0,   0, 0, 1,   2, 0, 0,       0, 2, 0,   0, 0, 2,   1, 1, 0,   1, 0, 1,   0, 1, 1,  &
    !             11         12         13         14         15             16         17         18         19         20
            3, 0, 0,   0, 3, 0,   0, 0, 3,   1, 2, 0,   2, 1, 0,       2, 0, 1,   1, 0, 2,   0, 1, 2,   0, 2, 1,   1, 1, 1,  &
    !             21         22         23         24         25             26         27         28         29         30
            4, 0, 0,   0, 4, 0,   0, 0, 4,   3, 1, 0,   3, 0, 1,       1, 3, 0,   0, 3, 1,   1, 0, 3,   0, 1, 3,   2, 2, 0,  &
    !             31         32         33         34         35             36         37         38         39         40
            2, 0, 2,   0, 2, 2,   2, 1, 1,   1, 2, 1,   1, 1, 2,       0, 0, 5,   0, 1, 4,   0, 2, 3,   0, 3, 2,   0, 4, 1,  &
    !             41         42         43         44         45             46         47         48         49         50
            0, 5, 0,   1, 0, 4,   1, 1, 3,   1, 2, 2,   1, 3, 1,       1, 4, 0,   2, 0, 3,   2, 1, 2,   2, 2, 1,   2, 3, 0,  &
    !             51         52         53         54         55             56
            3, 0, 2,   3, 1, 1,   3, 2, 0,   4, 0, 1,   4, 1, 0,       5, 0, 0/
    !<<< Gaussian (see also function fnorm_lmn)
    !             21         22         23         24         25             26         27         28         29         30
    !        0, 0, 4,   0, 1, 3,   0, 2, 2,   0, 3, 1,   0, 4, 0,       1, 0, 3,   1, 1, 2,   1, 2, 1,   1, 3, 0,   2, 0, 2,  &
    !             31         32         33         34         35
    !        2, 1, 1,   2, 2, 0,   3, 0, 1,   3, 1, 0,   4, 0, 0/
    !>>>
    save PATDAT

    if(itype < 1 .or. itype > 56)then
        write(*,*) 'assign_powers: itype out of range = ', itype
        error stop -1
    end if

    l=PATDAT(1,itype)
    m=PATDAT(2,itype)
    n=PATDAT(3,itype)

    return
  end subroutine assign_powers

  subroutine density_eval(R, rho, drho, ddrho)
    implicit none
    double precision, intent(in)  ::  R(3)
    double precision, intent(out) :: rho
    double precision, intent(out), optional :: drho(3)
    double precision, intent(out), optional :: ddrho(3,3)

    integer :: i, j, m,n, norb2, iorb, npow2
    double precision :: accum

    neval = neval + 1

    ! So the trick here is to only do those integrals that have a sufficiently
    ! high exponent. We want to do this with minimal branching logic. The way we
    ! set this up is by thresholding at the very beginning so that we do most of
    ! the calculations we need only on the kept orbitals, saved in densed
    ! arrays.

    iorb = 1               ! this is the index for orbitals that are not skipped
    do i=1,norb
       ! first powers
       dx(iorb,1) = (x(i) - R(1))
       dy(iorb,1) = (y(i) - R(2))
       dz(iorb,1) = (z(i) - R(3))

       ! second powers
       dx(iorb,2) = dx(iorb,1) ** 2
       dy(iorb,2) = dy(iorb,1) ** 2
       dz(iorb,2) = dz(iorb,1) ** 2

       alphaR2(iorb) = -alpha(i) * (dx(iorb,2) + dy(iorb,2) + dz(iorb,2))

       if (alphaR2(iorb) < lcutoff) cycle ! we skip this orbital

       ! pack coefficient
       packed_alpha(iorb) = alpha(i)

       ! pack cartesian powers
       packed_ix(iorb) = ix(i)
       packed_iy(iorb) = iy(i)
       packed_iz(iorb) = iz(i)

       ! this array will be used to repack densmat
       indices(iorb) = i
       iorb = iorb + 1
    end do

    ! new number of orbitals
    norb2 = iorb - 1

    ! new max power, in case we removed all integrals with high powers
    npow2 = max(maxval(packed_ix(1:norb2)), maxval(packed_iy(1:norb2)), maxval(packed_iz(1:norb2)))
    npow2 = max(2, npow2)

    ! repack densmat
    do i=1, norb2
       do j=1, norb2
          packed_densmat(i,j) = densmat(indices(i), indices(j))
       end do
    end do

    ! powers of dx, dy, dz
    do i=2,npow2-1
       call vdMul(norb2, dx(:,1), dx(:,i), dx(:,i+1))
       call vdMul(norb2, dy(:,1), dy(:,i), dy(:,i+1))
       call vdMul(norb2, dz(:,1), dz(:,i), dz(:,i+1))
    end do

    if (present(drho) .or. present(ddrho)) then
       ! -1 powers
       call vdInv(norb2, dx(:,1), dx(:,-1))
       call vdInv(norb2, dy(:,1), dy(:,-1))
       call vdInv(norb2, dz(:,1), dz(:,-1))
    end if

    if (present(ddrho)) then
       ! -2 powers
       call vdInv(norb2, dx(:,2), dx(:,-2))
       call vdInv(norb2, dy(:,2), dy(:,-2))
       call vdInv(norb2, dz(:,2), dz(:,-2))
    end if

    ! exp(-alpha R^2)
    call vdExp(norb2, alphaR2, exp_alphaR2)

    ! C(i) = exp(-alpha(i) R(i)^2) x^ix(i) y^iy(i) z^iz(i)
    do i=1, norb2
       ao_C(i) = exp_alphaR2(i) * dx(i, packed_ix(i)) * dy(i, packed_iy(i))  * dz(i, packed_iz(i))
    end do

    if (present(drho) .or. present(ddrho)) then
       ! dC(i)/dx = d/dx exp(-alpha(i) (x^2 + y^2 + z^2)) x^ix(i) y^iy(i) z^iz(i)
       !          =  exp(-alpha(i) R(i)^2) x^ix(i) y^iy(i) z^iz(i) (-2 alpha(i) * x(i))
       !          +  exp(-alpha(i) R(i)^2) ix(i) * x^(ix(i) - 1) y^iy(i) z^iz(i) if ix(i) neq 0
       !          = C(i) * (ix / x - 2 alpha x)
       do i=1,norb2
          dao_C(i,1) = ao_C(i) * (-2d0 * packed_alpha(i) * dx(i,1) + packed_ix(i) * dx(i,-1))
          dao_C(i,2) = ao_C(i) * (-2d0 * packed_alpha(i) * dy(i,1) + packed_iy(i) * dy(i,-1))
          dao_C(i,3) = ao_C(i) * (-2d0 * packed_alpha(i) * dz(i,1) + packed_iz(i) * dz(i,-1))
       end do
    end if


    if ( present(ddrho) ) then
       ! diagonal elements of hessian
       ! d^2 C(i)/dx^2 = d/dx C(i) * (ix/x - 2 alpha x)
       !               = dC(i)/dx * (ix/x - 2 alpha x)
       !               +  C(i) * (-ix/x^2 - 2 alpha)
       do i=1,norb2
          ddao_C(i,1,1) = dao_C(i,1) * (-2d0 * packed_alpha(i) * dx(i,1) + packed_ix(i) * dx(i,-1)) &
               + ao_C(i) * (-2d0 * packed_alpha(i) - packed_ix(i) * dx(i,-2))

          ddao_C(i,2,2) = dao_C(i,2) * (-2d0 * packed_alpha(i) * dy(i,1) + packed_iy(i) * dy(i,-1)) &
               + ao_C(i) * (-2d0 * packed_alpha(i) - packed_iy(i) * dy(i,-2))

          ddao_C(i,3,3) = dao_C(i,3) * (-2d0 * packed_alpha(i) * dz(i,1) + packed_iz(i) * dz(i,-1)) &
               + ao_C(i) * (-2d0 * packed_alpha(i) - packed_iz(i) * dz(i,-2))

          ! off diagonal elements of hessian
          ! d^2 C(i)/dx dy = d/dy C(i) * (ix/x - 2 alpha x)
          !               = dC(i)/dy * (ix/x - 2 alpha x)
          ddao_C(i,1,2) = dao_C(i,2) * (-2d0 * packed_alpha(i) * dx(i,1) + packed_ix(i) * dx(i,-1))
          ddao_C(i,1,3) = dao_C(i,3) * (-2d0 * packed_alpha(i) * dx(i,1) + packed_ix(i) * dx(i,-1))
          ddao_C(i,2,3) = dao_C(i,3) * (-2d0 * packed_alpha(i) * dy(i,1) + packed_iy(i) * dy(i,-1))
       end do
    end if

    ! Calculate density
    rho = 0d0
    do j=1, norb2
       do i=1,j
          rho = rho + packed_densmat(i,j) * ao_C(i) * ao_C(j)
       end do
    end do

    if ( present(drho) ) then
       ! Calculate first order derivatives
       do m=1,3
          accum = 0d0
          do j=1, norb2
             do i=1,j
                accum = accum + packed_densmat(i,j) * (&
                     + dao_C(i,m) * ao_C(j) &
                     + ao_C(i) * dao_C(j, m))
             end do
          end do
          drho(m) = accum
       end do
    end if

    if ( present(ddrho) ) then
       ! Calculate second order derivatives
       do m=1,3
          do n=m,3
             accum = 0d0
             do j=1, norb2
                do i=1,j
                   accum = accum + packed_densmat(i,j) * (&
                        + ddao_C(i,m, n) * ao_C(j) &
                        + dao_C(i,m) * dao_C(j, n) &
                        + dao_C(i,n) * dao_C(j, m) &
                        + ao_C(i) * ddao_C(j, m, n))
                end do
             end do

             ddrho(m,n) = accum
             ddrho(n,m) = accum
          end do
       end do
    end if

  end subroutine density_eval

subroutine density_eval_test(R)
  implicit none
  double precision,intent(in) :: R(3)
  double precision :: Rdx(3), Rdy(3), Rdz(3)
  double precision :: r0, r1(3), r2(3,3)
  double precision :: r0x, r1x(3), r2x(3,3)
  integer :: i

  Rdx = R
  Rdy = R
  Rdz = R

  Rdx(1) = Rdx(1) + 1d-5
  Rdy(2) = Rdy(2) + 1d-5
  Rdz(3) = Rdz(3) + 1d-5

  ! density
  write(*,*) "rho0"
  call density_eval(R, r0, r1, r2)
  write(*,*) r0
  write(*,*) ""

  ! 1st derivatives
  write(*,*) "d rho0 / dx, numerical vs analytical"
  call density_eval(Rdx, r0x, r1x, r2x)
  write(*,*) "x", (r0 - r0x)/1d-5, r1(1)
  call density_eval(Rdy, r0x, r1x, r2x)
  write(*,*) "y", (r0 - r0x)/1d-5, r1(2)
  call density_eval(Rdz, r0x, r1x, r2x)
  write(*,*) "z", (r0 - r0x)/1d-5, r1(3)
  write(*,*) ""

  ! 2nd derivatives
  write(*,*) "d^2 rho0 / dx dy, numerical vs analytical, diagonal"
  call density_eval(Rdx, r0x, r1x, r2x)
  write(*,*) "xx", (r1(1) - r1x(1))/1d-5, r2(1,1)

  call density_eval(Rdy, r0x, r1x, r2x)
  write(*,*) "yy", (r1(2) - r1x(2))/1d-5, r2(2,2)

  call density_eval(Rdz, r0x, r1x, r2x)
  write(*,*) "zz", (r1(3) - r1x(3))/1d-5, r2(3,3)
  write(*,*) ""

  write(*,*) "d^2 rho0 / dx dy, numerical vs analytical, off-diagonal"
  call density_eval(Rdx, r0x, r1x, r2x)
  write(*,*) "xy", (r1(2) - r1x(2))/1d-5, r2(1,2)

  call density_eval(Rdx, r0x, r1x, r2x)
  write(*,*) "xz", (r1(3) - r1x(3))/1d-5, r2(1,3)

  call density_eval(Rdy, r0x, r1x, r2x)
  write(*,*) "yz", (r1(3) - r1x(3))/1d-5, r2(2,3)
  write(*,*) ""

end subroutine density_eval_test

subroutine density_write_cube(iel, origin, displacements, npts)
   use log, only: verbosity, log_program_substep, log_progress_bar
   implicit none
   integer, intent(in) :: iel
   double precision, intent(in) :: origin(3), displacements(3,3)
   integer, intent(in) :: npts(3)

   integer :: i,j,k,count
   double precision :: x(3), rho

    if (verbosity .ge. 0) then 
       write(*,'(a,i3,a,i3,a,i3)') '    -> cube  nx=',npts(1),'  ny=',npts(2),'  nz=',npts(3)
       write(*,'(a,i10)') '       tot. number of points = ', product(npts)
    end if
   count = 1


   call log_program_substep('electronic density')
   do i=0,npts(1)-1
      call log_progress_bar(i+1, npts(1))
      do j=0,npts(2)-1
         do k=0,npts(3)-1
            ! center of the cube element
            x = origin &
                    + i * displacements(1,:) &
                    + j * displacements(2,:) &
                    + k * displacements(3,:)

            ! do electron
            call density_eval(x, rho)
            write(iel) rho
            count = count + 1
         end do
      end do 
   end do
   if (verbosity .ge. 0) write(*,*)

end subroutine density_write_cube

end module density
