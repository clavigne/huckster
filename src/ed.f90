module ed
  use log
  use integrals, only: UnrolledMOs
  use mkl, only: ddot, vdAdd, vdPackV, vdExp, vdMul, vdAdd, dsymm

  private

  type(UnrolledMOs) :: orbitals

  public :: ed_initialize, ed_close, ed_eval, ed_eval_test, ed_eval_bench

  integer, public :: neval

  integer :: ncenter, norb
  double precision, allocatable :: centers(:, :), densmat(:, :)
  double precision, allocatable :: malpha(:) ! -alpha
  integer, allocatable :: icnt(:)

  ! code is simpler if gx,gy,gz etc are contiguous
  integer, parameter :: g = 1, gx = 2, gy = 3, gz = 4
  integer, parameter :: gxx = 5, gyy = 6, gzz = 7
  integer, parameter :: gxy = 8, gxz = 9, gyz = 10

  ! polys in final form
  integer :: nterms, npow
  double precision, allocatable :: terms_coeff(:)
  integer, allocatable :: terms_powers(:, :), terms_for(:), terms_alpha(:), terms_centers(:)

  ! workspace, here to avoid allocating in the hot loop
  double precision, allocatable :: exps1(:), exps2(:), exps(:), dr(:, :, :), rho_orb(:, :), flat(:)
  double precision, allocatable :: dmat_rho(:, :)

contains
  subroutine ed_initialize(orbs)
    implicit none
    type(UnrolledMOs), intent(in) :: orbs

    ! workspace
    integer :: i, j

    orbitals = orbs

    ncenter = orbitals%natm
    norb = orbitals%nprim

    allocate (centers(ncenter, 3))
    allocate (malpha(norb))
    allocate (icnt(norb))
    centers(:, :) = orbitals%xyz
    malpha(:) = -1d0*orbitals%exps

    ! NOTE: 0-indexed for vdpacki
    icnt(:) = orbitals%icnt - 1

    ! build density matrix
    allocate (densmat(norb, norb))
    densmat = 0d0
    do i = 1, orbitals%nmos
      do j = 1, norb
        densmat(j, :) = densmat(j, :) + orbitals%coeffs(j, i)*orbitals%coeffs(:, i) &
                        *orbitals%occ(i)
      end do
    end do

    call ed_compute_terms(norb)

    ! create our temporaries
    allocate (dr(ncenter, 3, 0:npow))
    dr(:, :, 0) = 1d0 ! powers of zero are always 1
    allocate (exps(norb), exps1(norb), exps2(norb))
    allocate (rho_orb(norb, g:gyz), flat(g:gyz))
    allocate (dmat_rho(norb, g:gyz))
  end subroutine

  subroutine ed_compute_terms(n)
    implicit none
    integer, intent(in) :: n
    integer :: i, k, prev, which, t0, t1
    integer:: count0, count1, count2, count
    ! How does this thing work?
    !
    ! We have to compute up to 10 different terms of rho (rho, rho_x, rho_y, rho_z, rho_xx etc.), each built from orbital coefficients multiplied
    ! by exp(-alpha R**2) and a polynomial in x, y and z.
    !
    ! For each orbital we setup a bunch of Terms, these are temporary and they only exist to make our life easier right here. This is basically a
    ! handwritten, artisanal computer algebra system.
    type Term
      logical :: used = .false.
      double precision :: c
      integer :: center
      integer :: pow(3)
    end type Term

    ! The number of terms for each contribution is always the same, as will be seen below
    !
    ! note that the assign method relies on the first dimension being the orbital
    type(Term) :: rho0(n, g:g, 1)
    type(Term) :: rho1(n, gx:gz, 2)
    type(Term) :: rho2d(n, gxx:gzz, 6)
    type(Term) :: rho2o(n, gxy:gyz, 4)
    if (logging()) then
      call log_program_substep("Generating polynomials...")
    end if

    ! rho calculation is straightforward!
    do i = 1, n
      ! The mo coefficient will be in the norb x norb density matrix and the exp(-alpha R^2) contribution is calculated separately, so we start
      ! with coefficient 1 for the base terms:
      rho0(i, g, 1)%used = .true.
      rho0(i, g, 1)%c = 1d0
      rho0(i, g, 1)%center = orbitals%icnt(i)
      call assign_powers(orbitals%ityp(i), rho0(i, g, 1)%pow)
    end do

    do i = 1, n
      ! first derivatives (C(i) below is exp(-alpha R^2)) for orbital i
      !
      ! dC(i)/dx = d/dx exp(-alpha(i) (x^2 + y^2 + z^2)) x^ix(i) y^iy(i) z^iz(i)
      !          =  exp(-alpha(i) R(i)^2) x^ix(i) y^iy(i) z^iz(i) (-2 alpha(i) * x(i))
      !          +  exp(-alpha(i) R(i)^2) ix(i) * x^(ix(i) - 1) y^iy(i) z^iz(i) if ix(i) neq 0
      !          = C(i) * (ix / x - 2 alpha x)
      !
      ! Note that each of the contribution is two terms
      do k = 1, 3
        which = gx + k - 1
        call modify(rho0(i, g, 1), rho1(i, which, 1), +1d0*rho0(i, g, 1)%pow(k), k, -1)
        call modify(rho0(i, g, 1), rho1(i, which, 2), 2d0*malpha(i), k, 1)
      end do
    end do

    ! Second derivatives (diagonal elements)
    !
    ! d^2 C(i)/dx^2 = d/dx C(i) * (ix/x - 2 alpha x)
    !               = dC(i)/dx * (ix/x - 2 alpha x)
    !               +  C(i) * (-ix/x^2 - 2 alpha)
    !
    ! now the fun really starts!
    do i = 1, n
      do k = 1, 3
        prev = gx + k - 1
        which = gxx + k - 1
        t1 = 1
        do t0 = 1, 2
          ! this is dC(i)/dx * (ix/x - 2 alpha x)
          call modify(rho1(i, prev, t0), rho2d(i, which, t1), +1d0*rho0(i, g, 1)%pow(k), k, -1)
          call modify(rho1(i, prev, t0), rho2d(i, which, t1 + 1), 2d0*malpha(i), k, 1)
          t1 = t1 + 2
        end do
        ! and this is C(i) * (-ix/x^2 - 2 alpha)
        call modify(rho0(i, g, 1), rho2d(i, which, t1), -1d0*rho0(i, g, 1)%pow(k), k, -2)
        call modify(rho0(i, g, 1), rho2d(i, which, t1 + 1), 2d0*malpha(i), k, 0)
      end do
    end do

    ! Second derivatives (off-diagonal elements)
    !
    ! d^2 C(i)/dx dy = d/dy C(i) * (ix/x - 2 alpha x)
    !               = dC(i)/dy * (ix/x - 2 alpha x)
    do i = 1, n
      t1 = 1
      do t0 = 1, 2
        call modify(rho1(i, gy, t0), rho2o(i, gxy, t1), +1d0*rho0(i, g, 1)%pow(1), 1, -1)
        call modify(rho1(i, gy, t0), rho2o(i, gxy, t1 + 1), 2d0*malpha(i), 1, 1)

        call modify(rho1(i, gz, t0), rho2o(i, gxz, t1), +1d0*rho0(i, g, 1)%pow(1), 1, -1)
        call modify(rho1(i, gz, t0), rho2o(i, gxz, t1 + 1), 2d0*malpha(i), 1, 1)

        call modify(rho1(i, gz, t0), rho2o(i, gyz, t1), +1d0*rho0(i, g, 1)%pow(2), 2, -1)
        call modify(rho1(i, gz, t0), rho2o(i, gyz, t1 + 1), 2d0*malpha(i), 2, 1)
        t1 = t1 + 2
      end do
    end do

    ! Now we put the polynomials in a convenient form for the calculator
    if (logging()) then
      write (*, "('        Terms')")
    end if
    count0 = count_used_terms(rho0)
    if (logging()) then
      write (*, "('             density:         ', I6)") count0
    end if

    count1 = count_used_terms(rho1)
    if (logging()) then
      write (*, "('           + 1st derivatives: ', I6)") count1
    end if

    count2 = count_used_terms(rho2d)
    count2 = count_used_terms(rho2o) + count2
    count = count0 + count1 + count2
    if (logging()) then
      write (*, "('           + 2nd derivatives: ', I6)") count2
      write (*, "('                            ----------')")
      write (*, "('           =                ', I6, ' polynomial terms')") count
    end if

    nterms = count
    allocate (terms_coeff(nterms))
    allocate (terms_centers(nterms))
    allocate (terms_powers(nterms, 1:3))
    allocate (terms_alpha(nterms))
    allocate (terms_for(nterms))

    k = 1
    call assign(k, rho0(:, g, :), g)

    do i = gx, gz
      call assign(k, rho1(:, i, :), i)
    end do

    do i = gxx, gzz
      call assign(k, rho2d(:, i, :), i)
    end do

    do i = gxy, gyz
      call assign(k, rho2o(:, i, :), i)
    end do

    npow = max(maxval(terms_powers), 2)
  contains

    ! Copy a term, multiplying its coefficient by mult and modifying its powers by power
    !
    ! If the output would have a power less than zero, it is dropped, as this is used to implement derivatives.
    subroutine modify(from, to, mult, dim, power)
      type(Term), intent(in) :: from
      type(Term), intent(out) :: to
      double precision, intent(in) :: mult
      integer, intent(in) :: dim, power

      if (to%used) then
        call log_err("ed_compute_terms", "attempted to modify already created term")
        error stop - 1
      end if

      if (.not. from%used) then
        to%used = .false.
        return
      end if

      to%pow = from%pow
      to%pow(dim) = to%pow(dim) + power

      if (to%pow(dim) .lt. 0) then
        to%used = .false.
        return
      end if

      to%c = from%c*mult
      to%center = from%center
      to%used = .true.

    end subroutine

    pure integer function count_used_terms(terms)
      type(Term), intent(in) :: terms(:, :, :)
      integer :: lmn(3)
      integer :: i, j, k

      lmn = shape(terms)

      count_used_terms = 0
      do i = 1, lmn(1)
        do j = 1, lmn(2)
          do k = 1, lmn(3)
            if (terms(i, j, k)%used) count_used_terms = count_used_terms + 1
          end do
        end do
      end do

    end function

    subroutine assign(count, terms, is_for)
      integer, intent(in) :: is_for
      integer, intent(inout) :: count
      type(Term), intent(in) :: terms(:, :)
      integer :: lmn(2)
      integer :: i, k

      lmn = shape(terms)

      do i = 1, lmn(1)
        do k = 1, lmn(2)
          if (.not. terms(i, k)%used) cycle

          terms_for(count) = is_for
          terms_coeff(count) = terms(i, k)%c
          terms_centers(count) = terms(i, k)%center
          terms_powers(count, 1:3) = terms(i, k)%pow
          terms_alpha(count) = i
          count = count + 1
        end do
      end do

    end subroutine

  end subroutine

  subroutine ed_close
    implicit none

    ! deallocate( &
    !   densmat, &
    !   ix, iy, iz, &
    !   dr, dr2, exp_factor, &
    ! )
  end subroutine

  subroutine ed_eval(R, rho, drho, ddrho)
    implicit none
    double precision, intent(in)  ::  R(3)
    double precision, intent(out) :: rho
    double precision, intent(out), optional :: drho(3)
    double precision, intent(out), optional :: ddrho(3, 3)

    double precision :: term
    integer :: i

    neval = neval + 1

    do i = 1, ncenter
      dr(i, :, 1) = centers(i, :) - R
    end do

    ! ! computer powers of 2 up
    do i = 1, npow - 1
      call vdMul(ncenter*3, dr(:, :, 1), dr(:, :, i), dr(:, :, i + 1))
    end do

    ! vdpackv does the equivalent of dr(icnt, 1, 2) but fAsTeR
    !
    ! this whole bit of code here just sums the squares for all the orbitals
    call vdPackV(norb, dr(:, 1, 2), icnt, exps)
    call vdPackV(norb, dr(:, 2, 2), icnt, exps1)
    call vdAdd(norb, exps1, exps, exps2)
    call vdPackV(norb, dr(:, 3, 2), icnt, exps1)
    call vdAdd(norb, exps1, exps2, exps)

    ! and then computes e^-alpha R2
    call vdMul(norb, exps, malpha, exps1)
    call vdExp(norb, exps1, exps) ! now exp contains (-alpha R2)

    rho_orb = 0d0
    do i = 1, nterms
      term = dr(terms_centers(i), 1, terms_powers(i, 1)) &
             *dr(terms_centers(i), 2, terms_powers(i, 2)) &
             *dr(terms_centers(i), 3, terms_powers(i, 3))*exps(terms_alpha(i))*terms_coeff(i)

      rho_orb(terms_alpha(i), terms_for(i)) = rho_orb(terms_alpha(i), terms_for(i)) + term
    end do

    if (present(ddrho)) then
      ! this takes most of the cycles, but can't really be optimized!
      call dsymm('L', 'U', norb, gyz, 1d0, densmat, norb, rho_orb, norb, 0d0, dmat_rho, norb)
    else if (present(drho)) then
      call dsymm('L', 'U', norb, gz, 1d0, densmat, norb, rho_orb(:, g:gz), norb, 0d0, dmat_rho(:, g:gz), norb)
    else
      call dsymm('L', 'U', norb, g, 1d0, densmat, norb, rho_orb(:, g:g), norb, 0d0, dmat_rho(:, g:g), norb)
    end if

    rho = dorb(rho_orb(:, g), dmat_rho(:, g))

    if (present(drho)) then
      drho(1) = dorb(rho_orb(:, gx), dmat_rho(:, g)) + dorb(rho_orb(:, g), dmat_rho(:, gx))
      drho(2) = dorb(rho_orb(:, gy), dmat_rho(:, g)) + dorb(rho_orb(:, g), dmat_rho(:, gy))
      drho(3) = dorb(rho_orb(:, gz), dmat_rho(:, g)) + dorb(rho_orb(:, g), dmat_rho(:, gz))
    end if

    if (present(ddrho)) then

      ddrho(1, 1) = dorb(rho_orb(:, gxx), dmat_rho(:, g)) &
                    + dorb(rho_orb(:, gx), dmat_rho(:, gx))*2d0 &
                    + dorb(rho_orb(:, g), dmat_rho(:, gxx))

      ddrho(2, 2) = dorb(rho_orb(:, gyy), dmat_rho(:, g)) &
                    + dorb(rho_orb(:, gy), dmat_rho(:, gy))*2d0 &
                    + dorb(rho_orb(:, g), dmat_rho(:, gyy))

      ddrho(3, 3) = dorb(rho_orb(:, gzz), dmat_rho(:, g)) &
                    + dorb(rho_orb(:, gz), dmat_rho(:, gz))*2d0 &
                    + dorb(rho_orb(:, g), dmat_rho(:, gzz))

      ddrho(1, 2) = dorb(rho_orb(:, gxy), dmat_rho(:, g)) &
                    + dorb(rho_orb(:, gx), dmat_rho(:, gy)) &
                    + dorb(rho_orb(:, gy), dmat_rho(:, gx)) &
                    + dorb(rho_orb(:, g), dmat_rho(:, gxy))
      ddrho(2, 1) = ddrho(1, 2)

      ddrho(1, 3) = dorb(rho_orb(:, gxz), dmat_rho(:, g)) &
                    + dorb(rho_orb(:, gx), dmat_rho(:, gz)) &
                    + dorb(rho_orb(:, gz), dmat_rho(:, gx)) &
                    + dorb(rho_orb(:, g), dmat_rho(:, gxz))
      ddrho(3, 1) = ddrho(1, 3)

      ddrho(2, 3) = dorb(rho_orb(:, gyz), dmat_rho(:, g)) &
                    + dorb(rho_orb(:, gy), dmat_rho(:, gz)) &
                    + dorb(rho_orb(:, gz), dmat_rho(:, gy)) &
                    + dorb(rho_orb(:, g), dmat_rho(:, gyz))
      ddrho(3, 2) = ddrho(2, 3)
    end if

  contains
    double precision function dorb(v1, v2)
      double precision, intent(in) :: v1(norb)
      double precision, intent(in) :: v2(norb)
      dorb = ddot(norb, v1, 1, v2, 1)
    end function

  end subroutine

  subroutine assign_powers(itype, lmn)
    implicit none
    integer, intent(in) :: itype
    integer, intent(out) :: lmn(3)

    ! From molden2aim
    integer  :: PATDAT(3, 56)
    data PATDAT/ &
      !              1          2          3          4          5              6          7          8          9         10
      0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 2, 0, 0, 0, 2, 0, 0, 0, 2, 1, 1, 0, 1, 0, 1, 0, 1, 1, &
      !             11         12         13         14         15             16         17         18         19         20
      3, 0, 0, 0, 3, 0, 0, 0, 3, 1, 2, 0, 2, 1, 0, 2, 0, 1, 1, 0, 2, 0, 1, 2, 0, 2, 1, 1, 1, 1, &
      !             21         22         23         24         25             26         27         28         29         30
      4, 0, 0, 0, 4, 0, 0, 0, 4, 3, 1, 0, 3, 0, 1, 1, 3, 0, 0, 3, 1, 1, 0, 3, 0, 1, 3, 2, 2, 0, &
      !             31         32         33         34         35             36         37         38         39         40
      2, 0, 2, 0, 2, 2, 2, 1, 1, 1, 2, 1, 1, 1, 2, 0, 0, 5, 0, 1, 4, 0, 2, 3, 0, 3, 2, 0, 4, 1, &
      !             41         42         43         44         45             46         47         48         49         50
      0, 5, 0, 1, 0, 4, 1, 1, 3, 1, 2, 2, 1, 3, 1, 1, 4, 0, 2, 0, 3, 2, 1, 2, 2, 2, 1, 2, 3, 0, &
      !             51         52         53         54         55             56
      3, 0, 2, 3, 1, 1, 3, 2, 0, 4, 0, 1, 4, 1, 0, 5, 0, 0/
    !<<< Gaussian (see also function fnorm_lmn)
    !             21         22         23         24         25             26         27         28         29         30
    !        0, 0, 4,   0, 1, 3,   0, 2, 2,   0, 3, 1,   0, 4, 0,       1, 0, 3,   1, 1, 2,   1, 2, 1,   1, 3, 0,   2, 0, 2,  &
    !             31         32         33         34         35
    !        2, 1, 1,   2, 2, 0,   3, 0, 1,   3, 1, 0,   4, 0, 0/
    !>>>
    save PATDAT

    if (itype < 1 .or. itype > 56) then
      write (*, *) 'assign_powers: itype out of range = ', itype
      error stop - 1
    end if

    lmn(1) = PATDAT(1, itype)
    lmn(2) = PATDAT(2, itype)
    lmn(3) = PATDAT(3, itype)

    return
  end subroutine assign_powers

  subroutine ed_eval_test(R)
    implicit none
    double precision, intent(in) :: R(3)
    double precision :: Rdx(3), Rdy(3), Rdz(3)
    double precision :: r0, r1(3), r2(3, 3)
    double precision :: r0x, r1x(3), r2x(3, 3)

    Rdx = R
    Rdy = R
    Rdz = R

    Rdx(1) = Rdx(1) + 1d-5
    Rdy(2) = Rdy(2) + 1d-5
    Rdz(3) = Rdz(3) + 1d-5

    ! density
    write (*, *) "rho0"
    call ed_eval(R, r0, r1, r2)
    write (*, *) r0
    write (*, *) ""

    ! 1st derivatives
    write (*, *) "d rho0 / dx, numerical vs analytical"
    call ed_eval(Rdx, r0x, r1x, r2x)
    write (*, *) "x", (r0 - r0x)/1d-5, r1(1)
    call ed_eval(Rdy, r0x, r1x, r2x)
    write (*, *) "y", (r0 - r0x)/1d-5, r1(2)
    call ed_eval(Rdz, r0x, r1x, r2x)
    write (*, *) "z", (r0 - r0x)/1d-5, r1(3)
    write (*, *) ""

    ! 2nd derivatives
    write (*, *) "d^2 rho0 / dx dy, numerical vs analytical, diagonal"
    call ed_eval(Rdx, r0x, r1x, r2x)
    write (*, *) "xx", (r1(1) - r1x(1))/1d-5, r2(1, 1)

    call ed_eval(Rdy, r0x, r1x, r2x)
    write (*, *) "yy", (r1(2) - r1x(2))/1d-5, r2(2, 2)

    call ed_eval(Rdz, r0x, r1x, r2x)
    write (*, *) "zz", (r1(3) - r1x(3))/1d-5, r2(3, 3)
    write (*, *) ""

    write (*, *) "d^2 rho0 / dx dy, numerical vs analytical, off-diagonal"
    call ed_eval(Rdx, r0x, r1x, r2x)
    write (*, *) "xy", (r1(2) - r1x(2))/1d-5, r2(1, 2)

    call ed_eval(Rdx, r0x, r1x, r2x)
    write (*, *) "xz", (r1(3) - r1x(3))/1d-5, r2(1, 3)

    call ed_eval(Rdy, r0x, r1x, r2x)
    write (*, *) "yz", (r1(3) - r1x(3))/1d-5, r2(2, 3)
    write (*, *) ""

  end subroutine ed_eval_test

  subroutine ed_eval_bench()
    implicit none
    integer, parameter :: n = 100
    double precision :: r0, r1(3), r2(3, 3)
    double precision :: r(3, n)

    real :: t1, t2
    integer :: i

    call random_number(r)

    call cpu_time(t1)

    do i = 1, n
      call ed_eval(r(:, i), r0, r1, r2)
    end do

    call cpu_time(t2)

    write (*, "('total time: ', F9.2, ' ms | per iter: ', F6.2, ' ms')") (t2 - t1)*1000d0, (t2 - t1)*1000d0/n

  end subroutine
end module ed
