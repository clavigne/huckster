module crits
  use log
  use ed
  use mkl, only: dsyev, dsysv
  use params, only: CritNR, CritPath, CritSearch, Parameters
  private

  type(CritPath) :: ps_path
  type(CritNR) :: ps_nr
  type(CritSearch) :: ps_search

  double precision, allocatable :: atoms(:, :)
  double precision, allocatable :: radii(:)
  integer :: natm

  ! crit points for output
  type CritPoint
    double precision :: position(3)
    integer :: rank
    integer :: curvature
    integer :: atom_id
    double precision :: electron_density
    double precision :: ellipticity(3)
  end type CritPoint

  ! Internal list of critical points
  type CNode
    ! Data for critical point
    double precision :: R(3)
    integer :: rank
    integer :: curv
    double precision :: rho
    double precision :: w(3)
    double precision :: V(3, 3)
    integer :: niter

    type(CNode), allocatable :: next
    type(CNode), pointer :: prev
  end type CNode

  type(CNode), target :: cps
  type(CNode), pointer :: cps_next => cps
  integer :: ncrits

  ! Statistics
  integer :: crits_stat_searches
  integer :: crits_stat_maxiter
  ! statistics are computed for iteration 1 (first element) and every other iteration (second element)
  ! so that we can tell which grid points were immediately eliminated, for example
  integer :: crits_stat_nodensity(2)
  integer :: crits_stat_noatoms(2)
  integer :: crits_stat_near_atom_cp(2)
  integer :: crits_stat_prev_cp(2)

  character(len=255) :: str

  public :: crits_search, crits_initialize
  public :: crits_close, crits_print
  public :: crits_do_atoms, crits_do_bonds, crits_do_grid
  public :: crits_perceive_graph, crits_get_CPs

  public :: CritPoint

contains
  subroutine crits_initialize(MOs, params)
    use constants
    use integrals
    implicit none
    type(UnrolledMOs), intent(in) :: MOs
    type(Parameters), intent(in) :: params
    ! workspace
    double precision, allocatable :: my_radii(:)
    integer :: i

    ps_nr = params%nr
    ps_path = params%path
    ps_search = params%search

    ! clear any previous data
    call crits_close
    ncrits = 0
    nullify (cps%prev)
    cps_next => cps

    ! Allocate atoms
    natm = MOs%natm
    atoms = MOs%xyz

    ! atomic radii for atomic densities
    allocate (my_radii(natm))
    do i = 1, natm
      my_radii(i) = ATOMIC_RADII(int(MOs%Z(i)))
    end do
    radii = my_radii

    crits_stat_searches = 0
    crits_stat_maxiter = 0
    crits_stat_nodensity = 0
    crits_stat_noatoms = 0
    crits_stat_near_atom_cp = 0
    crits_stat_prev_cp = 0

    ! ready to roll
  end subroutine crits_initialize

  ! Search for CPs midpoint between every atom pair that are closer than MAX_BOND_DISTANCE from eachother
  !
  ! This is not usually necessary because the grid search will usually find those anyway.
  subroutine crits_do_atoms
    implicit none
    integer :: ia
    integer :: info, found

    found = 0
    do ia = 1, natm
      call crits_search(atoms(ia, :), info)
      if (info .eq. 1) then
        found = found + 1
      end if
    end do

    if (logging()) then
      write (*, "('          found ', I3, ' critical points')") found
    end if
  end subroutine crits_do_atoms

  ! Search for CPs midpoint between every atom pair that are closer than MAX_BOND_DISTANCE from eachother
  !
  ! This is not usually necessary because the grid search will usually find those anyway.
  subroutine crits_do_bonds
    implicit none
    integer :: ia, ib
    double precision :: x0(3)
    integer :: info, found

    found = 0
    do ia = 1, natm
      do ib = ia + 1, natm
        if (sqrt(sum((atoms(ia, :) - atoms(ib, :))**2)) < ps_search%max_bond) then
          x0 = (atoms(ia, :) + atoms(ib, :))/2.0
          call crits_search(x0, info)
          if (info .eq. 1) then
            found = found + 1
          end if
        end if
      end do
    end do

    if (logging()) then
      write (*, "('          found ', I3, ' critical points')") found
    end if
  end subroutine crits_do_bonds

  ! Search for CPs everywhere on a grid
  subroutine crits_do_grid
    implicit none
    integer :: i, j, k, n(3)
    integer :: info
    integer :: count
    double precision :: ul(3), ll(3), d(3)
    double precision :: x(3)

    ! Compute grid boundaries
    do i = 1, 3
      ul(i) = maxval(atoms(:, i)) + ps_search%grid*1.5
      ll(i) = minval(atoms(:, i)) - ps_search%grid*1.5
      n(i) = int(ceiling((ul(i) - ll(i))/ps_search%grid))
      d(i) = (ul(i) - ll(i))/n(i)
    end do

    if (logging()) then
      write (*, '(a,i3,a,i3,a,i3)') '  -> grid search, nx=', n(1), '  ny=', n(2), '  nz=', n(3)
      write (*, '(a,i8)') '     tot. number of points = ', product(n)
    end if

    count = 1
    do i = 0, n(1) - 1
      do j = 0, n(2) - 1
        call log_progress_bar(count, product(n))
        do k = 0, n(3) - 1
          x = ll
          x(1) = ll(1) + i*d(1)
          x(2) = ll(2) + j*d(2)
          x(3) = ll(3) + k*d(3)

          call crits_search(x, info)
          if (info .eq. 1 .and. very_verbose()) then
            write (*, '(i3,a,i3,i3,i3)') ncrits, 'g:', i, j, k
          end if
          count = count + 1
        end do
      end do
    end do

    if (logging()) write (*, *) ''
  end subroutine crits_do_grid

  subroutine crits_close
    implicit none
    type(CNode), pointer :: current

    current => cps

    ! go to the end of the list
    do while (.true.)
      if (allocated(current%next)) then
        current => current%next
      else
        exit
      end if
    end do

    ! deallocate moving backwards
    do while (.true.)
      if (associated(current%prev)) then
        current => current%prev
        deallocate (current%next)
      else
        exit
      end if
    end do
  end subroutine crits_close

  subroutine crits_search(x0, info)
    implicit none
    ! Argument
    double precision, intent(in) :: x0(3) ! starting point for search
    integer, intent(out) :: info

    ! Work parameters
    double precision :: rho0, rho1(3), rho2(3, 3)
    double precision :: x(3)
    double precision :: dx_last(3), dx(3)
    integer :: ipiv(3), linfo, i
    integer, parameter :: lwork = 256
    double precision :: work(lwork)

    crits_stat_searches = crits_stat_searches + 1
    x = x0
    info = 0
    dx_last = 0d0

    ! newton raphson
    do i = 1, ps_nr%max_neval
      ! Checks are done in order of copmutational expense :

      ! if we are far enough from every atoms we just don't bother computing anything
      ! this is very cheap
      if (atom_density(x) < ps_nr%min_atom_density) then
        if (i .eq. 1) then
          crits_stat_noatoms(1) = crits_stat_noatoms(1) + 1
        else
          crits_stat_noatoms(2) = crits_stat_noatoms(2) + 1
        end if

        return
      end if

      if (is_colliding(x) .ne. 0) then
        if (i .eq. 1) then
          crits_stat_prev_cp(1) = crits_stat_prev_cp(1) + 1
        else
          crits_stat_prev_cp(2) = crits_stat_prev_cp(2) + 1
        end if

        return
      end if

      call ed_eval(x, rho0, rho1, rho2)
      dx = dx_last*ps_nr%mixing + rho1*(1d0 - ps_nr%mixing)
      dx_last = rho1

      if (is_below_threshold(rho0, rho1)) then
        if (i .eq. 1) then
          crits_stat_nodensity(1) = crits_stat_nodensity(1) + 1
        else
          crits_stat_nodensity(2) = crits_stat_nodensity(2) + 1
        end if

        return
      end if

      ! solve for rho2 B = rho1 to obtain B = rho2^-1 rho1
      call dsysv('L', 3, 1, rho2, 3, ipiv, dx, 3, work, lwork, linfo)

      if (linfo .ne. 0) then
        write (str, "('error in lapack dsysv, info = ', I5)") linfo
        call log_err("crits_search", trim(str))
        info = -1
        return
      end if

      if (sqrt(dot_product(rho1, rho1)) < ps_nr%tol_root) then
        ! found a critical point.
        call crits_push_cp(x, i)
        info = 1
        return
      else

        x = x + dx
      end if
    end do

    crits_stat_maxiter = crits_stat_maxiter + 1
  contains
    double precision function atom_density(x)
      implicit none
      double precision, intent(in) :: x(3)
      double precision :: d
      integer :: i

      atom_density = 0d0
      do i = 1, natm
        d = sum((x - atoms(i, :))**2)
        atom_density = atom_density + exp(-sqrt(d)/radii(i))
      end do

    end function atom_density
  end subroutine crits_search

  logical function is_below_threshold(rho0, rho1)
    double precision, intent(in) :: rho0, rho1(3)

    is_below_threshold = .false.
    if (rho0 < ps_nr%tol_density_abandon) then
      if (abs(rho1(1)) > ps_nr%tol_grad_abandon) return
      if (abs(rho1(2)) > ps_nr%tol_grad_abandon) return
      if (abs(rho1(3)) > ps_nr%tol_grad_abandon) return

      ! rho0 < threshold and so are all three components of rho1
      is_below_threshold = .true.
    end if
  end function is_below_threshold

  integer function is_colliding(x)
    implicit none
    double precision, intent(in) :: x(3)
    integer :: i
    double precision :: d
    type(CNode), pointer :: n
    is_colliding = 0

    ! Check for collision against CPs
    n => cps
    do i = 1, ncrits
      d = sum((x - n%R)**2)
      if (sqrt(d) < ps_nr%rtrust) then
        is_colliding = i
        return
      end if

      n => n%next
    end do
  end function is_colliding

  subroutine crits_push_cp(x, niter)
    implicit none
    double precision, intent(in) :: x(3)
    integer, intent(in) :: niter

    ! Work parameters
    double precision :: rho0, rho1(3), rho2(3, 3)
    double precision :: w(3)
    integer ::  rank, curv, linfo, j
    integer, parameter :: lwork = 256
    double precision :: work(lwork)

    call ed_eval(x, rho0, rho1, rho2)

    ! Diagonalize Hessian
    call dsyev('V', 'U', 3, rho2, 3, w, work, lwork, linfo)

    if (linfo .ne. 0) then
      write (str, "('fatal error in lapack dsyev, info = ', I5)") linfo
      call log_err("crits_push_cp", trim(str))
      error stop - 1
    end if

    rank = 0
    curv = 0
    do j = 1, 3
      if (abs(w(j)) > ps_nr%tol_rank) then
        rank = rank + 1
        curv = curv + int(sign(1d0, w(j)))
      end if
    end do

    ! Push the acquired crit point
    cps_next%R = x
    cps_next%rho = rho0
    cps_next%V = rho2
    cps_next%w = w
    cps_next%niter = niter
    cps_next%rank = rank
    cps_next%curv = curv

    ! Move pointer to the next element
    allocate (cps_next%next)
    cps_next => cps_next%next
    ncrits = ncrits + 1

  end subroutine crits_push_cp

  subroutine crits_get_CPs(OutputCPs, np)
    implicit none
    ! Initialize the adjacency matrix and the list of critical points
    integer, intent(out) :: np
    type(CritPoint), allocatable, intent(out) :: OutputCPs(:)

    ! outputs
    type(CritPoint), allocatable :: C(:)

    ! workspace
    type(CNode), pointer :: ni
    double precision :: dist(natm)
    integer :: k, i

    np = ncrits
    allocate (C(np))

    ! Loop through cps
    ni => cps
    do k = 1, ncrits
      C(k)%position = ni%R
      C(k)%rank = ni%rank
      C(k)%curvature = ni%curv

      if (is_atom(ni)) then
        do i = 1, natm
          dist(i) = sum((ni%R - atoms(i, :))**2)
          C(k)%atom_id = minloc(dist, 1)
        end do
      else
        C(k)%atom_id = -1
      end if

      C(k)%electron_density = ni%rho
      C(k)%ellipticity = ni%w
      ni => ni%next
    end do

    OutputCPs = C
  end subroutine crits_get_CPs

  subroutine crits_perceive_graph(adjacency)
    implicit none
    integer, allocatable, intent(out) :: adjacency(:, :)

    ! adjacency matrix
    integer, allocatable :: A(:, :)

    integer :: i, n, np
    double precision :: x0(3)
    type(CNode), pointer :: ni
    integer :: endpoints(6)

    integer :: counts(3)

    ! Allocate empty adjacency matrix
    counts = 0
    np = ncrits
    allocate (A(np, np))
    A = 0

    ! Loop through cps
    ni => cps
    crits: do i = 1, ncrits
      call log_progress_bar(i, ncrits)
      ! if i is a bond
      if ((ni%curv .eq. -1) .and. (ni%rank .eq. 3)) then
        ! start R, but with a small nudge in the correct direction
        x0 = ni%R + ni%V(:, 3)*ps_path%nudge
        call crits_trace_path(i, x0, endpoints(1), -1d0)

        ! do opposite direction
        x0 = ni%R - ni%V(:, 3)*ps_path%nudge
        call crits_trace_path(i, x0, endpoints(2), -1d0)

        if (very_verbose()) then
          write (*, *) ''
          write (*, '(a,i3,a)', advance='no') 'CP ', i
        end if

        if (endpoints(1) > 0 .and. endpoints(2) > 0) then
          if (very_verbose()) &
            write (*, '(a, i2, a, i2)') ' is a bond between atom ', &
            endpoints(1), ' and atom ', endpoints(2)
          A(i, endpoints(1)) = 1
          A(i, endpoints(2)) = 1
          A(endpoints(1), i) = 1
          A(endpoints(2), i) = 1
          counts(1) = counts(1) + 2
        else
          if (very_verbose()) write (*, '(a, i2, a, i2)') ' is a bond CP but did not converge!!'
          cycle crits
        end if

        ! Test if it connects to rings...
        x0 = ni%R + ni%V(:, 2)*ps_path%nudge
        call crits_trace_path(i, x0, endpoints(3), 1d0)

        x0 = ni%R - ni%V(:, 2)*ps_path%nudge
        call crits_trace_path(i, x0, endpoints(4), 1d0)

        x0 = ni%R + ni%V(:, 1)*ps_path%nudge
        call crits_trace_path(i, x0, endpoints(5), 1d0)

        x0 = ni%R - ni%V(:, 1)*ps_path%nudge
        call crits_trace_path(i, x0, endpoints(6), 1d0)

        do n = 3, 6
          if (endpoints(n) > 0 .and. is_ring(endpoints(n))) then
            if (very_verbose()) write (*, '(a, i2)') '           it also connects to ring CP#', endpoints(n)
            A(i, endpoints(n)) = 1
            A(endpoints(n), i) = 1
            counts(2) = counts(2) + 1
          end if
        end do

      else if ((ni%curv .eq. 1) .and. (ni%rank .eq. 3)) then
        ! start R, but with a small nudge in the direction of the negative
        ! eigenvalue
        x0 = ni%R + ni%V(:, 1)*ps_path%nudge
        call crits_trace_path(i, x0, endpoints(1), 1d0)

        ! do opposite direction
        x0 = ni%R - ni%V(:, 1)*ps_path%nudge
        call crits_trace_path(i, x0, endpoints(2), 1d0)

        if (very_verbose()) write (*, '(a,i3,a)') 'CP ', i, ' is a ring CP'
        do n = 1, 2
          if (endpoints(n) > 0 .and. is_cage(endpoints(n))) then
            if (very_verbose()) write (*, '(a, i2)') '           that connects to cage CP#', endpoints(n)
            A(i, endpoints(n)) = 1
            A(endpoints(n), i) = 1
            counts(3) = counts(3) + 1
          end if
        end do
      end if

      ni => ni%next
    end do crits

    ! output
    adjacency = A

    if (logging()) then
      write (*, *) ''
      write (*, *) ''
      write (*, *) '    number of graph connections'
      write (*, '(a,i4)') '             | atoms to bonds:', counts(1)
      write (*, '(a,i4)') '             | bonds to rings:', counts(2)
      write (*, '(a,i4)') '             | rings to cages:', counts(3)
    end if
  end subroutine crits_perceive_graph

  logical function is_atom(nj)
    type(CNode), pointer :: nj

    if ((nj%curv .eq. -3) .and. (nj%rank .eq. 3)) then
      is_atom = .true.
    else
      is_atom = .false.
    end if
  end function is_atom

  logical function is_ring(i)
    ! true if CP i (including atoms) is a ring
    integer, intent(in) :: i
    type(CNode), pointer :: nj
    integer :: k

    nj => cps
    do k = 2, i
      nj => nj%next
    end do

    if ((nj%curv .eq. 1) .and. (nj%rank .eq. 3)) then
      is_ring = .true.
    else
      is_ring = .false.
    end if
  end function is_ring

  logical function is_cage(i)
    ! true if CP i (including atoms) is a cage
    integer, intent(in) :: i
    type(CNode), pointer :: nj
    integer :: k

    nj => cps
    do k = 2, i
      nj => nj%next
    end do

    if ((nj%curv .eq. 3) .and. (nj%rank .eq. 3)) then
      is_cage = .true.
    else
      is_cage = .false.
    end if
  end function is_cage

  subroutine crits_trace_path(icrit, x, endpoint, constant)
    implicit none

    ! multiplicative constant for the electron density. If it is negative then
    ! we move towards maxima of the electron density, and if it is positive we
    ! move towards minima of the same.
    double precision, intent(in) :: constant

    double precision, intent(in) :: x(3)
    integer, intent(in) :: icrit
    integer, intent(out) :: endpoint
    ! endpoint > 0    => successful, index of final CP
    !          = -1   => threshold electron density
    !          = -10  => ran out of segments
    !          = -20  => error in dopri5

    integer :: i

    ! work
    double precision :: t0, t1
    integer :: itol, iout
    integer, parameter :: lwork = 200, liwork = 200
    double precision :: work(lwork)
    integer :: iwork(liwork), idid
    integer :: ipar(2)
    double precision :: rpar

    ! tolerances
    itol = 0

    ! parameters for dopri5
    iout = 1
    work = 0d0
    iwork = 0
    iwork(3) = -1        ! mute
    ipar(1) = icrit
    ipar(2) = -1
    rpar = constant

    do i = 1, ps_path%max_segments
      ! times...
      t0 = 0.0
      t1 = ps_path%time

      call dopri5(3, &
                  gradrho, t0, x, t1, &
                  ps_path%dopri5_rtol, ps_path%dopri5_atol, itol, &
                  gradrho_out, iout, &
                  work, lwork, iwork, liwork, rpar, ipar, idid)

      if (idid .eq. 2) then
        ! interrupted by gradrho_out: either hit a CP, or went off in the
        ! distance
        endpoint = ipar(2)
        return
        !
      elseif (idid .eq. 1) then
        ! Ran out of time, keep going
        t0 = 0.0
        t1 = 1.0
      else
        write (str, "('error in dopri5, info = ', I5)") idid
        call log_err("crits_trace_path", trim(str))
        endpoint = -20
      end if

    end do
    endpoint = -10

  contains

    subroutine gradrho(n, x, y, f, rpar, ipar)
      implicit none
      integer, intent(in) :: n
      double precision, intent(in) :: x, y(n)
      double precision, intent(out) :: f(n)
      integer, intent(inout) :: ipar(2)
      double precision, intent(in) :: rpar

      ! e density
      double precision :: rho0, rho1(3), rho2(3, 3)

      call ed_eval(y, rho0, rho1, rho2)
      f = rpar*rho1
    end subroutine gradrho

    subroutine gradrho_out(NR, xold, x, y, n, con, icomp, nd, rpar, ipar, irtrn)
      implicit none
      integer, intent(in) :: n, nd, nr
      integer, intent(out) :: irtrn
      double precision, intent(in) :: x, xold, y(n), con(5*nd)
      integer, intent(inout) :: ipar(2)
      integer :: icomp(nd)
      double precision :: rpar

      ! work
      integer :: collision

      if (is_below_threshold(x, y)) then
        ! Not enough electron density, give up
        ipar(2) = -1
        irtrn = -1
        return
      end if

      collision = is_colliding(y)

      if (collision < 1) return
      if (collision .eq. ipar(1)) return

      if (collision > 0) then
        irtrn = -1
        ipar(2) = collision
      end if

    end subroutine gradrho_out
  end subroutine crits_trace_path

  subroutine crits_print
    implicit none
    integer :: i
    type(CNode), pointer :: n
    integer :: nbp, nrp, ncp, nnp
    integer :: hopf

    if (.not. logging()) return

    write (*, '(a,i4,a)') 'Search uncovered ', ncrits, ' critical points...'
    write (*, *) ''

    nnp = 0
    nbp = 0
    nrp = 0
    ncp = 0
    n => cps

    do i = 1, ncrits
      if (verbose()) then
        write (*, '(a, i3)', advance='no') 'CP# ', i
        write (*, '(a,i1,a,i2,a)') '  type: (', n%rank, ',', n%curv, ')'
        write (*, '(a, 3f10.5)') 'coord: ', n%R
        write (*, '(a,  e14.5)') '   rho:  ', n%rho
        write (*, '(a,  3e14.5)') '   ellip:', n%w
        write (*, *) ''
      end if

      if (n%rank .eq. 3) then
        select case (n%curv)
        case (-3)
          nnp = nnp + 1
        case (-1)
          nbp = nbp + 1
        case (1)
          nrp = nrp + 1
        case (3)
          ncp = ncp + 1
        end select
      end if

      n => n%next
    end do
    if (verbose()) then
      write (*, *) 'Numerical analysis ---'
      write (*, '(a,a)') '                   |', '  Number of occurences'
      write (*, '(a,i5)') '           attempts|', crits_stat_searches
      write (*, '(a,i5)') '          successes|', ncrits
      write (*, '(a,i5)') '           failures|', crits_stat_searches - ncrits
      write (*, '(a)') '-------------------| tot.    iter=1'
      write (*, '(a,i5,a, i8)') '    No atom density|', sum(crits_stat_noatoms), '  ', crits_stat_noatoms(1)
      write (*, '(a,i5,a, i8)') 'No electron density|', sum(crits_stat_nodensity), '  ', crits_stat_nodensity(1)
      write (*, '(a,i5,a, i8)') '    Reached an atom|', sum(crits_stat_near_atom_cp), '  ', crits_stat_near_atom_cp(1)
      write (*, '(a,i5,a, i8)') '   Reached other CP|', sum(crits_stat_prev_cp), '  ', crits_stat_prev_cp(1)
      write (*, '(a,i5)') ' Reached max. iter.|', crits_stat_maxiter
    end if

    write (*, *) '    Summary ---'
    write (*, '(a,i3)') '     atom crit. points:', nnp
    write (*, '(a,i3)') '     bond crit. points:', nbp
    write (*, '(a,i3)') '     ring crit. points:', nrp
    write (*, '(a,i3)') '     cage crit. points:', ncp
    write (*, *) ''

    hopf = nnp - nbp + nrp - ncp
    write (*, '(a)') '    Poincare-Hopf number (should = 1) '
    write (*, '(a,a,a,a,a,a,a,a)') '    ', 'nnp', '   ', 'nbp', '   ', 'nrp', '   ', 'ncp'
    write (*, '(a,i3,a,i3,a,i3,a,i3,a,i3)') '    ', nnp, ' - ', nbp, ' + ', nrp, ' - ', ncp, ' = ', hopf
    if (hopf == 1) then
      write (*, '(a)') ':D All expected CPs were found!'
    else
      write (*, '(a)') ':( Some expected CPs were not found!'
    end if

  end subroutine crits_print
end module crits
