program tester
  use log
  implicit none

  character(len=150)           :: arg
  integer :: i

  call set_verbosity(2)

  i = 1
  print *, command_argument_count()
  do while (i .le. command_argument_count())
    call get_command_argument(i, arg)
    print *, "test: " // arg
    select case (arg)

    case ("e2e_h2o")
      call e2e("test-data/h2o.xyz")

    case default
      call log_err('tester', 'unrecognized command-line option '//trim(arg))
    end select
    i = i + 1
  end do

  contains
  subroutine e2e(input_file)
    use integrals
    character(len=*), intent(in) :: input_file
    type(ElectronicSystem) :: electrons
    integer :: info

    call integrals_initialize

    open (unit=2, file=input_file, action='READ', iostat=info)
    if (info .ne. 0) then
      call log_err('tester', 'could not open input geometry file: '//input_file)
      error stop - 1
    end if
    call integrals_init_from_file(2, electrons, 0)
    close (unit=2)

    ! call integrals_build_basis(electrons)
    ! call integrals_build_atomic_orbitals(electrons)
    ! call integrals_overlaps(electrons, S_prim)

  end subroutine
end program
