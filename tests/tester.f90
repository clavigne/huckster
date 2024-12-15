program tester
  use log
  implicit none

  character(len=150)           :: arg
  integer :: i

  call set_verbosity(2)

  do while (i .le. command_argument_count())
    call get_command_argument(i, arg)
    select case (arg)
    case ("e2e_benzene")
      
    case default
      call log_err('tester', 'unrecognized command-line option '//trim(arg))
    end select
    i = i + 1
  end do
end program