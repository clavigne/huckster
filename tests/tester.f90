program tester
  use log
  use params
  use models
  use integrals, only: UnrolledMOs
  implicit none
  character(len=150)           :: arg
  integer :: i
  type(Parameters) :: ps
  type (UnrolledMOs) :: umos

  call set_verbosity(2)

  i = 1
  print *, command_argument_count()
  do while (i .le. command_argument_count())
    call get_command_argument(i, arg)
    print *, "test: "//arg
    select case (arg)

    case ("huckel_h2o")
      call set_verbosity(-1)
      ps = params_default()
      ps%route%calc_type = calc_huckel
      ps%route%find_crits = .false.
      ps%route%path_crits = .false.
      call models_from_params("test-data/h2o.xyz", ps, umos)

      
      
      write(*,*) umos%E

    case default
      call log_err('tester', 'unrecognized command-line option '//trim(arg))
      error stop -1
    end select
    i = i + 1
  end do

  contains 
    subroutine assert(cond, explain)
      logical, intent(in) :: cond
      character(len=*), intent(in) :: explain
    end subroutine
end program
