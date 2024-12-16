module log
  implicit none
  include "current-commit.f90"
  integer, private :: verbosity = 0
  integer, private :: step_counter = 0

contains

  subroutine set_verbosity(level)
    integer, intent(in) :: level

    verbosity = level
  end subroutine

  pure logical function logging()
    logging = verbosity .ge. 0
  end function

  pure logical function verbose()
    verbose = verbosity .ge. 1
  end function

  pure logical function very_verbose()
    very_verbose = verbosity .ge. 2
  end function

  subroutine log_out
  end subroutine

  subroutine log_err(routine, string)
    implicit none
    character(len=*), intent(in) :: routine, string
    write (*, '(a,a,a)') routine, ' : ', string
  end subroutine log_err

  subroutine log_banner
    if (.not. logging()) return
    write (*, *) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
    write (*, *) '+======================================================+'
    write (*, *) '|                **   HUCKSTER   **                    |'
    write (*, *) '|                                                      |'
    write (*, *) '|        a program for Extended Huckel Theory          |'
    write (*, *) '|           & promolecule density estimation           |'
    write (*, *) '|        -------------- by ------------------          |'
    write (*, *) '|               -- Cyrille Lavigne --                  |'
    write (*, *) '|                                                      |'
    write (*, "(' |                                    (commit: ', a, ') |')") git_hash
    write (*, *) '+======================================================+'
  end subroutine log_banner

  subroutine log_program_step(step)
    implicit none
    character(len=*) :: step

    if (.not. logging()) return

    write (*, *) '+-------------------------------------------------------'
    write (*, '(a,i3,a)') ' ', step_counter, ' - '//step

  end subroutine log_program_step

  subroutine log_program_substep(step)
    implicit none
    character(len=*) :: step
    if (.not. logging()) return

    write (*, '(a)') '      |'
    write (*, '(3a)') '      ', step, ' ...'
  end subroutine log_program_substep

  subroutine log_program_step_end
    step_counter = step_counter + 1
    if (.not. logging()) return

    write (*, *) '     done!'
    write (*, *) ''
  end subroutine log_program_step_end

  !************************************************************
  ! Maciej Żok, 2010 MIT License
  ! https://github.com/macie/fortran-libs
  subroutine log_progress_bar(iteration, maximum)
    !
    ! Prints progress bar.
    !
    ! Args:
    !     iteration - iteration number
    !     maximum - total iterations
    !
    implicit none
    integer :: iteration, maximum
    integer :: counter
    integer :: step, done

    ! only display at normal and verbose levels, we don't need it at very verbose level because
    ! we already have a pile of output
    if (.not. logging()) return
    if (very_verbose()) return

    step = nint(iteration*100/(1.0*maximum))
    done = floor(step/5.0)  ! mark every 5%

    do counter = 1, 36                    ! clear whole line - 36 chars
      write (6, '(a)', advance='no') '\b'  ! (\b - backslash)
    end do

    write (6, '(a)', advance='no') '     ['
    if (done .LE. 0) then
      do counter = 1, 20
        write (6, '(a)', advance='no') '='
      end do
    else if ((done .GT. 0) .and. (done .LT. 20)) then
      do counter = 1, done
        write (6, '(a)', advance='no') '#'
      end do
      do counter = done + 1, 20
        write (6, '(a)', advance='no') '='
      end do
    else
      do counter = 1, 20
        write (6, '(a)', advance='no') '#'
      end do
    end if
    write (6, '(a)', advance='no') '] '
    write (6, '(I3.1)', advance='no') step
    write (6, '(a)', advance='no') '%'
  end subroutine log_progress_bar
end module log
