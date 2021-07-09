include 'mkl_vsl.f90'

program dancer
  use integrals
  use constants
  use density
  use log

  ! random numbers
  use mkl_vsl_type
  use mkl_vsl

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

  ! TODO: move
  double precision :: com(3)
  double precision :: rotx(3,3), rotz(3,3), theta(2)
  double precision, parameter :: pi = 3.14159265358979323846264338327
  ! rng
  integer(kind=4) errcode
  integer brng,method,seed,n
  type(vsl_stream_state) :: stream
  ! cube
  double precision :: cubed
  integer :: cuben
  integer :: npts(3)
  double precision :: origin(3)
  double precision :: displ(3,3)
  integer, parameter :: ielectrons = 12
  integer, parameter :: iatoms = 13

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

  ! ---------------------------------------------------------------------------------
  ! Geometry adjustment (TODO: move)
  call log_program_step('Molecular geometry')
  call log_program_substep('centering system')
  com = sum(MOs%xyz,1)/MOs%natm
  do i=1, MOs%natm
    MOs%xyz(i,:) = MOs%xyz(i,:) - com
  end do

  call log_program_substep('performing random rotations')

  brng=VSL_BRNG_MT19937
  method=VSL_RNG_METHOD_UNIFORM_STD 
  seed = int(sum(com) * 1000000)

  errcode=vslnewstream( stream, brng,  seed)
  errcode=vdrnguniform( method, stream, 2, theta, -2 * pi , 2 * pi)
  errcode=vsldeletestream( stream )

  if (verbosity.ge.0) write(*, '(a,i8)') '         seed: ', seed
  if (verbosity.ge.0) write(*, '(a,2f8.2)') '         random angles: ', theta

  rotx(1,:) = (/ 1d0,        0d0,         0d0 /)
  rotx(2,:) = (/ 0d0, cos(theta(1)), -sin(theta(1)) /)
  rotx(3,:) = (/ 0d0, sin(theta(1)), -cos(theta(1)) /)

  rotz(1,:) = (/ cos(theta(2)), -sin(theta(2)), 0d0 /)
  rotz(2,:) = (/ sin(theta(2)), -cos(theta(2)), 0d0 /)
  rotz(3,:) = (/ 0d0,        0d0,         1d0 /)

  ! rotate
  MOs%xyz = matmul(MOs%xyz, rotx)
  MOs%xyz = matmul(MOs%xyz, rotz)
  call log_program_step_end

  ! ---------------------------------------------------------------------------------
  ! Make cube file
  call log_program_step('Generating cube')

  call density_initialize(MOs)

  ! todo
  cuben = 32
  cubed = 0.25d0 / conv_bohr ! to bohr
  origin = -cubed * cuben/2
  displ = 0
  do i=1,3
     displ(i,i) = cubed
  end do
  npts = cuben

  call log_program_substep('opening file: ' // output_name // '.electrons')
  open(unit=ielectrons, file= output_name // '.electrons', access="stream")
  call density_write_cube( &
       ielectrons,&
       origin,&
       displ,&
       npts)
  close(ielectrons)

  call log_program_substep('opening file: ' // output_name // '.atoms')
  open(unit=iatoms, file= output_name // '.atoms', access="stream")
  call density_write_atoms( &
       iatoms,&
       origin,&
       displ,&
       npts)
  close(iatoms)
  call log_program_step_end



  if (verbosity .ge. 0) write(*,*) '~~~~~END~~~~END~~~~END~~~~END~~~~END~~~~END~~~~END~~~~~~'

  stop 0

contains
  subroutine print_help
    write(*,'(a)') 'USAGE: dancer [options] -- WFN [OUTPUT]'
    write(*,'(a)') ''
    write(*,'(a)') 'Dancer is a program to produce electron density cubes from wfn files.'
  end subroutine print_help


end program dancer
