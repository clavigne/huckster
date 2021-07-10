program critter
  use integrals
  use constants
  use density
  use log
  use crits

  implicit none

  ! Command line arguments
  character(len=150)           :: arg
  character(len=150):: input_file, output_file
  character(len=:), allocatable :: output_name
  integer :: info

  logical :: do_graph

  ! wfn output
  integer,parameter :: iwfn = 10
  integer,parameter :: icrit = 11

  ! Electronic density
  type(UnrolledMOs) :: MOs

  ! Critical points
  type(CritPoint), allocatable :: CPs(:)
  integer :: ncp
  integer,parameter :: ifile=15

  ! Graph
  integer, allocatable :: adjacency(:,:)

  ! workspace variables
  integer :: i, j

  ! Load command line arguments
  i = 1
  verbosity = 0
  do_graph = .false.
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

     case ('-g', '--graph')
        do_graph = .true.

     case ('-v', '--verbose')
        verbosity = 1

     case ('-vv', '--very-verbose')
        verbosity = 2

     case default
        call log_err('critter', 'unrecognized command-line option ' // trim(arg))
     end select

     i = i + 1
  end do


  ! These next few sections are identical to dancer
  ! ---------------------------------------------------------------------------------
  ! Banner and file loading
  call log_banner_critter

  if ((i .eq. command_argument_count()) .or. (command_argument_count() .eq. 0)) then
     call log_err('critter', 'no wfn file passed')
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
     call log_err('critter', 'could not open input geometry file: ' // input_file)
     error stop -1
  end if

  call integrals_read_wfn(2, MOs)
  close(unit=2)
  call log_program_step_end

  ! --------------------------------------------------------------------------------
  call log_program_step('Finding critical points')

  ! initialize density evaluator
  call density_initialize(MOs)

  ! evaluate critical points
  call crits_initialize(MOs)
  call crits_do_grid
  call crits_print
  ! todo: if some are missing, maybe try do_bonds etc.
  call crits_get_CPs(CPs, ncp)
  if (verbosity .ge. 0) write(*,*) '         writing CPs to csv file'

  open(unit=ifile, file= output_name // '.csv')
  write(ifile, *) 'index,atom,rank,curv,x,y,z,rho,ellip_x,ellip_y,ellip_z'
  do i=1, ncp
  write(ifile,'(i4,a)', advance='no') i, ','
  if (CPs(i)%atom_id > 0) then
      write(ifile,'(a,a)', advance='no') mos%atoms(CPs(i)%atom_id), ','
  else
      write(ifile,'(a)', advance='no') ','
  end if
  write(ifile,'(i2,a)', advance='no') CPs(i)%rank, ','
  write(ifile,'(i2,a)', advance='no') CPs(i)%curvature, ','
  write(ifile,'(E23.15,a)', advance='no') CPs(i)%position(1), ','
  write(ifile,'(E23.15,a)', advance='no') CPs(i)%position(2), ','
  write(ifile,'(E23.15,a)', advance='no') CPs(i)%position(3), ','
  write(ifile,'(E23.15,a)', advance='no') CPs(i)%electron_density, ','
  write(ifile,'(E23.15,a)', advance='no') CPs(i)%ellipticity(1), ','
  write(ifile,'(E23.15,a)', advance='no') CPs(i)%ellipticity(2), ','
  write(ifile,'(E23.15)', advance='no') CPs(i)%ellipticity(3)
  write(ifile,*) ''
  end do
  close(unit=ifile)

  call log_program_step_end

  if (do_graph) then
    call log_program_step('Connecting critical points into graph')
    call crits_perceive_graph(adjacency)
    call log_program_step_end
    if (verbosity .ge. 0) write(*,*) '         writing adjacency matrix to .mat file'

    open(unit=ifile, file=output_name // '.mat')
    do i=1,ncp
        do j=1,ncp
            write(ifile,'(i2)',advance='no') adjacency(j, i)
        end do
        write(ifile,*) ''
    end do

    close(unit=ifile)
    if (verbosity .ge. 0) write(*,*) 
  end if



  if (verbosity .ge. 0) write(*,*) '~~~~~END~~~~END~~~~END~~~~END~~~~END~~~~END~~~~END~~~~~~'

  stop 0

contains
  subroutine print_help
    write(*,'(a)') 'USAGE: critter [options] -- WFN [OUTPUT]'
    write(*,'(a)') ''
    write(*,'(a)') 'Critter is a program that performs QTAIM calculations. Critter is specifically'
    write(*,'(a)') 'designed to work well with huckster semi-empirical electron densities. It takes'
    write(*,'(a)') 'a single mandatory argument (a WFN file) and generates at OUTPUT.csv data for all'
    write(*,'(a)') 'obtained critical points.'
    write(*,'(a)') ''
    write(*,'(a)') 'The -g flag can be used to produce a molecular graph, which is exported as an'
    write(*,'(a)') 'adjacency matrix.'
    write(*,'(a)') ''
    write(*,'(a)') 'Command-line options'
    write(*,'(a)') '   -g, --graph                Link CPs into a molecular graph.'
    write(*,'(a)') ''
    write(*,'(a)') '   -h, --help                 Print these instructions.'
    write(*,'(a)') '   -v, --verbose              Print more information.'
    write(*,'(a)') '   -vv, --very-verbose        Print more information.'
    write(*,'(a)') '   -q, --quiet                Silence output'
  end subroutine print_help


end program critter
