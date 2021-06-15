
! ---------------------------------------------------------------------------------- !
! This is a machine-generated data file.
! It was generated using atomic-data.py
! todo: add some identifiers, like time etc.
!
!
! This is not really made for human consumption, but to ensure the final
! program can be statically linked with absolutely no dependencies. It should
! not be modified by hand! Read at your own risk!
!
!                                                            CYRILLE LAVIGNE
! ---------------------------------------------------------------------------------- !

integer, parameter :: SIZEOF_BAS = 393 
integer, parameter :: SIZEOF_AOS = 110 
integer, parameter :: SIZEOF_ENV = 25387 
integer, parameter :: ENV_COORD_PTR = 22387 
integer, parameter :: MAX_NATOMS = 1000 
integer :: pbas(BAS_SLOTS, SIZEOF_BAS)
integer :: paos(AOS_SLOTS, SIZEOF_AOS)
double precision :: ENV(SIZEOF_ENV)
