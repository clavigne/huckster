module hamiltonians
  
contains
  subroutine hamiltonians_diag_on_basis(H, S, w, C)
    implicit none
    double precision, intent(in) :: H(1:,1:), S(1:,1:)
    double precision, allocatable, intent(out) :: w(:), C(:,:)

    ! workspace
    integer :: naos, info, lwork 
    double precision, allocatable :: work(:), Scopy(:,:)
    naos = size(H, 1)
    allocate(w(naos))
    Scopy = S
    C = H

   ! workspace query
   lwork = -1
   info = 0
   allocate(work(1))
   call dsygv(1, 'V', 'U',  naos,&
        C, naos, Scopy, naos, w, &
        work, lwork, info)

   if (info .ne. 0) then
      write(*,*) 'hamiltonians_diag_on_basis: error in dsygv work query', info
      error stop -1
   end if

   lwork = int(work(1))

   ! allocate work array
   deallocate(work)
   allocate(work(lwork))

   ! Diagonalizing...
   call dsygv(1, 'V', 'U',  naos,&
        C, naos, Scopy, naos, w, &
        work, lwork, info)

   if (info .ne. 0) then
      write(*,*) 'hamiltonians_diag_on_basis: error in dsygv diagonalization', info
      error stop -1
   end if
end subroutine hamiltonians_diag_on_basis
end module hamiltonians
