module history

   implicit none

   integer, private :: max_events = 0, n_events = 0

   type event 
      integer, dimension(:), allocatable :: part1,  part2   ! части кластера
      integer :: n1, n2		!размеры частей n1 >= n2
      integer :: time		!время столкновения или время разваливания
      integer :: fusion		!if grow == 1 - fusion, if grow == -1 - dissociation
   end type event

   type evptr
       type(event), pointer :: p =>null()
   end type evptr  
   
   type(evptr), dimension(:), allocatable :: all_events
 
   private :: grow_all_events

contains
!--------------------------------------------
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
subroutine add_fusion(new_p,new_cl,old_cl,time)

   implicit none
   integer :: new_cl(*),old_cl(*), time
   integer :: old1_cl(1000), i,j,k
   type(event), pointer :: new_p
    
!   write(*,*) 'hist: add fusion procedure', time, new_cl(1:5),old_cl(1:5)
   i = 1
   j = 1
   k = 1
   do while (new_cl(i) > 0)
     if (new_cl(i) == old_cl(j)) then 
       j = j + 1
     else
       old1_cl(k) = new_cl(i)
       k = k + 1
     endif
     i = i + 1
   enddo
    write(*,*) 'hist: fusion: size1 = ', j-1, old_cl(1:j-1),' size2 = ', k-1, old1_cl(1:k-1)

   if (n_events+1 .ge. max_events) then 
      call grow_all_events(100)
   endif
   n_events = n_events + 1
   allocate(all_events(n_events)%p)
   new_p => all_events(n_events)%p
   
   if (j .gt. k) then 
      allocate(new_p%part1(j-1))
      allocate(new_p%part2(k-1))
      new_p%part1(:) = old_cl(1:j-1)
      new_p%part2(:) = old1_cl(1:k-1)
      new_p%n1 = j-1
      new_p%n2 = k-1
   else
      allocate(new_p%part1(k-1))
      allocate(new_p%part2(j-1))
      new_p%part2(:) = old_cl(1:j-1)
      new_p%part1(:) = old1_cl(1:k-1)
      new_p%n2 = j-1
      new_p%n1 = k-1
   endif
   new_p%time = time
   new_p%fusion = 1

end subroutine add_fusion
!----------------------------------------------
subroutine add_diss(new_p,new_cl,old_cl,time)

   implicit none
   integer :: new_cl(*),old_cl(*), time
   integer :: new1_cl(1000), i,j,k
   type(event), pointer :: new_p
!   write(*,*) 'hist: add dissociation procedure', time, new_cl(1:5),old_cl(1:5)

   i = 1
   j = 1
   k = 1
   do while (old_cl(i) > 0)
     if (old_cl(i) == new_cl(j)) then 
       j = j + 1
     else
       new1_cl(k) = old_cl(i)
       k = k + 1
     endif
     i = i + 1
   enddo
   write(*,*) 'hist: dissociation:  size1 = ', j-1, new_cl(1:j-1),' size2 = ', k-1, new1_cl(1:k-1)
   
   if (n_events+1 .ge. max_events) then 
      call grow_all_events(100)
   endif
   n_events = n_events + 1
   allocate(all_events(n_events)%p)
   new_p => all_events(n_events)%p

   if (j .gt. k) then 
      allocate(new_p%part1(j-1))
      allocate(new_p%part2(k-1))
      new_p%part1(:) = new_cl(1:j-1)
      new_p%part2(:) = new1_cl(1:k-1)
      new_p%n1 = j-1
      new_p%n2 = k-1
   else
      allocate(new_p%part1(k-1))
      allocate(new_p%part2(j-1))
      new_p%part2(:) = new_cl(1:j-1)
      new_p%part1(:) = new1_cl(1:k-1)
      new_p%n2 = j-1
      new_p%n1 = k-1
   endif
   new_p%time = time
   new_p%fusion = -1 
   
end subroutine add_diss
!----------------------------------------------
subroutine rm_event(ev)

   implicit none
   type(event), pointer :: ev
   integer :: i
   
   if (.not. associated(ev)) then
     write(*,*) 'hist: warning try delete bad event' 
     return
   endif
   write(*,*) 'hist: remove event',n_events
   
   i = 1
   do while (i .le. n_events)
      if ( associated(ev, all_events(i)%p) ) then
!         write(*,*) 'hist: ', ev%part1,ev%part2
         write(*,*) 'hist: del event nimber', i
         n_events = n_events -1
         do while (i .le. n_events)
            all_events(i)%p => all_events(i+1)%p
            i = i+1
         enddo
         deallocate(all_events(n_events+1)%p%part1)
         deallocate(all_events(n_events+1)%p%part2)
         deallocate(all_events(n_events+1)%p)
         nullify(all_events(n_events+1)%p)
         
         deallocate(ev%part1)
         deallocate(ev%part2)
         deallocate(ev)
         nullify(ev)
         return
      endif 
      i = i+1
   enddo
   write(*,*) 'hist: bad pointer. No event'

end subroutine rm_event
!----------------------------------------------

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
integer function loop0(ev1,ev2)

   implicit none
   type(event), pointer :: ev1,ev2
   integer :: i
!   write(*,*) 'hist: Check for zero loop'

   loop0 = 0
   if (.not. associated(ev1)) return
   if (.not. associated(ev2)) return
   if (ev1%part1(1) .eq. ev2%part1(1)) then
        if (size(ev1%part1(:)) .ne. size(ev2%part1(:))) return
        do i=1, size(ev1%part1(:))
          if (ev1%part1(i) .ne. ev2%part1(i)) return
        enddo 
        
        if (size(ev1%part2(:)) .ne. size(ev2%part2(:))) return
        do i=1, size(ev1%part2(:))
          if (ev1%part2(i) .ne. ev2%part2(i)) return
        enddo 
    else if (ev1%part1(1) .eq. ev2%part2(1)) then
        if (size(ev1%part1(:)) .ne. size(ev2%part2(:))) return
        do i=1, size(ev1%part1(:))
          if (ev1%part1(i) .ne. ev2%part2(i)) return
        enddo 
        
        if (size(ev1%part2(:)) .ne. size(ev2%part1(:))) return
        do i=1, size(ev1%part2(:))
          if (ev1%part2(i) .ne. ev2%part1(i)) return
        enddo 
    
    else 
        return
    endif
    loop0 = 1
end function loop0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
integer function loop1(ev1,ev2,ev3)

   implicit none
   type(event), pointer :: ev1,ev2,ev3
   integer :: i

   loop1 = 0
   if (.not. associated(ev1)) return
   if (.not. associated(ev2)) return
   if (.not. associated(ev3)) return


end function loop1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

subroutine grow_all_events(n)

   implicit none
   integer n,i
   type(evptr),  dimension(:), allocatable :: buf
  
   write(*,*) 'hist: start grow all_events array', max_events, n_events
   if (allocated(all_events)) then 
      allocate(buf(size(all_events(:))))
      do i = 1, size(all_events(:))
        buf(i)%p => all_events(i)%p
      enddo
      deallocate(all_events)
   endif
   
   max_events = max_events + n
   allocate(all_events(max_events))
     
   if (allocated(buf)) then 
      do i = 1, size(buf(:))
        all_events(i)%p => buf(i)%p 
      enddo
      deallocate(buf)
   endif
   write(*,*) 'hist: finalize grow all_events array', max_events, n_events

end subroutine grow_all_events
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11111





end module history