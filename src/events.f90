module events

   implicit none

   integer, private :: max_events = 0, n_events = 0

   type event 
!      integer, dimension(:), allocatable :: part1,  part2   ! части кластера
      integer :: n1, n2		!размеры частей n1 >= n2
      integer :: time = 0		!время столкновения или время разваливания
      integer :: fusion = 0		!if grow == 1 - fusion, if grow == -1 - dissociation
      integer :: t_next = 0
      integer :: written = 0
      integer, dimension(:), pointer :: atoms => null()   ! все атомы      
   end type event

   type evptr			! ссылка на событие
       type(event), pointer :: p =>null()
   end type evptr  
   
   type(evptr), dimension(:), allocatable :: all_events
 
   private :: grow_all_events

contains

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
subroutine add_fusion(new_p,new_cl,old_cl,time)

   implicit none
   integer, intent(in) :: new_cl(*),old_cl(*), time
   type(event), pointer :: new_p
   

   integer :: old1_cl(1000), i,j,k
    
!   write(*,*) 'event: add fusion procedure', time, new_cl(1:5),old_cl(1:5)
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

   if (n_events+1 .ge. max_events) then 
      call grow_all_events(100)
   endif
   n_events = n_events + 1
   allocate(all_events(n_events)%p)
   allocate(all_events(n_events)%p%atoms(j-1+k-1))

   all_events(n_events)%p%n1 = j-1
   all_events(n_events)%p%n2 = k-1
   all_events(n_events)%p%atoms(1:j-1) = old_cl(1:j-1)
   all_events(n_events)%p%atoms(j:j-1+k-1) = old1_cl(1:k-1)
   all_events(n_events)%p%time = time
   all_events(n_events)%p%fusion = 1

   new_p => all_events(n_events)%p

   write(*,*) 'event:',n_events,' fusion: size1 = ', j-1, old_cl(1:j-1),' size2 = ', k-1, old1_cl(1:k-1)

end subroutine add_fusion
!----------------------------------------------
subroutine add_diss(new_p,new_cl,old_cl,time)

   implicit none
   integer,intent(in) :: new_cl(*),old_cl(*), time
   type(event), pointer :: new_p

   integer :: new1_cl(1000), i,j,k
!   write(*,*) 'event: add dissociation procedure', time, new_cl(1:5),old_cl(1:5)

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

   if (n_events+1 .ge. max_events) then 
      call grow_all_events(100)
   endif
   
   n_events = n_events + 1
   allocate(all_events(n_events)%p)
   allocate(all_events(n_events)%p%atoms(j-1+k-1))
 
   all_events(n_events)%p%n1 = j-1
   all_events(n_events)%p%n2 = k-1
   all_events(n_events)%p%atoms(1:j-1) = new_cl(1:j-1)
   all_events(n_events)%p%atoms(j:j-1+k-1) = new1_cl(1:k-1)
   all_events(n_events)%p%time = time
   all_events(n_events)%p%fusion = 1

   new_p => all_events(n_events)%p

   write(*,*) 'event:',n_events,' dissociation:  size1 = ', j-1, new_cl(1:j-1),' size2 = ', k-1, new1_cl(1:k-1)

end subroutine add_diss
!----------------------------------------------
subroutine rm_event(ev_target)

   implicit none
   type(event), pointer :: ev,ev_target
   integer :: i

   ev => ev_target
   if (.not. associated(ev)) then
     write(*,*) 'event: warning try delete bad event' 
     return
   endif
!   write(*,*) 'event: remove event',n_events

   if ( associated(ev, all_events(n_events)%p) ) then
      deallocate(all_events(n_events)%p%atoms)
      deallocate(all_events(n_events)%p)
      nullify(all_events(n_events)%p)
      n_events = n_events -1
      return
   endif
!     
   i = 1
   do while (i .le. n_events)
      if ( associated(ev, all_events(i)%p) ) then
!         write(*,*) 'event: ', ev%part1,ev%part2
         deallocate(ev%atoms)
         deallocate(ev)
         nullify(ev)

         write(*,*) 'event: del event nimber', i,n_events
         do while (i .lt. n_events)
            all_events(i)%p => all_events(i+1)%p
            i = i+1
         enddo
         write(*,*) 'event:',n_events-1,all_events(n_events-1)%p%atoms, all_events(n_events-1)%p%time
         write(*,*) 'event:',n_events,all_events(n_events)%p%n1,all_events(n_events)%p%fusion, all_events(n_events)%p%time
         nullify(all_events(n_events)%p)
         n_events = n_events -1

         return
      endif 
      i = i+1
   enddo
   write(*,*) 'event: bad pointer. No event'

end subroutine rm_event

!----------------------------------------------
integer function ev_number(ev)

   implicit none
   type(event), pointer :: ev
   integer :: i

   ev_number = 0
   if (.not. associated(ev)) return

   i = 1
   do while (i .le. n_events)
      if ( associated(ev, all_events(i)%p) ) then
        ev_number = i
        return
      endif 
      i = i+1
   enddo
   ev_number = -1 

end function ev_number
!----------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
subroutine grow_all_events(n)

   implicit none
   integer n,i
   type(evptr),  dimension(:), allocatable :: buf
  
   write(*,*) 'event: start grow all_events array', max_events, n_events
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
   write(*,*) 'event: finalize grow all_events array', max_events, n_events

end subroutine grow_all_events
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11111

end module events