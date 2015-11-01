module events

   implicit none

   integer, private :: max_events = 0, n_events = 0

   type event 
!      integer, dimension(:), allocatable :: part1,  part2   ! части кластера
      integer :: n1, n2		!размеры частей n1 >= n2
      integer :: time = 0		!время столкновения или время разваливания
      integer :: fusion = 0		!if grow == 1 - fusion, if grow == -1 - dissociation
      integer :: t_next = 0
      integer :: ref_number = 0
      real    :: e_part1, e_part2, e_tot
      integer, dimension(:), pointer :: atoms => null()   ! все атомы      
   end type event

   type evptr			! ссылка на событие
       type(event), pointer :: p =>null()
   end type evptr  
   
   type(evptr), dimension(:), allocatable :: all_events
 
   private :: grow_all_events

contains

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
subroutine add_event(new_p,new_cl,new_size,old_cl,old_size,time)
!   make an event and add it to events list  
!   new_cl, new_size  - list of atoms in new cluster and it's size
!   old_cl, old_size  - list of atoms in old cluster and it's size
!   time  - time of event
!   new_p - pointer to the new event

   implicit none
   integer,intent(in) :: new_size,old_size
   integer,intent(in), target :: new_cl(new_size),old_cl(old_size), time
   integer, pointer, dimension(:) :: big
   type(event), pointer :: new_p

   integer ::  i,j,k,n1
!   write(*,*) 'event: add dissociation procedure', time, new_cl(1:5),old_cl(1:5)

   if (n_events+1 .ge. max_events) then 
      call grow_all_events(1000)
   endif
   n_events = n_events + 1
   allocate(all_events(n_events)%p)

   if (new_size .gt. old_size) then
     big   => new_cl(:)
     all_events(n_events)%p%fusion = 1
     allocate(all_events(n_events)%p%atoms(new_size))
     n1 = old_size
     all_events(n_events)%p%atoms(1:n1) = old_cl(:)
   else
     big   => old_cl(:)
     all_events(n_events)%p%fusion = -1
     allocate(all_events(n_events)%p%atoms(old_size))
     n1 = new_size
     all_events(n_events)%p%atoms(1:n1) = new_cl(:)
   endif
 !  write(*,*) big,'+',part1
   
   k = n1+1
   j = 1
   do i = 1,size(big)
     if (big(i) == all_events(n_events)%p%atoms(j)) then 
       j = j + 1
     else
       all_events(n_events)%p%atoms(k) = big(i)
       k = k + 1
     endif
   enddo

   all_events(n_events)%p%n1 = n1
   all_events(n_events)%p%n2 = size(big)-n1
   all_events(n_events)%p%time = time
   all_events(n_events)%p%ref_number = size(big)

   new_p => all_events(n_events)%p

   write(*,*) 'event:',n_events, new_p%fusion, ' size1 = ', new_p%n1, ' size2 = ', new_p%n2, new_p%atoms

end subroutine add_event

!----------------------------------------------
subroutine rm_event(ev_target)
!   remove  an event from the events list
!   ev_target - pointer to the event

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
!   return event position in events list
!   ev - pointer to the event

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
subroutine grow_all_events(n)
!   grow all_events array 
!   n - increment 

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
!=================================
subroutine  write_event(p)
! write an event to file and delete it from the list (if there are no references)
! p - pointer to the event
   implicit none
   type(event), pointer :: p

    if (.not. associated(p)) return

! if number of references to the event is 0 then write it and delete it from the list
    p%ref_number = p%ref_number - 1
    if (p%ref_number .gt. 0) return

!    write(*,*) 'hist: write event to file!!!'
    open(40,file='hist.dat',access='append')
    write (40,'(i10, i7, i3, "(",i4," + ", i4,")",3f10.5)') p%time,p%t_next-p%time, p%fusion,p%n1,p%n2,p%e_tot,p%e_part1,p%e_part2
    write (40,*) p%atoms(:)
    close(40)
    call rm_event(p)

end subroutine write_event

end module events