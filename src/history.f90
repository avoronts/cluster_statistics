module history
! History of events for atoms/
! for each atom there is an array of pointers to last events
!
!     subroutine create_history(nat)             - create history array
!     subroutine update_history(ev)              - write event to histories of atoms
!     subroutine rm_from_history(ev_target)      - remove event from all histories of atoms and from events array
!     subroutine update_history_check(ev_target) - check for loops and write an event to history
!     integer function loop1(ev)                 - check for complex loops and correct the history

 

   use events

   implicit none
! number of levels in history
   integer, private :: max_hist = 5
! history array for iach atom 
   type (evptr), dimension(:,:), allocatable ::  hist

contains
!--------------------------------------------
subroutine create_history(nat)
! make history array 
! nat - number of atoms
  implicit none 
  integer nat

  allocate(hist(nat,max_hist))

end subroutine create_history

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
subroutine update_history(ev)
! add an event to histories of each atom included in it
! ev - pointer to event in events array 
   implicit none
   integer :: i,n1,iat
   type(event), pointer :: ev 

   if (.not. associated(ev)) return
!    write(*,*) 'hist: Update history' 

   do i = 1, size(ev%atoms)
     iat = ev%atoms(i)
     call write_event(hist(iat,max_hist)%p)

     n1 = max_hist
     do while (n1 .gt. 1) 
        hist(iat,n1)%p => hist(iat,n1-1)%p
        n1 =n1-1
     enddo
     hist(iat,n1)%p => ev
     if (associated(hist(iat,2)%p))  hist(iat,2)%p%t_next = ev%time
   enddo 

end subroutine update_history

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
subroutine rm_from_history(ev_target)
! remove an event from from events array and from histories of each atom included in it
! ev_target - pointer to event in events array 

   implicit none
   type(event), pointer :: ev,ev_target
   integer :: jat,i,iat,n,n1,j

   ev => ev_target      ! make copy
   if (.not. associated(ev)) return
!   write(*,*) 'hist: remove history ',ev%fusion, ev%n1, ev%n2, ev%atoms(:),ev_number(ev)
!   ev => hist(jat,n)%p

   do i = 1, ev%n1 + ev%n2
     iat = ev%atoms(i)
     n1 = 1
!     write(*,*) 'hist: remove history from atom ',iat,(ev_number(hist(iat,j)%p),j=1,5)
     do while (n1 .le. max_hist) 
        if (associated(ev,hist(iat,n1)%p)) then
          do while (n1 .lt. max_hist) 
            hist(iat,n1)%p => hist(iat,n1+1)%p
            n1 =n1+1
          enddo
          hist(iat,n1)%p => NULL()
          exit
        endif
        n1 = n1+1
     enddo
   enddo 
   call rm_event(ev)
!  write(*,*) 'hist: done remove history '
!
end subroutine rm_from_history

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>
subroutine update_history_check(ev_target)
! check a new event for loops and update the hisory array 
! ev_target - pointer to event in events array 

    implicit none
    integer :: iat1,iat2
    type(event), pointer :: ev,ev_target

    ev => ev_target
    if (ev%fusion .eq. -1) then
        call update_history(ev)
        return
    endif

!   write(*,*) 'hist: add dissociation procedure', time, new_cl(1:5),old_cl(1:5)
    iat1 = ev%atoms(1)
    iat2 = ev%atoms(ev%n1+1)

!------ check for loops ---------
    if ( associated(hist(iat1,1)%p,hist(iat2,1)%p) ) then 
       write(*,*) 'stat: Simple loop. Remove last step (hist(j))'
        !   ev => hist(iat1,1)%p
       call rm_from_history(hist(iat1,1)%p)
       call rm_event(ev)
       return
    endif

    if ( loop1(ev) .eq. 1 ) then 
       return
    endif

    call update_history(ev)

end subroutine update_history_check

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
integer function loop1(ev)
! check a new event for complex loop in history
! (if two parts of a cluster were separated and having a collision at that time)
! ev - pointer to event in events array 

   implicit none
   type(event), pointer :: ev, fake_ev, true_ev
   integer :: i,iat, iat1,iat2 , n0,n1,n2
   integer, allocatable :: all(:), part(:)
   
   integer, allocatable :: buf1(:), buf2(:), add_atoms(:)
   

   loop1 = 0
   if (.not. associated(ev)) return

   iat1 = ev%atoms(1)
   iat2 = ev%atoms(ev%n1+1)
   
   if (.not. associated(hist(iat1,1)%p) ) return
   if (.not. associated(hist(iat2,1)%p) ) return
   if (.not. associated(hist(iat1,2)%p) ) return
   if (.not. associated(hist(iat2,2)%p) ) return

   if ( associated(hist(iat1,1)%p,hist(iat2,2)%p) ) then 
       fake_ev => hist(iat1,1)%p
       allocate(add_atoms(ev%n1))
       add_atoms = ev%atoms(1:ev%n1)
       true_ev => hist(iat2,1)%p
   endif

   if ( associated(hist(iat1,2)%p,hist(iat2,1)%p) ) then 
       fake_ev => hist(iat2,1)%p
       allocate(add_atoms(ev%n2))
       add_atoms = ev%atoms(ev%n1+1:)
       true_ev => hist(iat1,1)%p
   endif

   if (.not. allocated(add_atoms)) return

   write(*,*) 'hist: Complex loop. reconstruction of history',add_atoms,'+', true_ev%atoms
   return
   
   loop1 = 1
!   read(*,*)

   n1 = true_ev%n1
   n2 = true_ev%n2
   n0 = size(add_atoms(:))
   write(*,*) 'hist: size',n1, n2, n0

   allocate( buf1(n1), buf2(n2))
   buf1(:) = true_ev%atoms(1:n1)
   buf2(:) = true_ev%atoms(n1+1 : )
   write(*,*) 'hist: buf',buf1(:),buf2(:)

   deallocate(true_ev%atoms)
   allocate(true_ev%atoms(n1 + n2 + n0 ))

   if ( associated( hist(true_ev%atoms(1),2)%p, fake_ev ) ) then 
       ! add to part 1
       write(*,*) 'hist: add to part 1'
       true_ev%atoms(1 : n1) = buf1(:)
       true_ev%atoms(n1+1 : ) = add_atoms(:)
       true_ev%n1 = n1 + n0
       true_ev%atoms(n1 + n0 + 1 : ) = buf2(:)
   else 
      ! add to part 2
       write(*,*) 'hist: add to part 2'
       true_ev%atoms(1:n1) = buf1(:)
       true_ev%atoms(n1+1 : ) = buf2(:)
       true_ev%atoms(n1 + n2 + 1 : ) = add_atoms(:)
       true_ev%n2 = n2 + n0
   endif
   deallocate(buf1,buf2)
   
   call rm_from_history(fake_ev)

!  update history for new atoms
   write(*,*) 'hist: Complex loop. reconstruction of history'

   do i=1,size(add_atoms)

      iat = add_atoms(i)
      n1 = max_hist
      do while (n1 .gt. 1) 
         hist(iat,n1)%p => hist(iat,n1-1)%p
         n1 =n1-1
      enddo
      hist(iat,n1)%p => true_ev
      if (associated(hist(iat,2)%p))  hist(iat,2)%p%t_next = true_ev%time

   enddo
   write(*,*) 'hist: Complex loop. reconstruction of history'

   call rm_event(ev)
   write(*,*) 'hist: Complex loop. reconstruction of history'

   deallocate(add_atoms)
end function loop1
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
end module history