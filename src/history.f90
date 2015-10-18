module history

   use events
 
   implicit none

   integer, private :: max_hist = 5
   type (evptr), dimension(:,:), allocatable ::  hist

contains
!--------------------------------------------
subroutine create_history(nat)

  implicit none 
  integer nat

  allocate(hist(nat,max_hist))
  

end subroutine create_history

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
subroutine update_history(ev)

   implicit none
   integer :: i,n1,iat
   type(event), pointer :: ev 
   
   if (.not. associated(ev)) return
!    write(*,*) 'hist: Update history' 

   do i = 1, size(ev%atoms)
     iat = ev%atoms(i)
     call write_hist(hist(iat,max_hist)%p)
!   call rm_old_history(iat,max_hist)

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
recursive subroutine rm_old_history(ev)

   implicit none
   type(event), pointer, intent(in) :: ev
   integer :: jat,i,iat,n,n1
   
   if (.not. associated(hist(jat,n)%p)) return
   write(*,*) 'hist: remove history ',n,' atom ',jat
   
   do i = 1, ev%n1+ev%n2
     iat = ev%atoms(i)
     n1 = 1
     do while (n1 .le. max_hist) 
        if (hist(iat,n1)%p%time .lt. ev%time) then
          call rm_old_history(hist(iat,n1)%p)
        endif
     enddo
   enddo 

   call rm_from_history(hist(iat,n1)%p)

end subroutine rm_old_history

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>
subroutine update_history_check(ev_target)

    implicit none
    integer :: iat1,iat2
    type(event), pointer :: ev,ev_target

    ev => ev_target
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

   write(*,*) 'hist: Complex loop. reconstruction of history',add_atoms
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

!=================================
subroutine  write_hist(p)
   
   implicit none
   type(event), pointer :: p
       
    if (.not. associated(p)) return
    if (p%written .eq. 1) return
    
!    write(*,*) 'hist: write event to file!!!'
    open(40,file='hist.dat',access='append')
    write (40,'(i10, i7, i3, "(",i4," + ", i4,")",3f10.5)') p%time,p%t_next-p%time, p%fusion,p%n1,p%n2,p%e_tot,p%e_part1,p%e_part2
    write (40,*) p%atoms(:)
    close(40)
    p%written = 1
       
end subroutine write_hist
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

end module history