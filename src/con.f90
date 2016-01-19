module con
! calculate concentration of speaces/
!
   implicit none
   integer, private :: max_n = 0, max_dt = 0
   integer, private :: tbegin, tend, dt = 10 !steps

   integer, dimension(:,:), allocatable :: times

contains
!--------------------------------------------
subroutine create_con(n,t,ddt)
! make history array 
! nat - number of atoms
  implicit none 
  integer n,t,ddt

  max_n = n
  max_dt = t
  dt = ddt
  
  allocate(times(max_dt,max_n))
  times(:,:) = 0
  tbegin = - (max_dt-2) * dt
  tend = 0 

end subroutine create_con

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
subroutine update_con(n,nsize)
   implicit none
   integer :: i, nsize
   integer :: n(nsize)

   do i = 1, max_n
     times(1,i) = times(1,i) + n(i)
   enddo
   
   tend = tend + 1
   if (int(tend/dt) .eq. 1.*tend/dt) then 
      
     do i = 1, max_dt-1
       times(max_dt+1-i,:) = times(max_dt-i,:) 
     enddo
     tbegin = tbegin + dt

   endif

end subroutine update_con

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
subroutine get_con(ccon)

   implicit none
   integer :: t, i
   real :: ccon(*)

   t = (tend-tbegin)
   if (tbegin.le.0) t = tend
   do i = 1, max_n
     ccon(i) = (times(1,i)-times(max_dt,i))
     ccon(i) = ccon(i)/t
   enddo
!   write(80,*) tend,tbegin,dt,ccon(1),times(:,1)
!
end subroutine get_con

end module con