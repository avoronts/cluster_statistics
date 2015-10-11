!!!!!!!!!!!!my programm !!!!!!!!!!!!!
program simple

  use LAMMPS_gether
  use history
  
  implicit none
  

integer nstat,max_atom,nclu,n_entre
parameter(nstat=125,max_atom=25001,nclu=2500,n_entre=15)

!this data for exchange with LAMMPS program

    !integer me - from LAMMPS_gether module
    !integer natoms - from LAMMPS_gether module
    !integer int_buf - from LAMMPS_gether module
    !integer dbl_buf - from LAMMPS_gether module
    integer  typ(max_atom), num_atom, num_vecino4(max_atom), num_vecino5(max_atom)
    double precision pot(max_atom), kin(max_atom), vx(max_atom),vy(max_atom),vz(max_atom), & 
                     x(max_atom),y(max_atom),z(max_atom)

!this piece of text is devoted to statistical relations betw...
    integer*4 statist(300), npromediar /300/,delta/5/, delta_tau, nstatist /10/
    real*8 mpro,dpro,a
    integer curnum1 /2/
    logical :: sexist

!
    character*200 sfile 
    character*150 sline

! this data for cluster storage
    integer*4 ncluster(max_atom),ncluster_o(max_atom),mcluster(0:nclu), & 
              mcluster_o(0:nclu)  ! 1- только металл, 2 - есть аргон
    integer cluster(0:nclu,1:nstat+1),cluster_o(0:nclu,1:nstat+1)


!    integer*4 numj_all(1000), numt_t(1000),numj_o(1000),numjt(1000),numt_o(1000),numt_all(1000)
!    integer   knumj_all, knumt_t, knumj_o, knumjt, knumt_o, knumt_all
    
    integer cyp(max_atom)

    integer*4 i, j, jj, i1, k, nt, cur_clu , icl, iat, iat1, iat2
    integer t, p, n, itwo,ipp, itres, ityp, ikk,ip,iexit,ifile,ik,n1,n2,n3
        
    integer*4  nmax /10000/,natraso /5/,sostav(500)
    !integer cluster7_o(1:nclu,1:nstat),cluster7(1:nclu,1:nstat),mcluster7(0:nclu,2),ncluster7(1:natom),mcluster7_o(0:nclu,2),ncluster7_o(1:natom)

    type (evptr), dimension(max_atom,5) ::  hist
    type (event), pointer :: nev
     
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!--------------------------------------------------------
CALL init_mpi_lammps('in.1') !!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  
nt=0
mcluster(:) = 0        !numero de atoms en cada cluster
ncluster(:) = 0        !numero de cluster que mantiene este atom 
cluster(:,:) = 0          !numero del atom metallico

!!me=1 ! DON'T forget to borrar it!
  ! do something
do i=1,5000
  typ(i)=1
enddo    
do i=5001,max_atom
   typ(i)=2
enddo

open(150,file='stat2.dat')
open(44000,File='r1.dat')      

do while(.true.)
! entro y esta pensando
!								if (me.eq.1) then 

!!!!!!!!!!!!!!!!!!!!!!! read _data !!!!!!!!!!!!!!!!!!!!!!!!!
!  do j=1,m
!    num_atom,typ(num_atom),num_vecino4(num_atom),num_vecino5(num_atom),kin(num_atom), pot(num_atom)  !x(num_atom),y(num_atom),z(num_atom),
!  enddo
    call one_step_data(5)
    if (me .ne. 0) cycle

        j=1
        jj=1
        do i=1,natoms
          num_atom=int_buf(j)
          typ(num_atom) = int_buf(j+1) 
          num_vecino4(num_atom) = int_buf(j+2) 
          num_vecino5(num_atom) = int_buf(j+3)
          j=j+4
          
          kin(num_atom) = dbl_buf(jj)
          pot(num_atom) = dbl_buf(jj+1)
	  x(num_atom) = dbl_buf(jj+2)
	  y(num_atom) = dbl_buf(jj+3)
	  z(num_atom) = dbl_buf(jj+4)
	  vx(num_atom) = dbl_buf(jj+5)
	  vy(num_atom) = dbl_buf(jj+6)
	  vz(num_atom) = dbl_buf(jj+7)
!          print  *, num_atom,x(num_atom),y(num_atom),z(num_atom),vx(num_atom),vy(num_atom),vz(num_atom)
       
          jj=jj+8

        enddo
        !print  *, '<<<<<<<<<<< read done  >>>>>>>>>>>>>>>>>'
        nt=nt+1

!!!!!!!!!!!!!!!!!!!!! check 1. my code 21.05.2015 !!!!!!!!!!!!!!!!1
   open(5050,File='qqq')    !'//trim(sfile)//'
   write(5050,*) nt			!real*8 kin(1: natom)
   close(5050)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11


!!!!!!!!!!!!!!!!!!!!! check 1. my code 21.05.2015 !!!!!!!!!!!!!!!!1
!   write(sfile,'(i0)') (nt+1)*natraso
!
!   sfile=trim('dump/check_in'//trim(sfile)//'.bin')
!   write(*,*) sfile
!   open(5050,File=sfile,form='unformatted')    !'//trim(sfile)//'
!
!   write(5050) typ 			!integer typ(1:natom)
!   write(5050) num_vecino5		!integer*4 num_vecino5(1:natom)
!   write(5050) kin			!real*8 kin(1: natom)
!   close(5050)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if(int(1.*nt/nmax).eq. 1.*nt/nmax) then
   close(44000)

   !<<<<<<<<<<<<<<< statistics >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   if((int(1.*(nt)/nstatist/nmax).eq. 1.*(nt)/nstatist/nmax) .and. nt .ne. 0) then
   !if((int(1.*(nt-1)/nstatist/nmax).eq. 1.*(nt-1)/nstatist/nmax) .and. nt .ne. 1) then

     statist=0 
     a=0.
     j=int(1.*nt/nmax)
     do   ifile=j-nstatist + 1, j, 1
       write(sfile,'(i0)') ifile
       sfile=('r'//trim(sfile)//'.dat')
       inquire (file=sfile,exist=sexist)
       if(sexist) then
         open(18,File=sfile,err=12) 

         do while(.true.)    

        !   read(18,('(2(i12,1x),a12,1x,a12,1x,3(i12,1x),125(i5,1x))'),end=11 ) ityp,n,sline,sline, n3, n1,n2, & 
        !             (numjt(ik), ik=1,12)
           if(sline.ne.'*******') then
               read(sline,'(i12)') delta_tau
           endif

           if (curnum1==n .and. (ityp==1 .or. ityp==-1) .and. sline.ne.'*******') then    !ityp==1 .or. ityp==-1)
               a=a+1.
               do i=1,npromediar
                  if (delta_tau>=(i-1.)*delta .and. delta_tau < i*delta) then
                     statist(i)=statist(i)+1
                     exit
                  endif 
               enddo 
           endif

         enddo
     11  continue
     
       endif  ! if(sexist) then
     enddo   !from  do   ifile=j-nstatist + 1, j, 1
  12 continue

     Mpro=0.
     dpro=0.
     do i=1,npromediar
        Mpro=Mpro+statist(i)*(i*delta+delta/2.)/a
        dpro=dpro+statist(i)*(i*delta+delta/2.)**2/a
     enddo
     dpro=dpro-mpro**2
      ! write(14, '( 425(i5,1x))' ) (sargon(i),i=1,120)
     write(150,('(i5,1x,3(F10.3, 1x), 425(i5,1x))')) ifile,a,mpro,dpro,(statist(i),i=1,npromediar)

   endif    !from if(int(1.*nt/nstatist).eq. 1.*nt/nstatist) then
!<<<<<<<<<<<<<<<<<<<<<<<<<<<< end statistics >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

   write(sfile,'(i0)') int(nt/nmax)+1
   !open(41000,File='osc'//trim(sfile)//'.dat') 
   !open(15,File='t'//trim(sfile)//'.dat')      
   open(44000,File='r'//trim(sfile)//'.dat')      
   !open(42000,File='pro_argon'//trim(sfile)//'.dat')
   ! open(43000,File='pro_argon_menos'//trim(sfile)//'.dat')      
   ! open(47000,File='prot'//trim(sfile)//'.dat')      
endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!                do ifile=1,nmax
!										    if (me.eq.1) then  
 mcluster_o(:) = mcluster(:)
 ncluster_o(:) = ncluster(:)
 cluster_o(:,:) = cluster(:,:)

 mcluster(:) = 0
 ncluster(:) = 0
 cluster(:,:) = 0
 mcluster(0) = 1     ! нулевой кластер описывает все мономеры (размер 1) 
 cur_clu=0


 do iat = 1,natoms
     if (typ(iat) .ne. 1) cycle
     if (iat .eq. num_vecino5(iat)) cycle

     icl = Ncluster(num_vecino5(iat))
!     write(*,*) iat, num_vecino5(iat),typ(iat),icl
     if (icl .eq. 0) then
            cur_clu = cur_clu+1    !nuevo numero
            icl = cur_clu
            if (typ(num_vecino5(iat)) .eq. 1) then   ! add base atom
               Ncluster(num_vecino5(iat)) = icl
               mcluster(icl) = mcluster(icl)+1
               cluster(icl,mcluster(icl)) = num_vecino5(iat)
            endif
     endif

     Ncluster(iat) = icl
     if (mcluster(icl) .lt. nstat) then
        mcluster(icl) = mcluster(icl)+1
        cluster(icl,mcluster(icl)) = iat
     endif
    ! write(*,*) iat, num_vecino5(iat),icl,Ncluster(iat),mcluster(icl)
     
 enddo    !!num_atom==1,m

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
if(int(1.*nt/1000).eq. 1.*nt/1000) write(*,*) 'stat: time = ', nt

cyp(:) = 0
do j=1,5000  ! todos los atoms

  if (nt .lt. 10) cycle
  if (cyp(j) .eq. 1) cycle      ! skip atoms if they are in known clusters
  if (typ(j) .ne. 1) cycle
  if (mcluster(ncluster(j)) .eq. mcluster_o(ncluster_o(j))) cycle
  
  if (mcluster(ncluster(j)) < mcluster_o(ncluster_o(j))) then 
!        write(*,*) 'stat: ! add dissociation'
        if (ncluster(j) .eq. 0) cluster(ncluster(j),1) = j  ! мономер
        
        call add_diss(nev,cluster(ncluster(j),:), cluster_o(ncluster_o(j),:), nt)
        do i=1,mcluster_o(ncluster_o(j))
          iat = cluster_o(ncluster_o(j),i)
          if (associated(hist(iat,5)%p)) then
             write(*,*) 'stat: diss: write event to file!!!'
             open(40,file='hist.dat',access='append')
             write (40,*) 'diss: time = ', hist(j,5)%p%time,', (',hist(iat,5)%p%n1,'+',hist(iat,5)%p%n2,')', & 
                hist(iat,5)%p%part1(:),  hist(iat,5)%p%part2(:)
             close(40)
             call rm_event(hist(iat,5)%p)
          endif  
          hist(iat,5)%p => hist(iat,4)%p
          hist(iat,4)%p => hist(iat,3)%p
          hist(iat,3)%p => hist(iat,2)%p
          hist(iat,2)%p => hist(iat,1)%p
          hist(iat,1)%p => nev
          cyp(iat) = 1
        enddo
        
  endif

  if (mcluster(ncluster(j)) > mcluster_o(ncluster_o(j))) then 
!        write(*,*) 'stat: Add fusion',j
        if (ncluster_o(j) .eq. 0) cluster_o(ncluster_o(j),1) = j  ! мономер
        call add_fusion(nev,cluster(ncluster(j),:), cluster_o(ncluster_o(j),:), nt)
    
!------ check for loops ---------
        if (loop0(nev,hist(j,1)%p) == 1) then 
           write(*,*) 'stat: Simple loop. Remove last step (hist(j))'
           call rm_event(hist(j,1)%p)
           do i = 1, mcluster(ncluster(j))
              iat = cluster(ncluster(j),i)
              hist(iat,1)%p => hist(iat,2)%p
              hist(iat,2)%p => hist(iat,3)%p
              hist(iat,3)%p => hist(iat,4)%p
              hist(iat,4)%p => hist(iat,5)%p
              hist(iat,5)%p => NULL()
           enddo
        else
           iat1 = nev%part1(1)
           iat2 = nev%part2(1)
           if (loop1(nev,hist(iat1,1)%p,hist(iat2,1)%p) == 1) then
               write(*,*) 'stat: Complex loop. reconstruction of history'
           
           else
              write(*,*) 'stat: No loops. Write new history'
              do i = 1,mcluster(ncluster(j))
                 iat = cluster(ncluster(j),i)
                 if (associated(hist(iat,5)%p)) then
                    write (*,*) 'stat: fusion: write event to file!!!'
                    open(40,file='hist.dat',access='append')
                    write (40,*) 'fusion: time = ', hist(j,5)%p%time,', (',hist(iat,5)%p%n1,'+',hist(iat,5)%p%n2,')', & 
                        hist(iat,5)%p%part1(:),  hist(iat,5)%p%part2(:)
                    close(40)
                    call rm_event(hist(iat,5)%p)
                 endif  
                 hist(iat,5)%p => hist(iat,4)%p
                 hist(iat,4)%p => hist(iat,3)%p
                 hist(iat,3)%p => hist(iat,2)%p
                 hist(iat,2)%p => hist(iat,1)%p
                 hist(iat,1)%p => nev
              enddo
           endif
        endif
 
        do i = 1,mcluster(ncluster(j))
           iat = cluster(ncluster(j),i)
           cyp(iat) = 1
        enddo
  endif

enddo   !form j=1,5000

enddo    !from do while .true. <<<<<< end of main circle


call finalize_mpi_lammps()
10 continue

end program 

