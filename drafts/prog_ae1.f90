program simple

   use MPI
   use LAMMPS
   use LAMMPS_gether
   use, intrinsic :: ISO_C_binding, only : C_double, C_ptr, C_int

   implicit none

   integer :: ierr, me, nprocs    ! for MPI

   type (C_ptr) :: lmp		! for lammps communication
   integer :: lammps_fl, nprocs_main, nprocs_lammps, comm_lammps
   integer, dimension(:), allocatable :: buf
   integer :: buf_step, natoms, loop, nloop, mv

interface   one_step_data  ! interface  for geting info from lammps nodes
  subroutine one_step_data_simple (data, natoms, lmp, me,lammps_fl,nprocs)
     use, intrinsic :: ISO_C_binding, only : C_double, C_ptr, C_int
     implicit none
     integer , dimension(:), allocatable,  intent(out) :: data
     type (C_ptr), intent(in) :: lmp
     integer, intent(in) :: me, lammps_fl,nprocs
     integer, intent(out) :: natoms
 end subroutine one_step_data
end interface one_step_data


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!<<<<<<<<<<<<<<<< put initial block here >>>>>>>>>>>>>>>>>>


!   include 'mpif.h'
integer nstat,nclu,natom, n_entre
parameter(nstat=125,natom=25001,nclu=2500,n_entre=15)
!!!   integer :: ierr, me, nprocs,natoms
!!!  type (lammps_instance) :: lmp
!!!   double precision, dimension(:), allocatable :: compute_v, mass, r,xx
!!!   double precision, dimension(:,:), allocatable :: comp
   real, dimension(:,:), allocatable :: x_r

integer ik,ip,ikk,ikp, itwo, itres, nttfile,ifile
integer  knumjt, knumj_o,  knumj_all, knumt_o, knumt_all, knumt_t, ipp, iexit
double precision en_arg, t,p
          

logical exists
    real*4  statinum(nstat), para_esc(nstat),pott,kinn,kin_entre(n_entre),l_entre(n_entre)   !x(natom),y(natom),z(natom),
    integer eof /0/,typ(1:natom),cyp(natom),i1,i2,i3,k,Code,cur_rec ,iras,np(12),lstat(nstat),iii
    character,allocatable:: struu(:)
    character*5 s1
    character*6 s2,c2
    character*200 sfile 
    character*150 sline
    integer*4 num_vecino4(1:natom),num_vecino5(1:natom),ncluster(1:natom),mcluster(0:nclu,2),ncluster_o(1:natom),mcluster_o(0:nclu,2),lcluster(0:nclu)  ! 1- ������ ������, 2 - ���� �����
integer*8 n1, n2, iarg
    integer*4 i,j,num_atom,nt, M, cur_num ,cur_clu,cur_clu7, numj_all(1000), numt_t(1000),numj_o(1000),numjt(1000),numt_o(1000),numt_all(1000)
    integer*4  nmax /10000/,natraso /5/,sostav(500)
  
    real*8 pot(1 :natom),kin(1: natom),l2
character*17 sdir/'c:\1500_1\for_me\'/
    integer cluster(1:nclu,1:nstat),cluster_o(1:nclu,1:nstat)
    integer cluster7_o(1:nclu,1:nstat),cluster7(1:nclu,1:nstat),mcluster7(0:nclu,2),ncluster7(1:natom),mcluster7_o(0:nclu,2),ncluster7_o(1:natom)
        type hist
            integer*4 ntime  !����� ������������ � �.�. ! ��� ����� ������������
            integer*4 n1  !� ������ ���� ������ ���� J   ! �� ��� ��������
            integer*4 n2 !� ����� �������� ����������
            integer*4 narg ! ������� ������������ � ������� ����, ������� � ntime
real*4 stat(nstat)   ! ��� �� �������������� � �������
integer cluster(nstat) !� ����� ������� ������ ������
        end type
    type (hist) hist_dec(25000), hist_uni(25000)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
mcluster_o=0        !numero de atoms en cada cluster
ncluster_o=0        !numero de cluster que mantiene este atom 
cluster_o =0          !numero del atom metallico

  nttfile=1

!--------------------------------------------------------
  call init_all('in.1',lmp,me,nprocs,lammps_fl)
!--------------------------------------------------------
nt=1
!!me=1 ! DON'T forget to borrar it!
  ! do something
          do i=1,5000
          typ(i)=1
          enddo
      do i=5001,natom
      typ(i)=2
      enddo
 ikk=0
ikp=0
      
										do iii =1,10  !while(.true.)
! entro y esta pensando
										if (me.eq.0) then 
    do i=0,20
 if (int(nttfile/(10.**i)) .eq.0) then
 write(sfile,'(i1)') nttfile
 exit
 endif
 enddo
 if (i==1) then
 write(sfile,'(i1)') nttfile
 endif
  if (i==2) then
 write(sfile,'(i2)') nttfile

 endif
  if (i==3) then
 write(sfile,'(i3)') nttfile
 endif
  if (i==4) then
 write(sfile,'(i4)') nttfile
 endif
  if (i==5) then
 write(sfile,'(i5)') nttfile
 endif
  if (i==6) then
 write(sfile,'(i6)') nttfile
 endif
  if (i==7) then
 write(sfile,'(i7)') nttfile
 endif
  if (i==8) then
 write(sfile,'(i8)') nttfile
 endif
 open(41000,File='osc'//trim(sfile)//'.dat') 
 !open(15,File='t'//trim(sfile)//'.dat')      
 open(44000,File='r'//trim(sfile)//'.dat')      
open(42000,File='pro_argon'//trim(sfile)//'.dat')
! open(43000,File='pro_argon_menos'//trim(sfile)//'.dat')      
! open(47000,File='prot'//trim(sfile)//'.dat')      
 
                do ifile=1,nmax
  
  do i=0,20
 if (int(nt*natraso/(10.**i)) .eq.0) then
 exit
 endif
 enddo
 if (i==1) then
 write(sfile,'(i1)') nt*natraso
 endif
  if (i==2) then
 write(sfile,'(i2)') nt*natraso

 endif
  if (i==3) then
 write(sfile,'(i3)') nt*natraso
 endif
  if (i==4) then
 write(sfile,'(i4)') nt*natraso
 endif
  if (i==5) then
 write(sfile,'(i5)') nt*natraso
 endif
  if (i==6) then
 write(sfile,'(i6)') nt*natraso
 endif
  if (i==7) then
 write(sfile,'(i7)') nt*natraso
 endif
  if (i==8) then
 write(sfile,'(i8)') nt*natraso
 endif
   if (i==9) then
 write(sfile,'(i9)') nt*natraso
 endif
   if (i==10) then
 write(sfile,'(i10)') nt*natraso
 endif
   if (i==11) then
 write(sfile,'(i11)') nt*natraso
 endif
   if (i==12) then
 write(sfile,'(i12)') nt*natraso
 endif
   if (i==13) then
 write(sfile,'(i13)') nt*natraso
 endif
   if (i==14) then
 write(sfile,'(i14)') nt*natraso
 endif
   if (i==15) then
 write(sfile,'(i15)') nt*natraso
 endif
   if (i==16) then
 write(sfile,'(i16)') nt*natraso
 endif
   if (i==17) then
 write(sfile,'(i17)') nt*natraso
 endif 
  if (i==18) then
 write(sfile,'(i18)') nt*natraso
 endif
   if (i==19) then
 write(sfile,'(i19)') nt*natraso
 endif
   if (i==20) then
 write(sfile,'(i20)') nt*natraso
 endif
   if (i==21) then 
 write(sfile,'(i21)') nt*natraso
 endif
!! sfile=trim('dmp_clust'//trim(sfile)//'.txt')
 sfile=trim('dmp_clust'//trim(sfile)//'.txt')
    mcluster_o= mcluster
    ncluster_o=ncluster
    cluster_o=cluster
      mcluster7_o= mcluster7
    ncluster7_o=ncluster7
    cluster7_o=cluster7          

mcluster=0
mcluster7=0



iras=0

   open(50000,File=sfile,err=10)    !'//trim(sfile)//'
 do i =1,3
 read(50000,*) sline     !lea las demas lineas no necesarias
 enddo
  Read(50000,*) M  !����� ������
  do i =5,9
 read(50000,*) sline     !lea las demas lineas no necesarias
 enddo
         !lea la primera linea
      ! do while (.not. EOF(1))
          do j=1,m
       Read(50000,*) num_atom,typ(num_atom),num_vecino4(num_atom),num_vecino5(num_atom),kin(num_atom), pot(num_atom)  !x(num_atom),y(num_atom),z(num_atom),
       enddo
       close(50000,status='delete')
! first we have to determine the kinetic energy        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!call next_data()
!num_atom,typ(num_atom),num_vecino4(num_atom),num_vecino5(num_atom),kin(num_atom), pot(num_atom)

      ncluster=0
      cluster=0
      cluster7=0
      ncluster7=0
      cur_clu=1
      cur_clu7=1
      ik=0
cyp=0

do num_atom=1,m
    if (num_vecino4(num_atom).ne.num_atom ) then   
            
             if (Ncluster7(num_vecino4(num_atom))==0) then
                    Ncluster7(num_vecino4(num_atom))=cur_clu7
                   
                    mcluster7(Ncluster7(num_vecino4(num_atom)),1)=1
                    !if (typ(num_atom).ne. typ(num_vecino4(num_atom)) )  then 
                    !mcluster(ncluster(num_vecino4(num_atom)),2)=2
                    !endif
                    cur_clu7=cur_clu7+1    !nuevo numero
                   
                    if (cluster7(ncluster7(num_vecino4(num_atom)),1) .eq. 0) cluster7(ncluster7(num_vecino4(num_atom)),1)=num_vecino4(num_atom)
              endif
                    Ncluster7(num_atom)= Ncluster7(num_vecino4(num_atom))
                   
                    i=ncluster7(num_vecino4(num_atom))
                    
                        if (mcluster7(Ncluster7(num_vecino4(num_atom)),1)<nstat) then

                        mcluster7(Ncluster7(num_vecino4(num_atom)),1)=mcluster7(Ncluster7(num_vecino4(num_atom)),1)+1
                        !if (typ(num_atom).ne. typ(num_vecino4(num_atom)) )  then 
                        !mcluster(ncluster(num_vecino4(num_atom)),2)=2
                        !else
                        !mcluster(Ncluster(num_vecino4(num_atom)),2)=1 
                        !endif

                        cluster7(i,mcluster7(Ncluster7(num_vecino4(num_atom)),1))=num_atom
                        endif
    endif   
 ! el chiquito
            
         if (num_vecino5(num_atom).ne.num_atom ) then   !quitamos the clusteres of argon solo .and. num_vecino5(num_atom)<5000 
            
             if (Ncluster(num_vecino5(num_atom))==0) then
                    Ncluster(num_vecino5(num_atom))=cur_clu
                    lcluster(cur_clu)=1
                    mcluster(Ncluster(num_vecino5(num_atom)),1)=1
                    !if (typ(num_atom).ne. typ(num_vecino5(num_atom)) )  then 
                    !mcluster(ncluster(num_vecino5(num_atom)),2)=2
                    !endif
                    cur_clu=cur_clu+1    !nuevo numero
                    if (cluster(ncluster(num_vecino5(num_atom)),1) .eq. 0) cluster(ncluster(num_vecino5(num_atom)),1)=num_vecino5(num_atom)
              endif
                    Ncluster(num_atom)= Ncluster(num_vecino5(num_atom))
                    cyp(num_vecino5(num_atom))=10
                    i=ncluster(num_vecino5(num_atom))
if (typ(num_atom) .ne.2) lcluster(i)=lcluster(i)+1
                if (mcluster(Ncluster(num_vecino5(num_atom)),1)<nstat) then
                mcluster(Ncluster(num_vecino5(num_atom)),1)=mcluster(Ncluster(num_vecino5(num_atom)),1)+1
                !if (typ(num_atom).ne. typ(num_vecino5(num_atom)) )  then 
                !mcluster(ncluster(num_vecino5(num_atom)),2)=2
                !else
                !mcluster(Ncluster(num_vecino5(num_atom)),2)=1 
                !endif
                cluster(i,mcluster(Ncluster(num_vecino5(num_atom)),1))=num_atom
                endif
          endif
 ! else 
        
          !if (Code <> 0)
        enddo    !!num_atom==1,m
        lstat=0
    do i=1,cur_clu-1   ! ��� ����� ��������� �� ������ ������
    lstat(lcluster(i))=lstat(lcluster(i))+1
    enddo         
         
   
do i=1,natom
if (cyp(i) .ne. 10)   cyp(i)=0
enddo

        do j=1,5000  ! todos los atoms

   
  numjt=0
  numj_o=0
  numj_all=0
  numt_o=0
  numt_all=0
  knumjt=0
  knumj_o=0
  knumj_all=0
  knumt_o=0
  knumt_all=0 
numt_t=0
knumt_t=0
  
        if (cyp(j).ne.1 ) then     !.and.cyp(j).ne.1
        !comparar las lineas
  
    if (ncluster7(j) .ne. 0) then 
                    ip=1
                    do k=1,mcluster7(ncluster7(j),1)
                    if (typ(cluster7(ncluster7(j),k)) .ne.2) then
                    numjt(ip)=cluster7(ncluster7(j),k)         
                    ip=ip+1
                    endif
                    enddo   
                        else
                        numjt(1)=j
                        ip=2
                        endif
            knumjt=ip-1
!� numjt ����� ��� ����� �� cluser(ncluster7(j)), ����� ������ ������
        if (ncluster7_o(j) .ne. 0) then 
                    ip=1
                    do k=1,mcluster7_o(ncluster7_o(j),1)
                    if (typ(cluster7_o(ncluster7_o(j),k)) .ne.2) then
                    numj_o(ip)=cluster7_o(ncluster7_o(j),k)         
                    ip=ip+1
                    endif
                    enddo   
        else
                        numj_o(1)=j
                        ip=2
        endif
            knumj_o=ip-1 
            !!!!!!��� ����������� � �������!!!\!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       if ( knumjt == knumj_o .and. mcluster7(ncluster7(j),1) > mcluster7_o(ncluster7_o(j),1) ) then   
       !����� ����� ������, ��� ����� ������� ��.

         do i=1,mcluster7(ncluster7(j),1)   !�����
      hist_uni(cluster7(ncluster7(j),i))%stat(1)=kin(cluster7(ncluster7(j),i))
      hist_uni(cluster7(ncluster7(j),i))%stat(2)=pot(cluster7(ncluster7(j),i))
         cyp( cluster7(ncluster7(j),i))=1  
        enddo   

            do i=1,knumjt
            hist_uni(numjt(i))%narg= hist_uni(numjt(i))%narg+1
            enddo
            
       endif     !��� �����������  � �������!!!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !��� ���������� � �������!!!\!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

    if ( knumjt == knumj_o .and. knumjt.ne.1 .and. mcluster7(ncluster7(j),1) < mcluster7_o(ncluster7_o(j),1) ) then   !��� ���������� � �������!!!\
en_arg=0.
         do i=1,mcluster7_o(ncluster7_o(j),1)   !������
!if (typ(cluster7_o(ncluster7_o(j),i))==2) then
         iexit=0
           
                 do ik=1, mcluster7(ncluster7(j),1)
            if (cluster7_o(ncluster7_o(j),i)==cluster7(ncluster7(j),ik)) then 
            iexit=iexit+1 ! encontramos
            exit
            endif
                enddo
    if (iexit == 0) then
    en_arg= en_arg+pot(cluster7_o(ncluster7_o(j),i))+kin(cluster7_o(ncluster7_o(j),i))-(hist_uni(cluster7_o(ncluster7_o(j),i))%stat(1)+hist_uni(cluster7_o(ncluster7_o(j),i))%stat(2))
    endif
!endif               
    
            enddo
                         if (ncluster(j) .ne. 0) then 
                    ip=1
                    do k=1,mcluster(ncluster(j),1)
                    if (typ(cluster(ncluster(j),k)) .ne.2) then
                           
                    ip=ip+1
                    endif
                    enddo   
                        else
                        
                        ip=2
                        endif
            knumjt=ip-1
  if ( en_arg .ne. 0 .and. knumjt>1 )   write(42000,('(5(i10,1x),f10.7,1x,115(i7,1x))') ) nt,cluster7_o(ncluster7_o(j),1), knumjt,mcluster7_o(ncluster7_o(j),1),mcluster7(ncluster7(j),1),en_arg, (lstat(ip),ip=1,15)

        do i=1,mcluster7_o(ncluster7_o(j),1)   
        hist_uni(cluster7_o(ncluster7_o(j),i))%stat(1)=kin(cluster7_o(ncluster7_o(j),i))
        hist_uni(cluster7_o(ncluster7_o(j),i))%stat(2)=pot(cluster7_o(ncluster7_o(j),i))
        enddo 
                    do i=1,mcluster7_o(ncluster7_o(j),1)   !�����
                    cyp(cluster7_o(ncluster7_o(j),i))=1
                    enddo                
  endif     !��� ���������� � �������!!! 
  
  numjt=0
  numj_o=0
    knumjt=0
  knumj_o=0  
           ! if (mcluster(ncluster(j),2)==2) then   ! cluster tiene argon para a dentro
            !���������� ���-�� ������ � ���
                        if (ncluster(j) .ne. 0) then 
                    ip=1
                    do k=1,mcluster(ncluster(j),1)
                    if (typ(cluster(ncluster(j),k)) .ne.2) then
                    numjt(ip)=cluster(ncluster(j),k)         
                    ip=ip+1
                    endif
                    enddo   
                        else
                        numjt(1)=j
                        ip=2
                        endif
            knumjt=ip-1
!� numjt ����� ��� ����� �� cluser(ncluster(j)), ����� ������ ������
        if (ncluster_o(j) .ne. 0) then 
                    ip=1
                    do k=1,mcluster_o(ncluster_o(j),1)
                    if (typ(cluster_o(ncluster_o(j),k)) .ne.2) then
                    numj_o(ip)=cluster_o(ncluster_o(j),k)         
                    ip=ip+1
                    endif
                    enddo   
        else
                        numj_o(1)=j
                        ip=2
        endif
            knumj_o=ip-1 
!IN THIS LOCATION STOPPED YESTerDAY WE HAVE 
         if ( knumjt > knumj_o ) then   !��� �����������!!!

        ! aqui biene la prueba que este es nuevo encuentro y no regreso de prodigo pedazo

        do i=1,knumjt   !�����
        iexit=1
        do ik=1,knumj_o    !������
            if (numjt(i)>=numj_o(ik)) then
        if (numjt(i) .eq. numj_o(ik)) then 
        iexit=0 ! encontramos
        exit
        endif             
            endif
                 !numj_o(k) - ����� �� ������� �������� � j-�� ������
        ! �������
        enddo
      
                     if (iexit.eq.1) then
                     itwo=numjt(i)
                     exit 
                     endif
     
        enddo

           if (ncluster_o(itwo) .ne. 0) then 
                    ip=1
                    do k=1,mcluster_o(ncluster_o(itwo),1)
                    if (typ(cluster_o(ncluster_o(itwo),k)) .ne.2) then
                    numt_o(ip)=cluster_o(ncluster_o(itwo),k)         
                    ip=ip+1
                    endif
                    enddo   
           else
                        numt_o(1)=itwo
                        ip=2
           endif
            
        knumt_o=ip-1

!numj_all
 if (hist_uni(j).cluster(1) .ne. 0) then ! the biggest cluster with j
            if (hist_uni(j).cluster(1) .ne. 0) then ! the biggest cluster with j
                    ip=1
                    do k=1,nstat
                    
                    if ( hist_uni(j).cluster(k) .ne. 0  ) then
                    if (typ(hist_uni(j).cluster(k)) .ne.2) then
                    numj_all(ip)=hist_uni(j).cluster(k)         
                    ip=ip+1
                    endif
                    endif
                    enddo   
           else
                        numj_all(1)=j
                        ip=2
           endif
knumj_all=ip-1
             
      ik=0
         iexit=0
            do k=1,knumj_all
                do ik=1,knumt_o
            if (numt_o(ik)==numj_all(k)) then 
            iexit=iexit+1 ! encontramos
            exit
            endif
                enddo
            enddo
                 if (iexit.eq.0)  then   !��� ��� �� ���!
                 ik=2
                 else
                 ik=-4
                 endif
  else
ik=0    !j ������� �� ��� � ��������. ��� ����� �������
endif  
 if (hist_uni(itwo).cluster(1) .ne. 0) then  !the biggest cluster with itwo
            if (hist_uni(itwo).cluster(1) .ne. 0) then  !the biggest cluster with itwo
                    ip=1
                    do k=1,nstat
                    if ( hist_uni(itwo).cluster(k) .ne. 0 ) then
                    if (typ(hist_uni(itwo).cluster(k)) .ne.2) then
                    numt_all(ip)=hist_uni(itwo).cluster(k)         
                    ip=ip+1
                    endif
                    endif
                    enddo   
           else
                        numt_all(1)=j
                        ip=2
           endif
          knumt_all=ip-1
          iexit=0
            do k=1,knumt_all
                do i=1,knumj_o
            if (numj_o(i)==numt_all(k)) then 
            iexit=iexit+1 ! encontramos
            exit
            endif
                enddo
            enddo
            ip=0
                 if (iexit.eq.0)  then   !��� ��� �� ���!
                 ip=2
                 else
                 ip=-4
                 endif
   
else
ip=0    !j ������� �� ��� � ��������. ��� ����� �������
endif 


!!!!!!!
            if (ik==-4 .or. ip==-4) then
            ikk=1 !
                 write(41000,'(25(i5,1x))' ) ikk,nt,(hist_uni(j).cluster(ip),ip=1,12)

if( knumjt==knumj_all .or. knumjt==knumt_all) then   ! el cluster se recupero hasta el tama?o inicial
do i=1,knumjt   !�����
hist_dec(numjt(i)).cluster = 0
enddo
endif

                   do i=1,knumjt   !�����
                 hist_dec(numjt(i)).ntime = 0
                 enddo
            else    !if (ik==-4 .or. ip==-4) then
   
        !�������� �����
        
               if (ik==2) then
               ! escribimos la historia vieja de decay:
 
        if (hist_dec(j).ntime .ne. 0) then
    !��� ��������� ������� j �����������? �� � ������� ������������ ������
    i=1  ! ���� ���
    t=hist_dec(j)%ntime-hist_uni(j)%ntime

    write(44000,('(7(i7,1x),125(i5,1x))') ) i,knumj_all, nt, hist_dec(j)%ntime-hist_uni(j)%ntime, hist_uni(j)%narg, &
    hist_dec(j)%n1, hist_dec(j)%n2, (numj_all(i),i=1,12),(lstat(i),i=1,15)
    if (knumj_o.ne.1) then
 i=12  ! ����� �������� �������� �����������
 i1=0
    write(44000,('(7(i7,1x),125(i5,1x))') ) i,knumj_o, nt, nt-hist_dec(j)%ntime, i1, &
    i1, i1, (numj_all(i),i=1,12),(lstat(i),i=1,15)
    endif
    t= hist_uni(j)%ntime
     do ik=1,knumj_all    ! para evitar la doble marka
     if (hist_uni(numj_all(ik))%ntime==t) then
            hist_uni(numj_all(ik))%n1=0
            hist_uni(numj_all(ik))%n2=0
            hist_uni(numj_all(ik))%cluster=0 
            hist_uni(numj_all(ik))%ntime=hist_dec(j)%ntime
     endif
     enddo
       
        
        else  !if (hist_dec(j).ntime .ne. 0) then then ����������� ��������� ����������
            if (hist_uni(j)%ntime .ne. 0) then
          i=12  ! 
    write(44000,('(7(i7,1x),125(i5,1x))') ) i,knumj_all, nt, nt-hist_uni(j)%ntime, hist_uni(j)%narg, &
    hist_dec(j)%n1, hist_dec(j)%n2, (numj_all(i),i=1,12),(lstat(i),i=1,15)
            endif
        endif
        endif   !if (ik==2) then
     
        if(ip==2) then
        
   if (hist_dec(itwo).ntime .ne. 0) then
    !��� ��������� ������� itwo �����������? �� � ������� ������������ ������
    i=1  ! ���� ���
    t=hist_dec(itwo)%ntime-hist_uni(itwo)%ntime

    write(44000,('(7(i7,1x),125(i5,1x))') ) i,knumt_all,nt, hist_dec(itwo)%ntime-hist_uni(itwo)%ntime,  & 
    hist_uni(itwo)%narg, hist_dec(itwo)%n1, hist_dec(itwo)%n2,(numt_all(i),i=1,12),(lstat(i),i=1,15)
            if (knumt_o.ne.1) then
            i=12
            i1=0
            write(44000,('(7(i7,1x),125(i5,1x))') ) i,knumt_o, nt, nt-hist_dec(itwo)%ntime, i1, &
            i1, i1, (numt_all(i),i=1,12),(lstat(i),i=1,15)
            endif
      t= hist_uni(itwo)%ntime
            do ik=1,knumt_all    ! para evitar la doble marka
            if (hist_uni(numt_all(ik))%ntime==t) then
            hist_uni(numt_all(ik))%n1=0
            hist_uni(numt_all(ik))%n2=0
            hist_uni(numt_all(ik))%cluster=0 
            hist_uni(numt_all(ik))%ntime=hist_dec(itwo)%ntime
            endif
            enddo
       
        
    else  !if  (hist_dec(itwo).ntime .ne. 0) then����������� ��������� ����������
     if (hist_uni(itwo)%ntime .ne. 0) then
    i=12  ! 
    write(44000,('(7(i7,1x),125(i5,1x))') ) i,knumt_all,nt, nt-hist_uni(itwo)%ntime,  & 
    hist_uni(itwo)%narg, hist_dec(itwo)%n1, hist_dec(itwo)%n2,(numt_all(i),i=1,12),(lstat(i),i=1,15)
    endif   
    endif      
    endif   !        if(ip==2) then
            do ik=1,knumjt
            hist_uni(numjt(ik))%n1=knumj_o
            hist_uni(numjt(ik))%n2=knumt_o
            hist_uni(numjt(ik))%stat=0. 
            hist_uni(numjt(ik))%ntime=nt
            hist_uni(numjt(ik))%cluster=0
                do ip= 1,knumjt
                hist_uni(numjt(ik))%cluster(ip)=numjt(ip)
                enddo           
             hist_uni(numjt(ik))%narg=0
             
            hist_dec(numjt(ik))%cluster=0
            hist_dec(numjt(ik))%n1=0
            hist_dec(numjt(ik))%n2=0
            hist_dec(numjt(ik))%stat=0. 
            hist_dec(numjt(ik))%ntime=0
            
             hist_dec(numjt(ik))%narg=0
       
           ! cyp(numjt(ik))=1
            enddo

  endif   !if (ik==-4 .or. ip==-4) then
          
               do ik=1,mcluster(ncluster(j),1)
            cyp(cluster(ncluster(j),ik))=1
            enddo
    
            endif !from if (mcluster(ncluster(j),1) > mcluster_o(ncluster_o(j),1)) then   !��� �����������!!!
    !!!!!
           
    if ( knumjt < knumj_o ) then   !claster_o �������������?!!! ������� ���, �������, �������!
 
! ������� ��� �� ���� ��������� ! ��� ������ �� ���������!!!
  do i=1,knumj_o   !������ �������
        iexit=1
        do ik=1,knumjt    !����� ���������
            if (numjt(ik) .eq. numj_o(i)) then 
            iexit=0 ! encontramos
            exit
            endif             
        enddo
                     if (iexit.eq.1) then
                     itwo=numj_o(i)
                     exit 
                     endif
 enddo

           if (ncluster(itwo) .ne. 0) then 
                    ip=1
                    do k=1,mcluster(ncluster(itwo),1)
                    if (typ(cluster(ncluster(itwo),k)) .ne.2) then
                    numt_t(ip)=cluster(ncluster(itwo),k)         
                    ip=ip+1
                    endif
                    enddo   
           else
                        numt_t(1)=itwo
                        ip=2
           endif
           knumt_t=ip-1

 ! a �� �������� �� ���, ��� ����� ����� ����� �� ������?           
    if (hist_dec(j)%cluster(1) .ne. 0) then ! the biggest cluster from which j ��� come off
            !if (hist_dec(j).cluster(1) .ne. 0) then 
                    ip=1
                    do k=1,nstat
                    
                    if ( hist_dec(j)%cluster(k) .ne. 0  ) then
                    if (typ(hist_dec(j)%cluster(k)) .ne.2) then
                    numj_all(ip)=hist_dec(j)%cluster(k)         
                    ip=ip+1
                    endif
                    endif
                    enddo   
         
knumj_all=ip-1
        itres= hist_dec(j)%cluster(1)     ! ��, ��� �� ���� ���������� �����
      ikk=0
         iexit=0
            do k=1,knumj_all  ! from what come off earlier
                do ik=1,knumt_t
            if (numt_t(ik)==numj_all(k)) then 
            iexit=iexit+1 ! encontramos
            exit
            endif
                enddo
            enddo
                 if (iexit.eq.0)  then   !��� ��� �� ���!
                 ikk=2
                 else
                 ikk=-4
                 endif
    else
    ikk=0    !j ������� �� ��� � ������������� ��������. ��� ������ ������ despues de unirse
    endif  
          
             if (hist_dec(itwo)%cluster(1).ne.0) then

                    ik=1
                    do k=1,nstat
                    
                    if ( hist_dec(itwo).cluster(k) .ne. 0  ) then
                    if (typ(hist_dec(itwo).cluster(k)) .ne.2) then
                    numj_all(ik)=hist_dec(itwo).cluster(k)         
                    ik=ik+1
                    endif
                    endif
                    enddo   
 
knumj_all=ik-1
      !!!!!!!la correcci?n es necesaria!       
      ipp=0
         iexit=0
            do k=1,knumj_all  ! from what come off earlier
                do ip=1,knumt_t
            if (numjt(ip)==numj_all(k)) then 
            iexit=iexit+1 ! encontramos
            exit
            endif
                enddo
            enddo
                 if (iexit.eq.0)  then   !��� ��� �� ���!
                 ipp=2
                 else
                 ipp=-4
                 endif
  else
ipp=0    !itwo ������� �� ��� � ������������� ��������. ��� ������ ������
endif    
           
    !!!!!!!!!!!!!!!!!!!!!!!!!!       

        if (ikk==-4 .or. ipp==-4) then
            ikk=-1 !
                 write(41000,'(25(i5,1x))' ) ikk,nt,(hist_uni(j).cluster(ip),ip=1,15)
                   do i=1,knumj_o   !�����
                 hist_dec(numj_o(i))%ntime = nt
                 enddo
        else    !if (ik==-4 .or. ip==-4) then
        
         
       if( ikk == 2 .and. ipp==2) then    !������ �������
        
        if (hist_uni(j).cluster(1) .ne. 0) then ! the biggest cluster from which itwo ��� come off
          
                           ip=1
                    do k=1,nstat
                    
                    if ( hist_uni(j).cluster(k) .ne. 0  ) then
                    if (typ(hist_uni(j).cluster(k)) .ne.2) then
                    numj_all(ip)=hist_uni(j).cluster(k)         
                    ip=ip+1
                    endif
                    endif
                    enddo   
           else
                        numj_all(1)=j
                        ip=2
           endif
knumj_all=ip-1 
 i=-1 
 t=hist_dec(j)%ntime-hist_uni(j)%ntime
 if(t<0) then
      ik=1
 endif
 write(44000,('(7(i7,1x),125(i5,1x))'))  i,knumj_all,nt, hist_dec(j)%ntime-hist_uni(j)%ntime, &
 hist_uni(j)%narg, hist_dec(j)%n1, hist_dec(j)%n2,(numj_all(ik),ik=1,12),(lstat(i),i=1,15)

     
      ! para asegurar que la doble notacion no se ocurre  
       t= hist_uni(j)%ntime
   
    ! ����� ������ �����
       
         if (ncluster_o(itres) .ne. 0) then    !  �� ���������� ����, �.�. � ����� ������������ ���-������?
                    ip=1
                    do k=1,mcluster_o(ncluster_o(itres),1)
                    if (typ(cluster_o(ncluster_o(itres),k)) .ne.2) then
                    numt_o(ip)=cluster_o(ncluster_o(itres),k)         
                    ip=ip+1
                    endif
                    enddo   
           else
                        numt_o(1)=itres
                        ip=2
           endif
           knumt_o=ip-1
           
          if (hist_uni(itres)%ntime==t .and. knumt_o.ne.1) then
                                   
                do k=1,knumt_o
                hist_uni( numt_o(k) )%cluster=0
                do i=1,knumt_o
                hist_uni( numt_o(k) )%cluster(i)=numt_o(i)! �����!!!!!
                enddo
                enddo 
          endif    ! if (hist_uni(itwo)%ntime==t) then
  p=hist_dec(j)%ntime
      do k=1,knumj_all
            if (hist_uni( numj_all(k) )%ntime==t) then
            hist_uni( numj_all(k) )%ntime=p
            hist_uni( numj_all(k) )%n1=0
            hist_uni( numj_all(k) )%n2=0          
            hist_dec( numj_all(k) )%ntime=0
            hist_dec( numj_all(k) )%cluster=0  ! no estoy muy segura

            endif
       enddo   
              do k=1,knumj_o
         hist_uni( numj_o(k) )%cluster=0     
             do i=1,knumj_o
      hist_uni( numj_o(k) )%cluster(i)=numj_o(i)! �����!!!!!
      enddo
      enddo 
          
    if (knumj_o==1)  hist_uni( numj_o(1) )%ntime=0   ! si esta solito!
      endif   !if( ik == 2 .or. ip==2) then    !������ ����
      
            if (mcluster(ncluster(j),1).ne.0) then 
        do i=1,mcluster(ncluster(j),1)

        hist_dec( cluster(ncluster(j),i) )%ntime=nt
        hist_dec(cluster(ncluster(j),i))%n1=knumjt
        hist_dec(cluster(ncluster(j),i))%n2=knumt_t
           hist_dec(cluster(ncluster(j),i))%cluster=0
             do ip= 1,knumt_t
                hist_dec(cluster(ncluster(j),i))%cluster(ip)=numt_t(ip)
              enddo
  
            hist_dec(cluster(ncluster(j),i))%narg=0
        enddo
      else ! � ����� ������!
   
        hist_dec(j)%ntime=nt
        hist_dec(j)%n1=knumjt
        hist_dec(j)%n2=knumt_t
     hist_dec(j)%cluster=0
              do ip= 1,knumt_t
                hist_dec(j)%cluster(ip)=numt_t(ip)
              enddo         
        hist_dec(j)%narg=0
      endif   
      
       if (mcluster(ncluster(itwo),1).ne.0) then   
       do i=1,mcluster(ncluster(itwo),1)
            hist_dec(cluster(ncluster(itwo),i) )%ntime=nt
            hist_dec(cluster(ncluster(itwo),i))%n1=knumt_t
            hist_dec(cluster(ncluster(itwo),i))%n2=knumjt
                hist_dec(cluster(ncluster(itwo),i))%cluster=0
             do ip= 1,knumjt
                hist_dec(cluster(ncluster(itwo),i))%cluster(ip)=numjt(ip)
              enddo
            hist_dec(cluster(ncluster(itwo),i))%narg=0
       enddo   
        else
         hist_dec(itwo)%ntime=nt
         hist_dec(itwo)%n1=knumjt
         hist_dec(itwo)%n2=knumt_t
         hist_dec(itwo)%cluster=0
             do ip= 1,knumjt
                hist_dec(itwo)%cluster(ip)=numjt(ip)
              enddo        
         hist_dec(itwo)%narg=0
        endif 
                    endif !(if(ip<.-4 and ik<>-4)            endif     
           do i=1,knumj_o
                cyp(numj_o(i))=1
                enddo  
            do ik=1,mcluster_o(ncluster_o(j),1)
            cyp(cluster_o(ncluster_o(j),ik))=1
            enddo
            
        endif !from (.not.(mc    
!!!!!!!           
    
    endif   ! if (cyp(j).ne.1
    
        enddo   !form j=1,5000
 
!end of and now the ar-ar collisions will be counted
nt=nt+1
 !change it after the prove!

      enddo    !from do ifile=1,nmax para leerlos
												endif   !from if it is the first process							
   
! esta borrando  

  nttfile=nttfile+1
																
 call lammps_command (lmp, 'run 5000 pre no post no')
									    enddo    !from do while .true.

   CALL lammps_close(lmp);
    CALL mpi_finalize(ierr)
10 continue
end program 

