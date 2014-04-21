program simple

   use MPI
   use LAMMPS
   use, intrinsic :: ISO_C_binding, only : C_double, C_ptr, C_int

   implicit none

!-----------------------------------------------------------
interface one_step_data    ! interface  for geting info from lammps nodes
  subroutine one_step_data (data, natoms, lmp, me,lammps_fl,nprocs)
     use, intrinsic :: ISO_C_binding, only : C_double, C_ptr, C_int
     implicit none
     integer , dimension(:), allocatable,  intent(out) :: data
     type (C_ptr), intent(in) :: lmp
     integer, intent(in) :: me, lammps_fl,nprocs
     integer, intent(out) :: natoms
  end subroutine one_step_data     
end interface one_step_data
!----------------------------------------------------------

   integer :: ierr, me, nprocs    ! for MPI

   type (C_ptr) :: lmp		! for lammps communication
   integer :: lammps_fl, nprocs_main, nprocs_lammps, comm_lammps
   integer, dimension(:), allocatable :: buf
   integer :: buf_step, natoms, loop, nloop, mv


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!<<<<<<<<<<<<<<<< put initial block here >>>>>>>>>>>>>>>>>>
integer nstat,nclu,natom
parameter(nstat=125,natom=15001,nclu=2500)
  ! integer :: ierr, me, nprocs,natoms
  ! type (lammps_instance) :: lmp
  ! double precision :: compute, fix, fix2
  ! double precision, dimension(:), allocatable :: compute_v, mass, r,xx
  ! double precision, dimension(:,:), allocatable :: comp
 

logical exists
    real*4  statinum(nstat), para_esc(nstat),pott,kinn
    integer eof /0/,gist(nclu),typ(1:natom),cyp(natom),z,y,icolor,k,Code,cur_rec,ip,n,iexit,itwo,im,ntimej,istat,knum_sta 
    integer iras,ikk,ikp,ik,knumjt,knumj_o,knumj_all,knumt_o,knumt_all,np(12),nttfile
 !   character,allocatable:: struu(:)
    character*5 s1
    character*6 s2,c2
    character*50 sfile 
    
    integer*4 num_vecino(1:5),ncluster(1:natom),mcluster(0:nclu,2),ncluster_o(1:natom),mcluster_o(0:nclu,2)  ! 1- только металл, 2 - есть аргон
integer*8 ibeg(1:natom), n1, n2
    integer*4 i,j,num_atom,nt, M, cur_num ,cur_clu, numj_all(1000), numj_o(1000),numjt(1000),numt_o(1000),numt_all(1000)
    integer*4  nmax /10000/,natraso /5/
  
    real*8 pot(1 :natom), kin_all, kin(nclu),kin_clu(nclu)     !,ke1(1: natom) !Сашин массив
!kin(j) is the summ of all k.es, j is the number of atoms
!kin_clu(j) is k.e., j is el cluster's number
!character*17 sdir/'c:\1500_1\for_me\'/
    character*2250  cluster(0:nclu),cluster_o(0:nclu),clupro
        type hist
            integer*4 ntime  !время столкновения в н.у.
            integer*4 n1  !в состав кого входил атом J
            integer*4 n2 !с каким размером столкнулся
real*4 stat(nstat)
character*2250 cluster
        end type
    type (hist) hist_dec(15000), hist_uni(15000)
     

  nttfile=1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  ! setup MPI and various communicators
  ! driver runs on all procs in MPI_COMM_WORLD
  ! comm_lammps only has 1st P procs (could be all or any subset)

  CALL mpi_init(ierr)
  CALL mpi_comm_rank(MPI_COMM_WORLD,me,ierr);
  CALL mpi_comm_size(MPI_COMM_WORLD,nprocs,ierr);

  nprocs_main = 1
  nprocs_lammps = nprocs-nprocs_main

  IF ((nprocs_lammps < 1).or.(nprocs_main < 1)) THEN
     IF (me == 0) THEN
        PRINT *, 'ERROR: LAMMPS cannot use more procs than available'
        CALL mpi_abort(MPI_COMM_WORLD,2,ierr)
     END IF
  END IF

  lammps_fl = 0
  IF (me < nprocs_main) THEN
     lammps_fl = MPI_UNDEFINED
  ELSE
     lammps_fl = 1
  END IF

  CALL mpi_comm_split(MPI_COMM_WORLD,lammps_fl,0,comm_lammps,ierr)

  ! open LAMMPS input script on rank zero

  if (me == 0) print  *, 'RUN LAMMPS on ',nprocs_lammps,' processos'
  if (lammps_fl == 1)  then
    CALL lammps_open('',comm_lammps,lmp)
    call lammps_file (lmp, 'in.1')
  end if

!!-------  AE data-------
!<<<<<<<<<<<<<<<< put everythig befor main loop here >>>>>>>>>>>>>>>>>>
! open(15,File='t'//trim(sfile)//'.dat')      
! do ifile=1,nmax
!   open(1,File=sfile)    !'//trim(sfile)//'

!mcluster_o  = mcluster        !numero de atoms en cada cluster
!ncluster_o  = ncluster        !numero de cluster que mantiene este atom 
!cluster_o   = cluster           !numero del atom metallico
!-------------------------- 
! end of AE block   >>>>>>>>>>>>>>>>>>>>>>

!!!!!!!!! main loop begin here !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
nloop = 0  
nt=1
open(1,File='gist.dat') 
open(2,File='kinet.dat') 
 open(15,File='t1.dat') 
!if (me == 0) print  *, 'make ', nloop, ' steps'
        do while(.true.)
        
        if (int(nt/nmax) .eq. 1.*nt/(1.*nmax)) then
        nttfile=int(nt/nmax)+1
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
 close(15)
 open(15,File='t'//trim(sfile)//'.dat')  
 endif 
        do loop=1,nmax
!!! nt = loop*timestep for AE

  call one_step_data (buf, natoms, lmp, me,lammps_fl,nprocs)
!! end communication let's do separate work,
!!  meet at MPI_Reduce point (number of atoms)

  if (me .ne.  0)  cycle  ! Only the me=0 processor make statistics
!! dont foget to deallocale "buf, s_list" on the me=0 processor!!!

!!! take data from bufer and deallocate it
  buf_step=size(buf)/natoms
  mv = maxval(buf(1::buf_step)) ! max atom id

!<<<<<<<<<<<<<<<< put everythig here >>>>>>>>>>>>>>>>>>
!----  done reading cluster data from bufer ---------------
   eof=0
mcluster_o=0        !numero de atoms en cada cluster
ncluster_o=0        !numero de cluster que mantiene este atom 
cluster_o   =''           !numero del atom metallico
          !  do while (eof==0)
            
    mcluster_o= mcluster
    ncluster_o=ncluster
    cluster_o=cluster
 

10 format(i5) 
mcluster=0
ncluster=0
cluster=''
ikk=0
ikp=0
iras=0

     ncluster=0
      cluster=''
      cur_clu=1
      ik=0
          do j=1,natoms
      ! Read(1,*) num_atom,k,num_vecino(5)   !,kin(num_atom), pot(num_atom)
     num_atom=buf((j-1)*4+1)
     k=buf((j-1)*4+2)
     num_vecino(5)=buf((j-1)*4+3)
            if ((k.ne.2) .and. (k.ne.0)) then       !k=2  - es Argon, vamos de aqui
            !estoy leendo la infn sobre el atom de metall
            if (nt.eq.1) typ(num_atom)=1
            if (num_vecino(5).ne.num_atom) then
            
                    if (Ncluster(num_vecino(5))==0) then
                    Ncluster(num_vecino(5))=cur_clu
                    mcluster(Ncluster(num_vecino(5)),1)=1
                    if (typ(num_atom).ne. typ(num_vecino(5)) )  then 
                    mcluster(ncluster(num_vecino(5)),2)=2
                    else 
                    mcluster(Ncluster(num_vecino(5)),2)=1 
                    endif
                    cur_clu=cur_clu+1    !nuevo numero
                        if (len_trim(cluster(Ncluster(num_vecino(5)))).eq.0) then
                        write(s1,10)  num_vecino(5)
                        i=ncluster(num_vecino(5))    
                      cluster(i)=cluster(i)(:len_trim(cluster(i)))//' '//s1//' '   !(:len_trim(s1))
                        endif
                    endif
                    Ncluster(num_atom)= Ncluster(num_vecino(5))
                    i=ncluster(num_vecino(5))

                    mcluster(Ncluster(num_vecino(5)),1)=mcluster(Ncluster(num_vecino(5)),1)+1
                    if (typ(num_atom).ne. typ(num_vecino(5)) )  then 
                    mcluster(ncluster(num_vecino(5)),2)=2
                    else 
                    mcluster(Ncluster(num_vecino(5)),2)=1 
                    endif
                    write(s1,10)  num_atom
                    cluster(i)=cluster(i)(:len_trim(cluster(i)))//' '//s1//' '
     !       else
     !       mcluster(Ncluster(num_vecino(5)),1)=1 
     !       mcluster(Ncluster(num_vecino(5)),2)=1                     
            endif
       else 
            !from if k=2  - аргон, vamos de aqui por ahora
            if (nt.eq.1) typ(num_atom)=2
            if (num_vecino(5).ne.num_atom) then
                    if (Ncluster(num_vecino(5))==0) then
                    Ncluster(num_vecino(5))=cur_clu
                    mcluster(Ncluster(num_vecino(5)),1)=1
                   if (typ(num_atom).ne. typ(num_vecino(5)) )  then 
                    mcluster(ncluster(num_vecino(5)),2)=2
                    else 
                    mcluster(Ncluster(num_vecino(5)),2)=1 
                    endif
                    cur_clu=cur_clu+1    !nuevo numero
                        if (len_trim(cluster(Ncluster(num_vecino(5)))).eq.0) then
                        write(s1,10)  num_vecino(5)
                        i=ncluster(num_vecino(5))    
                         cluster(i)=cluster(i)(:len_trim(cluster(i)))//' '//s1//' '
                        endif
                    endif
                    Ncluster(num_atom)= Ncluster(num_vecino(5))
                    i=ncluster(num_vecino(5))

                   ! mcluster(Ncluster(num_vecino(5)),1)=mcluster(Ncluster(num_vecino(5)),1)+1
                    if (typ(num_atom).ne. typ(num_vecino(5)))  then
                    mcluster(Ncluster(num_vecino(5)),2)=2
                    else
                    mcluster(Ncluster(num_vecino(5)),2)=1
                    endif
                   ! write(s1,10)  num_atom
                   ! cluster(i)=cluster(i)(:len_trim(cluster(i)))//' '//s1//' '
            !no escribe numero de atom en cluster si es el argon!!!
            endif  
       endif
         
          !if (Code <> 0)
        enddo    !!j=1,m
                    do j=1,0
                   
                        if (mcluster(ncluster(j),1)==0 ) then   ! под ВОПРОСОМ!
                        mcluster(ncluster(j),1)=1
                        write(s1,10)  j
                        cluster(j)=' '//s1(:len_trim(s1))//' '
                        endif
                    enddo
                    cur_clu=cur_clu-1
        do i=1,cur_clu
        if (mcluster(i,2)==2) ik=ik+1
        enddo   
   
    if (nt==1) then
    mcluster_o= mcluster
    ncluster_o=ncluster
    cluster_o=cluster
    ibeg=0
  
     do j=1,natoms
     s2=cluster(ncluster(j))(1:6)
     if (mcluster(ncluster(j),1) > 1 .and. mcluster(ncluster(j),2)==1) then
        read(s2,*) k
                if (j==k) then
                write(s1,'(i5)') mcluster(ncluster(j),1)*6
                clupro= trim(s1)//' '//trim(cluster(ncluster(j)))
                endif
      hist_uni(j)%ntime=1 
      if (mcluster(ncluster(j),1)==2) then
      hist_uni(j)%n1=1
      hist_uni(j)%n2=1
      hist_uni(j)%stat=0.
      endif
      if (mcluster(ncluster(j),1)==3) then
      hist_uni(j)%n1=2
      hist_uni(j)%n2=1
      hist_uni(j)%stat=0.
      endif  
      if (mcluster(ncluster(j),1)==4) then
      hist_uni(j)%n1=2
      hist_uni(j)%n2=2
      hist_uni(j)%stat=0.
      endif
     hist_uni(j)%cluster=trim(clupro)
     ibeg(j)=1
     endif
     enddo
      
    else  !начальная конфигурация
    ! para haces el analisis ooo!
    
  
    cyp=0
        do j=1,natoms  ! todos los atoms

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
        if (typ(j).ne.2 .and.cyp(j).ne.1 ) then     !.and.cyp(j).ne.1
        !comparar las lineas

           ! if (mcluster(ncluster(j),2)==2) then   ! cluster tiene argon para a dentro
            !подсчитаем кол-во аргона в нем
                    clupro=''
                    numjt=0.
                        if (ncluster(j) .ne. 0) then 
                    clupro=cluster(ncluster(j))
                    z=1 
                    ip=1
                    do k=1,mcluster(ncluster(j),1)
                    y=z+5
                    s2=clupro(z:y)
                    z=y+1
                    
                    read(s2,*) numjt(ip)
                    if (typ(numjt(ip)) .ne.2) then
                    
                    ip=ip+1
                    endif
                    enddo   
                        else
                        numjt(1)=j
                        ip=2
                        endif
           ! endif   !from if (mcluster(ncluster(j),2)==2) then   ! cluster tiene argon para adentro
            knumjt=ip-1
!в numjt лежат все атомы из cluser(ncluster(j)), кроме атомов аргона

numj_o=0 
        if (ncluster_o(j) .ne. 0) then 
        clupro=cluster_o(ncluster_o(j))
        z=1 
        ip=1
        do k=1,mcluster_o(ncluster_o(j),1)
        y=z+5
        s2=clupro(z:y)
        z=y+1
        read(s2,*) numj_o(ip)
        if (typ(numj_o(ip)).ne.2)  then
        
        ip=ip+1
        endif
        enddo   ! now numj_o(1:N) -номера атомов, содержащихся в старом cluster_o, кроме атомов Ar
        else
        numj_o(1)=j
        ip=2
        endif
        knumj_o=ip-1
!if  (j==k) then  
        if ( knumjt < knumj_o ) then   !claster_o разваливается?!!! запишем его, запишем, запишем!
     
        s2=cluster_o(ncluster_o(j))(1:6)     !ik-это первый атом в viejo кластере
                    if (len_trim(s2) .ne. 0) then
                    read(s2,*) ik
                    endif

                    if(ik==j)  iras=iras+1 !это тупо число распадов
                    !и когда распался и на что
        z=1
                do ik=1,mcluster_o(ncluster_o(j),1)   !копаемся в старом, том, который распался
                y=z+5
                s2=cluster_o(ncluster_o(j))(z:y) !это из нового 
                z=y+1
                read(s2,*) n 

                 hist_dec(n)%ntime=nt    ! pos1 - это номер записи самого большого кластера, который был у атома j
            if (mcluster(ncluster(n),1)==0) then 
            hist_dec(n)%n2=1  
            else
            hist_dec(n)%n2=mcluster(ncluster(n),1)  !Какой стал
            endif
           
                hist_dec(n)%n1=mcluster_o(ncluster_o(n),1)        !Какой был
                enddo
        endif !from (.not.(mc
!!!!!!!           
         if ( knumjt > knumj_o ) then   !они встретились!!!
! а если это встреча атома металла с аргоном?
        ! aqui biene la prueba que este es nuevo encuentro y no regreso de prodigo pedazo
!!      read(12,*) strii !прочел запись атомов, составляющих самый большой кластер, в который входил этот атом.
!с кем же он столкнулся? найти атом, который не принадлежит _о
!write(14,'(4(i5,1x))')  nt,j,mcluster(ncluster(j),1),mcluster_o(ncluster_o(j),1)


        do i=1,knumjt   !новый
        iexit=1
        do ik=1,knumj_o    !старый
        if (numjt(i) .eq. numj_o(ik)) then 
        iexit=0 ! encontramos
        exit
        endif                                   !numj_o(k) - взято из старого кластера с j-ым атомом
        ! счетчик
        enddo
      
                     if (iexit.eq.1) then
                     itwo=numjt(i)
                     exit 
                     endif
     
        enddo
        if (ncluster_o(itwo) .ne. 0) then
        z=1
        ip=1
            do ik=1,mcluster_o(ncluster_o(itwo),1)   !копаемся в itwo перед столкновением и надо ли искать их все?
            y=z+5
            s2=cluster_o(ncluster_o(itwo))(z:y) !это из нового 
            z=y+1
             read(s2,*) numt_o(ip)    !это из itwo на предыдущем шаге. Его примеряют к numj_all
            if(typ(numt_o(ip)) .ne.2) then
          
            ip=ip+1
            endif
            enddo
        else
        numt_o(1)=itwo
        ip=2
        endif
        knumt_o=ip-1

                        if (ibeg(j).ne.0) then   ! иначе он никогда не был в кластере
                     
       !древняя запись наибольшего j, наблюдавшегося ранee
       !сoдержатся ли атомы ione в ней?

       !  read(12,POS=ibeg(j)) s2
       s2=hist_uni(j)%cluster(1:5)
          read(s2,*) ik 
        z=6
      
        knumj_all=ik/6
        ip=1
            do k=1,knumj_all
           ! read(12,POS=ibeg(j)+k*6) s2
           y=z+5 
             s2=hist_uni(j)%cluster(z:y)
               z=y+1
                
            read(s2,*) numj_all(ip)
          !  if (typ(numj_all(ip)).ne.2) then
             ip=ip+1
          !   endif
            enddo   ! now numj_all(1:N) -номера атомов, содержащихся в старом j из файла 
    knumj_all=ip-1
   

        
         iexit=0
            do k=1,knumj_all
                do ik=1,knumt_o
            if (numt_o(ik)==numj_all(k)) then 
            iexit=iexit+1 ! encontramos
            exit
            endif
                enddo
            enddo
  
  
  
           
                 if (iexit.eq.0)  then   !они там не ФСЕ!
                 ik=0
                 else
                 ik=-4
                 endif
   
else
ik=0    !j никогда не был в кластере. это новая встреча
endif 

if (ibeg(itwo).ne.0) then   ! иначе он никогда не был в кластере

       !древняя запись наибольшего itwo, наблюдавшегося ранee
       !сoдержатся ли атомы ione в ней?

       !  read(12,POS=ibeg(itwo)) s2
         s2=hist_uni(itwo)%cluster(1:6)
         read(s2,*) im 
        z=6
        
        knumt_all=im/6
        ip=1
            do k=1,knumt_all
             y=z+5
            !read(12,POS=ibeg(itwo)+k*6) s2
             s2=hist_uni(itwo)%cluster(z:y)
                z=y+1
               
            read(s2,*) numt_all(ip)
            !if (typ(numt_all(ip)).ne.2) then
            ip=ip+1
            !endif
            enddo   ! now numj_all(1:N) -номера атомов, содержащихся в старом j из файла 
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
  
           
                 if (iexit.eq.0)  then   !они там не ФСЕ!
                 ip=0
                 else
                 ip=-4
                 endif
   
else
ip=0    !j никогда не был в кластере. это новая встреча
endif 


            if (ik==-4 .or. ip==-4) then
            ikk=ikk+1 !
            endif
       
        !подводим итоги
        if (.not.(ik==-4 .or. ip==-4)) then  !новая запись не забудь ее сделать!!!
        ikp=ikp+1
        !numpro - все там которые были. calculemos estos atomes 
 !предыдущая жизнь и ее описание
 !подсчитаем, сколько разных записей в all задействовано при этом соударении: их 2, короче, и проверим их...

! despues de todo esto tenemos:
!- line_all(ik) - las notas en all.dat
!- num_all(ik) - uno de los atomos que estaban juntos muy lejos del momento
!- iline_all - quantas grupas задействовано
!do i=1,iline_all   j:
ntimej=hist_uni(j)%ntime
statinum=hist_uni(j)%stat
    if (hist_uni(j)%ntime.ne.0 ) then 
    ! если это=0, то описывать раннее поведение нельзя? он просто атом?
    !посчитаем сколько каждого размера. strii это то, от чего откололся кусок c j-ым атомом
    !mcluster_o(ncluster_o(numj_o(ik)),1)-сколько атомов теперь на предыдущем шаге, считающимся состоявшимся было у атома numj_o(ik)

  !          read(12,POS=ibeg(num_all(i))) s2
  !          read(s2,*) ik 
  !          knum=ik/6
  !          do k=1,knum
  !          read(12,POS=ibeg(num_all(i))+k*6) s2
  !          read(s2,*) numj_all(k)
  !          enddo   ! какой состав был у группы i в далеком прошлом 
  !  write(17,17) nt,i,(numpro_old(k),k=1,knum)
 17 format(i7,2x,i3,2x,35(i5,1x))
 
!numj_all
        do k=1,knumj_all 
        z=1 
        ip=1
        ! в те, что не участвуют в numjt - ничего не пишем здесь
        iexit=0
            do ik=1,knumj_o
            if (numj_all(k)==numj_o(ik)) then
            iexit=1
            exit
            endif
            enddo
                    if (iexit==1) then
                    knum_sta=knumj_o   ! в таком составе он был перед последним столкновением    
                    statinum(knum_sta)=statinum(knum_sta)+1./knum_sta
                  ! statinum(knumj_o)=statinum(knumj_o)+1.
                    hist_uni(numj_all(k))%ntime=nt
                  !  exit
                    endif
                          !  endif  !from  if ( cyp(numj_all(k) ==0) then
        enddo
  
        istat= 0
        do k=1,knumj_all   !y ahora elaboremos los otros atomes de la viejita
        if (ntimej .eq. hist_uni(numj_all(k))%ntime ) istat=1
        enddo
        
        
   if  (istat==0) then   ! в группе нет больше ничего, что могло бы столкнуться.    
    write(15,('(5(i7,1x),155(f4.1,1x))') ) knumj_all,nt,nt-ntimej, hist_uni(j)%n1, hist_uni(j)%n2,   &
    (statinum (ik), ik=1,nstat)
   else
        !запишем историю столкновений statinum в оставшиеся
        do k=1,knumj_all   !y ahora elaboremos los otros atomes de la viejita
        if (ntimej .eq. hist_uni(numj_all(k))%ntime )then
        hist_uni(numj_all(k))%stat=statinum
        endif
        enddo
   endif
 endif   ! if hist_uni&time=0 ! если это=0, то описывать раннее поведение нельзя

ntimej=hist_uni(itwo)%ntime
statinum=hist_uni(itwo)%stat
    if (hist_uni(itwo)%ntime.ne.0 ) then 

 
!numt_all
        do k=1,knumt_all 
        z=1 
        ip=1
        ! в те, что не участвуют в numtt - ничего не пишем здесь
        iexit=0
            do ik=1,knumjt
            if (numt_all(k)==numjt(ik)) then
            iexit=1
            exit
            endif
            enddo
                    if (iexit==1) then
                   knum_sta=knumt_o   ! в таком составе он был перед последним столкновением    
                   statinum(knum_sta)=statinum(knum_sta)+1./knum_sta
                  ! statinum(knumt_o)=statinum(knumt_o)+1.
                    hist_uni(numt_all(k))%ntime=nt
                   
                    endif
                          !  endif  !from  if ( cyp(numt_all(k) ==0) then
        enddo
  
        istat= 0
        do k=1,knumt_all   !y ahora elaboremos los otros atomes de la viejita
        if (ntimej .eq. hist_uni(numt_all(k))%ntime ) istat=1
        enddo
        
        
   if  (istat==0) then   ! в группе нет больше ничего, что могло бы столкнуться.    
    write(15,('(i4,1x,i12,1x,3(i7,1x),155(f4.1,1x),4(i7,1x))') ) knumt_all,nt,nt-ntimej, hist_uni(itwo)%n1, hist_uni(itwo)%n2,   &
    (statinum (ik), ik=1,nstat)
   else
        !запишем историю столкновений statinum в оставшиеся
        do k=1,knumt_all   !y ahora elaboremos los otros atomes de la viejita
        if (ntimej .eq. hist_uni(numt_all(k))%ntime )then
        hist_uni(numt_all(k))%stat=statinum
        endif
        enddo
   endif
 endif   ! if hist_uni&time=0 ! если это=0, то описывать раннее поведение нельзя


 !escribimos este bueno nuevo grande cluster 
 !   Записывается только один раз, когда первый атом, остальные только отмечают себе адрес записи в 12
    
             !   ibeg(j)=cur_rec !номер строки в 12, где записаны самые толстые кластеры
             
    
                write(s1,'(i5)') knumjt*6
               clupro= s1
             
                    do ik=1,knumjt
                    write(s1,'(i5)')  numjt(ik)
                    clupro= (trim(clupro)//' '//s1//' ')
                    enddo
            do ik=1,knumjt
         
            hist_uni(numjt(ik))%n1=knumj_o
            hist_uni(numjt(ik))%n2=knumt_o
            hist_uni(numjt(ik))%stat=0. 
            hist_uni(numjt(ik))%ntime=nt  
            hist_uni(numjt(ik))%cluster=trim(clupro)
            ibeg(numjt(ik))=1
cyp(numjt(ik))=1
            enddo
            ke1=0.
            !the cluster kinetic energy calc. (en cada atom)
            do ik=1,knumjt
            kin_clu(ncluster(ik))= kin_clu(ncluster(ik))+ke1(numjt(ik))   !/knumjt
            enddo
            
  endif   !откуда он? ik<>-4 
            endif !from if (mcluster(ncluster(j),1) > mcluster_o(ncluster_o(j),1)) then   !они встретились!!!
    endif   !from если это не аргон 
    
        enddo   !form j=1,m
  
    endif   !if nt==1

        gist=0
        do j=ik,cur_clu
        gist(mcluster(ik,1))=gist(mcluster(ik,1))+1
        kin(mcluster(ik,1))=kin(mcluster(ik,1))+kin_clu(ik)
        enddo
        kin_all=0.
        do j=1,natoms
        kin_all=kin_all+ke1(j)
        enddo
        write(1,('(i12,1x,f12.4,1x,155(i5,1x))') ) nt, kin_all,  (gist(i), i=1,nstat)
        write(2,('(i12,1x,f12.4,1x,155(i5,1x))') ) nt,kin_all,(kin(i), i=1,nstat)
nt=nt+1 

!close(1)
          
!close(1,status='delete')
                            enddo    !from do ifile=1,nmax para leerlos
! end of AE block   >>>>>>>>>>>>>>>>>>>>>>
  print  *, '<<<<<<<<<<< step ', loop, ' done >>>>>>>>>>>>>>>>>'
enddo	! end of main loop "do loop=1,nloop" (true) ahora

!-----------------------------------------------------
  if (allocated(buf)) deallocate(buf)
  IF (lammps_fl == 1) CALL lammps_close(lmp);
  ! close down MPI
  CALL mpi_finalize(ierr)

end program simple

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine one_step_data (buf, natoms, lmp, me,lammps_fl,nprocs)
   use MPI
   use LAMMPS
   use, intrinsic :: ISO_C_binding, only : C_double, C_ptr, C_int
   implicit none

    integer :: status(MPI_STATUS_SIZE)
      integer , dimension(:), allocatable, intent(out) :: buf
      type (C_ptr), intent(in) :: lmp
      integer, intent(in) :: me, lammps_fl,nprocs
      integer, intent(out) :: natoms
   integer :: i,j,iproc,ierr
 !  integer, dimension(:), allocatable,target :: buf
   real (C_double), dimension(:,:), pointer :: comp => NULL()
   integer (C_int), dimension(:), pointer :: tag => NULL()
   integer (C_int), dimension(:), pointer :: type => NULL()

! every processor has number of atoms. me=0 takes total number of atoms
  natoms = 0
  if (lammps_fl == 1)  then
    call lammps_command(lmp,'run 5')
    call lammps_extract_atom (tag, lmp, 'id')
    call lammps_extract_atom (type, lmp, 'type')
    call lammps_extract_compute (comp, lmp, 'clu', 1, 2)
    natoms = size(tag)
  endif
  call MPI_Reduce(natoms, i, 1, MPI_INT,MPI_SUM, 0, MPI_COMM_WORLD ,ierr)
  if (me == 0)  natoms = i

!  print *, 'done',me, natoms
!  call MPI_BARRIER(MPI_COMM_WORLD, ierr)

! send all information to me=0 node
  if (allocated(buf)) deallocate(buf)
  allocate (buf(4*natoms))

  if (me .ne. 0)  then
    j=1
    do i=1, natoms 
!      print *, tag(i),type(i),comp(:,i)
      buf(j)= tag(i)
      buf(j+1)= type(i)
      buf(j+2)= idnint(comp(1,i))
      buf(j+3)= idnint(comp(2,i))
      j=j+4
    enddo

!!send inforation to prosess with rank 0
    call MPI_Send(buf, 4*natoms, MPI_INT, 0, me, MPI_COMM_WORLD, ierr)
    deallocate(buf)
  else		! I am rank 0 processor!!!
!! receive information from other nodes in one bufer
    j=1
    do iproc = 1,nprocs-1 
      call MPI_Recv(buf(j:), 4*natoms, MPI_INT,  iproc, iproc, MPI_COMM_WORLD, status, ierr)
      call MPI_Get_count( status,  MPI_INT, i, ierr )
!      print *, me, iproc, i,j,(j+3),(j+i+3)
      j=j+i
    enddo
  endif
contains
end subroutine one_step_data

