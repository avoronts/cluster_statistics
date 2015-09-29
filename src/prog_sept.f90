!!!!!!!!!!!!my programm !!!!!!!!!!!!!
program simple

  use LAMMPS_gether
  implicit none

integer nstat,natom,nclu,n_entre
!parameter(nstat=125,natom=25001,nclu=2500)
parameter(nstat=625,natom=15001,nclu=15500)   !this data for exchange with LAMMPS program

    !integer me - from LAMMPS_gether module
    !integer natos - from LAMMPS_gether module
    !integer int_buf - from LAMMPS_gether module
    !integer dbl_buf - from LAMMPS_gether module
    integer  typ(1:natom), num_atom, num_vecino4(1:natom), num_vecino5(1:natom),ncu
    double precision pot(1:natom), kin(1:natom), vx(1:natom),vy(1:natom),vz(1:natom), & 
                     x(1:natom),y(1:natom),z(1:natom)

!this piece of text is devoted to statistical relations betw...
    integer*4 statist(300), npromediar/300/,delta/5/, delta_tau, nstatist/10/
    real*8 mpro,dpro,a,en_cu
    integer curnum1 /2/,lstat(nstat)
        logical :: sexist

!
    character*200 sfile 
    character*150 sline

! this data for cluster storage
    integer*4 ncluster(1:natom),ncluster_o(1:natom),mcluster(0:nclu,2), & 
              mcluster_o(0:nclu,2),lcluster(0:nclu), pos_cluster  ! 1- cколько всех , 2 - сколько аргон
    integer cluster(1:nclu,1:nstat),cluster_o(1:nclu,1:nstat)


    integer*4 numj_all(1000), numt_t(1000),numj_o(1000),numjt(1000),numt_o(1000),numt_all(1000)
    integer   ipos_cluster, ithree,ipt,jpos_cluster, knumj_all, knumt_t, knumj_o, knumjt, knumt_o, knumt_all
    
    integer cyp(natom)

    integer*4 i, j, jj, i1, k, nt, cur_clu , cur_clu7
    integer t, p, n, itwo,ipp, itres, ityp, ikk,ip,iexit,ifile,ik,n1,n2,n3
        
    integer*4  nmax /10000/,natraso /5/,sostav(500)
    integer cluster7_o(1:nclu,1:nstat),cluster7(1:nclu,1:nstat),mcluster7(0:nclu,2),ncluster7(1:natom),mcluster7_o(0:nclu,2),ncluster7_o(1:natom)

        type hist
            integer*4 ntime  !время столкновения в н.у. ! или время разваливания
            integer*4 n1  !в состав кого входил атом J   ! на что распался
            integer*4 n2 !с каким размером столкнулся
            integer*4 narg ! сколько столкновений с аргоном было, начиная с ntime
	    real*4 stat(nstat)   ! это ек задействованно в распаде
	    integer cluster(nstat) !В какой кл,астер теперь входит
        end type
    type (hist) hist_dec(25000), hist_uni(25000)
     
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

!--------------------------------------------------------
CALL init_mpi_lammps('in.1') !!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  
nt=0
mcluster_o(:,:) = 0        !numero de atoms en cada cluster
ncluster_o(:) = 0        !numero de cluster que mantiene este atom 
cluster_o(:,:) = 0          !numero del atom metallico

mcluster7_o(:,:) = 0        !numero de atoms en cada cluster
ncluster7_o(:) = 0        !numero de cluster que mantiene este atom 
cluster7_o(:,:) = 0          !numero del atom metallico



open(150,file='stat2.dat')
open(41000,file='s1231.dat')
   open(44000,File='r1.dat')      
   open(42000,File='en_cu1.dat')
   open(43000,File='collis1.dat')      
do i1=1,1000000
    call one_step_data(5)
          enddo

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
nt =nt+1
i=0
if(nt.eq.1) then
do ik=1,natoms
if (typ(ik) .eq.1) i=i+1
enddo
ncu=i
endif



!write(41000,('(i12,1x,2(e20.13,1x))')) nt,pot(12),kin(12)
!write(41000,*) nt,pot(12),kin(12)
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Here we read 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if(int(1.*nt/nmax).eq. 1.*nt/nmax) then
   close(44000)
close(42000)
close(43000) 
   write(sfile,'(i0)') int(nt/nmax)+1
   !open(41000,File='osc'//trim(sfile)//'.dat') 
   !open(15,File='t'//trim(sfile)//'.dat')      
   open(44000,File='r'//trim(sfile)//'.dat')      
   open(42000,File='en_cu'//trim(sfile)//'.dat')
   open(43000,File='collis'//trim(sfile)//'.dat')      
   ! open(47000,File='prot'//trim(sfile)//'.dat')      
endif



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!                do ifile=1,nmax
!										    if (me.eq.1) then  
 mcluster7_o(:,:) = mcluster7(:,:)
 ncluster7_o(:) = ncluster7(:)
 cluster7_o(:,:) = cluster7(:,:)

 mcluster7(:,:)=0
 ncluster7(:)=0
 cluster7(:,:)=0
 
  mcluster_o(:,:) = mcluster(:,:)
 ncluster_o(:) = ncluster(:)
 cluster_o(:,:) = cluster(:,:)

 mcluster(:,:)=0
 ncluster(:)=0
 cluster(:,:)=0
 
 
 cur_clu=1
 cur_clu7=1
 ! el chiquito
 do num_atom=1,natoms
    ipos_cluster = Ncluster(num_vecino5(num_atom))
   
   if (num_vecino5(num_atom).ne.num_atom ) then   !quitamos the clusteres of argon solo .and. num_vecino5(num_atom)<5000 
            cyp(num_atom)=10 
             if (ipos_cluster == 0) then
                    Ncluster(num_vecino5(num_atom))=cur_clu
                    ipos_cluster = cur_clu
                   
                    mcluster(ipos_cluster,1)=1
                     if (typ(num_vecino5(num_atom)).ne. 1 )  then 
                    mcluster(ncluster(num_vecino5(num_atom)),2)=1
                    endif
                    !if (typ(num_atom).ne. typ(num_vecino5(num_atom)) )  then 
                    !mcluster(ncluster(num_vecino5(num_atom)),2)=2
                    !endif
                    cur_clu=cur_clu+1    !nuevo numero
                    if (cluster(ipos_cluster,1) .eq. 0) cluster(ipos_cluster,1) = num_vecino5(num_atom)
             endif
             Ncluster(num_atom) = ipos_cluster 

             !if (typ(num_atom) .ne.2) lcluster(ipos_cluster) = lcluster(ipos_cluster)+1
               if (typ(num_atom).ne. 1 )  then 
                    mcluster(ncluster(num_atom),2)=mcluster(ncluster(num_atom),2)+1
                    endif
             if (mcluster(ipos_cluster,1)<nstat) then
                mcluster(ipos_cluster,1) = mcluster(ipos_cluster,1)+1
                !if (typ(num_atom).ne. typ(num_vecino5(num_atom)) )  then 
                !mcluster(ncluster(num_vecino5(num_atom)),2)=2
                !else
                !mcluster(Ncluster(num_vecino5(num_atom)),2)=1 
                !endif
                cluster(ipos_cluster,mcluster(ipos_cluster,1)) = num_atom
            endif
    endif

 
   jpos_cluster= Ncluster7(num_vecino4(num_atom))
   
   if (num_vecino4(num_atom).ne.num_atom ) then   !quitamos the cluster7es of argon solo .and. num_vecino4(num_atom)<5000 
         !  cyp(num_atom)=10 
             if (jpos_cluster== 0) then
                    Ncluster7(num_vecino4(num_atom))=cur_clu7
                    jpos_cluster= cur_clu7
                   
                    mcluster7(jpos_cluster,1)=1
                     if (typ(num_vecino4(num_atom)).ne. 1 )  then 
                    mcluster7(ncluster7(num_vecino4(num_atom)),2)=1
                    endif
                    !if (typ(num_atom).ne. typ(num_vecino4(num_atom)) )  then 
                    !mcluster7(ncluster7(num_vecino4(num_atom)),2)=2
                    !endif
                    cur_clu7=cur_clu7+1    !nuevo numero
                    if (cluster7(jpos_cluster,1) .eq. 0) cluster7(jpos_cluster,1) = num_vecino4(num_atom)
             endif
             Ncluster7(num_atom) = jpos_cluster

             !if (typ(num_atom) .ne.2) lcluster7(jpos_cluster) = lcluster7(jpos_cluster)+1
               if (typ(num_atom).ne. 1 )  then 
                    mcluster7(ncluster7(num_atom),2)=mcluster7(ncluster7(num_atom),2)+1
                    endif
             if (mcluster7(jpos_cluster,1)<nstat) then
                mcluster7(jpos_cluster,1) = mcluster7(jpos_cluster,1)+1
                !if (typ(num_atom).ne. typ(num_vecino4(num_atom)) )  then 
                !mcluster7(ncluster7(num_vecino4(num_atom)),2)=2
                !else
                !mcluster7(Ncluster7(num_vecino4(num_atom)),2)=1 
                !endif
                cluster7(jpos_cluster,mcluster7(jpos_cluster,1)) = num_atom
            endif
    endif
 enddo    !!num_atom==1,natoms
 if (nt.eq.1) then
   mcluster_o(:,:) = mcluster(:,:)
 ncluster_o(:) = ncluster(:)
 cluster_o(:,:) = cluster(:,:)
    mcluster7_o(:,:) = mcluster7(:,:)
 ncluster7_o(:) = ncluster7(:)
 cluster7_o(:,:) = cluster7(:,:)
 endif
          lstat=0
    do i=1,cur_clu-1   ! это набор кластеров на данный момент
    lstat(mcluster(i,1)-mcluster(i,2))=lstat(mcluster(i,1)-mcluster(i,2))+1
    enddo 
    lstat(1)=0
        do i=2,natoms !
            if (typ(i)==1) then
                if(ncluster(i)==0) then
                lstat(1)=lstat(1)+1
                endif
            endif
        enddo    

do i=1,0
if (cyp(i) .eq. 10)  then
 if (i.eq.4916 .or. i.eq.2535 .or. i.eq.758.or. i.eq.7134.or. i.eq.3140.or. i.eq.1654 .or.& 
 i.eq.790.or. i.eq.693.or. i.eq.29.or. i.eq.651.or. i.eq.46.or. i.eq.3886) then

write(i,'(i10,1x,120(i10,1x,8(F10.5,1x)))') nt, (cluster(ncluster(i),j),x(cluster(ncluster(i),j)),y(cluster(ncluster(i),j)), &
z(cluster(ncluster(i),j)),vx(cluster(ncluster(i),j)),vy(cluster(ncluster(i),j)), &
vz(cluster(ncluster(i),j)),kin(cluster(ncluster(i),j)),pot(cluster(ncluster(i),j)),j=1,mcluster(ncluster(i),1))
do j=1,mcluster(ncluster(i),1)
 cyp(cluster(ncluster(i),j))=0
 enddo
endif
 if(i.eq.552.or. i.eq.995.or. i.eq.2522.or. i.eq.730 .or. i.eq.4314.or. i.eq.2480.or. i.eq.398.or. i.eq.2907.or. i.eq.415.or. i.eq.2653) then
 write(i,'(i10,1x,120(i10,1x,8(F10.5,1x)))') nt, (cluster(ncluster(i),j),x(cluster(ncluster(i),j)),y(cluster(ncluster(i),j)), &
z(cluster(ncluster(i),j)),vx(cluster(ncluster(i),j)),vy(cluster(ncluster(i),j)), &
vz(cluster(ncluster(i),j)),kin(cluster(ncluster(i),j)),pot(cluster(ncluster(i),j)),j=1,mcluster(ncluster(i),1))
do j=1,mcluster(ncluster(i),1)
 cyp(cluster(ncluster(i),j))=0
 enddo
endif


 if(i.eq.92.or. i.eq.737.or. i.eq.1213.or. i.eq.1325 .or. i.eq.553.or. i.eq.1081.or. i.eq.1604) then
 write(i,'(i10,1x,120(i10,1x,8(F10.5,1x)))') nt, (cluster(ncluster(i),j),x(cluster(ncluster(i),j)),y(cluster(ncluster(i),j)), &
z(cluster(ncluster(i),j)),vx(cluster(ncluster(i),j)),vy(cluster(ncluster(i),j)), &
vz(cluster(ncluster(i),j)),kin(cluster(ncluster(i),j)),pot(cluster(ncluster(i),j)),j=1,mcluster(ncluster(i),1))
do j=1,mcluster(ncluster(i),1)
 cyp(cluster(ncluster(i),j))=0
 enddo
endif

endif

enddo

!!!!!!!!!!!!!!!!!!!!! check2. my code 21.05.2015 !!!!!!!!!!!!!!!!1
!  write(sfile,'(i0)') (nt+1)*natraso
!  sfile=trim('dump/check_out'//trim(sfile)//'.bin')
!  open(5050,File=sfile,form='unformatted')    !'//trim(sfile)//'
!
!  write(5050) Ncluster		!    integer*4 ncluster(1:natom),mcluster(0:nclu,2),lcluster(0:nclu)
!  write(5050) mcluster
!  write(5050) lcluster
!  write(5050) cluster		! integer cluster(1:nclu,1:nstat)
!  close(5050)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11

cyp(:)=0
do j=1,natoms  ! todos los atoms

   if( typ(j).eq.1) then
  if (cyp(j).ne.1 ) then     !.and.cyp(j).ne.1

      !!!!!!!el pedazo de energia!!!!start
       numjt(:)=0
    knumjt=0
    
        numj_o(:)=0
    knumj_o=0
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
!в numjt лежат все атомы из cluser(ncluster7(j)), кроме атомов аргона
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
            
         !  busca el numero del atomo solitario
            
            
            !!!!!!они встретились с Cu!!!\!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       if ( knumjt - knumj_o==1 .and. mcluster7(ncluster7(j),2) == mcluster7_o(ncluster7_o(j),2) ) then   
       !найти энергию прибывшего атома
         do i=1,knumjt   !новый
        iexit=1
        do ik=1,knumj_o    !старый
         !   if (numjt(i)>=numj_o(ik)) then
        if (numjt(i) .eq. numj_o(ik)) then 
        iexit=0 ! encontramos
        exit
        endif             
          !  endif
                 !numj_o(k) - взято из старого кластера с j-ым атомом
        ! счетчик
        enddo
      
                     if (iexit.eq.1) then
                     itwo=numjt(i)   !это его номеро
                     exit 
                     endif
     
        enddo
         do i=1,mcluster7(ncluster7(j),1)   !новый
      hist_uni(cluster7(ncluster7(j),i))%stat(1)=kin(itwo)
      hist_uni(cluster7(ncluster7(j),i))%stat(2)=pot(itwo)
         cyp( cluster7(ncluster7(j),i))=1  
        enddo   

        endif     !они встретились  с Cuprumom !!!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !они расстались с аргоном!!!\!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

    if ( knumjt-knumj_o==-1 .and. knumjt.ne.1 .and. mcluster7(ncluster7(j),2)== mcluster7_o(ncluster7_o(j),2) .and. hist_uni(j)%stat(1) .ne. 0. ) then   !они расстались с аргоном!!!\

 !найти номер Cu, его кинет энергию after. !найти энергию убывшего атома
         do i=1,knumj_o  !старый gordo  NO PASO PRUEBA 10.07
        iexit=1
        do ik=1,knumjt    !новый flaco
            !if (numjt(ik)>=numj_o(i)) then
        if (numjt(ik) .eq. numj_o(i)) then 
        iexit=0 ! encontramos
        exit
        endif             
            !endif
                 !numj_o(k) - взято из старого кластер а с j-ым атомом
        ! счетчик
        enddo
      
                     if (iexit.eq.1) then
                     itwo=numj_o(i)
                     exit 
                     endif
     
        enddo

if (hist_uni(j)%stat(1).ne.0) then
   en_cu = kin(itwo)-hist_uni(j)%stat(1)
                          
  if ( abs(en_cu) .gt. 1.e-7)   write(42000,('(4(i10,1x),2(f14.9,1x),115(i7,1x))') ) nt,cluster7_o(ncluster7_o(j),1), knumjt,knumj_o, &
  en_cu, hist_uni(j)%stat(2),(lstat(ip),ip=1,50)
endif
        do i=1,mcluster7_o(ncluster7_o(j),1)   
        hist_uni(cluster7_o(ncluster7_o(j),i))%stat(1)=0.
        hist_uni(cluster7_o(ncluster7_o(j),i))%stat(2)=0.
         cyp(cluster7_o(ncluster7_o(j),i))=1
        enddo 
                               
  endif     !они расстались с Cu!!! 
  endif   !if typ==1
  endif   !if cyp=1
  enddo   !from do j=1,natoms  ! todos los atoms
  
        !comparar las lineas
cyp(:)=0
!<<<<<<<<<<<<<< present cluster >>>>>>>>>>>>>>> 
    numjt(:)=0
    knumjt=0
     do j=1,natoms  ! todos los atoms II
    
     
        if( typ(j).eq.1) then
  if (cyp(j).ne.1 ) then     !.and.cyp(j).ne.1
 numjt(:)=0
    knumjt=0
     
    ! if (mcluster(ncluster(j),2)==2) then   ! cluster tiene argon para a dentro
    !подсчитаем кол-во аргона в нем
    if (ncluster(j) .ne. 0) then 
        ip=1
        do k=1,mcluster(ncluster(j),1)
          if (typ(cluster(ncluster(j),k)) .ne.2) then
            numjt(ip)=cluster(ncluster(j),k)  
                if (numjt(ip).eq.10282) then
                i1=0
                endif        
            ip=ip+1
          endif
        enddo   
    else
        numjt(1)=j
        
        ip=2
    endif
    knumjt=ip-1
      if (nt==1 .and. knumjt.ne.1) then
     
      
         do ik=1,knumjt    ! para evitar la doble marka
                   hist_uni(numjt(ik))%ntime=nt
                     hist_uni(numjt(ik))%n1=knumjt-1
                     hist_uni(numjt(ik))%n2=1
                      hist_uni(numjt(ik))%cluster=0
                     do i=1,knumjt
                     hist_uni(numjt(ik))%cluster(i)=numjt(i) 
                     enddo
               cyp(numjt(ik))=1 
                 enddo
      
      
      
      
      
      endif ! from nt==1      
!<<<<<<<<<<<<<< old cluster >>>>>>>>>>>>>>>
    numj_o(:)=0
    knumj_o=0 
    !в numjt лежат все атомы из cluser(ncluster(j)), кроме атомов аргона
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
    if (nt==24.or.nt==848.or. nt==896 .or. nt==390 .or. nt==638) then
     i1=0
    endif
 

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!<<<<<<<<<<<<<<<<< checking>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!IN THIS LOCATION STOPPED YESTerDAY WE HAVE 
    if ( knumjt > knumj_o ) then   !они встретились!!!
   
   if (j.eq.1711 .and. nt>45000) then
                i1=0
                endif  
        ! aqui biene la prueba que este es nuevo encuentro y no regreso de prodigo pedazo

        do i=1,knumjt   !новый 
          iexit=1
          do ik=1,knumj_o    !старый
            if (numjt(i)>=numj_o(ik)) then
                if (numjt(i) .eq. numj_o(ik)) then 
                   iexit=0 ! encontramos
                   exit
                endif
            endif
          !numj_o(k) - взято из старого кластера с j-ым атомом
          ! счетчик
          enddo

          if (iexit.eq.1) then
             itwo=numjt(i)
             exit 
          endif
        enddo


        numt_o(:)=0
        knumt_o=0
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

!^^^^^^^^^^^^^^^^^^^^^ What is it?????????? 
        !numj_all
        numj_all(:)=0
        knumj_all=0
        if (hist_uni(j)%cluster(1) .ne. 0) then        ! the biggest cluster with j
            if (hist_uni(j)%cluster(1) .ne. 0) then  ! the biggest cluster with j
                ip=1
                do k=1,nstat
                    
                    if ( hist_uni(j)%cluster(k) .ne. 0  ) then
                       if (typ(hist_uni(j)%cluster(k)) .ne.2) then
                          numj_all(ip)=hist_uni(j)%cluster(k)         
                         ip=ip+1
                       endif
                    else
                    exit
                    endif
                enddo   
            else
                numj_all(1)=j
                ip=2
            endif
            knumj_all=ip-1
                 if (knumj_all==1) then
               i1=0
               endif
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
             if (iexit.eq.0)  then   !они там не ФСЕ!
                 ik=2
             else
                 ik=-4
             endif
        else
          ik=0    !j никогда не был в кластере. это новая встреча
        endif 


!^^^^^^^^^^^^^^^^^^^^^ What is it?????????? 
        numt_all(:)=0
        knumt_all=0 
        if (hist_uni(itwo)%cluster(1) .ne. 0) then  !the biggest cluster with itwo
            if (hist_uni(itwo)%cluster(1) .ne. 0) then  !the biggest cluster with itwo
                ip=1
                do k=1,nstat
                  if ( hist_uni(itwo)%cluster(k) .ne. 0 ) then
                    if (typ(hist_uni(itwo)%cluster(k)) .ne.2) then
                      numt_all(ip)=hist_uni(itwo)%cluster(k)         
                      ip=ip+1
                    endif
                  else
                    exit
                  endif
                enddo   
           else
                numt_all(1)=j
                ip=2
           endif
           knumt_all=ip-1
               if (knumt_all==1) then
               i1=0
               endif
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
            if (iexit.eq.0)  then   !они там не ФСЕ!
                 ip=2
             else
                 ip=-4
             endif
   
        else
        ip=0    !j никогда не был в кластере. это новая встреча
        endif 



!<<<<<<<<<<<<<<<< start save grow >>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        if (ik==-4 .or. ip==-4) then
            ikk=1 !
!                 write(41000,'(25(i5,1x))' ) ikk,nt,(hist_uni(j)%cluster(ip),ip=1,12)

           ! if( knumjt==knumj_all .or. knumjt==knumt_all) then   ! el cluster se recupero hasta el tama?o inicial
               do i=1,knumjt   !новый
                  hist_dec(numjt(i))%cluster = 0
               enddo
           ! endif

            do i=1,knumjt   !новый
                 hist_dec(numjt(i))%ntime = 0
            enddo
        else    !if (ik==-4 .or. ip==-4) then
   write(43000,('(3(i12,1x),125(i5,1x))') ) knumjt, knumj_o, nt, numjt(1),(lstat(i),i=1,50)   
        !подводим итоги
!!!!!!!!!!!!!!!!!!!!!
            if (ik==2) then
               ! escribimos la historia vieja de decay:
 
              if (hist_dec(j)%cluster(1) .ne. 0) then
                 !это результат развала j соединяется? да и запишем состоявшийся развал
                 i=1  ! сдох сам
                 t=hist_dec(j)%ntime-hist_uni(j)%ntime
                 if(t<0)then
                   i1=1
                 endif
                 write(44000,('(7(i12,1x),125(i5,1x))') ) i,knumj_all, nt, hist_dec(j)%ntime-hist_uni(j)%ntime, & 
                              numj_all(1), hist_dec(j)%n1, hist_dec(j)%n2, (lstat(i),i=1,50)
     if (knumj_o.ne.1) then
       i=13  ! жизнь продукта прервана соединением, но он был развалившимся до этого
       i1=0
       write(44000,('(7(i12,1x),125(i5,1x))') ) i,knumj_o, nt, & 
                   nt-hist_dec(j)%ntime, numj_o(1), i1, i1, (lstat(i),i=1,50)
     endif
          ! y esto es para todos
                 t= hist_uni(j)%ntime
                 !!!!!
                      ithree=hist_dec(j)%cluster(1)  !кто-нибудь отсоединившийся.hist_dec(j)%ntime .ne. 0)!!!
               ! развал закончился.запишем остатки
                
                      numt_t(:)=0
        knumt_t=0
        if (ncluster(ithree) .ne. 0) then 
            ip=1
            do k=1,mcluster(ncluster(ithree),1)
                if (typ(cluster(ncluster(ithree),k)) .ne.2) then
                    numt_t(ip)=cluster(ncluster(ithree),k)         
                    ip=ip+1
                endif
            enddo   
        else
            numt_t(1)=ithree
            ip=2
        endif
        knumt_t=ip-1

             ! but only if nothing passed before!!!     
             if (knumt_t.ne.1) then
                  t= hist_uni(j)%ntime
                  do ik=1,knumt_t    ! para evitar la doble marka
                     if (hist_uni(numt_t(ik))%ntime==t) then
                     hist_uni(numt_t(ik))%ntime=0
                        hist_uni(numt_t(ik))%n1=0
                        hist_uni(numt_t(ik))%n2=0
                       hist_uni(numt_t(ik))%cluster=0 !это случалось здесь!!!
                       do i=1,knumt_t
                        hist_uni(numt_t(ik))%cluster(i)=numt_t(i)
                       enddo
                        hist_uni(numt_t(ik))%ntime=hist_dec(itwo)%ntime
                     endif
                  enddo
        t= hist_dec(j)%ntime
                  do ik=1,knumt_t    ! para evitar la doble marka
                     if (hist_dec(numt_t(ik))%ntime==t) then
                     hist_dec(numt_t(ik))%ntime=0
                        hist_dec(numt_t(ik))%n1=0
                        hist_dec(numt_t(ik))%n2=0
                       hist_dec(numt_t(ik))%cluster=0 
                     endif
                  enddo
                 
                    else  !!from(knumt_t.ne.1) then
                    t= hist_uni(j)%ntime
                   ik=1
                     if (hist_uni(numt_t(ik))%ntime==t) then
                     hist_uni(numt_t(ik))%ntime=0
                        hist_uni(numt_t(ik))%n1=0
                        hist_uni(numt_t(ik))%n2=0
                       hist_uni(numt_t(ik))%cluster=0 !это случалось здесь!!!
                        hist_uni(numt_t(ik))%ntime=hist_dec(itwo)%ntime
                     endif
               
        t= hist_dec(j)%ntime
                 ik=1
                     if (hist_dec(numt_t(ik))%ntime==t) then
                     hist_dec(numt_t(ik))%ntime=0
                        hist_dec(numt_t(ik))%n1=0
                        hist_dec(numt_t(ik))%n2=0
                       hist_dec(numt_t(ik))%cluster=0 
                     endif
               
                 
                 endif   !!from(knumt_t.ne.1) then
            
                 do ik=1,knumj_all    ! para evitar la doble marka
                   if (hist_uni(numj_all(ik))%ntime==t) then
                     hist_uni(numj_all(ik))%n1=0
                     hist_uni(numj_all(ik))%n2=0
                    ! hist_uni(numj_all(ik))%cluster=0 
                     hist_uni(numj_all(ik))%ntime=hist_dec(j)%ntime
                   endif
                 enddo
               else  !if (hist_dec(j)%ntime .ne. 0) then then соединяется результат соединения
                 if (hist_uni(j)%ntime .ne. 0) then 
                    i=12  ! 
                    if (knumj_all==2) then
                    write(10,('(125(i8,1x))')) i,j,nt,numj_all(1)
                    endif
                    write(44000,('(7(i12,1x),125(i5,1x))') ) i,knumj_all, nt, & 
                    nt-hist_uni(j)%ntime, numj_all(1), hist_dec(j)%n1, hist_dec(j)%n2, (lstat(i),i=1,50)
                 endif
               endif
            endif   !if (ik==2) then
!!!!!!!!!!!!!!!!!!!!!!!!!!1
            if(ip==2) then
        
               if (hist_dec(itwo)%cluster(1) .ne. 0) then
                  !это результат развала itwo соединяется? да и запишем состоявшийся развал
                  i=1  ! сдох сам
              
                  t=hist_dec(itwo)%ntime-hist_uni(itwo)%ntime
                  if(t<0)then
                     i1=1
                  endif
                  write(44000,('(7(i12,1x),125(i5,1x))') ) i,knumt_all,nt, hist_dec(itwo)%ntime-hist_uni(itwo)%ntime,  &	
                         numt_all(1), hist_dec(itwo)%n1, hist_dec(itwo)%n2,(lstat(i),i=1,50)
                  if (knumt_o.ne.1) then   
                     ! но он был развалившимся до этого
                     i=13
                     i1=0
                     write(44000,('(7(i12,1x),125(i5,1x))') ) i,knumt_o, nt, nt-hist_dec(itwo)%ntime, itwo, & 
                           i1, i1, (lstat(i),i=1,50)
                  endif
               ithree=hist_dec(itwo)%cluster(1)  !кто-нибудь отсоединившийся.
               ! развал закончился.запишем остатки
                  t= hist_uni(itwo)%ntime
                      numt_t(:)=0
        knumt_t=0
        if (ncluster(ithree) .ne. 0) then 
            ip=1
            do k=1,mcluster(ncluster(ithree),1)
                if (typ(cluster(ncluster(ithree),k)) .ne.2) then
                    numt_t(ip)=cluster(ncluster(ithree),k)         
                    ip=ip+1
                endif
            enddo   
        else
            numt_t(1)=ithree
            ip=2
        endif
        knumt_t=ip-1
! тот, кто отделился ранее и ни с кем не встретился
             ! but only if nothing passed before!!!
             if (knumt_t .ne.1) then     
                   t= hist_uni(itwo)%ntime
                  do ik=1,knumt_t    ! para evitar la doble marka
                     if (hist_uni(numt_t(ik))%ntime==t) then
                        hist_uni(numt_t(ik))%n1=0
                        hist_uni(numt_t(ik))%n2=0
                       hist_uni(numt_t(ik))%cluster=0 !это случалось здесь!!!
                       do i=1,knumt_t
                        hist_uni(numt_t(ik))%cluster(i)=numt_t(i)
                       enddo
                        hist_uni(numt_t(ik))%ntime=hist_dec(itwo)%ntime
                     endif
                  enddo
       t= hist_dec(itwo)%ntime
                  do ik=1,knumt_t    ! para evitar la doble marka
                     if (hist_dec(numt_t(ik))%ntime==t) then
                        hist_dec(numt_t(ik))%n1=0
                        hist_dec(numt_t(ik))%n2=0
                       hist_dec(numt_t(ik))%cluster=0 !это случалось здесь!!!
                      endif
                  enddo
                  
                  else   !from if (knumt_t .ne.1) then 
                        t= hist_uni(itwo)%ntime
                  ik=1
                     if (hist_uni(numt_t(ik))%ntime==t) then
                        hist_uni(numt_t(ik))%n1=0
                        hist_uni(numt_t(ik))%n2=0
                       hist_uni(numt_t(ik))%cluster=0 !это случалось здесь!!!
                       hist_uni(numt_t(ik))%ntime=0
                     endif
               
       t= hist_dec(itwo)%ntime
                  ik=1
                     if (hist_dec(numt_t(ik))%ntime==t) then
                        hist_dec(numt_t(ik))%n1=0
                        hist_dec(numt_t(ik))%n2=0
                       hist_dec(numt_t(ik))%cluster=0 
                         hist_dec(numt_t(ik))%ntime=0
                      endif
                                  
                  endif     !from if (knumt_t .ne.1) then 
               else  !if  (hist_dec(itwo)%ntime .ne. 0) thenсоединяется результат соединения
                  if (hist_uni(itwo)%ntime .ne. 0) then
                     i=12  ! 
                     if (knumj_all==2) then
                    write(10,('(125(i8,1x))')) i+2,j,nt,numj_all(1)
                    endif
                     write(44000,('(7(i12,1x),125(i5,1x))') ) i,knumt_all,nt, nt-hist_uni(itwo)%ntime,  & 
                              itwo, hist_uni(itwo)%n1, hist_uni(itwo)%n2,(lstat(i),i=1,50)
                  endif   
               endif      
            endif   !        if(ip==2) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
          if (hist_dec(j)%cluster(1)==0  .and. hist_dec(j)%ntime.ne.0) then
          ip=0
          write(10,('(12(i10,1x))') )   nt,ip,j,hist_dec(j)%ntime
          endif
             if (hist_uni(j)%cluster(1)==0  .and. hist_uni(j)%ntime.ne.0) then
             ip=1
          write(10,('(12(i10,1x))') )   nt,ip,j,hist_uni(j)%ntime
          endif
       do ik=1,mcluster(ncluster(j),1)
         cyp(cluster(ncluster(j),ik))=1
       enddo
    
    endif ! if ( knumjt > knumj_o ) then    !они встретились!!!
!<<<<<<<<<<<<<  end save grow >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!<<<<<<<<<<<<<<< start save decay >>>>>>>>>>>>>>>>>>>>>>>>>>
    if ( knumjt < knumj_o ) then   !claster_o разваливается?!!! запишем его, запишем, запишем!
      
   if (j.eq.1711 .and. nt>45000) then
                i1=0
                endif   
       ! выясним кто от него отделился ! нет защиты от колебаний!!!
       do i=1,knumj_o   !старый больной
          iexit=1
          do ik=1,knumjt    !новый маленький
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


       numt_t(:)=0
       knumt_t=0
       if (ncluster(itwo) .ne. 0) then 
          ip=1
          do k=1,mcluster(ncluster(itwo),1)
            if (typ(cluster(ncluster(itwo),k)) .ne.2) then
                numt_t(ip)=cluster(ncluster(itwo),k)         
                ! тот, кто оторвался сейчас, в настоящем состоянии
                ip=ip+1
            endif
          enddo   
       else
          numt_t(1)=itwo
          ip=2
       endif
       knumt_t=ip-1


       numj_all(:)=0   !??????????????
       knumj_all = 0   !?????????????
       ! a не случится ли так, что отрыв будет похож на старый?           
       if (hist_dec(j)%cluster(1) .ne. 0) then ! the biggest cluster from which j цше come off
          !if (hist_dec(j)%cluster(1) .ne. 0) then 
          ip=1
          do k=1,nstat
             if ( hist_dec(j)%cluster(k) .ne. 0  ) then
                if (typ(hist_dec(j)%cluster(k)) .ne.2) then
                    numj_all(ip)=hist_dec(j)%cluster(k)    ! и кто уже удалялся от j      
                    ip=ip+1
                endif
             else
                exit
             endif
          enddo   
          knumj_all=ip-1
         
          itres= hist_dec(j)%cluster(1)     ! то, что от него отвалилось тогда
          ikk=0
          iexit=0
          do k=1,knumj_all  ! from what come off earlier
             do ik=1,knumt_t   ! сравним тех, кто уже, с теми кто сейчас
                if (numt_t(ik)==numj_all(k)) then 
                   iexit=iexit+1 ! encontramos
                   exit
                endif
             enddo
          enddo
          if (iexit.eq.0)  then   !они там не ФСЕ!
             ikk=2
          else
             ikk=-4
          endif
       else
          ikk=0    !j никогда не был в распадающемся кластере. это первый распад despues de unirse
       endif  

       numj_all(:)=0   !??????????????
       knumj_all = 0   !?????????????
      
       if (hist_dec(itwo)%cluster(1).ne.0) then
          ik=1
          do k=1,nstat
             if ( hist_dec(itwo)%cluster(k) .ne. 0  ) then
                 if (typ(hist_dec(itwo)%cluster(k)) .ne.2) then
                    numj_all(ik)=hist_dec(itwo)%cluster(k)         
                    ik=ik+1
                 endif
             else
                 exit
             endif
          enddo   
          knumj_all=ik-1

          !!!!!!!la correcciоn es necesaria!       
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
          if (iexit.eq.0)  then   !они там не ФСЕ!
             ipp=2
          else
             ipp=-4
          endif
       else
          ipp=0    !itwo никогда не был в распадающемся кластере. это первый распад
       endif    
           
    !!!!!!!!!!!!!!!!!!!!!!!!!!       

       if (ikk==-4 .or. ipp==-4) then
            ikk=-1 !
            !  write(41000,'(25(i5,1x))' ) ikk,nt,(hist_uni(j)%cluster(ip),ip=1,15)
            do i=1,knumj_o   !новый
               hist_dec(numj_o(i))%ntime = nt
                hist_dec(itwo)%n1=knumjt
                hist_dec(itwo)%n2=knumt_t
            enddo
       else    !if (ik==-4 .or. ip==-4) then
            !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
         
            if (ikk == 2 .and. ipp==2) then    !распад распада
               numj_all(:)=0   !??????????????
               knumj_all = 0   !?????????????
               if (hist_uni(j)%cluster(1) .ne. 0) then ! the biggest cluster from which itwo цше come off
                  ip=1
                  do k=1,nstat
                    if ( hist_uni(j)%cluster(k) .ne. 0  ) then
                      if (typ(hist_uni(j)%cluster(k)) .ne.2) then
                        numj_all(ip)=hist_uni(j)%cluster(k)         
                        ip=ip+1
                      endif
                    else
                      exit
                    endif
                  enddo   
               else
                  numj_all(1)=j
                  ip=2
               endif
               knumj_all=ip-1 
               if (knumj_all==1) then
               i1=0
               endif
               i=1 
               t=hist_dec(j)%ntime-hist_uni(j)%ntime
               write(44000,('(7(i12,1x),125(i5,1x))'))  i,knumj_all,nt, hist_dec(j) %ntime-hist_uni(j)%ntime, & 
                          numj_all(1), hist_dec(j)%n1, hist_dec(j)%n2,(lstat(ik),ik=1,50)
               ! para asegurar que la doble notacion no se ocurre  
               t= hist_uni(j)%ntime
   
               numt_o(:)=0   !??????????????
               knumt_o = 0   !?????????????
               ! найди второй кусок
               if (ncluster_o(itres) .ne. 0) then    !  на предыдущем шаге, т.л. а вдруг одновременно что-нибудь?
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
               
               
               ipt=hist_dec(j)%ntime
               ! второй свободен как птица. узаконим это.
               do k=1,knumt_o
                  if (hist_uni(numt_o(k))%ntime==t) then
                     hist_uni(numt_o(k))%ntime=0
                     hist_uni( numt_o(k) )%cluster=0
                     hist_uni( numt_o(k) )%n1=0
                     hist_uni( numt_o(k) )%n2=0  
                     hist_dec( numt_o(k) )%ntime=0 
                     hist_dec( numt_o(k) )%cluster=0 
                      hist_dec( numt_o(k) )%n1=0
                     hist_dec( numt_o(k) )%n2=0 
                     if (knumt_o.ne.1) then
                        do i=1,knumt_o
                           hist_uni( numt_o(k) )%cluster(i)=numt_o(i)! 
                        enddo
                        hist_uni(numt_o(k))%ntime=ipt
                     endif
                  endif    
               enddo 
               
               !el primero es libre tambien
               do k=1,knumj_o
                  if (hist_uni( numj_o(k) )%ntime==t) then
                      hist_uni( numj_o(k) )%ntime=ipt
                      hist_uni( numj_o(k) )%n1=0
                      hist_uni( numj_o(k) )%n2=0          
                      hist_dec( numj_o(k) )%ntime=0
                      hist_dec( numj_o(k) )%cluster=0 
                      hist_uni( numj_o(k) )%cluster=0  ! no estoy muy segura
                      do i=1,knumj_o
                         hist_uni( numj_o(k) )%cluster(i)=numj_o(i)! 
                      enddo
                  endif
               enddo 

            endif   !if( ik == 2 .or. ip==2) then    !распад расп
            !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
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
            else ! я снова одинок!
   
                 hist_dec(j)%ntime=nt
                 hist_dec(j)%n1=knumjt
                 hist_dec(j)%n2=knumt_t
                 hist_dec(j)%cluster=0
                 do ip= 1,knumt_t
                    hist_dec(j)%cluster(ip)=numt_t(ip)
                 enddo         
                 hist_dec(j)%narg=0
            endif   
            !<<<<<<<<<<<<<<<<<<<<<<<y porque aqui se hace un error?
            ! он унес с собой аргон и я в этот аргон тоже пишу. оно мне не нада!!! потерпи!
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
            !<<<<<<<<<<<<<<<<<<<<<<<<<
        endif !(if(ip<.-4 and ik<>-4)            endif     
      if (hist_dec(j)%cluster(1)==0  .and. hist_dec(j)%ntime.ne.0) then
          ip=3
          write(10,('(12(i10,1x))') )   nt,ip,j,hist_dec(j)%ntime
          endif
             if (hist_uni(j)%cluster(1)==0  .and. hist_uni(j)%ntime.ne.0) then
             ip=4
          write(10,('(12(i10,1x))') )   nt,ip,j,hist_uni(j)%ntime
          endif
        do i=1,knumj_o
            cyp(numj_o(i))=1
        enddo  
     
            
    endif !if ( knumjt < knumj_o ) then   !claster_o разваливается?!!!   
!<<<<<<<<<<<<<<< end saving decay >>>>>>>>>>>>>>>>>>>>>>>>>>

  endif   ! if (cyp(j).ne.1
  endif    !from typ.eq.1
enddo   !form j=1,natoms


!      <<<<<< end of main circle
!  .nd now the ar-ar collisions will be counted
!      enddo    !from do ifile=1,nmax para leerlos
enddo    !from do while .true.

call finalize_mpi_lammps()
10 continue
end program 
 