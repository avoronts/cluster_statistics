!!!!!!!!!!!!my programm !!!!!!!!!!!!!
program simple

  use LAMMPS_gether
   !implicit none

integer nstat,nclu,natom
parameter(nstat=125,natom=25001,nclu=2500,n_entre=15)

!this piece of text is devoted to statistical relations betw...
integer*4 statist(300), npromediar /300/,delta/5/, delta_tau, nstatist /10/
real*8 mpro,dpro,a
integer curnum1 /2/
!endik
logical :: sexist

    real*4  statinum(nstat), para_esc(nstat),pott,kinn,kin_entre(n_entre),l_entre(n_entre) ,  x(natom),y(natom),z(natom)
    integer eof /0/,typ(1:natom),cyp(natom),i1,i2,i3,k,Code,cur_rec ,iras,np(12),lstat(nstat)
    character,allocatable:: struu(:)
    character*5 s1
    character*6 s2,c2
    character*200 sfile 
    character*150 sline
    integer*4 num_vecino4(1:natom),num_vecino5(1:natom),ncluster(1:natom),mcluster(0:nclu,2),ncluster_o(1:natom),mcluster_o(0:nclu,2),lcluster(0:nclu)  ! 1- только металл, 2 - есть аргон
integer*8 n1, n2, iarg
    integer*4 i,j,jj,num_atom,nt, M, cur_num ,cur_clu,cur_clu7, numj_all(1000), numt_t(1000),numj_o(1000),numjt(1000),numt_o(1000),numt_all(1000)
    integer*4  nmax /10000/,natraso /5/,sostav(500)
  
    real*8 pot(1 :natom),kin(1: natom),l2, vx(1:natom),vy(1:natom),vz(1:natom)
    integer cluster(1:nclu,1:nstat),cluster_o(1:nclu,1:nstat)
    !integer cluster7_o(1:nclu,1:nstat),cluster7(1:nclu,1:nstat),mcluster7(0:nclu,2),ncluster7(1:natom),mcluster7_o(0:nclu,2),ncluster7_o(1:natom)
        type hist
            integer*4 ntime  !время столкновения в н.у. ! или время разваливания
            integer*4 n1  !в состав кого входил атом J   ! на что распался
            integer*4 n2 !с каким размером столкнулся
            integer*4 narg ! сколько столкновений с аргоном было, начиная с ntime
real*4 stat(nstat)   ! это ек задействованно в распаде
integer cluster(nstat) !В какой кластер теперь входит
        end type
    type (hist) hist_dec(25000), hist_uni(25000)
     
mcluster_o=0        !numero de atoms en cada cluster
ncluster_o=0        !numero de cluster que mantiene este atom 
cluster_o =0          !numero del atom metallico

  nttfile=1

!--------------------------------------------------------
  CALL init_mpi_lammps('in.1') !!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  
  
nt=0
!!me=1 ! DON'T forget to borrar it!
  ! do something
 	  do i=1,5000
          typ(i)=1
       !   write(sfile,'(i0)') i
       !   open(i,file=trim(sfile)//'.dat')
          enddo    
      do i=5001,natom
      typ(i)=2
      enddo
 ikk=0
ikp=0
open(150,file='stat2.dat')
      
										do while(.true.)
! entro y esta pensando
!										if (me.eq.1) then 

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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if(int(1.*nt/nmax).eq. 1.*nt/nmax) then
close(44000)

         
if((int(1.*(nt)/nstatist/nmax).eq. 1.*(nt)/nstatist/nmax) .and. nt .ne. 0) then
!if((int(1.*(nt-1)/nstatist/nmax).eq. 1.*(nt-1)/nstatist/nmax) .and. nt .ne. 1) then
statist=0 
cyp=0
a=0.
 j=int(1.*nt/nmax)
 do   ifile=j-nstatist + 1, j, 1

 write(sfile,'(i0)') ifile
 sfile=('r'//trim(sfile)//'.dat')
   inquire (file=sfile,exist=sexist)
   if(sexist) then
    open(18,File=sfile,err=12) 

        do while(.true.)    

    read(18,('(2(i12,1x),a12,1x,a12,1x,3(i12,1x),125(i5,1x))'),end=11 ) ityp,n,sline,sline, n3, n1,n2, (numjt(ik), ik=1,12)
    if(sline.ne.'*******') then
    read(sline,'(i12)') delta_tau
    endif

          if   (curnum1==n .and. (ityp==1 .or. ityp==-1) .and. sline.ne.'*******') then    !ityp==1 .or. ityp==-1)
             a=a+1.
               do i=1,npromediar
               if (delta_tau>=(i-1.)*delta .and. delta_tau < i*delta) then
               statist(i)=statist(i)+1
               exit
               endif 
               enddo 
   
          endif

        enddo
         11 continue


         
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


 write(sfile,'(i0)') nt/nmax+1
 !open(41000,File='osc'//trim(sfile)//'.dat') 
 !open(15,File='t'//trim(sfile)//'.dat')      
 open(44000,File='r'//trim(sfile)//'.dat')      
!open(42000,File='pro_argon'//trim(sfile)//'.dat')
! open(43000,File='pro_argon_menos'//trim(sfile)//'.dat')      
! open(47000,File='prot'//trim(sfile)//'.dat')      
endif

!                do ifile=1,nmax
!										    if (me.eq.1) then  
    mcluster_o= mcluster
    ncluster_o=ncluster
    cluster_o=cluster
      
mcluster=0
mcluster7=0
	ncluster=0
      cluster=0
      
      cur_clu=1

 ! el chiquito
 do num_atom=1,natoms
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
                  !  cyp(num_vecino5(num_atom))=10
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
           
!write(112,('(40(i7,1x))')) (lstat(i),i=1,15)
 !do i=1,natom
!if (cyp(i) .ne. 10)   cyp(i)=0
!enddo

lstat(1)=0.
do i=1,0
if (cyp(i) .eq. 10)  then
 if (i.eq.255 .or. i.eq.8 .or. i.eq.1651) then
write(i,'(i10,1x,120(i10,1x,8(F10.5,1x)))') nt, (cluster(ncluster(i),j),x(cluster(ncluster(i),j)),y(cluster(ncluster(i),j)), &
z(cluster(ncluster(i),j)),vx(cluster(ncluster(i),j)),vy(cluster(ncluster(i),j)), &
vz(cluster(ncluster(i),j)),kin(cluster(ncluster(i),j)),pot(cluster(ncluster(i),j)),j=1,mcluster(ncluster(i),1))
do j=1,mcluster(ncluster(i),1)
 cyp(cluster(ncluster(i),j))=0
 enddo
endif
endif
if (ncluster(i).eq.0) then
 lstat(1)=lstat(1)+1.
 endif
enddo


cyp=0
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
  
 
  numjt=0
  numj_o=0
    knumjt=0
  knumj_o=0  
           ! if (mcluster(ncluster(j),2)==2) then   ! cluster tiene argon para a dentro
            !подсчитаем кол-во аргона в нем
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
!IN THIS LOCATION STOPPED YESTerDAY WE HAVE 
         if ( knumjt > knumj_o ) then   !они встретились!!!

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
                    else
                    exit
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
                 if (iexit.eq.0)  then   !они там не ФСЕ!
                 ik=2
                 else
                 ik=-4
                 endif
  else
ik=0    !j никогда не был в кластере. это новая встреча
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
                    else
                    exit
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
                 if (iexit.eq.0)  then   !они там не ФСЕ!
                 ip=2
                 else
                 ip=-4
                 endif
   
else
ip=0    !j никогда не был в кластере. это новая встреча
endif 



!!!!!!!
            if (ik==-4 .or. ip==-4) then
            ikk=1 !
!                 write(41000,'(25(i5,1x))' ) ikk,nt,(hist_uni(j).cluster(ip),ip=1,12)

if( knumjt==knumj_all .or. knumjt==knumt_all) then   ! el cluster se recupero hasta el tama?o inicial
do i=1,knumjt   !новый
hist_dec(numjt(i)).cluster = 0
enddo
endif

                   do i=1,knumjt   !новый
                 hist_dec(numjt(i)).ntime = 0
                 enddo
            else    !if (ik==-4 .or. ip==-4) then
   
        !подводим итоги
        
               if (ik==2) then
               ! escribimos la historia vieja de decay:
 
        if (hist_dec(j).ntime .ne. 0) then
    !это результат развала j соединяется? да и запишем состоявшийся развал
    i=1  ! сдох сам
    t=hist_dec(j)%ntime-hist_uni(j)%ntime
    if(t<0)then
    i1=1
    endif
    write(44000,('(7(i12,1x),125(i5,1x))') ) i,knumj_all, nt, hist_dec(j)%ntime-hist_uni(j)%ntime, hist_uni(j)%narg, &
    hist_dec(j)%n1, hist_dec(j)%n2, (numj_all(i),i=1,20)
    if (knumj_o.ne.1) then
 i=13  ! жизнь продукта прервана соединением, но он был развалившимся до этого
 i1=0
    write(44000,('(7(i12,1x),125(i5,1x))') ) i,knumj_o, nt, nt-hist_dec(j)%ntime, i1, &
    i1, i1, (numj_o(i),i=1,20)
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
       
        
        else  !if (hist_dec(j).ntime .ne. 0) then then соединяется результат соединения
            if (hist_uni(j)%ntime .ne. 0) then 
          i=12  ! 
    write(44000,('(7(i12,1x),125(i5,1x))') ) i,knumj_all, nt, nt-hist_uni(j)%ntime, hist_uni(j)%narg, &
    hist_dec(j)%n1, hist_dec(j)%n2, (numj_all(i),i=1,20)
            endif
        endif
        endif   !if (ik==2) then
     
        if(ip==2) then
        
   if (hist_dec(itwo).ntime .ne. 0) then
    !это результат развала itwo соединяется? да и запишем состоявшийся развал
    i=1  ! сдох сам
    t=hist_dec(itwo)%ntime-hist_uni(itwo)%ntime
     if(t<0)then
    i1=1
    endif
    write(44000,('(7(i12,1x),125(i5,1x))') ) i,knumt_all,nt, hist_dec(itwo)%ntime-hist_uni(itwo)%ntime,  & 
    hist_uni(itwo)%narg, hist_dec(itwo)%n1, hist_dec(itwo)%n2,(numt_all(i),i=1,20)
            if (knumt_o.ne.1) then   
            ! но он был развалившимся до этого
            i=13
            i1=0
            write(44000,('(7(i12,1x),125(i5,1x))') ) i,knumt_o, nt, nt-hist_dec(itwo)%ntime, i1, &
            i1, i1, (numt_o(i),i=1,20)
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
       
        
    else  !if  (hist_dec(itwo).ntime .ne. 0) thenсоединяется результат соединения
     if (hist_uni(itwo)%ntime .ne. 0) then
    i=12  ! 
    write(44000,('(7(i12,1x),125(i5,1x))') ) i,knumt_all,nt, nt-hist_uni(itwo)%ntime,  & 
    hist_uni(itwo)%narg, hist_dec(itwo)%n1, hist_dec(itwo)%n2,(numt_all(i),i=1,20)
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
    
            endif ! if ( knumjt > knumj_o ) then    !они встретились!!!
    !!!!!
           
    if ( knumjt < knumj_o ) then   !claster_o разваливается?!!! запишем его, запишем, запишем!

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

 ! a не случится ли так, что отрыв будет похож на старый?           
    if (hist_dec(j)%cluster(1) .ne. 0) then ! the biggest cluster from which j цше come off
            !if (hist_dec(j).cluster(1) .ne. 0) then 
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
          
             if (hist_dec(itwo)%cluster(1).ne.0) then

                    ik=1
                    do k=1,nstat
                    
                    if ( hist_dec(itwo).cluster(k) .ne. 0  ) then
                    if (typ(hist_dec(itwo).cluster(k)) .ne.2) then
                    numj_all(ik)=hist_dec(itwo).cluster(k)         
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
               !  write(41000,'(25(i5,1x))' ) ikk,nt,(hist_uni(j).cluster(ip),ip=1,15)
                   do i=1,knumj_o   !новый
                 hist_dec(numj_o(i))%ntime = nt
                 enddo
        else    !if (ik==-4 .or. ip==-4) then
        
         
       if( ikk == 2 .and. ipp==2) then    !распад распада
        
        if (hist_uni(j).cluster(1) .ne. 0) then ! the biggest cluster from which itwo цше come off
          
                           ip=1
                    do k=1,nstat
                    
                    if ( hist_uni(j).cluster(k) .ne. 0  ) then
                    if (typ(hist_uni(j).cluster(k)) .ne.2) then
                    numj_all(ip)=hist_uni(j).cluster(k)         
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
 i=1 
 t=hist_dec(j)%ntime-hist_uni(j)%ntime

 write(44000,('(7(i12,1x),125(i5,1x))'))  i,knumj_all,nt, hist_dec(j) %ntime-hist_uni(j)%ntime, &
 hist_uni(j)%narg, hist_dec(j)%n1, hist_dec(j)%n2,(numj_all(ik),ik=1,20)

     
      ! para asegurar que la doble notacion no se ocurre  
       t= hist_uni(j)%ntime
   
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
            p=hist_dec(j)%ntime
           ! второй свободен как птица. узаконим это.
           do k=1,knumt_o
          if (hist_uni(numt_o(k))%ntime==t) then
            hist_uni(numt_o(k))%ntime=0
            hist_uni( numt_o(k) )%cluster=0
            hist_uni( numt_o(k) )%n1=0
            hist_uni( numt_o(k) )%n2=0  
            hist_dec( numt_o(k) )%ntime=0
            hist_dec( numt_o(k) )%cluster=0 
                if (knumt_o.ne.1) then
                do i=1,knumt_o
                hist_uni( numt_o(i) )%cluster(i)=numt_o(i)! 
                enddo
                hist_uni(numt_o(k))%ntime=p
                endif
           endif    
           enddo 
     !el primero es libre tambien
      do k=1,knumj_o
            if (hist_uni( numj_o(k) )%ntime==t) then
            hist_uni( numj_o(k) )%ntime=p
            hist_uni( numj_o(k) )%n1=0
            hist_uni( numj_o(k) )%n2=0          
            hist_dec( numj_o(k) )%ntime=0
            hist_uni( numj_o(k) )%cluster=0  ! no estoy muy segura
        do i=1,knumj_o
      hist_uni( numj_o(k) )%cluster(i)=numj_o(i)! 
      enddo
            endif
         
    
      enddo 
          

      endif   !if( ik == 2 .or. ip==2) then    !распад расп
      
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
            
        endif !if ( knumjt < knumj_o ) then   !claster_o разваливается?!!!   
!!!!!!!           
    
    endif   ! if (cyp(j).ne.1
    
        enddo   !form j=1,5000
 
!  .nd now the ar-ar collisions will be counted
nt=nt+1
 !change it after the pro!
!  close(50000)   
!!  close(1,status='delete')
!      enddo    !from do ifile=1,nmax para leerlos
			!	endif   !from if it is the first process							
   
! esta borrando  

!  nttfile=nttfile+1	
! call lammps_command (lmp, 'run 50000 pre no post no')
									    enddo    !from do while .true.

   call finalize_mpi_lammps()
!   CALL lammps_close(lmp);
!    CALL mpi_finalize(ierr)
10 continue
end program 
