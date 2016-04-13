!!!!!!!!!!!!my programm !!!!!!!!!!!!!
program simple

  use LAMMPS_gether
  use history
  use events
  use sort
  use con
  
  implicit none
  

integer nstat,max_atom,nclu,n_entre
parameter(nstat=125,max_atom=25001,nclu=2500,n_entre=15)

double precision mat
parameter(mat = 63.55e-3/6.02/1.6)

!this data for exchange with LAMMPS program

    !integer me - from LAMMPS_gether module
    !integer natoms - from LAMMPS_gether module
    !integer int_buf - from LAMMPS_gether module
    !integer dbl_buf - from LAMMPS_gether module
    integer  typ(max_atom), num_atom, num_vecino4(max_atom), num_vecino5(max_atom)
    double precision  x(max_atom),y(max_atom),z(max_atom), vx(max_atom),vy(max_atom),vz(max_atom)

    double precision, target :: pot1(max_atom), kin1(max_atom),kin2(max_atom), pot2(max_atom) 

    double precision, pointer, dimension (:) :: pot, pot_o, kin, kin_o, tmp_ptr

!this piece of text is devoted to statistical relations betw...
    integer statist(300), npromediar /300/,delta/5/, delta_tau, nstatist /10/
    real*8 mpro,dpro,a
    integer curnum1 /2/
    logical :: sexist


    character*200 sfile 
    character*150 sline

! this data for cluster storage
    integer ncluster(max_atom),ncluster_o(max_atom),mcluster(0:nclu), mcluster_o(0:nclu)
    integer cluster(1:nstat+1,0:nclu),cluster_o(1:nstat+1,0:nclu)
    integer conts(1:nstat+1)
    real conts1(1:nstat+1)
    double precision  energ(max_atom), energ_o(max_atom),vx_cl(max_atom), vy_cl(max_atom), vz_cl(max_atom)

!    integer*4 numj_all(1000), numt_t(1000),numj_o(1000),numjt(1000),numt_o(1000),numt_all(1000)
!    integer   knumj_all, knumt_t, knumj_o, knumjt, knumt_o, knumt_all

    integer cyp(max_atom)

    integer i, j, jj, i1, k, nt, cur_clu , icl, iat, iat1, iat2, icl1, icl2
    integer t, p, n, itwo,ipp, itres, ityp, ikk,ip,iexit,ifile,ik,n1,n2,n3

    real    :: e_part1, e_part2, e_tot


    integer*4  nmax /10000/,natraso /5/,sostav(500)
    !integer cluster7_o(1:nclu,1:nstat),cluster7(1:nclu,1:nstat),mcluster7(0:nclu,2),ncluster7(1:natom),mcluster7_o(0:nclu,2),ncluster7_o(1:natom)

    type (event), pointer :: nev

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!--------------------------------------------------------
CALL init_mpi_lammps('in.1') !!!!!!!!!!!!!!!!!!!!!!!!!!!!1

if (me .eq. 0) then	! init statistics
  nt=0
  mcluster(:) = 0        !numero de atoms en cada cluster
  ncluster(:) = 0        !numero de cluster que mantiene este atom 
  cluster(:,:) = 0          !numero del atom metallico
  vx_cl(:) = 0
  vy_cl(:) = 0
  vz_cl(:) = 0
  energ(:) = 0

  kin1(:) = 0.
  kin2(:) = 0.
  pot1(:) = 0.
  pot2(:) = 0.

  kin => kin1(:)
  pot => pot1(:)
  kin_o => kin2(:)
  pot_o => pot2(:)

  do i=1,5000
    typ(i)=1
  enddo    
  do i=5001,max_atom
     typ(i)=2
  enddo

  call create_history(max_atom)
  call create_con(nstat, 100, 100)
endif

!open(150,file='stat2.dat')
!open(44000,File='r1.dat')      

do while(.true.)
! entro y esta pensando
!								if (me.eq.1) then 

!!!!!!!!!!!!!!!!!!!!!!! read _data !!!!!!!!!!!!!!!!!!!!!!!!!
!  do j=1,m
!    num_atom,typ(num_atom),num_vecino4(num_atom),num_vecino5(num_atom),kin(num_atom), pot(num_atom)  !x(num_atom),y(num_atom),z(num_atom),
!  enddo


    call one_step_data(5)
    if (me .ne. 0) cycle

        tmp_ptr => kin_o
        kin_o => kin
        kin => tmp_ptr
        
        tmp_ptr => pot_o
        pot_o => pot
        pot => tmp_ptr

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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11

   open(5050,File='qqq')    !'//trim(sfile)//'
   write(5050,*) nt			!real*8 kin(1: natom)
   close(5050)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11

 mcluster_o(:) = mcluster(:)
 ncluster_o(:) = ncluster(:)
 cluster_o(:,:) = cluster(:,:)
 energ_o(:)   = energ(:)

 mcluster(:) = 0
 ncluster(:) = 0
 cluster(:,:) = 0
 vx_cl(:) = 0.
 vy_cl(:) = 0.
 vz_cl(:) = 0.
 energ(:) = 0.

 mcluster(0) = 1     ! нулевой кластер описывает все мономеры (размер 1) 
 cur_clu=0

 conts(:) = 0  ! consentration
 do iat = 1,natoms
     if (typ(iat) .ne. 1) cycle
     if (iat .eq. num_vecino5(iat)) then
       conts(1) = conts(1) + 1
       cycle
     endif

     icl = Ncluster(num_vecino5(iat))
!     write(*,*) iat, num_vecino5(iat),typ(iat),icl
     if (icl .eq. 0) then		! add new cluster
            cur_clu = cur_clu+1    ! number of clusters
            icl = cur_clu
            if (typ(num_vecino5(iat)) .eq. 1) then   ! add base atom
               Ncluster(num_vecino5(iat)) = icl
               mcluster(icl) = mcluster(icl)+1
               cluster(mcluster(icl),icl) = num_vecino5(iat)
               vx_cl(icl) = vx_cl(icl) + vx(num_vecino5(iat))
               vy_cl(icl) = vy_cl(icl) + vy(num_vecino5(iat))
               vz_cl(icl) = vz_cl(icl) + vz(num_vecino5(iat))
               energ(icl) = energ(icl) + kin(num_vecino5(iat)) + pot(num_vecino5(iat))
               conts(1) = conts(1) - 1
            endif
     endif

     Ncluster(iat) = icl
     if (mcluster(icl) .lt. nstat) then		! add atom to cluster
        mcluster(icl) = mcluster(icl)+1
        cluster(mcluster(icl),icl) = iat
        vx_cl(icl) = vx_cl(icl) + vx(iat)
        vy_cl(icl) = vy_cl(icl) + vy(iat)
        vz_cl(icl) = vz_cl(icl) + vz(iat)
        energ(icl) = energ(icl) + kin(iat) + pot(iat)
     endif
    ! write(*,*) iat, num_vecino5(iat),icl,Ncluster(iat),mcluster(icl)
     
 enddo    !!num_atom==1,m

do icl = 1, cur_clu
   call QsortC(cluster(1:mcluster(icl),icl))
!   write(*,*) icl, energ(icl), mat*(vx_cl(icl)**2. + vy_cl(icl)**2. + vz_cl(icl)**2. )/mcluster(icl)/2., mcluster(icl)
   energ(icl) = (energ(icl) - mat*(vx_cl(icl)**2. + vy_cl(icl)**2. + vz_cl(icl)**2. )/mcluster(icl)/2.)/mcluster(icl)
!   write(*,*) icl, energ(icl), mcluster(icl)
   if (mcluster(icl) .le. nstat) then
      conts(mcluster(icl)) = conts(mcluster(icl)) + 1
   endif
enddo    !!num_atom==1,m

call update_con(conts,nstat)
!  write (80,*) conts(1:10)

call get_con(conts1)
!      write (80,*) conts1(1:10)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
if(int(1.*nt/1000).eq. 1.*nt/1000) write(*,*) 'stat: time = ', nt

cyp(:) = 0
do j=1,5000  ! todos los atoms

  if (nt .lt. 5) cycle
  if (cyp(j) .eq. 1) cycle      ! skip atoms if they are in known clusters
  if (typ(j) .ne. 1) cycle
  if (mcluster(ncluster(j)) .eq. mcluster_o(ncluster_o(j))) cycle

!---------------------
  if (ncluster(j) .eq. 0) cluster(1,ncluster(j)) = j  ! мономер
  if (ncluster_o(j) .eq. 0) cluster_o(1,ncluster_o(j)) = j  ! мономер
  
  call add_event(nev,cluster(:,ncluster(j)), mcluster(ncluster(j)),cluster_o(:,ncluster_o(j)),mcluster_o(ncluster_o(j)),nt)

  do i=1,size(nev%atoms)	!mark atoms as seen
      iat = nev%atoms(i)   !
      cyp(iat) = 1
  enddo

  if (nev%fusion.eq. 1) then
    nev%e_tot = energ(ncluster(nev%atoms(1)))
    nev%e_part1 = energ_o(ncluster_o(nev%atoms(1)))
    nev%e_part2 = energ_o(ncluster_o(nev%atoms(nev%n1+1)))
  else
    nev%e_tot = energ_o(ncluster_o(nev%atoms(1)))
    nev%e_part1 = energ(ncluster(nev%atoms(1)))
    nev%e_part2 = energ(ncluster(nev%atoms(nev%n1+1)))
  endif

  nev%c1 = conts1(nev%n1)
  nev%c2 = conts1(nev%n2)
  call update_history_check(nev) !------ check for loops and update history---------


enddo   !form j=1,5000

enddo    !from do while .true. <<<<<< end of main circle


call finalize_mpi_lammps()
10 continue

end program 

