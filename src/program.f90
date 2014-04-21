program cluster_count

   use MPI
   use LAMMPS
   use, intrinsic :: ISO_C_binding, only : C_double, C_ptr, C_int

   implicit none

!-----------------------------------------------------------
interface one_step_data    ! interface  for geting info from lammps nodes
  subroutine one_step_data (data, s_list,natoms, lmp, me,lammps_fl,nprocs)
     use, intrinsic :: ISO_C_binding, only : C_double, C_ptr, C_int
     implicit none
     integer , dimension(:), allocatable,  intent(out) :: data, s_list 
     type (C_ptr), intent(in) :: lmp
     integer, intent(in) :: me, lammps_fl,nprocs
     integer, intent(out) :: natoms
  end subroutine one_step_data     
end interface one_step_data

interface part_my    ! ! interface for find common part of two ordered set of numbers (clusters)
  subroutine part_my(clust1,clust2,comm,part1,part2)
    implicit none
    integer, dimension(:), intent(in):: clust1,clust2
    integer, dimension(:), allocatable, intent(out):: comm, part1,part2
!    integer, intent(in) :: s1,s2
  end subroutine part_my
end interface part_my

!----------------------------------------------------------
   logical :: in_my
   integer :: ierr, me, nprocs    ! for MPI

   type (C_ptr) :: lmp		! for lammps communication
   integer :: lammps_fl, nprocs_main, nprocs_lammps, comm_lammps
   integer, dimension(:), allocatable :: buf, s_list
   integer :: buf_step, ibuf, natoms, nclusters, mv

!--------------------------------------------------------
   integer :: ind, first, color_tmp , n_cl, icl
   integer, dimension(:), pointer :: clust_all, color_at, size_at,type_at
   integer, dimension(:), pointer :: clust_all_o, color_at_o, size_at_o, type_at_o

   double precision :: di

integer :: i,j,loop, nloop

   integer :: ibeg, iend, ibeg_o,iend_o, iat, jat,nat,kat, time, flag
   integer, dimension(:), pointer :: tmp
   integer, dimension(:), allocatable :: comm,part1,part2

integer nstat,nclu,natom
parameter(nstat=125,natom=25001,nclu=2500)

    type hist
            integer*4 time  
            integer*4 n1  
            integer*4 n2 
            real*4 stat(nstat)
            integer, dimension(:), allocatable :: cluster 
    end type
    type (hist), dimension(:), pointer :: hist_all => NULL(), hist_tmp => NULL()

 
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!! AE data!!

logical exists
    real*4  statinum(nstat), para_esc(nstat),pott,kinn
    integer eof /0/,typ(1:natom),cyp(natom),z,y,icolor,k,Code,cur_rec ,iras,np(12)
    character,allocatable:: struu(:)
    character*5 s1
    character*6 s2,c2
    character*200 sfile 
    character*150 sline
    integer*4 num_vecino(1 :5),ncluster(1:natom),mcluster(0:nclu,2),ncluster_o(1:natom),mcluster_o(0:nclu,2)  ! 1- только металл, 2 - есть аргон
!integer*8 ibeg(1:natom), n1, n2
    integer*4 num_atom,nt, M, cur_num ,cur_clu, numj_all(1000), numj_o(1000),numjt(1000),numt_o(1000),numt_all(1000)
    integer*4  nmax /10000/,natraso /5/
    integer*4 ik,ikk,ikp
  
    real*8 pot(1 :natom),kin(1: natom)
character*17 sdir/'c:\1500_1\for_me\'/
    character*2250  cluster(0:nclu),cluster_o(0:nclu),clupro
    
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  ! setup MPI and various communicators
  ! driver runs on all procs in MPI_COMM_WORLD
  ! comm_lammps only has 1st P procs (could be all or any subset)

  CALL mpi_init(ierr)
  CALL mpi_comm_rank(MPI_COMM_WORLD,me,ierr);
  CALL mpi_comm_size(MPI_COMM_WORLD,nprocs,ierr);

!  narg = command_argument_count()
!  IF (narg /= 2) THEN
!     PRINT *, 'Syntax: f_driver P in.lammps'
!     CALL mpi_abort(MPI_COMM_WORLD,1,ierr)
!  END IF

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


open(15, file='out.dat')
!!!!!!!!! main loop begin here !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
nloop = 100
if (me == 0) print  *, 'make ', nloop, ' steps'
do loop = 1,nloop
!!! nt = loop*timestep for AE

  call one_step_data (buf, s_list,natoms, lmp, me,lammps_fl,nprocs)
!! end communication let's do separate work,
!!  meet at MPI_Reduce point (number of atoms)

  if (me .ne.  0)  cycle  ! Only the me=0 processor make statistics
!! dont foget to deallocale "buf, s_list" on the me=0 processor!!!

!!! take data from bufer and deallocate it
     print *, 'ok',loop,natoms, nclusters

  buf_step=size(buf)/natoms
  mv = maxval(buf(1::buf_step)) ! max atom id
  
  allocate(color_at(mv))
  allocate(size_at(mv))
  allocate(type_at(mv))
  allocate(clust_all(natoms))
  
  if (.not.associated(hist_all)) then	! increase history size (if necessary)
    allocate(hist_all(mv))
    do i = 1,mv
      hist_all(i)%time = 0
      hist_all(i)%stat(:) = 0
    enddo
  elseif (size(hist_all).lt.mv) then
    hist_tmp => hist_all
    allocate(hist_all(mv))
    do i = 1,mv
      hist_all(i)%time = 0
      hist_all(i)%stat(:) = 0
    enddo
    hist_all(:size(hist_tmp)) = hist_tmp(:)
    deallocate(hist_tmp)
  endif

  ind = 1
  n_cl=0
  di=0
  do i =1, natoms
    j=(s_list(i)-1)*buf_step
    clust_all(ind) = buf(j+1)
    if (color_tmp .ne. buf(j+3) .or. (i .eq. 1) ) then
      first = ind
      color_tmp = buf(j+3)
      n_cl = n_cl+1		! increace when new cluster start
    endif
    color_at(buf(j+1)) = first   ! atom is in the cluster started "fom fist"
    size_at(buf(j+1)) = buf(j+4) ! size of cluster
    type_at(buf(j+1)) = buf(j+2) ! type of atom
    if (buf(j+4) .ne. 0) di = di +1.d0/buf(j+4)
    ind = ind + 1
  enddo 
!    do i=1,natoms
!      j=(s_list(i)-1)*4
!      j=i !clust_all(i)
!      print *, me,color_at(j),size_at(j),type_at(j)
!    enddo
  deallocate(buf,s_list)
  print *,n_cl,di
!----  done reading cluster data from bufer ---------------
!print *, 'go'

!<<<<<<<<<<<<<<<<<<< AE data <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

!call seg_argon(nclu,nstat,natoms,ncluster7,cluster7,mcluster7,
!                                 ncluster7_o,cluster7_o,mcluster7_o,
!                                 ncluster,mcluster,cluster,kin,pot)
!call seg_met(nclu,nstat,natoms,ncluster,cluster,mcluster,
!                               ncluster_o,cluster_o,mcluster_o)
!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

enddo	! end of main loop "do loop=1,nloop"
close(15)
!-----------------------------------------------------
  if (allocated(buf)) deallocate(buf)
  IF (lammps_fl == 1) CALL lammps_close(lmp);
  ! close down MPI
  CALL mpi_finalize(ierr)

end program cluster_count

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! make one step of LAMMPS simulation 
! buf (out) - output data
! list - 
! natoms (out) - number of atoms
! lmp (in) - LAMMPS pointer
! me (in) - processor rank
! lammps_fl (in) - if LAMMPS runing on this processor
! nprocs (in) - number of processors
subroutine one_step_data (buf, list, natoms, lmp, me,lammps_fl,nprocs)
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
   integer , dimension(:), allocatable  :: list

! every processor has number of atoms. me=0 takes total number  of atoms
  natoms = 0
  if (lammps_fl == 1)  then
    call lammps_command(lmp,'run 5 pre no post no')
    call lammps_extract_atom (tag, lmp, 'id')
    call lammps_extract_atom (type, lmp, 'type')
!    print *, me, 'comp' 
    call lammps_extract_compute (comp, lmp, 'clu', 1, 2)
!    print *, me, 'comp ok ' 
    natoms = size(tag)
  endif
  call MPI_Reduce(natoms, i, 1, MPI_INT,MPI_SUM, 0, MPI_COMM_WORLD ,ierr)
  if (me == 0)  natoms = i
!  print *, 'done',me, natoms

! send all information to me=0 node
  if (allocated(buf)) deallocate(buf)
  allocate (buf(4*natoms))
  allocate (list(natoms))
  do i=1,natoms
      list(i)=i
  enddo

  if (me .ne. 0)  then
    !!!! let's sort local data by cluster ID!!!!
    do i=1,natoms
      buf(i)=idnint(comp(1,i))
    enddo
    call sort_my(list,buf,tag,1,natoms+1)
!    print *, me, 'sort ok ' 
    !!!!!  pack data to buffer !!!!!!!!!!!!!!
    j=1
    do i=1, natoms 
!      print *, tag(i),type(i),comp(:,i)
      buf(j)= tag(list(i))
      buf(j+1)= type(list(i))
      buf(j+2)= idnint(comp(1,list(i)))
      buf(j+3)= idnint(comp(2,list(i)))
      j=j+4
    enddo

!!send inforation to prosess with rank 0
    call MPI_Send(buf, 4*natoms, MPI_INT, 0, me, MPI_COMM_WORLD, ierr)
    deallocate(buf)
    deallocate(list)
  else		! I am rank 0 processor!!!
!! receive information from other nodes in one bufer
    j=1
    do iproc = 1,nprocs-1 
      call MPI_Recv(buf(j:), 4*natoms, MPI_INT,  iproc, iproc, MPI_COMM_WORLD, status, ierr)
      call MPI_Get_count( status,  MPI_INT, i, ierr )
!      print *, me, iproc, i,j,(j+3),(j+i+3)
      !!!! let's sort data from different processors!!!!
      call merge_my(list,buf(3::4),buf(1::4),1,(j+3)/4,(j+3)/4+i/4)
!      print *, me, iproc,'sort done'
      j=j+i
    enddo
!    do i=1,natoms
!      j=(list(i)-1)*4
!      print *, me,buf(j+1),buf(j+2),buf(j+3),buf(j+4)
!    enddo
!    deallocate(list)
  endif
contains
end subroutine one_step_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! sort pat of bufer 'buf' between elements 'ibeg' and 'iend'
! according comparison parameter 'comp1' and 'comp2'
! buf - array for sorting
! comp1 - first sorting criterion
! comp2 - second sorting criterion
! ibeg, iend - rnge of sorting
 
recursive subroutine sort_my(buf,comp1,comp2,ibeg,iend)
implicit none

integer, dimension(*), intent(inout):: buf
integer, dimension(*), intent(in):: comp1,comp2
integer, intent(in):: ibeg, iend
integer :: level
integer :: imid

if (iend-ibeg .le. 1) then
  return
else
   imid = (iend+ibeg)/2
   call sort_my(buf,comp1,comp2,ibeg,imid)
   call sort_my(buf,comp1,comp2,imid,iend)
   call merge_my(buf,comp1,comp2,ibeg,imid,iend)
end if
end subroutine sort_my

!------------------------------------------------------
! merge with sorting of two parts of indexed array
! buf  - indexes of elements
! comp1 (comp2) - sorting  criterion
! ibeg - index of the begin element
! imid - index of the separator
! iend - index of the end  element

subroutine merge_my(buf,comp1,comp2,ibeg,imid,iend)
implicit none

integer, dimension(*), intent(inout):: buf
integer, dimension(*), intent(in):: comp1, comp2
integer, intent(in) :: ibeg, imid, iend
integer :: l_beg, l_end, r_beg, r_end, j,i
integer, dimension(iend-ibeg) :: tmp

l_beg=ibeg
l_end=imid
r_beg=imid
r_end=iend

j=1
do while ((l_end .gt. l_beg) .and. (r_end .gt. r_beg))
  if (comp1(buf(l_beg)) .lt. comp1(buf(r_beg))) then    !first condition of sorting (comp1)
    tmp(j)= buf(l_beg)
    j=j+1
    l_beg =l_beg+1
  elseif (comp1(buf(l_beg)) .gt. comp1(buf(r_beg))) then !first condition of sorting (comp1)
    tmp(j)= buf(r_beg)
    j=j+1
    r_beg = r_beg+1
  elseif (comp2(buf(l_beg)) .le. comp2(buf(r_beg))) then  !first condition the same, second condition (comp2)
    tmp(j)= buf(l_beg)
    j=j+1
    l_beg =l_beg+1
  else						 !first condition the same, second condition (comp2)
    tmp(j)= buf(r_beg)
    j=j+1
    r_beg = r_beg+1
  end if
enddo

do while (l_end .gt. l_beg)
   tmp(j) = buf(l_beg)
   j=j+1
   l_beg=l_beg+1
enddo

do while (r_end .gt. r_beg) 
   tmp(j) = buf(r_beg)
   j=j+1
   r_beg=r_beg+1
enddo

j=1
do i=ibeg,iend-1
  buf(i)=tmp(j)
  j=j+1
enddo
end subroutine merge_my

!------------------------------------------------------
! find common part of two ordered set of numbers (clusters)
! clust1, clust2 - odered sets (clusters)
! comm - array of common elements (sorted)
! part1 - unique part of first cluster (sorted)
! part2 - unique part of second cluster (sorted)

subroutine part_my(clust1,clust2,comm,part1,part2)
implicit none

integer, dimension(:), intent(in):: clust1,clust2
integer, dimension(:), allocatable, intent(out):: comm, part1,part2
integer  :: s1,s2

integer, dimension(size(clust1)) :: comm_tmp, part1_tmp
integer, dimension(size(clust2)) :: part2_tmp
integer :: l_beg, l_end, r_beg, r_end, j,j1,j2

s1=size(clust1)
s2=size(clust2)

l_beg=1
l_end=s1
r_beg=1
r_end=s2

j=1
j1=1
j2=1
do while ((l_end .gt. l_beg) .and. (r_end .gt. r_beg))
  if (clust1(l_beg) .eq. clust2(r_beg)) then    ! atom in both clusters
    comm_tmp(j)= clust1(l_beg)
    j=j+1
    l_beg =l_beg+1
    r_beg =r_beg+1
  elseif (clust1(l_beg) .lt. clust2(r_beg)) then ! atom only in fist cluster
    part1_tmp(j1)= clust1(l_beg)
    j1=j1+1
    l_beg =l_beg+1
  else 				! atom only in second cluster
    part2_tmp(j2)= clust2(r_beg)
    j2=j2+1
    r_beg =r_beg+1
  end if
enddo

do while (l_end .gt. l_beg)  ! tail of fist cluster
  part1_tmp(j1)= clust1(l_beg)
  j1=j1+1
  l_beg =l_beg+1
enddo

do while (r_end .gt. r_beg)  ! tail of second cluster
  part2_tmp(j2)= clust2(r_beg)
  j2=j2+1
  r_beg =r_beg+1
enddo

allocate(comm(j-1))
comm=comm_tmp(:j-1)
allocate(part1(j1-1))
part1=part1_tmp(:j1-1)
allocate(part2(j2-1))
part2=part2_tmp(:j2-1)

end subroutine part_my

!------------------------------------------------------
! compare two clusters
integer function bigger_my(clust1,clust2)
implicit none
integer, dimension(:), intent(in):: clust1,clust2
integer, dimension(:), allocatable :: comm, part1,part2

call part_my(clust1,clust2,comm,part1,part2)

if ((size(part1) .ne. 0) .and. (size(part2) .ne. 0)) then
  bigger_my = -2			! wery diffeent
elseif (size(part1) .ne. 0) then
  bigger_my = 1				! clust1 is bigger
elseif (size(part2) .ne. 0) then
  bigger_my = -1			! clust2 is bigger
else
  bigger_my = 0			! are equal
endif
deallocate(comm, part1,part2)

end function bigger_my

!------------------------------------------------
! find by bisection method 
! if the elemnt 'at' in sorted set 'clust'
! return 'true' if 'at' in 'clust'
!
logical function in_my(at,clust)
implicit none
integer, dimension(:), intent(in):: clust
integer, intent(in) :: at
integer :: ibeg, iend, imid

ibeg=1
iend=size(clust)
in_my = .false.

do while (iend - ibeg .gt. 1)
  imid = (iend - ibeg)/2
  if (clust(imid) .le. at) then
    ibeg = imid
  else 
    iend = imid
  endif
enddo
if (clust(ibeg) .eq. at) in_my = .true. 
if (clust(iend) .eq. at) in_my = .true. 

end function in_my

