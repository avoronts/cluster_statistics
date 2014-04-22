!test simple mpi interface for LAMMPS
! read data from dummp files and from LAMMPS code
! and compare  them
 
program simple_test

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
! end of AE block   >>>>>>>>>>>>>>>>>>>>>>
  print  *, '<<<<<<<<<<< step ', loop, ' done >>>>>>>>>>>>>>>>>'
enddo	! end of main loop "do loop=1,nloop" (true) ahora

!-----------------------------------------------------
  if (allocated(buf)) deallocate(buf)
  IF (lammps_fl == 1) CALL lammps_close(lmp);
  ! close down MPI
  CALL mpi_finalize(ierr)

end program simple_test

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

