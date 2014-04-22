module LAMMPS_gether


   use MPI
   use LAMMPS
   use, intrinsic :: ISO_C_binding, only : C_double, C_ptr, C_int

   implicit none

   private


!interface   one_step_data  ! interface  for geting info from lammps nodes
!  subroutine one_step_data_sorted (data, s_list,natoms, lmp, me,lammps_fl,nprocs)
!     use, intrinsic :: ISO_C_binding, only : C_double, C_ptr, C_int
!     implicit none
!     integer , dimension(:), allocatable,  intent(out) :: data, s_list 
!     type (C_ptr), intent(in) :: lmp
!     integer, intent(in) :: me, lammps_fl,nprocs
!     integer, intent(out) :: natoms
 ! end subroutine one_step_data_sorted
  
!  subroutine one_step_data_simple (data, natoms, lmp, me,lammps_fl,nprocs)
!     use, intrinsic :: ISO_C_binding, only : C_double, C_ptr, C_int
!     implicit none
!     integer , dimension(:), allocatable,  intent(out) :: data
!     type (C_ptr), intent(in) :: lmp
!     integer, intent(in) :: me, lammps_fl,nprocs
!     integer, intent(out) :: natoms
 ! end subroutine one_step_data_simple
!end interface one_step_data


   public :: one_step_data
!   public :: lammps_instance, C_ptr, C_double, C_int


contains
!--------------------------------------------
!interface part_my    ! ! interface for find common part of two ordered set of numbers (clusters)
!  subroutine part_my(clust1,clust2,comm,part1,part2)
!    implicit none
!    integer, dimension(:), intent(in):: clust1,clust2
!    integer, dimension(:), allocatable, intent(out):: comm, part1,part2
!!    integer, intent(in) :: s1,s2
!  end subroutine part_my
!end interface part_my

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   logical :: in_my
!   integer :: ierr, me, nprocs    ! for MPI

!   type (C_ptr) :: lmp		! for lammps communication
!   integer :: lammps_fl, nprocs_main, nprocs_lammps, comm_lammps
!   integer, dimension(:), allocatable :: buf, s_list
!   integer :: buf_step, ibuf, natoms, nclusters, mv
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Init MPI, split processors to main and LAMMPS
! Init LAMMPS
subroutine init_all (in_file, lmp, me,lammps_fl,nprocs)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   use MPI
   use LAMMPS
   use, intrinsic :: ISO_C_binding, only : C_double, C_ptr, C_int

   implicit none
   integer :: ierr, me, nprocs    ! for MPI

   type (C_ptr) :: lmp		! for lammps communication
   integer :: lammps_fl, nprocs_main, nprocs_lammps, comm_lammps
   integer, dimension(:), allocatable :: buf
   integer :: buf_step, natoms, loop, nloop, mv

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
end subroutine init_all

subroutine final_all()
   if (allocated(buf)) deallocate(buf)
  IF (lammps_fl == 1) CALL lammps_close(lmp);
  ! close down MPI
  CALL mpi_finalize(ierr)
end subroutine
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
end subroutine one_step_data




end module LAMMPS_gether