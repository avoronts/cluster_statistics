module LAMMPS_gether


   use MPI
   use LAMMPS
   use, intrinsic :: ISO_C_binding, only : C_double, C_ptr, C_int

   implicit none

   private
     integer :: ierr, nprocs,lammps_fl     ! for MPI

     type (C_ptr) :: lmp		! for lammps communication
     integer :: nprocs_main, nprocs_lammps, comm_lammps
     
    integer :: natoms

   integer, public :: me    ! for MPI

   integer, dimension(:), allocatable,public :: buf
   integer :: int_peratom = 4
   integer :: dbl_peratom = 2
   
     
   public ::  init_mpi_lammps, finalize_mpi_lammps, one_step_data


contains
!--------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Init MPI
! split processors to main and LAMMPS
! Init LAMMPS

subroutine init_mpi_lammps (file_name)
!   use MPI
!   use LAMMPS
!   use, intrinsic :: ISO_C_binding, only : C_double, C_ptr, C_int

   implicit none
   
   character(*), intent(in) :: file_name

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
     call lammps_file (lmp, file_name)
   end if
end subroutine init_mpi_lammps

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11111
! finalize lammps
! deallocate arrays
! finalize MPI
subroutine finalize_mpi_lammps()

   if (allocated(buf)) deallocate(buf)
   IF (lammps_fl == 1) CALL lammps_close(lmp);
  ! close down MPI
   CALL mpi_finalize(ierr)
end subroutine finalize_mpi_lammps

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! make one step of LAMMPS simulation 
! buf (out) - output data
! list - 
! natoms (out) - number of atoms
! lmp (in) - LAMMPS pointer
! me (in) - processor rank
! lammps_fl (in) - if LAMMPS runing on this processor
! nprocs (in) - number of processors

subroutine one_step_data (natoms)
!    use MPI
!    use LAMMPS
!    use, intrinsic :: ISO_C_binding, only : C_double, C_ptr, C_int
    implicit none

    integer :: status(MPI_STATUS_SIZE)
    integer, intent(out) :: natoms
!    integer, dimension(:), allocatable, intent(inout) :: buf


    integer :: i,j,iproc,ierr
 !  integer, dimension(:), allocatable,target :: buf
    real (C_double), dimension(:), pointer :: comp => NULL()
    real (C_double), dimension(:), pointer :: comp1 => NULL()
    
    integer (C_int), dimension(:), pointer :: tag => NULL()
    integer (C_int), dimension(:), pointer :: type => NULL()

! every processor has number of atoms. me=0 takes total number of atoms
    natoms = 0
    if (lammps_fl == 1)  then
       call lammps_command(lmp,'run 5 pre no post no')
       call lammps_extract_atom (tag, lmp, 'id')
       call lammps_extract_atom (type, lmp, 'type')
       call lammps_extract_compute (comp, lmp, 'clu8', 1, 1)
       call lammps_extract_compute (comp1, lmp, 'clu5', 1, 1)
!       call lammps_extract_compute (comp, lmp, 'ke1', 1, 1)
!       call lammps_extract_compute (comp, lmp, 'pe1', 1, 1)
!       call lammps_extract_compute (comp, lmp, 'clu', 1, 2) ! 2d array
       natoms = size(tag)
    endif
    call MPI_Reduce(natoms, i, 1, MPI_INT,MPI_SUM, 0, MPI_COMM_WORLD ,ierr)
    if (me == 0)  natoms = i

!    print *, 'done',me, natoms
!  call MPI_BARRIER(MPI_COMM_WORLD, ierr)

! send all information to me=0 node
!    buf_step=size(buf)/natoms

    
    if (.not.allocated(buf)) allocate (buf(int_peratom*natoms)) 
    if (size(buf) .ne. int_peratom*natoms) then
       deallocate(buf)
       allocate (buf(int_peratom*natoms))
    endif

    if (me .ne. 0)  then
        j=1
        do i=1, natoms 
!      print *, tag(i),type(i),comp(:,i)
          buf(j)= tag(i)
          buf(j+1)= type(i)
          buf(j+2)= idnint(comp(i))
          buf(j+3)= idnint(comp1(i))
    !      buf(j+3)= idnint(comp(2,i))
          j=j+4
        enddo

!!send inforation to prosess with rank 0
        call MPI_Send(buf, int_peratom*natoms, MPI_INT, 0, me, MPI_COMM_WORLD, ierr)
        deallocate(buf)
    else		! I am rank 0 processor!!!
!! receive information from other nodes in one bufer
        j=1
        do iproc = 1,nprocs-1 
          call MPI_Recv(buf(j:), int_peratom*natoms, MPI_INT,  iproc, iproc, MPI_COMM_WORLD, status, ierr)
          call MPI_Get_count( status,  MPI_INT, i, ierr )
!      print *, me, iproc, i,j,(j+3),(j+i+3)
          j=j+i
        enddo
    endif
end subroutine one_step_data




end module LAMMPS_gether