module LAMMPS_gether


   use MPI
   use LAMMPS
   use, intrinsic :: ISO_C_binding, only : C_double, C_ptr, C_int

   implicit none

   private
     integer :: ierr, nprocs,lammps_fl     ! for MPI

     type (C_ptr) :: lmp		! for lammps communication
     integer :: nprocs_main, nprocs_lammps, comm_lammps
     integer :: icall = 0
     integer :: ncall = 10000       ! number of calling between saving information

   integer, public :: me    ! for MPI

   integer, public :: natoms = 0
   integer, dimension(:), allocatable,public :: int_buf
   double precision, dimension(:), allocatable,public :: dbl_buf
   
   integer :: int_peratom = 4
   integer :: dbl_peratom = 8
   
     
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
! deallocate arrays
! finalize lammps
! finalize MPI
subroutine finalize_mpi_lammps()

   if (allocated(int_buf)) deallocate(int_buf)
   if (allocated(dbl_buf)) deallocate(dbl_buf)
   if (lammps_fl == 1) call lammps_close(lmp);
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

subroutine one_step_data (step)

    implicit none

    integer :: step 
    integer :: status(MPI_STATUS_SIZE)
    integer :: i,j,j_int, j_dbl,iproc,ierr
    character*15 :: str_step

    real (C_double), dimension(:), pointer :: comp => NULL()
    real (C_double), dimension(:), pointer :: comp1 => NULL()
    real (C_double), dimension(:), pointer :: comp2 => NULL()
    real (C_double), dimension(:), pointer :: comp3 => NULL()

    real (C_double), dimension(:,:), pointer :: x => NULL()
    real (C_double), dimension(:,:), pointer :: v => NULL()

    integer (C_int), dimension(:), pointer :: tag => NULL()
    integer (C_int), dimension(:), pointer :: type => NULL()

! every processor has number of atoms. me=0 takes total number of atoms
    write(str_step,'(i0)') step

    natoms = 0
    if (lammps_fl == 1)  then
       if ((icall .eq. ncall) .and. (ncall .gt. 0)) then
           icall = 0
           call lammps_command(lmp,'write_restart restart_file1')
       endif
       icall = icall +1
!     call lammps_command(lmp,'fix fff all ave/atom 1 1 '//trim(str_step)//' c_pe1')

       call lammps_command(lmp,'run '//trim(str_step)//' pre no post no')
       call lammps_extract_atom (tag, lmp, 'id')
       call lammps_extract_atom (type, lmp, 'type')
       call lammps_extract_compute (comp, lmp, 'clu8', 1, 1)
       call lammps_extract_compute (comp1, lmp, 'clu5', 1, 1)
       call lammps_extract_compute (comp2, lmp, 'ke1', 1, 1)
       call lammps_extract_compute (comp3, lmp, 'pe1', 1, 1)
call lammps_extract_atom (x, lmp, 'x')
call lammps_extract_atom (v, lmp, 'v')

!       call lammps_extract_compute (comp, lmp, 'clu', 1, 2) ! 2d array
       natoms = size(tag)
    endif
    call MPI_Reduce(natoms, i, 1, MPI_INT,MPI_SUM, 0, MPI_COMM_WORLD ,ierr)
    if (me == 0)  natoms = i

!    print *, 'done',me, natoms
!  call MPI_BARRIER(MPI_COMM_WORLD, ierr)

!   send all information to me=0 node
    if (.not.allocated(int_buf)) allocate(int_buf(int_peratom*natoms)) 
    if (size(int_buf) .ne. int_peratom*natoms) then
       deallocate(int_buf)
       allocate (int_buf(int_peratom*natoms))
    endif
    
    if (.not.allocated(dbl_buf)) allocate(dbl_buf(dbl_peratom*natoms)) 
    if (size(dbl_buf) .ne. dbl_peratom*natoms) then
       deallocate(dbl_buf)
       allocate (dbl_buf(dbl_peratom*natoms))
    endif

    if (me .ne. 0)  then
        j_int=1
        j_dbl=1
        do i=1, natoms 
!      print *, tag(i),type(i),comp(:,i)
          int_buf(j_int)= tag(i)
          int_buf(j_int+1)= type(i)
          int_buf(j_int+2)= idnint(comp(i))
          int_buf(j_int+3)= idnint(comp1(i))
    !      buf(j+3)= idnint(comp(2,i))
          j_int=j_int+4

          dbl_buf(j_dbl)= comp2(i)
          dbl_buf(j_dbl+1)=comp3(i)
dbl_buf(j_dbl+2)=x(1,i)
dbl_buf(j_dbl+3)=x(2,i)
dbl_buf(j_dbl+4)=x(3,i)

dbl_buf(j_dbl+5)=v(1,i)
dbl_buf(j_dbl+6)=v(2,i)
dbl_buf(j_dbl+7)=v(3,i)
          j_dbl=j_dbl+dbl_peratom
        enddo

!!send inforation to prosess with rank 0
        call MPI_Send(int_buf, int_peratom*natoms, MPI_INT, 0, me, MPI_COMM_WORLD, ierr)
        deallocate(int_buf)

        call MPI_Send(dbl_buf, dbl_peratom*natoms, MPI_DOUBLE, 0, me, MPI_COMM_WORLD, ierr)
        deallocate(dbl_buf)
        
    else		! I am rank 0 processor!!!
!! receive information from other nodes in one bufer
        j=1
        do iproc = 1,nprocs-1 
          call MPI_Recv(int_buf(j:), int_peratom*natoms, MPI_INT,  iproc, iproc, MPI_COMM_WORLD, status, ierr)
          call MPI_Get_count( status,  MPI_INT, i, ierr )
!      print *, me, iproc, i,j,(j+3),(j+i+3)
          j=j+i
        enddo
!! receive information from other nodes in one bufer
        j=1
        do iproc = 1,nprocs-1 
          call MPI_Recv(dbl_buf(j:), dbl_peratom*natoms, MPI_DOUBLE,  iproc, iproc, MPI_COMM_WORLD, status, ierr)
          call MPI_Get_count( status,  MPI_DOUBLE, i, ierr )
!      print *, me, iproc, i,j,(j+3),(j+i+3)
          j=j+i
        enddo

    endif

end subroutine one_step_data




end module LAMMPS_gether