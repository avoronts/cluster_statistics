!test simple mpi interface for LAMMPS
! read data from dummp files and from LAMMPS code
! and compare  them
 
program simple_test

!   use MPI
!   use LAMMPS
    use LAMMPS_gether
!   use, intrinsic :: ISO_C_binding, only : C_double, C_ptr, C_int

    implicit none
    integer :: natoms
!    integer, dimension(:), allocatable :: buf
    integer nstat,nclu,natom, n_entre,nloop,iloop
    parameter(nstat=125,natom=25001,nclu=2500,n_entre=15)

    character*150 sline,sfile
    integer*4 typ(1:natom),num_vecino4(1:natom),num_vecino5(1:natom)
    integer*4 num_atom, M, i,j,id
    real*8 pot(1 :natom),kin(1: natom)


    call init_mpi_lammps('in.1')
    if (me == 0) then
       print  *, '<<<<<<<<<<< mpi_init done >>>>>>>>>>>>>>>>>'
    endif

!!-------  AE data-------
!<<<<<<<<<<<<<<<< put everythig befor main loop here >>>>>>>>>>>>>>>>>>
! open(15,File='t'//trim(sfile)//'.dat')      
nloop = 25
do iloop=5,nloop,5
!   open(1,File=sfile)    !'//trim(sfile)//'
!-------------------------- 
! end of AE block   >>>>>>>>>>>>>>>>>>>>>>
    call one_step_data(natoms)
!!!!!!!!! main loop begin here !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
    if (me == 0) then
        print  *, '<<<<<<<<<<< one step done >>>>>>>>>>>>>>>>>',me 
        write(sfile,'(i0)') iloop
        open(50000,File='dmp_clust'//trim(sfile)//'.txt')    !'//trim(sfile)//'
        do i =1,3
	  read(50000,*) sline     !lea las demas lineas no necesarias
	enddo
	Read(50000,*) m  !
	do i = 5,9
	  read(50000,*) sline     !lea las demas lineas no necesarias
	enddo
         !lea la primera linea
      ! do while (.not. EOF(1))
        do j=1,m
    	Read(50000,*) num_atom,typ(num_atom),num_vecino4(num_atom),num_vecino5(num_atom),kin(num_atom), pot(num_atom)
        enddo
        close(50000,status='delete')
        print  *, '<<<<<<<<<<< read dump done  >>>>>>>>>>>>>>>>>'
        
        j=1
        do i=1,m 
          id=buf(j)
          if (buf(j+1) .ne. typ(id)) then
            print  *, 'bad type in', j,id,buf(j+1),typ(id)
          endif
          if (buf(j+2) .ne. num_vecino4(id)) then
            print  *, 'bad vencino4', id,buf(j+2),num_vecino4(id)
          endif
          if (buf(j+3) .ne. num_vecino5(id)) then
            print  *, 'bad vencino5', id,buf(j+3),num_vecino5(id)
          endif
          j=j+4
        enddo
        print  *, '<<<<<<<<<<< check done  >>>>>>>>>>>>>>>>>'
        
    endif

enddo
!-----------------------------------------------------
    call finalize_mpi_lammps()

end program simple_test

