!test simple mpi interface for LAMMPS
! read data from dummp files and from LAMMPS code
! and compare  them
 
program simple_test

    use LAMMPS_gether

    implicit none

    integer natom
    parameter(natom=25001)
    integer :: nloop,iloop
    character*150 sline,sfile
    integer*4 typ(1:natom),num_vecino4(1:natom),num_vecino5(1:natom)
    real*8 pot(1 :natom),kin(1: natom)
    integer*4 num_atom, m, i,j,id,j_dbl

!<<<<<<<<<<<<<<<< program start here !  >>>>>>>>>>>>>>>>>>
    call init_mpi_lammps('test.in')
    if (me == 0) print *, '<<<<<<<<<<< mpi_init done >>>>>>>>>>>>>>>>>'

!<<<<<<<<<<<<<<<< main loop >>>>>>>>>>>>>>>>>>
nloop = 250000
do iloop = 5,nloop,5
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
    call one_step_data(5)
!    if (me == 0) print  *, '<<<<<<<<<<< one step done >>>>>>>>>>>>>>>>>'
    
!!!!!!!!!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
    if (me == 0) then
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
!        close(50000)
    
!        print  *, '<<<<<<<<<<< read dump done  >>>>>>>>>>>>>>>>>'
        
        j=1
        j_dbl=1
        do i=1,m 
          id=int_buf(j)
          if (int_buf(j+1) .ne. typ(id))   print  *, 'bad type in', j,id,int_buf(j+1),typ(id)
          if (int_buf(j+2) .ne. num_vecino4(id)) print  *, 'bad vencino4', id,int_buf(j+2),num_vecino4(id)
          if (int_buf(j+3) .ne. num_vecino5(id)) print  *, 'bad vencino5', id,int_buf(j+3),num_vecino5(id)
          j=j+4
          
          if (abs(dbl_buf(j_dbl)-kin(id)) .gt. 1e-5) print  *, 'bad kin', id,dbl_buf(j_dbl),kin(id)
          if (abs(dbl_buf(j_dbl+1)- pot(id)) .gt. 1e-5) print  *, 'bad pot', id,dbl_buf(j_dbl+1),pot(id)
          j_dbl=j_dbl+8
        enddo
        print  *, '<<<<<<<<<<< check done  >>>>>>>>>>>>>>>>>'
        
    endif   

enddo       !<<<<<<<<<<<<<<<< end of main loop >>>>>>>>>>>>>>>>>>
!-----------------------------------------------------
    call finalize_mpi_lammps()

end program simple_test

