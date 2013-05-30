      module obj_fun
     
      type :: ener_ref_type
          integer :: id
          character(80) :: mol
          character(80) :: mol_ref1
          character(80) :: mol_ref2
          character(80) :: mol_ref3
          character(80) :: mol_ref4
          integer   :: stoch1
          integer   :: stoch2
          integer   :: stoch3
          integer   :: stoch4
          integer   :: stoch5
          real(kind=8) :: ref  
      end type ener_ref_type     
        
      contains
       
      subroutine fun(x,n,f,group,mpi_master_worker_comm,mpi_worker_comm, mpi_worker_worker_comm, &
                     master_worker_group,input_files,ener_ref)
      use mpi
      use ls_rmsd
           
      implicit none
      integer,intent(in)   :: n
      integer              :: i,j,funct,env_id,ierr,iter_inp,group,my_id
      integer              :: comm,mpirank,master_worker_group,mpisize,mpierror
      integer              :: mpi_master_worker_comm,mpi_worker_comm,mpi_worker_worker_comm
      integer,parameter    :: h2o_id=0,h3o_id=15,ohm_id=21,h2o8s4_id=10,h2o_eaip_id=30
      real(kind=dp) :: f(:),dummy2,energy
      real(kind=dp),allocatable :: res(:),dummy(:),time(:)
      real(kind=dp) :: x(n),a,temp,rnd,h2o_ref,h2o_eaip_ref,h3op_ref,ohm_ref,h2o8s4_ref
      real(kind=dp) :: z(n),M(n,n),tmp
      real(kind=dp),dimension(3,3) :: dipole
      real(kind=dp)           :: dipole_res
      real(kind=dp),parameter :: debye=0.393430307_dp
      real(kind=dp),allocatable :: pos(:),pos1(:)
      real(kind=dp),parameter :: pi=2*asin(1._dp),kcalmol=6.27509468713739E+02_dp,ev=27.2116_dp
      character(30)           :: out_file,out_group
      character(30),dimension(:)    :: input_files
      type(ener_ref_type)  :: ener_ref(:)
      logical,allocatable :: converged(:),log_dummy(:)
      !--------rmsd-----------------
      real(dp),dimension(3,3) :: U
      real(dp),dimension(3) :: center1,center2
      logical :: calc_g = .false.
      real(dp),dimension(:,:),allocatable :: g 
      real(dp)  :: rmsd_error,ipea_error(2),ipea_ref(2) 
 !------------------timing----------------------------------------
      integer :: t1,t2,clock_rate,clock_max
      
      CALL MPI_Comm_size(mpi_comm_world, mpisize, mpierror)
      allocate(dummy(mpisize),log_dummy(mpisize))
      allocate(time(mpisize),converged(mpisize))
      converged=.true.
      log_dummy=.true.
      f=0._dp
      rmsd_error=0._dp
      ipea_error=0._dp
      time=0._dp
!----------------------------initialize fore env-------------------------      
      write(out_group,*) group 
      out_group=adjustl(out_group)
      out_file=trim(out_group)//"out"
      my_id = group+1
!      print *,my_id,input_files(my_id),group
      CALL cp_create_fenv_comm(env_id,trim(input_files(my_id)),out_file,mpi_worker_comm,ierr)
      IF (ierr.NE.0) THEN
        print *,"my bad",group
        STOP "cp_create_fenv"
      ENDIF 
!==========================worker caluculations========================     
      call system_clock ( t1, clock_rate, clock_max )
      
      CALL cp_get_natom(env_id, i, ierr)
      CALL cp_get_nparticle(env_id, j, ierr)
      allocate(pos(3*j))
      allocate(pos1(3*j))
      CALL cp_get_pos(env_id, pos, 3*j, ierr)   
!      print *,my_id,input_files(my_id),group
      if  (INDEX (trim(adjustl(input_files(my_id))), 'geo_opt') .eq. 1) then
          !  geo opt
         call cp_do_geo_opt(env_id,converged(my_id),ierr)
         allocate(g(2,3*j))
         CALL cp_get_pos(env_id, pos1, 3*j, ierr)
         call rmsd(j,reshape(pos,(/ 3, j /)),reshape(pos1,(/ 3, j /) ),0, U,center2,center1, f(my_id),calc_g,g)
         dummy2=-1.0_dp
         if (.not. converged(my_id)) then
             f(my_id) = sqrt(dummy2) 
         endif
      rmsd_error=0._dp
      else if (INDEX (trim(adjustl(input_files(my_id))), 'dipole') .eq. 1) then
        dipole = 0.0_dp
     !   CALL cp_calc_multipole(env_id,dipole,3,ierr)
     !   dipole_res=norm2(dipole(:,1)+dipole(:,2)+dipole(:,3))/debye
        dipole_res = 1.85_dp
      else
        CALL cp_calc_energy(env_id,pos,3*j,f(my_id),ierr,converged(my_id))
!        print *,f(my_id),'test'
        dummy2=-1.0_dp
        if (.not. converged(my_id)) then
            f(my_id) = sqrt(dummy2)
        endif
      endif
      IF (ierr.NE.0) THEN
        print *,"my bad",group
        STOP "cp_calc_energy"
      ENDIF
   
!==========================main iteration end==============================
!      IF (ierr.NE.0) STOP "cp_destroy_fenv"
!      ! TO BE CHANGED !!!!!!!!!!! STOICHMOETRIC CORRECTNESS
!      !-----------------ip/ea error--------
!          !first +1 and second -1 reference in eV
!          ipea_ref(1) = 12.621_dp
!          ipea_ref(2) = 0.16_dp
!         ! ipea_ref(2) = 7.16_dp
!          if (my_id-h2o_eaip_id .gt. 1) then
!          !-------this is also just kinda hack, should be improved------------------
!          res(2) = abs((f(1) - h2o_eaip_ref)*eV - ipea_ref(my_id-h2o_eaip_id-1)) 
!         ! print *,my_id,'ipea',f(1)*eV,h2o_eaip_ref*eV,abs((f(1) - h2o_eaip_ref)*eV - ipea_ref(my_id-h2o_eaip_id)),res(2)
!          endif
!      else if (INDEX (trim(adjustl(input_files(my_id))), 'dipole') .eq. 1) then
!      !--------------dipole error---------
!          !print *,'dipole',dipole_res
!          res(3)=abs(1.85_dp-dipole_res) 
!      else  if  (INDEX (trim(adjustl(input_files(my_id))), 'geo_opt') .eq. 1) then 
!      !--------------rmsd error-----------
!          res(4)=rmsd_error
!=========================================================================
            call MPI_Reduce(f,dummy,mpisize,MPI_DOUBLE_PRECISION, MPI_SUM,0,mpi_master_worker_comm,ierr)
            call MPI_Reduce(time,dummy,mpisize,MPI_DOUBLE_PRECISION, MPI_SUM,0,mpi_master_worker_comm,ierr)
            call MPI_Reduce(converged,log_dummy,mpisize,MPI_LOGICAL, MPI_LAND,0,mpi_master_worker_comm,ierr)
!=========================================================================

      call cp_destroy_fenv(env_id,ierr)
      call system_clock ( t2, clock_rate, clock_max )   
      time(my_id) =  real ( t2 - t1 ) / real ( clock_rate ) 
      deallocate(pos,pos1,time)
      deallocate(dummy,log_dummy,converged)
      end subroutine fun
!--------------------------------------------------------------------------
      subroutine init_cp2k(master)
      use mpi
      logical master 
      IF(.NOT. master) THEN
        call cp_init_cp2k(0,ierr)
        IF (ierr.NE.0) STOP "cp_init_cp2k"
      ENDIF
      end subroutine
       
      subroutine finalize_cp2k(master)
      use mpi
      logical master 
      IF(.NOT. master) THEN
        call cp_finalize_cp2k(0,ierr)
        IF (ierr.NE.0) STOP "cp_init_cp2k"
      ENDIF
      end subroutine
      
      subroutine var_trans(x,y,lbnd,ubnd)
      use random
      real(kind=dp),dimension(:) :: x,lbnd,ubnd
      real(kind=dp),dimension(:),intent(out) :: y
      y = ((x - lbnd) / (ubnd - lbnd))
      end subroutine

      subroutine var_back_trans(x,y,lbnd,ubnd)
      use random
      real(kind=dp),dimension(:) :: x,y,lbnd,ubnd
      y = x*(ubnd-lbnd) + lbnd
      end subroutine

end module
