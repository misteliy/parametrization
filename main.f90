PROGRAM main
    USE defs
    USE mpi
    USE obj_fun
    USE optim
    USE write_nc
    IMPLICIT NONE
    INTEGER :: mpierror, mpisize, mpirank 
    INTEGER :: mpi_top_level_comm, mpi_worker_comm, mpi_master_worker_comm
    INTEGER :: mpi_worker_worker_comm
    INTEGER :: group, group_rank, master_worker_rank=-1,ierr
    INTEGER :: master_worker_group, nr_energy_groups, nr_geom_groups
    INTEGER :: ene_group_size, geo_opt_group_size,ipea_group_size,dipole_group_size, &
               tmp_rank,tmp,env_id,my_id
    INTEGER :: test,i,file_id,file_id2,file_id4
    LOGICAL :: master=.FALSE., incomplete_group =.FALSE.
    logical,allocatable :: ener_converged(:),log_dummy(:)
    integer,parameter :: n=16,num_inp=30
    real(kind=dp),allocatable :: f(:),res_f(:),res(:),dummy(:)
    real(kind=dp) :: x(n),rand
    integer, parameter :: maxiter=1e6
    integer :: j,k,irandom,t,irand
    !--------------------------------------------------------------
    real(kind=dp),parameter :: pi=2*asin(1._dp),kcalmol=6.27509468713739E+02_dp,ev=27.2116_dp,ev_kcal=23.061
    INTEGER :: id,stoch1,stoch2,stoch3,stoch4,stoch5,s
    character(15) :: mol,mol_ref1,mol_ref2,mol_ref3,mol_ref4
    character(30) :: str,out_file,out_group
    real(kind=dp) :: ref,stepr,energy
    integer,parameter    :: h2o_id=1,h3o_id=16,ohm_id=22,h2o8s4_id=11,h2o_eaip_id=31 ! Id of the array not MPI rank
    real(kind=dp) :: h2o8s4_ref,h2o_eaip_ref,h2o_ref,h3op_ref,ohm_ref
    real(kind=dp),dimension(2),parameter :: ipea_ref = (/ 12.530_dp , 5.36_dp /)
    !--------------------------------------------------------------
    real(kind=dp),parameter  :: sampler=1._dp/(2._dp*sqrt(real(n))) 
    real(kind=dp) :: sampx(n),sr,tmpx(n),tmpx2(n)
    logical :: generate,restart
    !------------gaussian adapt-------------------------------------
    type(gauss_opt_type),pointer  :: gauss_opt
    character(30),pointer   :: task
    real(kind=dp),pointer   :: c_t,eta(:),m(:),Q(:,:),C(:,:),r
    real(kind=dp),pointer   :: l_bound(:),u_bound(:),xmin(:),fmin
    integer,pointer         :: file_id3,nf
    real(kind=dp),pointer   :: l_bound_trans(:),u_bound_trans(:)
    real(kind=dp),pointer   :: obj_f
    real(kind=dp)           :: random
    !----------------------------------------------------------------
    integer,pointer         :: opt_state,maxfun,iprint,nvar
    real(kind=dp),pointer   :: rhobeg,rhoend
    !----------------------------------------------------------------
    character(30),dimension(num_inp+5) :: input_files
    type(ener_ref_type),dimension(num_inp) :: ener_ref
    INTEGER :: total_cores,cores_ene,cores_geo_opt,cores_ipea,cores_dipole
    INTEGER :: nr_geo_opt_groups,nr_ene_groups,nr_ipea_groups,nr_dipole_groups
    !------------------timing----------------------------------------
    integer :: t1,t2,clock_rate,clock_max
    !----------------------------------------------------------------
    integer :: npt,np,ixn,ixo
    type(opt_state_type),pointer :: powell_opt
    real(kind=dp) :: rnum(n)
    logical :: file_ex
!===============================================================================
    !----------------------chose algirthm----------------------------
    !allocate(gauss_opt)
    allocate(powell_opt)    
    !----------------set pointers gauss adapt------------------------
    if (associated(gauss_opt)) then
    allocate(gauss_opt%eta(n),gauss_opt%m(n),gauss_opt%xmin(n))
    allocate(gauss_opt%l_bound(n),gauss_opt%u_bound(n))
    allocate(gauss_opt%l_bound_trans(n),gauss_opt%u_bound_trans(n))
    allocate(gauss_opt%x(n),gauss_opt%Q(n,n),gauss_opt%C(n,n))
    task => gauss_opt%task
    c_t => gauss_opt%c_t
    r => gauss_opt%r
    eta => gauss_opt%eta
    m => gauss_opt%m
    Q => gauss_opt%Q
    C => gauss_opt%C
    l_bound => gauss_opt%l_bound
    l_bound_trans => gauss_opt%l_bound_trans
    u_bound => gauss_opt%u_bound
    u_bound_trans => gauss_opt%u_bound_trans
    fmin => gauss_opt%fmin
    xmin => gauss_opt%xmin
    obj_f => gauss_opt%f
    file_id3 => gauss_opt%file_id
    file_id3 = 300
    nf => gauss_opt%nf
    endif
    !---------------set pointers powell opt--------------------------
    if (associated(powell_opt)) then
    np=n+1
    npt=2*n+1
    allocate(powell_opt%w((npt+13)*(npt+n)+3*n*(n+3)/2))
    allocate(powell_opt%xopt(n))
    allocate(l_bound(n),u_bound(n),l_bound_trans(n),u_bound_trans(n))
    allocate(task,file_id3)
    opt_state => powell_opt%state
    maxfun => powell_opt%maxfun
    iprint => powell_opt%iprint
    nvar => powell_opt%nvar
    rhobeg => powell_opt%rhobeg
    rhoend => powell_opt%rhoend
    obj_f => powell_opt%f
    fmin => powell_opt%fopt
    xmin => powell_opt%xopt
    nf => powell_opt%nf
    powell_opt%unit = 6
    file_id3 = 300
    endif
    !-------------------mpi stuff------------------------------------
    master_worker_group = MPI_UNDEFINED
    group = MPI_UNDEFINED
    group_rank = MPI_UNDEFINED
    mpi_worker_comm = MPI_UNDEFINED
    mpierror = 0
    CALL MPI_Init(mpierror)
    CALL MPI_COMM_DUP(MPI_COMM_WORLD, mpi_top_level_comm, mpierror)
    CALL MPI_Comm_size(mpi_top_level_comm, mpisize, mpierror)
    ! -- needs at least one worker (ant the master)
    IF(mpisize.LE.1) THEN
       WRITE(*,*)"error: minimum 2 cores needed, one master AND one worker."
       !CALL write_help()
       CALL MPI_Finalize(mpierror)
       !RETURN
    END IF   
    ! get individual comunication rank
    CALL MPI_Comm_rank(mpi_top_level_comm, mpirank, mpierror)
    IF(mpirank.EQ.mpisize-1) master = .TRUE. ! set master node
!-----------------------------------------------------------------------
   allocate(f(mpisize),res_f(mpisize),res(mpisize),dummy(mpisize))
   allocate(ener_converged(mpisize),log_dummy(mpisize))
   f=0.0_dp    
   res_f=0._dp 
   res=0._dp
   ener_converged=.true.
   log_dummy=.true. 
   dummy=0._dp
!---------------------input file(s)--------------------------------------
      open(88,file='WATER27structures/mol_list.txt')
      do i=1,num_inp
          read(88,*) str
          input_files(i)= trim(str) // ".inp"
      enddo
      close(88)
      input_files(num_inp+1) = 'eaip_H2O.inp'
      input_files(num_inp+2) = 'eaip_H2Op.inp'
      input_files(num_inp+3) = 'eaip_H2On.inp'
      input_files(num_inp+4) = 'dipole.inp'
      input_files(num_inp+5) = 'geo_opt.inp'
    !----only master needs ref data------------------
    if (master) then
      print *,'reading reference data...'
      open(99,file='WATER27.dat')
      do i=1,num_inp
          read(99,*) id,mol,mol_ref1,mol_ref2,mol_ref3,mol_ref4, &
          stoch1,stoch2,stoch3,stoch4,stoch5,ref
          ener_ref(i)%id = id
          ener_ref(i)%mol = trim(mol)
          ener_ref(i)%mol_ref1 = trim(mol_ref1)
          ener_ref(i)%mol_ref2 = trim(mol_ref2)
          ener_ref(i)%mol_ref3 = trim(mol_ref3)
          ener_ref(i)%mol_ref4 = trim(mol_ref4)
          ener_ref(i)%stoch1 = stoch1
          ener_ref(i)%stoch2 = stoch2
          ener_ref(i)%stoch3 = stoch3
          ener_ref(i)%stoch4 = stoch4
          ener_ref(i)%stoch5 = stoch5
          ener_ref(i)%ref = ref
          !if (master) print *,ener_ref(i)
      enddo
      close(99)
      print *,'end reading reference data'
    endif
!--------------------------------------------------------------------------
! this probably needs to be adjusted 
! here the simplest case - can be omitted since it's simpler as I thought :)
!---------------------------------------------------------------------------
   ! initialization
   nr_ene_groups=num_inp
   nr_ipea_groups=3
   nr_geo_opt_groups=1
   nr_dipole_groups=1
   ! IMPLEMENT SOME KIND OF ITERATIVE DISTRIBUTION OF RESOURCES WITH TARGET INPUT
   cores_ene = num_inp
   ene_group_size = 1
   cores_geo_opt = 1
   geo_opt_group_size = 1
   cores_ipea = 3
   ipea_group_size = 1
   cores_dipole = 1
   dipole_group_size = 1
   total_cores = cores_geo_opt + cores_ene + cores_ipea + cores_dipole
   ! check input and set cores and group size
   if (nr_ene_groups.le.0) then
       print *,"energy group at least size 1 please"
       CALL MPI_Finalize(mpierror)
       STOP
   endif 
   IF (total_cores.gt.mpisize-1) THEN
          IF(master) WRITE(*,*)"not enough resources for calculations"
          CALL MPI_Finalize(mpierror)
          STOP
   ENDIF
   ! --------------------------------------------------------------------------
   ! dividing availible cores into certain groups
   IF(master) THEN ! master 
      master_worker_group = 1                      ! belong to master_worker_comm
      master_worker_rank = 0                                ! rank in m_w_comm
      group_rank = 0                                        ! rank in worker group
   ELSE ! workers
      ! energy calculation groups
      IF(mpirank .LT. cores_ene) THEN
         group = INT(mpirank/ene_group_size)                  ! assign to groups
         ! master of worker group
         tmp = MODULO(mpirank,ene_group_size)
         IF(tmp.EQ.0) THEN           ! group masters
            master_worker_group = 1             ! belong to master_worker_comm
            master_worker_rank = group +1          ! rank in m_w_comm
         END IF
         group_rank = tmp          ! rank in worker group
      ! configurational change groups
      ELSE IF(mpirank .LT.(cores_ene + cores_ipea)) THEN
            !group = nr_ene_groups + INT((tmp_rank)/cores_geo_opt)+1 !group in geo opt calculation groups
            group = INT((mpirank-cores_ene)/ipea_group_size) + nr_ene_groups
            ! master of worker group
            tmp = MODULO(mpirank,ipea_group_size)
            IF(tmp.EQ.0)THEN ! group masters
               master_worker_group = 1                      ! belong to master_nmc_worker_comm
               master_worker_rank = group +1   ! rank in m_w_comm
            END IF
            group_rank = tmp 
      ELSE IF(mpirank .LT.(cores_ene + cores_ipea + cores_dipole)) THEN
            !group = nr_ene_groups + INT((tmp_rank)/cores_geo_opt)+1 !group in geo opt calculation groups
            group = INT((mpirank-(cores_ene+cores_ipea))/dipole_group_size) + nr_ene_groups + nr_ipea_groups
            ! master of worker group
            tmp = MODULO(mpirank,ipea_group_size)
            IF(tmp.EQ.0)THEN ! group masters
               master_worker_group = 1                      ! belong to master_nmc_worker_comm
               master_worker_rank = group +1   ! rank in m_w_comm
            END IF
            group_rank = tmp 
      ELSE IF(mpirank .LT.(total_cores)) THEN
            !group = nr_ene_groups + INT((tmp_rank)/cores_geo_opt)+1 !group in geo opt calculation groups
            group = INT((mpirank-(cores_ene+cores_ipea+cores_dipole))/geo_opt_group_size) + nr_ene_groups  &
            + nr_ipea_groups + nr_dipole_groups
            tmp = MODULO(mpirank,geo_opt_group_size)
            IF(tmp.EQ.0)THEN ! group masters
               master_worker_group = 1                      ! belong to master_nmc_worker_comm
               master_worker_rank = group +1   ! rank in m_w_comm
            END IF
            group_rank = tmp 
      ELSE
         !print *,'yeah it is me',mpirank
         ! not used cores
         incomplete_group = .TRUE.
      END IF
   !   print *,mpirank,group,group_rank,master_worker_group,master_worker_rank,input_files(mpirank+1)
   END IF

   ! -- spiltting communicators
   ! worker intern communication
   CALL MPI_COMM_SPLIT(mpi_top_level_comm, group, group_rank, mpi_worker_comm, mpierror)
   ! worker master communication
   CALL MPI_COMM_SPLIT(mpi_top_level_comm, master_worker_group, &
            master_worker_rank, mpi_master_worker_comm, mpierror)
   ! worker - worker communcation
   IF(master) master_worker_group = MPI_UNDEFINED
   CALL MPI_COMM_SPLIT(mpi_top_level_comm, master_worker_group, &
            master_worker_rank, mpi_worker_worker_comm, mpierror)
   ! wait till all processes are up and running
   CALL MPI_BARRIER(mpi_top_level_comm, mpierror)
   
   ! --------------------------------------------------------------------------
   ! START PROGRAM
   ! --------------------------------------------------------------------------
  
   IF ( (.NOT. master) .and. group.ge.0) call init_cp2k(master)
   IF (master) THEN
!              l_bound=(/-2._dp,0.05_dp,1._dp,-1.5_dp,0._dp,0._dp,-50._dp,&
!        -30._dp,0.05_dp,0.05_dp,-3._dp,0.05_dp,-5._dp,0._dp,0._dp,-50._dp/)
!        u_bound=(/-0.5_dp,0.4_dp,1.5_dp,-0.5_dp,10._dp,300._dp,100._dp,&
!        -20._dp,0.4_dp,3._dp,-2._dp,3._dp,-2._dp,10._dp,300._dp,100._dp/)
           l_bound=(/-2._dp,0.05_dp,1._dp,-1.5_dp,0._dp,0._dp,-50._dp,&
        -40._dp,0.01_dp,0.05_dp,-3._dp,0.05_dp,-5._dp,0._dp,0._dp,-50._dp/)
        u_bound=(/-0.5_dp,1.0_dp,1.5_dp,-0.5_dp,10._dp,300._dp,300._dp,&
        -20._dp,0.4_dp,3._dp,-2._dp,5._dp,-2._dp,10._dp,400._dp,300._dp/)
    l_bound_trans=0._dp
    u_bound_trans=1._dp
    call read4file(n,x,restart)
    call para2file2(n,x,"para.xml",2222)
    file_id = 100
    file_id2 = 200
    file_id4 = 400
    call open_file("test.out",file_id)
    call open_file("test_scaled.out",file_id2)
   !--------------restart files  -----------------------------------------------
   if (associated(gauss_opt)) then
         inquire(file="gauss.restart",exist=file_ex)
           if (file_ex) then
              print *,'Gauss parameters from Restart file'
              call read_gauss_file(eta,m,C,Q,r,c_t,xmin,fmin,nf,file_id3)
              task="new_x"
           else
              task="start"
           endif
   else if (associated(powell_opt)) then
         inquire(file="powell.restart",exist=file_ex) 
         if (file_ex) then
             print *,'Powell parameters from Restart file'
             call read_powell_file(powell_opt)
             task='new_x'
             if (powell_opt%state .eq. -1) opt_state = 0
             if (powell_opt%state .eq. 7) opt_state = 0
         else
            print *,'Powell new run :) '
            opt_state = 0
            maxfun = maxiter
            rhobeg = 1._dp
            rhoend = 0.001_dp
            iprint = 1
            nvar = n
            task='start'
         endif
   endif
   endif
   ! --------------Iteration heeeerre-------------------------------------------
   do i=1,maxiter
        call MPI_BCAST(task,20,MPI_CHARACTER,mpisize-1,mpi_top_level_comm,ierr)
        call MPI_BCAST(x,n,MPI_DOUBLE_PRECISION,mpisize-1,mpi_top_level_comm,ierr)
        IF (master) THEN
            res=0._dp
            if (task(1:4).eq.'stop') STOP !and terminate mpi...balblalba
            if (task(1:4).eq.'skip') res(1) = 1e9_dp 
            if (task(1:5).eq.'new_x') then
                call MPI_Reduce(dummy,res_f,mpisize,MPI_DOUBLE_PRECISION, MPI_SUM,0,mpi_master_worker_comm,ierr)
                call MPI_Reduce(log_dummy,ener_converged,mpisize,MPI_LOGICAL, MPI_LAND,0,mpi_master_worker_comm,ierr)
            !---------------energy references-------------------------------------------
            h2o_ref = res_f(h2o_id)
            h2o_eaip_ref = res_f(h2o_eaip_id)
            h3op_ref = res_f(h3o_id)
            ohm_ref = res_f(ohm_id) 
            h2o8s4_ref = res_f(h2o8s4_id)
            !-----------------calculate obj_fun--------------------
            write(*,*) '           !!! energy in kcalmol'
            write(*,*) '           mol ','s1 ','s2 ','s3 ','s4 ','s5 ','ener       ','ener diff ','ref     ','error   ','converged'
                    !---calcualte error difference----------------
                    do s=1,30
                          energy = (ener_ref(s)%stoch1*res_f(s) + &
                          ener_ref(s)%stoch2*h3op_ref + &
                          ener_ref(s)%stoch3*ohm_ref + &
                          ener_ref(s)%stoch4*h2o_ref + & 
                          ener_ref(s)%stoch5*h2o8s4_ref)*kcalmol
                          res(s) = abs(energy-ener_ref(s)%ref)
                          write(*,'(A15,1X,I2,1X,I1,1X,I1,1X,I2,1X,I2,F10.2,F10.2,F10.2,1X,F10.2,1X,L1)') &
                          trim(ener_ref(s)%mol),ener_ref(s)%stoch1,ener_ref(s)%stoch2, &
                          ener_ref(s)%stoch3,ener_ref(s)%stoch4,ener_ref(s)%stoch5,res_f(s)*kcalmol, &
                          energy,ener_ref(s)%ref,res(s),ener_converged(s)
                    end do
                    !print *,'res',res
                    !---calculate ea/ip---------------------------- 
                    print *,'h2o  ',res_f(31)
                    res(32) = (abs(res_f(32)-h2o_eaip_ref)*eV-ipea_ref(2))*ev_kcal
                    print *,'h2o_n',res_f(32),res(32),'kcalmol'
                    res(33) = (abs(res_f(33)-h2o_eaip_ref)*eV-ipea_ref(1))*ev_kcal
                    print *,'h2o_p',res_f(33),res(33),'kcalmol'
                    !---calculate rmsd---------------------------- 
                    !---calculate dipole--------------------------
            endif
            !------------------------------------------------------
            call system_clock ( t2, clock_rate, clock_max )            
            print *,'-----------------------------'
            print *,'iteration =',i,'Elapsed real time = ', real ( t2 - t1 ) / real ( clock_rate )
            obj_f=sum(res(1:30))
            print *,'-----------------------------'
            print *, 'task  ',trim(task),'  obj_f',obj_f,'geo rmsd',res_f(35),'nf',nf,'fmin',fmin,'state',opt_state
            print *,'-----------------------------'
!            print *,'error in:            energy,   ip/ea,    dipole,      rmsd   '
!            print *,'function =',f(1),f(2),f(3),f(4)
            print *,'-----------------------------'
            call system_clock ( t1, clock_rate, clock_max )
            tmpx = 0.0_dp
            call var_trans(x,tmpx,l_bound,u_bound)
            !-------------optimization algorithm-------------------
            if (associated(powell_opt)) call cp_do_powell_optimize(n,tmpx,powell_opt)
            !------------gauss adapt in normalized variables-------
            if (associated(gauss_opt)) call GaussAdapt(n,tmpx,gauss_opt)
            !-------------------------------------------------------
            !call sample(n,x,l_bound,u_bound)
            !-------------------------------------------------------
            call var_back_trans(tmpx,x,l_bound,u_bound)
            !-----------simple check for boundaries ----------------
            !--------should be replaced by bobyqa------------------- 
            task = "new_x"
            do k=1,n
                if (x(k) .lt. l_bound(k)) then
                x(k) = l_bound(k)
                task(1:4) = 'skip'
                endif
                if (x(k) .gt. u_bound(k)) then
                x(k) = u_bound(k)
                task(1:4) = 'skip'
                endif        
            enddo
            !-------------------------------------------------------
            if ( (i .gt. 1) .and. (.not. isnan(obj_f)) ) then
            call write_file(x,obj_f,file_id)
            call write_file(tmpx,obj_f,file_id2)
            if (associated(gauss_opt)) call write_gauss_file(eta,m,C,Q,r,c_t,xmin,fmin,nf,file_id3)
            if (associated(powell_opt)) call write_powell_file(powell_opt,file_id3)
            call write_optim_file(xmin,fmin,file_id4)
            endif
            if (task(1:4) .ne. 'skip') call para2file2(n,x,"para.xml",2222)
            !-------------------------------------------
        ELSE IF (group.ge.0)  THEN
            ! call fun_eval_routine
            IF (task(1:5).eq.'new_x') THEN
            call fun(x,n,f,group,mpi_master_worker_comm,mpi_worker_comm, & 
                     mpi_worker_worker_comm,master_worker_group,input_files,ener_ref)
            ENDIF
        ENDIF
   enddo
   ! --------------------------------------------------------------------------
   IF ( (.NOT. master) .and. group.ge.0) call finalize_cp2k(master) 
   ! --------------------------------------------------------------------------
   ! -- deallocate stuff (initialized in tmc_command_line_opt)
   ! CALL finalize_mv_types()
   ! wait till all processes are here and quit
   !CALL MPI_BARRIER(mpi_top_level_comm, mpierror)
 
   CALL MPI_Finalize(mpierror)
!write(6,*)"node ",mpirank, "finalized and ready to stop."
   IF(master) call close_file(file_id) 
   IF(master) call close_file(file_id2) 
   IF(master) WRITE(6,*)"all work done."
   IF(master) WRITE(6,*) "#geogroups",nr_geo_opt_groups
   IF(master) WRITE(6,*) "#enegroups",nr_ene_groups
   IF(master) WRITE(6,*) "ene groups size",ene_group_size
   IF(master) WRITE(6,*) "geo opt groups size",geo_opt_group_size
   IF(master) WRITE(6,*) "cores_ene",cores_ene
   IF(master) WRITE(6,*) "cores_geo_opt",cores_geo_opt
   IF(master) WRITE(6,*) "incomplete group",incomplete_group
FLUSH(6)
END PROGRAM

subroutine sample(n,x,l_bound,u_bound)
    use obj_fun
    !--------------------------------------------------------------
    integer,intent(in) :: n
    integer,parameter :: dp=8
    integer :: irand,t
    real(kind=dp) :: x(n)
    real(kind=dp) :: sampler,rnum,step,temp
    real(kind=dp) :: sampx(n),sr,tmpx(n),tmpx2(n),l_bound(n),u_bound(n)
    logical :: generate
    
    CALL init_random_seed() 
    sampler=1._dp/(2._dp*sqrt(real(n)))
    generate=.TRUE.
    step=1._dp/2._dp
    t=1
    irand=7
            do while (generate .and. t .lt. 10)
!--------------------------random only one variable---------------------------
                print *,'generation loop',t
                do while (irand .eq. 7)
                 call random_number(rnum)
                 irand=int(rnum*n)+1
                enddo
                print *,'change var',irand
                call var_trans(x,tmpx,l_bound,u_bound)
                print *,'after var trans',tmpx(irand)
                call random_number(rnum)
                if (rnum .gt. 0.5_dp) then
                    step=step
                else
                    step=-step
                endif
                tmpx(irand) = tmpx(irand) + step
                print *,'before check',tmpx(irand) 
                if (tmpx(irand) .gt. 1.0_dp) then 
                    tmpx(irand) = tmpx(irand) - 2*step
                endif
                
                if (tmpx(irand) .lt. 0._dp) then
                    tmpx(irand) = tmpx(irand) + 2*abs(step)
                endif
               
                print *,'test',tmpx(irand)
                print *,'old var',x(irand), &
                'new var',(tmpx(irand)*((u_bound(irand)) - (l_bound(irand)))) + (l_bound(irand))
!--------------------------random ball walk-------------------------------------  
!                call random_number(sampx)   
!                call var_trans(x,tmpx,l_bound,u_bound)
!                sampx = sampx-tmpx
!                sampx(7) = 0._dp
!                sr = sqrt(sum(sampx**2))
!                stepr = 1._dp/(2._dp*sqrt(real(n)))
!                sampx = stepr*sampx/sr
!                tmpx = tmpx+sampx
!--------------------------random sample---------------------------------------
!                call random_number(sampx)
!                tmpx=sampx
!                tmpx(7)=1._dp
!              !  if (tmpx(4) .gt. tmpx(5)) then 
!              !      temp=tmpx(4)
!              !      tmpx(4)=tmpx(5)
!              !      tmpx(5)=temp 
!              !  endif
!              !  if (tmpx(21) .gt. tmpx(22)) then 
!              !      temp=tmpx(21)
!              !      tmpx(21)=tmpx(22)
!              !      tmpx(22)=temp 
!              !  endif
!--------------------------change one random variable---------------------------
!                print *,'generation loop',t
!                do while (irand .eq. 7)
!                 call random_number(rnum)
!                 irand=int(rnum*n)+1
!                enddo
!                call random_number(rnum)
!                call var_trans(x,tmpx,l_bound,u_bound)
!                tmpx(irand) = rnum
!                print *,'change var',irand
!                print *,'test',tmpx(irand)
!                print *,'old var',x(irand), &
!                'new var',(tmpx(irand)*((u_bound(irand)) - (l_bound(irand)))) + (l_bound(irand))
!-------------------------------------------------------------------------------
                call var_back_trans(tmpx,tmpx2,l_bound,u_bound)
                tmpx2(7) = u_bound(7)
                

                do j=1,n
                    if (tmpx2(j) .lt. l_bound(j)) then
                        generate = .TRUE. 
                        print *,j,tmpx2(j)
                        exit
                    else
                        generate = .FALSE.
                    endif
                    if (tmpx2(j) .gt. u_bound(j)) then
                        generate = .TRUE.
                        print *,j,tmpx2(j)
                        exit
                    else
                        generate = .FALSE.
                    endif   
                    !    call random_number(random)
                    !    j = random*n
                    !    if (j .eq. 0) j = 1
                    !    call random_number(random)
                    !    x(j) = random*abs(u_bound(j)-l_bound(j)) + l_bound(j)
                enddo
                t = t +1
            enddo
            x = tmpx2

end subroutine

SUBROUTINE init_random_seed()
  INTEGER :: i, n, clock
  INTEGER, DIMENSION(:), ALLOCATABLE :: seed

  CALL RANDOM_SEED(size = n)
  ALLOCATE(seed(n))

  CALL SYSTEM_CLOCK(COUNT=clock)

  seed = clock + 37 * (/ (i - 1, i = 1, n) /)
  CALL RANDOM_SEED(PUT = seed)

  DEALLOCATE(seed)
END SUBROUTINE
