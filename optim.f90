MODULE optim
    use defs 
  TYPE opt_state_type
    INTEGER            :: state
    INTEGER            :: nvar
    INTEGER            :: iprint
    INTEGER            :: unit
    INTEGER            :: maxfun
    REAL(dp)           :: rhobeg, rhoend
    REAL(dp), DIMENSION(:), POINTER  :: w
    REAL(dp), DIMENSION(:), POINTER  :: xopt
    ! local variables
    INTEGER            :: np, nh, nptm, nftest, idz, itest, nf, nfm, nfmm, &
                          nfsav, knew, kopt, ksave, ktemp
    REAL(dp)           :: rhosq, recip, reciq, fbeg, fopt, diffa, xoptsq, &
                          rho, delta, dsq, dnorm, ratio, temp, tempq, beta, &
                          dx, vquad, diff, diffc, diffb, fsave, detrat, hdiag, &
                          distsq, gisq, gqsq, f, bstep, alpha, dstep
  END TYPE opt_state_type

  TYPE gauss_opt_type
    integer                                  :: n,nf,file_id
    real(kind=dp),dimension(:),allocatable   :: eta,m,x,xmin
    real(kind=dp),dimension(:),allocatable   :: l_bound,u_bound,l_bound_trans,u_bound_trans
    real(kind=dp),dimension(:,:),allocatable :: Q,C
    real(kind=dp)                            :: c_t,f,fmin,r
    real(kind=dp)                            :: wm,wc,wt,beta,p,fe,fc
    character(30)                            :: task

  END TYPE gauss_opt_type

CONTAINS

!-----------------------Gaussian Adaption Algorithm-----------------------------
	subroutine GaussAdapt(n,x,gauss_opt)
    use defs
	
    implicit none
    type(gauss_opt_type),target        :: gauss_opt
	character(30),pointer              :: task
	integer,intent(in)                 :: n
    real(kind=dp),parameter            :: sto_eps=1e-2_dp
    real(kind=dp)                      :: x(:)
	real(kind=dp),pointer              :: f,c_t,r,eta(:),m(:)
	real(kind=dp),pointer              :: Q(:,:),C(:,:)
	real(kind=dp),pointer              :: l_bound(:),u_bound(:)
	real(kind=dp),pointer              :: fmin,xmin(:)
	real(kind=dp)                      :: detQ,dQ(n,n),delta(n)
	real(kind=dp)                      :: Id(n,n),eig_val(n),work(3*n-1),test,temp(n,n),x_old(n)
	real(kind=dp)                      :: wm,wc,wt,beta,p,fe,fc
	integer                            :: i,j,info,algorithm
    integer,pointer                    :: file_id,nf
    logical                            :: u_check(n),l_check(n),file_ex
!----------------------------set pointers to type struct-----------------------
    f => gauss_opt%f
    c_t => gauss_opt%c_t
    r => gauss_opt%r
    eta => gauss_opt%eta
    m => gauss_opt%m
    Q => gauss_opt%Q
    C => gauss_opt%C
    l_bound => gauss_opt%l_bound_trans
    u_bound => gauss_opt%u_bound_trans
    fmin => gauss_opt%fmin
    xmin => gauss_opt%xmin
    task => gauss_opt%task
    file_id => gauss_opt%file_id
    nf => gauss_opt%nf
!----------------------------set GauAdapt prameter------------------------------
    print *,'in gauss adapt    ',trim(task),r,c_t
    algorithm=1
    x_old=x
	Id=0._dp
    do i=1,n; Id(i,i)=1._dp; enddo
    wm=exp(1._dp)*n
!    wm=1._dp
	wc=(n+1._dp)**2/(log(n+1._dp))
!	wt=exp(1._dp)*n
    wt=wc/2._dp
!	wt=10._dp*n**2
!    wt=1._dp
	beta=1._dp/wc
	p=1._dp/exp(1._dp)
    fe=1._dp + beta*(1._dp-p)
     fc=1._dp-beta*p
!-------------------------------------------------------------------------------
	if (task(1:5).eq.'new_x') then
		if (f .lt. c_t) then
            delta = x - m
			task='new_x_accepted'
!----------------------------save current best----------------------------------
			if (f .lt. fmin) then
			fmin = f 
			xmin = x 
			endif
!------------------------------adapt parameters---------------------------------
			r=fe*r
			c_T = (1._dp-1._dp/wt)*c_T+f/wt
			m = (1._dp-1._dp/wm)*m + x/wm
!-----------------------------update convariance--------------------------------

!-----------------------------rank 1 update ------------------------------------
     select case (algorithm)
        case(1)
            dQ=(1._dp-1._dp/wc)*Id
            call DSYR('l',n,1._dp/wc,eta,int(1),dQ,n)
        case(2)
            dQ=(1._dp-1._dp/wc)*C
!            call DSYR('l',n,1._dp/wc,delta,int(1),dQ,n)
			call DGER(n,n,1._dp/wc,delta,int(1),delta,int(1),dQ,n)
            C=dQ
     end select
!-----------------------------eig.val decomposition-----------------------------
			call DSYEV( 'V', 'l', n, dQ, n, eig_val, WORK, 3*n-1, INFO )
			if (info.ne.0) task='stop_eigendecomp'
            if (any(isnan(eig_val))) task='stop'
            if (any(eig_val.lt.0)) then
                task='stop eig_val'
                return 
            endif
!-----------------------------normalization of Q and C--------------------------
	    detQ=product(eig_val)
            detQ=detQ**(1.d0/n)
            C=C/detQ
            eig_val=eig_val/detQ
            eig_val=sqrt(eig_val)
!-----------------------------update Q such that Q=Q*D^1/2------------------------
      select case (algorithm)
         case(1)
            do i=1,n; dQ(:,i)=dQ(:,i)*sqrt(eig_val(i)); enddo
            call dgemm('n','t',n,n,n,1._dp,dQ,n,dQ,n,0._dp,temp,n)
            dQ=temp
            call dsymm('r','l',n,n,1._dp,dQ,n,Q,n,0._dp,temp,n)
            Q=temp
            call dgemm('n','t',n,n,n,1._dp,Q,n,Q,n,0._dp,temp,n)
            C=temp
          case(2)
            do i=1,n; dQ(:,i)=dQ(:,i)*(eig_val(i)); enddo
            Q=dQ
!            call dgemm('n','t',n,n,n,1._dp,Q,n,Q,n,0._dp,temp,n)
!            C=temp
!            print *,temp
!            print *,'---------------------'
!            print *,C
!            if(any(temp.ne.C)) print *,'something wroing'
      end select
!--------------------------------update ---------------------------------------
		else if (isnan(f)) then
            task='new_x_rejected'
            r=fc*r
        else
!---------------rejected lower step size dont adopt cov. + mean-----------------
			task='new_x_rejected'
			r=fc*r
		endif
!------------------------------------sample-------------------------------------
		if (r.lt.1e-9) then
            task='stop_r'
            return
        endif
            do i=1,n; eta(i)=random_normal(); enddo
            !CALL reset_to_next_rng_substream(rng_stream)
            x=m
            call DGEMV('N',n,n,r,Q,n,eta,int(1),1._dp,x,int(1))
            
            if (any(isnan(x))) then 
            task='stop_sampling x'
            ! write nex x to se_para_file
            endif

            do i=1,n
                if (x(i) .lt. l_bound(i)) then
                x(i) = l_bound(i)
                endif
                if (x(i) .gt. u_bound(i)) then
                x(i) = u_bound(i)
                endif        
            enddo
            !-------------end boundary check-----------------
                        
	endif
	
	if (task(1:5).eq.'start') then
     print *,'This is the gauss start - yeah'
     inquire(file="gauss.restart",exist=file_ex)
     if (file_ex) then
        print *,'Gauss parameters from Restart file'
        ! restart file should be read in main!!!!!
        !call read_gauss_file(eta,m,C,Q,r,c_t,xmin,fmin,file_id)
     else    
		Q = Id 
        C = Q
        x = x_old
		m = x
		task='new_x'
		c_t=2e4
        r=0.1_dp
		fmin=c_t
		xmin=x
        file_id=300
        nf = 0
      endif
	endif

END SUBROUTINE



    FUNCTION random_normal() RESULT(fn_val)
    
    integer,parameter :: dp=8
    REAL(kind=8) :: fn_val,half=0.5_dp
    
    REAL(kind=8)     :: s = 0.449871, t = -0.386595, a = 0.19600, b = 0.25472,           &
                r1 = 0.27597, r2 = 0.27846, u, v, x, y, q
    
    !     Generate P = (u,v) uniform in rectangle enclosing acceptance region
    DO
      CALL RANDOM_NUMBER(u)
      CALL RANDOM_NUMBER(v)
      v = 1.7156 * (v - half)
    
    !     Evaluate the quadratic form
      x = u - s
      y = ABS(v) - t
      q = x**2 + y*(a*y - b*x)
    
    !     Accept P if inside inner ellipse
      IF (q < r1) EXIT
    !     Reject P if outside outer ellipse
      IF (q > r2) CYCLE
    !     Reject P if outside acceptance region
      IF (v**2 < -4.0*LOG(u)*u**2) EXIT
    END DO
    
    !     Return ratio of P's coordinates as the normal deviate
    fn_val = v/u
    RETURN
    END FUNCTION random_normal

END MODULE

