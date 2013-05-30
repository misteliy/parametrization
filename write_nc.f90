module write_nc

   use optim

contains

subroutine open_file(file_name,file_id)
    character (len= *), intent(in) :: file_name
    integer,intent(out) :: file_id
    OPEN ( unit = file_id, file=file_name, access=trim("append"))
end subroutine

subroutine write_file(x,f,file_id)
    integer,intent(in) :: file_id
    real(kind=8),DIMENSION(:)   :: x
    real(kind=8)   :: f !(:)
    write(file_id,*) x,f
    flush(file_id)
end subroutine

subroutine write_gauss_file(eta,m,C,Q,r,c_t,xmin,fmin,nf,file_id)
    integer,intent(in) :: file_id,nf
    real(kind=8)   :: eta(:),m(:),C(:,:),Q(:,:),r,c_t,xmin(:),fmin
    open(unit=file_id,file='gauss.restart')
    write(file_id,*) eta
    write(file_id,*) m
    write(file_id,*) C
    write(file_id,*) Q
    write(file_id,*) r
    write(file_id,*) c_t
    write(file_id,*) xmin,fmin
    write(file_id,*) nf
    flush(file_id)
    close(file_id)
end subroutine
subroutine read_gauss_file(eta,m,C,Q,r,c_t,xmin,fmin,nf,file_id)
    integer :: file_id,nf
    real(kind=8)   :: eta(:),m(:),C(:,:),Q(:,:),r,c_t,xmin(:),fmin
    open(unit=file_id, file='gauss.restart')
    read(file_id,*) eta
    read(file_id,*) m
    read(file_id,*) C
    read(file_id,*) Q
    read(file_id,*) r
    read(file_id,*) c_t
    read(file_id,*) xmin,fmin
    read(file_id,*) nf
    close(file_id)
end subroutine
subroutine write_powell_file(opt,file_id)
    integer :: file_id
    type(opt_state_type) :: opt
    open(unit=file_id,file='powell.restart')
    write(file_id,*) opt%state 
    write(file_id,*) opt%xopt 
    write(file_id,*) opt%maxfun
    write(file_id,*) opt%rhobeg
    write(file_id,*) opt%rhoend
    write(file_id,*) opt%w     
    write(file_id,*) opt%nvar
    write(file_id,*) opt%iprint 
    write(file_id,*) opt%np       
    write(file_id,*) opt%nh       
    write(file_id,*) opt%nptm     
    write(file_id,*) opt%nftest   
    write(file_id,*) opt%idz      
    write(file_id,*) opt%itest    
    write(file_id,*) opt%nf       
    write(file_id,*) opt%nfm      
    write(file_id,*) opt%nfmm     
    write(file_id,*) opt%nfsav    
    write(file_id,*) opt%knew     
    write(file_id,*) opt%kopt     
    write(file_id,*) opt%ksave    
    write(file_id,*) opt%ktemp    
    write(file_id,*) opt%rhosq    
    write(file_id,*) opt%recip    
    write(file_id,*) opt%reciq    
    write(file_id,*) opt%fbeg     
    write(file_id,*) opt%fopt     
    write(file_id,*) opt%diffa    
    write(file_id,*) opt%xoptsq   
    write(file_id,*) opt%rho      
    write(file_id,*) opt%delta    
    write(file_id,*) opt%dsq      
    write(file_id,*) opt%dnorm    
    write(file_id,*) opt%ratio    
    write(file_id,*) opt%temp     
    write(file_id,*) opt%tempq    
    write(file_id,*) opt%beta     
    write(file_id,*) opt%dx       
    write(file_id,*) opt%vquad    
    write(file_id,*) opt%diff     
    write(file_id,*) opt%diffc    
    write(file_id,*) opt%diffb    
    write(file_id,*) opt%fsave    
    write(file_id,*) opt%detrat   
    write(file_id,*) opt%hdiag    
    write(file_id,*) opt%distsq   
    write(file_id,*) opt%gisq     
    write(file_id,*) opt%gqsq     
    write(file_id,*) opt%f        
    write(file_id,*) opt%bstep    
    write(file_id,*) opt%alpha    
    write(file_id,*) opt%dstep    
    flush(file_id)
    close(file_id)
end subroutine
subroutine read_powell_file(opt)
    integer :: file_id=439857
    type(opt_state_type) :: opt
    open(unit=file_id,file='powell.restart')
    read(file_id,*) opt%state 
    read(file_id,*) opt%xopt 
    read(file_id,*) opt%maxfun
    read(file_id,*) opt%rhobeg
    read(file_id,*) opt%rhoend
    read(file_id,*) opt%w     
    read(file_id,*) opt%nvar
    read(file_id,*) opt%iprint 
    read(file_id,*) opt%np       
    read(file_id,*) opt%nh       
    read(file_id,*) opt%nptm     
    read(file_id,*) opt%nftest   
    read(file_id,*) opt%idz      
    read(file_id,*) opt%itest    
    read(file_id,*) opt%nf       
    read(file_id,*) opt%nfm      
    read(file_id,*) opt%nfmm     
    read(file_id,*) opt%nfsav    
    read(file_id,*) opt%knew     
    read(file_id,*) opt%kopt     
    read(file_id,*) opt%ksave    
    read(file_id,*) opt%ktemp    
    read(file_id,*) opt%rhosq    
    read(file_id,*) opt%recip    
    read(file_id,*) opt%reciq    
    read(file_id,*) opt%fbeg     
    read(file_id,*) opt%fopt     
    read(file_id,*) opt%diffa    
    read(file_id,*) opt%xoptsq   
    read(file_id,*) opt%rho      
    read(file_id,*) opt%delta    
    read(file_id,*) opt%dsq      
    read(file_id,*) opt%dnorm    
    read(file_id,*) opt%ratio    
    read(file_id,*) opt%temp     
    read(file_id,*) opt%tempq    
    read(file_id,*) opt%beta     
    read(file_id,*) opt%dx       
    read(file_id,*) opt%vquad    
    read(file_id,*) opt%diff     
    read(file_id,*) opt%diffc    
    read(file_id,*) opt%diffb    
    read(file_id,*) opt%fsave    
    read(file_id,*) opt%detrat   
    read(file_id,*) opt%hdiag    
    read(file_id,*) opt%distsq   
    read(file_id,*) opt%gisq     
    read(file_id,*) opt%gqsq     
    read(file_id,*) opt%f        
    read(file_id,*) opt%bstep    
    read(file_id,*) opt%alpha    
    read(file_id,*) opt%dstep  
    close(file_id)
end subroutine

subroutine write_optim_file(xmin,fmin,file_id)
    integer,intent(in) :: file_id
    real(kind=8)   :: xmin(:),fmin
    open(unit=file_id,file='optim.out',access=trim("append"))
    write(file_id,*) xmin,fmin
    close(file_id)
end subroutine

subroutine close_file(file_id)
    integer, intent(in) :: file_id
    close(file_id)
end subroutine
!----------------------------parameter function---------------------------------
subroutine para2file2(n,x,inputstr,fileid2)
    use random
    use mpi
    implicit none
    integer,intent(in)          :: n,fileid2
    integer                     :: i,rank,ierr,fileid,k
    real(kind=dp),intent(in)    :: x(n)
    real(kind=dp)               :: dummy=0._dp
    character(len=*),intent(in) :: inputstr
    character(200)              :: line
    character(24)              :: str,str2
  
    fileid=fileid2+100 
    open(fileid,file=adjustl(trim(inputstr)))
    open(fileid2,file='temp')
    k=1 
!   print *,'x to file',k,x
    do i=1,59
!    print *,'k',k
    select case (i)
    case default
     read(fileid,'(A)') line
!     str=adjustl(trim(str))
     write(fileid2,'(A)') line
    case (17:18)
     read(fileid,'(A24,F13.7,A)') str,dummy,str2
!     str=adjustl(trim(str))
     str2=adjustl(trim(str2))
     write(fileid2,'(A24,F13.7,A1,A20)') str,x(k),' ',str2
     k=k+1
    case (20:21)        
     read(fileid,'(A24,F13.7,A)') str,dummy,str2
!     str=adjustl(trim(str))d
     str2=adjustl(trim(str2))
     write(fileid2,'(A24,F13.7,A20)') str,x(k),str2
     k=k+1
    case (28:30)
     read(fileid,'(A16,F13.7,A)') str,dummy,str2
!     str=adjustl(trim(str))
     str2=adjustl(trim(str2))
     write(fileid2,'(A16,F13.7,A,A16)') str,x(k),' ',str2
        k=k+1
    case (36:37)
     read(fileid,'(A24,F13.7,A)') str,dummy,str2
!     str=adjustl(trim(str))
     str2=adjustl(trim(str2))
     write(fileid2,'(A24,F13.7,A,A20)') str,x(k),' ',str2
        k=k+1
    case (39:40)
     read(fileid,'(A24,F13.7,A)') str,dummy,str2
!     str=adjustl(trim(str))
     str2=adjustl(trim(str2))
     write(fileid2,'(A24,F13.7,A24)') str,x(k),str2
        k=k+1
    case (45:46)
     read(fileid,'(A24,F13.7,A)') str,dummy,str2
!     str=adjustl(trim(str))
     str2=adjustl(trim(str2))
     write(fileid2,'(A24,F13.7,A20)') str,x(k),str2
        k=k+1
    case (54:56) 
     read(fileid,'(A16,F13.7,A)') str,dummy,str2
!     str=adjustl(trim(str))
     str2=adjustl(trim(str2))
     write(fileid2,'(A16,F13.7,A,A16)') str,x(k),' ',str2
        k=k+1
     end select 
    enddo

    print *,'file written'
    call RENAME(inputstr,"para.xml.restart",ierr)
    if (ierr .ne. 0) print*,'shit...moving file failed'
    call RENAME('temp',inputstr,ierr)
    if (ierr .ne. 0) print*,'shit...moving file failed'
   
    close(fileid)
    close(fileid2)
end subroutine


subroutine read4file(n,x,ex)
    use random
    implicit none
    integer,intent(in)       :: n
    integer                  :: i,fileid,k
    real(kind=dp),intent(out) :: x(n)
    character(200)            :: str,str2,filename
    logical,intent(out)      :: ex

    inquire(file="para.xml.restart",exist=ex)
    k=1
    if (ex) then
        print *,'read from restart file'
        filename="para.xml.restart"
    else
        filename="para.xml"
    endif


    fileid=100
    open(fileid,file=filename)
    
    do i=1,59 
    select case (i)
    case default
     read(fileid,*) str
    case (17:18)
     read(fileid,*) str,x(k),str2
     k=k+1
    case (20:21)        
     read(fileid,*) str,x(k),str2
     k=k+1
     
    case (28:30)
     read(fileid,*) str,x(k),str2
     k=k+1

    case (36:37)
     read(fileid,*) str,x(k),str2
     k=k+1

    case (39:40)
     read(fileid,*) str,x(k),str2
     k=k+1

    case (45:46)
     read(fileid,*) str,x(k),str2
     k=k+1

    case (54:56) 
     read(fileid,*) str,x(k),str2
     k=k+1
    
    end select
    end do
    if (k-1 .ne. n) write(*,*) "something is very wrong"
    end subroutine

!subroutine pres_temp_4D_wr
!  use netcdf
!  implicit none
!
!  ! This is the name of the data file we will create.
!  character (len = *), parameter :: FILE_NAME = "pres_temp_4D.nc"
!  integer :: ncid
!
!  ! We are writing 4D data, a 2 x 6 x 12 lvl-lat-lon grid, with 2
!  ! timesteps of data.
!  integer, parameter :: NDIMS = 4, NRECS = 2
!  integer, parameter :: NLVLS = 2, NLATS = 6, NLONS = 12
!  character (len = *), parameter :: LVL_NAME = "level"
!  character (len = *), parameter :: LAT_NAME = "latitude"
!  character (len = *), parameter :: LON_NAME = "longitude"
!  character (len = *), parameter :: REC_NAME = "time"
!  integer :: lvl_dimid, lon_dimid, lat_dimid, rec_dimid
!
!  ! The start and count arrays will tell the netCDF library where to
!  ! write our data.
!  integer :: start(NDIMS), count(NDIMS)
!
!  ! These program variables hold the latitudes and longitudes.
!  real :: lats(NLATS), lons(NLONS)
!  integer :: lon_varid, lat_varid
!
!  ! We will create two netCDF variables, one each for temperature and
!  ! pressure fields.
!  character (len = *), parameter :: PRES_NAME="pressure"
!  character (len = *), parameter :: TEMP_NAME="temperature"
!  integer :: pres_varid, temp_varid
!  integer :: dimids(NDIMS)
!
!  ! We recommend that each variable carry a "units" attribute.
!  character (len = *), parameter :: UNITS = "units"
!  character (len = *), parameter :: PRES_UNITS = "hPa"
!  character (len = *), parameter :: TEMP_UNITS = "celsius"
!  character (len = *), parameter :: LAT_UNITS = "degrees_north"
!  character (len = *), parameter :: LON_UNITS = "degrees_east"
!
!  ! Program variables to hold the data we will write out. We will only
!  ! need enough space to hold one timestep of data; one record.
!  real :: pres_out(NLONS, NLATS, NLVLS)
!  real :: temp_out(NLONS, NLATS, NLVLS)
!  real, parameter :: SAMPLE_PRESSURE = 900.0
!  real, parameter :: SAMPLE_TEMP = 9.0
!
!  ! Use these to construct some latitude and longitude data for this
!  ! example.
!  real, parameter :: START_LAT = 25.0, START_LON = -125.0
!
!  ! Loop indices
!  integer :: lvl, lat, lon, rec, i
!
!  ! Create pretend data. If this wasn't an example program, we would
!  ! have some real data to write, for example, model output.
!  do lat = 1, NLATS
!     lats(lat) = START_LAT + (lat - 1) * 5.0
!  end do
!  do lon = 1, NLONS
!     lons(lon) = START_LON + (lon - 1) * 5.0
!  end do
!  i = 0
!  do lvl = 1, NLVLS
!     do lat = 1, NLATS
!        do lon = 1, NLONS
!           pres_out(lon, lat, lvl) = SAMPLE_PRESSURE + i
!           temp_out(lon, lat, lvl) = SAMPLE_TEMP + i
!           i = i + 1
!        end do
!     end do
!  end do
!
!  ! Create the file. 
!  call check( nf90_create(FILE_NAME, nf90_clobber, ncid) )
!  
!  ! Define the dimensions. The record dimension is defined to have
!  ! unlimited length - it can grow as needed. In this example it is
!  ! the time dimension.
!  call check( nf90_def_dim(ncid, LVL_NAME, NLVLS, lvl_dimid) )
!  call check( nf90_def_dim(ncid, LAT_NAME, NLATS, lat_dimid) )
!  call check( nf90_def_dim(ncid, LON_NAME, NLONS, lon_dimid) )
!  call check( nf90_def_dim(ncid, REC_NAME, NF90_UNLIMITED, rec_dimid) )
!
!  ! Define the coordinate variables. We will only define coordinate
!  ! variables for lat and lon.  Ordinarily we would need to provide
!  ! an array of dimension IDs for each variable's dimensions, but
!  ! since coordinate variables only have one dimension, we can
!  ! simply provide the address of that dimension ID (lat_dimid) and
!  ! similarly for (lon_dimid).
!  call check( nf90_def_var(ncid, LAT_NAME, NF90_REAL, lat_dimid, lat_varid) )
!  call check( nf90_def_var(ncid, LON_NAME, NF90_REAL, lon_dimid, lon_varid) )
!
!  ! Assign units attributes to coordinate variables.
!  call check( nf90_put_att(ncid, lat_varid, UNITS, LAT_UNITS) )
!  call check( nf90_put_att(ncid, lon_varid, UNITS, LON_UNITS) )
!
!  ! The dimids array is used to pass the dimids of the dimensions of
!  ! the netCDF variables. Both of the netCDF variables we are creating
!  ! share the same four dimensions. In Fortran, the unlimited
!  ! dimension must come last on the list of dimids.
!  dimids = (/ lon_dimid, lat_dimid, lvl_dimid, rec_dimid /)
!
!  ! Define the netCDF variables for the pressure and temperature data.
!  call check( nf90_def_var(ncid, PRES_NAME, NF90_REAL, dimids, pres_varid) )
!  call check( nf90_def_var(ncid, TEMP_NAME, NF90_REAL, dimids, temp_varid) )
!
!  ! Assign units attributes to the netCDF variables.
!  call check( nf90_put_att(ncid, pres_varid, UNITS, PRES_UNITS) )
!  call check( nf90_put_att(ncid, temp_varid, UNITS, TEMP_UNITS) )
!  
!  ! End define mode.
!  call check( nf90_enddef(ncid) )
!  
!  ! Write the coordinate variable data. This will put the latitudes
!  ! and longitudes of our data grid into the netCDF file.
!  call check( nf90_put_var(ncid, lat_varid, lats) )
!  call check( nf90_put_var(ncid, lon_varid, lons) )
!  
!  ! These settings tell netcdf to write one timestep of data. (The
!  ! setting of start(4) inside the loop below tells netCDF which
!  ! timestep to write.)
!  count = (/ NLONS, NLATS, NLVLS, 1 /)
!  start = (/ 1, 1, 1, 1 /)
!
!  ! Write the pretend data. This will write our surface pressure and
!  ! surface temperature data. The arrays only hold one timestep worth
!  ! of data. We will just rewrite the same data for each timestep. In
!  ! a real :: application, the data would change between timesteps.
!  do rec = 1, NRECS
!     start(4) = rec
!     call check( nf90_put_var(ncid, pres_varid, pres_out, start = start, &
!                              count = count) )
!     call check( nf90_put_var(ncid, temp_varid, temp_out, start = start, &
!                              count = count) )
!  end do
!  
!  ! Close the file. This causes netCDF to flush all buffers and make
!  ! sure your data are really written to disk.
!  call check( nf90_close(ncid) )
!  
!  print *,"*** SUCCESS writing example file ", FILE_NAME, "!"
!  end subroutine
!
!  subroutine check(status)
!    integer, intent ( in) :: status
!    
!    if(status /= nf90_noerr) then 
!      !print *, trim(nf90_strerror(status))
!      stop "Stopped"
!    end if
!  end subroutine check 
   
end module 

