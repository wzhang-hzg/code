    program main
    use netcdf
	parameter(im=3000,imax=10726,tm=4000,pi=3.14159265,time_steps1=20000,years_number=68,zmax=31)   ! imax=number of grid cells; time_steps1=initial iteration steps; years_number= number of years to simulate; zmax=depth range (cm) in sediment
	intrinsic random_seed, random_number 
	integer*4 i,j,m,n,k,t_step,loop,z0,time_steps,ii,year_loop,tag_index_lon(imax),tag_index_lat(imax),current_step,seed,om,hypo_tag(imax),count_num,count_step,t_day
	real*4 w(imax,zmax),kv(imax,zmax),R1,R2,R3,adv_1(zmax),adv_2(zmax),adv_3(zmax),diff_1(zmax),diff_2(zmax),diff_3(zmax),poro(imax,zmax),delta_z,suspension_feeder(imax)
    real*4 c1(imax,zmax,2),c2(imax,zmax,2),c3(imax,zmax,2),oc_pool_1,oc_pool_2,oc_pool_3,gamma(imax,3,2),pool1,pool2,pool1_degradation,x,y,bathy(imax)
	real*4 bio_sum,bio_z0,erro1,erro2,erro3,erro4,oc_sum,toc_sum,bio_sum0,oc_sum0,dis,hypoxia_count(imax),hypo_area,hypo_area_2,deadzone_area,deadzone_area_2
    real*4 xgrid(imax),ygrid(imax),biomass(imax,zmax,2),TOC(imax,zmax,2),toc0(imax),total_advection,total_degradation,de_1(30),de_2(30),de_3(30),total_uptake
    real*4 f1(imax),f2(imax),f3(imax),total_biomass,smin,total_toc,total_oc1,total_toc_2,total_oc1_2,total_biomass_2,total_diffusion,dif_1(30),dif_2(30),dif_3(30)
    real*4 biomass_shallow_s,biomass_shallow_w,toc_shallow_s,toc_shallow_w,oc1_shallow_s,oc1_shallow_w,sus_max,biomass_norderei(4)
    real*4 diffusion_shallow,advection_shallow,degradation_shallow,hypo_shallow,uptake_shallow,deadzone_shallow,x_low(1634),y_low(1634),dis2
    real*4 biomass_apr,toc_apr,oc1_apr,biomass_50m_apr,toc_50m_apr,oc1_50m_apr,biomass_oct,toc_oct,oc1_oct,biomass_50m_oct,toc_50m_oct,oc1_50m_oct
    real*4 uptake(4),diffusion(4),advection(4),degradation(4),uptake_50m(4),diffusion_50m(4),advection_50m(4),degradation_50m(4)
    real*4 check_point_x(4),check_point_y(4),kv_15  ! to print vertical profiles at four field points for paper
    CHARACTER*12 netcdf_file(68)
    character*14 result_toc_w(68),result_toc_s(68)
    character*13 result_z0_w(68),result_z0_s(68)
    character*13 result_Kv_w(68),result_Kv_s(68)
    character*15 result_susp_w(68),result_susp_s(68)
    character*15 result_hypo_w(68),result_hypo_s(68)
    character*14 result_oc1_s(68),result_oc1_w(68)
    character*18 result_bio_s(68),result_bio_w(68)
!-------------------------------------------------------------------------------------
!The TOC-Macrobenthos model is based on the following study:
!Zhang,W., and K.,Wirtz, 2018. Mutual dependence between sedimentary organic carbon and macrobenthos resolved by mechanistic modelling. JGR-Biogeosciences, 122. DOI:10.1002/2017JG003909
!-------------------------------------------------------------------------------------
! This is the name of the data file we will read.
  integer :: ncid

  ! We are reading 3D data from hydrodynamic-biogeochemicam model, a 177 x 207 lat-lon grid for 365 or 366 DAYS.
  integer, parameter :: NDIMS = 3
  integer, parameter :: NLATS = 177, NLONS = 207, NDAYS=366    ! 366 is pre-defined, and re-directed to the exact length of current processing netcdf file
  character (len = *), parameter :: LAT_NAME = "lat"
  character (len = *), parameter :: LON_NAME = "lon"
  character (len = *), parameter :: DAY_NAME = "time"
  integer :: lat_dimid, lon_dimid, day_dimid
  
  integer :: day_len

  ! For the lat lon coordinate netCDF variables.
  real :: lats(NLATS), lons(NLONS)
  real, allocatable :: days(:)    ! dynamic array for number of days 
  integer :: lat_varid, lon_varid, day_varid

  ! We will read bottom phyto-detritus,oxygen and taub fields.
  character (len = *), parameter :: OC1_NAME = "phybot"     
  character (len = *), parameter :: OXY_NAME = "o2b"
  character (len = *), parameter :: TAU_NAME = "tau"
  integer :: oc1_varid, oxy_varid,tau_varid
  integer :: dimids(NDIMS)

  ! To check the units attributes.
  character (len = *), parameter :: UNITS = "units"
  character (len = *), parameter :: OC1_UNITS = "mgC/m3"   ! phyto biomass bottom
  character (len = *), parameter :: OXY_UNITS = "ml/l" 
  character (len = *), parameter :: TAU_UNITS = "N/m2"  ! undefined unit in ECOSMO output 
  character (len = *), parameter :: LAT_UNITS = "degrees_north"
  character (len = *), parameter :: LON_UNITS = "degrees_east"
  character (len = *), parameter :: DAY_UNITS = "days since 1948-1-1"   ! this must be different from each netcdf file
  integer, parameter :: MAX_ATT_LEN = 80
  integer :: att_len
  character*(MAX_ATT_LEN) :: oc1_units_in, oxy_units_in, tau_units_in
  character*(MAX_ATT_LEN) :: lat_units_in, lon_units_in, day_units_in

  ! Read the data into these arrays.
  ! Very important note: The output from ncdump has the dimensions in "C" order, so in fortran we
  ! need to reverse the order of the array dimension!!!! e.g. ncdump -h ecosmo nc file (nc) shows tau(time,lon,lat)!  
  
  real, allocatable :: oc1_in(:,:,:), oxy_in(:,:,:), tau_in(:,:,:)    ! dynamic array for number of days 

  ! We will learn about the data file and store results in these
  ! program variables.
  integer :: ndims_in, nvars_in, ngatts_in, unlimdimid_in

  ! Loop indices
  integer :: lat, lon, day   
!  end defining netcdf variables   
!-------------------------------------------
   
! read initial TOC information!    
   open(10,file='SNS_TOC_map.dat',action='read')   ! read intial TOC map of study area

   do i=1,imax
   read(10,'(2f10.5,f10.3)')xgrid(i),ygrid(i),toc0(i)   ! TOC in %
   enddo
 
!------------model initialization------  
    DO i=1,zmax   ! depth from 1 to zmax cm, 1cm interval
	 do j=1,imax
	 w(j,i)=(0.001*toc0(j)+0.0006)*exp((1-i)*0.012)    ! site dependent sedimentation rate
     kv(j,i)=0.5*exp((1-i)*0.5)  ! initial bioturbation diffusivity 
     poro(j,i)=(0.275*toc0(j)+0.296)*exp((1-i)*0.11)    ! site dependant porosity on mud content
     gamma(j,1,1)=0.  ! macrobenthic uptake coefficient for Pool_1
	 gamma(j,2,1)=0.  ! macrobenthic uptake coefficient for Pool_2
	 gamma(j,3,1)=0.  ! macrobenthic uptake coefficient for Pool_3
     enddo  
	 R1=20/365.        ! first-order mineralization rate constant of Pool_1 per day at z=0 cm   
	 R2=2/365.         ! first-order mineralization rate constant of Pool_2 per day at z=0 cm    
	 R3=0.002/365.     ! first-order mineralization rate constant of Pool_3 per day at z=0 cm
	
	 sed_thickness=0.  ! accumulative sedimentation/erosion thickness (cm) 
	enddo  

!--------------------------------------  
!-----computation for equilibrium distribution of TOC and biomass based on intial toc content, assuming constant rate and dominant Pool_3 !
   do j=1,imax
         f1(j)=toc0(j)/1.6*(0.05/365/0.14)   ! 
         f2(j)=toc0(j)/1.6*(0.05/365/0.14)   !
         f3(j)=toc0(j)/1.6*(1.30/365/0.14)	 !                      
    enddo
!---------------------------------------------  
	
	total_biomass=0.
	total_toc=0.
	total_oc1=0.
	
	delta_z=1. ! spatial step = 1 cm
	
  do j=1,imax    ! 
	
	time_steps=time_steps1
		
    DO i=1,zmax   ! from 1 to 31 cm depth, initialization of parameters for every grid cell
	 c1(j,i,1)=0.005  !initialization at t_step=1
	 c2(j,i,1)=0.005
	 c3(j,i,1)=0.005
	 biomass(j,i,1)=0.
	enddo 
	
	count_step=1
	z0=0
	
    do t_step=2,time_steps1  ! time step loop

        c1(j,1,2)=f1(j)    ! read input at current time step
        c2(j,1,2)=f2(j)
        c3(j,1,2)=f3(j)  
!

!--------------solution of C1 at time step i---     
     do i=2,zmax-1  
       adv_1(i)=0.0022*(1-poro(j,i))*(c1(j,i,1)-c1(j,i-1,1))/delta_z+(0.0022*(1-poro(j,i))-0.0022*(1-poro(j,i-1)))*c1(j,i,1)/delta_z   !advection      
       
       diff_1(i)=kv(j,i)*(1-poro(j,i))*(c1(j,i+1,1)+c1(j,i-1,1)-2*c1(j,i,1))/delta_z/delta_z+(kv(j,i)*(1-poro(j,i))-kv(j,i-1)*(1-poro(j,i-1)))*(c1(j,i,1)-c1(j,i-1,1))/delta_z/delta_z   ! diffusion
 
       c1(j,i,2)=c1(j,i,1)+(diff_1(i)-adv_1(i)-(1-poro(j,i))*c1(j,i,1)*(R1*exp((2-i)*0.3)+biomass(j,i,1)*gamma(j,1,1)))/(1-poro(j,i))    !update c1
       
       if(c1(j,i,2).le.0.)c1(j,i,2)=0.
       
     enddo
     c1(j,zmax,2)=c1(j,zmax-1,2)
 
!---------------------------------------------         
     do i=2,zmax-1  
       adv_2(i)=0.0022*(1-poro(j,i))*(c2(j,i,1)-c2(j,i-1,1))/delta_z+(0.0022*(1-poro(j,i))-0.0022*(1-poro(j,i-1)))*c2(j,i,1)/delta_z   !advection    
       
       diff_2(i)=kv(j,i)*(1-poro(j,i))*(c2(j,i+1,1)+c2(j,i-1,1)-2*c2(j,i,1))/delta_z/delta_z+(kv(j,i)*(1-poro(j,i))-kv(j,i-1)*(1-poro(j,i-1)))*(c2(j,i,1)-c2(j,i-1,1))/delta_z/delta_z   ! diffusion

       c2(j,i,2)=c2(j,i,1)+(diff_2(i)-adv_2(i)-(1-poro(j,i))*c2(j,i,1)*(R2*exp((2-i)*0.3)+biomass(j,i,1)*gamma(j,2,1)))/(1-poro(j,i))    !update c1
       
       if(c2(j,i,2).le.0.)c2(j,i,2)=0.
       
     enddo
     c2(j,zmax,2)=c2(j,zmax-1,2)
 
!---------------------------------------------     
     
!--------------solution of C3 at time step i--         
     do i=2,zmax-1  
       adv_3(i)=0.0022*(1-poro(j,i))*(c3(j,i,1)-c3(j,i-1,1))/delta_z+(0.0022*(1-poro(j,i))-0.0022*(1-poro(j,i-1)))*c3(j,i,1)/delta_z   !advection      
       
       diff_3(i)=kv(j,i)*(1-poro(j,i))*(c3(j,i+1,1)+c3(j,i-1,1)-2*c3(j,i,1))/delta_z/delta_z+(kv(j,i)*(1-poro(j,i))-kv(j,i-1)*(1-poro(j,i-1)))*(c3(j,i,1)-c3(j,i-1,1))/delta_z/delta_z   ! diffusion
 
       c3(j,i,2)=c3(j,i,1)+(diff_3(i)-adv_3(i)-(1-poro(j,i))*c3(j,i,1)*(R3*exp((2-i)*0.3)+biomass(j,i,1)*gamma(j,3,1)))/(1-poro(j,i))    !update c1
       
       if(c3(j,i,2).le.0.)c3(j,i,2)=0.
       
     enddo
     c3(j,zmax,2)=c3(j,zmax-1,2)     
     
!   
    do i=1,zmax           ! after iteration at the current time step, update the values for next time step
       c1(j,i,1)=c1(j,i,2)
       c2(j,i,1)=c2(j,i,2)
       c3(j,i,1)=c3(j,i,2)
    enddo

!     
!----------------updating Kv and biomass every few (=user defined) time steps 
    if(t_step==2.or.mod(t_step,10).eq.1)then ! updating kv and biomass
     bio_sum=0.
     toc_sum=0.
     oc_sum=0.
     OC_Pool_1=0.
     OC_Pool_2=0.
     OC_Pool_3=0.
     count_step=t_step
     do i=1,zmax
       oc_sum=oc_sum+1000.*(c1(j,i,2)/1.5+c2(j,i,2)/4.+c3(j,i,2)/25.)
       OC_Pool_1=OC_Pool_1+1000.*c1(j,i,2)
       OC_Pool_2=OC_Pool_2+1000.*c2(j,i,2)
       OC_Pool_3=OC_Pool_3+1000.*c3(j,i,2)
     enddo
     
      
     do i=1,zmax      
       biomass(j,i,2)=1000.*(c1(j,i,2)/1.5+c2(j,i,2)/4.+c3(j,i,2)/25.)*exp((1-i)*0.11)*oc_sum/zmax/(0.05+toc0(j)*0.5+(0.3+toc0(j)*0.3)*exp((1-i)*0.22))
       kv(j,i)=0.22*(biomass(j,i,2)**1.334)/(1000.*(c1(j,i,2)/1.5+c2(j,i,2)/4.+c3(j,i,2)/25.))
       bio_sum=bio_sum+biomass(j,i,2)  

       TOC(j,i,1)=1000.*(c1(j,i,2)+c2(j,i,2)+c3(j,i,2))
       toc_sum=toc_sum+TOC(j,i,1)
     enddo 
       gamma(j,1,2)=OC_Pool_1/(205.*zmax)  
       gamma(j,2,2)=OC_Pool_2/(1600.*zmax)
       gamma(j,3,2)=OC_Pool_3/(250000.*zmax)
     
     bio_z0=0.
     z0=0
     do i=1,zmax
      bio_z0=bio_z0+biomass(j,i,2)
      if(bio_z0.le.(bio_sum/2.))z0=z0+1
      
      biomass(j,i,1)=biomass(j,i,2)
     enddo   
     
      gamma(j,1,1)=gamma(j,1,2)
      gamma(j,2,1)=gamma(j,2,2)
      gamma(j,3,1)=gamma(j,3,2)   
     
    endif   ! end updating kv and biomass

!------end updating Kv and Biomass          

     enddo   ! end iteration for time_steps
     
      
   enddo   ! end of looping all study area grid points here
   
   
!------end calculating the initial (equilibrium) status of toc and biomass based on pool3


!----start annual calculation of the variation

!------read netcdf file for the first year 1948 to derive the nearest point for interpolation with the TOC map--------------------

  call check( nf90_open(netcdf_file(1), nf90_nowrite, ncid) ) ! read the first netcdf file

  ! Get the varids of the latitude and longitude coordinate and time variables.
  call check( nf90_inq_varid(ncid, LAT_NAME, lat_varid) )
  call check( nf90_inq_varid(ncid, LON_NAME, lon_varid) )
  call check( nf90_inq_varid(ncid, DAY_NAME, day_varid) )
  
!--- Get the exact number of days in the current netcdf file
  call check( nf90_inq_dimid(ncid, DAY_NAME, day_dimid) )
  call check( nf90_inquire_dimension(ncid, day_dimid, len = day_len) )
!-----derive the exact number of days = da_len  


!---allocate the exact length of four arrays: days(:); oc1_in(:,:,:),oxy_in(:,:,:), tau_in(:,:,:)
    if (.not.allocated(days)) allocate ( days(day_len) )
    if (.not.allocated(oc1_in)) allocate ( oc1_in(NLATS,NLONS,day_len) )
    if (.not.allocated(oxy_in)) allocate ( oxy_in(NLATS,NLONS,day_len) )
    if (.not.allocated(tau_in)) allocate ( tau_in(NLATS,NLONS,day_len) )
!----------------

  ! Read the latitude and longitude, and time data.
  call check( nf90_get_var(ncid, lat_varid, lats) )
  call check( nf90_get_var(ncid, lon_varid, lons) )
  call check( nf90_get_var(ncid, day_varid, days) )

  ! Get the varids of the oc1, oxy and taub netCDF variables.
  call check( nf90_inq_varid(ncid, OC1_NAME, oc1_varid) )
  call check( nf90_inq_varid(ncid, OXY_NAME, oxy_varid) )
  call check( nf90_inq_varid(ncid, TAU_NAME, tau_varid) )
  
  ! Read the OC1,OXY and TAU data from the file.
  ! Since we know the contents of the file we know that the data
  ! arrays in this program are the correct size to hold all the data.
  call check( nf90_get_var(ncid, oc1_varid, oc1_in) )
  call check( nf90_get_var(ncid, oxy_varid, oxy_in) )
  call check( nf90_get_var(ncid, tau_varid, tau_in) )


  ! Close the file. This frees up any internal netCDF resources
  ! associated with the file.
  call check( nf90_close(ncid) )

  ! If we got this far, everything worked as expected. Yipee! 
  print *,"number of days in",netcdf_file(1), day_len
  print *,"*** SUCCESS reading",netcdf_file(1)," for interpolation purpose"



!------------end reading the 1948 netcdf file

!--------------------------------
! find the nearest points from two input files
     
     do j=1,imax 
! find the nearest OC point 
      smin=sqrt((xgrid(j)-lons(1))**2.+(ygrid(j)-lats(1))**2.)*110000.	
      tag_index_lon(j)=1
      tag_index_lat(j)=1 
      do lon = 1, NLONS
       do lat = 1, NLATS     
      dis=sqrt((xgrid(j)-lons(lon))**2.+(ygrid(j)-lats(lat))**2.)*110000.	
	  if(dis.lt.smin)then
	  smin=dis
	  tag_index_lon(j)=lon
	  tag_index_lat(j)=lat 
	  endif
       enddo
      enddo 
! end finding the nearest OC point for point j
    enddo   
    
!----deallocate the four arrays (days(:); oc1_in(:,:,:),oxy_in(:,:,:), tau_in(:,:,:)) after finding the nearest point  
  deallocate ( days )
  deallocate ( oc1_in )
  deallocate ( oxy_in )
  deallocate ( tau_in )
!--------------------------- 

! end finding the nearest points and close the files
!-------------------------------

do year_loop=1,years_number 

! Open the netcdf file for the current computation loop.
 
  call check( nf90_open(netcdf_file(year_loop), nf90_nowrite, ncid) )

  ! In this case we know that there are 3 netCDF dimensions, 3 netCDF
  ! variables to read

  ! Get the varids of the latitude and longitude coordinate and time variables.
  call check( nf90_inq_varid(ncid, LAT_NAME, lat_varid) )
  call check( nf90_inq_varid(ncid, LON_NAME, lon_varid) )
  call check( nf90_inq_varid(ncid, DAY_NAME, day_varid) )
  
!--- Get the exact number of days in the current netcdf file
  call check( nf90_inq_dimid(ncid, DAY_NAME, day_dimid) )
  call check( nf90_inquire_dimension(ncid, day_dimid, len = day_len) )
!-----derive the exact number of days = da_len  

!---allocate the exact length of four arrays: days(:); oc1_in(:,:,:),oxy_in(:,:,:), tau_in(:,:,:)
    if (.not.allocated(days)) allocate ( days(day_len) )
    if (.not.allocated(oc1_in)) allocate ( oc1_in(NLATS,NLONS,day_len) )
    if (.not.allocated(oxy_in)) allocate ( oxy_in(NLATS,NLONS,day_len) )
    if (.not.allocated(tau_in)) allocate ( tau_in(NLATS,NLONS,day_len) )
!----------------

  ! Read the latitude and longitude, and time data.
  call check( nf90_get_var(ncid, lat_varid, lats) )
  call check( nf90_get_var(ncid, lon_varid, lons) )
  call check( nf90_get_var(ncid, day_varid, days) )

  ! Get the varids of the oc1, oxy and taub netCDF variables.
  call check( nf90_inq_varid(ncid, OC1_NAME, oc1_varid) )
  call check( nf90_inq_varid(ncid, OXY_NAME, oxy_varid) )
  call check( nf90_inq_varid(ncid, TAU_NAME, tau_varid) )
  
  ! Read the OC1,OXY and TAU data from the file.
  ! Since we know the contents of the file we know that the data
  ! arrays in this program are the correct size to hold all the data.
  call check( nf90_get_var(ncid, oc1_varid, oc1_in) )
  call check( nf90_get_var(ncid, oxy_varid, oxy_in) )
  call check( nf90_get_var(ncid, tau_varid, tau_in) )       

! end reading the netcdf file, the values of oc1, oxy and tau are ready to use!        
       
       
     open(150+year_loop,file='SNS_results/'//result_oc1_s(year_loop),action='write')
     open(250+year_loop,file='SNS_results/'//result_toc_s(year_loop),action='write')
     open(350+year_loop,file='SNS_results/'//result_bio_s(year_loop),action='write')
     open(450+year_loop,file='SNS_results/'//result_z0_s(year_loop),action='write')
     open(550+year_loop,file='SNS_results/'//result_hypo_s(year_loop),action='write')
     open(650+year_loop,file='SNS_results/'//result_susp_s(year_loop),action='write')
     open(750+year_loop,file='SNS_results/'//result_Kv_s(year_loop),action='write')

     open(100+year_loop,file='SNS_results/'//result_oc1_w(year_loop),action='write')
     open(200+year_loop,file='SNS_results/'//result_toc_w(year_loop),action='write')
     open(300+year_loop,file='SNS_results/'//result_bio_w(year_loop),action='write')
     open(400+year_loop,file='SNS_results/'//result_z0_w(year_loop),action='write')
     open(500+year_loop,file='SNS_results/'//result_hypo_w(year_loop),action='write')
     open(600+year_loop,file='SNS_results/'//result_susp_w(year_loop),action='write')
     open(700+year_loop,file='SNS_results/'//result_Kv_w(year_loop),action='write')
    
     seed=760015+year_loop*1234
        
   do t_day=1,365   
           
        current_step=time_steps1+(year_loop-1)*365+t_day
!  start calculation of new status of TOC and macrobenthos at t_day        
        do j=1,imax
        
         count_step=1
         
         f1(j)=toc0(j)/1.6*(0.01/365/0.14)   
         f2(j)=toc0(j)/1.6*(0.01/365/0.14)
         f3(j)=toc0(j)/1.6*(1.30/365/0.14)	          
          
         if(oc1_in(tag_index_lat(j),tag_index_lon(j),t_day).ge.0.)then   ! only for wet cells with organic detritus
         
         f1(j)=toc0(j)/1.6*(0.01/365/0.14)+(oc1_in(tag_index_lat(j),tag_index_lon(j),t_day)/300.)*(0.15/365/0.14)
         f2(j)=toc0(j)/1.6*(0.01/365/0.14)+(oc1_in(tag_index_lat(j),tag_index_lon(j),t_day)/300.)*(0.2/365/0.14)  
!----check whether deposition/erosion occurs according to bottom shear stress
         if(tau_in(tag_index_lat(j),tag_index_lon(j),t_day).gt.0.03.and.tau_in(tag_index_lat(j),tag_index_lon(j),t_day).lt.0.1)then    ! larger than threshold for deposition of organic detritus (0.03 N/m2) and less than threshold for erosion of sediment(0.1 N/m2)
         f1(j)=toc0(j)/1.6*(0.01/365/0.14)+(f1(j)-toc0(j)/1.6*(0.01/365/0.14))*(0.1-tau_in(tag_index_lat(j),tag_index_lon(j),t_day))/0.07
         f2(j)=toc0(j)/1.6*(0.01/365/0.14)+(f2(j)-toc0(j)/1.6*(0.01/365/0.14))*(0.1-tau_in(tag_index_lat(j),tag_index_lon(j),t_day))/0.07
         endif
          
         if(tau_in(tag_index_lat(j),tag_index_lon(j),t_day).ge.0.1)then    ! larger than threshold of erosion 
          f1(j)=0.     
         endif
        
         endif        
                  
!-------calculate advection and diffusion

        c1(j,1,2)=f1(j)    ! read inputs at current step
        c2(j,1,2)=f2(j)
        c3(j,1,2)=f3(j)
!     

!--------------solution of C1 at time step i---     
     do i=2,zmax-1  
       adv_1(i)=w(j,i)*(1-poro(j,i))*(c1(j,i,1)-c1(j,i-1,1))/delta_z+(w(j,i)*(1-poro(j,i))-w(j,i-1)*(1-poro(j,i-1)))*c1(j,i,1)/delta_z   !advection      
       
       diff_1(i)=kv(j,i)*(1-poro(j,i))*(c1(j,i+1,1)+c1(j,i-1,1)-2*c1(j,i,1))/delta_z/delta_z+(kv(j,i)*(1-poro(j,i))-kv(j,i-1)*(1-poro(j,i-1)))*(c1(j,i,1)-c1(j,i-1,1))/delta_z/delta_z   ! diffusion

       c1(j,i,2)=c1(j,i,1)+(diff_1(i)-adv_1(i)-(1-poro(j,i))*c1(j,i,1)*(R1*exp((2-i)*0.3)+biomass(j,i,1)*gamma(j,1,1)))/(1-poro(j,i))    !update c1
       
       if(c1(j,i,2).le.0.)c1(j,i,2)=0.
       
     enddo
     c1(j,zmax,2)=c1(j,zmax-1,2)
 
!---------------------------------------------         
     do i=2,zmax-1  
       adv_2(i)=w(j,i)*(1-poro(j,i))*(c2(j,i,1)-c2(j,i-1,1))/delta_z+(w(j,i)*(1-poro(j,i))-w(j,i-1)*(1-poro(j,i-1)))*c2(j,i,1)/delta_z   !advection      
       
       diff_2(i)=kv(j,i)*(1-poro(j,i))*(c2(j,i+1,1)+c2(j,i-1,1)-2*c2(j,i,1))/delta_z/delta_z+(kv(j,i)*(1-poro(j,i))-kv(j,i-1)*(1-poro(j,i-1)))*(c2(j,i,1)-c2(j,i-1,1))/delta_z/delta_z   ! diffusion
 
       c2(j,i,2)=c2(j,i,1)+(diff_2(i)-adv_2(i)-(1-poro(j,i))*c2(j,i,1)*(R2*exp((2-i)*0.3)+biomass(j,i,1)*gamma(j,2,1)))/(1-poro(j,i))    !update c1
       
       if(c2(j,i,2).le.0.)c2(j,i,2)=0.
       
     enddo
     c2(j,zmax,2)=c2(j,zmax-1,2)
 
!---------------------------------------------     
     
!--------------solution of C3 at time step i--         
     do i=2,zmax-1  
       adv_3(i)=w(j,i)*(1-poro(j,i))*(c3(j,i,1)-c3(j,i-1,1))/delta_z+(w(j,i)*(1-poro(j,i))-w(j,i-1)*(1-poro(j,i-1)))*c3(j,i,1)/delta_z   !advection      
       
       diff_3(i)=kv(j,i)*(1-poro(j,i))*(c3(j,i+1,1)+c3(j,i-1,1)-2*c3(j,i,1))/delta_z/delta_z+(kv(j,i)*(1-poro(j,i))-kv(j,i-1)*(1-poro(j,i-1)))*(c3(j,i,1)-c3(j,i-1,1))/delta_z/delta_z   ! diffusion
 
       c3(j,i,2)=c3(j,i,1)+(diff_3(i)-adv_3(i)-(1-poro(j,i))*c3(j,i,1)*(R3*exp((2-i)*0.3)+biomass(j,i,1)*gamma(j,3,1)))/(1-poro(j,i))    !update c1
       
       if(c3(j,i,2).le.0.)c3(j,i,2)=0.
       
     enddo
     c3(j,zmax,2)=c3(j,zmax-1,2)
     
     
!   
    do i=1,zmax           ! after iteration at the current time step, update the values for next time step
       c1(j,i,1)=c1(j,i,2)
       c2(j,i,1)=c2(j,i,2)
       c3(j,i,1)=c3(j,i,2)
    enddo
!-------------------------    

!----------------updating Kv and biomass every few (=user defined) time steps
    if(t_day==1.or.mod(t_day,10).eq.1)then 
     bio_sum=0.
     toc_sum=0.
     oc_sum=0.
     OC_Pool_1=0.
     OC_Pool_2=0.
     OC_Pool_3=0.
     count_step=t_day
     do i=1,zmax
       oc_sum=oc_sum+1000.*(c1(j,i,2)/1.5+c2(j,i,2)/4.+c3(j,i,2)/25.)
       OC_Pool_1=OC_Pool_1+1000.*c1(j,i,2)
       OC_Pool_2=OC_Pool_2+1000.*c2(j,i,2)
       OC_Pool_3=OC_Pool_3+1000.*c3(j,i,2)
     enddo
     
      
     do i=1,zmax   
       biomass(j,i,2)=1000.*(c1(j,i,2)/1.5+c2(j,i,2)/4.+c3(j,i,2)/25.)*exp((1-i)*0.11)*oc_sum/zmax/(0.05+toc0(j)*0.5+(0.3+toc0(j)*0.3)*exp((1-i)*0.22))
       kv(j,i)=0.22*(biomass(j,i,2)**1.334)/(1000.*(c1(j,i,2)/1.5+c2(j,i,2)/4.+c3(j,i,2)/25.))
      
       TOC(j,i,1)=1000.*(c1(j,i,2)+c2(j,i,2)+c3(j,i,2))
       toc_sum=toc_sum+TOC(j,i,1)
     enddo
       
       gamma(j,1,1)=OC_Pool_1/(205.*zmax)  
       gamma(j,2,1)=OC_Pool_2/(1600.*zmax)
       gamma(j,3,1)=OC_Pool_3/(250000.*zmax)
     
     bio_z0=0.
     z0=0
     do i=1,zmax
!-------maximum growth rate method 2------scaling with biomass based on Cammen (1980)      

      if(biomass(j,i,1).gt.0.and.biomass(j,i,2).gt.(biomass(j,i,1)+0.25*(biomass(j,i,1)**0.742)))then
      biomass(j,i,1)=biomass(j,i,1)+0.25*(biomass(j,i,1)**0.742)
      endif
      if(biomass(j,i,1).gt.0.and.biomass(j,i,2).le.(biomass(j,i,1)+0.25*(biomass(j,i,1)**0.742)).and.biomass(j,i,2).ge.biomass(j,i,1))then
      biomass(j,i,1)=biomass(j,i,2) 
      endif
      
      if(biomass(j,i,1).gt.0.and.biomass(j,i,2).lt.biomass(j,i,1))then
       if(biomass(j,i,2).le.(biomass(j,i,1)-0.15*(biomass(j,i,1)**0.742)))then
           biomass(j,i,1)=biomass(j,i,1)-0.15*(biomass(j,i,1)**0.742)        ! maximum decline rate per day
        else
        biomass(j,i,1)=biomass(j,i,2)
       endif 
      endif
!-----------------------------------                 
      
      if(biomass(j,i,1).le.0.and.biomass(j,i,2).le.0)then
      biomass(j,i,1)=0. 
      endif
      
      if(biomass(j,i,1).le.0.and.biomass(j,i,2).gt.0)then
      biomass(j,i,1)=biomass(j,i,2)/2. 
      endif

      bio_sum=bio_sum+biomass(j,i,1)
       
     enddo   
     
     do i=1,zmax
      bio_z0=bio_z0+biomass(j,i,1)
      if(bio_z0.le.(bio_sum/2.))z0=z0+1
     enddo
     
    endif

!------end updating Kv and Biomass 
  
       enddo ! end looping the entire grid cells	     
! end looping the entire grid cells	   
      
   enddo ! end days loop
   
 ! Close the current netcdf file. This frees up any internal netCDF resources
  ! associated with the file.
  call check( nf90_close(ncid) )

  ! If we got this far, everything worked as expected. Yipee! 
  print *,"number of days in",netcdf_file(year_loop), day_len
  print *,"*** SUCCESS finishing simulation of year",1947+year_loop

!----deallocate the four arrays (days(:); oc1_in(:,:,:),oxy_in(:,:,:), tau_in(:,:,:)) at the end  
  deallocate ( days )
  deallocate ( oc1_in )
  deallocate ( oxy_in )
  deallocate ( tau_in )
!--------------------------- 
      
enddo  ! end year_loop  

!----netcdf function

contains
  subroutine check(status)
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
  end subroutine check  

!----------------------
  
	end


   