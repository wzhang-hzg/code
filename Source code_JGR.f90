    program Macrobenthos_TOC_Interaction
    implicit none
	integer:: i,j,m,n,k,t_step,loop,min_loop,z0(100),cold_start,io_input,time_steps,zmax,count_num,count_step
	parameter(time_steps=20000,zmax=31)    ! maximum depth = zmax
	character(len=2):: site
	real:: R1,R2,R3,delta_z,f1,f2,f3,oc_pool_1, oc_pool_2,oc_pool_3,rmse_biomass,rmse_toc
	real:: c1(zmax,time_steps),c2(zmax,time_steps),c3(zmax,time_steps),biomass(zmax,time_steps),TOC(zmax,time_steps)
	real:: w(zmax),kv(zmax),adv_1(zmax),adv_2(zmax),adv_3(zmax),diff_1(zmax),diff_2(zmax),diff_3(zmax),poro(zmax)
	real:: bio_sum,toc_sum,bio_z0,oc_sum,bio_check(100,zmax),oc_check(100,zmax),gamma(3,time_steps),sed_thickness
	
!-------------------------------------------------------------------------------------
!Please cite Zhang,W., and K.,Wirtz, subm. Mutual dependence between sedimentary 
!organic carbon and macrobenthos resolved by mechanistic modelling. JGR-Biogeosciences
!-------------------------------------------------------------------------------------
  
    print*,'Please specify study site(SK,GB,FF,BF):'
    read(*,*)site    
    print*,'Please specify boundary input'
    print*,'Please specify sedimentation rate w (cm/day):'
    read(*,*)w(1)
	print*,'Please specify input rate of OC Pool_1 (g/g/day) into seabed:'
    read(*,*)f1                                                          ! example: 0.002(SK),0.015(GB),0.0008(FF),0.0002(BF)
    print*,'Please specify input rate of OC Pool_2 (g/g/day) into seabed:'
    read(*,*)f2                                                          ! example: 0.001(SK),0.02(GB),0.0058(FF),0.0001(BF)
    print*,'Please specify input rate of OC Pool_3 (g/g/day) into seabed:'
    read(*,*)f3                                                          ! example: 0.023(SK),0.008(GB),0.0031(FF),0.00035(BF)

!----Note that the sum of input values (f1+f2+f3) equals approximately measured TOC at depth z=0 cm--------------

  cold_start=1  ! specification of the starting type (1=from user defined initial conditions; 0=from external data of initial conditions)

 if(cold_start==1)then   !simulation starts from user-specified initial conditions
   
!------------model initialization------  
    DO i=1,zmax   ! depth from 1 to zmax cm, 1cm interval
	 c1(i,1)=0.01   !volumeric concentration of Pool 1 at time step=1,  unit: g/g  
	 c2(i,1)=0.01   !volumeric concentration of Pool 2 at time step=1,  unit: g/g  
	 c3(i,1)=0.01   !volumeric concentration of Pool 3  at time step=1,  unit: g/g  
	 biomass(i,1)=0.  !biomass at time step=1
	 w(i)=w(1)    ! sedimentation rate
     kv(i)=0.5*exp((1-i)*0.5)  ! initial bioturbation diffusivity 
     poro(i)=0.7*exp((1-i)*0.11)    ! porosity  
	 R1=20/365.        ! first-order mineralization rate constant of Pool_1 per day at z=0 cm   
	 R2=2/365.  ! first-order mineralization rate constant of Pool_2 per day at z=0 cm    
	 R3=0.02/365.  ! first-order mineralization rate constant of Pool_3 per day at z=0 cm
	 gamma(1,1)=0. ! macrobenthic uptake coefficient for Pool_1
	 gamma(2,1)=0.  ! macrobenthic uptake coefficient for Pool_2
	 gamma(3,1)=0.  ! macrobenthic uptake coefficient for Pool_3
	 sed_thickness=0.  ! accumulative sedimentation/erosion thickness (cm) 
	enddo
	
	delta_z=1. ! vertical spatial resolution = 1 cm

!-----end model initialization 

 else  !simulation starts from reading external initial fields
   io_input=100
   open(io_input,file='Initial TOC and biomass profile.dat',action='read')
   do j=1,50
   read(io_input,'(i5,4f10.5)')i,c1(i,1),c2(i,1),c3(i,1),biomass(i,1)
   enddo
   
   w(i)=w(1)    ! sedimentation rate
   kv(i)=0.5*exp((1-i)*0.5)  ! initial bioturbation diffusivity
   poro(i)=0.7*exp((1-i)*0.11)    ! porosity    
   R1=20/365.        ! first-order mineralization rate of Pool_1 per day   
   R2=2/365.  ! first-order mineralization rate of Pool_2 per day   
   R3=0.02/365.  ! first-order mineralization rate of Pool_3 per day   
   delta_z=1. ! vertical spatial resolution = 1 cm
   gamma(1,1)=0. ! macrobenthic uptake coefficient for Pool_1
   gamma(2,1)=0.  ! macrobenthic uptake coefficient for Pool_2
   gamma(3,1)=0.  ! macrobenthic uptake coefficient for Pool_3
   sed_thickness=0.  ! accumulative sedimentation/erosion thickness (cm) 
 endif
     
   
!!------Start interaction between TOC and macrobenthos 

     select case(site)  ! loop for model convergence. The minimum value depends on whether an equilibrium is achieved. Larger value is needed in BF due to a very low sedimentation rate.   
       case("FF")
       min_loop=20 
       case("GB")
       min_loop=20
       case("SK")
       min_loop=20
       case("BF")
       min_loop=100
     end select  
   
	do loop=1,min_loop  
	
    do t_step=2,time_steps  ! time step iterations for achieving equilibrium status

! Note that concentration at the first cell point, namely c1(1,t_step),c2(1,t_step)and c3(1,t_step), represent boundary input and not involved in iteration!
! Model result output for TOC and each OC pool should start from the second cell point, namely c1(2,t_step),c2(2,t_step)and c3(2,t_step)
    
     c1(1,t_step)=f1     
	 c2(1,t_step)=f2  
	 c3(1,t_step)=f3
	   
     do i=2,zmax-1  
       adv_1(i)=w(i)*(1-poro(i))*(c1(i,t_step-1)-c1(i-1,t_step-1))/delta_z+(w(i)*(1-poro(i))-w(i-1)*(1-poro(i-1)))*c1(i,t_step-1)/delta_z   !advection       
       diff_1(i)=kv(i)*(1-poro(i))*(c1(i+1,t_step-1)+c1(i-1,t_step-1)-2*c1(i,t_step-1))/delta_z/delta_z+(kv(i)*(1-poro(i))-kv(i-1)*(1-poro(i-1)))*(c1(i,t_step-1)-c1(i-1,t_step-1))/delta_z/delta_z   ! diffusion
      
       c1(i,t_step)=c1(i,t_step-1)+(diff_1(i)-adv_1(i)-(1-poro(i))*c1(i,t_step-1)*(R1*exp((2-i)*0.3)+biomass(i,t_step-1)*gamma(1,t_step-1)))/(1-poro(i))    !update c1 
       
       if(c1(i,t_step).le.0.)c1(i,t_step)=0.
       
     enddo
     c1(zmax,t_step)=min(c1(zmax-1,t_step),c1(1,t_step))

        
     do i=2,zmax-1    
       adv_2(i)=w(i)*(1-poro(i))*(c2(i,t_step-1)-c2(i-1,t_step-1))/delta_z+(w(i)*(1-poro(i))-w(i-1)*(1-poro(i-1)))*c2(i,t_step-1)/delta_z   !advection         
       diff_2(i)=kv(i)*(1-poro(i))*(c2(i+1,t_step-1)+c2(i-1,t_step-1)-2*c2(i,t_step-1))/delta_z/delta_z+(kv(i)*(1-poro(i))-kv(i-1)*(1-poro(i-1)))*(c2(i,t_step-1)-c2(i-1,t_step-1))/delta_z/delta_z    ! diffusion
       c2(i,t_step)=c2(i,t_step-1)+(diff_2(i)-adv_2(i)-(1-poro(i))*c2(i,t_step-1)*(R2*exp((2-i)*0.3)+biomass(i,t_step-1)*gamma(2,t_step-1)))/(1-poro(i))    !update c2
       
       if(c2(i,t_step).le.0.)c2(i,t_step)=0.
       
     enddo
     c2(zmax,t_step)=min(c2(zmax-1,t_step),c2(1,t_step))
    
     
     do i=2,zmax-1    
       adv_3(i)=w(i)*(1-poro(i))*(c3(i,t_step-1)-c3(i-1,t_step-1))/delta_z+(w(i)*(1-poro(i))-w(i-1)*(1-poro(i-1)))*c3(i,t_step-1)/delta_z   !advection        
       diff_3(i)=kv(i)*(1-poro(i))*(c3(i+1,t_step-1)+c3(i-1,t_step-1)-2*c3(i,t_step-1))/delta_z/delta_z+(kv(i)*(1-poro(i))-kv(i-1)*(1-poro(i-1)))*(c3(i,t_step-1)-c3(i-1,t_step-1))/delta_z/delta_z     ! diffusion
       c3(i,t_step)=c3(i,t_step-1)+(diff_3(i)-adv_3(i)-(1-poro(i))*c3(i,t_step-1)*(R3*exp((2-i)*0.3)+biomass(i,t_step-1)*gamma(3,t_step-1)))/(1-poro(i))   !update c3
       
       if(c3(i,t_step).le.0.)c3(i,t_step)=0.
       
     enddo
     c3(zmax,t_step)=min(c3(zmax-1,t_step),c3(1,t_step))
          
     
     if(sed_thickness.ge.1)then    ! check whether an upward shift of grid cell index due to deposition is needed
       
       do i=3,zmax
         c1(i,t_step)=(c1(i-1,t_step)+c1(i,t_step))/2.   ! Smoothing is implemented to avoid computational instability
         c2(i,t_step)=(c2(i-1,t_step)+c2(i,t_step))/2.  
         c3(i,t_step)=(c3(i-1,t_step)+c3(i,t_step))/2.    
       enddo          
         
       
       sed_thickness=0.
       
     endif
     
      if(sed_thickness.lt.-1)then    ! check whether a downward shift of grid cell index due to erosion is needed
       
       do i=2,zmax-1
         c1(i,t_step)=(c1(i+1,t_step)+c1(i,t_step))/2.   ! Smoothing is implemented to avoid computational instability
         c2(i,t_step)=(c2(i+1,t_step)+c2(i,t_step))/2.   
         c3(i,t_step)=(c3(i+1,t_step)+c3(i,t_step))/2.  
       enddo
       c1(zmax,t_step)=(c1(zmax-1,t_step)+c1(zmax,t_step))/2.
       c2(zmax,t_step)=(c2(zmax-1,t_step)+c2(zmax,t_step))/2.
       c3(zmax,t_step)=(c3(zmax-1,t_step)+c3(zmax,t_step))/2.
      
       sed_thickness=0.
    
     endif

 
    if(t_step==2.or.mod(t_step,10).eq.1)then ! updating kv and biomass every 10 time steps
     bio_sum=0.
     toc_sum=0.
     oc_sum=0.
     OC_Pool_1=0.
     OC_Pool_2=0.
     OC_Pool_3=0.
     count_step=t_step
     sed_thickness=sed_thickness+w(1)*10
     do i=1,zmax
       oc_sum=oc_sum+1000.*(c1(i,t_step)/1.5+c2(i,t_step)/4.+c3(i,t_step)/25.)
       OC_Pool_1=OC_Pool_1+1000.*c1(i,t_step)
       OC_Pool_2=OC_Pool_2+1000.*c2(i,t_step)
       OC_Pool_3=OC_Pool_3+1000.*c3(i,t_step)
     enddo
     
      
     do i=1,zmax
     select case(site)
       case("FF")
       biomass(i,t_step)=1000.*(c1(i,t_step)/1.5+c2(i,t_step)/4.+c3(i,t_step)/25.)*exp((1-i)*0.11)*oc_sum/30./(0.21+0.64*exp((1-i)*0.22))  !Station FF  
       case("GB")
       biomass(i,t_step)=1000.*(c1(i,t_step)/1.5+c2(i,t_step)/4.+c3(i,t_step)/25.)*exp((1-i)*0.11)*oc_sum/30./(1.3+0.75*exp((1-i)*0.22))  !Station GB
       case("SK")
        biomass(i,t_step)=1000.*(c1(i,t_step)/1.5+c2(i,t_step)/4.+c3(i,t_step)/25.)*exp((1-i)*0.11)*oc_sum/30./(0.7+0.5*exp((1-i)*0.22))  !Station SK 
       case("BF")
       biomass(i,t_step)=1000.*(c1(i,t_step)/1.5+c2(i,t_step)/4.+c3(i,t_step)/25.)*exp((1-i)*0.11)*oc_sum/30./(0.017+0.02*exp((1-i)*0.22)) !Station BF
     end select
       bio_sum=bio_sum+biomass(i,t_step)
       kv(i)=0.22*(biomass(i,t_step)**1.334)/(1000.*(c1(i,t_step)/1.5+c2(i,t_step)/4.+c3(i,t_step)/25.))   
       kv(i)=min(kv(i),0.52)   ! An upper limit is introduced to stabalize computation
       TOC(i,t_step)=1000.*(c1(i,t_step)+c2(i,t_step)+c3(i,t_step))
       toc_sum=toc_sum+TOC(i,t_step)
     enddo 
       gamma(1,t_step)=OC_Pool_1/(205.*zmax)  
       gamma(2,t_step)=OC_Pool_2/(1600.*zmax)
       gamma(3,t_step)=OC_Pool_3/(250000.*zmax)
     
     bio_z0=0.
     z0(1)=0
     do i=2,zmax
      bio_z0=bio_z0+biomass(i,t_step)
      if(bio_z0.le.(bio_sum/2.))z0(1)=z0(1)+1
     enddo   
     
    if(loop.ge.20)then
    rmse_biomass=0.
    rmse_toc=0.
    count_num=0
    do i=1,zmax
    rmse_biomass=rmse_biomass+(biomass(i,t_step)-biomass(i,t_step-1))**2./zmax
    rmse_toc=rmse_toc+(TOC(i,t_step)-TOC(i,t_step-1))**2./zmax
    enddo
    rmse_biomass=sqrt(rmse_biomass)/(bio_sum/zmax)
    rmse_toc=sqrt(rmse_toc)/(toc_sum/zmax)
         if(rmse_biomass.le.0.01.and.rmse_toc.le.0.01)count_num=count_num+1
         
         if(count_num.ge.2000)then
         print*,'Equilibrium status achieved'
         exit
         endif
     endif      
  
  else
    do i=1,zmax
     biomass(i,t_step)=biomass(i,count_step)
    enddo
     gamma(1,t_step)=gamma(1,count_step)
     gamma(2,t_step)=gamma(2,count_step)
     gamma(3,t_step)=gamma(3,count_step)     
    endif   ! end updating kv and biomass

                  
     enddo   ! end iteration for time_step calculation
     
     enddo  ! end loop
     
 
!----Examples for model results output-----------------
       open(1,file='OC_Pool_1.dat',action='write')
	   open(2,file='OC_Pool_2.dat',action='write')
	   open(3,file='OC_Pool_3.dat',action='write')
	   open(4,file='TOC.dat',action='write')
	   open(5,file='Macrobenthic_biomass.dat',action='write')
	   open(6,file='Bioturbation_diffusivity.dat',action='write')

    
      do j=2,zmax
      write(1,'(f10.5)')1000.*c1(j,time_steps)     ! verticial distribution of OC_Pool 1, unit: mg/g
      write(2,'(f10.5)')1000.*c2(j,time_steps)     ! verticial distribution of OC_Pool 2, unit: mg/g
      write(3,'(f10.5)')1000.*c3(j,time_steps)     ! verticial distribution of OC_Pool 3, unit: mg/g
      write(4,'(f10.5)')1000.*(c1(j,time_steps)+c2(j,time_steps)+c3(j,time_steps))  ! verticial distribution of TOC, unit: mg/g
      write(5,'(f10.5)')biomass(j,time_steps)      ! verticial distribution of macrobenthic biomass, unit: *100 g/m**3  
      write(6,'(f10.5)')kv(j)                 ! verticial distribution of bioturbation diffusivity, unit: cm**2/day
        
      enddo

    do i=1,6
    close(i)
    enddo
!----end model results output-----------------    
    
     
	end  ! end program


   
