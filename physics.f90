Module mp_physics

! --------------------------------------------------------------------------
! MPI Version of the Cartesian Solar Wave Simulator.
! Copyright 2008, Shravan Hanasoge
                                                                                                                                                         
! W. W. Hansen Experimental Physics Laboratory
! Stanford University, Stanford, CA 94305, USA
! Email: shravan@solar.stanford.edu
! --------------------------------------------------------------------------

! SUBROUTINES IN THIS MODULE:
!
! INITIALIZES THINGS REQUIRED FOR COMPUTING RHS
!
! SPONGE BCs - COMPUTE RHS FOR QUIET, MHD, FLOWS
!
! FILTERING
!
! SIXTH-ORDER LAGRANGE INTERPOLATION FOR FORCING FUNC.
!


  use initialize
  use all_modules
  use derivatives
  use damping
  implicit none

Contains

!============================================================================

 SUBROUTINE mp_Initialize_RHS()

   implicit none
   integer i,j,k,l,flag(2),trans,counter
   real*8 temp_r,temp_t,temp_p,in(nz), out(nz),advect0(nx,ny,nz),alphaz(nz), cmax
   real*8 temp,output,temp2,start_damping,cchar, finish_damping
   real*8, allocatable, dimension(:,:,:) :: tempin, tempout, dampz, functiondecay, c2a

   ! Determination of excitation location.

   temp2 = (1.0 + excitdep/Rsun)
!   temp2 = 0.9997
   temp = 1.0
   do k=1,nz
    if (abs(z(k) - temp2) .LE. temp) then
     temp = abs(z(k) - temp2)
     e_rad = k
    endif
   enddo


   ! Determination of the output height 
   temp_r = 1.0 + obsheight/rsun
   temp = 1.0
   do k=1,nz
    if (abs(z(k) - temp_r) .LE. temp) then
     temp = abs(z(k) - temp_r)
     o_rad = k
    endif
   enddo


   ! COMPUTING EFFECT OF WILSON DEPRESSION
   if (MAGNETIC) then 
    do j=1,dim2(rank)
     do i=1,nx
     
      temp = 1.0
      do k=e_rad-50,e_rad + 10
       if ( abs(rhoq(e_rad) - rho0(i,j,k)) .le. temp) then
 	temp = abs(rhoq(e_rad) - rho0(i,j,k))
 	erad_2d(i,j) = k
       endif
      enddo

      temp = 1.0
      do k=o_rad-50,o_rad + 10
       if ( abs(rhoq(o_rad) - rho0(i,j,k)) .le. temp) then
	temp = abs(rhoq(o_rad) - rho0(i,j,k))
	orad_2d(i,j) = k
       endif
      enddo

      if (compute_adjoint) delta_width_2d(i,j) = &
		2.0/(z(orad_2d(i,j)+1) - z(orad_2d(i,j)-1))
    
      if (compute_forward) delta_width_2d(i,j) = &
		2.0/(z(erad_2d(i,j)+1) - z(erad_2d(i,j)-1))


     enddo
    enddo
   endif


   if (compute_adjoint) delta_width = 2.0d0/(z(o_rad+1) - z(o_rad-1))
   if (compute_forward) delta_width = 2.0d0/(z(e_rad+1) - z(e_rad-1))

   ! Computing the absorbent boundary adjoining sponge.
 
   do k=1,nz
    ! Sponge in the z-direction
    spongez(k) =  200.0/(1.0 + exp((1.002 - z(k))/0.0001)) & 
		+ 300.0/(1.0 + exp((z(k) - (z(1) + 0.005))/0.0005))
! Source depth is where the sources are located
    source_dep(k) = exp(-10.0**(7.0)*4.0 *(z(k) - z(e_rad))**2.0)

   enddo

   if (magnetic) then
 
    call READ_IN_FIELD   ! IN ALL_MODULES.f90
  
    do k=1,nz
     reduction(:,:,k) = (boz(:,:,k)**2.+ box(:,:,k)**2.+boy(:,:,k)**2.)/(rho0(:,:,k)*cq(k)**2.)
    enddo


    reduction = 8./(8. + reduction)

    do k=nz-10, nz
     reduction(:,:,k) = reduction(:,:,k) * ((k - nz)/10.0)**4.
    enddo

    box = box * reduction**0.5
    boy = boy * reduction**0.5
    boz = boz * reduction**0.5

    call ddx(boy, curlboz, 1)
    call ddy(box, reduction, 1)
    curlboz = curlboz - reduction 
  
    call ddy(boz, curlbox, 1)
    call ddz(boy, reduction, 1)
     do k=nz-40,nz 
     reduction(:,:,k) = reduction(:,:,k)/(1. + exp(k-(nz-35.0)) )
    enddo

    do k=1,40
     reduction(:,:,k) = reduction(:,:,k)/(1. + exp(k-(35.0)) )
    enddo

    reduction = 0.0
    curlbox = curlbox - reduction 

    call ddz(box, curlboy, 1)
    call ddx(boz, reduction, 1)

    do k=1,40
     curlboy(:,:,k) = curlboy(:,:,k)/(1. + exp(k-(35.0)) )
    enddo

    curlboy = 0.0
    curlboy = curlboy - reduction

    deallocate(reduction)
 
   endif
    
 
   ! Damping length from the horizontal boundaries
   start_damping = 10.0*10.0**(8.0)/xlength
   finish_damping = 1. -  10.0*10.0**(8.0)/xlength

  spongex = 0.0
  spongey = 0.0

  if (.NOT. PERIODIC) then
   ! Sponge in the x-direction
    do i=1,nx
     spongex(i) = (1.0 - 1.0/( (1.0 + exp((start_damping-x(i))/0.008))*&
				(1.0 + exp((x(i)-finish_damping)/0.008)) ))
    enddo 

   ! Sponge in the distributed y-direction

    do i=1,dim2(rank)
     spongey(i) = (1.0 - 1.0/( (1.0 + exp((start_damping-y(i+ystart(rank)-1))/0.008))&
			*(1.0 + exp((y(i+ystart(rank)-1)-finish_damping)/0.008)) ))
    enddo 

   endif

! Total sponge
   spongex = 0.0
   spongey = 0.0
   if (PERIODIC) then 
     do k=1,nz
     do j=1,dim2(rank)
      spongexyz(:,j,k) = spongez(k)
     enddo
    enddo
   else
     do k=1,nz
     do j=1,dim2(rank)
      spongexyz(:,j,k) = (spongez(k) + 100.*spongex(:)*spongey(j))
     enddo
    enddo
   endif 

   rhoinv = 1.0/rho0
   c2rho0 = c2 * rho0
   
   call ddz(rho0, gradrho0_z, 1)
   call ddz(p0, gradp0_z, 1)


! Stabilize the background model

  if (STABILIZE_MODEL) then

   allocate(c2a(nx,dim2(rank),nz))
!    FOR CSM_B
   do k=1,nz
    do j=1,dim2(rank)
     do i=1,nx

     c2a(i,j,k) = c2(i,j,k) *(1. + 0.1 * exp(-((Rsun-Rsun*z(k))*10.0**(-8.) )**2. ) )
     c2(i,j,k) = c2a(i,j,k)* (1. - 0.03 * exp(-((Rsun-Rsun*z(k))/(5*10.0**(8.)))**2.))

     c_speed(i,j,k) = c2(i,j,k)**0.5
     c2rho0(i,j,k) = c2(i,j,k)*rho0(i,j,k)

     enddo
     enddo
    enddo
    deallocate(c2a)

    do k=1,nz
     do j=1,dim2(rank)
      do i=1,nx

      ! FOR NEAR SURFACE LAYERS ONLY, (-0.15<z<0.1)                                                     
       gradp0_z(i,j,k) =max(c2(i,j,k)*gradrho0_z(i,j,k)+10.**(-4.)*Rsun/(dimrho*dimc**2.),-rho0(i,j,k)*g(k))

      ! AND FOR THE REST OF THE BOX z<= -0.15 and z>=0.1 Mm
      if (z(k) .LE. 0.9998) then
       gradp0_z(i,j,k) =max(c2(i,j,k)*gradrho0_z(i,j,k),-rho0(i,j,k)*g(k))
      endif

      if (z(k) .GE. 1.0001) then
       gradp0_z(i,j,k) =max(c2(i,j,k)*gradrho0_z(i,j,k),-rho0(i,j,k)*g(k))
      endif

     enddo
    enddo
   enddo


  endif

   if (use_pml) then
    allocate(functiondecay(nx,dim2(rank),nz),dampz(nx,dim2(rank),nz),tempin(nx,dim2(rank),nz))
!    spongexyz = 0.
 
    if (magnetic) tempin= (c2 + (box**2. + boy**2. + boz**2.)/rho0)**0.5
    if (.not. magnetic) tempin = c2**0.5
    dampz = 0.
    alphaz = pi*0.005*diml/dimc
    functiondecay = 1.

    do k=1,nz
     if (k < npmlbot+1)  then
      alphaz(k) = pi * 0.005 * diml/dimc*(z(k)-z(1))/(z(npmlbot)-z(1))
      dampz(:,:,k) = 2.*3.*3.*tempin(:,:,1)/(z(npmlbot) - z(1))*((z(npmlbot)-z(k))/(z(npmlbot)-z(1)))**2.

      functiondecay(:,:,k) = alphaz(k)/(dampz(:,:,k) + alphaz(k)) 

      gradp0_z(:,:,k)  =gradp0_z(:,:,k) * functiondecay(:,:,k)
      gradrho0_z(:,:,k)  =gradrho0_z(:,:,k) * functiondecay(:,:,k)

      g(k)  =g(k) * functiondecay(1,1,k)
     endif

     if (k .gt. nz-npmltop) then
      alphaz(k) = pi * 0.005 * diml/dimc * (z(nz) -z(k))/(z(nz)-z(nz-npmltop+1))
      dampz(:,:,k) = 3.*3.*tempin(:,:,nz)/(z(nz) - z(nz-npmltop+1))*((z(nz-npmltop +1)-z(k))/(z(nz)-z(nz-npmltop+1)))**2.
 
      functiondecay(:,:,k) = alphaz(k)/(dampz(:,:,k) + alphaz(k)) 
      gradp0_z(:,:,k)  =gradp0_z(:,:,k) * functiondecay(:,:,k)
      gradrho0_z(:,:,k)  =gradrho0_z(:,:,k) * functiondecay(:,:,k)
      g(k)  =g(k) * functiondecay(1,1,k)
     endif

     bzpml(:,:,k) = - alphaz(k) - dampz(:,:,k)
     az(:,:,k) = -dampz(:,:,k)


     ! NEW PML STUFF - ONLY DEAL WITH THE LOWER BOUNDARY
     spongexyz(:,:,1:90) = 0.0

    enddo
    c_speed = c2**0.5

    deallocate(tempin, dampz, functiondecay)

    scrpml = 0.0 
    pmlvars = 0.0
    psipml = 0.0
   endif

   if (HORIZONTAL_PMLS) THEN

    do i=1,npmlhor
     axpml(i,:,:) = - 0.25*3. * 3. * c2(1,:,:)**0.5/&
     (x(npmlhor) - x(1)) * ( (x(i) - x(npmlhor))/(x(npmlhor)-x(1)) )**2.
     bxpml(i,:,:) = axpml(i,:,:) - 0.005 * &
     diml/dimc * pi * (x(i) - x(1))/(x(npmlhor) - x(1)) 
    enddo
  
    counter = npmlhor + 1 
    do i=nx-npmlhor+1,nx
     axpml(counter,:,:) = - 0.25*3. * 3. * c2(1,:,:)**0.5/&
     (x(nx) - x(nx-npmlhor+1)) * ( (x(nx-npmlhor+1) - x(i))/(x(nx-npmlhor+1)-x(nx)) )**2.

     bxpml(counter,:,:) = axpml(counter,:,:) - 0.005 * &
     diml/dimc * pi * (x(nx) - x(i))/(x(nx) - x(nx-npmlhor+1)) 
     counter = counter + 1
    enddo
 

    if (PROC_HAS_PML)  then
       aypml = 0.0
       bypml = 0.0
    endif

    if (ny > npmlhor) then
     do j=1,dim2(rank)

      if (j + ystart(rank) -1 .le. npmlhor) then 
       aypml(:,j,:) = - 0.5*3. * 3.* c2(:,j,:)**0.5/&
       (y(npmlhor) - y(1)) * ((y(j+ystart(rank)-1) - y(npmlhor))/(y(npmlhor) - y(1)))**2.

       bypml(:,j,:) = - 0.5*3. * 3.* c2(:,j,:)**0.5/&
       (y(npmlhor) - y(1)) * ((y(j+ystart(rank)-1) - y(npmlhor))/(y(npmlhor) - y(1)))**2. &
	- 0.005 * diml/dimc * pi * (y(j+ystart(rank)-1) -y(1))/(y(npmlhor) - y(1))

      elseif (j + ystart(rank)-1 .ge. ny - npmlhor + 1) then 
       aypml(:,j,:) = - 0.5*3. * 3.* c2(:,j,:)**0.5/&
       (y(ny) - y(ny-npmlhor+1)) * ((y(ny-npmlhor+1) - y(j+ystart(rank)-1))/&
       (y(ny) - y(ny-npmlhor+1)))**2.
       bypml(:,j,:) = - 0.5*3. * 3.* c2(:,j,:)**0.5/(y(ny) - y(ny-npmlhor+1)) * &
       ((y(ny-npmlhor+1) - y(j+ystart(rank)-1))/(y(ny) - y(ny-npmlhor+1)))**2. &
 	- 0.005 * diml/dimc * pi * (y(ny) - y(j+ystart(rank)-1))/(y(ny) - y(ny-npmlhor+1))

      endif

     enddo
    endif 

   endif  
 

   if (magnetic) then

    call ddx(rho0, gradrho0_x, 1)
    call ddy(rho0, gradrho0_y, 1)

    call ddx(p0, gradp0_x, 1)
    call ddy(p0, gradp0_y, 1)

   endif
  
   if (FLOWS) then

!    call READ_IN_FLOWS  ! IN ALL_MODULES.f90
   
    if (.not. displ) then
     allocate(tempin(nx, dim2(rank), nz), tempout(nx, dim2(rank), nz))

     call curl(v0_x, v0_y, v0_z, omega0_x, omega0_y, omega0_z)
     call cross(v0_x, v0_y, v0_z, omega0_x, omega0_y, omega0_z, advect0_x, advect0_y, advect0_z)

     call ddx(v0_x, tempin, 1)  
     advect0_x = advect0_x - tempin * v0_x

     call ddy(v0_y, tempout, 1) 
     advect0_y = advect0_y - tempout * v0_y

     call ddz(v0_z, div_v0, 1) 
     advect0_z = advect0_z - div_v0 * v0_z
  
     div_v0 = div_v0 + tempin + tempout

     c2div_v0 = c2 * div_v0
     advect0_x = advect0_x/rho0
     advect0_y = advect0_y/rho0
     advect0_z = advect0_z/rho0
     deallocate(tempin, tempout)
    endif
   endif


   if (rank == 0) then 
     print *,'Source excitation at radius', z(e_rad)
     print *,'The corresponding radial gridpoint', e_rad
   endif

 if (DAMP_WAVES) call INIT_DAMPING

end SUBROUTINE mp_Initialize_RHS



!==================================================================================


SUBROUTINE MP_QUIET_SPONGE_3D

 ! CONSERVATIVE FORM OF THE FLUXES.

 ! I can't really do the conservative form - doesn't work for the incredibly small density 
 ! and forcing functions anyway. I don't know what the error is anyway
 implicit none
 integer k,i,j, bc

 bc  = 1
 call ddxyz(v_x, dvxdx, v_y, dvydy, v_z, dvzdz, bc)

 div = dvxdx + dvydy + dvzdz

 do j=1,dim2(rank)
  do i=1,nx
   gradp_z(i,j,1)  =   (- c_speed(i,j,1)*rho0(i,j,1)*dvzdz(i,j,1)   &
		  	    & - rho(i,j,1)*g(1))*unstretch(1)
   gradp_z(i,j,nz)  =  (c_speed(i,j,nz)*rho0(i,j,nz)*dvzdz(i,j,nz)     &
      			    & - rho(i,j,nz)*g(nz))*unstretch(nz)
  enddo
 enddo


 bc = 4

 call ddxyz(p, gradp_x, p, gradp_y, p, gradp_z, bc)

 RHScont = - gradrho0_z * v_z - rho0 * div - spongexyz*rho            

 RHSv_x = - rhoinv * gradp_x - spongexyz * v_x

 RHSv_y = - rhoinv * gradp_y - spongexyz * v_y

 RHSv_z = - rhoinv * gradp_z - spongexyz * v_z

 do k= 1,nz
   RHSv_z(:,:,k) = RHSv_z(:,:,k) - rho(:,:,k)*g(k)*rhoinv(:,:,k)
   if (ABS(k-e_rad) .LE. 4) &
    RHSv_z(:,:,k) = RHSv_z(:,:,k) +  forcing*source_dep(k)
 enddo

 RHSp = - c2rho0 * div - v_z * gradp0_z  - spongexyz * p


  if (DAMP_WAVES) call DAMP_VELOCITY

 END SUBROUTINE MP_QUIET_SPONGE_3D

!================================================================================================

SUBROUTINE MP_MHD_SPONGE_3D

 ! CONSERVATIVE FORM OF THE FLUXES.

 ! I can't really do the conservative form - doesn't work for the incredibly small density 
 ! and forcing functions anyway. I don't know what the error is anyway
 implicit none
 integer i,j,k, bc

 bc = 1

 call ddxyz(v_x, dvxdx, v_y, dvydy, v_z, dvzdz, bc)

 flux1 = - (v_z * boy - v_y * boz)
 flux2 = - (v_x * boz - v_z * box)
 flux3 = - (v_y * box - v_x * boy)

! call ddz(v_z, dvzdz, 1)
! call ddz(p, gradp_z, 1)

! call ddz(flux2, dzflux2, 1)
! call ddz(flux1, dzflux1, 1)

! call ddz(bx, dzbx, 1)
! call ddz(by, dzby, 1)


 ! CALLING OVERLAP DERIVATIVE ROUTINE
 
 call ddxyz(flux2, RHSb_z, bz, curlbx, bx, dzbx, bc)
 call ddxyz(p, gradp_x, flux3, RHSb_x, flux2, dzflux2, bc)
 call ddxyz(flux3, RHSb_y, p, gradp_y, flux1, dzflux1, bc)

 ! FLUX2, FLUX3 ARE FINISHED - I.E., ALL RELEVANT
 ! DERIVATIVES HAVE BEEN COMPUTED. WILL NOW USE
 ! THEM AS SCRATCH ARRAYS

 call ddxyz(by, curlbz, bx, flux2, p, gradp_z, bc)
 call ddxyz(bz, curlby, flux1, flux3, by, dzby, bc)

 ! NOTE THAT flux3 is overwritten in the final ddxyz call
 
! call ddx(flux3, RHSb_y, 1)
! call ddy(flux3, RHSb_x, 1)

! call ddx(flux2, RHSb_z,1)
! call ddy(flux1, flux2, 1)

 RHSb_x = RHSb_x - dzflux2
 RHSb_y =-RHSb_y + dzflux1
 RHSb_z = RHSb_z - flux3

! call ddx(v_x, dvxdx, 1)
! call ddx(p, gradp_x, 1)

! call ddy(v_y, dvydy, 1)
! call ddy(p, gradp_y, 1)

 div = dvxdx + dvydy + dvzdz

! call ddy(bz, curlbx, 1)
 curlbx = curlbx - dzby

! call ddx(bz, curlby, 1)
 curlby = -curlby + dzbx

! call ddx(by, curlbz, 1)
! call ddy(bx, flux1, 1)
 curlbz = curlbz - flux2


 ! CONTIUNUITY 
 ! --------------------------------------

 RHScont = - gradrho0_x * v_x - gradrho0_y * v_y &
           - gradrho0_z * v_z - rho0 * div - spongexyz * rho

 ! PRESSURE

 RHSp = - c2rho0 * div - v_x * gradp0_x - v_y * gradp0_y &
	- v_z * gradp0_z  - spongexyz * p 


 
 call cross(curlbx, curlby, curlbz, box, boy, boz, flux1, flux2, flux3)

 RHSv_x = -rhoinv*(gradp_x - flux1)
 RHSv_y = -rhoinv*(gradp_y - flux2)
 RHSv_z = -rhoinv*(gradp_z - flux3)


 call cross(curlbox, curlboy, curlboz, bx, by, bz, flux1, flux2, flux3)


 RHSv_x = RHSv_x + rhoinv*flux1 - spongexyz * v_x
 RHSv_y = RHSv_y + rhoinv*flux2 - spongexyz * v_y
 RHSv_z = RHSv_z + rhoinv*flux3 - spongexyz * v_z


 do k= 1,nz
   RHSv_z(:,:,k) = RHSv_z(:,:,k) - rhoinv(:,:,k) * rho(:,:,k)*g(k)

   if (abs(k-e_rad) .LE. 4) &
    RHSv_z(:,:,k) = RHSv_z(:,:,k) + forcing*source_dep(k) 

 enddo 

 do k=1,8
  scr(:,:,nz,k) = 0.
  scr(:,:,1,k) = 0.
 enddo


  if (DAMP_WAVES) call DAMP_VELOCITY

! deallocate(flux1, flux2, flux3)
 END SUBROUTINE MP_MHD_SPONGE_3D


!================================================================================================


 Subroutine MP_FLOWS_SPONGE_3D ()

  ! Right hand side computation

  implicit none
  integer i,j,k, bc
  real*8, allocatable, dimension(:,:,:) :: temp_x, temp_y, temp_z, advect

!  call derivatives_all ()

  allocate(temp_x(nx,dim2(rank),nz),temp_y(nx,dim2(rank),nz),&
   temp_z(nx,dim2(rank),nz),advect(nx,dim2(rank),nz))

  bc =1
  call ddxyz(v_x, dvxdx, v_y, dvydy, v_z, dvzdz, bc)
 
  div = dvxdx + dvydy + dvzdz

  do j=1,dim2(rank)
   do i=1,nx
    gradp_z(i,j,1)  =   (- c_speed(i,j,1)*rho0(i,j,1)*dvzdz(i,j,1)   &
 		  	    & - rho(i,j,1)*g(1))*unstretch(1)
    gradp_z(i,j,nz)  =  (c_speed(i,j,nz)*rho0(i,j,nz)*dvzdz(i,j,nz)     &
      			    & - rho(i,j,nz)*g(nz))*unstretch(nz)
   enddo
  enddo



 bc = 4

 call ddxyz(p, gradp_x, p, gradp_y, p, gradp_z, bc)

  call curl(v_x, v_y, v_z, omega_x, omega_y, omega_z)
  call cross(v0_x, v0_y, v0_z, omega_x, omega_y, omega_z, flow_x, flow_y, flow_z)  
  call cross(v_x, v_y, v_z, omega0_x, omega0_y, omega0_z, temp_x, temp_y, temp_z)  


  ! Using the vector identity that v_0 dot nalbla v' + v' dot nabla v_0 =
  ! nalbla(v_0 dot v') - v_0 X omega' - v' X omega_0 

  flow_x = flow_x + temp_x
  flow_y = flow_y + temp_y
  flow_z = flow_z + temp_z

  advect = v0_x * v_x + v0_y * v_y + v0_z * v_z

  call ddx(advect, temp_x, 1)
  call ddy(advect, temp_y, 1)
  call ddz(advect, temp_z, 1)

  flow_x = flow_x - temp_x
  flow_y = flow_y - temp_y
  flow_z = flow_z - temp_z

  deallocate(temp_x, temp_y, temp_z, advect)

  call ddx(rho, gradrho_x, 1)
  call ddy(rho, gradrho_y, 1)
  call ddz(rho, gradrho_z, 1)


 ! Gravity is assumed to be directed inwards
  ! \vec{g}  = - g \vec{e_r}

   do k = 1,nz

     RHSv_z(:,:,k)  = - rhoinv(:,:,k) *  rho(:,:,k) * g(k)

     if (ABS(k-e_rad) .LE. 4) then
        RHSv_z(:,:,k) = RHSv_z(:,:,k) +  forcing*source_dep(k)
     endif

   enddo

   RHScont = -rho0*div - v_z * gradrho0_z - v0_x * gradrho_x - v0_y * gradrho_y &
	  - v0_z * gradrho_z - rho * (div_v0 + spongexyz)

   RHSv_x = -rhoinv * gradp_x + flow_x + rho * advect0_x - spongexyz*v_x

   RHSv_y = -rhoinv * gradp_y + flow_y + rho * advect0_y - spongexyz*v_y

   RHSv_z = RHSv_z - rhoinv * gradp_z + flow_z + rho * advect0_z - spongexyz*v_z

   RHSp   = - v_z * gradp0_z - rho * c2div_v0 - c2rho0 * div - v0_x * gradp_x - v0_y * gradp_y & 
	 & - v0_z * gradp_z - spongexyz*p  
	 

  if (DAMP_WAVES) call DAMP_VELOCITY

 end SUBROUTINE MP_FLOWS_SPONGE_3D

!================================================================================================

SUBROUTINE LAGRANGE_INTERP()

  implicit none
  real*8 xt

  xt = DBLE(time)*timestep/cadforcing
  time1 = FLOOR(xt)+1
 
  xt = time*timestep 

  if (time_old .LT. time1) then
   x0 = (time1-4)*cadforcing
   x1 = (time1 -3)*cadforcing
   x2 = (time1-2)*cadforcing
   x3 = (time1-1)*cadforcing
   x4 = (time1)*cadforcing
   x5 = (time1+1)*cadforcing
   x6 = (time1+2)*cadforcing  
   
   LC0 = 0.0
   LC1 = 0.0
   LC2 = 0.0
   LC3 = 0.0
   LC4 = 0.0
   LC5 = 0.0
   LC6 = 0.0
 
   if (x0 .GE. 0) & 
    LC0 = vr(:,:,(time1-3))/((x0-x1) * (x0-x2) * (x0-x3) * (x0-x4) * (x0-x5) *(x0-x6))

   if (x1 .GE. 0) &
    LC1 = vr(:,:,(time1-2))/((x1-x0) * (x1-x2) * (x1-x3) * (x1-x4) * (x1-x5) *(x1-x6))

   if (x2 .GE. 0) &
    LC2 = vr(:,:,(time1-1))/((x2-x0) * (x2-x1) * (x2-x3) * (x2-x4) * (x2-x5) *(x2-x6))

   if (x3 .GE. 0) &
    LC3 = vr(:,:,(time1))/((x3-x0) * (x3-x1) * (x3-x2) * (x3-x4) * (x3-x5) *(x3-x6))

   LC4 = vr(:,:,(time1+1))/((x4-x0) * (x4-x1) * (x4-x2) * (x4-x3) * (x4-x5) *(x4-x6))
   
   if (time1 .LT. num_steps) &
    LC5 = vr(:,:,(time1+2))/((x5-x0) * (x5-x1) * (x5-x2) * (x5-x3) * (x5-x4) *(x5-x6))

   if (time1 .LT. (num_steps-1)) &
    LC6 = vr(:,:,(time1+3))/((x6-x0) * (x6-x1) * (x6-x2) * (x6-x3) * (x6-x4) *(x6-x5))

   time_old = time1
  endif

  forcing = LC0 * (xt-x1) * (xt-x2) * (xt-x3) * (xt-x4) * (xt-x5) *(xt-x6)+ &
  	    LC1 * (xt-x0) * (xt-x2) * (xt-x3) * (xt-x4) * (xt-x5) *(xt-x6)+ &
	    LC2 * (xt-x0) * (xt-x1) * (xt-x3) * (xt-x4) * (xt-x5) *(xt-x6)+ &
	    LC3 * (xt-x0) * (xt-x1) * (xt-x2) * (xt-x4) * (xt-x5) *(xt-x6)+ &
	    LC4 * (xt-x0) * (xt-x1) * (xt-x2) * (xt-x3) * (xt-x5) *(xt-x6)+ &
	    LC5 * (xt-x0) * (xt-x1) * (xt-x2) * (xt-x3) * (xt-x4) *(xt-x6)+ &
	    LC6 * (xt-x0) * (xt-x1) * (xt-x2) * (xt-x3) * (xt-x4) *(xt-x5)


END SUBROUTINE LAGRANGE_INTERP


!================================================================================================
SUBROUTINE FILTER_VARS_3D
  implicit none
  integer i,k,l,kst
  real*8 coeffs(0:5)
  real*8, allocatable, dimension(:,:,:) :: tempfilt,temp
  complex*16, allocatable, dimension(:,:,:) :: coeff 
  
  coeffs(0) = 0.75390625000000
  coeffs(1) = 0.41015625000000
  coeffs(2) = -0.23437500000000
  coeffs(3) = 0.08789062500000
  coeffs(4) = -0.01953125000000
  coeffs(5) = 0.00195312500000
  
  kst = 6
  allocate(temp(nx,dim2(rank),nz))

  do l=1,dimfive

   !if (USE_PML) then
    temp(:,:,2) = 0.5 * (a(:,:,3,l) + a(:,:,1,l))
    temp(:,:,3) = 0.5 * (a(:,:,4,l) + a(:,:,2,l))
    temp(:,:,4) = 0.5 * (a(:,:,5,l) + a(:,:,3,l))
    temp(:,:,5) = 0.5 * (a(:,:,6,l) + a(:,:,4,l))
    kst = 2
!
!   endif

   do k=6,nz-5
    temp(:,:,k) = coeffs(0)*a(:,:,k,l) + 0.5*coeffs(1)*(a(:,:,k-1,l) + a(:,:,k+1,l)) &
                +  0.5*coeffs(2)*(a(:,:,k-2,l) + a(:,:,k+2,l))  &
                +  0.5*coeffs(3)*(a(:,:,k-3,l) + a(:,:,k+3,l))  &
                +  0.5*coeffs(4)*(a(:,:,k-4,l) + a(:,:,k+4,l))  &
                +  0.5*coeffs(5)*(a(:,:,k-5,l) + a(:,:,k+5,l))
   enddo

   temp(:,:,nz-4) = 0.5 * (a(:,:,nz-3,l) + a(:,:,nz-5,l))
   temp(:,:,nz-3) = 0.5 * (a(:,:,nz-2,l) + a(:,:,nz-4,l))
   temp(:,:,nz-2) = 0.5 * (a(:,:,nz-1,l) + a(:,:,nz-3,l))
   temp(:,:,nz-1) = 0.5 * (a(:,:,nz,l) + a(:,:,nz-2,l))

   do k=kst,nz-1
    a(:,:,k,l) = temp(:,:,k)
   enddo 

  enddo

  deallocate(temp)


  if (mod(time,freq_filtering) .EQ. 0) then

   allocate(tempfilt(nx,dim2(rank),nz),coeff(nx/2+1,dim2(rank),nz))

   do k=1,5
    call dfftw_execute_dft_r2c(fftw_plan_fwd_x, a(:,:,:,k), coeff)

    do i=1,nx/2+1
     coeff(i,:,:) = coeff(i,:,:)*decay_coeff_x(i)
    enddo

    call dfftw_execute_dft_c2r(fftw_plan_inv_x, coeff, a(:,:,:,k))

   enddo
  
   deallocate(tempfilt, coeff)

   allocate(tempfilt(ny,dim1(rank),nz),coeff(ny/2+1,dim2(rank),nz))

!  ---- X-filtering done, now doing the y-filtering
!
  
   do k=1,5

    call transpose_3D_y(a(:,:,:,k),tempfilt)

    call dfftw_execute_dft_r2c(fftw_plan_fwd_y, tempfilt, coeff)

    do i=1,ny/2+1
     coeff(i,:,:) = coeff(i,:,:)*decay_coeff_y(i)
    enddo

    call dfftw_execute_dft_c2r(fftw_plan_inv_y, coeff, tempfilt)


    call inv_transpose_3D_y(tempfilt,a(:,:,:,k))
!-----
   enddo
   deallocate(tempfilt,coeff)
  
  endif 

END SUBROUTINE FILTER_VARS_3D

!================================================================================================
SUBROUTINE FILTER_VARS_3D_HORIZ_PML
  implicit none
  integer i,k,l,j
  real*8 coeffs(0:5)
  real*8, allocatable, dimension(:,:,:) :: temp
  
  coeffs(0) = 0.75390625000000
  coeffs(1) = 0.41015625000000
  coeffs(2) = -0.23437500000000
  coeffs(3) = 0.08789062500000
  coeffs(4) = -0.01953125000000
  coeffs(5) = 0.00195312500000
  
  allocate(temp(nx,dim2(rank),nz))

  do l=1,dimfive

   do k=npmlbot+1,nz-npmltop
    do j=1,dim2(rank) 
!     if (j+ystart(rank)-1 .gt. npmlhor+5 .and. (j+ystart(rank)-1) .lt. ny-npmlhor-4) then
      do i=1,nx
       temp(i,j,k) = coeffs(0)*a(i,j,k,l) + 0.5*coeffs(1)*(a(i,j,k-1,l) + a(i,j,k+1,l)) &
                +  0.5*coeffs(2)*(a(i,j,k-2,l) + a(i,j,k+2,l))  &
                +  0.5*coeffs(3)*(a(i,j,k-3,l) + a(i,j,k+3,l))  &
                +  0.5*coeffs(4)*(a(i,j,k-4,l) + a(i,j,k+4,l))  &
                +  0.5*coeffs(5)*(a(i,j,k-5,l) + a(i,j,k+5,l))
      enddo
     !endif 
    enddo
   enddo

   do k=npmlbot+1,nz-npmltop
    do j=1,dim2(rank) 
     !if (j+ystart(rank)-1 .gt. npmlhor+5 .and. (j+ystart(rank)-1) .lt. ny-npmlhor-4) then
      do i=1,nx!-npmlhor-4
        a(i,j,k,l) = temp(i,j,k)
      enddo
     !endif 
    enddo
   enddo

  enddo

  deallocate(temp)


END SUBROUTINE FILTER_VARS_3D_HORIZ_PML
!==================================================================================


END MODULE mp_physics
