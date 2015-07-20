Module mp_physics_2d

! --------------------------------------------------------------------------
! MPI Version of the Cartesian Solar Wave Simulator.
! Copyright 2008, Shravan Hanasoge
                                                                                                                                                         
! W. W. Hansen Experimental Physics Laboratory
! Stanford University, Stanford, CA 94305, USA
! Email: shravan@solar.stanford.edu
! --------------------------------------------------------------------------

  use initialize
  use all_modules
  use derivatives
  use damping
  implicit none

Contains


!============================================================================

 SUBROUTINE mp_Initialize_RHS_2d()

   implicit none
   integer i,j,k,l,flag(2),trans, counter
   real*8 temp_r,temp_t,temp_p,in(nz), out(nz),advect0(nx,ny,nz),alphaz(nz), cmax
   real*8 temp,output,temp2,start_damping,cchar, finish_damping, tempxy
   real*8, allocatable, dimension(:,:,:) :: tempin, tempout, dampz, functiondecay, c2a, n2
   logical iteration

   ! Determination of excitation location.

   temp2 = (1.0 + excitdep/Rsun)
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


   ! ALLOCATE FORCING FUNCTION STUFF



  if (.not. kernel_mode)  then
   num_steps = FLOOR(DBLE(maxtime)*timestep/cadforcing) + 2
   num_steps = 2*(num_steps/2)+2

   allocate(tempin(nx, nx, num_steps))
   allocate(vr(nx,dim2(rank),num_steps))

   call readfits_local(forcingfunc,tempin,nx,nx,num_steps)

   vr(:,1,:) = tempin(:,20,:)
 
   j=1
   do k=1,num_steps
     do i=1,nx
      tempxy=(((x(i)-0.5)*xlength*10.0**(-8.0))**2.0)**0.5
      vr(i,j,k) = vr(i,j,k) /(1.0+ exp((5.-k))) * &
	     1.0/(1.0+exp((10.0-tempxy)*1.5))
 
     enddo
   enddo

   deallocate(tempin)
  !!!!!!!!!!!!!!!!
 endif

 if (damp_waves) call init_damping_2d

    if (magnetic) then
    !    call create_2d_slices ()
        
        if ((.not. kernel_mode) .or. (kernel_mode .and. COMPUTE_DATA) ) then

            call readfits(dirbacktrue//'B_x2D.fits',box,nz)
            call readfits(dirbacktrue//'B_z2D.fits',boz,nz)
            call readfits(dirbacktrue//'pressure2D.fits',p0,nz)
            call readfits(dirbacktrue//'density2D.fits',rho0,nz)
            call readfits(dirbacktrue//'soundspeed2D.fits',c_speed,nz)

        elseif (kernel_mode .and. (.not. compute_data)) then

            call readfits(dirbackmag//'pressure2D.fits',p0,nz)
            call readfits(dirbackmag//'density2D.fits',rho0,nz)

            inquire(file=directory//'model_vectorpsi_ls'//jobno//'.fits', exist = iteration)
            if (iteration) then
                allocate(psivar(nx,dim2(rank),nz))
                call readfits(directory//'model_vectorpsi_ls'//jobno//'.fits',psivar,nz)

                call ddz(psivar, box, 0)
                !call dbyd2(box,psivar,nx,dim2(rank)*nz,1)
                box = -box
                !box(:,1,nz-1) = 0.5*(box(:,1,nz) + box(:,1,nz-2))
                !box(:,1,nz) = box(:,1,nz-1)

                !       box(:,:,180:230) = 0.0

                !       call ddx(psivar, boz, 1)
                call dbyd1(boz,psivar,nx,dim2(rank)*nz,0)
                boz = boz*stretchx
                deallocate(psivar)
            endif
       
            !call readfits(dirbackmag//'B_x2D.fits',box,nz)
            !call readfits(dirbackmag//'B_z2D.fits',boz,nz)
            inquire(file=directory//'model_c_ls'//jobno//'.fits', exist = iteration)
            if (iteration) &
                call readfits(directory//'model_c_ls'//jobno//'.fits',c_speed,nz)
            if (contrib == '01') then
                call writefits_3d('Bx_ls'//jobno//'.fits',box,nz)
                call writefits_3d('Bz_ls'//jobno//'.fits',boz,nz)
            endif
        endif

     
        box = box/dimb
        boz = boz/dimb
        p0 = p0/(dimrho * dimc**2.)
        rho0 = rho0/dimrho
        c_speed = c_speed/dimc

        !print *,maxval((box**2.+boz**2.)**0.5/rho0**0.5*dimc*10.**(-5.))
        c2 = c_speed**2.
        cq = c_speed(1,1,:)
        rhoq = rho0(1,1,:)
        
        do k=1,nz
            reduction(:,:,k) = (boz(:,:,k)**2.+ box(:,:,k)**2.)/(rho0(:,:,k)*cq(k)**2.) 
        enddo
        reduction = 8./(8. + reduction)

        call ddz(box, curlboy, 1)
        call ddx(boz, reduction,1)
        curlboy = curlboy-reduction


        reduction = 1.
        do k=nz-10, nz
            reduction(:,:,k) = reduction(:,:,k) * ((k - nz)/10.0)**4.
        enddo
        reduction=1.0  
        
    !    do k=nz-30,nz 
    !     curlboy(:,:,k) = curlboy(:,:,k)/(1. + exp(k-(nz-20.0)) )
     !   enddo
    !    curlboy(:,:,nz-30:nz) = 0.
    !    curlboy = 0.
        
        box = box * reduction**0.5
        boz = boz * reduction**0.5
        curlboy = curlboy * reduction**0.5


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
      spongexyz(:,j,k) = (spongez(k) + 100.*spongex(:))
     enddo
    enddo
   endif 

   rhoinv = 1.0/rho0
   c2rho0 = c2 * rho0
   
   call ddz(rho0, gradrho0_z, 1)
   call ddz(p0, gradp0_z, 1)

   allocate(n2(nx,dim2(rank),nz))
   n2 = gradp0_z/c2 - gradrho0_z
   do k=180,nz
    do i=1,nx
!     if (n2(i,1,k) < 0) gradp0_z(i,1,k) = c2(i,1,k)*gradrho0_z(i,1,k)
  !c2(i,1,k) = abs(gradp0_z(i,1,k)/gradrho0_z(i,1,k))
    enddo
   enddo
   deallocate(n2)        

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
    spongexyz = 0.
 
    if (magnetic) then
     do k=1,nz
      tempin(:,:,k) = (c2(:,:,k) + (box(:,:,k)**2.+ boz(:,:,k)**2.)/rho0(:,:,k))**0.5
     enddo
    endif
    if (.not. magnetic) tempin = c2**0.5

    dampz = 0.
    alphaz = pi*0.005*diml/dimc
    functiondecay = 1.

    do k=1,nz
     if (k < npmlbot+1)  then
      alphaz(k) = pi * 0.005 * diml/dimc*(z(k)-z(1))/(z(npmlbot)-z(1))
      dampz(:,:,k) = 3.*3.*tempin(:,:,1)/(z(npmlbot) - z(1))*&
           ((z(npmlbot)-z(k))/(z(npmlbot)-z(1)))**2. !tempin(:,:,1)
      functiondecay(:,:,k) = alphaz(k)/(dampz(:,:,k) + alphaz(k)) 
      gradp0_z(:,:,k)  =gradp0_z(:,:,k) * functiondecay(:,:,k)
      gradrho0_z(:,:,k)  =gradrho0_z(:,:,k) * functiondecay(:,:,k)
      g(k)  =g(k) * functiondecay(1,1,k)
     endif

     if (k .gt. nz-npmltop) then
      alphaz(k) = pi * 0.005 * diml/dimc * (z(nz) -z(k))/(z(nz)-z(nz-npmltop+1))
      dampz(:,:,k) = 3.*3.*tempin(:,:,nz-5)/(z(nz) - z(nz-npmltop+1))*&
        ((z(nz-npmltop +1)-z(k))/(z(nz)-z(nz-npmltop+1)))**2. ! tempin(:,:,nz)
 
      functiondecay(:,:,k) = alphaz(k)/(dampz(:,:,k) + alphaz(k)) 
      gradp0_z(:,:,k)  =gradp0_z(:,:,k) * functiondecay(:,:,k)
      gradrho0_z(:,:,k)  =gradrho0_z(:,:,k) * functiondecay(:,:,k)
      g(k)  =g(k) * functiondecay(1,1,k)
     endif

     bzpml(:,:,k) = - alphaz(k) - dampz(:,:,k)
     az(:,:,k) = -dampz(:,:,k)


    enddo
    c_speed = c2**0.5

    deallocate(tempin, dampz, functiondecay)

    scrpml = 0.0 
    pmlvars = 0.0
    psipml = 0.0
   endif

   if (HORIZONTAL_PMLS) THEN

    do i=1,npmlhor
     axpml(i,:,:) = - 0.25*3. * 3. * &
	c2(1,:,:)**0.5/(x(npmlhor) - x(1)) * ( (x(i) - x(npmlhor))/(x(npmlhor)-x(1)) )**2.
     bxpml(i,:,:) = axpml(i,:,:) - &
	0.005 * diml/dimc * pi * (x(i) - x(1))/(x(npmlhor) - x(1)) 
    enddo
  
    counter = npmlhor + 1 
    do i=nx-npmlhor+1,nx
     axpml(counter,:,:) = - 0.25*3. * 3. * c2(1,:,:)**0.5/(x(nx) - x(nx-npmlhor+1)) * &
		( (x(nx-npmlhor+1) - x(i))/(x(nx-npmlhor+1)-x(nx)) )**2.
     bxpml(counter,:,:) = axpml(counter,:,:) - 0.005 * diml/dimc * pi * &
		(x(nx) - x(i))/(x(nx) - x(nx-npmlhor+1)) 
     counter = counter + 1
    enddo
  
   endif

   if (magnetic) then
 
    call ddx(rho0, gradrho0_x ,1)
    call ddx(p0, gradp0_x ,1)

!gradp0_x = boz*curlboy
   ! COMPUTING EFFECT OF WILSON DEPRESSION
    do j=1,dim2(rank)
     do i=1,nx
     
      temp = 1.0
      do k=e_rad-5,e_rad + 5
       if ( abs(rhoq(e_rad) - rho0(i,j,k)) .le. temp) then
 	temp = abs(rhoq(e_rad) - rho0(i,j,k))
 	erad_2d(i,j) = k
       endif
      enddo

      temp = 1.0
      do k=o_rad-50,o_rad !+ 10
       if ( abs(rhoq(o_rad) - rho0(i,j,k)) .le. temp) then
	temp = abs(rhoq(o_rad) - rho0(i,j,k))
	orad_2d(i,j) = k
       endif
      enddo

      erad_2d = e_rad
      orad_2d = o_rad

      if (compute_adjoint) delta_width_2d(i,j) = &
		2.0/(z(orad_2d(i,j)+1) - z(orad_2d(i,j)-1))
    
 !     if (compute_forward) delta_width_2d(i,j) = &
!		2.0/(z(erad_2d(i,j)+1) - z(erad_2d(i,j)-1))


     enddo
    enddo

   endif
  
   if (rank == 0) then 
     print *,'Source excitation at radius', z(e_rad)
     print *,'The corresponding radial gridpoint', e_rad
   endif

end SUBROUTINE mp_Initialize_RHS_2d


!============================================================================


SUBROUTINE MP_MHD_PML_2d


 ! NON-CONSERVATIVE FORM OF THE FLUXES.

 implicit none
 integer i,j,k,pmlindex
 real*8 tempxy, blip
 real*8, dimension(nx,dim2(rank),nz) :: temp1, temp, temp3

 temp = - (v_x * boz - v_z * box)

 ! COMPUTE ALL z-derivatives and account for the PML
 
 call ddz(v_z, dvzdz, 1)
 call ddz(p, gradp_z, 1)

 call ddz(bx, dzbx, 1)
 call ddz(temp, dzflux2, 1)

 call PML_BOUNDARY_ACCOUNTING_2D ()

 ! ALL x-derivatives
 
 call ddx(v_x, dvxdx, 1)
 call ddx(p, gradp_x, 1)

 call ddx(temp, RHSb_z, 1)
 call ddx(bz, temp, 1)

 div = dvxdx + dvzdz


 ! CONTIUNUITY FLUX
 ! --------------------------------------
 RHScont = - gradrho0_z * v_z - rho0 * div &
	   - gradrho0_x * v_x 
           
 RHSp = - c2rho0 * div - v_x * gradp0_x &
	- v_z * gradp0_z 


 temp1 =   boz * (dzbx - temp) + bz * curlboy
 temp3 = - box * (dzbx - temp) - bx * curlboy
 
 RHSv_x = rhoinv *(temp1 - gradp_x) 
 RHSv_z = rhoinv *(temp3 - gradp_z) 

 do k=1,nz
  RHSv_z(:,:,k) = RHSv_z(:,:,k)  - rhoinv(:,:,k) * rho(:,:,k) * g(k) 

  if (abs(k-e_rad) .LE. 4) &
    RHSv_z(:,:,k) = RHSv_z(:,:,k) + forcing*source_dep(k)
  
 enddo

 RHSb_x =  - dzflux2

 do k=1,8
  scr(:,:,1,k) = 0.
  scr(:,:,nz,k) = 0.
 enddo

 END SUBROUTINE MP_MHD_PML_2d

!================================================================================================

SUBROUTINE CREATE_2D_SLICES ()

 implicit none
 integer i, j, k
 real*8 tempxy
 real*8, allocatable, dimension(:,:,:) :: temp

 allocate(temp(nx, nx, nz))
 call readfits_local(dirbackmag//'density2D.fits',temp,nx,1,nz)
 rho0(:,1,:) = temp(:,1,:)/dimrho
 
 call readfits_local(dirbackmag//'pressure2D.fits',temp,nx,1,nz)
 p0(:,1,:) = temp(:,1,:)/(dimrho * dimc**2.0)

 call readfits_local(dirbackmag//'B_x2D.fits',temp,nx,1,nz)
 box(:,1,:) = temp(:,1,:)/dimb

 call readfits_local(dirbackmag//'B_z2D.fits',temp,nx,1,nz)
 boz(:,1,:) =  temp(:,1,:)/dimb

 do j=1,nz
  c2(:,:,j) = gamma(j)*p0(:,:,j)/rho0(:,:,j)
 enddo
 c_speed= c2**(0.5)

 deallocate(temp)

  
END SUBROUTINE CREATE_2D_SLICES
!================================================================================================

SUBROUTINE FILTER_VARS_2D
  implicit none
  integer i,k,l,j,kst
  real*8 coeffs(0:5)
  real*8, allocatable, dimension(:,:,:) :: temp
  real*8, allocatable, dimension(:,:,:) :: tempfilt
  complex*16, allocatable, dimension(:,:,:) :: coeff 
  
  coeffs(0) = 0.75390625000000
  coeffs(1) = 0.41015625000000
  coeffs(2) = -0.23437500000000
  coeffs(3) = 0.08789062500000
  coeffs(4) = -0.01953125000000
  coeffs(5) = 0.00195312500000
  
  allocate(temp(nx,dim2(rank),nz))
  kst = 6
  do l=1,dimfive
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

   if (USE_PML) then
    temp(:,:,5) = 0.5 * (a(:,:,6,l) + a(:,:,4,l))
    temp(:,:,4) = 0.5 * (a(:,:,5,l) + a(:,:,3,l))
    temp(:,:,3) = 0.5 * (a(:,:,4,l) + a(:,:,2,l))
    temp(:,:,2) = 0.5 * (a(:,:,3,l) + a(:,:,1,l))
!    kst = 2
   endif

   do k=kst,nz-1
    if (.not. horizontal_pmls) & 
      a(:,:,k,l) = temp(:,:,k)


    if (horizontal_pmls) & 
      a(npmlhor:nx-npmlhor + 1,:,k,l) = temp(npmlhor:nx-npmlhor+1,:,k)

   enddo
  enddo

  deallocate(temp)

  if (mod(time+1,freq_filtering) .EQ. 0 .and. (.not. HORIZONTAL_PMLS)) then
   allocate(coeff(nx/2+1,dim2(rank),nz))

   do l=1,dimfive
    call dfftw_execute_dft_r2c(fftw_plan_fwd_x, a(:,:,:,l), coeff)

    do k=1,nz
     do j=1,dim2(rank)
      do i=1,nx/2+1
       coeff(i,j,k) = coeff(i,j,k)*decay_coeff_x(i)
       enddo
     enddo
    enddo

    call dfftw_execute_dft_c2r(fftw_plan_inv_x, coeff, a(:,:,:,l))

   enddo
  
   deallocate(coeff)

  endif 

END SUBROUTINE FILTER_VARS_2D

!================================================================================================


SUBROUTINE MP_MHD_SPONGE_2D

 ! CONSERVATIVE FORM OF THE FLUXES.

 ! I can't really do the conservative form - doesn't work for the incredibly small density 
 ! and forcing functions anyway. I don't know what the error is anyway
 implicit none
 integer i,j,k
 real*8 tempxy, blip, alf(nx,dim2(rank))
 real*8 temp1(nx,dim2(rank),nz), temp2(nx,dim2(rank),nz), temp3(nx,dim2(rank),nz)
 real*8 flux1(nx,dim2(rank),nz), flux2(nx,dim2(rank),nz), flux3(nx,dim2(rank),nz)

 call ddx(v_x, dvxdx, 1)
 call ddz(v_z, dvzdz, 1)
 div = dvxdx + dvzdz

 call ddx(p, gradp_x, 1)

 do j=1,dim2(rank)
  do i=1,nx
   gradp_z(i,j,1)  =   (- c_speed(i,j,1)*rho0(i,j,1)*dvzdz(i,j,1)   &
  	    & - rho(i,j,1)*g(1))*unstretch(1)
   gradp_z(i,j,nz)  =  (c_speed(i,j,nz)*rho0(i,j,nz)*dvzdz(i,j,nz)     &
	    & - rho(i,j,nz)*g(nz))*unstretch(nz)
  enddo
 enddo

 call ddz(p, gradp_z, 4)

 ! CONTIUNUITY FLUX
 ! --------------------------------------
 RHScont = - gradrho0_z * v_z - rho0 * div &
	   - spongexyz*rho - gradrho0_x * v_x 
           

 flux3 = bz * box + boz * bx

 call ddx(bz, flux1, 1)
 call ddz(bx, flux2, 0)

 temp1 =   boz * (flux2 - flux1) + bz * curlboy
 temp3 = - box * (flux2 - flux1) - bx * curlboy

 RHSv_x = rhoinv *(temp1 - gradp_x) - spongexyz * v_x
 RHSv_z = rhoinv *(temp3 - gradp_z) - spongexyz * v_z

 do k=1,nz
  RHSv_z(:,:,k) = RHSv_z(:,:,k)  - rhoinv(:,:,k) * rho(:,:,k) * g(k) 

  if (abs(k-e_rad) .LE. 4) &
    RHSv_z(:,:,k) = RHSv_z(:,:,k) + forcing*source_dep(k)

 enddo

 RHSp = - c2rho0 * div - v_x * gradp0_x &
	- v_z * gradp0_z - spongexyz * p 


 flux2 = - (v_x * boz - v_z * box)

 call ddz(flux2, temp2, 1)
 RHSb_x = -temp2

 call ddx(flux2, RHSb_z, 1)

 do k=1,8
  scr(:,:,1,k) = 0.0
  scr(:,:,nz,k) = 0.0
 enddo

END SUBROUTINE MP_MHD_SPONGE_2D

!================================================================================================

SUBROUTINE MP_QUIET_PML_2D

 ! CONSERVATIVE FORM OF THE FLUXES.

 ! I can't really do the conservative form - doesn't work for the incredibly small density 
 ! and forcing functions anyway. I don't know what the error is anyway
 implicit none
 integer k,i,j,pmlindex


 call ddz(v_z, dvzdz, 1)
 call ddz(p, gradp_z, 1)

 ! PML Accounting for dvzdz and gradp_z
 call PML_BOUNDARY_ACCOUNTING_2D ()

 call ddx(v_x, dvxdx, 1)
 call ddx(p, gradp_x, 1)

 ! Compute divergence

 div = dvxdx + dvzdz

 ! RHS continuity and pressure

 RHScont =  - rho0 * div - gradrho0_z * v_z

 RHSp = - c2rho0 * div - v_z * gradp0_z

! RHS v_x

 RHSv_x = - rhoinv * gradp_x 

! RHS v_z and PML accounting for gradp_z

 do k= 1,nz
  RHSv_z(:,:,k) =  - rhoinv(:,:,k) * (gradp_z(:,:,k) + rho(:,:,k)*g(k))

  if (ABS(k-e_rad) .LE. 4) &
    RHSv_z(:,:,k) = RHSv_z(:,:,k) +  forcing*source_dep(k)
 enddo

 RHSv_x(:,:,nz) = 0.
 RHSv_z(:,:,nz) = 0.
 RHSv_x(:,:,1) = 0.
 RHSv_z(:,:,1) = 0.
 RHSp(:,:,1) = 0.
 RHSp(:,:,nz) = 0.
 RHScont(:,:,1) = 0.
 RHScont(:,:,nz) = 0.

 END SUBROUTINE MP_QUIET_PML_2D


!================================================================================================


SUBROUTINE MP_MHD_PML_2D_DISPL

 ! CONSERVATIVE FORM OF THE FLUXES.

 ! I can't really do the conservative form - doesn't work for the incredibly small density 
 ! and forcing functions anyway. I don't know what the error is anyway
 implicit none
 integer i,j,k, pmlindex,bc
! real*8, dimension(nx, dim2(rank), nz) :: flux1, flux2, flux3


 bc =1
 call ddx(xi_x, dxixdx, 1) 
 dxizdz(:,:,nz) = -dxixdx(:,:,nz)*unstretch(nz)
 call ddz(xi_z, dxizdz, 3)

 pmlindex = 0
 do k=1,nz

  if ((k .le. npmlbot) .or. (k .ge. (nz-npmltop+1))) then
   pmlindex = pmlindex + 1

   RHSpsivz(:,:,pmlindex) = az(:,:,k)*dxizdz(:,:,k) + bzpml(:,:,k)*psivz(:,:,pmlindex)    
   dxizdz(:,:,k) = dxizdz(:,:,k) + psivz(:,:,pmlindex)
  endif 

 enddo

 if (HORIZONTAL_PMLS) THEN
  
  do k=1,nz
   do j=1,dim2(rank)
     do i=1,npmlhor

      dxixdx(i,j,k) = dxixdx(i,j,k) + psivx(i,j,k)
      RHSpsivx(i,j,k) = axpml(i,j,k) * dxixdx(i,j,k) + bxpml(i,j,k) * psivx(i,j,k)

      dxixdx(nx-npmlhor + i,j,k) = dxixdx(nx-npmlhor + i,j,k) + psivx(i+npmlhor,j,k)
      RHSpsivx(i+npmlhor,j,k) = axpml(i+npmlhor,j,k) * dxixdx(nx-npmlhor + i,j,k) + bxpml(i+npmlhor,j,k) * psivx(i+npmlhor,j,k)

     enddo
    enddo
   enddo 

  endif

 div = dxixdx  + dxizdz
! div(:,:,nz) = 0.0

 p = - c2rho0 * div - xi_z * gradp0_z - xi_x * gradp0_x 
 rho = - rho0 * div - xi_z * gradrho0_z - xi_x * gradrho0_x 

 !p = gradp0_z/c2 - gradrho0_z
 !call writefits_3d('n2.fits',p,nz)
 !stop
 bc =0
 call ddx(p, gradp_x, 0) 
 call ddz(p, gradp_z, 0)


 pmlindex = 0
 do k=1,nz

  if ((k .le. npmlbot) .or. (k .ge. (nz-npmltop+1))) then
   pmlindex = pmlindex + 1

   RHSpsip(:,:,pmlindex) = az(:,:,k)*gradp_z(:,:,k) + bzpml(:,:,k)*psip(:,:,pmlindex)    
   gradp_z(:,:,k) =  gradp_z(:,:,k) + psip(:,:,pmlindex)
  endif 

 enddo

 if (HORIZONTAL_PMLS) THEN
  
  do k=1,nz
   do j=1,dim2(rank)
     do i=1,npmlhor

      gradp_x(i,j,k) = gradp_x(i,j,k) + psi_gradp_x(i,j,k)
      RHS_psi_gradp_x(i,j,k) = axpml(i,j,k) * gradp_x(i,j,k) + bxpml(i,j,k) * psi_gradp_x(i,j,k)

      gradp_x(nx-npmlhor + i,j,k) = gradp_x(nx-npmlhor+ i,j,k) + psi_gradp_x(i+npmlhor,j,k)
      RHS_psi_gradp_x(i+npmlhor,j,k) = axpml(i+npmlhor,j,k) * &
      gradp_x(nx-npmlhor+i,j,k) + bxpml(i+npmlhor,j,k) * psi_gradp_x(i+npmlhor,j,k)

     enddo
    enddo
   enddo 
  
 endif



 flux2 = - (xi_x * boz - xi_z * box)
 bc =1
 call ddz(flux2, bx, 0)
 
 call ddx(flux2, bz, 1)

 ! NOTE THAT flux3 is overwritten in the third ddxyz call

 ! THE PML STUFF HAS TO COME HERE UNFORTUNATELY

 pmlindex=0
 do k=1,nz

  if ((k .le. npmlbot) .or. (k .ge. (nz-npmltop+1))) then
   pmlindex = pmlindex + 1

   RHSpsiinductionbx(:,:,pmlindex) = az(:,:,k)*bx(:,:,k) &
			+ bzpml(:,:,k)*psiinductionbx(:,:,pmlindex)    

   bx(:,:,k) = bx(:,:,k) + psiinductionbx(:,:,pmlindex)


!   RHSpsiinductionby(:,:,pmlindex) = az(:,:,k)*dzflux1(:,:,k) &
!			+ bzpml(:,:,k)*psiinductionby(:,:,pmlindex)    

 !  dzflux1(:,:,k) = dzflux1(:,:,k) + psiinductionby(:,:,pmlindex)

  endif 
 
 enddo

 bx = -bx

! bx =   bx - dzflux2
! bz =   bz - flux3

 bc =1
 call ddz(bx, dzbx, 0)

 call ddx(bz, curlby, 1)


 ! THE PML STUFF HAS TO COME HERE UNFORTUNATELY

 pmlindex=0
 do k=1,nz
  if ((k .le. npmlbot) .or. (k .ge. (nz-npmltop+1))) then
    pmlindex = pmlindex + 1

    RHSpsidzbx(:,:,pmlindex) = az(:,:,k)*dzbx(:,:,k) &
			+ bzpml(:,:,k)*psidzbx(:,:,pmlindex)    

    dzbx(:,:,k) = dzbx(:,:,k) + psidzbx(:,:,pmlindex)

   endif
 enddo


 curlby = -curlby + dzbx


 flux1 = curlby * boz
 flux3 = - curlby * box 

 RHSv_x = rhoinv * (- gradp_x + flux1)
 RHSv_z = rhoinv * (- gradp_z + flux3)


 flux1 = curlboy * bz
 flux3 = - curlboy * bx


 RHSv_x = RHSv_x + rhoinv*flux1 
 RHSv_z = RHSv_z + rhoinv*flux3

 do k= 1,nz
   RHSv_z(:,:,k) = RHSv_z(:,:,k) - rhoinv(:,:,k) * rho(:,:,k)*g(k)
 enddo 

 RHSxi_x = v_x
 RHSxi_z = v_z


 if (DAMP_WAVES) call DAMP_VELOCITY



 if (time .le. forcing_length*cadforcing_step) then

   ! WILSON DEPRESSION TAKEN INTO ACCOUNT
   ! ASSUME SOURCES ARE UNIFORMLY DISTRIBUTED HORIZONTALLY
   ! BUT NOT AT THE SAME GEOMETRICAL HEIGHT
   ! SAME OPTICAL DEPTH (i.e., SAME DENSITY)

   if (compute_adjoint) then

    do j=1,dim2(rank)
     do i=1,nx
      RHSv_z(i,j,orad_2d(i,j)) = RHSv_z(i,j,orad_2d(i,j)) +  &
	forcing(i,j) * delta_width_2d(i,j)/rho0(i,j,orad_2d(i,j)) 
     enddo
    enddo

   elseif (compute_forward) then

     do j=1,dim2(rank)
      do i=1,nx
       RHSv_z(i,j,erad_2d(i,j)) = RHSv_z(i,j,erad_2d(i,j)) +  &
	  vr(i,j,1) * fwdsource(time+1)/rho0(i,j,erad_2d(i,j)) * delta_width_2d(i,j)
      enddo
     enddo

   endif
  endif


 do k=1,6
!  scr(:,:,nz,k) = 0.
  scr(:,:,1,k) = 0.
 enddo

 if (HORIZONTAL_PMLS) then
  scr(1,:,:,:) = 0.0
  scr(nx,:,:,:) = 0.0
 endif

 END SUBROUTINE MP_MHD_PML_2D_DISPL


!================================================================================================

SUBROUTINE MP_QUIET_PML_2D_DISPL

 ! CONSERVATIVE FORM OF THE FLUXES.

 ! I can't really do the conservative form - doesn't work for the incredibly small density 
 ! and forcing functions anyway. I don't know what the error is anyway
 implicit none
 integer k,i,j,pmlindex
 real*8 f(nx,dim2(rank))


 call ddx(xi_x, dxixdx, 1)
 dxizdz(:,:,nz) = -dxixdx(:,:,nz)*unstretch(nz)
 call ddz(xi_z, dxizdz, 3)

 pmlindex = 0
 do k=1,nz

  if ((k .le. npmlbot) .or. (k .ge. (nz-npmltop+1))) then
   pmlindex = pmlindex + 1

   RHSpsivz(:,:,pmlindex) = az(:,:,k)*dxizdz(:,:,k) + bzpml(:,:,k)*psivz(:,:,pmlindex)    
   dxizdz(:,:,k) = dxizdz(:,:,k) + psivz(:,:,pmlindex)
  endif 

 enddo


 if (HORIZONTAL_PMLS) THEN
  
  do k=1,nz
   do j=1,dim2(rank)
     do i=1,npmlhor

      dxixdx(i,j,k) = dxixdx(i,j,k) + psivx(i,j,k)
      RHSpsivx(i,j,k) = axpml(i,j,k) * dxixdx(i,j,k) + bxpml(i,j,k) * psivx(i,j,k)

      dxixdx(nx-npmlhor + i,j,k) = dxixdx(nx-npmlhor + i,j,k) + psivx(i+npmlhor,j,k)
      RHSpsivx(i+npmlhor,j,k) = axpml(i+npmlhor,j,k) * dxixdx(nx-npmlhor + i,j,k) + bxpml(i+npmlhor,j,k) * psivx(i+npmlhor,j,k)

     enddo
    enddo
   enddo 

  endif


 ! Compute divergence

 div = dxixdx + dxizdz

 ! RHS continuity and pressure

 !g = 0
 !gradp0_z = 0

! div(:,:,nz) = 0.0
 p = - c2rho0 * div - gradp0_z *xi_z 
 
 call ddx(p, gradp_x, 0)


 if (HORIZONTAL_PMLS) THEN
  
  do k=1,nz
   do j=1,dim2(rank)
     do i=1,npmlhor

      gradp_x(i,j,k) = gradp_x(i,j,k) + psi_gradp_x(i,j,k)
      RHS_psi_gradp_x(i,j,k) = axpml(i,j,k) * gradp_x(i,j,k) + bxpml(i,j,k) * psi_gradp_x(i,j,k)

      gradp_x(nx-npmlhor + i,j,k) = gradp_x(nx-npmlhor+ i,j,k) + psi_gradp_x(i+npmlhor,j,k)
      RHS_psi_gradp_x(i+npmlhor,j,k) = axpml(i+npmlhor,j,k) * &
      gradp_x(nx-npmlhor+i,j,k) + bxpml(i+npmlhor,j,k) * psi_gradp_x(i+npmlhor,j,k)

     enddo
    enddo
   enddo 
  
 endif

 call ddz(p, gradp_z, 2)

 pmlindex = 0
 do k=1,nz

  if ((k .le. npmlbot) .or. (k .ge. (nz-npmltop+1))) then
   pmlindex = pmlindex + 1

   RHSpsip(:,:,pmlindex) = az(:,:,k)*gradp_z(:,:,k) + bzpml(:,:,k)*psip(:,:,pmlindex)    
   gradp_z(:,:,k) =  gradp_z(:,:,k) + psip(:,:,pmlindex)
  endif 

 enddo


! RHS v_x

 RHSv_x = - rhoinv * gradp_x 

! RHS v_z and PML accounting for gradp_z

 rho =  - rho0 * div - gradrho0_z * xi_z
 do k= 1,nz
  RHSv_z(:,:,k) =  - rhoinv(:,:,k) * (gradp_z(:,:,k) + rho(:,:,k)*g(k))
 enddo


 RHSxi_x = v_x
 RHSxi_z = v_z

 if (time .le. (forcing_length)*cadforcing_step) then 
  if (compute_adjoint) then

    RHSv_z(:,:,o_rad) = RHSv_z(:,:,o_rad) + forcing * delta_width/rho0(:,:,o_rad)

  elseif (compute_forward) then
 
    RHSv_z(:,:,e_rad) = RHSv_z(:,:,e_rad) + vr(:,:,1) * fwdsource(time+1) * &
			  delta_width/rho0(:,:,e_rad)
  endif
 endif

! if (damp_waves) call damp_velocity

 if (BACKGROUND_FLOWS_EXIST) then

  call ddx(v_x,dxixdx,1)
  RHSv_x = RHSv_x - 2. * v0_x * dxixdx

  call ddz(v_x,dxixdx,1)
  RHSv_x = RHSv_x - 2. * v0_z * dxixdx
 
  call ddx(v_z,dxixdx,1)
  RHSv_z = RHSv_z - 2. * v0_x * dxixdx

  call ddz(v_z,dxixdx,1)
  RHSv_z = RHSv_z - 2. * v0_z * dxixdx

 endif


! RHSv_x(:,:,nz) = 0.
! RHSv_z(:,:,nz) = 0.
 RHSv_x(:,:,1) = 0.
 RHSv_z(:,:,1) = 0.
 RHSxi_x(:,:,1) = 0.
!! RHSxi_x(:,:,nz) = 0.
 RHSxi_z(:,:,1) = 0.
! RHSxi_z(:,:,nz) = 0.


 if (HORIZONTAL_PMLS) then
  scr(1,:,:,:) = 0.0
  scr(nx,:,:,:) = 0.0
 endif

 END SUBROUTINE MP_QUIET_PML_2D_DISPL



!================================================================================================
SUBROUTINE MP_QUIET_SPONGE_2D

 ! CONSERVATIVE FORM OF THE FLUXES.

 ! I can't really do the conservative form - doesn't work for the incredibly small density 
 ! and forcing functions anyway. I don't know what the error is anyway
 implicit none
 integer k,i

 call ddx(v_x, dvxdx, 1)
 call ddz(v_z, dvzdz, 1)

 call ddx(p, gradp_x, 1)
 call ddz(p, gradp_z, 1)

 div = dvxdx + dvzdz
 
 RHScont = - gradrho0_z * v_z - rho0 * div - spongexyz*rho            

 RHSv_x = - rhoinv * gradp_x - spongexyz * v_x

 RHSv_z = - rhoinv * gradp_z - spongexyz * v_z

 do k= 1,nz
   RHSv_z(:,:,k) = RHSv_z(:,:,k) - rho(:,:,k)*g(k)*rhoinv(:,:,k)
   if (ABS(k-e_rad) .LE. 4) &
    RHSv_z(:,:,k) = RHSv_z(:,:,k) +  forcing*source_dep(k)
 enddo

 RHSp = - c2rho0 * div - v_z * gradp0_z  - spongexyz * p

 RHSv_x(:,:,nz) = 0.
 RHSv_z(:,:,nz) = 0.
 RHSv_x(:,:,1) = 0.
 RHSv_z(:,:,1) = 0.
 RHSp(:,:,1) = 0.
 RHSp(:,:,nz) = 0.

 END SUBROUTINE MP_QUIET_SPONGE_2D


!================================================================================================

SUBROUTINE PML_BOUNDARY_ACCOUNTING_2D
 implicit none
 integer pmlindex, k
 real*8 tempfilt(nx, dim2(rank), npmltop)
 
 pmlindex = 0
 do k=1,nz

  if ((k .le. npmlbot) .or. (k .ge. (nz-npmltop+1))) then
   pmlindex = pmlindex + 1

   RHSpsip(:,:,pmlindex) = az(:,:,k)*gradp_z(:,:,k) + bzpml(:,:,k)*psip(:,:,pmlindex)    
   gradp_z(:,:,k) =  gradp_z(:,:,k) + psip(:,:,pmlindex)

    RHSpsivz(:,:,pmlindex) = az(:,:,k)*dvzdz(:,:,k) + bzpml(:,:,k)*psivz(:,:,pmlindex)    
    dvzdz(:,:,k) = dvzdz(:,:,k) + psivz(:,:,pmlindex)
  endif 

 enddo

 if (magnetic) then

  pmlindex = 0 
  do k=2,npmlbot
   pmlindex = pmlindex + 1
   tempfilt(:,:,pmlindex) = 0.5 * (dvzdz(:,:,k-1) + dvzdz(:,:,k+1))
  enddo

  dvzdz(:,:,2:npmlbot) = tempfilt(:,:,1:npmlbot-1)
 
  pmlindex = 0 
  do k=nz-npmltop+1,nz-1
   pmlindex = pmlindex + 1
   tempfilt(:,:,pmlindex) = 0.5 * (dvzdz(:,:,k-1) + dvzdz(:,:,k+1))
  enddo

  dvzdz(:,:,nz-npmltop+1:nz-1) = tempfilt(:,:,1:npmltop-1)
 

  pmlindex=0
  do k=1,nz

   if ((k .le. npmlbot) .or. (k .ge. (nz-npmltop+1))) then
    pmlindex = pmlindex + 1
    RHSpsiinductionbx(:,:,pmlindex) = az(:,:,k)*dzflux2(:,:,k) &
			+ bzpml(:,:,k)*psiinductionbx(:,:,pmlindex)    

    dzflux2(:,:,k) = dzflux2(:,:,k) + psiinductionbx(:,:,pmlindex)

    RHSpsidzbx(:,:,pmlindex) = az(:,:,k)*dzbx(:,:,k) &
			+ bzpml(:,:,k)*psidzbx(:,:,pmlindex)    

    dzbx(:,:,k) = dzbx(:,:,k) + psidzbx(:,:,pmlindex)

   endif 
 
  enddo

 endif

END SUBROUTINE PML_BOUNDARY_ACCOUNTING_2D

!================================================================================================

END MODULE MP_PHYSICS_2D

