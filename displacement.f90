Module DISPLACEMENT 

! --------------------------------------------------------------------------
! MPI Version of the Cartesian Solar Wave Simulator.
! Copyright 2008, Shravan Hanasoge
                                                                                                                                                         
! W. W. Hansen Experimental Physics Laboratory
! Stanford University, Stanford, CA 94305, USA
! Email: shravan@solar.stanford.edu
! --------------------------------------------------------------------------

! SUBROUTINES IN THIS MODULE:
!

  use initialize
  use derivatives
  use mp_physics
  use damping
  implicit none

Contains

!============================================================================


SUBROUTINE MP_QUIET_SPONGE_3D_DISPL

 ! CONSERVATIVE FORM OF THE FLUXES.

 ! I can't really do the conservative form - doesn't work for the incredibly small density 
 ! and forcing functions anyway. I don't know what the error is anyway
 implicit none
 integer k,i,j,bc 

 bc = 1

 call ddxyz(xi_x, dxixdx, xi_y, dxiydy, xi_z, dxizdz, bc)
 div = dxixdx + dxiydy + dxizdz

 p = - c2rho0 * div - xi_z * gradp0_z
 rho = - rho0 * div - xi_z * gradrho0_z 

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

 RHSv_x = - rhoinv * gradp_x - spongexyz * v_x

 RHSv_y = - rhoinv * gradp_y - spongexyz * v_y

 RHSv_z = - rhoinv * gradp_z - spongexyz * v_z

 do k= 1,nz
   RHSv_z(:,:,k) = RHSv_z(:,:,k) - rho(:,:,k)*g(k)*rhoinv(:,:,k)
 enddo

 RHSxi_x = v_x
 RHSxi_y = v_y
 RHSxi_z = v_z


 if (.not. kernel_mode) then

  do k=e_rad-3,e_rad+3
    RHSv_z(:,:,k) = RHSv_z(:,:,k) +  forcing*source_dep(k)
  enddo

 elseif (kernel_mode) then

   if (compute_adjoint) then

     do j=1,dim2(rank)
      RHSv_z(:,j,o_rad) = RHSv_z(:,j,o_rad) + cos(4.*pi*(time*timestep-650.)/700.)  * &
      exp(-powerspec_fac*(time*timestep-650.)**2./(2.*200.**2.)) &
      * exp(-(x - loc_x)**2./(2. * (x(2)-x(1)))**2.)* &
      exp(-(y(j+ystart(rank)-1) - loc_y)**2./(2. * (x(2)-x(1)))**2.)
     enddo

   elseif (compute_forward .and. (time .le. maxtime/2)) then

     RHSv_z(:,:,e_rad) = RHSv_z(:,:,e_rad) +  forcing * delta_width
   endif

  endif


 do k=1,6
  scr(:,:,nz,k) = 0.0
  scr(:,:,1,k) = 0.0
 enddo


 if (DAMP_WAVES) call DAMP_VELOCITY

 END SUBROUTINE MP_QUIET_SPONGE_3D_DISPL

!================================================================================================

SUBROUTINE MP_MHD_SPONGE_3D_DISPL

 ! CONSERVATIVE FORM OF THE FLUXES.

 ! I can't really do the conservative form - doesn't work for the incredibly small density 
 ! and forcing functions anyway. I don't know what the error is anyway
 implicit none
 integer i,j,k, bc
 real(kind=real64), dimension(nx, dim2(rank), nz) :: temp1, temp2, temp3, flux1, flux2, flux3

  bc = 1
 call ddxyz(xi_x, dxixdx, xi_y, dxiydy, xi_z, dxizdz, bc)
 div = dxixdx + dxiydy + dxizdz

 p = - c2rho0 * div - xi_z * gradp0_z - xi_x * gradp0_x - xi_y * gradp0_y
 rho = - rho0 * div - xi_z * gradrho0_z - xi_x * gradrho0_x - xi_y * gradrho0_y

 flux1 = - (xi_z * boy - xi_y * boz)
 flux2 = - (xi_x * boz - xi_z * box)
 flux3 = - (xi_y * box - xi_x * boy)

! call ddz(v_z, dvzdz, 1)
! call ddz(p, gradp_z, 1)

! call ddz(flux2, dzflux2, 1)
! call ddz(flux1, dzflux1, 1)

! call ddz(bx, dzbx, 1)
! call ddz(by, dzby, 1)


 ! CALLING OVERLAP DERIVATIVE ROUTINE
 
 ! NOTE THAT THE ORDER IS SENSITIVE.
 ! MUST FULLY COMPUTE bx, by, bz BEFORe
 ! COMPUTING DERIVATIVES THEREOF
 do j=1,dim2(rank)
  do i=1,nx
   gradp_z(i,j,1)  =   (- c_speed(i,j,1)*rho0(i,j,1)*dvzdz(i,j,1)   &
            & - rho(i,j,1)*g(1))*unstretch(1)
   gradp_z(i,j,nz)  =  (c_speed(i,j,nz)*rho0(i,j,nz)*dvzdz(i,j,nz)     &
                & - rho(i,j,nz)*g(nz))*unstretch(nz)
  enddo
 enddo


 bc =1
 call ddxyz(flux3, by, p, gradp_y, flux2, dzflux2, bc)
 call ddxyz(p, gradp_x, flux3, bx, flux1, dzflux1, bc)
 bc =4
 call ddxyz(flux2, bz, flux1, flux3, p, gradp_z, bc)

 bx =   bx - dzflux2
 by = - by + dzflux1
 bz =   bz - flux3

 ! NOTE THAT flux3 is overwritten in the third ddxyz call

 bc =1
 call ddxyz(by, curlbz, bz, curlbx,  bx, dzbx, bc)
 call ddxyz(bz, curlby, bx, flux1, by, dzby,bc)


 
! call ddx(flux3, RHSb_y, 1)
! call ddy(flux3, RHSb_x, 1)

! call ddx(flux2, RHSb_z,1)
! call ddy(flux1, flux2, 1)


! call ddx(v_x, dvxdx, 1)
! call ddx(p, gradp_x, 1)

! call ddy(v_y, dvydy, 1)
! call ddy(p, gradp_y, 1)

! call ddy(bz, curlbx, 1)
 curlbx = curlbx - dzby

! call ddx(bz, curlby, 1)
 curlby = -curlby + dzbx

! call ddx(by, curlbz, 1)
! call ddy(bx, flux1, 1)
 curlbz = curlbz - flux1


 call cross(curlbx, curlby, curlbz, box, boy, boz, flux1, flux2, flux3)

 RHSv_x = rhoinv * (- gradp_x + flux1)
 RHSv_y = rhoinv * (- gradp_y + flux2)
 RHSv_z = rhoinv * (- gradp_z + flux3)


 call cross(curlbox, curlboy, curlboz, bx, by, bz, flux1, flux2, flux3)


 RHSv_x = RHSv_x + rhoinv*flux1 - spongexyz * v_x
 RHSv_y = RHSv_y + rhoinv*flux2 - spongexyz * v_y
 RHSv_z = RHSv_z + rhoinv*flux3 - spongexyz * v_z


 do k= 1,nz
   RHSv_z(:,:,k) = RHSv_z(:,:,k) - rhoinv(:,:,k) * rho(:,:,k)*g(k)

!   if (abs(k-e_rad) .LE. 4) &
!    RHSv_z(:,:,k) = RHSv_z(:,:,k) + forcing*source_dep(k) 

 enddo 

 RHSxi_x = v_x
 RHSxi_y = v_y
 RHSxi_z = v_z

 if (.not. kernel_mode) then

  do k=e_rad-3,e_rad+3
    RHSv_z(:,:,k) = RHSv_z(:,:,k) +  forcing*source_dep(k)
  enddo

 elseif (kernel_mode) then

   if (compute_adjoint) then

     do j=1,dim2(rank)
      RHSv_z(:,j,o_rad) = RHSv_z(:,j,o_rad) + cos(4.*pi*(time*timestep-650.)/700.) *&
       source_dep(k) * exp(-powerspec_fac*(time*timestep-650.)**2./(2.*200.**2.)) &
      * exp(-(x - loc_x)**2./(2. * (x(2)-x(1)))**2.)  !* exp(-(y(j+ystart(rank)-1) - loc_y)**2./(2. * (x(2)-x(1)))**2.)
     enddo

!
     !RHSv_z(:,:,o_rad) = RHSv_z(:,:,o_rad) +  forcing/rho0(:,:,o_rad) * delta_width
   elseif (compute_forward .and. (time .le. maxtime/2)) then

     RHSv_z(:,:,e_rad) = RHSv_z(:,:,e_rad) +  forcing/rho0(:,:,e_rad) * delta_width

   endif

  endif

 do k=1,6
  scr(:,:,nz,k) = 0.
  scr(:,:,1,k) = 0.
 enddo


  if (DAMP_WAVES) call DAMP_VELOCITY


 END SUBROUTINE MP_MHD_SPONGE_3D_DISPL


!================================================================================================


 Subroutine MP_FLOWS_SPONGE_3D_DISPL ()

  ! Right hand side computation

  implicit none
  integer i,j,k, bc
  real(kind=real64), allocatable, dimension(:,:,:) :: temp_x, temp_y, temp_z

  allocate(temp_x(nx,dim2(rank),nz),temp_y(nx,dim2(rank),nz),&
   temp_z(nx,dim2(rank),nz))


  ! ONLY THE LINEAR ADVECTION TERM IS ACCOUNTED FOR. MORE TERMS
  ! NEED TO BE ADDED TO ENSURE CONSISTENCY BUT FOR ALL PRACTICAL
  ! PURPOSES, I DEEM THIS APPROXIMATION SUFFICIENT

  bc =1
  call ddxyz(xi_x, dxixdx, xi_y, dxiydy, xi_z, dxizdz, bc)
  div = dxixdx + dxiydy + dxizdz

  p = - c2rho0 * div - xi_z * gradp0_z
  rho = - rho0 * div - xi_z * gradrho0_z 

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


! call ddxyz(p, gradp_x, p, gradp_y, p, gradp_z)

 bc =1
 call ddxyz(v_x, temp_x, v_x, temp_y, v_x, temp_z, bc)

 RHSv_x = - 2.0 * (v0_x * temp_x + v0_y * temp_y + v0_z * temp_z) - &
    rhoinv * gradp_x - spongexyz * v_x


 call ddxyz(v_y, temp_x, v_y, temp_y, v_y, temp_z,bc)

 RHSv_y = - 2.0 * (v0_x * temp_x + v0_y * temp_y + v0_z * temp_z) - &
    rhoinv * gradp_y - spongexyz * v_y


 call ddxyz(v_z, temp_x, v_z, temp_y, v_z, temp_z,bc)

 RHSv_z = - 2.0 * (v0_x * temp_x + v0_y * temp_y + v0_z * temp_z) - &
    rhoinv * gradp_z - spongexyz * v_z


 deallocate(temp_x, temp_y, temp_z)

 do k = 1,nz

     RHSv_z(:,:,k)  = RHSv_z(:,:,k) - rhoinv(:,:,k) *  rho(:,:,k) * g(k)

!     if (ABS(k-e_rad) .LE. 4) then
 !       RHSv_z(:,:,k) = RHSv_z(:,:,k) +  forcing*source_dep(k)
  !   endif

  enddo
  
  RHSxi_x = v_x
  RHSxi_y = v_y
  RHSxi_z = v_z


 if (.not. kernel_mode) then

  do k=e_rad-3,e_rad+3
    RHSv_z(:,:,k) = RHSv_z(:,:,k) +  forcing*source_dep(k)
  enddo

 elseif (kernel_mode) then

   if (compute_adjoint) then
     do j=1,dim2(rank)
      RHSv_z(:,j,o_rad) = RHSv_z(:,j,o_rad) + cos(4.*pi*(time*timestep-650.)/700.)  * &
      exp(-powerspec_fac*(time*timestep-650.)**2./(2.*200.**2.)) &
      * exp(-(x - loc_x)**2./(2. * (x(2)-x(1)))**2.)* &
      exp(-(y(j+ystart(rank)-1) - loc_y)**2./(2. * (x(2)-x(1)))**2.)
!*source_dep(k)

     !RHSv_z(:,:,o_rad) = RHSv_z(:,:,o_rad) +  forcing/rho0(:,:,o_rad) * delta_width
     enddo

   elseif (compute_forward .and. (time .le. maxtime/2)) then
     RHSv_z(:,:,e_rad) = RHSv_z(:,:,e_rad) +  forcing/rho0(:,:,e_rad) * delta_width
   endif
  endif


 do k=1,6
  scr(:,:,nz,k) = 0.
  scr(:,:,1,k) = 0.
 enddo

  if (DAMP_WAVES) call DAMP_VELOCITY

 end SUBROUTINE MP_FLOWS_SPONGE_3D_DISPL

!================================================================================================

 Subroutine MP_FLOWS_PML_3D_DISPL ()

  ! Right hand side computation

  implicit none
  integer i,j,k,pmlindex, bc
  real(kind=real64), allocatable, dimension(:,:,:) :: temp_x, temp_y, temp_z

  allocate(temp_x(nx,dim2(rank),nz),temp_y(nx,dim2(rank),nz),&
   temp_z(nx,dim2(rank),nz))


  ! ONLY THE LINEAR ADVECTION TERM IS ACCOUNTED FOR. MORE TERMS
  ! NEED TO BE ADDED TO ENSURE CONSISTENCY BUT FOR ALL PRACTICAL
  ! PURPOSES, I DEEM THIS APPROXIMATION SUFFICIENT

  ! NOTE THAT I ASSUME THAT v0 = 0.0 in the PML ZONE!


 bc =1
 call ddxyz(xi_x, dxixdx, xi_y, dxiydy, xi_z, dxizdz, bc)

 pmlindex = 0
 do k=1,nz

  if ((k .le. npmlbot) .or. (k .ge. (nz-npmltop+1))) then
   pmlindex = pmlindex + 1

   RHSpsivz(:,:,pmlindex) = az(:,:,k)*dxizdz(:,:,k) + bzpml(:,:,k)*psivz(:,:,pmlindex)    
   dxizdz(:,:,k) = dxizdz(:,:,k) + psivz(:,:,pmlindex)
  endif 

 enddo

 div = dxixdx + dxiydy + dxizdz

 p = - c2rho0 * div - xi_z * gradp0_z
 rho = - rho0 * div - xi_z * gradrho0_z 


 bc = 0

 call ddxyz(p, gradp_x, p, gradp_y, p, gradp_z, bc)

! call ddxyz(p, gradp_x, p, gradp_y, p, gradp_z)

 pmlindex = 0
 do k=1,nz

  if ((k .le. npmlbot) .or. (k .ge. (nz-npmltop+1))) then
   pmlindex = pmlindex + 1

   RHSpsip(:,:,pmlindex) = az(:,:,k)*gradp_z(:,:,k) + bzpml(:,:,k)*psip(:,:,pmlindex)    
   gradp_z(:,:,k) =  gradp_z(:,:,k) + psip(:,:,pmlindex)
  endif 

 enddo
 
 bc =1
 call ddxyz(v_x, temp_x, v_x, temp_y, v_x, temp_z, bc)

 RHSv_x = - 2.0 * (v0_x * temp_x + v0_y * temp_y + v0_z * temp_z) - &
    rhoinv * gradp_x 

 bc =1
 call ddxyz(v_y, temp_x, v_y, temp_y, v_y, temp_z, bc)

 RHSv_y = - 2.0 * (v0_x * temp_x + v0_y * temp_y + v0_z * temp_z) - &
    rhoinv * gradp_y 

 bc =1
 call ddxyz(v_z, temp_x, v_z, temp_y, v_z, temp_z, bc)

 RHSv_z = - 2.0 * (v0_x * temp_x + v0_y * temp_y + v0_z * temp_z) - &
    rhoinv * gradp_z 

 deallocate(temp_x, temp_y, temp_z)

 do k = 1,nz

     RHSv_z(:,:,k)  = RHSv_z(:,:,k) - rhoinv(:,:,k) *  rho(:,:,k) * g(k)

     if (ABS(k-e_rad) .LE. 4) then
        RHSv_z(:,:,k) = RHSv_z(:,:,k) +  forcing*source_dep(k)
     endif

  enddo
  
  RHSxi_x = v_x
  RHSxi_y = v_y
  RHSxi_z = v_z

 do k=1,6
  scr(:,:,nz,k) = 0.
  scr(:,:,1,k) = 0.
 enddo


  if (DAMP_WAVES) call DAMP_VELOCITY

 end SUBROUTINE MP_FLOWS_PML_3D_DISPL

!================================================================================================


SUBROUTINE MP_QUIET_PML_3D_DISPL

 ! CONSERVATIVE FORM OF THE FLUXES.

 ! I can't really do the conservative form - doesn't work for the incredibly small density 
 ! and forcing functions anyway. I don't know what the error is anyway
 implicit none
 integer k,i,j,pmlindex, bc

 bc =1
 call ddxyz(xi_x, dxixdx, xi_y, dxiydy, xi_z, dxizdz, bc)


 pmlindex = 0
 do k=1,nz

  if ((k .le. npmlbot) .or. (k .ge. (nz-npmltop+1))) then
  !if (k .le. npmlbot) then
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

   if (PROC_HAS_PML) then
      dxiydy = dxiydy + psivy
      RHSpsivy = aypml * dxiydy + bypml * psivy
   endif
 
  endif


 div = dxixdx + dxiydy + dxizdz

 p = - c2rho0 * div - xi_z * gradp0_z
 rho = - rho0 * div - xi_z * gradrho0_z 

! if (horizontal_pmls) then
!  p(1,:,:) =0.0
!  p(nx,:,:) = 0.0

!  if (ystart(rank) .eq. 1) p(:,1,:) = 0.0
!  if (ystart(rank) + dim2(rank)-1 .eq. ny) p(:,dim2(rank), :)=0.0

! endif

 bc = 0
 call ddxyz(p, gradp_x, p, gradp_y, p, gradp_z, bc)

 pmlindex = 0
 do k=1,nz

  if ((k .le. npmlbot) .or. (k .ge. (nz-npmltop+1))) then
  !if (k .le. npmlbot) then
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
      RHS_psi_gradp_x(i+npmlhor,j,k) = axpml(i+npmlhor,j,k) * gradp_x(nx-npmlhor+i,j,k) &
      + bxpml(i+npmlhor,j,k) * psi_gradp_x(i+npmlhor,j,k)

     enddo
    enddo
   enddo 

   if (PROC_HAS_PML) then

      gradp_y = gradp_y + psi_gradp_y
      RHS_psi_gradp_y= aypml * gradp_y + bypml * psi_gradp_y

   endif

 endif 
 
 RHSv_x = - rhoinv * gradp_x  

 RHSv_y = - rhoinv * gradp_y 

 RHSv_z = - rhoinv * gradp_z 

 do k= 1,nz
   RHSv_z(:,:,k) = RHSv_z(:,:,k) - rho(:,:,k)*g(k)*rhoinv(:,:,k)
 enddo

! do k=nz-40,nz
!   RHSv_x(:,:,k) = RHSv_x(:,:,k) - spongexyz(:,:,k) * v_x(:,:,k)
!   RHSv_y(:,:,k) = RHSv_y(:,:,k) - spongexyz(:,:,k) * v_y(:,:,k)
!   RHSv_z(:,:,k) = RHSv_z(:,:,k) - spongexyz(:,:,k) * v_z(:,:,k)
! enddo

 RHSxi_x = v_x
 RHSxi_y = v_y
 RHSxi_z = v_z

 if (.not. kernel_mode) then

  do k=e_rad-3,e_rad+3
    RHSv_z(:,:,k) = RHSv_z(:,:,k) +  forcing*source_dep(k)
  enddo

 elseif (kernel_mode) then

   if (compute_adjoint) then

   !  do j=1,dim2(rank)
   !   RHSv_z(:,j,o_rad) = RHSv_z(:,j,o_rad) + cos(4.*pi*(time*timestep-450.)/700.)  * &
!       exp(-powerspec_fac*(time*timestep-450.)**2./(2.*200.**2.))  * &
!                          exp(-(x - loc_x)**2./(2. * (x(2)-x(1)))**2.) * &
!       exp(-(y(j+ystart(rank)-1) - loc_y)**2./(2. * (x(2)-x(1)))**2.) * &
!       delta_width/rho0(:,j,e_rad)
!     enddo

     RHSv_z(:,:,o_rad) = RHSv_z(:,:,o_rad) +  forcing/rho0(:,:,o_rad) * delta_width

   elseif (compute_forward) then

     RHSv_z(:,:,e_rad) = RHSv_z(:,:,e_rad) +  forcing/rho0(:,:,e_rad) * delta_width

   endif
  endif


 if (DAMP_WAVES) call DAMP_VELOCITY

 do k=1,dimfive
  scr(:,:,nz,k) = 0.0
  scr(:,:,1,k) = 0.0
 enddo


 if (HORIZONTAL_PMLS) then

  scr(1,:,:,:) = 0.0
  scr(nx,:,:,:) = 0.0

  if (ystart(rank) .eq. 1) &
    scr(:,1,:,:) = 0.0
  if (ystart(rank) + dim2(rank) -1 .eq. ny) &
    scr(:,dim2(rank),:,:) = 0.0

  endif


 END SUBROUTINE MP_QUIET_PML_3D_DISPL

!================================================================================================


SUBROUTINE MP_MHD_PML_3D_DISPL

 ! CONSERVATIVE FORM OF THE FLUXES.

 ! I can't really do the conservative form - doesn't work for the incredibly small density 
 ! and forcing functions anyway. I don't know what the error is anyway
 implicit none
 integer i,j,k, pmlindex,bc
 real(kind=real64), dimension(nx, dim2(rank), nz) :: flux1, flux2, flux3

 flux1 = - (xi_z * boy - xi_y * boz)
 flux2 = - (xi_x * boz - xi_z * box)
 flux3 = - (xi_y * box - xi_x * boy)

 bc =1
 call ddxyz(xi_x, dxixdx, xi_y, dxiydy, xi_z, dxizdz, bc)

 pmlindex = 0
 do k=1,nz

  !if ((k .le. npmlbot) .or. (k .ge. (nz-npmltop+1))) then
  if (k .le. npmlbot) then
   pmlindex = pmlindex + 1

   RHSpsivz(:,:,pmlindex) = az(:,:,k)*dxizdz(:,:,k) + bzpml(:,:,k)*psivz(:,:,pmlindex)    
   dxizdz(:,:,k) = dxizdz(:,:,k) + psivz(:,:,pmlindex)
  endif 

 enddo

 div = dxixdx + dxiydy + dxizdz

 p = - c2rho0 * div - xi_z * gradp0_z - xi_x * gradp0_x - xi_y * gradp0_y
 rho = - rho0 * div - xi_z * gradrho0_z - xi_x * gradrho0_x - xi_y * gradrho0_y

 bc =0
 call ddxyz(p, gradp_x, p, gradp_y, p, gradp_z, 0)

 pmlindex = 0
 do k=1,nz

  !if ((k .le. npmlbot) .or. (k .ge. (nz-npmltop+1))) then
  if (k .le. npmlbot) then
   pmlindex = pmlindex + 1

   RHSpsip(:,:,pmlindex) = az(:,:,k)*gradp_z(:,:,k) + bzpml(:,:,k)*psip(:,:,pmlindex)    
   gradp_z(:,:,k) =  gradp_z(:,:,k) + psip(:,:,pmlindex)
  endif 

 enddo

 bc =1
 call ddxyz(flux3, by, flux3, bx, flux2, dzflux2, bc)
 call ddxyz(flux2, bz, flux1, flux3, flux1, dzflux1, bc)


 ! NOTE THAT flux3 is overwritten in the third ddxyz call

 ! THE PML STUFF HAS TO COME HERE UNFORTUNATELY

 pmlindex=0
 do k=1,nz

  !if ((k .le. npmlbot) .or. (k .ge. (nz-npmltop+1))) then
  if (k .le. npmlbot) then
   pmlindex = pmlindex + 1

   RHSpsiinductionbx(:,:,pmlindex) = az(:,:,k)*dzflux2(:,:,k) &
      + bzpml(:,:,k)*psiinductionbx(:,:,pmlindex)    

   dzflux2(:,:,k) = dzflux2(:,:,k) + psiinductionbx(:,:,pmlindex)


   RHSpsiinductionby(:,:,pmlindex) = az(:,:,k)*dzflux1(:,:,k) &
      + bzpml(:,:,k)*psiinductionby(:,:,pmlindex)    

   dzflux1(:,:,k) = dzflux1(:,:,k) + psiinductionby(:,:,pmlindex)

  endif 
 
 enddo

 bx =   bx - dzflux2
 by = - by + dzflux1
 bz =   bz - flux3

 bc =1
 call ddxyz(by, curlbz, bz, curlbx,  bx, dzbx, bc)
 call ddxyz(bz, curlby, bx, flux1, by, dzby, bc)

 ! NOTE THAT FLUX1 IS OVERWRITTEN IN THE FINAL CALL

 ! THE PML STUFF HAS TO COME HERE UNFORTUNATELY

 pmlindex=0
 do k=1,nz
  !if ((k .le. npmlbot) .or. (k .ge. (nz-npmltop+1))) then
  if (k .le. npmlbot) then
    pmlindex = pmlindex + 1

    RHSpsidzbx(:,:,pmlindex) = az(:,:,k)*dzbx(:,:,k) &
      + bzpml(:,:,k)*psidzbx(:,:,pmlindex)    

    dzbx(:,:,k) = dzbx(:,:,k) + psidzbx(:,:,pmlindex)

    RHSpsidzby(:,:,pmlindex) = az(:,:,k)*dzby(:,:,k) &
      + bzpml(:,:,k)*psidzby(:,:,pmlindex)    

    dzby(:,:,k) = dzby(:,:,k) + psidzby(:,:,pmlindex)
   endif
 enddo


 curlbx = curlbx - dzby
 curlby = -curlby + dzbx
 curlbz = curlbz - flux1

 call cross(curlbx, curlby, curlbz, box, boy, boz, flux1, flux2, flux3)

 RHSv_x = rhoinv * (- gradp_x + flux1)
 RHSv_y = rhoinv * (- gradp_y + flux2)
 RHSv_z = rhoinv * (- gradp_z + flux3)


 call cross(curlbox, curlboy, curlboz, bx, by, bz, flux1, flux2, flux3)


 RHSv_x = RHSv_x + rhoinv*flux1 
 RHSv_y = RHSv_y + rhoinv*flux2 
 RHSv_z = RHSv_z + rhoinv*flux3

 do k=nz-35,nz
  RHSv_x(:,:,k) = RHSv_x(:,:,k) - spongexyz(:,:,k) * v_x(:,:,k) 
  RHSv_y(:,:,k) = RHSv_y(:,:,k) - spongexyz(:,:,k) * v_y(:,:,k) 
  RHSv_z(:,:,k) = RHSv_z(:,:,k) - spongexyz(:,:,k) * v_z(:,:,k) 
 enddo

 do k= 1,nz
   RHSv_z(:,:,k) = RHSv_z(:,:,k) - rhoinv(:,:,k) * rho(:,:,k)*g(k)
 enddo 

 RHSxi_x = v_x
 RHSxi_y = v_y
 RHSxi_z = v_z


 if (DAMP_WAVES) call DAMP_VELOCITY

 if (.not. kernel_mode) then

 elseif (kernel_mode .and. (time .le. forcing_length*cadforcing_step)) then

   ! WILSON DEPRESSION TAKEN INTO ACCOUNT
   ! ASSUME SOURCES ARE UNIFORMLY DISTRIBUTED HORIZONTALLY
   ! BUT NOT AT THE SAME GEOMETRICAL HEIGHT
   ! SAME OPTICAL DEPTH (i.e., SAME DENSITY)

   if (compute_adjoint) then

    do j=1,dim2(rank)
     do i=1,nx
      RHSv_z(i,j,orad_2d(i,j)) = RHSv_z(i,j,orad_2d(i,j)) +  &
  forcing(i,j)/rho0(i,j,orad_2d(i,j)) * delta_width_2d(i,j) 
     enddo
    enddo

   elseif (compute_forward) then

     do j=1,dim2(rank)
      do i=1,nx
       RHSv_z(i,j,erad_2d(i,j)) = RHSv_z(i,j,erad_2d(i,j)) +  &
    forcing(i,j)/rho0(i,j,erad_2d(i,j)) * delta_width_2d(i,j)
      enddo
     enddo

   endif
  endif

 do k=1,6
  scr(:,:,nz,k) = 0.
  scr(:,:,1,k) = 0.
 enddo




 END SUBROUTINE MP_MHD_PML_3D_DISPL


!================================================================================================

END MODULE DISPLACEMENT
