MODULE PML

 use initialize
 use derivatives
 use all_modules
 use damping
 implicit none

Contains

!==================================================================================
SUBROUTINE MP_QUIET_PML_3D

 ! CONSERVATIVE FORM OF THE FLUXES.

 ! I can't really do the conservative form - doesn't work for the incredibly small density 
 ! and forcing functions anyway. I don't know what the error is anyway
 implicit none
 integer k,i,j,pmlindex


 ! CALL THE OVERLAP-DERIVATIVE ROUTINE

 call ddxyz(v_x, dvxdx, v_y, dvydy, v_z, dvzdz, 1)

 call ddxyz(p, gradp_x, p, gradp_y, p, gradp_z, 0)



 !call ddx(p, gradp_x, 1)
 !call ddy(p, gradp_y, 1)
 !call ddz(p, gradp_z, 1)

 ! NOTE THAT flux3 is overwritten in the final ddxyz call
 
 call PML_BOUNDARY_ACCOUNTING_3D ()
  
 div = dvxdx + dvydy + dvzdz

 RHScont = - gradrho0_z * v_z - rho0 * div

 RHSv_x = - rhoinv * gradp_x 

 RHSv_y = - rhoinv * gradp_y

 do k= 1,nz
    RHSv_z(:,:,k) =  - rhoinv(:,:,k) * (gradp_z(:,:,k) + rho(:,:,k)*g(k))
 
   if (ABS(k-e_rad) .LE. 4) &
    RHSv_z(:,:,k) = RHSv_z(:,:,k) +  forcing*source_dep(k)
 enddo

 RHSp = - c2rho0 * div - v_z * gradp0_z
 

 do k=1,5
  scr(:,:,1,k) = 0.
  scr(:,:,nz,k) = 0.
 enddo


 if (HORIZONTAL_PMLS) then

  scr(1,:,:,:) = 0.0
  scr(nx,:,:,:) = 0.0

  if (ystart(rank) .eq. 1) &
    scr(:,1,:,:) = 0.0
  if (ystart(rank) + dim2(rank) -1 .eq. nx) &
    scr(:,dim2(rank),:,:) = 0.0

 endif

  if (DAMP_WAVES) call DAMP_VELOCITY

 END SUBROUTINE MP_QUIET_PML_3D

!================================================================================================

SUBROUTINE MP_MHD_PML_3D

 ! CONSERVATIVE FORM OF THE FLUXES.

 ! I can't really do the conservative form - doesn't work for the incredibly small density 
 ! and forcing functions anyway. I don't know what the error is anyway
 implicit none

 integer i,j,k
 real(kind=real64), dimension(nx,dim2(rank),nz) ::  flux1, flux2, flux3
 
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

 call ddxyz(v_x, dvxdx, v_y, dvydy, v_z, dvzdz, 1)
 call ddxyz(flux3, RHSb_y, p, gradp_y, flux1, dzflux1, 1)
 call ddxyz(p, gradp_x, flux3, RHSb_x, flux2, dzflux2, 1)
 call ddxyz(flux2, RHSb_z, bz, curlbx,bx, dzbx, 1)

 ! FLUX2, FLUX3 ARE FINISHED - I.E., ALL RELEVANT
 ! DERIVATIVES HAVE BEEN COMPUTED. WILL NOW USE
 ! THEM AS SCRATCH ARRAYS

 call ddxyz(by, curlbz, bx, flux2, p, gradp_z, 1)
 call ddxyz(bz, curlby, flux1, flux3, by, dzby,1)

 ! NOTE THAT flux3 and flux2 overwritten in the final ddxyz calls
 
 call PML_BOUNDARY_ACCOUNTING_3D ()

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
! call ddy(bx, flux2, 1)
 curlbz = curlbz - flux2


 ! CONTIUNUITY 
 ! --------------------------------------

 RHScont = - gradrho0_x * v_x - gradrho0_y * v_y &
           - gradrho0_z * v_z - rho0 * div 

 ! PRESSURE

 RHSp = - c2rho0 * div - v_x * gradp0_x - v_y * gradp0_y &
  - v_z * gradp0_z 


 
 call cross(curlbx, curlby, curlbz, box, boy, boz, flux1, flux2, flux3)

 RHSv_x = -rhoinv*(gradp_x - flux1)
 RHSv_y = -rhoinv*(gradp_y - flux2)
 RHSv_z = -rhoinv*(gradp_z - flux3)


 call cross(curlbox, curlboy, curlboz, bx, by, bz, flux1, flux2, flux3)


 RHSv_x = RHSv_x + rhoinv*flux1
 RHSv_y = RHSv_y + rhoinv*flux2
 RHSv_z = RHSv_z + rhoinv*flux3


 do k= 1,nz
   RHSv_z(:,:,k) = RHSv_z(:,:,k) - rhoinv(:,:,k) * rho(:,:,k)*g(k)

   if (abs(k-e_rad) .LE. 4) &
    RHSv_z(:,:,k) = RHSv_z(:,:,k) + forcing*source_dep(k) 

 enddo 

 do k=1,8
  scr(:,:,nz,k) = 0.
  scr(:,:,1,k) = 0.
 enddo


 if (HORIZONTAL_PMLS) then

  scr(1,:,:,:) = 0.0
  scr(nx,:,:,:) = 0.0

  if (ystart(rank) .eq. 1) &
    scr(:,1,:,:) = 0.0
  if (ystart(rank) + dim2(rank) -1 .eq. nx) &
    scr(:,dim2(rank),:,:) = 0.0

 endif


  if (DAMP_WAVES) call DAMP_VELOCITY

 END SUBROUTINE MP_MHD_PML_3D

!================================================================================================


SUBROUTINE PML_BOUNDARY_ACCOUNTING_3D
 implicit none

 integer pmlindex, i, j, k
 real(kind=real64) tempfilt(nx, dim2(rank), npmltop)
 
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

!  pmlindex = 0 
!  do k=2,npmlbot
!   pmlindex = pmlindex + 1
!   tempfilt(:,:,pmlindex) = 0.5 * (dvzdz(:,:,k-1) + dvzdz(:,:,k+1))
!  enddo

!  dvzdz(:,:,2:npmlbot) = tempfilt(:,:,1:npmlbot-1)
 
!  pmlindex = 0 
!  do k=nz-npmltop+1,nz-1
!   pmlindex = pmlindex + 1
!   tempfilt(:,:,pmlindex) = 0.5 * (dvzdz(:,:,k-1) + dvzdz(:,:,k+1))
!  enddo

!  dvzdz(:,:,nz-npmltop+1:nz-1) = tempfilt(:,:,1:npmltop-1)

  pmlindex=0
  do k=1,nz

   if ((k .le. npmlbot) .or. (k .ge. (nz-npmltop+1))) then
    pmlindex = pmlindex + 1

    RHSpsiinductionbx(:,:,pmlindex) = az(:,:,k)*dzflux2(:,:,k) &
      + bzpml(:,:,k)*psiinductionbx(:,:,pmlindex)    

    dzflux2(:,:,k) = dzflux2(:,:,k) + psiinductionbx(:,:,pmlindex)


    RHSpsiinductionby(:,:,pmlindex) = az(:,:,k)*dzflux1(:,:,k) &
      + bzpml(:,:,k)*psiinductionby(:,:,pmlindex)    

    dzflux1(:,:,k) = dzflux1(:,:,k) + psiinductionby(:,:,pmlindex)


    RHSpsidzbx(:,:,pmlindex) = az(:,:,k)*dzbx(:,:,k) &
      + bzpml(:,:,k)*psidzbx(:,:,pmlindex)    

    dzbx(:,:,k) = dzbx(:,:,k) + psidzbx(:,:,pmlindex)

    RHSpsidzby(:,:,pmlindex) = az(:,:,k)*dzby(:,:,k) &
      + bzpml(:,:,k)*psidzby(:,:,pmlindex)    

    dzby(:,:,k) = dzby(:,:,k) + psidzby(:,:,pmlindex)


   endif 
 
  enddo

 endif


 if (HORIZONTAL_PMLS) THEN
  
  do k=1,nz
   do j=1,dim2(rank)
     do i=1,npmlhor

      dvxdx(i,j,k) = dvxdx(i,j,k) + psivx(i,j,k)
      RHSpsivx(i,j,k) = axpml(i,j,k) * dvxdx(i,j,k) + bxpml(i,j,k) * psivx(i,j,k)

      dvxdx(nx-npmlhor + i,j,k) = dvxdx(nx-npmlhor + i,j,k) + psivx(i+npmlhor,j,k)
      RHSpsivx(i+npmlhor,j,k) = axpml(i+npmlhor,j,k) * dvxdx(nx-npmlhor + i,j,k) &
      + bxpml(i+npmlhor,j,k) * psivx(i+npmlhor,j,k)


      gradp_x(i,j,k) = gradp_x(i,j,k) + psi_gradp_x(i,j,k)
      RHS_psi_gradp_x(i,j,k) = axpml(i,j,k) * gradp_x(i,j,k) + bxpml(i,j,k) * psi_gradp_x(i,j,k)

      gradp_x(nx-npmlhor + i,j,k) = gradp_x(nx-npmlhor+ i,j,k) + psi_gradp_x(i+npmlhor,j,k)
      RHS_psi_gradp_x(i+npmlhor,j,k) = axpml(i+npmlhor,j,k) * &
      gradp_x(nx-npmlhor+i,j,k) + bxpml(i+npmlhor,j,k) * psi_gradp_x(i+npmlhor,j,k)


     enddo
    enddo
   enddo 

   if (PROC_HAS_PML) then

      dvydy = dvydy + psivy
      RHSpsivy = aypml * dvydy + bypml * psivy
 
      gradp_y = gradp_y + psi_gradp_y
      RHS_psi_gradp_y= aypml * gradp_y + bypml * psi_gradp_y
 
   endif
 
  endif

END SUBROUTINE PML_BOUNDARY_ACCOUNTING_3D

!================================================================================================



 Subroutine MP_FLOWS_PML_3D ()

  ! Right hand side computation

  implicit none
  
  integer i,j,k
  real(kind=real64), allocatable, dimension(:,:,:) :: temp_x, temp_y, temp_z, advect

  allocate(temp_x(nx,dim2(rank),nz),temp_y(nx,dim2(rank),nz),&
   temp_z(nx,dim2(rank),nz),advect(nx,dim2(rank),nz))


  call ddxyz(v_x, dvxdx, v_y, dvydy, v_z, dvzdz,1)
 
  call ddxyz(p, gradp_x, p, gradp_y, p, gradp_z,0)

!  call ddx(p, gradp_x, 1)
!  call ddy(p, gradp_y, 1)
!  call ddz(p, gradp_z, 4)
 
  call PML_BOUNDARY_ACCOUNTING_3D ()

  div = dvxdx + dvydy + dvzdz

  call curl(v_x, v_y, v_z, omega_x, omega_y, omega_z)
  call cross(v0_x, v0_y, v0_z, omega_x, omega_y, omega_z, flow_x, flow_y, flow_z)  
  call cross(v_x, v_y, v_z, omega0_x, omega0_y, omega0_z, temp_x, temp_y, temp_z)  


  ! Using the vector identity that v_0 dot nalbla v' + v' dot nabla v_0 =
  ! nalbla(v_0 dot v') - v_0 X omega' - v' X omega_0 

  flow_x = flow_x + temp_x
  flow_y = flow_y + temp_y
  flow_z = flow_z + temp_z

  advect = v0_x * v_x + v0_y * v_y + v0_z * v_z

  call ddxyz(advect, temp_x, advect, temp_y, advect, temp_z,1)
!  call ddx(advect, temp_x, 1)
!  call ddy(advect, temp_y, 1)
!  call ddz(advect, temp_z, 1)

  flow_x = flow_x - temp_x
  flow_y = flow_y - temp_y
  flow_z = flow_z - temp_z

  deallocate(temp_x, temp_y, temp_z, advect)

  call ddxyz(rho, gradrho_x, rho, gradrho_y, rho, gradrho_z,1)

!  call ddx(rho, gradrho_x, 1)
!  call ddy(rho, gradrho_y, 1)
!  call ddz(rho, gradrho_z, 1)


 ! Gravity is assumed to be directed inwards
  ! \vec{g}  = - g \vec{e_r}

   do k = 1,nz

     RHSv_z(:,:,k)  = - rhoinv(:,:,k) *  rho(:,:,k) * g(k)

     if (ABS(k-e_rad) .LE. 4) then
        RHSv_z(:,:,k) = RHSv_z(:,:,k) +  forcing*source_dep(k)
     endif

   enddo

   RHScont = -rho0*div - v_z * gradrho0_z - v0_x * gradrho_x - v0_y * gradrho_y &
    - v0_z * gradrho_z - rho * (div_v0)

   RHSv_x = -rhoinv * gradp_x + flow_x + rho * advect0_x

   RHSv_y = -rhoinv * gradp_y + flow_y + rho * advect0_y 

   RHSv_z = RHSv_z - rhoinv * gradp_z + flow_z + rho * advect0_z

   RHSp   = - v_z * gradp0_z - rho * c2div_v0 - c2rho0 * div - v0_x * gradp_x - v0_y * gradp_y & 
   & - v0_z * gradp_z  
   

  if (HORIZONTAL_PMLS) then

   scr(1,:,:,:) = 0.0
   scr(nx,:,:,:) = 0.0

   if (ystart(rank) .eq. 1) &
     scr(:,1,:,:) = 0.0
   if (ystart(rank) + dim2(rank) -1 .eq. nx) &
     scr(:,dim2(rank),:,:) = 0.0

   endif

  if (DAMP_WAVES) call DAMP_VELOCITY

 end SUBROUTINE MP_FLOWS_PML_3D

!==================================================================================

END MODULE PML


