MODULE KERNELS

 use initialize
 use all_modules
 implicit none
 integer fwd_start_stamp, adj_end_stamp
 Contains



!==================================================================================

SUBROUTINE PRODUCE_KERNELS

 implicit none
 integer i,j,k, ierr, iter
 real*8 time1, time2 , tot_time, dimen,tempnor, ltime, UNKNOWN_FACTOR,Lregular
 real*8, allocatable, dimension(:,:,:) :: temp1, temp2, temp3, tempx, &
	tempy, tempz, fbx, fby, fbz, abx, aby, abz, for_figures, temp, kernelpsi,c2a
 real*8, allocatable, dimension(:,:,:,:,:) :: deriv_ij
 character* (timestamp_size) charnum
 character*80 directory_rel
 logical lexist

 directory_rel = directory//'kernel/'
 UNKNOWN_FACTOR = 1.!2./pi

 allocate( temp1(nx,dim2(rank),nz_kern), temp2(nx,dim2(rank), nz_kern), &
	temp3(nx,dim2(rank),nz_kern), tempx(nx,dim2(rank),nz_kern), &
	deriv_ij(nx, dim2(rank), nz_kern, 3, 3),kernelpsi(nx,dim2(rank),nz))

  test_in_2d = .true.


  inquire(file=directory//'model_c_ls'//jobno//'.fits',exist = lexist)
  if (lexist) then
   call readfits(directory//'model_c_ls'//jobno//'.fits', temp1, nz)
   if (sum(temp1) .ne. 0) c2 = temp1
   c2 = (c2/dimc)**2.0
   c_speed = c2**0.5
   cq = c2(1,1,st_z:fi_z)**0.5
  endif
  

  if (flows) then
   inquire(file=directory//'model_psi_ls'//jobno//'.fits', exist = lexist)
   if (lexist) then
    allocate(psivar(nx,dim2(rank),nz_kern))
    call readfits(directory//'model_psi_ls'//jobno//'.fits',psivar,nz)
!    deallocate(psivar)
   endif
  endif



  if (magnetic) then

   allocate(fbx(nx,dim2(rank),nz_kern), fby(nx,dim2(rank),nz_kern), &
        fbz(nx,dim2(rank),nz_kern), abx(nx,dim2(rank),nz_kern), &
        aby(nx,dim2(rank),nz_kern), abz(nx,dim2(rank),nz_kern), &
        tempy(nx,dim2(rank),nz_kern), tempz(nx,dim2(rank),nz_kern),&
	temp(nx,dim2(rank),nz))


   !call readfits(dirbackmag//'densityclassic.fits', temp, nz)
!   call readfits(dirbackmag//'density2D.fits', temp, nz)
 !  rho0 = temp(:,:,st_z:fi_z)/dimrho
  
   !call readfits(dirbackmag//'pressureclassic.fits', temp, nz)
!   call readfits(dirbackmag//'pressure2D.fits', temp, nz)
!   p0 = temp(:,:,st_z:fi_z)/(dimrho * dimc**2.)

   !call readfits(dirbackmag//'B_xclassic.fits', temp, nz)
!   call readfits(dirbackmag//'B_x2D.fits', temp, nz)
!   box = temp(:,:,st_z:fi_z)/dimb

!   boy = 0.0
!   if (.not. test_in_2D) then 
!    call readfits(dirbackmag//'B_yclassic.fits', temp, nz)
!    boy = temp(:,:,st_z:fi_z)/dimb
!   endif

   !call readfits(dirbackmag//'B_zclassic.fits', temp, nz)
!   call readfits(dirbackmag//'B_z2D.fits', temp, nz)
!   boz = temp(:,:,st_z:fi_z)/dimb


   !call readfits(dirbackmag//'soundspeedclassic.fits', temp, nz)
!   call readfits(dirbackmag//'soundspeed2D.fits', temp, nz)
 !  c2 = temp(:,:,st_z:fi_z)**2./dimc**2.
  ! cq = c2(1,1,st_z:fi_z)**0.5

!   allocate(reduction(nx,dim2(rank),nz_kern))
!    do k=1,nz_kern
!     reduction(:,:,k) = (boz(:,:,k)**2.+ box(:,:,k)**2.+boy(:,:,k)**2.)/(rho0(:,:,k)*cq(k)**2.)
!    enddo

!    reduction = 8./(8. + reduction)
!    reduction = 1.

!    do k=nz-10, nz
!     reduction(:,:,k) = reduction(:,:,k) * ((k - nz)/10.0)**4.
!    enddo
!    reduction = 1.0
    boy=0.0
  
    inquire(file=directory//'model_psi_ls'//jobno//'.fits', exist = lexist)
    lexist = .true.
    if (lexist) then
     allocate(psivar(nx,dim2(rank),nz_kern))
     call readfits(directory//'model_psi_ls'//jobno//'.fits',psivar,nz)
     call ddzkern(psivar, box, 0)
     box = -box
!     box(:,:,180:230) = 0.0
     print *,"Max value of Bx",maxval(abs(box))
!~      call writefits_3d(directory//'Bx.fits', box, nz)
     !call ddxkern(psivar, boz, 1)
     call dbyd1(boz,psivar,nx,dim2(rank)*nz,0)
     boz = boz*stretchx
     print *,"Max value of Bz",maxval(abs(boz))
     deallocate(psivar)
   

     call readfits(dirbacktrue//'pressure2D.fits',p0,nz)
     call readfits(dirbacktrue//'density2D.fits',rho0,nz)

!     call readfits(dirbackmag//'B_x2D.fits', temp, nz)
!     call readfits(dirbackmag//'B_z2D.fits', temp, nz)

     box = box/dimb
     boz = boz/dimb
     p0 = p0/(dimrho * dimc**2.)
     rho0 = rho0/dimrho



     call curl_kern(box, boy, boz, curlbox, curlboy, curlboz)
     !curlboy = curlboy * reduction**0.5
!    box = box * reduction**0.5
!    if (.not. test_in_2D) boy = boy * reduction**0.5
!    boz = boz * reduction**0.5

!    deallocate(reduction)
    else
     print *,'No psi!'
     stop
    endif

    deallocate(temp)
  endif



  call ddxkern(rho0, gradrho0_x, 1)
  call ddykern(rho0, gradrho0_y, 1)
  call ddzkern(rho0, gradrho0_z, 1)

  if (STABILIZE_MODEL) then
   allocate(c2a(nx,dim2(rank),nz_kern))
!    FOR CSM_B
   do k=1,nz_kern
    do j=1,dim2(rank)
     do i=1,nx

      c2a(i,j,k) = c2(i,j,k) *(1. + 0.1 * exp(-((Rsun-Rsun*zkern(k))*10.0**(-8.) )**2. ) )
      c2(i,j,k) = c2a(i,j,k)* (1. - 0.03 * exp(-((Rsun-Rsun*zkern(k))/(5*10.0**(8.)))**2.))

     enddo
     enddo
    enddo
    deallocate(c2a)
  endif



!  if (background_flows_exist) then
!   call readfits(dirbackmag//'v0_x.fits', v0_x, nz_kern)
!   call readfits(dirbackmag//'v0_y.fits', v0_y, nz_kern)
!   call readfits(dirbackmag//'v0_z.fits', v0_z, nz_kern)
!  endif

 hessian = 0.0
 kernelv = 0.0
 kernelc2 = 0.0
 kernelp = 0.0
 kernelrho = 0.0



 sound_speed_kernels = .true.
 magnetic_kernels = .false.
 flow_kernels = .false.
 if (FLOWS) flow_kernels = .true.
 if (magnetic) magnetic_kernels = .true.
 density_kernels = .true.


! grav =0.
 
 if (magnetic) kernelb = 0.0

 stepskern = floor(outputcad/timestep)

 call read_parameters_kernel

! nt_kern = floor((final_time - local_time)*timestep/out	putcad) + 1
 if (rank==0) print *,'Total number of temporal steps:', nt_kern,totkern

 do k=1,nt_kern

  if (rank==0)  print *,'At time step:', k, ' out of ', nt_kern


  time1 = MPI_WTIME ()
  ltime = (k-1.)*stepskern + fwd_start_stamp
  call convert_to_string(int(abs(ltime)),charnum,timestamp_size)

  !if (ltime .lt. 0) charnum = '-'//tempct
  !if (ltime .ge. 0) charnum = '+'//tempct


  ! FORWARD-TIME FORWARD 

  call readfits(directory//'forward_src'//contrib//'_ls'//jobno//'/xiz_'//charnum//'_partial.fits', &
		f_xi_z, nz_kern)
  f_xi_y = 0.
  if (.not. test_in_2d) call readfits(directory//'forward_src'//contrib//'_ls'//jobno//'/xiy_'//charnum//'_partial.fits', &
		f_xi_y, nz_kern)
  call readfits(directory//'forward_src'//contrib//'_ls'//jobno//'/xix_'//charnum//'_partial.fits', &
		f_xi_x, nz_kern)

  call readfits(directory//'forward_src'//contrib//'_ls'//jobno//'/vz_'//charnum//'_partial.fits', &
		f_vel_z, nz_kern)

  f_vel_y = 0.
  if (.not. test_in_2d)  &
  call readfits(directory//'forward_src'//contrib//'_ls'//jobno//'/vy_'//charnum//'_partial.fits', &
		f_vel_y, nz_kern)

  call readfits(directory//'forward_src'//contrib//'_ls'//jobno//'/vx_'//charnum//'_partial.fits', &
		f_vel_x, nz_kern)



  call readfits(directory//'forward_src'//contrib//'_ls'//jobno//'/acc_z_'//charnum//'_partial.fits', &
		f_acc_z, nz_kern)
  f_acc_y = 0.
  if (.not. test_in_2d) call readfits(directory//'forward_src'//contrib//'_ls'//jobno//'/acc_y_'//charnum//'_partial.fits', &
		f_acc_y, nz_kern)
  call readfits(directory//'forward_src'//contrib//'_ls'//jobno//'/acc_x_'//charnum//'_partial.fits', & 
	f_acc_x, nz_kern)



  ! REVERSED-TIME ADJOINT 
  ltime = adj_end_stamp - (k-1)*stepskern !*timestep + timestep
  call convert_to_string(int(abs(ltime)),charnum,timestamp_size)


!  call readfits(adjoint_directory//'xiz_'//charnum//'_partial.fits', a_xi_z, nz_kern) 
!  call readfits(adjoint_directory//'xiy_'//charnum//'_partial.fits', a_xi_y, nz_kern)
!  call readfits(adjoint_directory//'xix_'//charnum//'_partial.fits', a_xi_x, nz_kern)

  call readfits(directory//'adjoint_src'//contrib//'/vz_'//charnum//'_partial.fits', & !'_'//contrib//
	a_vel_z, nz_kern) 
  a_vel_y = 0.
  if (.not. test_in_2d) call readfits(directory//'adjoint_src'//contrib//'/vy_'//charnum//'_partial.fits', & !'_'//contrib//
	a_vel_y, nz_kern)
  call readfits(directory//'adjoint_src'//contrib//'/vx_'//charnum//'_partial.fits', &!//'_'//contrib
	a_vel_x, nz_kern)



  call readfits(directory//'adjoint_src'//contrib//'/acc_z_'//charnum//'_partial.fits', &
		a_acc_z, nz_kern)
  a_acc_y = 0.
  if (.not. test_in_2d) call readfits(directory//'adjoint_src'//contrib//'/acc_y_'//charnum//'_partial.fits', &
		a_acc_y, nz_kern)
  call readfits(directory//'adjoint_src'//contrib//'/acc_x_'//charnum//'_partial.fits', & 
	a_acc_x, nz_kern)


  time2 = MPI_WTIME ()

  if (rank==0) print *,'Read in time:', (time2 - time1) 

  time1 = MPI_WTIME ()


!  hessian = hessian + f_acc_z * a_acc_z + f_acc_x * a_acc_x + &
!	f_acc_y * a_acc_y

  hessian = hessian + f_acc_z * f_acc_z + f_acc_x * f_acc_x + &
	f_acc_y * f_acc_y

 
  !print *,f_acc_x 
  !print *,f_acc_y
  !print *,f_acc_z

  if (sound_speed_kernels) then

   call compute_div(f_xi_x,f_xi_y,f_xi_z,temp2)
   call compute_div(a_vel_x,a_vel_y,a_vel_z,temp3)
   temp1 = temp2 * temp3 
   kernelc2 = kernelc2 + temp1

  endif

  tempnor = normkern(kernelc2)
  if (rank==0) print*,'Norm is:', tempnor,maxval(abs(a_vel_z))


  if(density_kernels) then

    kernelrho = kernelrho + temp1 * c2 + f_acc_z * a_vel_z  + &
		f_acc_x * a_vel_x + &
		f_acc_y * a_vel_y
 


   do j=1,nz_kern
    temp2(:,:,j) = a_vel_z(:,:,j) * grav(j)
   enddo
    call ddxkern(temp2, deriv_ij(:,:,:,3,1), 1)
    call ddykern(temp2, deriv_ij(:,:,:,3,2), 1)
    call ddzkern(temp2, deriv_ij(:,:,:,3,3), 1)
 ! DONT KNOW ABOUT THE NEGATIVE/POSITIVE SIGN FOR THE PRODUCT AT THE END
   do j=1,nz_kern
    kernelrho(:,:,j) = kernelrho(:,:,j) - grav(j) * f_xi_z(:,:,j) * temp3(:,:,j) &
		+ (f_xi_x(:,:,j) * deriv_ij(:,:,j,3,1) + f_xi_y(:,:,j) * deriv_ij(:,:,j,3,2) + &
		f_xi_z(:,:,j) * deriv_ij(:,:,j,3,3)) 

   enddo

  endif

  if (FLOW_KERNELS) THEN

   call compute_del(f_vel_x,f_vel_y,f_vel_z,deriv_ij)

   kernelv(:,:,:,1) = deriv_ij(:,:,:,1,1) * a_vel_x + &
		     deriv_ij(:,:,:,1,2) * a_vel_y + &
		     deriv_ij(:,:,:,1,3) * a_vel_z + &
		     kernelv(:,:,:,1)

   kernelv(:,:,:,2) = deriv_ij(:,:,:,2,1) * a_vel_x + &
		     deriv_ij(:,:,:,2,2) * a_vel_y + &
		     deriv_ij(:,:,:,2,3) * a_vel_z + &
		     kernelv(:,:,:,2)

   kernelv(:,:,:,3) = deriv_ij(:,:,:,3,1) * a_vel_x + &
		     deriv_ij(:,:,:,3,2) * a_vel_y + &
		     deriv_ij(:,:,:,3,3) * a_vel_z + &
		     kernelv(:,:,:,3)


  endif


  if (MAGNETIC_KERNELS) then

   call cross_kern(f_xi_x,f_xi_y,f_xi_z, &
		box, boy, boz, tempx, tempy, tempz)
 
   call curl_kern(tempx, tempy, tempz, fbx, fby, fbz)
 
   call curl_kern(fbx, fby, fbz, tempx, tempy, tempz)

   call cross_kern(a_vel_x,a_vel_y,a_vel_z, & 
	tempx, tempy, tempz, temp1, temp2, temp3)

   kernelb(:,:,:,1) = kernelb(:,:,:,1) + temp1



   kernelb(:,:,:,2) = kernelb(:,:,:,2) + temp2
   kernelb(:,:,:,3) = kernelb(:,:,:,3) + temp3

  ! THE FIRST TERM

   call cross_kern(fbx, fby, fbz, a_vel_x,a_vel_y,a_vel_z, &
	tempx, tempy, tempz)

   call curl_kern(tempx, tempy, tempz, temp1, temp2, temp3)


   kernelb(:,:,:,1) = kernelb(:,:,:,1) + temp1
   kernelb(:,:,:,2) = kernelb(:,:,:,2) + temp2
   kernelb(:,:,:,3) = kernelb(:,:,:,3) + temp3


  ! THE THIRD TERM 



   call cross_kern(a_vel_x,a_vel_y,a_vel_z, &
		box, boy, boz, tempx, tempy, tempz)

   call curl_kern(tempx, tempy, tempz, abx, aby, abz)

   call curl_kern(abx, aby, abz, tempx, tempy, tempz)

   call cross_kern(f_xi_x, f_xi_y, f_xi_z, & 
	tempx, tempy, tempz, temp1, temp2, temp3)


   kernelb(:,:,:,1) = kernelb(:,:,:,1) + temp1
   kernelb(:,:,:,2) = kernelb(:,:,:,2) + temp2
   kernelb(:,:,:,3) = kernelb(:,:,:,3) + temp3

 ! THE SECOND TERM (SYMMETRIC)

   call cross_kern(a_vel_x,a_vel_y,a_vel_z, &
		curlbox, curlboy, curlboz, temp1, temp2, temp3)



   call curl_kern(temp1, temp2, temp3, tempx, tempy, tempz)

   call cross_kern(tempx, tempy, tempz, f_xi_x, f_xi_y, &
 	f_xi_z, temp1, temp2, temp3)
 
 
   kernelb(:,:,:,1) = kernelb(:,:,:,1) + temp1
   kernelb(:,:,:,2) = kernelb(:,:,:,2) + temp2
   kernelb(:,:,:,3) = kernelb(:,:,:,3) + temp3
 


  ! THE FOURTH TERM  (FINAL TERM)

   call compute_div(a_vel_x,a_vel_y,a_vel_z,temp3)
   

   call cross_kern(f_xi_x, f_xi_y, f_xi_z, curlbox, curlboy, curlboz,  &
	tempx, tempy, tempz) 

   kernelb(:,:,:,1) = kernelb(:,:,:,1) + tempx * temp3
   kernelb(:,:,:,2) = kernelb(:,:,:,2) + tempy * temp3
   kernelb(:,:,:,3) = kernelb(:,:,:,3) + tempz * temp3

   ! EQUILIBRIUM TERM 2

   tempx = f_xi_x * temp3
   tempy = f_xi_y * temp3
   tempz = f_xi_z * temp3

   call cross_kern(box,boy, boz, tempx, tempy, tempz, &
	temp1, temp2, temp3) 

   call curl_kern(temp1, temp2, temp3, tempx, tempy, tempz)

   kernelb(:,:,:,1) = kernelb(:,:,:,1) + tempx
   kernelb(:,:,:,2) = kernelb(:,:,:,2) + tempy
   kernelb(:,:,:,3) = kernelb(:,:,:,3) + tempz



  ! EQUILIBRIUM CONTRIBUTION 1


  endif


!  if (mod(k,49) == 0 ) then
   
  ! call convert_to_string(k,charnum,timestamp_size)

 !  allocate(for_figures(nx, dim2(rank), nz))

!   for_figures = deriv_ij(:,:,:,1,1) * a_vel_x + &
!		     deriv_ij(:,:,:,1,2) * a_vel_y + &
!		     deriv_ij(:,:,:,1,3) * a_vel_z 
!
!   call writefits_3d(KERNEL_DIRECTORY//'interaction_x_'//charnum//'.fits', for_figures, nz_kern)
!   call writefits_3d(KERNEL_DIRECTORY//'adjoint_z_'//charnum//'.fits', a_vel_x, nz_kern)
!   call writefits_3d(KERNEL_DIRECTORY//'forward_z_'//charnum//'.fits', f_vel_x, nz_kern)

!   call writefits_3d(KERNEL_DIRECTORY//'kernel_v0_x_'//charnum//'.fits', kernelv(:,:,:,1), nz_kern)
!   call writefits_3d(KERNEL_DIRECTORY//'kernel_v0_y_'//charnum//'.fits', kernelv(:,:,:,2), nz_kern)
!   call writefits_3d(KERNEL_DIRECTORY//'kernel_v0_z_'//charnum//'.fits', kernelv(:,:,:,3), nz_kern)

!   deallocate(for_figures)
!  endif

  time2 = MPI_WTIME ()
  if (rank==0) print *,'Compute time:', (time2 - time1) 


 enddo


 print *,"Writing out hessian"
 call writefits_3d(adjustl(trim(directory_rel))//'hessian_'//contrib//'.fits', hessian, nz_kern)

 print *,"Flow kernels"
 if (FLOW_KERNELS) then

  ! Mm^-3 km^-1 
  !dimen = (diml * 10.0**(-8.))**(-4.) * 0.001  !* (dimc * 10.0**(-5.))**(-1.) 
  dimen = (diml * 10.0**(-8.))**(-3.) * (dimc * 10.0**(-5.))**(-1.)  
  if (test_in_2d) dimen = (diml * 10.0**(-8.) * xlength * 10.0**(-8) * dimc * 10.0**(-5.))**(-1.) * UNKNOWN_FACTOR
           ! Mm^(-3.) bit !                 ! s/km bit !  The other delta-t^2. already incorporated  

  Lregular = 30.*10.**8/diml
  temp2= 0.0


  if (enf_cont .and. psi_cont) then
  call curl_kern(kernelv(:,:,:,1), temp2, kernelv(:,:,:,3), tempx, kernelpsi, temp1)
  kernelpsi = 2.0*kernelpsi * dimen * UNKNOWN_FACTOR * &
                   (rho0*c2**0.5*psivar)*stepskern*timestep * dimen
!~  without a leading minus sign, the kernel is gradient in param space.
!~  with a leading minus sign, it is the negative of the gradient.
  
  end if

  print *,"Velocity kernels"
  do i=1,3
   kernelv(:,:,:,i) = - 2.0 * kernelv(:,:,:,i) * rho0 * stepskern * timestep * dimen ! the other time-unit
  enddo



  if (rank==0) print *,'DIMENSIONS OF VELOCITY KERNELS ARE Mm^-3 km^-1 s^3'
  call writefits_3d(adjustl(trim(directory_rel))//'kernel_vx_'//contrib//'.fits', kernelv(:,:,:,1), nz_kern)
  if (.not. TEST_IN_2D) &
     call writefits_3d(adjustl(trim(directory_rel))//'kernel_vy_'//contrib//'.fits', kernelv(:,:,:,2), nz_kern)
  call writefits_3d(adjustl(trim(directory_rel))//'kernel_vz_'//contrib//'.fits', kernelv(:,:,:,3), nz_kern)
  if (enf_cont .and. psi_cont) then
  call writefits_3d(adjustl(trim(directory_rel))//'kernel_psi_'//contrib//'.fits', kernelpsi, nz_kern)
  endif

 endif

 print *,"Magnetic kernels"
 if (MAGNETIC_KERNELS) then

  dimen =  (diml * 10.0**(-8.))**(-3.) * dimb**(-1.) * 1000.

  if (test_in_2d) dimen = (diml * 10.0**(-8.) * xlength * 10.0**(-8.) * dimb)**(-1.) * UNKNOWN_FACTOR !* 1000.

  do i=1,3
   kernelb(:,:,:,i) = kernelb(:,:,:,i) * stepskern * timestep * dimen 
  enddo

  temp2 = 0.0
  call curl_kern(kernelb(:,:,:,1), temp2, kernelb(:,:,:,3), tempx, kernelpsi, temp1)
!  kernelpsi = kernelpsi/(diml*10.0**(-8.))

  if (rank==0) print*,'DIMENSIONS OF B-FIELD KERNELS ARE Mm^(-3) * kG^(-1.) * s^2'
  call writefits_3d(adjustl(trim(directory_rel))//'kernel_bx_'//contrib//'.fits', kernelb(:,:,:,1), nz_kern)
  if (.not. test_in_2d) &
    call writefits_3d(adjustl(trim(directory_rel))//'kernel_by_'//contrib//'.fits', kernelb(:,:,:,2), nz_kern)
  call writefits_3d(adjustl(trim(directory_rel))//'kernel_bz_'//contrib//'.fits', kernelb(:,:,:,3), nz_kern)
  call writefits_3d(adjustl(trim(directory_rel))//'kernel_psi_'//contrib//'.fits', kernelpsi, nz_kern)


 endif
 
! if (rank==0) print *,'DIMENSIONS OF NORMALIZED PRESSURE KERNELS ARE Mm^-3 s^2'
 dimen = (diml * 10.0**(-8.))**(-3.)
! kernelp = - kernelp * p0 * stepskern * timestep * dimen
!
 if (test_in_2d) dimen = (diml * 10.0**(-8.) * xlength * 10.0**(-8.))**(-1.) *UNKNOWN_FACTOR

 print *,"Density kernels"
 if (DENSITY_KERNELS) then
  if (rank==0) print *,'DIMENSIONS OF NORMALIZED DENSITY KERNELS ARE Mm^-3 s^2'
  kernelrho = - rho0 * kernelrho * stepskern * timestep * dimen
  if (flows .and. psi_cont) kernelrho = kernelrho + kernelpsi*(1.-psivar(1,1,1)/psivar)
  call writefits_3d(adjustl(trim(directory_rel))//'kernel_rho_'//contrib//'.fits', kernelrho, nz_kern)
 endif

 print *,"Sound speed kernels"
 if (SOUND_SPEED_KERNELS) then
  if (rank==0) print *,'DIMENSIONS OF NORMALIZED SOUND-SPEED KERNELS ARE s^2 Mm^-3'
  kernelc2 = - kernelc2 * rho0 * c2 * stepskern * timestep * dimen * 2.0
  if (flows .and. psi_cont) kernelc2 = kernelc2 + kernelpsi*(1.-psivar(1,1,1)/psivar)*2.0
  call writefits_3d(adjustl(trim(directory_rel))//'kernel_c_'//contrib//'.fits', kernelc2, nz_kern)

!~   if (rank==0) then
!~   open(443,file=directory//'status/kernel_src'//contrib,status='unknown')
!~   close(443)
!~   endif

 endif

 
 if (rank ==0) then
  print *,'Removing forward and adjoint data'
  call system('rm -rf '//directory//'forward_src'//contrib//'_ls'//jobno//'/*partial*')
!  call system('rm -rf '//directory//'forward_src'//contrib//'_ls'//jobno//'/*full*')
  call system('rm -rf '//directory//'adjoint_src'//contrib//'/*partial*')
 ! call system('rm -rf '//directory//'adjoint'//contrib//'/*full*')
 endif

 deallocate(temp1, temp2, temp3, tempx, deriv_ij,kernelpsi)

 if (magnetic) deallocate(tempy, tempz, fbx, fby, fbz, abx, aby, abz)


 if (rank==0) print *,''
 if (rank==0) print *,'FINISHED COMPUTING KERNELS'
 if (rank==0) print *,''



!  call writefits_3d(KERNEL_DIRECTORY//'kernel_p0.fits', kernelp, nz_kern)


  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  call MPI_FINALIZE(ierr)

END SUBROUTINE PRODUCE_KERNELS 

!================================================================================

SUBROUTINE COMPUTE_FREQS(nt, deltatime, nus)
 implicit none
 integer nt,i
 real*8 deltatime, nus(nt)

 if ((nt/2) *2 .eq. nt) then

  do i=1,nt/2+1
   nus(i) = 1./(nt*deltatime) * (i - 1.0)
  enddo

  do i=2,nt/2
   nus(nt-i +2) = -nus(i)
  enddo

 else

  do i=1,nt/2+1
   nus(i) = 1./(nt*deltatime) * (i - 1.0)
  enddo

  do i=2,nt/2+1
   nus(nt-i +2) = - nus(i)
  enddo

 endif

END SUBROUTINE COMPUTE_FREQS

!================================================================================


SUBROUTINE ddzkern(var, dvar, bc)
 implicit none
 real*8 var(nx,dim2(rank),nz_kern), dvar(nx,dim2(rank),nz_kern)
 integer i,j,k,bc


  call dbyd2(dvar,var,nx*dim2(rank),nz_kern,bc)
  do k=1,nz_kern
   dvar(:,:,k) = dvar(:,:,k)*stretchkern(k)
  enddo
  
  call filter_var_z(dvar)
 
END SUBROUTINE ddzkern

!================================================================================================ 

SUBROUTINE ddxkern(var, dvar,bc)
 implicit none
 real*8 var(nx,dim2(rank),nz_kern), dvar(nx,dim2(rank),nz_kern)
 integer i,j,k,bc
 complex*16, allocatable, dimension(:,:,:) :: temp

 if ((PERIODIC) .AND. (USE_FFT)) then
  allocate(temp(nx/2+1,dim2(rank), nz_kern))
  call dfftw_execute_dft_r2c(fftw_plan_fwd_x, var, temp)

  do i=1,nx/2+1
   temp(i,:,:) = temp(i,:,:)*eyekx(i)
  enddo

  call dfftw_execute_dft_c2r(fftw_plan_inv_x, temp, dvar)
  deallocate(temp)

 elseif (PERIODIC .AND. (.NOT. USE_FFT)) then
  call dbyd1(dvar,var,nx,dim2(rank)*nz_kern,5)

 elseif (.NOT. PERIODIC) then
  call dbyd1(dvar,var,nx,dim2(rank)*nz_kern,bc)

 endif
 dvar = dvar*stretchx

END SUBROUTINE ddxkern

!================================================================================================ 


SUBROUTINE ddykern(var, dvar,bc)
 implicit none
 integer i,j,k,ierr,sendtag,recvtag,stat(MPI_STATUS_SIZE,2)
 integer n, statt(MPI_STATUS_SIZE), req(2), bc, reqq, tag, nelements
 real*8 var(nx,dim2(rank),nz_kern), dvar(nx,dim2(rank),nz_kern)
 real*8, dimension(:,:,:), allocatable :: trans, temp
 real*8, dimension(:,:), allocatable :: rectemp
 complex*16, allocatable, dimension(:,:,:) :: tempcomp

 allocate(temp(nx,dim1(rank),nz_kern) , trans(ny, dim1(rank), nz_kern))
 call transpose_3D_y_kern(var, temp)
 


 if ((PERIODIC) .AND. (USE_FFT)) then
  allocate(tempcomp(ny/2+1,dim1(rank), nz_kern))
  call dfftw_execute_dft_r2c(fftw_plan_fwd_y, temp, tempcomp)

  do i=1,ny/2+1
   tempcomp(i,:,:) = tempcomp(i,:,:)*eyeky(i)
  enddo

  call dfftw_execute_dft_c2r(fftw_plan_inv_y, tempcomp, trans)
  deallocate(tempcomp)
 
 elseif (PERIODIC .AND. (.NOT. USE_FFT)) then
  call dbyd1(trans,temp,ny,dim1(rank)*nz_kern,5)

 elseif (.NOT. PERIODIC) then
  call dbyd1(trans,temp,ny,dim1(rank)*nz_kern,bc)
 endif


 call inv_transpose_3D_y_kern(trans, dvar)
   
 deallocate(trans, temp)

 dvar = dvar*stretchy

END SUBROUTINE ddykern

!================================================================================================ 


SUBROUTINE TRANSPOSE_3D_Y_KERN(input, output)

 implicit none
 integer i, j, k, sendtag, recvtag, req(2), ierr, stat(MPI_STATUS_SIZE, 2)
 real*8 input(nx, dim2(rank), nz_kern), output(ny, dim1(rank), nz_kern)

 ! Non - communicative transpose (local transpose)
 do j=1,dim2(rank)
   do i=1,dim1(rank)
    output(j+fwd_rdispls(rank)-1,i,:) = input(i+fwd_sdispls(rank)-1,j,:)
   enddo
 enddo

 ! Non-blocking send-recv transpose
 do k=1,numtasks-1
   i = modulo(rank+k,numtasks)
   j = modulo(rank-k,numtasks)
   sendtag = 0
   recvtag = 0
   call MPI_ISEND(input(fwd_sdispls(i),1,1),1,fwd_sendtypes(i),i,sendtag,MPI_COMM_WORLD,req(1),ierr)
   call MPI_IRECV(output(fwd_rdispls(j),1,1),1,fwd_recvtypes(j),j,recvtag,MPI_COMM_WORLD,req(2),ierr)
   call MPI_WAITALL(2,req,stat,ierr)
 enddo

END SUBROUTINE TRANSPOSE_3D_Y_KERN


!================================================================================================ 


SUBROUTINE INV_TRANSPOSE_3D_Y_KERN(input, output)

 implicit none
 integer i, j, k, sendtag, recvtag, req(2), ierr, stat(MPI_STATUS_SIZE, 2)
 real*8 input(ny, dim1(rank), nz_kern), output(nx, dim2(rank), nz_kern)

 ! Non - communicative transpose


 do j=1,dim1(rank)
  do i=1,dim2(rank)
   output(j+inv_rdispls(rank)-1,i,:) = input(i+inv_sdispls(rank)-1,j,:)
  enddo
 enddo

  ! Non-blocking send-recv transpose

 do k=1,numtasks-1
   i = modulo(rank+k,numtasks)
   j = modulo(rank-k,numtasks)
   sendtag = 0
   recvtag = 0
   call MPI_ISEND(input(inv_sdispls(i),1,1),1,inv_sendtypes(i),i,sendtag,MPI_COMM_WORLD,req(1),ierr)
   call MPI_IRECV(output(inv_rdispls(j),1,1),1,inv_recvtypes(j),j,recvtag,MPI_COMM_WORLD,req(2),ierr)
   call MPI_WAITALL(2,req,stat,ierr)
 enddo

END SUBROUTINE INV_TRANSPOSE_3D_Y_KERN

!================================================================================================
                                                                                                                                                            
SUBROUTINE CROSS_KERN(a_x, a_y, a_z, b_x, b_y, b_z, c_x, c_y, c_z)
                                                                                                                                                            
  implicit none
  real*8, dimension(nx, dim2(rank), nz_kern) :: a_x, a_y, a_z, b_x, b_y, b_z, c_x, c_y, c_z


  c_x = a_y * b_z - a_z * b_y
  c_y = b_x * a_z - a_x * b_z
  c_z = a_x * b_y - a_y * b_x

END SUBROUTINE CROSS_KERN

!==================================================================================

SUBROUTINE CURL_KERN(f_x, f_y, f_z, curl_x, curl_y, curl_z)
                                                                                                                                                            
  implicit none
  real*8, dimension(nx, dim2(rank), nz_kern) :: curl_x, curl_y, curl_z
  real*8, dimension(nx, dim2(rank), nz_kern) :: f_x, f_y, f_z
  real*8, dimension(:,:,:), allocatable :: temp

  allocate(temp(nx, dim2(rank), nz_kern))
                                                                                                                            
  call ddxkern(f_y, curl_z, bcx)
  call ddykern(f_x, temp, bcy)
  curl_z = curl_z - temp
                                                                                                                                                            
  call ddykern(f_z, curl_x, bcy)
  call ddzkern(f_y, temp, bcz)
  curl_x = curl_x - temp
                                                                                                                                                            
  call ddzkern(f_x, curl_y,bcz)
  call ddxkern(f_z, temp, bcx)
  curl_y = curl_y - temp
                                                                                                                                                     
  deallocate(temp)
                                                                                                                                                            
END SUBROUTINE CURL_KERN


!================================================================================================

SUBROUTINE COMPUTE_DIV(fx, fy, fz, dive)
 implicit none
 real*8, dimension(nx,dim2(rank),nz_kern) :: fx, fy, fz, dive
 real*8, allocatable, dimension(:,:,:) :: temp
 allocate(temp(nx, dim2(rank),nz_kern))
 
 call ddxkern(fx,dive,1)


 if (ny .gt. 4) then
  call ddykern(fy,temp,1)
  dive = dive + temp
 endif

 call ddzkern(fz,temp,1)
 dive = dive + temp

 deallocate(temp)
END SUBROUTINE COMPUTE_DIV

!================================================================================================


SUBROUTINE COMPUTE_DEL(fx, fy, fz, deriv_ij)

 implicit none
 real*8, dimension(nx,dim2(rank),nz_kern) :: fx, fy, fz
 real*8, dimension(nx,dim2(rank),nz_kern,3,3) :: deriv_ij 


 call ddxkern(fx,deriv_ij(:,:,:,1,1),1)
 call ddxkern(fy,deriv_ij(:,:,:,1,2),1)
 call ddxkern(fz,deriv_ij(:,:,:,1,3),1)

 if (ny .gt. 4) then
  call ddykern(fx,deriv_ij(:,:,:,2,1),1)
  call ddykern(fy,deriv_ij(:,:,:,2,2),1)
  call ddykern(fz,deriv_ij(:,:,:,2,3),1)
 endif

 call ddzkern(fx,deriv_ij(:,:,:,3,1),1)
 call ddzkern(fy,deriv_ij(:,:,:,3,2),1)
 call ddzkern(fz,deriv_ij(:,:,:,3,3),1)

END SUBROUTINE COMPUTE_DEL

!================================================================================================

SUBROUTINE FILTER_VAR_Z(var)
  implicit none
  integer i,k,l,kst
  real*8 coeffs(0:5)
  real*8, dimension(nx,dim2(rank),nz_kern) :: temp, var 
  
  coeffs(0) = 0.75390625000000
  coeffs(1) = 0.41015625000000
  coeffs(2) = -0.23437500000000
  coeffs(3) = 0.08789062500000
  coeffs(4) = -0.01953125000000
  coeffs(5) = 0.00195312500000


  temp(:,:,2) = 0.5 * (var(:,:,3) + var(:,:,1))
  temp(:,:,3) = 0.5 * (var(:,:,4) + var(:,:,2))
  temp(:,:,4) = 0.5 * (var(:,:,5) + var(:,:,3))
  temp(:,:,5) = 0.5 * (var(:,:,6) + var(:,:,4))


  do k=6,nz_kern-5
   temp(:,:,k) = coeffs(0)*var(:,:,k) + 0.5*coeffs(1)*(var(:,:,k-1) + var(:,:,k+1)) &
                +  0.5*coeffs(2)*(var(:,:,k-2) + var(:,:,k+2))  &
                +  0.5*coeffs(3)*(var(:,:,k-3) + var(:,:,k+3))  &
                +  0.5*coeffs(4)*(var(:,:,k-4) + var(:,:,k+4))  &
                +  0.5*coeffs(5)*(var(:,:,k-5) + var(:,:,k+5))
  enddo

  temp(:,:,nz_kern-4) = 0.5 * (var(:,:,nz_kern-3) + var(:,:,nz_kern-5))
  temp(:,:,nz_kern-3) = 0.5 * (var(:,:,nz_kern-2) + var(:,:,nz_kern-4))
  temp(:,:,nz_kern-2) = 0.5 * (var(:,:,nz_kern-1) + var(:,:,nz_kern-3))
  temp(:,:,nz_kern-1) = 0.5 * (var(:,:,nz_kern) + var(:,:,nz_kern-2))

  do k=2,nz_kern-1
   var(:,:,k) = temp(:,:,k)
  enddo 

END SUBROUTINE FILTER_VAR_Z

!================================================================================================

function normkern(matrix)

   implicit none
   integer i, j, k, ierr
   real*8 matrix(nx,dim2(rank),nz_kern)
   real*8 normkern, sum
  
   normkern  = 0.0
   sum = 0.0


   do k = 1,nz_kern
    do j =1,dim2(rank)
     do i =1,nx
       sum = sum + matrix(i,j,k)**2.0
     end do     
    end do     
   end do 
 
   call MPI_REDUCE(sum, normkern, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr) 

   normkern = (normkern/(DBLE(nx)*DBLE(ny)*DBLE(nz)))**0.5d0

   call MPI_BARRIER(MPI_COMM_WORLD, ierr)
end function normkern

!================================================================================

SUBROUTINE READ_PARAMETERS_KERNEL

 implicit none
 integer j
 character*80 calculation_type, directory_rel

 open(44,file=directory//'adjoint_src'//contrib//'/kernel_info',form='formatted',status='unknown')



  read(44,*) nt_kern
  read(44,*) adj_end_stamp

 close(44)


 open(44,file=directory//'forward_src'//contrib//'_ls'//jobno//'/kernel_info',form='formatted',status='unknown')

  read(44,*) fwd_start_stamp
  read(44,*) j



 close(44)

! fwd_start_stamp = 10800 - adj_end_stamp


END SUBROUTINE READ_PARAMETERS_KERNEL


!================================================================================



END MODULE KERNELS
