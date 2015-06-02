MODULE INITIALIZE

   implicit none
! -------READING IN BASIC PARAMETER FILES-------!
! --params.i : contains the simulation parameters
! --fftw3.f : The fftw parameter file
! --mpif.h : The MPI parameter file

   INCLUDE 'params.i'
   INCLUDE 'mpif.h'
   INCLUDE 'fftw3.f'
!   INCLUDE '/opt/users/apps/intel/composer_xe_2015.2.164/mkl/include/fftw/fftw3.f'
 !  INCLUDE '/opt/users/apps/intel/impi/5.0.3.048/intel64/include/mpif.h'

!   INCLUDE '/usr/lpp/ppe.poe/include/thread64/mpif.h'



!  INCLUDE THE LOCATION OF THE mpif.h FILE 
!  INCLUDE '/usr/lpp/ppe.poe/include/thread64/mpif.h'
   
! ----------------------------

   integer time,step_rk,dimfive,option, indexglob
   integer MAXTIME,e_rad,num_steps, steps, o_rad, freq_filtering

   integer*8 fftw_plan_fwd_x,fftw_plan_inv_x
   integer*8 fftw_plan_fwd_y,fftw_plan_inv_y
   integer*8 fftw_plan_fwd_2d,fftw_plan_inv_2d


   real*8  pi,deltat,rsun,z(nz), x(nx), y(ny), timeline
   real*8  dimc, diml, spongex(nx), dx, dy, rhoq(nz), &
	   dimrho, normx, normy,nu, spongey(ny), spongez(nz), &
	   height(nz), visc(nz), cq(nz)

   parameter(normx=1.0/DBLE(nx), normy = 1.0/DBLE(ny))
   parameter(rsun = 69598946770.0, pi = 3.14159265358979) 
   real*8  stretch(nz), unstretch(nz), stretchx, stretchy, decay_coeff_x(nx/2+1)
   real*8 decay_coeff_y(ny/2+1), delta_width
   logical initcond, generate_wavefield, linesearch

   ! MPI Related definitions  

   integer  numtasks, rank, status(MPI_STATUS_SIZE)

   ! Domain distribution 
   !
   integer, dimension(:), allocatable :: dim2, ystart, fwd_recvtypes, fwd_sendtypes
   integer, dimension(:), allocatable :: dim1, inv_recvtypes, inv_sendtypes, dimz
   integer, dimension(:), allocatable :: fwd_sdispls, inv_sdispls, fwd_z_displs
   integer, dimension(:), allocatable :: fwd_rdispls, inv_rdispls
   integer, dimension(:), allocatable :: fwd_z_send, inv_z_send, fwd_z_recv, inv_z_recv

   ! Definitions of rank and number of processes
!   integer  datatype, transp
 !  external datatype
  ! external transp

   ! END OF MPI DEFINITIONS
  
   ! Common to all calculations
   ! div - divergence of the velocity vector
   ! vr - forcing applied on the RHS of the radial momentum equation
   ! gamma - the adiabatic index of the background solar state
   ! p0, rho0, T0, c_speed, and g  - pressure, density, temperature, sound speed
   ! and gravity of the background solar state
   ! c2 = c_speed^2, rhoinv = 1/rho0
   ! gradp, gradvr - radial velocity and pressure gradients in 3-space 
   ! omega - vorticities vector
   ! a - array that contains the 5 variables in the calculation: density, radial
   ! velocity, latitudinal velocity, longitudinal velocity and pressure in that order
   ! temp_step, scr - scratch space arrays used in the time evolution algorithm
  
   real*8, allocatable, dimension(:,:,:) :: div, vr, c2, omega_r,spongexyz
   real*8, allocatable, dimension(:,:) :: forcing
   real*8, allocatable, dimension(:,:,:) :: p0, gradp0_z, c_speed, rho0, gradrho0_z, rhoinv, c2rho0, c2div_v0 ,reduction
   real*8, allocatable, dimension(:) :: source_dep, g, gamma
   real*8, allocatable, target, dimension(:,:,:,:) ::  a,temp_step,scr, gradp

   real*8, pointer, dimension(:,:,:) :: rho,RHSv_z,dvzdz,    &
	&	RHSv_y,RHScont,p,RHSp,gradp_z,gradp_y,gradp_x,RHSv_x

   real*8, allocatable, target, dimension(:,:,:) :: dvxdx, dvydy, gradvz

   ! For periodic horizontal boundaries
   complex*16, allocatable, dimension(:) :: eyekx, eyeky
   complex*16, parameter :: eye = (0.0d0,1.0d0)

   !-----Quadratic and Lagrange interpolation stuff------!
   integer time1, time_old
   real*8 x0, x1, x2, x3, x4, x5, x6
   real*8, dimension(:,:), allocatable :: z_i, z_iplus1, LC0, LC1, LC2, LC3, LC4, LC5, LC6

   !-----Flow stuf------! 
   real*8, allocatable, target, dimension(:,:,:,:) :: omega,v0,gradrho,omega0,advect0
   real*8, allocatable, dimension(:,:,:) :: flow_x, flow_y, flow_z, div_v0, psivar

   real*8, pointer, dimension(:,:,:) ::	v0_x, v0_y, v0_z, omega0_x, omega0_y, omega0_z, gradrho_x, &
		gradrho_y, gradrho_z, advect0_x, advect0_y, advect0_z, &
		omega_x, omega_y, omega_z

 
   !-----Magnetic field stuff-----!
  
   real*8 dimb!, reduction(nz)
   real*8, allocatable, dimension(:,:,:) :: box, boy, boz, gradp0_x, gradp0_y, dzbx, dzflux2, &
		curlbox, curlboy, curlboz, curlbx, curlby, curlbz, gradrho0_x, gradrho0_y,&
		dzby, dzflux1, flux1, flux2, flux3
   real*8, pointer, dimension(:,:,:) :: bx, by, bz, RHSb_x, RHSb_y, RHSb_z, v_x, v_y, v_z

   integer, dimension(:,:), allocatable :: erad_2d, orad_2d
   real*8, dimension(:,:), allocatable :: delta_width_2d


   !-----DISPLACEMENT--------!
   real*8, allocatable, target, dimension(:,:,:,:) :: scratch
   real*8, pointer, dimension(:,:,:) :: dxixdx, dxiydy, dxizdz, RHSxi_x, RHSxi_y, RHSxi_z, &
					xi_x, xi_y, xi_z

   !------PML STUFF------!
   real*8, allocatable, dimension(:,:,:) :: az, bzpml
   real*8, target, allocatable, dimension(:,:,:,:) :: scrpml, psipml, pmlvars
   real*8, pointer, dimension(:,:,:) :: psivz, psip, psidzbx, psidzby, RHSpsiinductionbx,&
					psiinductionbx, psiinductionby, RHSpsivz, RHSpsip, &
					RHSpsidzbx, RHSpsidzby, RHSpsiinductionby


   !------HORIZONTAL PML STUFF-----!
   real*8, allocatable, dimension(:,:,:) :: axpml, bxpml, aypml, bypml
   real*8, target, allocatable, dimension(:,:,:,:) :: scrpmlx, psipmlx, pmlvarsx,&
					scrpmly, psipmly, pmlvarsy
   real*8, pointer, dimension(:,:,:) :: psivx, psivy, RHSpsivx, RHSpsivy, psi_gradp_x,&
				RHS_psi_gradp_x, psi_gradp_y, RHS_psi_gradp_y
   logical :: PROC_HAS_PML


   !-----DAMPING STUFF---!
   real*8, allocatable, dimension(:,:,:) :: transfermatrix
   complex*16, allocatable, dimension(:,:,:) :: transfertemp
   real*8, allocatable, dimension(:,:) :: damping_rates,kay2d

   !------KERNEL STUFF-----!
 
   integer bcx, bcy, bcz, source_x, source_y, nmasters
   integer*8 :: fftw_time_fwd, fftw_time_inv
   logical :: PROC_HAS_ADJOINT, COMPUTE_DATA, COMPUTE_SYNTH
   real*8 powerspec_fac, loc_x, loc_y
   real*8, allocatable, dimension(:,:,:) :: f_xi_x, f_xi_y, f_xi_z, hessian, &
			 kernelc2, kernelrho, kernelp, f_vel_x, f_vel_y, f_vel_z, a_vel_x, &
			 a_vel_y, a_vel_z, f_acc_x, f_acc_y, f_acc_z, a_acc_x, a_acc_y, a_acc_z
   real*8, allocatable, dimension(:,:,:,:) :: kernelv, kernelb
   real*8, allocatable, dimension(:) :: sourcefunction, zkern, stretchkern, nus, grav, adjoint_source, fwdsource

   integer  st_adjoint, st_forward, end_forward, cadforcing_step
   real*8 local_time, final_time, fwd_start_time
   character*2 contrib,jobno


Contains

!==========================================================================


 Subroutine Initialize_all()

   implicit none

   integer i, j, k,fwd_sprev, inv_sprev, rem1, fwd_prev_z
   integer fwd_rprev, inv_rprev, ierr, rem2, rem3
   integer (KIND=MPI_ADDRESS_KIND) dum,sizeofdouble
   real*8, allocatable, dimension(:,:,:) :: in
   complex*16, allocatable, dimension(:,:,:) :: co1
   real*8 param1, param2, value1, value2, tempxyz, const, themax

   call MPI_INIT(ierr)
   call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
   call MPI_COMM_SIZE(MPI_COMM_WORLD, numtasks, ierr)
   call MPI_TYPE_GET_EXTENT(MPI_DOUBLE_PRECISION,dum,sizeofdouble,ierr)

   call PARSE_OPTIONS ()

   if (.NOT. (periodic)) then
      do i=1,nx
         x(i) = DBLE(i-1.0)/DBLE(nx - 1.0)
      enddo

      do i=1,ny
         y(i) = DBLE(i-1.0)/DBLE(ny - 1.0)
      enddo

   else
 
    do i=1,nx
     x(i) = DBLE(i-1.0)/DBLE(nx)
    enddo

    do i=1,ny
     y(i) = DBLE(i-1.0)/DBLE(ny)
    enddo

    allocate(eyekx(nx/2+1), eyeky(ny/2+1))

    do i=1,nx/2+1
     eyekx(i) = DBLE(i-1.0)*eye*normx*2.0*pi
    enddo

    do i=1,ny/2+1
     eyeky(i) = DBLE(i-1.0)*eye*normy*2.0*pi
    enddo

   endif

  
   do i=1,nx/2+1
      decay_coeff_x(i) = 1.0/(1.0 + exp((i-nx/3.0)*1.))*normx
   enddo
 
   do i=1,ny/2+1
      decay_coeff_y(i) = 1.0/(1.0 + exp((i-ny/3.0)*1.))*normy
   enddo

   stretchx = rsun/xlength
   stretchy = rsun/ylength

   allocate(dim1(0:numtasks-1), dim2(0:numtasks-1), dimz(0:numtasks-1),&
	fwd_sendtypes(0:numtasks-1), inv_sendtypes(0:numtasks-1), &
	fwd_recvtypes(0:numtasks-1), inv_recvtypes(0:numtasks-1), &
	ystart(0:numtasks-1), fwd_sdispls(0:numtasks-1), fwd_rdispls(0:numtasks-1),&
	inv_sdispls(0:numtasks-1), inv_rdispls(0:numtasks-1), fwd_z_send(0:numtasks-1),&
        inv_z_send(0:numtasks-1), fwd_z_recv(0:numtasks-1), inv_z_recv(0:numtasks-1),&
        fwd_z_displs(0:numtasks-1)) 

   dim1 = nx/numtasks
   dim2 = ny/numtasks
   dimz = nz/numtasks

   rem1 = nx - dim1(0)*numtasks
   rem2 = ny - dim2(0)*numtasks
   rem3 = nz - dimz(0)*numtasks

   do i=0,rem1-1
      dim1(i) = dim1(i) +1
   enddo

   do i=0,rem2-1
      dim2(i) = dim2(i) +1
   enddo

   do i=0,rem3-1
      dimz(i) = dimz(i) +1
   enddo

   ! Fwd transpose
   ! Sent packets are of dimensions dim1(receiver_rank), dim2(rank), nz
   ! Received packets are of dimensions dim2(sender_rank), dim1(rank), nz
   

   fwd_sprev = 1
   fwd_rprev = 1
   fwd_prev_z = 1
   ystart(0) = 1

   do j=0,numtasks-1

    fwd_sendtypes(j) = transp(nx,dim1(j),dim2(rank),nz)
    fwd_recvtypes(j) = datatype(ny,dim2(j),dim1(rank),nz)

    fwd_z_send(j) = datatype_z(nx,nx,dim2(rank),dimz(j))
    fwd_z_recv(j) = datatype_z(nx,nx,dim2(j),dimz(rank))

    fwd_sdispls(j) = fwd_sprev
    fwd_sprev = fwd_sprev + dim1(j)

    fwd_rdispls(j) = fwd_rprev
    fwd_rprev = fwd_rprev + dim2(j)

    fwd_z_displs(j) = fwd_prev_z
    fwd_prev_z = fwd_prev_z + dimz(j)

    if (j .NE. 0) ystart(j) = ystart(j-1) + dim2(j-1)
   enddo 


   ! Inv transpose
   ! Sent packets are of dimensions dim2(receiver_rank), dim1(rank), nz
   ! Received packets are of dimensions dim1(sender_rank), dim2(rank), nz
   
   inv_sprev = 1
   inv_rprev = 1

   do j=0,numtasks-1

      inv_sendtypes(j) = transp(ny,dim2(j),dim1(rank),nz)
      inv_recvtypes(j) = datatype(nx,dim1(j),dim2(rank),nz)

      inv_z_send(j) = datatype_z(nx,nx,dim2(j),dimz(rank))
      inv_z_recv(j) = datatype_z(nx,nx,dim2(rank),dimz(j))

      inv_sdispls(j) = inv_sprev
      inv_sprev = inv_sprev + dim2(j)

      inv_rdispls(j) = inv_rprev
      inv_rprev = inv_rprev + dim1(j)

   enddo 

   ! The size of the fourth dimension

   dimfive = 5
   if (magnetic) dimfive = 8
   if (displ) dimfive = 6
  

   if (kernel_mode) call read_parameters

   if (.not. CONSTRUCT_KERNELS) then

      !----COMPUTING FFT PLANS

      allocate(in(nx,dim2(rank),nz), co1(nx/2+1,dim2(rank),nz))

      call dfftw_plan_guru_dft_r2c(fftw_plan_fwd_x,1,nx,1,1,1,nz*dim2(rank),&
                   &            nx,(nx/2+1),in(1,1,1),co1(1,1,1),FFTW_MEASURE)
      call dfftw_plan_guru_dft_c2r(fftw_plan_inv_x,1,nx,1,1,1,nz*dim2(rank),&
                   &            (nx/2+1),nx,co1(1,1,1),in(1,1,1),FFTW_MEASURE)


      deallocate(in,co1)


      allocate(in(ny,dim1(rank),nz), co1(ny/2+1,dim1(rank),nz))

      call dfftw_plan_guru_dft_r2c(fftw_plan_fwd_y,1,ny,1,1,1,nz*dim1(rank),&
                   &            ny,(ny/2+1),in(1,1,1),co1(1,1,1),FFTW_MEASURE)
      call dfftw_plan_guru_dft_c2r(fftw_plan_inv_y,1,ny,1,1,1,nz*dim1(rank),&
                   &            (ny/2+1),ny,co1(1,1,1),in(1,1,1),FFTW_MEASURE)


      deallocate(in,co1)

      !-----   

      allocate(a(nx,dim2(rank),nz,dimfive), temp_step(nx,dim2(rank),nz,dimfive),&
                            scr(nx,dim2(rank),nz,dimfive))


      a = 0.0
      scr =0.0
      temp_step = 0.0

    

      if (magnetic) then 
   !        if (.not. displ) 
         if (displ) allocate(scratch(nx,dim2(rank),nz,5))

         allocate(box(nx,dim2(rank),nz), boy(nx, dim2(rank), nz), boz(nx, dim2(rank), nz),&
         curlbox(nx,dim2(rank),nz), curlboy(nx,dim2(rank),nz), curlboz(nx,dim2(rank),nz),&
         curlbx(nx,dim2(rank),nz), curlby(nx,dim2(rank),nz), curlbz(nx,dim2(rank),nz),&
         gradrho0_x(nx, dim2(rank), nz), gradrho0_y(nx, dim2(rank), nz), gradp0_x(nx, dim2(rank), nz),&
         gradp0_y(nx, dim2(rank), nz),reduction(nx,dim2(rank),nz),dzbx(nx,dim2(rank),nz), &
         dzflux2(nx,dim2(rank),nz), dzflux1(nx,dim2(rank),nz), dzby(nx,dim2(rank),nz),&
         flux1(nx,dim2(rank),nz), flux2(nx,dim2(rank),nz), flux3(nx,dim2(rank),nz), &
         erad_2d(nx,dim2(rank)), orad_2d(nx,dim2(rank)), delta_width_2d(nx, dim2(rank)))

         if (USE_PML) &
         allocate(scrpml(nx,dim2(rank),nzpml,6),pmlvars(nx,dim2(rank),nzpml,6),&
         psipml(nx,dim2(rank),nzpml,6))
    
      else
         !if (.not. displ) allocate(a(nx,dim2(rank),nz,5), temp_step(nx,dim2(rank),nz,5),&
         !!scr(nx,dim2(rank),nz,5), scrpml(nx,dim2(rank),nzpml,2), pmlvars(nx,dim2(rank),nzpml,2),&
         !psipml(nx,dim2(rank),nzpml,2))

         allocate(scrpml(nx,dim2(rank),nzpml,2), pmlvars(nx,dim2(rank),nzpml,2),&
         psipml(nx,dim2(rank),nzpml,2), scratch(nx, dim2(rank),nz,2))

      endif
      
      allocate(gradp(nx,dim2(rank),nz,3),g(nz), gamma(nz), source_dep(nz),&
           c2(nx,dim2(rank),nz), div(nx,dim2(rank),nz),&
      dvxdx(nx,dim2(rank),nz), forcing(nx,dim2(rank)), dvydy(nx,dim2(rank),nz),&
      spongexyz(nx, dim2(rank), nz), gradrho0_z(nx, dim2(rank), nz),&
      gradp0_z(nx, dim2(rank), nz), rhoinv(nx, dim2(rank), nz), rho0(nx, dim2(rank), nz), &
      p0(nx, dim2(rank), nz), c_speed(nx, dim2(rank), nz), c2rho0(nx,dim2(rank),nz),&
      z_i(nx,dim2(rank)), z_iplus1(nx,dim2(rank)), LC0(nx,dim2(rank)), LC1(nx,dim2(rank)), LC2(nx,dim2(rank)),&
      LC3(nx,dim2(rank)), LC4(nx,dim2(rank)), LC5(nx,dim2(rank)), LC6(nx,dim2(rank)),&
      az(nx,dim2(rank),nz),bzpml(nx,dim2(rank),nz),gradvz(nx,dim2(rank),nz))


      if (HORIZONTAL_PMLS) then

         allocate(axpml(2*npmlhor,dim2(rank),nz), bxpml(2*npmlhor,dim2(rank),nz),&
         psipmlx(2*npmlhor,dim2(rank),nz, 2), pmlvarsx(2*npmlhor,dim2(rank),nz,2),&
           scrpmlx(2*npmlhor, dim2(rank), nz, 2))

         psivx => psipmlx(:,:,:,1)
         psi_gradp_x => psipmlx(:,:,:,2)

         RHSpsivx => scrpmlx(:,:,:,1)
         RHS_psi_gradp_x => scrpmlx(:,:,:,2)



         PROC_HAS_PML = .FALSE.

         if ((ystart(rank) .LE. npmlhor) .OR. (ystart(rank) + dim2(rank)-1 .GE. ny-npmlhor+1 )) then

            PROC_HAS_PML = .TRUE.
            allocate(aypml(nx,dim2(rank),nz), bypml(nx,dim2(rank),nz),&
            psipmly(nx,dim2(rank),nz, 2), pmlvarsy(nx,dim2(rank),nz,2),&
            scrpmly(nx, dim2(rank), nz, 2))

            psivy => psipmly(:,:,:,1)
            psi_gradp_y => psipmly(:,:,:,2)

            RHSpsivy => scrpmly(:,:,:,1)
            RHS_psi_gradp_y => scrpmly(:,:,:,2)
            
         endif

      endif


      if (flows) then

         allocate(v0(nx,dim2(rank),nz,3))

         v0 = 0.0

         v0_x => v0(:,:,:,1)
         v0_y => v0(:,:,:,2)
         v0_z => v0(:,:,:,3)

         if (.not. displ) then
           allocate(omega(nx,dim2(rank),nz,3),omega0(nx,dim2(rank),nz,3),&
      div_v0(nx,dim2(rank),nz), advect0(nx,dim2(rank),nz,3), gradrho(nx,dim2(rank),nz,3),&
      flow_x(nx,dim2(rank),nz), flow_y(nx,dim2(rank),nz), flow_z(nx,dim2(rank),nz), c2div_v0(nx,dim2(rank),nz) )


           gradrho_x => gradrho(:,:,:,1)
           gradrho_y => gradrho(:,:,:,2)
           gradrho_z => gradrho(:,:,:,3)


           omega0_x => omega0(:,:,:,1)
           omega0_y => omega0(:,:,:,2)
           omega0_z => omega0(:,:,:,3)

           omega_x => omega(:,:,:,1)
           omega_y => omega(:,:,:,2)
           omega_z => omega(:,:,:,3)

           advect0_x => advect0(:,:,:,1)
           advect0_y => advect0(:,:,:,2)
           advect0_z => advect0(:,:,:,3)
          endif
       endif

      call solar_data ()
      ! Setup Radial derivative stuff

      call dbyd2(unstretch(1),z(1),1,nz,1)

      do k=1,nz
       stretch(k) = 1.0/unstretch(k)
      enddo

    
      ! Non-dimensional form, for the calculations
      deltat = timestep*dimc/diml

      if (rank == 0) &
         print *,'With cadence of 60 s, the timestep is', deltat*diml/dimc

      a = 0.0
      pmlvars = 0.0

      if (horizontal_pmls) then
         pmlvarsx = 0.0
         if (PROC_HAS_PML) pmlvarsy = 0.0
      endif


      psivz => psipml(:,:,:,1)
      psip => psipml(:,:,:,2)

      RHSpsivz => scrpml(:,:,:,1)
      RHSpsip => scrpml(:,:,:,2)

      gradp_x => gradp(:,:,:,1)
      gradp_y => gradp(:,:,:,2)
      gradp_z => gradp(:,:,:,3)

      if (.not. displ) then
    
         rho => temp_step(:,:,:,1)
         v_x => temp_step(:,:,:,2)
         v_y => temp_step(:,:,:,3)
         v_z => temp_step(:,:,:,4)
         p => temp_step(:,:,:,5)

         dvzdz => gradvz(:,:,:)

         RHScont => scr(:,:,:,1) ! Continuity eqn.
         RHSv_x  => scr(:,:,:,2) ! Radial momentum
         RHSv_y  => scr(:,:,:,3) ! Latitudinal momentum
         RHSv_z  => scr(:,:,:,4) ! Longitudinal momentum
         RHSp  => scr(:,:,:,5) ! Pressure eqn.
       
       
      elseif (displ) then

         xi_x => temp_step(:,:,:,1) ! Continuity eqn.
         xi_y  => temp_Step(:,:,:,2) ! Radial momentum
         xi_z  => temp_step(:,:,:,3) ! Latitudinal momentum
         v_x  => temp_step(:,:,:,4) ! Longitudinal momentum
         v_y  => temp_step(:,:,:,5) ! Pressure eqn.
         v_z  => temp_step(:,:,:,6) ! Pressure eqn.

         RHSxi_x => scr(:,:,:,1) ! Continuity eqn.
         RHSxi_y  => scr(:,:,:,2) ! Radial momentum
         RHSxi_z  => scr(:,:,:,3) ! Latitudinal momentum
         RHSv_x  => scr(:,:,:,4) ! Longitudinal momentum
         RHSv_y  => scr(:,:,:,5) ! Pressure eqn.
         RHSv_z  => scr(:,:,:,6) ! Pressure eqn.

         dxizdz => gradvz(:,:,:)
         dxixdx => dvxdx
         dxiydy => dvydy

         p => scratch(:,:,:,2)
         rho => scratch(:,:,:,1)
       
      endif

      if (magnetic) then 

      box = 0.0 
      boy = 0.0
      boz = 0.0

      if (.not. displ) then
         bx => temp_step(:,:,:,6)
         by => temp_step(:,:,:,7)
         bz => temp_step(:,:,:,8)

         RHSb_x => scr(:,:,:,6) ! equation for b_x
         RHSb_y => scr(:,:,:,7) ! equation for b_y
         RHSb_z => scr(:,:,:,8) ! equation for b_z
      endif

      if (USE_PML) then

         RHSpsidzbx => scrpml(:,:,:,3)
         RHSpsidzby => scrpml(:,:,:,4)
         RHSpsiinductionbx => scrpml(:,:,:,5)
         RHSpsiinductionby => scrpml(:,:,:,6)

         psidzbx => psipml(:,:,:,3)
         psidzby => psipml(:,:,:,4)
         psiinductionbx => psipml(:,:,:,5) 
         psiinductionby => psipml(:,:,:,6)

      endif

      if (displ) then
         bx => scratch(:,:,:,3) 
         by => scratch(:,:,:,4) 
         bz => scratch(:,:,:,5) 
      endif

      endif 

   else
      call INIT_KERNEL ()
   endif


   call MPI_BARRIER(MPI_COMM_WORLD, ierr)

end Subroutine Initialize_all


!==============================================================================

SUBROUTINE INIT_KERNEL
 
  implicit none
  integer i,j, rem1, rem2, rem3, fwd_sprev, fwd_rprev, inv_sprev, inv_rprev
  real*8, allocatable, dimension(:,:,:) :: in
  complex*16, allocatable, dimension(:,:,:) :: co1 
  real*8, allocatable, dimension(:,:,:,:) ::inn4
  complex*16, allocatable, dimension(:,:,:,:) :: fft_data
  real*8 unstretchkern(nz_kern)

   ! Fwd transpose
   ! Sent packets are of dimensions dim1(receiver_rank), dim2(rank), nz
   ! Received packets are of dimensions dim2(sender_rank), dim1(rank), nz
   
   fwd_sprev = 1
   fwd_rprev = 1
!   fwd_prev_z = 1
   ystart(0) = 1

   do j=0,numtasks-1

    fwd_sendtypes(j) = transp(nx,dim1(j),dim2(rank),nz_kern)
    fwd_recvtypes(j) = datatype(ny,dim2(j),dim1(rank),nz_kern)

 !   fwd_z_send(j) = datatype_z(nx,nx,dim2(rank),dimz(j))
  !  fwd_z_recv(j) = datatype_z(nx,nx,dim2(j),dimz(rank))

    fwd_sdispls(j) = fwd_sprev
    fwd_sprev = fwd_sprev + dim1(j)

    fwd_rdispls(j) = fwd_rprev
    fwd_rprev = fwd_rprev + dim2(j)

!    fwd_z_displs(j) = fwd_prev_z
 !   fwd_prev_z = fwd_prev_z + dimz(j)

    if (j .NE. 0) ystart(j) = ystart(j-1) + dim2(j-1)
   enddo 


   ! Inv transpose
   ! Sent packets are of dimensions dim2(receiver_rank), dim1(rank), nz_kern
   ! Received packets are of dimensions dim1(sender_rank), dim2(rank), nz_kern
   
   inv_sprev = 1
   inv_rprev = 1

   do j=0,numtasks-1

    inv_sendtypes(j) = transp(ny,dim2(j),dim1(rank),nz_kern)
    inv_recvtypes(j) = datatype(nx,dim1(j),dim2(rank),nz_kern)
  
!    inv_z_send(j) = datatype_z(nx,nx,dim2(j),dimz(rank))
 !   inv_z_recv(j) = datatype_z(nx,nx,dim2(rank),dimz(j))

    inv_sdispls(j) = inv_sprev
    inv_sprev = inv_sprev + dim2(j)

    inv_rdispls(j) = inv_rprev
    inv_rprev = inv_rprev + dim1(j)

  enddo 

 
  
  allocate(f_xi_x(nx,dim2(rank),nz_kern), f_xi_y(nx,dim2(rank),nz_kern), &
	   f_xi_z(nx,dim2(rank),nz_kern),  &
	   kernelc2(nx,dim2(rank),nz_kern), kernelrho(nx,dim2(rank),nz_kern), kernelp(nx,dim2(rank),nz_kern), &
	   kernelv(nx,dim2(rank),nz_kern,3), nus(nt_kern), &
	   sourcefunction(nt_kern/2+1), zkern(nz_kern), f_vel_x(nx,dim2(rank),nz_kern),&
	   f_vel_y(nx,dim2(rank),nz_kern), f_vel_z(nx,dim2(rank),nz_kern), &
	   a_vel_x(nx,dim2(rank),nz_kern),&
	   a_vel_y(nx,dim2(rank),nz_kern), a_vel_z(nx,dim2(rank),nz_kern), &
	   stretchkern(nz_kern), grav(nz_kern), gamma(nz_kern),&
	   gradrho0_x(nx,dim2(rank),nz_kern), gradrho0_y(nx,dim2(rank),nz_kern), &
	   gradrho0_z(nx,dim2(rank),nz_kern), f_acc_x(nx,dim2(rank),nz_kern), &
	   f_acc_y(nx,dim2(rank),nz_kern),  f_acc_z(nx,dim2(rank),nz_kern),&
	   a_acc_x(nx,dim2(rank),nz_kern), a_acc_y(nx,dim2(rank),nz_kern), &
	   a_acc_z(nx,dim2(rank),nz_kern), hessian(nx,dim2(rank),nz_kern))


  allocate(p0(nx, dim2(rank), nz_kern), rho0(nx, dim2(rank), nz_kern), c2(nx, dim2(rank), nz_kern),c_speed(nx,dim2(rank),nz_kern))

  if (magnetic) allocate(curlbox(nx, dim2(rank), nz_kern), curlboy(nx, dim2(rank), nz_kern), &
        curlboz(nx,dim2(rank),nz_kern), box(nx, dim2(rank), nz_kern), boy(nx, dim2(rank), nz_kern), &
	boz(nx, dim2(rank), nz_kern), kernelb(nx,dim2(rank),nz_kern,3))

  if (background_flows_exist) allocate(v0_x(nx, dim2(rank), nz_kern), v0_y(nx, dim2(rank), nz_kern), &
	v0_z(nx, dim2(rank), nz_kern))


   !----COMPUTING FFT PLANS



   allocate(in(nx,dim2(rank),nz_kern), co1(nx/2+1,dim2(rank),nz_kern))

   call dfftw_plan_guru_dft_r2c(fftw_plan_fwd_x,1,nx,1,1,1,nz_kern*dim2(rank),&
                &            nx,(nx/2+1),in(1,1,1),co1(1,1,1),FFTW_ESTIMATE)
   call dfftw_plan_guru_dft_c2r(fftw_plan_inv_x,1,nx,1,1,1,nz_kern*dim2(rank),&
                &            (nx/2+1),nx,co1(1,1,1),in(1,1,1),FFTW_ESTIMATE)


   deallocate(in,co1)

   allocate(in(ny,dim1(rank),nz_kern), co1(ny/2+1,dim1(rank),nz_kern))

   call dfftw_plan_guru_dft_r2c(fftw_plan_fwd_y,1,ny,1,1,1,nz_kern*dim1(rank),&
                &            ny,(ny/2+1),in(1,1,1),co1(1,1,1),FFTW_ESTIMATE)
   call dfftw_plan_guru_dft_c2r(fftw_plan_inv_y,1,ny,1,1,1,nz_kern*dim1(rank),&
                &            (ny/2+1),ny,co1(1,1,1),in(1,1,1),FFTW_ESTIMATE)


   deallocate(in,co1)

    !-----   

!   allocate(inn4(nx, dim2(rank), nz_kern, nt_kern))
!  call dfftw_plan_guru_dft_r2c(fftw_time_fwd, 1, nt_kern, nx*dim2(rank)*nz_kern,  nx*dim2(rank)*nz_kern, 1, &
!			nx*dim2(rank)*nz_kern,1,1,inn4, fft_data, FFTW_ESTIMATE)
                                
!  call dfftw_plan_guru_dft_c2r(fftw_time_inv, 1, nt_kern, nx*dim2(rank)*nz_kern,  nx*dim2(rank)*nz_kern, 1, &
!			nx*dim2(rank)*nz_kern,1,1, fft_data, inn4,  FFTW_ESTIMATE)

!  deallocate(fft_data, inn4)

  f_xi_x = 0.0

  call kernel_back ()
  call dbyd2(unstretchkern(1),zkern(1),1,nz_kern,1)
  stretchkern = 1.0/unstretchkern


END SUBROUTINE INIT_KERNEL



!==============================================================================


Subroutine solar_data

  implicit none
  integer k,q
  real*8 data(nz,6),temp


  ! Data in the file is arranged as 
  ! Non-dimensional Solar radius, sound speed, density, pressure, gravity, gamma_1
  ! In CGS Units

  ! solar mass: 1.9891d33 g
 
  open(7,file = file_data,form = 'formatted',status = 'old')
  do k = 1,nz
    read(7,*) data(k,:)
  enddo
  close(7)

  q = 1
  
! dimc = 10 km/s

  dimc = 10.0**6. !data(q,2)

! dimrho = 10^-3 g/cc
 
  dimrho = 10.0**(-3.0) !data(q,3)
  diml = rsun
  dimb = (4.0*pi*dimc**2.0*dimrho)**0.5

  do k =1,nz
    z(k) = data(k,1)
    c_speed(:,:,k) = data(k,2)/dimc  !*0.01
    c2(:,:,k) = (data(k,2)/dimc)**2.0 !c_speed(4,5,k)**2.0
    rho0(:,:,k) = data(k,3)/dimrho
    p0(:,:,k) = data(k,4)/(dimrho*dimc**2.0)
    g(k) = data(k,5)*diml/dimc**2.0
    gamma(k) = data(k,6)
    cq(k) =data(k,2)/dimc
    rhoq(k) = data(k,3)/dimrho
  enddo

  height = (z-1.)*695.98994

 
end Subroutine solar_data

!==========================================================================

SUBROUTINE PARSE_OPTIONS
 implicit none
 

 ! USE_PML OPTIONS ARE ALL EVEN
 ! 

 if ((.not. magnetic) .and. (.not. flows)) then
   if ((.not. TEST_IN_2D) .and. (.not. USE_PML) .and. (.not. displ)) option = 1
   if ((.not. TEST_IN_2D) .and. (USE_PML) .and. (.not. displ)) option = 2

   if ((.not. TEST_IN_2D) .and. (.not. USE_PML) .and. displ) option = 3
   if ((.not. TEST_IN_2D) .and. (USE_PML) .and. displ) option = 4

   if ((TEST_IN_2D) .and. (.not. USE_PML)) option = 5
   if ((TEST_IN_2D) .and. (USE_PML)) option = 6


   if ((TEST_IN_2D) .and. (USE_PML) .and. (KERNEL_MODE)) option = 17
 endif

 if (magnetic) then
   if ((.not. TEST_IN_2D) .and. (.not. USE_PML) .and. (.not. displ)) option = 7
   if ((.not. TEST_IN_2D) .and. (USE_PML) .and. (.not. displ)) option = 8

   if ((.not. TEST_IN_2D) .and. (.not. USE_PML) .and. displ) option = 9
   if ((.not. TEST_IN_2D) .and. (USE_PML) .and. displ) option = 10

   if ((TEST_IN_2D) .and. (.not. USE_PML)) option = 11
   if ((TEST_IN_2D) .and. (USE_PML)) option = 12

   if ((TEST_IN_2D) .and. (USE_PML) .and. (KERNEL_MODE)) option = 18
 endif

 if (flows) then
   if (.not. USE_PML .and. (.not. displ)) option = 13
   if (USE_PML .and. (.not. displ)) option = 14

   if (.not. USE_PML .and. displ) option = 15
   if (USE_PML .and. displ) option = 16

 endif




END SUBROUTINE PARSE_OPTIONS

!==========================================================================

FUNCTION TRANSP_COMP(sizeofcol,m,n,p)
   implicit none
   integer transp_comp,ierr,row,m,n,sizeofcol,p,transp2d
   integer(KIND=MPI_ADDRESS_KIND) extent2d, dum, extent

! Datatype for transpose of an m X n X p matrix
   call MPI_TYPE_GET_EXTENT(MPI_DOUBLE_COMPLEX,dum, extent,ierr)
   extent2d = sizeofcol*n*extent
   call MPI_TYPE_VECTOR(n,1,sizeofcol,MPI_DOUBLE_COMPLEX,row,ierr)
   call MPI_TYPE_CREATE_HVECTOR(m,1,extent,row,transp2d,ierr)
   call MPI_TYPE_CREATE_HVECTOR(p,1,extent2d,transp2d,transp_comp,ierr)
   call MPI_TYPE_COMMIT(transp_comp,ierr)

END FUNCTION TRANSP_COMP
!==========================================================================

FUNCTION DATATYPE_COMP(sizeofcol,m,n,p)
   implicit none
   integer datatype_comp,ierr,newtype,m,n,sizeofcol,p,datatype2d
   integer (KIND=MPI_ADDRESS_KIND) extentofcol, dum, extent2d, extent

! Data type for a submatrix m X n X p in a supermatrix M X n X p
   call MPI_TYPE_GET_EXTENT(MPI_DOUBLE_COMPLEX,dum,extent,ierr)
   extentofcol = sizeofcol*extent
   extent2d = sizeofcol*n*extent

   call MPI_TYPE_VECTOR(m,1,1,MPI_DOUBLE_COMPLEX,newtype,ierr)
   call MPI_TYPE_CREATE_HVECTOR(n,1,extentofcol,newtype,datatype2d,ierr)
   call MPI_TYPE_CREATE_HVECTOR(p,1,extent2d,datatype2d,datatype_comp,ierr)
   call MPI_TYPE_COMMIT(datatype_comp,ierr)

END FUNCTION DATATYPE_COMP
!==========================================================================

FUNCTION DATATYPE_Z(sizeofcol,m,n,p)
   implicit none
   integer datatype_z,ierr,newtype,m,n,sizeofcol,p,datatype2d
   integer (KIND=MPI_ADDRESS_KIND) extentofcol, dum, extent2d, extent

! Data type for a submatrix m X n X p in a supermatrix m X n X P
   call MPI_TYPE_GET_EXTENT(MPI_DOUBLE_PRECISION,dum,extent,ierr)
   extentofcol = sizeofcol*extent
   extent2d = sizeofcol*n*extent

   call MPI_TYPE_VECTOR(m,1,1,MPI_DOUBLE_PRECISION,newtype,ierr)
   call MPI_TYPE_CREATE_HVECTOR(n,1,extentofcol,newtype,datatype2d,ierr)
   call MPI_TYPE_CREATE_HVECTOR(p,1,extent2d,datatype2d,datatype_z,ierr)
   call MPI_TYPE_COMMIT(datatype_z,ierr)

END FUNCTION DATATYPE_Z
!==========================================================================

FUNCTION TRANSP(sizeofcol,m,n,p)
   implicit none
   integer transp,ierr,row,m,n,sizeofcol,p,transp2d
   integer(KIND=MPI_ADDRESS_KIND) extent2d, dum, extent

! Datatype for transpose of an m X n X p matrix
   call MPI_TYPE_GET_EXTENT(MPI_DOUBLE_PRECISION, dum, extent,ierr)
   extent2d = sizeofcol*n*extent
   call MPI_TYPE_VECTOR(n,1,sizeofcol,MPI_DOUBLE_PRECISION,row,ierr)
   call MPI_TYPE_CREATE_HVECTOR(m,1,extent,row,transp2d,ierr)
   call MPI_TYPE_CREATE_HVECTOR(p,1,extent2d,transp2d,transp,ierr)
   call MPI_TYPE_COMMIT(transp,ierr)

END FUNCTION TRANSP
!==========================================================================

FUNCTION DATATYPE(sizeofcol,m,n,p)
   implicit none
   integer datatype,ierr,newtype,m,n,sizeofcol,p,datatype2d
   integer (KIND=MPI_ADDRESS_KIND) extentofcol, dum, extent2d,extent

! Data type for a submatrix m X n X p in a supermatrix M X n X p
   call MPI_TYPE_GET_EXTENT(MPI_DOUBLE_PRECISION,dum, extent,ierr)
   extentofcol = sizeofcol*extent
   extent2d = sizeofcol*n*extent

   call MPI_TYPE_VECTOR(m,1,1,MPI_DOUBLE_PRECISION,newtype,ierr)
   call MPI_TYPE_CREATE_HVECTOR(n,1,extentofcol,newtype,datatype2d,ierr)
   call MPI_TYPE_CREATE_HVECTOR(p,1,extent2d,datatype2d,datatype,ierr)
   call MPI_TYPE_COMMIT(datatype,ierr)

END FUNCTION DATATYPE


!==========================================================================


Subroutine kernel_back ()

  implicit none
  integer k,q,kk
  real*8 data(nz,6),temp


  ! Data in the file is arranged as 
  ! Non-dimensional Solar radius, sound speed, density, pressure, gravity, gamma_1
  ! In CGS Units

  ! solar mass: 1.9891d33 g
 
  open(7,file = file_data,form = 'formatted',status = 'old')
  do k = 1,nz
    read(7,*) data(k,:)
  enddo
  close(7)

  q = 1
! dimc = 10 km/s

  dimc = 10.0**6. !data(q,2)

! dimrho = 10^-3 g/cc
 
  dimrho = 10.0**(-3.0) !data(q,3)
  diml = rsun
  dimb = (4.0*pi*dimc**2.0*dimrho)**0.5

  do k =1,nz_kern
    kk = k + st_z - 1
    zkern(k) = data(kk,1)
    !c_speed(:,:,k) = data(kk,2)/dimc  !*0.01
    c2(:,:,k) = (data(kk,2)/dimc)**2.0 !c_speed(4,5,k)**2.0 /dimc
    cq(k) = data(kk,2)/dimc
    rho0(:,:,k) = data(kk,3)/dimrho
    p0(:,:,k) = data(kk,4)/(dimrho*dimc**2.0)
    grav(k) = data(kk,5)*diml/dimc**2.0
    gamma(k) = data(kk,6)
  enddo

end Subroutine kernel_back

!==========================================================================


SUBROUTINE READ_PARAMETERS

 implicit none
 integer i
 character*80 calculation_type, directory_rel
 character*2 whether_2D, contribs
 integer ierr
 logical lexist1, lexist2, lexist
 logical, allocatable, dimension(:) :: forward, adjoint, kernels

 compute_forward = .false.
 compute_adjoint = .false.
 construct_kernels = .false.
 generate_wavefield = .false.
 cadforcing_step = FLOOR(cadforcing/timestep)
 maxtime = floor((solartime*3600. + 2.0/nupeak)/timestep) + 2
 maxtime = cadforcing_step * floor(maxtime/cadforcing_step*1.) 
 forcing_length = floor(maxtime/cadforcing_step*1.) + 1

 call getarg(1,contrib)
 call getarg(2,jobno)

 inquire(file='Instruction_src'//contrib//'_ls'//jobno,exist=lexist1)

 open(356,file=directory//'master.pixels',action='read')
 do i=1,100
  read(356,*,end=24432)
 enddo
 24432 close(356)
 nmasters = i-1

 allocate(forward(nmasters), adjoint(nmasters), kernels(nmasters))

! do i=1,nmasters
!  call convert_to_string(i,contribs,2) 
!  inquire(file=directory//'status',exist=lexist)
!  if (.not. lexist .and. (rank==0)) call system('mkdir '//directory//'status')
!  inquire(file=directory//'status/forward'//contribs,exist=forward(i))
!  inquire(file=directory//'status/adjoint'//contribs,exist=adjoint(i))
!  inquire(file=directory//'status/kernel'//contribs,exist=kernels(i))
! enddo
 
if (lexist1) then

 open(22, file = 'Instruction_src'//contrib//'_ls'//jobno, form='formatted', status='unknown', action='read')
  read(22,"(A)") calculation_type
!  read(22,"(A)") directory_rel
 ! read(22,"(A)") contrib 

  if (adjustl(trim(calculation_type)) == 'spectral') then
    COMPUTE_ADJOINT = .FALSE.
    GENERATE_WAVEFIELD = .TRUE.
    COMPUTE_FORWARD = .TRUE.
    !COMPUTE_SYNTH = .FALSE.
    inquire(file='compute_synth',exist=COMPUTE_SYNTH)
    !LINESEARCH = .FALSE.
    inquire(file='linesearch',exist=LINESEARCH)
    COMPUTE_DATA = .false.
    inquire(file='compute_data', exist = COMPUTE_DATA)
    !if (lexist) LINESEARCH = .TRUE.
!    i = 1
!    do while (forward(i) .and. (i .le. nmasters))
! 	i = i+1
!    enddo
!    if (i .le. nmasters) then 
!      call convert_to_string(i,contrib,2)
!      call getarg(1, contrib) 
      inquire(file=directory//'forward_src'//contrib//'_ls'//jobno,exist=lexist)
      if (.not. lexist .and. (rank==0)) call system('mkdir '&
                //directory//'forward_src'//contrib//'_ls'//jobno)
!    elseif (i .gt. nmasters) then
!      print *,'All forward calculations processed'
!      stop
!    endif
     
    directory_rel = directory//'forward_src'//contrib//'_ls'//jobno//'/'


    timeline = -solartime * 3600.!timeline
  endif



  if (adjustl(trim(calculation_type)) == 'adjoint') then
     COMPUTE_ADJOINT = .TRUE.
     !i = 1

!      call getarg(1, contrib) 
!      do while (adjoint(i) .and. (i .le. nmasters))
! 	i = i+1
!     enddo
!     if (i .le. nmasters) then
!      call convert_to_string(i,contrib,2)
      inquire(file=directory//'adjoint_src'//contrib,exist=lexist)
      if (.not. lexist .and. (rank==0)) call system('mkdir '//directory//'adjoint_src'//contrib)
!     elseif (i .gt. nmasters) then
!      print *,'All adjoint calculations processed'
!      stop
!     endif
     directory_rel = directory//'adjoint_src'//contrib//'/'

     final_time = solartime*3600.
     timeline = final_time - (forcing_length - 1.)*cadforcing_step*timestep

 endif

  read(22,*) whether_2D


 if (whether_2D == '2D') TEST_IN_2D = .true.
  if (generate_wavefield) read(22,*) initcond
 close(22)


elseif (.not. lexist1) then

  CONSTRUCT_KERNELS = .true.
 

!     i = 1
!      do while (kernels(i) .and. (i .le. nmasters))
! 	i = i+1
!     enddo
!     if (i .le. nmasters) then
!      call convert_to_string(i,contrib,2)
!      call getarg(1, contrib) 
      inquire(file=directory//'kernel',exist=lexist)
      if (.not. lexist .and. (rank==0)) call system('mkdir '//directory//'kernel')
!     elseif (i .gt. nmasters) then
!      print *,'All kernels computed'

!stop
 !    endif
     directory_rel = directory//'kernel/'

endif

!if (.not. CONSTRUCT_KERNELS) then 


!  if (rank==0) print *,'Reading information file: '//adjustl(trim(directory_rel))//'information'

!  open(44,file=adjustl(trim(directory_rel))//'information',form='formatted',status='unknown',position='rewind')

 
!  read(44,*) forcing_length
!  read(44,*) timeline 
!  read(44,*) final_time
  
!endif

 if (test_in_2d) then
   if (.not. magnetic) option = 17
   if (magnetic) option = 18
 endif


deallocate(kernels, forward, adjoint)

END SUBROUTINE READ_PARAMETERS


!================================================================================

 SUBROUTINE convert_to_string(numbe,sting,length_string)

  implicit none
  integer i,length_string,numbe,n(1:length_string),number_temp
  character*(length_string) sting
  character*1 charc(10)

  charc(1)  = '0'
  charc(2)  = '1'
  charc(3)  = '2'
  charc(4)  = '3'
  charc(5)  = '4'
  charc(6)  = '5'
  charc(7)  = '6'
  charc(8)  = '7'
  charc(9)  = '8'
  charc(10) = '9'


  number_temp = numbe
  do i=length_string,1,-1

    n(length_string-i+1) = floor(number_temp*10.0**(-(i-1.0)))
    number_temp = number_temp - n(length_string-i+1)*10**(i-1)
    sting(length_string-i+1:length_string-i+1) = charc(n(length_string-i+1)+1)

  enddo

  end SUBROUTINE convert_to_string

!================================================================================

SUBROUTINE distmat(n,m,f)
 
 implicit none
 integer m, n, i, j, i2, j2
 real*8 f(n,m), sig


 do j=1,m
   do i=1,n
       sig = 1.0
       i2 = min(i-1,n-i+1)
       j2 = min(j-1,m-j+1)
       if ((i-1) > (n-i+1)) sig = -1.0
       f(i,j) = (i2**2 + j2**2)**0.5 * sig
   enddo
 enddo
 
END SUBROUTINE distmat

!================================================================================

END MODULE INITIALIZE
