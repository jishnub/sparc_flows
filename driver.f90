Program driver

    ! --------------------------------------------------------------------------
    ! MPI Version of the Cartesian Acoustic Sun Simulator.
    ! Copyright 2006, Shravan Hanasoge

    ! Hansen Experimental Physics Laboratory
    ! 455 Via Palou way, Stanford
    ! CA 94305, USA
    ! Email: shravan@stanford.edu

    ! The driver routine. It utilizes the other subroutines to perform time-
    ! advancements. The time advancement is performed with the explicit Low
    ! Dissipation, Low Dispersion Runge Kutta 4th order method, described in
    ! Hu et al. (1996). Radial derivatives are computed using compact finite
    ! difference (Lele 1992). Horizontal variations are extracted by projecting
    ! the desired quantity onto Fourier harmonic space followed by an
    ! appropriate inverse fourier transform.


    !
    ! BASIC CHECKS:
    ! 1. If init is all right.
    ! 2. If the time-length of the simulation matches
    ! 3. If the directories are OK
    ! 4. If the forcing function is being read from the right directory
    ! 5. If the variables are initialized OK
    ! 6. If the flow profiles are allright
    ! 7. If output and status are deleted or appropriately moved
    !
    ! ---------------------------------------------------------------------------
    use initialize
    use all_modules
    use mp_physics
    use mp_physics_2d
    use time_step
    use kernels
    use bspline
    use integrals
    use derivatives

    implicit none

    integer i,j, init, ierr, t, k, index(1000,2), randinit(2), aind, nsteps_given,reclmax, loc
    real(kind=real64) start_time,end_time,tempf,T00,nor1,nor2,nor3,nor4,nor5,e, mp_time,Rchar
    real(kind=real64) total_time, start_mp_time, avg_time,zz, tempxy, rand1, rand2,con,kay,z0,sigmaz,sigmax
    real(kind=real64) bes(0:2), Lregular,signt, xcutoffpix, xcutoff
    logical saved, iteration, init_variables, tempbool
    character(len=1) ci
    real(kind=real64) dt

20  format(8f12.4)


    ! Initializing the computation
    ! --------------------------------------------

    call Initialize_all

!~     call adjoint_source_filt(520)
!~     stop

    Lregular = 30.0D8/diml

    if (COMPUTE_DATA) then

        if (FLOWS) then

           call readfits('true_vx.fits',v0_x,nz)
           v0_x = v0_x*1D2/dimc

           call readfits('true_vz.fits',v0_z,nz)
           v0_z = v0_z*1D2/dimc


            if (contrib=="01") then
                print *,"max vx",maxval(v0_x)*dimc*1D-2,"max vz",maxval(v0_z)*dimc*1D-2," m/s"
            endif

            !~             CONTINUITY
            call continuity_check(v0_x,v0_z)

        endif


    else

        if (FLOWS) then
            v0_x = 0.0
            v0_z = 0.0
        endif


        if (FLOWS) then
            if (psi_cont .and. enf_cont) then
                inquire(file=directory//'model_psi_ls'//jobno//'.fits', exist = iteration)
            
            endif

            if (iteration) then

!~                  These logical variables are defined in params.i
                if (psi_cont .and. enf_cont) then

                    allocate(psivar(nx,dim2(rank),nz))
                    call readfits(directory//'model_psi_ls'//jobno//'.fits',psivar,nz)

                    psivar = rho0*c2**0.5*(psivar/(diml/1.0D8))

                    if (.not. CONSTRUCT_KERNELS) then
                        call ddz(psivar, v0_x, 1)
                        v0_x = -v0_x/rho0

                        call ddx(psivar, v0_z, 1)
                        v0_z = v0_z/rho0
                    else
                        call ddzkern(psivar, v0_x, 1)
                        v0_x = -v0_x/rho0

                        call ddxkern(psivar, v0_z, 1)
                        v0_z = v0_z/rho0
                    endif

                    deallocate(psivar)

                endif

                print *,"Source ",contrib
                if (contrib=="01") then

                    print *, "vxmax",maxval(abs(v0_x)*dimc*10.**(-2.)),"m/s " &
                    ,"vzmax", maxval(abs(v0_z)*dimc*10.**(-2.)),"m/s"

                end if

                !~             CONTINUITY
                call continuity_check(v0_x,v0_z)



            endif


            if (compute_adjoint .and. FLOWS) then
                v0_x = -v0_x
                v0_z = -v0_z
            endif
        endif
    endif
!~      stop
    if (CONSTRUCT_KERNELS) then 
        call PRODUCE_KERNELS
    else


    start_mp_time = MPI_WTIME()
    start_time = MPI_WTIME()

    !!!!!!! SOUND-SPEED PERTURBATIONS HERE !!!!!!

    if (.not. kernel_mode) then
        maxtime = floor(solartime*3600.0/timestep) + 1
        maxtime = 2*(maxtime/2)
    endif

    if (.not. TEST_IN_2D) call mp_initialize_RHS
    if (TEST_IN_2D) call mp_initialize_RHS_2d

    call Initialize_step

    end_time = MPI_WTIME()

    if (rank == 0) &
        print *, 'Initialization took: ', (end_time -start_time)

    !---------------------------------------------


    if ((time0 .ne.0) .and. (.not. kernel_mode)) &
        call read_in_initial_condition(directory, time0)

    !c2 = (gradp0_x - boz*curlboy) !*diml*10.**(-8.) /( dimc**2. * dimrho) !
    !c2 = (box**2. + boz**2.)**0.5/rho0**0.5 * dimc * 10.**(-5.)
    !do k=1,nz

    !c2(:,1,k) = (gradp0_z(:,1,k) + rho0(:,1,k)*g(k) + boz(:,1,k)*curlboy(:,1,k))!/(rho0(1,1,k)*g(k))
    !c2(:,1,k) = g(k)*(gradp0_z(:,1,k)/c2(:,1,k) -gradrho0_z(:,1,k))/rho0(:,1,k) * (dimc/diml)**2.!+ rho0(:,1,k)*g(k) + boz(:,1,k)*curlboy(:,1,k))!/(rho0(1,1,k)*g(k))
    !enddo
    !c2 = (box**2. + boz**2.)**0.5/(rho0)**0.5 * 10.0
    !c2(:,:,1) = 0.
    !c2(:,:,nz) = 0.
    !call writefits_3d('ca.fits',c2,nz)!p0*dimc**2.*dimrho,nz)
    !stop
    !do k=1,nz
    ! c2(:,1,k) = c2(:,1,k) * 1.01!(1. + 0.01*exp(-(x -0.5 + 10./300.)**2./(64.*(x(2)-x(1))**2.)))
    !enddo
    !c2rho0 = c2*rho0
    !c_speed = c2**0.5
    ! c2 = c2*1.01 !* 1.01
    ! c2rho0 = c2rho0*1.01 !* 1.01
    ! c_speed = c2**0.5

    !!!! OTHER INITIAL CONDITIONS!!!!!!
    ! a(:,:,:,1) - density
    ! a(:,:,:,2) - v_x
    ! and so on - v_y, v_z, p, b_x, b_y, b_z
    ! a(:,:,1,1) = orad_2d
    ! call writefits_3d('orad_2d.fits',a(:,:,1,1),1)
    !stop
    ! Length of simulation

    total_time = wall_time * 3600.

    steps = floor(outputcad/timestep )
    freq_filtering = floor(60./timestep)

    if (rank ==0) then
        print *,'FREQUENCY OF OUTPUT AT', steps, 'TIMESTEPS'
        print *, 'TIMESTEP =', timestep
        !    print *,'XXXXXX ------ NOT WRITING OUT XI VARIABLES -------- XXXXXXXX'
        print *,'OBSERVATION GRID POINT: ', o_rad, 'AND RADIUS:', z(o_Rad)
    endif

    ! Excitation series computation -----------------------------
    num_steps = FLOOR(DBLE(maxtime)*timestep/cadforcing) + 2
    num_steps = 2*(num_steps/2) + 2
    saved = .false.
    cadforcing_step = FLOOR(cadforcing/timestep)

    !   if (kernel_mode) then
    !    maxtime = floor((final_time - timeline)/timestep)
    !    maxtime = cadforcing_step * floor(maxtime/cadforcing_step*1.)
    !    nsteps_given = floor(maxtime/cadforcing_step*1.) + 1
    !   endif

    if (time0 == 0) init = 0
    if (time0 .NE. 0) init = time0+1

    if (rank ==0)  &
        print *, 'NUMBER OF TIMESTEPS =', maxtime
    if (kernel_mode) then
        call DETERMINE_STATUS(init, maxtime)
        if (rank==0) call system('rm Instruction_src'//contrib//'_ls'//jobno)
        if (compute_forward) then
            indexglob = 0
            if (rank == 0) then
                open(19,file=directory//'forward_src'//contrib//'_ls'//jobno//'/rms_hist'&
                        ,status='replace',form='formatted',position='append')
                open(745, file=directory//'forward_src'//contrib//'_ls'//jobno//'/kernel_info'&
                        ,status='unknown', action='write')
                write(745,*) init
                write(745,*) maxtime
                close(745)
            endif

            delta_width = 2.0d0/(z(e_rad+1) - z(e_rad-1))

            if (magnetic) then
                do k=1,dim2(rank)
                    do j=1,nx
                        delta_width_2d(j,k) = 2.0d0/(z(erad_2d(j,k)+1) - z(erad_2d(j,k)-1))
                    enddo
                enddo
            endif
            call timestepping(init)
            if (rank==0) then
                close(28)
                close(1244)
            endif

            call read_binary_writefits(directory//'forward_src'//contrib//'_ls'//jobno//'/vz_cc.bin',&
              indexglob, directory//'forward_src'//contrib//'_ls'//jobno//'/vz_cc.fits')

            if (.not. compute_data) then
                call adjoint_source_filt(indexglob)
                call misfit_all(indexglob)
            endif
        endif
        if (compute_adjoint) then
            if (rank == 0) then
                open(19,file=directory//'adjoint_src'//contrib//'/rms_hist',&
                        status='replace',form='formatted',position='append')
                open(745, file=directory//'adjoint_src'//contrib//'/kernel_info',&
                        status='unknown', action='write')
                write(745,*) floor(maxtime/cadforcing_step*1.) + 1
            endif

            call timestepping(init)

            if (rank==0) then
                write(745,*) floor(time/cadforcing_step*1.)*cadforcing_step
                close(745)
                close(28)
                close(1244)
            endif

            open(44,file=directory//'adjoint_src'//contrib//'/kernel_info',form='formatted',status='unknown')
            read(44,*) init!nt_kern
            read(44,*) maxtime
            close(44)
        endif


    endif  ! CORRESPONDING TO "IF KERNEL_MODE"


    call MPI_BARRIER(MPI_COMM_WORLD, ierr)

    if (.not. kernel_mode) call write_out_full_state(directory)
    ! if (kernel_mode) then
    !!  if (compute_adjoint) call write_out_full_state(directory//'adjoint'//contrib//'/')
    !  if (compute_forward) call write_out_full_state(directory//'forward'//contrib//'/')
    ! endif

    call dfftw_destroy_plan(fftw_plan_fwd_x)
    call dfftw_destroy_plan(fftw_plan_inv_x)

    call dfftw_destroy_plan(fftw_plan_fwd_y)
    call dfftw_destroy_plan(fftw_plan_inv_y)

    if (rank == 0) then
        close(19)
    endif

    endif ! if kernel

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    call MPI_FINALIZE(ierr)


end program driver

! --------------------------------------------------------------------------------------------



SUBROUTINE TIMESTEPPING(init)

 use time_step
 use mp_physics
 use mp_physics_2d
 use initialize
 use all_modules

 implicit none

 integer i,j,k, init, nsteps_given
 real(kind=real64) start_time, end_time, T00, t1, nor1
 real(kind=real64) tempstore(nx, dim2(rank))

  do i = init,maxtime
    time = i
  ! code to be timed

    call cpu_time(t1)
    T00 = t1
    start_time = MPI_WTIME()
    if (.not. kernel_mode) call lagrange_interp
    if ( kernel_mode .and. compute_adjoint .and. &
    (i .le. forcing_length*cadforcing_step) ) call lagrange_interp
    call step()
    end_time = MPI_WTIME()
    call cpu_time(t1)
    if (.not. displ) nor1 = norm(a(:,:,:,4)) !* dimc
    if (displ .and. test_in_2d) nor1 = sum(abs(a(:,:,:,6)**2.))**0.5

  ! Timing and 'convergence' statistics
  ! The criterion used to see if the results are 'convergent' is the
  ! L2 norm of the radial velocity.

   if (rank ==0) then

    print *, '----------------------------------------------------------'
    print *, '  Iteration #,  Cpu Time and Real Time:'
    print *,       time,        (t1 - T00),        (time+1)*timestep
    print *
    print *,'   Wall clock time', (end_time - start_time)
    print *
    print *,nor1
    print *, '----------------------------------------------------------'
    print *,''

    write(19,*) nor1

   endif
   timeline = timeline + timestep

   if (mod(time,steps) == 0) then

     write(28, *) time, timeline
     flush(28)

       if (compute_adjoint) tempstore = a(:,:,o_rad,6) !* rho0(:,:,o_rad)**0.5
       if (compute_forward .and. (.not. magnetic)) tempstore = a(:,:,o_rad,6)


      if (compute_forward .and. magnetic) then
        do k=1,dim2(rank)
         do j=1,nx
      tempstore(j,k) = a(j,k,orad_2d(j,k),6)
     enddo
        enddo
      endif

      indexglob = indexglob +1
      if ( compute_forward) then
        call write_binary_seq_2D(1244, tempstore, 1, indexglob)
        flush(1244)
      endif
       if (rank == 0) print *, 'Writing output'

       if (compute_forward) then
         if (.not. (COMPUTE_SYNTH .OR. LINESEARCH .OR. COMPUTE_DATA)) &
        call write_out_partial_state(directory//'forward_src'//contrib//'_ls'//jobno//'/')
     !call write_out_slice(directory//'forward'//contrib//'/')
       endif
       if (compute_adjoint) call write_out_partial_state(directory//'adjoint_src'//contrib//'/')

   endif
   if (mod(time+1,200) == 0)  then

     if (.not. kernel_mode) call write_out_full_state(directory)
     if (kernel_mode) then

       if (rank==0) open(99, file=directory//'unfinished_calc_'//contrib//'_'//jobno,action='write', status='replace')

       if (compute_adjoint) then
          call write_out_full_state(directory//'adjoint_src'//contrib//'/')

          if (rank==0) then
            if (.not. generate_wavefield) write(99,"(A)") 'adjoint'
            if (generate_wavefield) write(99,"(A)") 'spectral'
            write(99,*) time
            write(99,"(A)") contrib
            write(99,*) indexglob
            close(99)
          endif

       else
          call write_out_full_state(directory//'forward_src'//contrib//'_ls'//jobno//'/')
          if (rank==0) then
            write(99,"(A)") 'forward'
            write(99,*) time
            write(99,"(A)") contrib
            write(99,*) indexglob
            close(99)
          endif
       endif

     endif


     if (rank ==0 ) print *,'Full output done, index:',time

   endif
 end do

 if (rank==0) then

  call system('rm '//directory//'unfinished_calc_'//contrib//'_'//jobno)

!~   if (compute_forward) then
!~    open(223, file=directory//'status/'//'forward_src'//contrib//'_ls'//jobno,status='unknown')
!~    close(223)
   !call system('rm -rf '//directory//'forward'//contrib//'/*full*')
!~   endif

!~   if (compute_adjoint) then
!~    open(223, file=directory//'status/'//'adjoint_src'//contrib,status='unknown')
!~    close(223)
   !call system('rm -rf '//directory//'adjoint'//contrib//'/*full*')
!~   endif

 endif
END SUBROUTINE TIMESTEPPING


! --------------------------------------------------------------------------------------------

SUBROUTINE DETERMINE_STATUS(init, nsteps_given)

 use time_step
 use mp_physics
 use mp_physics_2d
 use initialize
 use all_modules

 implicit none

 integer i,j,k, init, nsteps_given,reclmax, indexnum
 real(kind=real64) start_time, end_time, T00, t1, nor1,  xloc, tempoa(nx,dim2(rank))
 real(kind=real64) sincx(nx), xdist(nx)
 character(len=80) calctype
 character(len=7) keyword
 character(len=1) ci
 integer inde
 logical lexist
 inquire(iolength=reclmax) tempoa
 init =0

!  ci = '1'
 ! if (contrib == '1') ci = '2'

  if (generate_wavefield) then
    COMPUTE_ADJOINT = .FALSE.
    COMPUTE_FORWARD = .TRUE.
    generate_wavefield = .true.
    keyword = 'forward_'
  endif

  if (compute_adjoint) then
    COMPUTE_ADJOINT = .TRUE.
    COMPUTE_FORWARD = .FALSE.
    generate_wavefield = .false.
    keyword = 'adjoint'
  endif
  initcond = .false.

  if (compute_forward) then
   allocate(vr(nx,dim2(rank),1), fwdsource(nsteps_given+8))
   vr = 0.0
   call forward_source(nsteps_given + 8)



   !if (rank==0) print *,'Need a file of: ', forcing_length, ' temporal slices.'

!   call read_binary_reverse(directory//'forward'//contrib//'/vz_2D_data.bin', &
!   vr(:,:,1:nsteps_given),nsteps_given)

   !forcing_length = nsteps_given
!    call readfits(directory//'forward'//contrib//'/source.fits',& !'_'//contrib//
!   vr(:,:,1:forcing_length),forcing_length)

!   if (contrib=='1') then
     read(contrib,'(I2)') indexnum

     open(356,file=directory//'master.pixels',action='read', position='rewind')
     do i=1,nmasters
      if (indexnum .ne. i) read(356,*)
      if (indexnum .eq. i) read(356,*) xloc
     enddo
!~       print *,"x_location",xloc
     close(356)

    dx = (x(2) - x(1)) * xlength*10.**(-8)
    xdist = (x-0.5)*xlength*10.**(-8) - xloc
    !sincx = sin(0.1*pi*xdist)/(0.1*pi*xdist)
    sincx = sin(0.43*pi*xdist)/(0.43*pi*xdist)
    sincx(minloc(abs(xdist))) = 1.0
    do j=1,dim2(rank)
      vr(:,j,1) = exp(-xdist**2./(400.*dx**2.))*sincx
    enddo

    dx = dx/(xlength*10.**(-8))
!    timeline = timeline + (nsteps_given -1)*cadforcing_step*floor(timestep)
!    timeline = -timeline

!    maxtime = floor((abs(timeline))*2./timestep)
!    nsteps_given = floor(maxtime/cadforcing_step*1.) + 1

    if (rank==0) open(28, file=directory//'forward_src'//contrib//'_ls'//jobno//'/timeline',&
        status='replace', action='write',position='append')


    open(1244,file=directory//'forward_src'//contrib//'_ls'//jobno//'/vz_cc.bin',&
        form='unformatted',status='replace', action ='write',&
    access='direct',recl=reclmax)!recordtype='fixed',


    delta_width = 2.0d0/(z(e_rad+1) - z(e_rad-1))

    if (magnetic) then
     do k=1,dim2(rank)
      do j=1,nx
       delta_width_2d(j,k) = 2.0/(z(erad_2d(j,k)+1) - z(erad_2d(j,k)-1))
      enddo
     enddo
    endif


   endif

   if (compute_adjoint) then

     !forcing_length = floor(maxtime/cadforcing_step*1.) + 1
    print *,forcing_length,maxtime, cadforcing_step
    allocate(vr(nx,dim2(rank),forcing_length + 8))

    vr = 0.0

    if (rank==0) print *,'Need a file of: ', forcing_length, ' temporal slices.'
    call readfits(directory//'adjoint_src'//contrib//'/source.fits',& !'_'//contrib//
        vr(:,:,1:forcing_length),forcing_length)

    if (rank==0 .and. (init .gt. 1)) then
    if (keyword=='forward_') then
        open(28, file=directory//'forward_src'//contrib//'_ls'//jobno//'/timeline',&
            status='old', action='write', position='append')
    elseif (keyword=='adjoint') then
        open(28, file=directory//'adjoint_src'//contrib//'/timeline',&
            status='old', action='write', position='append')
    endif
    endif

    if (rank==0 .and. (init .eq. 1)) then
    if (keyword=='forward_') then
        open(28, file=directory//'forward_src'//contrib//'_ls'//jobno//'/timeline',&
            status='replace', action='write', position='rewind')
    elseif (keyword=='adjoint') then
        open(28, file=directory//'adjoint_src'//contrib//'/timeline',&
            status='replace', action='write', position='rewind')
    endif
    endif

    delta_width = 2.0d0/(z(o_rad+1) - z(o_rad-1))

    if (magnetic) then
     do k=1,dim2(rank)
      do j=1,nx
       delta_width_2d(j,k) = 2.0/(z(orad_2d(j,k)+1) - z(orad_2d(j,k)-1))
      enddo
     enddo
    endif


   endif
 indexglob = 0
 time_old = 0
 init=0
END SUBROUTINE DETERMINE_STATUS

!================================================================================

 SUBROUTINE convert_binary_to_fits(input_filename, nxs, nys, dim3, output_filename)
  use all_modules
  implicit none
  integer nxs, nys, dim3, i, j, k, reclmax
  real(kind=real64) array(nxs, nys, dim3), tempoa(nxs,nys)
  character*(*) input_filename, output_filename
  inquire(iolength=reclmax) tempoa

  open(44,file=input_filename, access='direct',&!,recordtype='fixed'
        form='unformatted', action='read', &
        recl = reclmax, status='old')

   do k=1,dim3
    read(44, rec=k) ((array(i,j,k),i=1,nxs),j=1,nys)
   enddo
   close(44)

   call writefits_local(output_filename, array, nxs, nys, dim3)

 END SUBROUTINE CONVERT_BINARY_TO_FITS

!================================================================================================

SUBROUTINE ADJOINT_SOURCE_COMP(nt)
 use initialize
 use all_modules
 implicit none
 integer(kind=8) fwdplantemp, invplantemp, invplantemp2, fwdplandata, invplandata
 logical lexist
 integer i, nt, indexnum, pord, timesmax, halftime!, jj
 integer loc, leng, lef, rig, j, nmeasurements, nmeas, ierr
 character(len=1) ord
 real(kind=real64) mindist, maxdist, window, distances(nx), x00, t(nt), dat(nx, dim2(rank),nt), tau
 real(kind=real64) pcoef(5,4), adj(nx,dim2(rank),nt), acc(nx,dim2(rank),nt), windows(nt), filt(nt),dt
 real(kind=real64) freqnu(nt), leftcorner, rightcorner, ccdot(nx,dim2(rank),nt), dnu, con, misfit,misfit_tot
 complex(kind=real64), dimension(nx/2+1,dim2(rank),nt) :: filtout, temp, filtdat

 filtout = 1.0
 pcoef = 0.0
 misfit = 0.0
 read(contrib,*) indexnum

 if (rank==0) print *,'Doing Adjoint Source'

 open(356,file=directory//'master.pixels',action='read', position='rewind')
 do i=1,nmasters
  if (indexnum .ne. i) read(356,*)
  if (indexnum .eq. i) read(356,*) x00
 enddo
 close(356)

 distances = abs((x-0.5)*xlength*10.0**(-8.) - x00)

 open(44,file=directory//'forward_src'//contrib//'_ls'//jobno//'/timeline',action='read')
 do i=1,nt
  read(44,*) loc, t(nt-i+1)
 enddo
 close(44)
 t = -t/60.0
 dt = t(2)-t(1)

 call distmat(nt,1,freqnu)
 freqnu = freqnu/(nt*dt*60.)
 dnu = freqnu(2) - freqnu(1)

 call dfftw_plan_dft_r2c_3d(fwdplantemp, nx, dim2(rank), nt, acc, filtout, FFTW_ESTIMATE)
 call dfftw_plan_dft_c2r_3d(invplantemp, nx, dim2(rank), nt, filtout, acc, FFTW_ESTIMATE)
 call dfftw_plan_dft_c2r_3d(invplantemp2, nx, dim2(rank), nt, filtout, ccdot, FFTW_ESTIMATE)

 call dfftw_plan_dft_r2c_3d(fwdplandata, nx, dim2(rank), nt, dat, filtdat, FFTW_ESTIMATE)
 call dfftw_plan_dft_c2r_3d(invplandata, nx, dim2(rank), nt, filtdat, dat, FFTW_ESTIMATE)


 call readfits(directory//'forward_src'//contrib//'_ls'//jobno//'/vz_cc.fits', acc, nt)
 call readfits(directory//'data/'//contrib//'.fits', dat, nt)

 call dfftw_execute(fwdplantemp)
 call dfftw_execute(fwdplandata)

 call dfftw_destroy_plan(fwdplantemp)
 call dfftw_destroy_plan(fwdplandata)

 inquire(file=directory//'filter.params', exist=lexist)
 if (lexist) then

  open(94,file=directory//'filter.params', action='read',position='rewind')
   read(94,*) leftcorner ! mHz
   read(94,*) rightcorner  ! mHz
  close(94)

  filt = 1./(1. + exp((leftcorner - abs(freqnu))*2./dnu)) &
       - 1./(1. + exp((rightcorner - abs(freqnu))*2./dnu))

  do i=1,nt
   filtout(:,1,i) = filtout(:,1,i) * filt(i)
   filtdat(:,1,i) = filtdat(:,1,i) * filt(i)
  enddo
  temp = filtout
  call dfftw_execute(invplantemp)
  call dfftw_execute(invplandata)
  filtout = temp

  call dfftw_destroy_plan(invplantemp)
  call dfftw_destroy_plan(invplandata)

 endif
 con = 2.0*pi/(dble(nt)*dble(nx))
 do i=1,nt
  filtout(:,1,i) = filtout(:,1,i) * eye * freqnu(i) * con
 enddo

 call dfftw_execute(invplantemp2)
 call dfftw_destroy_plan(invplantemp2)


 pcoef = 0.0
 pcoef(1,1) = 11.5414
 pcoef(2,1) = 0.742286
 pcoef(3,1) = -0.00599439
 pcoef(4,1) = 2.59565e-5
 pcoef(5,1) = -2.69437e-8

 pcoef(1,2) = 16.9286
 pcoef(2,2) = 0.840476
 pcoef(3,2) = -0.00238095

 pcoef(1,3) = 25.0
 pcoef(2,3) = 0.9
 pcoef(3,3) = -0.002222

 adj = 0.0
 nmeasurements = 0
 do pord=1,4
  call convert_to_string(pord, ord, 1)

  inquire(file=directory//'params.p'//ord, exist=lexist)
  if (lexist) then
    open(97, file = directory//'params.p'//ord, action = 'read', status='old')
    read(97,*) mindist
    read(97,*) maxdist
    read(97,*) window
    close(97)
!    open(993, file=directory//'forward'//contrib//'/times.p'//ord,status='unknown',position='rewind',action='write')
!    open(45, file='temp',status='unknown',position='rewind',action='write')
    halftime = nint(window/(2.*dt))
    leng = 2*halftime+1
    do i=1,nx

     if ((distances(i) > mindist) .and. (distances(i) < maxdist)) then

      timesmax = nint((pcoef(1,pord) - t(1) + pcoef(2,pord)*distances(i) + pcoef(3,pord)*distances(i)**2. &
     + pcoef(4,pord)*distances(i)**3. + pcoef(5,pord)*distances(i)**4.)/dt) + 1

      loc = maxloc(acc(i,1,(timesmax-6):(timesmax+6)),1) + timesmax - 7
      lef = loc - halftime
      rig = loc + halftime
      print *, '1'
      print *, 'lef', lef
      print *, 'rig', rig
      print *, 'loc', loc
      print *, 'halftime', halftime
      call compute_tt(acc(i,1,lef:rig),dat(i,1,lef:rig),tau,dt,leng)
!        open(45,file='temp',action='write')
!       do jj=1,nt !lef,rig
!       write(45,*) acc(i,1,jj), dat(i,1,jj)
!       enddo
!       close(45)
!   print *,lef, rig,tau, i
!   stop
!      write(993, *) (x(i)-0.5)*xlength*10.**(-8.),tau
      misfit = misfit + tau**2.
      nmeasurements = nmeasurements + 1
      con = -tau/(sum(ccdot(i,1,lef:rig)**2.)*dt)
      do j=lef,rig
       adj(i,1,nt-j+1) = ccdot(i,1,j) * con + adj(i,1,nt-j+1)
      enddo
     endif
    enddo
 !   close(993)
  endif
 enddo
 adj = adj * nx
 inquire(file=directory//'adjoint_src'//contrib, exist=lexist)
 if (.not. lexist .and. (rank==0)) &
    call system('mkdir '//directory//'adjoint_src'//contrib)

 call writefits_3d(directory//'adjoint_src'//contrib//'/source.fits',adj,nt)

 call MPI_REDUCE(misfit, misfit_tot, 1, MPI_DOUBLE_PRECISION, &
                MPI_SUM, 0, MPI_COMM_WORLD, ierr)

 call MPI_REDUCE(nmeasurements, nmeas, 1, MPI_INTEGER, &
                MPI_SUM, 0, MPI_COMM_WORLD, ierr)

! misfit_tot = misfit_tot * 0.5
! print *,'Total',misfit_tot, misfit_tot/nmeasurements
 if (rank==0) then
   open(34, file=directory//'adjoint_src'//contrib//'/information', &
        action='write', position='rewind', status='replace')
   write(34,*) nt
   write(34,*) t(1)*60.
   write(34,*) t(nt)*60.
   if (test_in_2d) write(34,*) '2D'
   close(34)

    open(88, file=directory//'kernel/misfit_'//contrib//'_'//jobno, status='unknown',&
               action='write')


   misfit_tot = misfit_tot * 0.5
   write(88,*) indexnum, nmeas, misfit_tot
   close(88)
   flush(88)

 endif

END SUBROUTINE ADJOINT_SOURCE_COMP

!================================================================================

SUBROUTINE COMPUTE_TT(u0, u, tau, dt, nt)

 implicit none
 integer, parameter :: real64 = kind(1.d0)

 integer nt,l1,l2, loc, i
 real(kind=real64) u0(nt), u(nt), cc(-nt+1:nt-1), t(-nt+1:nt-1)
 real(kind=real64) times(3), tau, invmat(3,3), p1, p2, dt, mat(3,3)
 do i=(-nt+1),0
  t(i) = i*dt
  l1 = -i + 1
  l2 = nt
  cc(i) = sum(u0(l1:l2)*u(1:(l2-l1+1))) * dt
 enddo

 do i=1,nt-1
  t(i) = i*dt
  l1 = i + 1
  l2 = nt-i
  cc(i) = sum(u0(1:l2)*u(l1:nt)) * dt
 enddo

 loc = maxloc(cc,1)-nt
 times(1) = t(loc-1)
 times(2) = t(loc)
 times(3) = t(loc+1)
 mat(:,1) = 1
 mat(:,2) = times
 mat(:,3) = times**2.
 call inverse(mat, invmat, 3)
 p1 = invmat(2,1)*cc(loc-1) + invmat(2,2) * cc(loc) + invmat(2,3) * cc(loc+1)
 p2 = invmat(3,1)*cc(loc-1) + invmat(3,2) * cc(loc) + invmat(3,3) * cc(loc+1)

 tau = -p1*0.5/p2
 !print *,loc,-nt+1, nt-1, tau
END SUBROUTINE COMPUTE_TT

!================================================================================

SUBROUTINE COMPUTE_TT_GIZONBIRCH(u0,u,tau,dt,nt, lef, rig)
    use initialize
    use integrals
    implicit none
    integer, intent(in) :: nt, lef, rig
    real(kind=real64), intent(in) :: u0(nt), u(nt)
    real(kind=real64) dt
    real(kind=real64), intent(inout) :: tau
    real(kind=real64) u0dot(nt), window(nt)
    integer k
    integer(kind=8) plan
    complex(kind=real64) u0w(nt/2+1)
    real(kind=real64) numerator,denominator


    window = 0.0
    window(lef:rig) = 1.0
    call dfftw_plan_dft_r2c_1d(plan,nt,u0,u0w,FFTW_ESTIMATE)
    call dfftw_execute_dft_r2c(plan, u0, u0w)
    call dfftw_destroy_plan(plan)

    do k=1,nt/2
        u0w(k)=u0w(k)*(eye*2*pi*(k-1.0)/(nt*dt))
    end do

    u0w(nt/2+1)=0

    call dfftw_plan_dft_c2r_1d(plan,nt,u0w,u0dot,FFTW_ESTIMATE)
    call dfftw_execute_dft_c2r(plan, u0w, u0dot)
    call dfftw_destroy_plan(plan)

    u0dot=u0dot/nt

!~     open(unit=395,file="u_udot",status='REPLACE')
!~     do k=1,nt
!~         write(395,'(F14.3,X,F14.3,X,F14.3)')  u0(k),u0dot(k),u(k)
!~     end do
!~     close(395)

    call integrate_time(window*u0dot*(u0-u),numerator,dt)
    call integrate_time(window*u0dot**2,denominator,dt)

    tau=numerator/denominator
    print *,tau*60.,dt*60.


END SUBROUTINE COMPUTE_TT_GIZONBIRCH

!================================================================================

SUBROUTINE VZ_FROM_VX_CONTINUITY(vx,vz)

    use initialize
    use derivatives
    use integrals
    use all_modules
    implicit none

    real(kind=real64), dimension(nx,1,nz), intent(in) :: vx
    real(kind=real64), dimension(nx,1,nz) :: dxrhovx
    real(kind=real64), dimension(nx,1,nz), intent(out) :: vz
!~     CHARACTER(len=255) :: cwd,homedir,pythoncmd
!~     integer e
!~     logical rhoex


    call ddx(rho0*vx,dxrhovx,1)
    call writefits_3d("dxrhovx.fits",dxrhovx,nz)
    call integrate_z(dxrhovx,vz)
    vz=-vz/rho0

!~     call writefits_3d("dxrhovx_"//contrib//"_"//jobno//".fits",dxrhovx,nz)
    !stop

!~     inquire(file="rho.fits",exist=rhoex)
!~     if (rhoex .eqv. .false.) then
!~     call writefits_3d("rho.fits",rho0,nz)
!~     endif
!~     CALL getcwd(cwd)
!~     CALL getenv("HOME", homedir)
!~     pythoncmd = trim(trim(homedir)//"/anaconda/bin/python "//trim(cwd)//&
!~     "/vz_from_vx_continuity.py "//contrib//" "//jobno)
!~     call system(pythoncmd,e)
!~     if (e /= 0) print *,"Python call exit code",e
!~     call readfits('vz_int_'//contrib//'_'//jobno//'.fits',vz,nz)
!~     call system('rm vz_int_'//contrib//'_'//jobno//'.fits')
!~     call system("rm dxrhovx_"//contrib//"_"//jobno//".fits")


END SUBROUTINE VZ_FROM_VX_CONTINUITY

!================================================================================

SUBROUTINE VX_FROM_VZ_CONTINUITY(vz,vx)

    use initialize
    use derivatives
    use integrals
    use all_modules
    implicit none

    real(kind=real64), dimension(nx,1,nz), intent(in) :: vz
    real(kind=real64), dimension(nx,1,nz) :: dzrhovz
    real(kind=real64), dimension(nx,1,nz), intent(out) :: vx

    call ddz(rho0*vz,dzrhovz,1)
!~     call writefits_3d("dzrhovz.fits",dzrhovz,nz)
    call integrate_x(dzrhovz,vx)
    vx=-vx/rho0

END SUBROUTINE

!================================================================================

SUBROUTINE CONTINUITY_CHECK(vx,vz)
    use initialize
    use derivatives
    use all_modules
    implicit none
    real(kind=real64), intent(in),dimension(nx,1,nz) :: vx,vz
    real(kind=real64), dimension(nx,1,nz) :: cont
    real(kind=real64), dimension(:,:,:), allocatable :: dxrhovx,dzrhovz,dzrho

    allocate(dxrhovx(nx,1,nz),dzrhovz(nx,1,nz),dzrho(nx,1,nz))

    call ddz(rho0*vz,dzrhovz,1)
    call ddx(rho0*vx,dxrhovx,1)
    call ddz(rho0,dzrho,1)


    cont=(dxrhovx + dzrhovz)/c_speed/dzrho

    deallocate(dzrhovz,dxrhovx,dzrho)

    print *,"Continuity check maxval",maxval(abs(cont)),", should ideally be zero"

!~     if (contrib=="01") call writefits_3d('contcheck_ls'//jobno//'.fits',cont,nz)

END SUBROUTINE CONTINUITY_CHECK


!================================================================================


  subroutine inverse(a,c,n)
!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed
! during the calculation
!===========================================================
implicit none
integer n
double precision a(n,n), c(n,n)
double precision L(n,n), U(n,n), b(n), d(n), x(n)
double precision coeff
integer i, j, k

! step 0: initialization for matrices L and U and b
! Fortran 90/95 aloows such operations on matrices
L=0.0
U=0.0
b=0.0

! step 1: forward elimination
do k=1, n-1
   do i=k+1,n
      coeff=a(i,k)/a(k,k)
      L(i,k) = coeff
      do j=k+1,n
         a(i,j) = a(i,j)-coeff*a(k,j)
      end do
   end do
end do

! Step 2: prepare L and U matrices
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
do i=1,n
  L(i,i) = 1.0
end do
! U matrix is the upper triangular part of A
do j=1,n
  do i=1,j
    U(i,j) = a(i,j)
  end do
end do

! Step 3: compute columns of the inverse matrix C
do k=1,n
  b(k)=1.0
  d(1) = b(1)
! Step 3a: Solve Ld=b using the forward substitution
  do i=2,n
    d(i)=b(i)
    do j=1,i-1
      d(i) = d(i) - L(i,j)*d(j)
    end do
  end do
! Step 3b: Solve Ux=d using the back substitution
  x(n)=d(n)/U(n,n)
  do i = n-1,1,-1
    x(i) = d(i)
    do j=n,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
    end do
    x(i) = x(i)/u(i,i)
  end do
! Step 3c: fill the solutions x(n) into column k of C
  do i=1,n
    c(i,k) = x(i)
  end do
  b(k)=0.0
end do
end subroutine inverse

!================================================================================

SUBROUTINE FORWARD_SOURCE(nt)
 use initialize
! use all_modules
 implicit none
 integer nt, k
 real(kind=real64) ttime!,nupea,nuwid

!nupeak = 0.0025
! nuwidth=0.002

!nupea = 0.0017
!nuwid=0.0055

 do k=1,nt
  ttime = (k- 1.)*timestep
  fwdsource(k) = cos(2.*pi*nupeak*ttime) * exp(-2.*(ttime-0.5/nupeak)**2.*nuwidth**2.)
  if (abs(fwdsource(k)) .lt. 0.00001) fwdsource(k) = 0.0
 enddo
 fwdsource = fwdsource * 10.**(-5.)
 open(2535,file='source.txt',action='write',status='unknown')
 do k=1,nt
   write(2535,*) fwdsource(k)
 enddo
 close(2535)

END SUBROUTINE FORWARD_SOURCE
!================================================================================

SUBROUTINE FMODE_FILTER(nt, fmode)
  use initialize
  use all_modules
  implicit none

  integer i,j,nt
  real(kind=real64) fmode(nx, dim2(rank), nt),f_low,df,k(nx),dt,f_mode_const
  real(kind=real64) Poly(0:2),Polylow(0:2), f0(nx),w(nt),f(nt),f1(nx),d,delta

  dt = outputcad

  f_mode_const=2.056565

  Poly(0)=1.1
    Poly(1)=1.9
    Poly(2)=-0.2

    Polylow(0)=0.7
    Polylow(1)=1.7
    Polylow(2)=-0.2

  f_low = 1.1
  df = 0.5

  call distmat(nx,1,k)
  call distmat(nt,1,w)
  k = abs(k) * 2.*pi/(xlength*10.**(-8.)*nx/(nx-1.))
  w = abs(w) * 2.*pi/(nt*dt)

  f0=Polylow(0) + Polylow(1)*k +Polylow(2)*k**2.
    f1=Poly(0) + Poly(1)*k +Poly(2)*k**2.
  f = w/(2.*pi)*1e3


  fmode = 0.0
  do i=1,nx
   delta = (f1(i) - f0(i))
    do j=1,nt
     d = f(j) - f0(i)
     if ((d .lt. delta) .and. (d>0)) &
        fmode(i,1,j) = 0.5*(1.+cos(pi*(d-delta/2.0)/(delta/2.0)))
    enddo
   enddo

   do j=1,nt
    if (f(j) .lt. f_low+df) &
      fmode(:,1,j) = fmode(:,1,j)*0.5*(1.+cos(pi*(f(j)-(f_low+df))/df) )
    if (f(j) .lt. f_low) fmode(:,1,j) = 0.
   enddo

END SUBROUTINE FMODE_FILTER
!================================================================================


SUBROUTINE HIGHPMODE_FILTER(nt, pmode)

  use initialize
  use all_modules
  implicit none
  integer i,j,nt
  real(kind=real64) pmode(nx, dim2(rank), nt),f_low,df,k(nx),dt
  real(kind=real64) Poly(0:2), f0(nx),w(nt),f(nt),f1(nx),d,delta

  dt = outputcad

  Poly(0)=0.0080513021
  Poly(1)=0.02369524
  Poly(2)=-0.00460597

  df = 0.5
  f_low = 1.6

  call distmat(nx,1,k)
  call distmat(nt,1,w)
  k = abs(k) * 2.*pi/(xlength*10.**(-8.)*nx/(nx-1.))
  w = abs(w) * 2.*pi/(nt*dt)

  f0=2.0*(254e-6*abs(k))**0.5/(2*pi)*1e3
  f1=1.85*(Poly(0) + Poly(1)*abs(k) +Poly(2)*abs(k)**2.)/(2*pi)*1e3
  f = w/(2.*pi)*1e3


  pmode = 0.0
  do i=1,nx
   delta = (f1(i) - f0(i))*0.5
    do j=1,nt
     d = f(j) - f0(i)
     if (abs(d) .gt. delta) pmode(i,1,j) = 0.
     if (abs(d) .lt. delta) pmode(i,1,j) = 1.
     if ((abs(d) .lt. delta) .and. (abs(d) .gt. delta*0.5)) &
        pmode(i,1,j) = 0.5*(1.+cos(pi*(abs(d)-abs(delta)*0.5)/(abs(delta)*0.5)))
    enddo
   enddo

   do j=1,nt
    if (f(j) .lt. f_low) pmode(:,1,j) = 0.
    if (f(j) .lt. f_low+df) &
      pmode(:,1,j) = pmode(:,1,j)*0.5*(1.+cos(pi*(f(j)-(f_low+df))/df) )
   enddo

END SUBROUTINE HIGHPMODE_FILTER
!================================================================================

SUBROUTINE P1MODE_FILTER(nt, pmode)

  use initialize
  use all_modules
  implicit none
  integer i,j,nt,nrow
  real(kind=real64) pmode(nx, dim2(rank), nt),f_low,df,k(nx),dt,f_mode_const
  real(kind=real64) Poly(0:2), f0(nx),w(nt),f(nt),f1(nx),d,delta,Polylow(0:2)

  open (unit=32,file="filter.txt",action="write",status="replace")
  dt = outputcad

!~   f_mode_const=3.096

  Poly(0)=1.4
    Poly(1)=3.0
    Poly(2)=-0.5

    Polylow(0)=1.1
    Polylow(1)=2.4
    Polylow(2)=-0.3

  f_low = 1.6
  df = 0.5

  call distmat(nx,1,k)
  call distmat(nt,1,w)
  k = abs(k) * 2.*pi/(xlength*10.**(-8.)*nx/(nx-1.))
  w = abs(w) * 2.*pi/(nt*dt)

!~   f0=f_mode_const*abs(k)**0.5
  f0=Polylow(0) + Polylow(1)*k +Polylow(2)*k**2.
  f1=Poly(0) + Poly(1)*k +Poly(2)*k**2.
  f = w/(2.*pi)*1e3

  pmode = 0.0
  do i=1,nx
   delta = (f1(i) - f0(i))
    do j=1,nt
     d = f(j) - f0(i)
     if ((d .lt. delta) .and. (d .gt. 0)) then
        pmode(i,1,j) = 0.5*(1.+cos(pi*(abs(d)-abs(delta)*0.5)/(abs(delta)*0.5)))
     end if
    enddo
   enddo

   do j=1,nt
    if (f(j) .lt. f_low) pmode(:,1,j) = 0.
    if (f(j) .lt. f_low+df) &
      pmode(:,1,j) = pmode(:,1,j)*0.5*(1.+cos(pi*(f(j)-(f_low+df))/df) )
   enddo

   close(32)

END SUBROUTINE P1MODE_FILTER
!================================================================================

SUBROUTINE P2MODE_FILTER(nt, pmode)

  use initialize
  use all_modules
  implicit none
  integer i,j,nt,nrow
  real(kind=real64) pmode(nx, dim2(rank), nt),f_low,df,k(nx),dt,f_mode_const
  real(kind=real64) Poly(0:2),Polylow(0:2), f0(nx),w(nt),f(nt),f1(nx),d,delta

  open (unit=32,file="filter.txt",action="write",status="replace")
  dt = outputcad

  f_mode_const=4.1004

  Poly(0)=1.6
    Poly(1)=3.8
    Poly(2)=-0.65

    Polylow(0)=1.4
    Polylow(1)=3.3
    Polylow(2)=-0.62

  f_low = 1.6
  df = 0.5

  call distmat(nx,1,k)
  call distmat(nt,1,w)
  k = abs(k) * 2.*pi/(xlength*10.**(-8.)*nx/(nx-1.))
  w = abs(w) * 2.*pi/(nt*dt)

  f0=f_mode_const*abs(k)**0.5
  f1=Poly(0) + Poly(1)*k +Poly(2)*k**2.
  f = w/(2.*pi)*1e3

  pmode = 0.0
  do i=1,nx
   delta = (f1(i) - f0(i))
    do j=1,nt
     d = f(j) - f0(i)
     if ((d .lt. delta) .and. (d .gt. 0)) then
        pmode(i,1,j) = 0.5*(1.+cos(pi*(abs(d)-abs(delta)*0.5)/(abs(delta)*0.5)))
     end if
    enddo
   enddo

   do j=1,nt
    if (f(j) .lt. f_low) pmode(:,1,j) = 0.
    if (f(j) .lt. f_low+df) &
      pmode(:,1,j) = pmode(:,1,j)*0.5*(1.+cos(pi*(f(j)-(f_low+df))/df) )
   enddo

   close(32)

END SUBROUTINE P2MODE_FILTER
!================================================================================

SUBROUTINE P3MODE_FILTER(nt, pmode)

  use initialize
  use all_modules
  implicit none
  integer i,j,nt,nrow
  real(kind=real64) pmode(nx, dim2(rank), nt),f_low,df,k(nx),dt!f_mode_const
  real(kind=real64) Poly(0:2),Polylow(0:2), f0(nx),w(nt),f(nt),f1(nx),d,delta

  open (unit=32,file="filter.txt",action="write",status="replace")
  dt = outputcad

!~   f_mode_const=4.7

    Poly(0)=2
    Poly(1)=4.1
    Poly(2)=-0.8

    Polylow(0)=2
    Polylow(1)=3.55
    Polylow(2)=-0.7

    f_low = 1.6
    df = 0.5

  call distmat(nx,1,k)
  call distmat(nt,1,w)
  k = abs(k) * 2.*pi/(xlength*10.**(-8.)*nx/(nx-1.))
  w = abs(w) * 2.*pi/(nt*dt)

  f0=Polylow(0) + Polylow(1)*k +Polylow(2)*k**2.
    f1=Poly(0) + Poly(1)*k +Poly(2)*k**2.
  f = w/(2.*pi)*1e3

  pmode = 0.0
  do i=1,nx
   delta = (f1(i) - f0(i))
    do j=1,nt
     d = f(j) - f0(i)
     if ((d .lt. delta) .and. (d .gt. 0)) then
        pmode(i,1,j) = 0.5*(1.+cos(pi*(abs(d)-abs(delta)*0.5)/(abs(delta)*0.5)))
     end if
    enddo
   enddo

   do j=1,nt
    if (f(j) .lt. f_low) pmode(:,1,j) = 0.
    if (f(j) .lt. f_low+df) &
      pmode(:,1,j) = pmode(:,1,j)*0.5*(1.+cos(pi*(f(j)-(f_low+df))/df) )
   enddo

   close(32)

END SUBROUTINE P3MODE_FILTER
!================================================================================

SUBROUTINE P4MODE_FILTER(nt, pmode)

  use initialize
  use all_modules
  implicit none
  integer i,j,nt,nrow
  real(kind=real64) pmode(nx, dim2(rank), nt),f_low,df,k(nx),dt,f_mode_const
  real(kind=real64) Poly(0:2),Polylow(0:2), f0(nx),w(nt),f(nt),f1(nx),d,delta

  open (unit=32,file="filter.txt",action="write",status="replace")
  dt = outputcad

  f_mode_const=5.4

    Poly(0)=2.3
    Poly(1)=4.4
    Poly(2)=-0.7

    Polylow(0)=2.15
    Polylow(1)=4
    Polylow(2)=-0.7

    f_low = 1.6
    df = 0.5

  call distmat(nx,1,k)
  call distmat(nt,1,w)
  k = abs(k) * 2.*pi/(xlength*10.**(-8.)*nx/(nx-1.))
  w = abs(w) * 2.*pi/(nt*dt)

  f0=Polylow(0) + Polylow(1)*k +Polylow(2)*k**2.
    f1=Poly(0) + Poly(1)*k +Poly(2)*k**2.
  f = w/(2.*pi)*1e3

  pmode = 0.0
  do i=1,nx
   delta = (f1(i) - f0(i))
    do j=1,nt
     d = f(j) - f0(i)
     if ((d .lt. delta) .and. (d .gt. 0)) then
        pmode(i,1,j) = 0.5*(1.+cos(pi*(abs(d)-abs(delta)*0.5)/(abs(delta)*0.5)))
     end if
    enddo
   enddo

   do j=1,nt
    if (f(j) .lt. f_low) pmode(:,1,j) = 0.
    if (f(j) .lt. f_low+df) &
      pmode(:,1,j) = pmode(:,1,j)*0.5*(1.+cos(pi*(f(j)-(f_low+df))/df) )
   enddo

   close(32)

END SUBROUTINE P4MODE_FILTER
!================================================================================

SUBROUTINE P5MODE_FILTER(nt, pmode)

  use initialize
  use all_modules
  implicit none
  integer i,j,nt,nrow
  real(kind=real64) pmode(nx, dim2(rank), nt),f_low,df,k(nx),dt,f_mode_const
  real(kind=real64) Poly(0:2),Polylow(0:2), f0(nx),w(nt),f(nt),f1(nx),d,delta

  open (unit=32,file="filter.txt",action="write",status="replace")
  dt = outputcad

!~   f_mode_const=4.1004

!~   f_mode_const=5.8

    Poly(0)=2.35
    Poly(1)=5.6
    Poly(2)=-1.1

    Polylow(0)=2.2
    Polylow(1)=4.7
    Polylow(2)=-1.0

    f_low = 1.6
    df = 0.5

  call distmat(nx,1,k)
  call distmat(nt,1,w)
  k = abs(k) * 2.*pi/(xlength*10.**(-8.)*nx/(nx-1.))
  w = abs(w) * 2.*pi/(nt*dt)

!~   f0=f_mode_const*abs(k)**0.5
  f0=Polylow(0) + Polylow(1)*k +Polylow(2)*k**2.
  f1=Poly(0) + Poly(1)*k +Poly(2)*k**2.
  f = w/(2.*pi)*1e3

  pmode = 0.0
  do i=1,nx
   delta = (f1(i) - f0(i))
    do j=1,nt
     d = f(j) - f0(i)
     if ((d .lt. delta) .and. (d .gt. 0)) then
        pmode(i,1,j) = 0.5*(1.+cos(pi*(abs(d)-abs(delta)*0.5)/(abs(delta)*0.5)))
     end if
    enddo
   enddo

   do j=1,nt
    if (f(j) .lt. f_low) pmode(:,1,j) = 0.
    if (f(j) .lt. f_low+df) &
      pmode(:,1,j) = pmode(:,1,j)*0.5*(1.+cos(pi*(f(j)-(f_low+df))/df) )
   enddo

   close(32)

END SUBROUTINE P5MODE_FILTER
!================================================================================

SUBROUTINE P6MODE_FILTER(nt, pmode)

  use initialize
  use all_modules
  implicit none
  integer i,j,nt,nrow
  real(kind=real64) pmode(nx, dim2(rank), nt),f_low,df,k(nx),dt,f_mode_const
  real(kind=real64) Poly(0:2),Polylow(0:2), f0(nx),w(nt),f(nt),f1(nx),d,delta

  open (unit=32,file="filter.txt",action="write",status="replace")
  dt = outputcad

!~   f_mode_const=4.1004

!~   f_mode_const=5.8

!~     Poly(0)=2.65
!~     Poly(1)=5.7
!~     Poly(2)=-1.2
!~
!~     Polylow(0)=2.4
!~     Polylow(1)=5.4
!~     Polylow(2)=-1.3
    Poly(0)=1.98567115899
    Poly(1)=8.09108986838
    Poly(2)=-3.20316331815

    Polylow(0)=1.80035417224
    Polylow(1)=7.42939105658
    Polylow(2)=-2.84595764385

    f_low = 1.6
    df = 0.5

  call distmat(nx,1,k)
  call distmat(nt,1,w)
  k = abs(k) * 2.*pi/(xlength*10.**(-8.)*nx/(nx-1.))
  w = abs(w) * 2.*pi/(nt*dt)

!~   f0=f_mode_const*abs(k)**0.5
  f0=Polylow(0) + Polylow(1)*k +Polylow(2)*k**2.
  f1=Poly(0) + Poly(1)*k +Poly(2)*k**2.
  f = w/(2.*pi)*1e3

  pmode = 0.0
  do i=1,nx
   delta = (f1(i) - f0(i))
    do j=1,nt
     d = f(j) - f0(i)
     if ((d .lt. delta) .and. (d .gt. 0)) then
        pmode(i,1,j) = 0.5*(1.+cos(pi*(abs(d)-abs(delta)*0.5)/(abs(delta)*0.5)))
     end if
    enddo
   enddo

   do j=1,nt
    if (f(j) .lt. f_low) pmode(:,1,j) = 0.
    if (f(j) .lt. f_low+df) &
      pmode(:,1,j) = pmode(:,1,j)*0.5*(1.+cos(pi*(f(j)-(f_low+df))/df) )
   enddo

   close(32)

END SUBROUTINE P6MODE_FILTER
!================================================================================

SUBROUTINE P7MODE_FILTER(nt, pmode)

  use initialize
  use all_modules
  implicit none
  integer i,j,nt,nrow
  real(kind=real64) pmode(nx, dim2(rank), nt),f_low,df,k(nx),dt,f_mode_const
  real(kind=real64) Poly(0:2),Polylow(0:2), f0(nx),w(nt),f(nt),f1(nx),d,delta

  open (unit=32,file="filter.txt",action="write",status="replace")
  dt = outputcad

!~   f_mode_const=4.1004

!~   f_mode_const=5.8

!~     Poly(0)=2.9
!~     Poly(1)=5.9
!~     Poly(2)=-1.3
!~
!~     Polylow(0)=2.65
!~     Polylow(1)=5.7
!~     Polylow(2)=-1.2
    Poly(0)=2.18544600032
    Poly(1)=8.68183289647
    Poly(2)=-3.84478880142

    Polylow(0)=1.98297334673
    Polylow(1)=7.94931076885
    Polylow(2)=-3.00725356897

    f_low = 1.6
    df = 0.5

  call distmat(nx,1,k)
  call distmat(nt,1,w)
  k = abs(k) * 2.*pi/(xlength*10.**(-8.)*nx/(nx-1.))
  w = abs(w) * 2.*pi/(nt*dt)

!~   f0=f_mode_const*abs(k)**0.5
  f0=Polylow(0) + Polylow(1)*k +Polylow(2)*k**2.
  f1=Poly(0) + Poly(1)*k +Poly(2)*k**2.
  f = w/(2.*pi)*1e3

  pmode = 0.0
  do i=1,nx
   delta = (f1(i) - f0(i))
    do j=1,nt
     d = f(j) - f0(i)
     if ((d .lt. delta) .and. (d .gt. 0)) then
        pmode(i,1,j) = 0.5*(1.+cos(pi*(abs(d)-abs(delta)*0.5)/(abs(delta)*0.5)))
     end if
    enddo
   enddo

   do j=1,nt
    if (f(j) .lt. f_low) pmode(:,1,j) = 0.
    if (f(j) .lt. f_low+df) &
      pmode(:,1,j) = pmode(:,1,j)*0.5*(1.+cos(pi*(f(j)-(f_low+df))/df) )
   enddo

   close(32)

END SUBROUTINE P7MODE_FILTER
!================================================================================

SUBROUTINE ALL_PMODE_FILTER(nt, pmode)

  use initialize
  use all_modules
  implicit none
  integer i,j,nt,nrow
  real(kind=real64) pmode(nx, dim2(rank), nt),f_low,df,k(nx),dt,f_mode_const
  real(kind=real64) Poly(0:2),Polylow(0:2), f0(nx),w(nt),f(nt),f1(nx),d,delta

  open (unit=32,file="filter.txt",action="write",status="replace")
  dt = outputcad

    Polylow(0)=1.1
    Polylow(1)=2.4
    Polylow(2)=-0.3

    Poly(0)=3.5
    Poly(1)=6.5
    Poly(2)=-1.3

    f_low = 1.6
    df = 0.5

  call distmat(nx,1,k)
  call distmat(nt,1,w)
  k = abs(k) * 2.*pi/(xlength*10.**(-8.)*nx/(nx-1.))
  w = abs(w) * 2.*pi/(nt*dt)

  f0=Polylow(0) + Polylow(1)*k +Polylow(2)*k**2.
  f1=Poly(0) + Poly(1)*k +Poly(2)*k**2.
  f = w/(2.*pi)*1e3

  pmode = 1.0
  do i=1,nx
   delta = (f1(i) - f0(i))
    do j=1,nt
     d = f(j) - f0(i)
     if ((d .lt. delta) .and. (d .gt. 0)) then
        pmode(i,1,j) = 1
     end if
    enddo
   enddo

   do j=1,nt
    if (f(j) .lt. f_low) pmode(:,1,j) = 0.
    if (f(j) .lt. f_low+df) &
      pmode(:,1,j) = pmode(:,1,j)*0.5*(1.+cos(pi*(f(j)-(f_low+df))/df) )
   enddo

   close(32)

END SUBROUTINE ALL_PMODE_FILTER
!================================================================================

SUBROUTINE FREQ_FILTER(f1, f2, nt, filt)
  use initialize
  implicit none
  integer nt, i
  real(kind=real64) f1,f2,filt(nt),w(nt),wid,dt
  ! f1, f2 are in mHz
  call distmat(nt,1,w)
  dt =outputcad
  w = abs(w)*1e3/(nt*dt)
  wid = 0.4*(w(2)-w(1))
  filt = 1./(1.+exp((f1-w)/wid)) - 1./(1. + exp((f2-w)/wid))

END SUBROUTINE FREQ_FILTER

!================================================================================

SUBROUTINE PHASE_FILTER(speed, var, nt, filt)
  use initialize
  implicit none
  integer nt, k
  real(kind=real64) speed, var,filt(nx,dim2(rank),nt),w(nt),kay(nx),dt,dw
  ! f1, f2 are in mHz
  call distmat(nt,1,w)
  dt =outputcad
  w = abs(w)*2.*pi/(nt*dt)
  dw = w(2) - w(1)
  call distmat(nx,1,kay)
  kay = abs(kay)*2.*pi/(xlength*10.**(-8.))
  do k=1,nt
   filt(:,1,k) = exp(-(w(k)*1e3/(kay+0.01) - speed)**2./(2.*var**2.)) * &
                1./(1. + exp((0.002*2.*pi - w(k))/(dw*0.5)))
  enddo

END SUBROUTINE PHASE_FILTER

!================================================================================

SUBROUTINE ADJOINT_SOURCE_FILT(nt)
    use initialize
    use all_modules
    use bspline
    implicit none
    integer(kind=8) fwdplantemp, invplantemp, invplantemp2!, fwdquiet, invquiet
    integer(kind=8)  invplantemp3, fwdplandata, invplandata, onedplan
    logical :: lexist,exist_win,file_open
    integer i, nt, indexnum, pord, timesmax, halftime, bounce,timefin, idiff,dumm(0:5)
    integer loc, leng, lef, rig, j, nmeasurements, nmeas, ierr,timest, filenum
    character(len=1) ord
    character(len=2) contribmatch
    real(kind=real64) mindist, maxdist, window, distances(nx), x00, t(nt),  tau, signed(nx),taus(nx)
    real(kind=real64) ampquiet, ampmag
    real(kind=real64),dimension(nx,dim2(rank),nt):: p1mode,p2mode,fmode,all_else,filter,phase,temparr,&
                                      all_pmode,p3mode,p4mode,p5mode,p6mode,p7mode
    real(kind=real64) pcoef(5,4), adj(nx,dim2(rank),nt), windows(nt), filt(nt),dt,xdim(nx), vel(0:10)
    real(kind=real64) freqnu(nt), leftcorner, rightcorner, dnu, con, misfit,misfit_tot
    complex(kind=real64), dimension(nx,dim2(rank),nt) :: filtout, tempout, tempdat,& !, filtquiet, tempquiet, &
                                 filtdat, filtex, filtemp,dat,acc,ccdot!, quiet, quietfilt
    complex(kind=real64), dimension(nx) :: eyekh
    complex(kind=real64), dimension(nt) :: oned, ccdotone
    complex(kind=real64) UNKNOWN
    integer kxord
    parameter(kxord=3)
    real(kind=real64) speed, var, param(6), offset, bcoef(nx),xknot(kxord+nx)
    real(kind=real64) ign1,ign2, dist


    UNKNOWN = 1.0/0.55
    filtout = 1.0
    pcoef = 0.0
    misfit = 0.0
    call distmat(nx,1,distances)

    dx = x(2)-x(1)
    eyekh= cmplx(0,1)*distances*2.*pi
    ! eyekh(nx/2+1) = 0.0
    read(contrib,*) indexnum

    if (rank==0) print *,'Doing Adjoint Source'

    open(356,file=directory//'master.pixels',action='read', position='rewind')
    do i=1,nmasters
        if (indexnum .ne. i) read(356,*)
        if (indexnum .eq. i) read(356,*) x00
    enddo
    close(356)

    distances = abs((x-0.5)*xlength*10.0**(-8.) - x00)
    signed= ((x-0.5)*xlength*10.0**(-8.) - x00)
    idiff = floor(x00/((x(2)-x(1))*xlength*10.0**(-8.)))
    print *,nt
    open(44,file=directory//'forward_src'//contrib//'_ls'//jobno//'/timeline',action='read')
    do i=1,nt
        read(44,*) loc, t(nt-i+1)
    enddo
    close(44)
    t = -t/60.0
    dt = t(2)-t(1)



    adj = 0.0
    nmeasurements = 0.0

    call distmat(nt,1,freqnu)
    freqnu = freqnu/(nt*dt*60.)
    dnu = freqnu(2) - freqnu(1)

    call dfftw_plan_dft_3d(fwdplantemp, nx, dim2(rank), nt, acc, filtout, -1, FFTW_ESTIMATE)
    call dfftw_plan_dft_3d(fwdplantemp, nx, dim2(rank), nt, acc, filtout, -1, FFTW_ESTIMATE)
    call dfftw_plan_dft_3d(invplantemp, nx, dim2(rank), nt, filtout, acc, 1, FFTW_ESTIMATE)
    call dfftw_plan_dft_3d(invplantemp2, nx, dim2(rank), nt, filtout, ccdot, 1, FFTW_ESTIMATE)
    call dfftw_plan_dft_3d(invplantemp3, nx, dim2(rank), nt, filtemp, filtex, 1, FFTW_ESTIMATE)

    call dfftw_plan_dft_3d(fwdplandata, nx, dim2(rank), nt, dat, filtdat, -1, FFTW_ESTIMATE)
    call dfftw_plan_dft_3d(invplandata, nx, dim2(rank), nt, filtdat, dat, 1, FFTW_ESTIMATE)

    call dfftw_plan_dft_1d(onedplan, nt, oned, ccdotone, -1, FFTW_ESTIMATE)

    call readfits(directory//'forward_src'//contrib//'_ls'//jobno//'/vz_cc.fits', temparr, nt)
    acc = cmplx(temparr)
    call readfits(directory//'data/'//contrib//'.fits', temparr, nt)
    dat = cmplx(temparr)

    call dfftw_execute(fwdplantemp)
    call dfftw_execute(fwdplandata)

    filt = 1.0
    call fmode_filter(nt, fmode)
    call p1mode_filter(nt, p1mode)
    call p2mode_filter(nt, p2mode)
    call p3mode_filter(nt, p3mode)
    call p4mode_filter(nt, p4mode)
    call p5mode_filter(nt, p5mode)
    call p6mode_filter(nt, p6mode)
    call p7mode_filter(nt, p7mode)
    call all_pmode_filter(nt, all_pmode)

    inquire(file=directory//'filter.params.allmodes', exist=lexist)
    if (lexist) then
    open(94,file=directory//'filter.params.allmodes', action='read',position='rewind')
        read(94,*) leftcorner ! mHz
        read(94,*) rightcorner  ! mHz
        close(94)

        call freq_filter(leftcorner, rightcorner, nt, filt)
    endif

    do i=1,nt
        fmode(:,1,i)     =     fmode(:,1,i) * filt(i)!* UNKNOWN
        p1mode(:,1,i)    =    p1mode(:,1,i) * filt(i)!* UNKNOWN
        p2mode(:,1,i)    =    p2mode(:,1,i) * filt(i)!* UNKNOWN
        p3mode(:,1,i)    =    p3mode(:,1,i) * filt(i)!* UNKNOWN
        p4mode(:,1,i)    =    p4mode(:,1,i) * filt(i)!* UNKNOWN
        p5mode(:,1,i)    =    p5mode(:,1,i) * filt(i)!* UNKNOWN
        p6mode(:,1,i)    =    p6mode(:,1,i) * filt(i)!* UNKNOWN
        p7mode(:,1,i)    =    p7mode(:,1,i) * filt(i)!* UNKNOWN
        all_pmode(:,1,i) = all_pmode(:,1,i) * filt(i)!* UNKNOWN
    end do

    inquire(file=directory//'filter.params.0', exist=lexist)
    if (lexist) then

        open(94,file=directory//'filter.params.0', action='read',position='rewind')
        read(94,*) leftcorner ! mHz
        read(94,*) rightcorner  ! mHz
        close(94)

        call freq_filter(leftcorner, rightcorner, nt, filt)
    endif

    do i=1,nt
        fmode(:,1,i) = fmode(:,1,i) * filt(i)!* UNKNOWN
    enddo

    call readfits(directory//'forward_src'//contrib//'_ls'//jobno//'/vz_cc.fits', temparr, nt)
    acc = cmplx(temparr)
    call readfits(directory//'data/'//contrib//'.fits', temparr, nt)
    dat = cmplx(temparr)

    inquire(file=directory//'filter.params.1', exist=lexist)
    if (lexist) then

        open(94,file=directory//'filter.params.1', action='read',position='rewind')
        read(94,*) leftcorner ! mHz
        read(94,*) rightcorner  ! mHz
        close(94)

        call freq_filter(leftcorner, rightcorner, nt, filt)
    endif

    do i=1,nt
        p1mode(:,1,i) = p1mode(:,1,i) * filt(i)!* UNKNOWN
    enddo


    inquire(file=directory//'filter.params.2', exist=lexist)
    if (lexist) then

        open(94,file=directory//'filter.params.2', action='read',position='rewind')
        read(94,*) leftcorner ! mHz
        read(94,*) rightcorner  ! mHz
        close(94)

        call freq_filter(leftcorner, rightcorner, nt, filt)
    endif

    do i=1,nt
        p2mode(:,1,i) = p2mode(:,1,i) * filt(i)!* UNKNOWN
    enddo

    tempout = filtout
    tempdat = filtdat

    vel = 0
    vel(0) = 0.44
    vel(1) = 0.60
    vel(2) = 0.75
    vel(3) = 0.9
    vel(4) = 1.2
    vel(5) = 1.4
    vel(6) = 1.7
    vel(7) = 1.9

    !RIDGE FILTERS
    do pord=0,8
        taus = 0.0
        call convert_to_string(pord, ord, 1)
        inquire(file=directory//'params.'//ord, exist=lexist)
        if (lexist) then
            open(97, file = directory//'params.'//ord, action = 'read', status='old')

            read(97,*) mindist
            read(97,*) maxdist
            read(97,*) window
            close(97)


            open(238, file = directory//'forward_src'//contrib//'_ls'//jobno//'/ttdiff.'//ord, action = 'write')

            halftime = nint(window/(2.*dt))
            leng = 2*halftime+1

            if (pord==0) then
                filter = fmode
                if (.not. linesearch) call writefits_3d('fmode_filter.fits',fmode,nt)
            else if (pord==1) then
                filter = p1mode
                call writefits_3d('p1mode_filter.fits',p1mode,nt)
            else if (pord==2) then
                call writefits_3d('p2mode_filter.fits',p2mode,nt)
                filter = p2mode
            else if (pord==3) then
                filter = p3mode
                call writefits_3d('p3mode_filter.fits',p3mode,nt)
            else if (pord==4) then
                filter = p4mode
                call writefits_3d('p4mode_filter.fits',p4mode,nt)
            else if (pord==5) then
                filter = p5mode
                call writefits_3d('p5mode_filter.fits',p5mode,nt)
            else if (pord==6) then
                filter = p6mode
                call writefits_3d('p6mode_filter.fits',p6mode,nt)
            else if (pord==7) then
                filter = p7mode
                call writefits_3d('p7mode_filter.fits',p7mode,nt)
            else if (pord==8) then
                filter = all_pmode
                call writefits_3d('first_bounce_pmode_filter.fits',all_pmode,nt)
            end if
            filtout = tempout * cmplx(filter)
            filtdat = tempdat * cmplx(filter)
            call dfftw_execute(invplantemp)
            call dfftw_execute(invplandata)

            con = 2.0*pi


!~             inquire(file='wavespeed',exist=ws_exist)
!~             if (ws_exist) then
!~                 open(3378,file="wavespeed",action="read")
!~                     do i=1,nmasters
!~                     call convert_to_string(i, contribmatch, 2)
!~                     if (contribmatch == contrib) then
!~                     read(3378,*) vel(0), vel(1), vel(2),vel(3),vel(4),&
!~                     vel(5),vel(6),vel(7) ,vel(8), vel(9),vel(10)
!~                     else
!~                         read(3378,*)
!~                     end if
!~                     end do
!~                 close(3378)
!~             end if

!~             if (pord==0) print *,"Using mode velocities",vel

            do i=1,nt
                filtout(:,1,i) = filtout(:,1,i) * eye * freqnu(i) * con
            enddo

            call dfftw_execute(invplantemp2)


            if (rank==0 .and. (.not. (linesearch))) then
                open(596,file=directory//'forward_src'//contrib//'_ls00/windows.'//ord,action='write',status='replace')
            elseif(rank==0 .and. (linesearch)) then
                inquire(file=directory//'forward_src'//contrib//'_ls00/windows.'//ord,exist=exist_win)
                if (exist_win) &
                    open(596,file=directory//'forward_src'//contrib//'_ls00/windows.'//ord,action='read',status='old')
            endif


            timest = 1
            timefin = nt
!~             RECEIVER PIXEL FLAG
            do i=1,nx
!~             RECEIVER PIXEL END FLAG
!~                 print *,"Using i =",i,"dist =",distances(i),"to compute misfits"
                if ((distances(i) > mindist) .and. (distances(i) < maxdist)) then

                    if (.not. linesearch) then
                        if (pord .ne. 8) then
                            timest = floor(distances(i)/vel(pord) * 1./dt)
                            timefin = timest + 40
                        else if (pord .eq. 8) then
                            dist = distances(i)
                            timest = floor((-0.000417147774671*dist**2&
                            +0.313350998096*dist+17.9609631186)* 1./(dt))
                            timefin = floor((-0.00028361249034*dist**2&
                            +0.29270337114*dist+36.4398734578)* 1./(dt))
                        end if

                        loc = maxloc(abs(real(acc(i,1,timest:timefin))),1)+timest-1
                        lef = max(loc - halftime,1)
                        rig = min(loc + halftime,nt)

                        inquire(unit=596,opened=file_open)
                        if (file_open) write(596,*) lef, rig

                    else
!~                       Read in windows
                        inquire(unit=596,opened=file_open)
                        if (file_open) read(596,*) lef, rig
                    endif

                    ! call compute_tt(real(acc(i,1,lef:rig)),real(dat(i,1,lef:rig)),tau,dt,leng)
                     call compute_tt_gizonbirch(real(acc(i,1,:)),real(dat(i,1,:)),tau,dt,nt, lef, rig)

138                 format (I3,1X,F14.8,1X,I4,1X,I4,1X,I4,1X,I4,1X,I4)

                    write(238,138) i,tau*60.,lef,rig,loc,timest,timefin

                    windows(:) = 0.0
                    windows(lef:rig) = 1.0
                    oned = ccdot(i,1,:) * cmplx(windows)

                    call dfftw_execute(onedplan) ! oned -> ccdotone
                    do j=1,nt
                        filtemp(:,1,j) = cmplx(filter(:,1,j))  * ccdotone(j) * exp(-eyekh*(x(i)))
                    enddo
                    call dfftw_execute(invplantemp3) ! filtemp -> filtex

                    misfit = misfit + tau**2. !* iwls_misfit_factor

                    nmeasurements = nmeasurements + 1
                    con = -tau/(sum(ccdot(i,1,lef:rig)**2.)*dt) !* iwls_misfit_factor !* sign(1.0,signed(i))
                    do j=1,nt
                        adj(:,1,nt-j+1) = real(filtex(:,1,j) * con) + adj(:,1,nt-j+1)
                    enddo
                endif ! if distances are within range
            enddo

            inquire(unit=238,opened=file_open)
            if (file_open) close(238)
            inquire(unit=237,opened=file_open)
            if (file_open) close(237)

            inquire(unit=596,opened=file_open)
            if (file_open) close(596)
        endif


    enddo

    pcoef = 0.0

    pcoef(1,1) =  22.5000
    pcoef(2,1) =   0.235317
    pcoef(3,1) =   0.00134167
    pcoef(4,1) = -1.32222e-05
    pcoef(5,1) = 3.00000e-08

    pcoef(1,2) = 22.6889
    pcoef(2,2) = 0.83
    pcoef(3,2) = -0.0032222

    pcoef(1,3) = 25.0
    pcoef(2,3) = 0.9
    pcoef(3,3) = -0.002222


!~     highpmode = 1.0             ! ???
!~
!~     do i=1,nt
!~         highpmode(:,1,i) = highpmode(:,1,i) * filt(i) !* UNKNOWN
!~     enddo

    !PHASE-SPEED FILTERS
!~     do pord=9,9
!~         taus=0.0
!~         call convert_to_string(pord, ord, 1)
!~         inquire(file=directory//'params.'//ord, exist=lexist)
!~         if (lexist) then
!~             open(97, file = directory//'params.'//ord, action = 'read', status='old')
!~             read(97,*) mindist
!~             read(97,*) maxdist
!~             read(97,*) window
!~             close(97)
!~
!~             halftime = nint(window/(2.*dt))
!~             leng = 2*halftime+1
!~
!~             filtout = tempout * cmplx(highpmode) ! * UNKNOWN
!~             filtdat = tempdat * cmplx(highpmode) ! * UNKNOWN
!~
!~             call dfftw_execute(invplantemp)
!~             call dfftw_execute(invplandata)
!~
!~             con = 2.0*pi!/(dble(nt)*dble(nx))
!~
!~             do i=1,nt
!~                 filtout(:,1,i) = filtout(:,1,i) * eye * freqnu(i) * con
!~             enddo
!~
!~             call dfftw_execute(invplantemp2)
!~
!~             do i=1,nx
!~                 if ((distances(i) > mindist) .and. (distances(i) < maxdist)) then
!~
!~                     if (distances(i) < 250.0) then
!~                         timesmax = nint((pcoef(1,pord-1) - t(1) + pcoef(2,pord-1)*distances(i) + pcoef(3,pord-1)*distances(i)**2. &
!~                         + pcoef(4,pord-1)*distances(i)**3. + pcoef(5,pord-1)*distances(i)**4.)/dt) + 1
!~
!~                         loc = maxloc(real(acc(i,1,(timesmax-6):(timesmax+6))),1) + timesmax - 7
!~                     else
!~                         loc = maxloc(real(acc(i,1,:)),1)
!~                     endif
!~                     lef = loc - halftime
!~                     rig = loc + halftime
!~                     call compute_tt(real(acc(i,1,lef:rig)),real(dat(i,1,lef:rig)),tau,dt,leng)
!~
!~                     windows(:) = 0.0
!~                     windows(lef:rig) = 1.0
!~                     oned = ccdot(i,1,:) * cmplx(windows)
!~                     call dfftw_execute(onedplan)
!~                     do j=1,nt
!~                         filtemp(:,1,j) = cmplx(highpmode(:,1,j))  * ccdotone(j) * exp(-eyekh*(x(i)))
!~                     enddo
!~                     call dfftw_execute(invplantemp3)
!~
!~                     misfit = misfit + tau**2.
!~                     nmeasurements = nmeasurements + 1
!~
!~                     con = -tau/(sum(ccdot(i,1,lef:rig)**2.)*dt) !* sign(1.0,signed(i))
!~                     do j=1,nt
!~                         adj(:,1,nt-j+1) = real(filtex(:,1,j) * con) + adj(:,1,nt-j+1)
!~                     enddo
!~                 endif !Mindist, maxdist
!~             enddo ! end of the do-loop for mindist,maxdist
!~         endif ! if params.p# exists
!~     enddo ! end of pord loop

    adj = adj *nt / 60.

    if (.not. (linesearch .or. compute_synth)) then
        inquire(file=directory//'adjoint_src'//contrib, exist=lexist)
        if (.not. lexist .and. (rank==0)) &
            call system('mkdir '//directory//'adjoint_src'//contrib)

        call writefits_3d(directory//'adjoint_src'//contrib//'/source.fits',adj,nt)
    endif

    call MPI_REDUCE(misfit, misfit_tot, 1, MPI_DOUBLE_PRECISION, &
                MPI_SUM, 0, MPI_COMM_WORLD, ierr)

    call MPI_REDUCE(nmeasurements, nmeas, 1, MPI_INTEGER, &
            MPI_SUM, 0, MPI_COMM_WORLD, ierr)


    if (rank==0) then
        if (.not. linesearch) then
            open(34, file=directory//'adjoint_src'//contrib//'/information', &
                action='write', position='rewind', status='replace')
            write(34,*) nt
            write(34,*) t(1)*60.
            write(34,*) t(nt)*60.
            if (test_in_2d) write(34,*) '2D'
            close(34)
        endif
        if (rank==0) open(88, file=directory//'kernel/misfit_'//contrib//'_'//jobno, status='unknown',&
            action='write')

        misfit_tot = misfit_tot * 0.5
        print *, misfit_tot
        write(88,*) indexnum, nmeas, misfit_tot
        close(88)

    endif

    call dfftw_destroy_plan(fwdplantemp)
    call dfftw_destroy_plan(fwdplandata)

    call dfftw_destroy_plan(invplantemp3)
    call dfftw_destroy_plan(invplantemp2)
    call dfftw_destroy_plan(invplantemp)
    call dfftw_destroy_plan(invplandata)
    call dfftw_destroy_plan(onedplan)

END SUBROUTINE ADJOINT_SOURCE_FILT

!================================================================================
! subroutine bessel(order,arg,outp)
!! use ifport
! implicit none
! integer, intent(in) :: order
! real(kind=real64), intent(in) :: arg
! real(kind=real64), intent(out) :: outp

! outp = BesJN(order,arg)

! end subroutine bessel

!================================================================================

SUBROUTINE FOURIER_SMOOTH_X(input_arr,nk,output_arr)
    use initialize
    implicit none
    integer, intent(in) :: nk
    real(kind=real64), dimension(nx,1,nz), intent(in) :: input_arr
    real(kind=real64), dimension(nx,1,nz), intent(out) :: output_arr
    real(kind=real64) psi_row(nx),sigmak_smooth,smoothing_function(nx/2+1),k_smooth
    complex(kind=real64) psi_fft(nx/2+1)
    integer zind,xind
    integer(kind=8) plan_fwd,plan_inv

    call dfftw_plan_dft_r2c_1d(plan_fwd,nx,psi_row,psi_fft,FFTW_ESTIMATE)
    call dfftw_plan_dft_c2r_1d(plan_inv,nx,psi_fft,psi_row,FFTW_ESTIMATE)

    sigmak_smooth = 2*pi*nk/(xlength/1D8)
    do xind=1,nx/2+1
        k_smooth = 2*pi*(xind-1)/(xlength/1D8)
        smoothing_function(xind) = exp(-k_smooth**2/(2*sigmak_smooth**2))
!~         print *,"k",k_smooth,"smoothing fn",exp(-k_smooth**2/(2*sigmak_smooth**2))
    end do

    do zind=1,nz
        psi_row = input_arr(:,1,zind)
        call dfftw_execute_dft_r2c(plan_fwd, psi_row, psi_fft)

        psi_fft=psi_fft*smoothing_function
        call dfftw_execute_dft_c2r(plan_inv,psi_fft,psi_row)
        output_arr(:,1,zind) = psi_row/nx
    end do

    call dfftw_destroy_plan(plan_fwd);
    call dfftw_destroy_plan(plan_inv);

END SUBROUTINE FOURIER_SMOOTH_X

!================================================================================




SUBROUTINE MISFIT_ALL(nt)
    use initialize
    use all_modules
    implicit none
    integer(kind=8) fwdplantemp, invplantemp, invplantemp2
    integer(kind=8)  invplantemp3, fwdplandata, invplandata, onedplan
    logical lexist,exist_win,file_open
    integer i, nt, indexnum, pord, timesmax, halftime, bounce, freqfilts,timest
    integer loc, leng, lef, rig, j, nmeasurements(0:10), nmeas(0:10), ierr,timefin
    character(len=1) ord,ffstr
    real(kind=real64) mindist, maxdist, window, distances(nx), x00, t(nt),  tau, vel(0:10)
    real(kind=real64),dimension(nx,dim2(rank),nt):: pmode,p2mode,fmode,all_else,filter,phase,temparr,&
                                      all_pmode,p3mode,p4mode,p5mode,p6mode,p7mode
    real(kind=real64) pcoef(5,4), adj(nx,dim2(rank),nt), windows(nt), filt(nt),dt,xdim(nx)
    real(kind=real64) freqnu(nt), leftcorner, rightcorner, dnu, con, misfit(0:10),misfit_tot(0:10)
    complex(kind=real64), dimension(nx,dim2(rank),nt) :: filtout, tempout, tempdat, &
                                 filtdat, filtex, filtemp,dat,acc,ccdot
    complex(kind=real64), dimension(nx) :: eyekh
    complex(kind=real64), dimension(nt) :: oned, ccdotone
    complex(kind=real64) UNKNOWN
    real(kind=real64) speed, var, param(6), d_i
    real(kind=real64) ign1,ign2
    character(len=2) fnum
    character(len=1) frqnum

    UNKNOWN = 1.0/0.55
    filtout = 1.0
    pcoef = 0.0
    misfit = 0.0
    call distmat(nx,1,distances)

    dx = x(2)-x(1)
    eyekh= cmplx(0,1)*distances*2.*pi
    ! eyekh(nx/2+1) = 0.0
    read(contrib,*) indexnum

    if (rank==0) print *,'Doing Adjoint Source'

    open(356,file=directory//'master.pixels',action='read', position='rewind')
    do i=1,nmasters
        if (indexnum .ne. i) read(356,*)
        if (indexnum .eq. i) read(356,*) x00
    enddo
    close(356)

    distances = abs((x-0.5)*xlength*10.0**(-8.) - x00)

    open(44,file=directory//'forward_src'//contrib//'_ls'//jobno//'/timeline',action='read')
    do i=1,nt
        read(44,*) loc, t(nt-i+1)
    enddo
    close(44)
    t = -t/60.0
    dt = t(2)-t(1)


    adj = 0.0
    nmeasurements = 0.0

    call distmat(nt,1,freqnu)
    freqnu = freqnu/(nt*dt*60.)
    dnu = freqnu(2) - freqnu(1)

    call dfftw_plan_dft_3d(fwdplantemp, nx, dim2(rank), nt, acc, filtout, -1, FFTW_ESTIMATE)
    call dfftw_plan_dft_3d(invplantemp, nx, dim2(rank), nt, filtout, acc, 1, FFTW_ESTIMATE)
    call dfftw_plan_dft_3d(invplantemp2, nx, dim2(rank), nt, filtout, ccdot, 1, FFTW_ESTIMATE)
    call dfftw_plan_dft_3d(invplantemp3, nx, dim2(rank), nt, filtemp, filtex, 1, FFTW_ESTIMATE)

    call dfftw_plan_dft_3d(fwdplandata, nx, dim2(rank), nt, dat, filtdat, -1, FFTW_ESTIMATE)
    call dfftw_plan_dft_3d(invplandata, nx, dim2(rank), nt, filtdat, dat, 1, FFTW_ESTIMATE)

    call dfftw_plan_dft_1d(onedplan, nt, oned, ccdotone, -1, FFTW_ESTIMATE)




    call readfits(directory//'forward_src'//contrib//'_ls'//jobno//'/vz_cc.fits', temparr, nt)
    acc = temparr

    call readfits(directory//'data/'//contrib//'.fits', temparr, nt)
    dat = temparr
    call dfftw_execute(fwdplantemp)
    call dfftw_execute(fwdplandata)

    call dfftw_destroy_plan(fwdplantemp)
    call dfftw_destroy_plan(fwdplandata)


    if (rank==0) open(543,file=directory//'kernel/misfit_all_'//contrib//'_'//jobno,&
                    action='write')

    tempout = filtout
    tempdat = filtdat

    vel = 0
    vel(0) = 0.44
    vel(1) = 0.60
    vel(2) = 0.75
    vel(3) = 0.9
    vel(4) = 1.2
    vel(5) = 1.4
    vel(6) = 1.7
    vel(7) = 1.9



    do freqfilts=1,3

        if (freqfilts==1) then
            leftcorner = 2.0
            rightcorner = 3.5
        elseif (freqfilts==2) then
            leftcorner = 3.5
            rightcorner = 4.5
        elseif (freqfilts==3) then
            leftcorner = 4.5
            rightcorner = 6.0
        endif

        call freq_filter(leftcorner, rightcorner, nt, filt)

        call fmode_filter(nt, fmode)
        call p1mode_filter(nt, pmode)
        call p2mode_filter(nt, p2mode)
        call p3mode_filter(nt, p3mode)
        call p4mode_filter(nt, p4mode)
        call p5mode_filter(nt, p5mode)
        call p6mode_filter(nt, p6mode)
        call p7mode_filter(nt, p7mode)
        call all_pmode_filter(nt, all_pmode)
!~         all_else = 1. - fmode - pmode          ! ???
!~         highpmode = 1.0

        do i=1,nt
            fmode(:,1,i) = fmode(:,1,i) * filt(i)!* UNKNOWN
            pmode(:,1,i) = pmode(:,1,i) * filt(i)!* UNKNOWN
            p2mode(:,1,i) = p2mode(:,1,i) * filt(i)!* UNKNOWN
            p3mode(:,1,i) = p3mode(:,1,i) * filt(i)!* UNKNOWN
            p4mode(:,1,i) = p4mode(:,1,i) * filt(i)!* UNKNOWN
            p5mode(:,1,i) = p5mode(:,1,i) * filt(i)!* UNKNOWN
            p6mode(:,1,i) = p6mode(:,1,i) * filt(i)!* UNKNOWN
            p7mode(:,1,i) = p7mode(:,1,i) * filt(i)!* UNKNOWN
            all_pmode(:,1,i) = all_pmode(:,1,i) * filt(i) !* UNKNOWN
        enddo

        misfit = 0.0
        nmeasurements = 0
        misfit_tot =0.0
        nmeas=0

        !RIDGE FILTERS, f, p1, ...
        do pord=0,8
            call convert_to_string(pord, ord, 1)
            inquire(file=directory//'params.'//ord, exist=lexist)
            if (lexist) then
                open(97, file = directory//'params.'//ord, action = 'read', status='old')
                read(97,*) mindist
                read(97,*) maxdist
                read(97,*) window
                close(97)

                halftime = nint(window/(2.*dt))
                leng = 2*halftime+1

                if (rank==0 .and. (.not. linesearch)) then
                    call convert_to_string(pord, fnum, 2)
                    call convert_to_string(freqfilts, frqnum, 1)
                    open(20380,file=directory//'forward_src'//contrib//'_ls00'//&
                    '/windows.all.'//fnum//'.'//frqnum,action='write',status='replace')
                elseif (rank==0 .and. linesearch) then

                    call convert_to_string(pord, fnum, 2)
                    call convert_to_string(freqfilts, frqnum, 1)

                    inquire(file=directory//'forward_src'//contrib//'_ls00'//&
                    '/windows.all.'//fnum//'.'//frqnum,exist=exist_win)
                    if (exist_win) &
                    open(20380,file=directory//'forward_src'//contrib//'_ls00'//&
                    '/windows.all.'//fnum//'.'//frqnum,action='read',status='old')

                endif

                if (pord==0) filter = fmode
                if (pord==1) filter = pmode
                if (pord==2) filter = p2mode
                if (pord==3) filter = p3mode
                if (pord==4) filter = p4mode
                if (pord==5) filter = p5mode
                if (pord==6) filter = p6mode
                if (pord==7) filter = p7mode
                if (pord==8) filter = all_pmode
                filtout = tempout * cmplx(filter)
                filtdat = tempdat * cmplx(filter)

                call dfftw_execute(invplantemp)
                call dfftw_execute(invplandata)

                con = 2.0*pi

                do i=1,nt
                    filtout(:,1,i) = filtout(:,1,i) * eye * freqnu(i) * con
                enddo
                call dfftw_execute(invplantemp2)

                do i=1,nx
                    if ((distances(i) > mindist) .and. (distances(i) < maxdist)) then

                        if (.not. linesearch) then
                            if (pord .ne. 8) then
                            timest = floor(distances(i)/vel(pord) * 1./dt)
                            timefin = timest + 40
                            else
                            d_i = distances(i)
                            timest = floor((-1.45511095e-04*(d_i)**2 &
                                    + 2.72140632e-01*d_i + 2.16331969e+01 )* 1./dt)

                            timefin = floor((3.85746516e-08*(d_i)**2 &
                                    + 2.21324786e-01*d_i + 4.36702155e+01 )* 1./dt)
                            end if

                            loc = maxloc(abs(real(acc(i,1,timest:timefin))),1)+timest-1
                            lef = loc - halftime
                            rig = loc + halftime

                            inquire(unit=20380,opened=file_open)
                            if (file_open) write(20380,*) lef, rig

                        elseif (linesearch) then
                            inquire(unit=20380,opened=file_open)
                            if (file_open) read(20380,*) lef, rig
                        endif

                        call compute_tt(real(acc(i,1,lef:rig)),real(dat(i,1,lef:rig)),tau,dt,leng)

                        misfit(pord) = misfit(pord) + tau**2.
                        nmeasurements(pord) = nmeasurements(pord) + 1

                    endif
                enddo

                if (rank==0) then
                    inquire(unit=20380,opened=file_open)
                    if (file_open) close(20380)

                end if
            endif

            call MPI_REDUCE(misfit(pord), misfit_tot(pord), 1, MPI_DOUBLE_PRECISION, &
                        MPI_SUM, 0, MPI_COMM_WORLD, ierr)

            call MPI_REDUCE(nmeasurements(pord), nmeas(pord), 1, MPI_INTEGER, &
                        MPI_SUM, 0, MPI_COMM_WORLD, ierr)

        enddo


!~         pcoef = 0.0
!~
!~         pcoef(1,1) =  22.5000
!~         pcoef(2,1) =   0.235317
!~         pcoef(3,1) =   0.00134167
!~         pcoef(4,1) = -1.32222e-05
!~         pcoef(5,1) = 3.00000e-08
!~
!~         pcoef(1,2) = 22.6889
!~         pcoef(2,2) = 0.83
!~         pcoef(3,2) = -0.0032222
!~
!~         pcoef(1,3) = 25.0
!~         pcoef(2,3) = 0.9
!~         pcoef(3,3) = -0.002222


!~         do pord=9,9
!~             call convert_to_string(pord, ord, 1)
!~             inquire(file=directory//'params.'//ord, exist=lexist)
!~             if (lexist) then
!~                 open(97, file = directory//'params.'//ord, action = 'read', status='old')
!~                 read(97,*) mindist
!~                 read(97,*) maxdist
!~                 read(97,*) window
!~                 close(97)
!~
!~                 halftime = nint(window/(2.*dt))
!~                 leng = 2*halftime+1
!~
!~                 filtout = tempout * cmplx(highpmode) ! * UNKNOWN
!~                 filtdat = tempdat * cmplx(highpmode) ! * UNKNOWN
!~
!~                 call dfftw_execute(invplantemp)
!~                 call dfftw_execute(invplandata)
!~
!~                 con = 2.0*pi!/(dble(nt)*dble(nx))
!~
!~                 do i=1,nt
!~                     filtout(:,1,i) = filtout(:,1,i) * eye * freqnu(i) * con
!~                 enddo
!~
!~                 call dfftw_execute(invplantemp2)
!~
!~                 do i=1,nx
!~                     if ((distances(i) > mindist) .and. (distances(i) < maxdist)) then
!~
!~                         if (distances(i) < 250.0) then
!~                             timesmax = nint((pcoef(1,pord-1) - t(1) + &
!~                              pcoef(2,pord-1)*distances(i) + &
!~                              pcoef(3,pord-1)*distances(i)**2. &
!~                             + pcoef(4,pord-1)*distances(i)**3. + &
!~                             pcoef(5,pord-1)*distances(i)**4.)/dt) + 1
!~
!~                             loc = maxloc(real(acc(i,1,(timesmax-6):(timesmax+6))),1) &
!~                             + timesmax - 7
!~                         else
!~                             loc = maxloc(real(acc(i,1,:)),1)
!~                         endif
!~                         lef = loc - halftime
!~                         rig = loc + halftime
!~
!~                         call compute_tt(real(acc(i,1,lef:rig)),real(dat(i,1,lef:rig)),tau,dt,leng)
!~
!~
!~                         misfit(pord) = misfit(pord) + tau**2.
!~                         nmeasurements(pord) = nmeasurements(pord) + 1
!~
!~                     endif !Mindist, maxdist
!~                 enddo ! end of the do-loop for mindist,maxdist
!~             endif ! if params.p# exists
!~
!~             call MPI_REDUCE(misfit(pord), misfit_tot(pord), 1, MPI_DOUBLE_PRECISION, &
!~                         MPI_SUM, 0, MPI_COMM_WORLD, ierr)
!~
!~             call MPI_REDUCE(nmeasurements(pord), nmeas(pord), 1, MPI_INTEGER, &
!~                         MPI_SUM, 0, MPI_COMM_WORLD, ierr)
!~
!~         enddo ! end of pord loop

        if (rank==0) then
            write(543,*) "#",indexnum,leftcorner, rightcorner
            write(543,*) "#",nmeas
            write(543,*) misfit_tot*0.5
            print *,misfit_tot*0.5

        endif
    enddo ! END OF FREQFILTS LOOP

    if (rank==0) then
        inquire(unit=543,opened=file_open)
        if (file_open) close(543)
    end if

    call dfftw_destroy_plan(invplantemp3)
    call dfftw_destroy_plan(invplantemp2)
    call dfftw_destroy_plan(invplantemp)
    call dfftw_destroy_plan(invplandata)
    call dfftw_destroy_plan(onedplan)


END SUBROUTINE MISFIT_ALL


!================================================================================
