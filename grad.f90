PROGRAM CONJGRAD

 
    implicit none
    include 'fftw3.f'
    include 'params.i'
    integer*8 fftw_x, fftw_z, fwdplan, invplan
    integer i, k, iteration, nmasters, iterati, start, itermax,stencilBFGS,minBFGS
    parameter(itermax = 25, stencilBFGS = 4)
    real*8 x(nx), smooth_x, smooth_z, ran, z0(nz),con,den,back(6,nz), eps
    parameter(smooth_x = 4, smooth_z = 4)
    real*8, dimension(:), allocatable :: alph
    real*8, dimension(nx,1,nz) :: totkern_c, totkern_psi, b, &
    kern, update_c, update, update0, update_psi, &
    lastupdate, lastmodel, model, gradc, lastgradc,  &
    lastgradpsi, gradpsi, hess, lastmodel_psi, lastmodel_c
    complex*16 filtx(nx/2+1), filtz(nz), temp(nx/2+1, 1, nz), z(nz)
    character*2 iter, itermin, pixnum, chartem, chartminone
    character*1 charnum
    real*8 data(nz,6)
    ! character*80 directory
    logical FORCE_STEEPEST, LBFGS, lexist


    ! directory = '/nobackup/shanasog/classic/'
    FORCE_STEEPEST = .true.
    LBFGS = .false.
    iteration = 1
    write(iter, '(I2.2)') iteration 
    inquire(file= adjustl(trim(directory))//'update/model_psi_'//iter//'.fits', exist = lexist)

    do while (lexist .and. (iteration < itermax))
        iteration = iteration + 1
        write(iter, '(I2.2)') iteration 
        inquire(file= adjustl(trim(directory))//'update/model_psi_'//iter//'.fits', exist = lexist)
    enddo

    if (LBFGS .and. (.not. FORCE_STEEPEST)) then
        minBFGS = iteration - 1 - stencilBFGS
        if (minBFGS < 1) then
            print *,'BFGS stencil too wide'
            stop
        endif
    endif 
    print *,iteration

    open(44,file='solar_deep',position='rewind',action='read')
    do i=1,nz
        read(44,*) back(:,i)
    enddo
    close(44)

    call dfftw_plan_dft_r2c_1d(fftw_x, nx, x, filtx, FFTW_ESTIMATE)
    call dfftw_plan_dft_1d(fftw_z, nz, z, filtz, FFTW_FORWARD, FFTW_ESTIMATE)

    call distmat(nx,1,x)
    call distmat(nz,1,z0)
    x = exp(-x**2./(2.*smooth_x**2.))
    z = exp(-z0**2./(2.*smooth_z**2.))

    call dfftw_execute(fftw_x)
    call dfftw_execute(fftw_z)

    call dfftw_destroy_plan(fftw_x)
    call dfftw_destroy_plan(fftw_z)


    write(iter, '(I2.2)') iteration
    write(itermin, '(I2.2)') iteration-1
 
    call dfftw_plan_dft_r2c_3d(fwdplan, nx, 1, nz, kern, temp, FFTW_ESTIMATE)
    call dfftw_plan_dft_c2r_3d(invplan, nx, 1, nz, temp, kern, FFTW_ESTIMATE)
    open(356,file=adjustl(trim(directory))//'master.pixels',action='read',position='rewind')
    do i=1,100
        read(356,*,end=24432)
    enddo
    24432 close(356)
    nmasters = i-1

    totkern_c = 0.0
    totkern_psi = 0.0
    ! totkern_vz = 0.0

    do i=1,nmasters

        write(pixnum, '(I2.2)') i 
        print *,pixnum,i

        call readfits(adjustl(trim(directory))//'kernel/kernel_c_'//pixnum//'.fits',kern, nx, 1, nz)
        call filter_z(kern)
        !kern(:,1,1:10) = 0.0
        !kern(:,1,nz-9:nz) = 0.0
       
        totkern_c = totkern_c + kern
        !call writefits_3d(adjustl(trim(directory))//'update/totkern_c_'//pixnum//'.fits',totkern_c,nx, 1,  nz)
        
        !call countzeros(totkern_c,'totkern_c')
        
        call readfits(adjustl(trim(directory))//'kernel/kernel_psi_'//pixnum//'.fits',kern, nx, 1, nz)
        call filter_z(kern)
        !kern(:,1,1:10) = 0.0
        !kern(:,1,nz-9:nz) = 0.0
        totkern_psi = totkern_psi + kern
        !call writefits_3d(adjustl(trim(directory))//'update/totkern_psi_'//pixnum//'.fits',totkern_psi,nx, 1,  nz)
	    !call countzeros(totkern_psi,'totkern_psi')
        
        call readfits(adjustl(trim(directory))//'kernel/hessian_'//pixnum//'.fits',kern, nx, 1, nz)
        call filter_z(kern)
        hess = hess + abs(kern)
        !call writefits_3d(adjustl(trim(directory))//'update/hess'//pixnum//'.fits',hess,nx, 1,  nz)
        !call countzeros(hess,'hess')
	
        !  call readfits(adjustl(trim(directory))//'kernel/kernel_vz_'//pixnum//'.fits',kern, nx, 1, nz)
        !  kern(:,1,1:10) = 0.0
        !  kern(:,1,nz-9:nz) = 0.0
        !  totkern_vz = totkern_vz + kern

    enddo

    do k=1,nz
        hess(:,:,k) = hess(:,:,k) * back(3,k)
    enddo
    !hess = abs(hess)/maxval(abs(hess)) + con

    !hess = hess/con

    con = 0.005
    hess = (hess)/maxval(abs(hess))
    !call countzeros(hess,'hess')
    !call writefits_3d(adjustl(trim(directory))//'update/hesss.fits',hess,nx, 1,  nz)

    do k=1,nz
        do i=1,nx
            !if (hess(i,1,k) < 0.01 .and. hess(i,1,k) .ge. 0.) hess(i,1,k) = 0.01
            if (hess(i,1,k) < con) hess(i,1,k) = con
            !if (hess(i,1,k) > -0.01 .and. hess(i,1,k) < 0.) hess(i,1,k) = -0.01
        enddo
    enddo
    hess = hess/maxval(abs(hess))
  
    !call writefits_3d(adjustl(trim(directory))//'update/hessss.fits',hess,nx, 1,  nz)

    !hess = 1.0
    kern = totkern_c/hess
    !call writefits_3d(adjustl(trim(directory))//'update/kernn.fits',kern,nx, 1,  nz)
    !call countzeros(kern,'kern')
    call dfftw_execute(fwdplan)
    do k=1,nz
        temp(:,1,k) = temp(:,1,k)*filtx!*filtz(k)
    enddo
    call dfftw_execute(invplan)
    call filter_z(kern)
    totkern_c = kern

    kern = totkern_psi/hess
    !call countzeros(kern,'kern')
    !call writefits_3d(adjustl(trim(directory))//'update/kernnn.fits',kern,nx, 1,  nz)
    
    !do k=1,nz
    !    kern(:,1,k) = sum(kern(:,1,k))
    !enddo
    !call writefits_3d(adjustl(trim(directory))//'update/kernnnn.fits',kern,nx, 1,  nz)

    !ckern,'kern')
    
    call dfftw_execute(fwdplan)
    do k=1,nz
        temp(:,1,k) = temp(:,1,k)*filtx*filtz(k)

        ! con = temp(1,1,k)
        !temp(:,1,k) = 0.0
        !temp(1,1,k) = con*filtz(k)

    enddo
    call dfftw_execute(invplan)
    totkern_psi = kern
    
    !call countzeros(totkern_psi,'totkern_psi')
    
    !call writefits_3d(adjustl(trim(directory))//'update/totkernpsi.fits', totkern_psi,nx, 1,  nz)
    !  kern = totkern_vz/hess
    !  call dfftw_execute(fwdplan)
    !  do k=1,nz
    !   temp(:,1,k) = temp(:,1,k)*filtx*filtz(k)
    !  enddo
    !  call dfftw_execute(invplan)
    !  totkern_vz = kern


    call dfftw_destroy_plan(fwdplan)
    call dfftw_destroy_plan(invplan)

    
    call antisymmetrize(totkern_psi)
    do i=1,5
        call filter_z(totkern_psi)
    enddo
    call filter_z(totkern_c)  
    !call countzeros(totkern_psi,'totkern_psi')
    
    call writefits_3d(adjustl(trim(directory))//'update/gradient_c_'//itermin//'.fits',totkern_c,nx, 1, nz)
    call writefits_3d(adjustl(trim(directory))//'update/gradient_psi_'//itermin//'.fits',totkern_psi,nx, 1, nz)
    ! call writefits_3d(adjustl(trim(directory))//'update/gradient_vx_'//itermin//'.fits',totkern_vx,nx, 1, nz)
    ! call writefits_3d(adjustl(trim(directory))//'update/gradient_vz_'//itermin//'.fits',totkern_vz,nx, 1, nz)


    if (iteration .eq. 1 .or. (FORCE_STEEPEST)) then 

        !  if (iteration == 1) then
        !   open(22,file='~/testclassic/solar_model',action='read',position='rewind')
        !   do i=1,nz
        !    read(22,*) ran, lastmodel(1,1,i), ran, ran, ran
        !    lastmodel(:,1,i) = lastmodel(1,1,i)
        !   enddo
        !   close(22)
        !  else
        print *,'Reading in previous model ', itermin
        call readfits(adjustl(trim(directory))//'update/model_c_'//itermin//'.fits', lastmodel_c,nx, 1,  nz)  
        call readfits(adjustl(trim(directory))//'update/model_psi_'//itermin//'.fits', lastmodel_psi,nx, 1,  nz)  
        
        call countzeros(lastmodel_psi,'lastmodel_psi')
        !   lastmodel_psi = 250e-4
        !  endif

        call writefits_3d(adjustl(trim(directory))//'update/update_c_'//itermin//'.fits',totkern_c,nx, 1, nz)
        call writefits_3d(adjustl(trim(directory))//'update/update_psi_'//itermin//'.fits',totkern_psi,nx, 1, nz)

        !  call writefits_3d(adjustl(trim(directory))//'update/update_vx_'//itermin//'.fits',totkern_vx,nx, 1, nz)
        !  call writefits_3d(adjustl(trim(directory))//'update/update_vz_'//itermin//'.fits',totkern_vz,nx, 1, nz)

        !  totkern_c = totkern_c/maxval(abs(totkern_c))
        !  totkern_psi = totkern_psi/maxval(abs(totkern_psi))
        !  totkern_vx = totkern_vx/maxval(abs(totkern_vx))
        !  totkern_vz = totkern_vz/maxval(abs(totkern_vz))

        update_c = totkern_c 
        update_psi = totkern_psi
        
        !  update_vx = totkern_vx
        !  update_vz = totkern_vz
        
        if (iteration > 1) print *,'FORCING STEEPEST DESCENT'

    elseif (iteration .gt. 1 .and. (.not. FORCE_STEEPEST) .and. (.not. LBFGS) ) then 
        print *,'CONJUGATE GRADIENT'
        write(itermin, '(I2.2)') iteration-1
        call readfits(adjustl(trim(directory))//'update/model_c_'//itermin//'.fits', lastmodel_c,nx, 1,  nz)  
        write(itermin, '(I2.2)') iteration-1
        call readfits(adjustl(trim(directory))//'update/model_psi_'//itermin//'.fits', lastmodel_psi,nx, 1,  nz)  

        write(itermin, '(I2.2)') iteration-1
        call readfits(adjustl(trim(directory))//'update/gradient_psi_'//itermin//'.fits', gradpsi,nx, 1,  nz)  

        write(itermin, '(I2.2)') iteration-2
        call readfits(adjustl(trim(directory))//'update/gradient_psi_'//itermin//'.fits', lastgradpsi,nx, 1,  nz)  

        con = sum(gradpsi*(gradpsi - lastgradpsi))
        den = sum(lastgradpsi**2.)
        !stop
        !---
        ! write(itermin, '(I2.2)') iteration-1
        ! call readfits(adjustl(trim(directory))//'update/gradient_vz_'//itermin//'.fits', gradvz,nx, 1,  nz)  

        ! write(itermin, '(I2.2)') iteration-2
        ! call readfits(adjustl(trim(directory))//'update/gradient_vz_'//itermin//'.fits', lastgradvz,nx, 1,  nz)  

        ! con = sum(gradvz*(gradvz - lastgradvz)) + con
        ! den = sum(lastgradvz**2.) + den

        !---

        write(itermin, '(I2.2)') iteration-1
        call readfits(adjustl(trim(directory))//'update/gradient_c_'//itermin//'.fits', gradc,nx, 1,  nz)  

        write(itermin, '(I2.2)') iteration-2
        call readfits(adjustl(trim(directory))//'update/gradient_c_'//itermin//'.fits', lastgradc,nx, 1,  nz)  

        ! UNDO THIS
        ! con = sum(gradc*(gradc - lastgradc)) + con
        ! den = sum(lastgradc**2.) + den

        !---
        con = con/den
        den = con
        if (con .lt. 0)  con = 0
        print *,con, den,'here'

        write(itermin, '(I2.2)') iteration-2
        call readfits(adjustl(trim(directory))//'update/update_c_'//itermin//'.fits', lastupdate,nx, 1,  nz)  
        update_c = totkern_c +  con* lastupdate

        write(itermin, '(I2.2)') iteration-1
        call writefits_3d(adjustl(trim(directory))//'update/update_c_'//itermin//'.fits', update_c,nx, 1,  nz)  

        !---

        write(itermin, '(I2.2)') iteration-2
        call readfits(adjustl(trim(directory))//'update/update_psi_'//itermin//'.fits', lastupdate,nx, 1,  nz)  
        update_psi = totkern_psi +  con*lastupdate 

        write(itermin, '(I2.2)') iteration-1
        call writefits_3d(adjustl(trim(directory))//'update/update_psi_'//itermin//'.fits', update_psi,nx, 1,  nz)  

        !---

        ! write(itermin, '(I2.2)') iteration-2
        ! call readfits(adjustl(trim(directory))//'update/update_vz_'//itermin//'.fits', lastupdate,nx, 1,  nz)  
        ! update_vz = totkern_vz +  con*lastupdate 

        ! write(itermin, '(I2.2)') iteration-1
        ! call writefits_3d(adjustl(trim(directory))//'update/update_vz_'//itermin//'.fits', update_vz,nx, 1,  nz)  


        !---


        ! update_c = update_c/maxval(abs(update_c))
        ! update_psi = update_psi/maxval(abs(update_psi))
        ! update_vx = update_vx/maxval(abs(update_vx))
        ! update_vz = update_vz/maxval(abs(update_vz))
        ! lastupdate = lastupdate/maxval(abs(lastupdate))

        ! update_vx = totkern_vx + sum(totkern_vx*(totkern_vx - lastupdate))/sum(lastupdate**2.) * lastupdate
        ! update_vz = totkern_vz + sum(totkern_vz*(totkern_vz - lastupdate))/sum(lastupdate**2.) * lastupdate

    elseif (iteration .gt. 1 .and. LBFGS) then

        allocate(alph(1:iteration))
        iterati = iteration - 1
        write(chartem, '(I2.2)') iterati
        call readfits(adjustl(trim(directory))//'update/model_psi_'//chartem//'.fits', model,nx, 1,  nz)  

        !  call readfits(adjustl(trim(directory))//'update/gradient_vx_'//chartem//'.fits', gradvx,nx, 1,  nz)  
        !  call readfits(adjustl(trim(directory))//'update/gradient_vz_'//chartem//'.fits', gradvy,nx, 1,  nz)  
        call readfits(adjustl(trim(directory))//'update/gradient_c_'//chartem//'.fits', gradc,nx, 1,  nz)  

        update_c = gradc
        do i=iterati-1, minBFGS, -1

            write(chartem, '(I2.2)') i

            call readfits(adjustl(trim(directory))//'update/model_psi_'//chartem//'.fits', lastmodel,nx, 1,  nz)  

            !  call readfits(adjustl(trim(directory))//'update/gradient_vx_'//chartem//'.fits', lastgradvx,nx, 1,  nz)  
            !  call readfits(adjustl(trim(directory))//'update/gradient_vz_'//chartem//'.fits', lastgradvy,nx, 1,  nz)  
            call readfits(adjustl(trim(directory))//'update/gradient_c_'//chartem//'.fits', lastgradc,nx, 1,  nz)  
            alph(i) = sum((model-lastmodel)*update_c) /sum((model-lastmodel)*(gradc-lastgradc))
            update_c = update_c - alph(i) * (gradc - lastgradc) 

            model = lastmodel
            gradc = lastgradc
        enddo

        write(chartem, '(I2.2)') i
        call readfits(adjustl(trim(directory))//'update/model_psi_'//chartem//'.fits', lastmodel,nx, 1,  nz)  

        !  call readfits(adjustl(trim(directory))//'update/gradient_vx_'//chartem//'.fits', gradvx,nx, 1,  nz)  
        !  call readfits(adjustl(trim(directory))//'update/gradient_vz_'//chartem//'.fits', gradvy,nx, 1,  nz)  
        call readfits(adjustl(trim(directory))//'update/gradient_c_'//chartem//'.fits', lastgradc,nx, 1,  nz)  
        start = i+1
        ! update_c = gradc
        do i=start, iterati-1
            print *,i, i-1 
            write(chartem, '(I2.2)') i

            call readfits(adjustl(trim(directory))//'update/model_psi_'//chartem//'.fits', model,nx, 1,  nz)  

            !  call readfits(adjustl(trim(directory))//'update/gradient_vx_'//chartem//'.fits', lastgradvx,nx, 1,  nz)  
            !  call readfits(adjustl(trim(directory))//'update/gradient_vz_'//chartem//'.fits', lastgradvy,nx, 1,  nz)  
            call readfits(adjustl(trim(directory))//'update/gradient_c_'//chartem//'.fits', gradc,nx, 1,  nz)  
            alph(i) = alph(i) - sum((gradc-lastgradc)*update_c) /sum((model-lastmodel)*(gradc-lastgradc))
            update_c = update_c + alph(i) * (model - lastmodel) 


            !  print *,sum((model-lastmodel)*(gradc-lastgradc)), sum((gradc - lastgradc))
            lastmodel = model
            lastgradc = gradc
        enddo
        !print *,alph!, bet
        ! print *,alph
        ! update_c = update_c/maxval(abs(update_c))


    endif
    call writefits_3d(adjustl(trim(directory))//'update/updatepsi.fits', update,nx, 1,  nz)
    con = maxval(abs(update_psi))
    den = maxval(abs(update_c))
    print *,con,den
    ! if (con .lt. den) con = den
    update_c = update_c/con
    update_psi = update_psi/con
    eps = 0.15

    open(7,file = file_data,form = 'formatted',status = 'old')
    do k = 1,nz
        read(7,*) data(k,:)
    enddo
    close(7)

    do k =1,nz
        z0(k) = (data(k,1)-1.)*695.9895
    enddo
    !update_psi=update_psi*8. 
    !print *,maxval(lastmodel_psi),minval(lastmodel_psi) ,maxval(update_psi)
    do i = 1,5
        write(charnum,'(I1.1)') i

        update = lastmodel_c*(1. + eps * (i) * update_c)
       ! call countzeros(update,'updatec')
        
        call writefits_3d(adjustl(trim(directory))//'update/test_c_'//charnum//'.fits', update,nx, 1,  nz)

        update = lastmodel_psi*(1. + eps * (i) * update_psi)
        !call countzeros(update,'updatepsi')
        
        !  do k=1,nz
        !    update(:,1,k) = (update(:,1,k) -update(1,1,1))/(1. + exp((z0(k) + 3.6)/0.3)) + update(1,1,1) 
        !  enddo
        call writefits_3d(adjustl(trim(directory))//'update/test_psi_'//charnum//'.fits', update,nx, 1,  nz)
    enddo


END PROGRAM CONJGRAD

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

SUBROUTINE readfits(filename,read,dim1,dim2,dim3)

    implicit none
    integer status,unit,readwrite,blocksize,naxes(3), tag, k, dim3 ,dim2,dim1
    integer group,firstpix, ierr
    integer disp, temptyp, sendtyp
    integer nelements
    real*8 nullval,read(dim1,dim2,dim3)
    real*8, dimension(:,:,:), allocatable :: temp
    logical anynull, lexist
    character*(*) filename

    status=0
    call ftgiou(unit,status)
    readwrite=0
    inquire(file=filename, exist=  lexist)
    if (.not. lexist) print *, 'File doesnt exist: '//filename
    if (lexist) then
        print *,'Now reading the file: '//filename
        call ftopen(unit,filename,readwrite,blocksize,status)

        naxes(1) = dim1
        naxes(2) = dim2
        naxes(3) = dim3
        nelements=naxes(1)*naxes(2)*naxes(3)
        group=1
        firstpix=1
        nullval=-999

        call ftgpvd(unit,group,firstpix,nelements,nullval, &
        &            read,anynull,status)


        call ftclos(unit, status)
        call ftfiou(unit, status)

    endif

end SUBROUTINE readfits

!================================================================================
SUBROUTINE writefits_3d(filename, dump_array, dim1, dim2, dim3)

    implicit none
    integer blocksize,bitpix,naxes(3),unit1,dim3, dim1, dim2
    integer status1,group,fpixel, ierr
    integer temptype, sendtyp
    integer*8 nelements
    real*8 dump_array(dim1,dim2,dim3)
    character*(*) filename
    logical simple,extend


    print *,'Writing file '//filename
    call system('rm -rf '//filename)
    status1 = 0
    call ftgiou(unit1,status1)
    blocksize=1
    call ftinit(unit1,filename,blocksize,status1)
    simple=.true.
    bitpix=-64
    naxes(1)=dim1
    naxes(2)=dim2
    naxes(3)=dim3
    nelements=naxes(1)*naxes(2)*naxes(3)
    extend=.false.
    group=1
    fpixel=1
                                                                                                                                                  
    call ftphpr(unit1,simple,bitpix,3,naxes,0,1,extend,status1)
    call ftpprd(unit1,group,fpixel,nelements,dump_array,status1)
    call ftclos(unit1, status1)
    call ftfiou(unit1, status1)


end SUBROUTINE writefits_3d


!================================================================================

SUBROUTINE FILTER_Z(a)

    implicit none
    include 'params.i'
    integer k,kst
    real*8 coeffs(0:5)
    real*8 a(nx,1,nz)
    real*8, allocatable, dimension(:,:,:) :: temp

    coeffs(0) = 0.75390625000000
    coeffs(1) = 0.41015625000000
    coeffs(2) = -0.23437500000000
    coeffs(3) = 0.08789062500000
    coeffs(4) = -0.01953125000000
    coeffs(5) = 0.00195312500000

    allocate(temp(nx,ny,nz))
    kst = 6
    do k=6,nz-5
        temp(:,:,k) = coeffs(0)*a(:,:,k) + 0.5*coeffs(1)*(a(:,:,k-1) + a(:,:,k+1)) &
                    +  0.5*coeffs(2)*(a(:,:,k-2) + a(:,:,k+2))  &
                    +  0.5*coeffs(3)*(a(:,:,k-3) + a(:,:,k+3))  &
                    +  0.5*coeffs(4)*(a(:,:,k-4) + a(:,:,k+4))  &
                    +  0.5*coeffs(5)*(a(:,:,k-5) + a(:,:,k+5))
    enddo

    if (kst==2) then

        temp(:,:,nz-4) = 0.5 * (a(:,:,nz-3) + a(:,:,nz-5))
        temp(:,:,nz-3) = 0.5 * (a(:,:,nz-2) + a(:,:,nz-4))
        temp(:,:,nz-2) = 0.5 * (a(:,:,nz-1) + a(:,:,nz-3))
        temp(:,:,nz-1) = 0.5 * (a(:,:,nz) + a(:,:,nz-2))

        temp(:,:,5) = 0.5 * (a(:,:,6) + a(:,:,4))
        temp(:,:,4) = 0.5 * (a(:,:,5) + a(:,:,3))
        temp(:,:,3) = 0.5 * (a(:,:,4) + a(:,:,2))
        temp(:,:,2) = 0.5 * (a(:,:,3) + a(:,:,1))

    endif

    do k=kst,nz-kst+1
        a(:,:,k) = temp(:,:,k)
    enddo

    a(:,:,1:kst-1) = 0.0
    a(:,:,nz-kst+2:nz) = 0.0

    deallocate(temp)


END SUBROUTINE FILTER_Z 

SUBROUTINE ANTISYMMETRIZE(a)
    implicit none
    include 'params.i'
    real*8 a(nx,1,nz),b(nx,1,nz)
    integer i,j,k
    print *, 'in antisym \n'
    !call countzeros(a,'a')
    do k=1,nz
        do i=1,nx
            b(i,1,k)=a(nx-i+1,1,k)
        end do   
    end do
    !call countzeros(b,'b')
    a=(a-b)*0.5
  
END SUBROUTINE ANTISYMMETRIZE

!================================================================================

subroutine countzeros(a,arrname)
    implicit none
    include 'params.i'
    real*8,intent(in) :: a(nx,1,nz)
    integer :: numzeros,i,k
    character(len=15),intent(in) :: arrname
    
    numzeros=0
    do k=1,nz
        do i=1,nx
            if (a(i,1,k)==0) numzeros=numzeros+1
        end do   
    end do
    
    print *,"Number of zeros in "//arrname,numzeros,"out of",nx*nz
    
end subroutine countzeros
    
