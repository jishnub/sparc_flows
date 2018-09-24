Module all_modules

! --------------------------------------------------------------------------
! MPI Version of the Spherical Acoustic Sun Simulator.
! Copyright 2006, Shravan Hanasoge
                                                                                                                                                        
! Hansen Experimental Physics Laboratory
! 455 Via Palou way, Stanford
! CA 94305, USA
! Email: shravan@stanford.edu
! --------------------------------------------------------------------------
!
! More Subroutines.
!

use initialize
implicit none

Contains
!==================================================================================

function norm2(matrix)
  
   implicit none
   integer i,j,ierr
   real*8 matrix(nx,dim2(rank))
   real*8 norm2, sum

   sum = 0.0  
   norm2  = 0.0
   do j =1,dim2(rank)
    do i =1,nx

      sum = sum + matrix(i,j)**2.0
    end do     
   end do     

   call MPI_REDUCE(sum, norm2, 1, MPI_DOUBLE_PRECISION, &
				MPI_SUM, 0, MPI_COMM_WORLD, ierr)
   

   norm2 = (norm2/(DBLE(nx)*DBLE(ny)))**0.5

end function norm2

!================================================================================

      SUBROUTINE printerror(status)

      ! See cookbook.f on FITSIO for details

      integer status
      character errtext*30,errmessage*80

      if (status .le. 0)return

      call ftgerr(status,errtext)
      print *,'FITSIO Error Status =',status,': ',errtext

      call ftgmsg(errmessage)
      do while (errmessage .ne. ' ')
          print *,errmessage
          call ftgmsg(errmessage)
      end do

      end SUBROUTINE printerror

!================================================================================
      SUBROUTINE deletefile(filename,status)

      integer status,unit,blocksize
      character*(*) filename

      if (status .gt. 0)return
      call ftgiou(unit,status)
      call ftopen(unit,filename,1,blocksize,status)
      

      if (status .eq. 0)then
          call ftdelt(unit,status)
      endif
      if (status .eq. 103)then
          status=0
          call ftcmsg
      endif 
      if ((status .NE. 0) .and. (status .ne. 103)) then 

          status=0
          call ftcmsg
          call ftdelt(unit,status)
      end if

      call ftfiou(unit, status)

      end SUBROUTINE deletefile

!================================================================================

   SUBROUTINE readfits(filename,read,dim3)

      implicit none
      integer status,unit,readwrite,blocksize,naxes(3), tag, k, dim3 
      integer group,firstpix, ierr, stat(MPI_STATUS_SIZE)
      integer disp, temptyp, sendtyp
      integer nelements
      real*8 nullval,read(nx,dim2(rank),dim3)
      real*8, dimension(:,:,:), allocatable :: temp
      logical anynull
      character*(*) filename
      integer (KIND=MPI_ADDRESS_KIND) ext, dum, spacing

      
      if (rank == 0) then
       allocate(temp(nx,ny,dim3))
       status=0
       call ftgiou(unit,status)
       readwrite=0
       print *,'Now reading the file: '//filename
       call ftopen(unit,filename,readwrite,blocksize,status)

       naxes(1) = nx
       naxes(2) = ny
       naxes(3) = dim3
       nelements=naxes(1)*naxes(2)*naxes(3)
       group=1
       firstpix=1
       nullval=-999

       call ftgpvd(unit,group,firstpix,nelements,nullval, &
       &            temp,anynull,status)


       call ftclos(unit, status)
       call ftfiou(unit, status)

       if (status .gt. 0)call printerror(status)

       read = temp(:,1:dim2(0),:)  

       call MPI_TYPE_GET_EXTENT(MPI_DOUBLE_PRECISION, dum, ext, ierr)
       spacing = nx*ny*ext
       disp = dim2(0) + 1
      endif

      call MPI_BARRIER(MPI_COMM_WORLD, ierr)       
      do k=1,numtasks-1
                     
        if (rank ==0) then
 	 call MPI_TYPE_VECTOR(nx*dim2(k),1,1,MPI_DOUBLE_PRECISION,temptyp,ierr)
	 call MPI_TYPE_CREATE_HVECTOR(dim3,1,spacing,temptyp,sendtyp,ierr)
	 call MPI_TYPE_COMMIT(sendtyp, ierr)
	 tag = 0
	 call MPI_SEND(temp(1,disp,1), 1, sendtyp, &
			k, tag, MPI_COMM_WORLD, ierr)
	 disp = disp + dim2(k)
        endif
        if (rank  == k) then
  	 nelements = nx * dim2(k) * dim3
	 tag = 0
	 call MPI_RECV(read, nelements, MPI_DOUBLE_PRECISION,&
			0, tag, MPI_COMM_WORLD, stat, ierr)
        endif

        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
 
      enddo

      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      if (rank == 0) deallocate(temp)

    end SUBROUTINE readfits

!================================================================================


function norm(matrix)

   implicit none
   integer i, j, k, ierr
   real*8 matrix(nx,dim2(rank),nz)
   real*8 norm, sums
  
   norm  = 0.0
   sums = 0.0

   do k = 1,nz
    do j =1,dim2(rank)
     do i =1,nx
       sums = sums + matrix(i,j,k)**2.0
     end do     
    end do     
   end do 
 
   call MPI_REDUCE(sums, norm, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr) 

   norm = (norm/(DBLE(nx)*DBLE(ny)*DBLE(nz)))**0.5d0

   call MPI_BARRIER(MPI_COMM_WORLD, ierr)
end function norm

!================================================================================

function norm_1D(matrix)

   use initialize   
   implicit none
   integer i, j, k, ierr
   real*8 matrix(nz)
   real*8 norm_1d, sum

  
   norm_1d  = 0.0
   sum = 0.0
   do k = 1,nz
     sum = sum + matrix(k)**2.0
   end do 

   norm_1d = (sum/DBLE(nz))**0.5


end function norm_1d

!================================================================================

     function date_time()
  
      implicit none 
      character*27 date_time
      integer i, j, one_day, day, month, hour, minutes
      parameter(one_day = 1440)
      character*2 monthchar, daychar, hourchar, minutechar, secondchar
      character*4 yearchar

      ! Assuming the timestep (in seconds) divides 60 seconds.

      yearchar = '2007'
      monthchar = '06'
      secondchar = '00'
      day = time/(steps*one_day) ! Number of days of computation
      hour = (time/steps - day*1440)/60 ! Remainder hours
      minutes =  (time/steps - day*1440 - hour*60) ! Remainder number of minutes
      day = day + 1 ! Starting day is first day of the month

      call convert_to_string(day, daychar, 2)
      call convert_to_string(hour, hourchar, 2)
      call convert_to_string(minutes, minutechar, 2)

!      date_time ='2015.05.11_'//t1//':'//t2//':'01':00.0_UT'
    
      date_time = yearchar//'.'//monthchar//'.'//daychar//'_'//hourchar//':'//minutechar//':'//secondchar//':00.0_UT'
     end function date_time

!================================================================================
     SUBROUTINE writefits_local(filename, dump_array, dim1, dime2, dim3)

      implicit none
      integer blocksize,bitpix,naxes(3),unit1,dim3
      integer status1,group,fpixel, ierr, dim1, dime2
      integer temptype, sendtyp
      integer*8 nelements
      real*8 dump_array(dim1,dime2,dim3)
      character*(*) filename
      logical simple,extend,lexist


      
      inquire(file=filename,exist=lexist)

      if (lexist) call system('rm '//filename)
      status1 = 0
      call ftgiou(unit1,status1)
      blocksize=1
      call ftinit(unit1,filename,blocksize,status1)
      simple=.true.
      bitpix=-64
      naxes(1)=dim1
      naxes(2)=dime2
      naxes(3)=dim3
      nelements=naxes(1)*naxes(2)*naxes(3)
      extend=.false.
      group=1
      fpixel=1
                                                                                                                                                      
      call ftphpr(unit1,simple,bitpix,3,naxes,0,1,extend,status1)
      call ftpprd(unit1,group,fpixel,nelements,dump_array,status1)
      call ftclos(unit1, status1)
      call ftfiou(unit1, status1)

      if (status1 .gt. 0) call printerror(status1)

!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

     end SUBROUTINE writefits_local
!===============================================================================

     SUBROUTINE writefits_3D(filename,dump_array,dim3)

      implicit none 
      integer blocksize,bitpix,naxes(3),unit1,k,disp, tag, dim3
      integer status1,group,fpixel,ierr,stat(MPI_STATUS_SIZE),flag
      integer temptype, sendtyp, nelements,req
      real*8 dump_array(nx,dim2(rank),dim3)
      real*8, allocatable, dimension(:,:,:) :: temp
      character*(*) filename
      logical simple,extend,lexist
      integer (KIND = MPI_ADDRESS_KIND) space, dum, ext
 
      if (rank == 0) then
	 allocate(temp(nx,ny,dim3))
         temp(:, 1:dim2(0), :) = dump_array
         disp = dim2(0) + 1
         inquire(file=filename,exist=lexist)

         if (lexist) then
          print *,'Deleting file '//filename
	  call system('rm -rf '//filename)
         endif


         print *,'Writing file: ', filename
         call MPI_TYPE_GET_EXTENT(MPI_DOUBLE_PRECISION,dum, ext,ierr)
      endif


       if (rank .ne. 0) then
	 tag = rank
         nelements = nx*dim3*dim2(rank)
         call MPI_ISEND(dump_array, nelements, MPI_DOUBLE_PRECISION, &
			0, tag, MPI_COMM_WORLD, req, ierr)
         !call MPI_TEST(req,flag,stat)
         !print *,flag,rank
       endif

       if (rank == 0) then
          do k=1,numtasks-1
	  space = nx*ny*ext
	  tag = k
	  call MPI_TYPE_VECTOR(nx*dim2(k),1,1,MPI_DOUBLE_PRECISION,temptype,ierr)
	  call MPI_TYPE_CREATE_HVECTOR(dim3,1,space,temptype,sendtyp,ierr)
	  call MPI_TYPE_COMMIT(sendtyp, ierr)
	  call MPI_RECV(temp(1,disp,1), 1, sendtyp, &
			k, tag, MPI_COMM_WORLD, stat, ierr)

	  disp = disp + dim2(k)
         enddo
        endif
 	call MPI_BARRIER(MPI_COMM_WORLD,ierr)

       if (rank == 0) then
         status1 = 0
         call ftgiou(unit1,status1)
         blocksize=1

         call ftinit(unit1,filename,blocksize,status1)
         simple=.true.
         bitpix=-64
         naxes(1)=nx
         naxes(2)=ny
         naxes(3)=dim3
         nelements=naxes(1)*naxes(2)*naxes(3)
         extend=.false.
         group=1
         fpixel=1
         
         call ftphpr(unit1,simple,bitpix,3,naxes,0,1,extend,status1)
         call ftpkyj(unit1,'TIME',time,'Also, init = time+1', status1)
         call ftpprd(unit1,group,fpixel,nelements,temp,status1)
         call ftclos(unit1, status1)
         call ftfiou(unit1, status1)

         deallocate(temp)
         if (status1 .gt. 0) call printerror(status1)

       endif

       call MPI_BARRIER(MPI_COMM_WORLD,ierr)

     end SUBROUTINE writefits_3D
!================================================================================
function sum2(matrix, dima, dimb, option)
  
   implicit none
   integer i,j,ierr, option, dima, dimb
   real*8 matrix(dima, dimb)
   real*8 sum2, summer

   summer = 0.0  
   sum2  = 0.0
   if (option == 0) then
    do j =1,dimb
     do i =1,dima
       summer = summer + matrix(i,j)
     end do     
    end do     
   else
    do j =1,dimb
     do i =1,dima
       summer = summer + abs(matrix(i,j))
     end do     
    end do     
   endif

   call MPI_REDUCE(summer, sum2, 1, MPI_DOUBLE_PRECISION, &
				MPI_SUM, 0, MPI_COMM_WORLD, ierr)
   

   sum2 = (sum2/(DBLE(dima)*DBLE(dimb)))

end function sum2

!================================================================================
function sum2_normal(matrix, dima, dimb, option)
  
   implicit none
   integer i,j,ierr, option, dima, dimb
   real*8 matrix(dima, dimb)
   real*8 sum2_normal, summer

   summer = 0.0  
   sum2_normal  = 0.0
   if (option == 0) then
    do j =1,dimb
     do i =1,dima
       summer = summer + matrix(i,j)
     end do     
    end do     
   else
    do j =1,dimb
     do i =1,dima
       summer = summer + abs(matrix(i,j))
     end do     
    end do     
   endif

   sum2_normal = (summer/(DBLE(dima)*DBLE(dimb)))

end function sum2_normal


!================================================================================

   SUBROUTINE readfits_local(filename,readarr,dim1,dime2,dim3)

      implicit none
      integer status,unit,readwrite,blocksize,naxes(3), tag, k, dim3 ,dime2,dim1
      integer group,firstpix, ierr, stat(MPI_STATUS_SIZE)
      integer disp, temptyp, sendtyp
      integer nelements
      real*8 nullval,readarr(dim1,dime2,dim3)
      real*8, dimension(:,:,:), allocatable :: temp
      logical anynull
      character*(*) filename
      integer (KIND=MPI_ADDRESS_KIND) ext, dum, spacing

      
      status=0
       call ftgiou(unit,status)
       readwrite=0
       print *,'Now reading the file: '//filename
       call ftopen(unit,filename,readwrite,blocksize,status)

       naxes(1) = dim1
       naxes(2) = dime2
       naxes(3) = dim3
       nelements=naxes(1)*naxes(2)*naxes(3)
       group=1
       firstpix=1
       nullval=-999

       call ftgpvd(unit,group,firstpix,nelements,nullval, &
       &            readarr,anynull,status)


       call ftclos(unit, status)
       call ftfiou(unit, status)


       if (status .gt. 0)call printerror(status)



    end SUBROUTINE readfits_local

!================================================================================
  
  SUBROUTINE READ_IN_FIELD ()

   implicit none
   integer j

   call readfits(dirbackmag//'density128.fits',rho0,nz)
   call readfits(dirbackmag//'soundspeed128.fits',c_speed,nz)
   call readfits(dirbackmag//'pressure128.fits',p0,nz)
   call readfits(dirbackmag//'B_x128.fits',box,nz)
   call readfits(dirbackmag//'B_y128.fits',boy,nz)
   call readfits(dirbackmag//'B_z128.fits',boz,nz)

   boz =  boz/dimb
   box = box/dimb
   boy = boy/dimb

   p0 = p0/(dimrho * dimc**2.0)
   rho0 = rho0/dimrho
   c_speed = c_speed/dimc
   c2 = c_speed**2.
  

!   do j=1,nz
!    c2(:,:,j) = gamma(j)*p0(:,:,j)/rho0(:,:,j)
!   enddo
!   c_speed= c2**(0.5)
  
  END SUBROUTINE READ_IN_FIELD
   
!================================================================================
  
  SUBROUTINE READ_IN_FLOWS ()
  
   implicit none

   call readfits(dirbackmag//'v0_x.fits',v0_x,nz)
   call readfits(dirbackmag//'v0_y.fits',v0_y,nz)
   call readfits(dirbackmag//'v0_z.fits',v0_z,nz)

   v0_x = v0_x/dimc
   v0_y = v0_y/dimc
   v0_z = v0_z/dimc

  END SUBROUTINE READ_IN_FLOWS

!================================================================================

  SUBROUTINE READ_IN_INITIAL_CONDITION(director, timest)
 
    implicit none
    character*(timestamp_size) tempc
    character*(*) director
    character*10 ci
    integer timest
  
    call convert_to_string(timest,tempc,timestamp_size)

    if (.not. displ .and. (.not. TEST_IN_2D)) then

     call readfits(director//'p_'//tempc//'_full.fits',a(:,:,:,5),nz)
     call readfits(director//'vy_'//tempc//'_full.fits',a(:,:,:,3),nz)
     call readfits(director//'rho_'//tempc//'_full.fits',a(:,:,:,1),nz)
     call readfits(director//'vx_'//tempc//'_full.fits',a(:,:,:,2),nz)
     call readfits(director//'vz_'//tempc//'_full.fits',a(:,:,:,4),nz)

     a(:,:,:,5) = a(:,:,:,5)/(dimc**2. * dimrho)
     a(:,:,:,1) = a(:,:,:,1)/dimrho
     a(:,:,:,2) = a(:,:,:,2)/dimc
     a(:,:,:,3) = a(:,:,:,3)/dimc
     a(:,:,:,4) = a(:,:,:,4)/dimc

     if (magnetic) then
      call readfits(director//'bx_'//tempc//'_full.fits',a(:,:,:,6),nz)
      call readfits(director//'by_'//tempc//'_full.fits',a(:,:,:,7),nz)
      call readfits(director//'bz_'//tempc//'_full.fits',a(:,:,:,8),nz)

      a(:,:,:,6) = a(:,:,:,6)/dimb
      a(:,:,:,7) = a(:,:,:,7)/dimb
      a(:,:,:,8) = a(:,:,:,8)/dimb

     endif

    elseif (displ) then

     call readfits(director//'xix_'//tempc//'_full.fits',a(:,:,:,1),nz)
     call readfits(director//'xiy_'//tempc//'_full.fits',a(:,:,:,2),nz)
     call readfits(director//'xiz_'//tempc//'_full.fits',a(:,:,:,3),nz)
     call readfits(director//'vx_'//tempc//'_full.fits',a(:,:,:,4),nz)
     call readfits(director//'vy_'//tempc//'_full.fits',a(:,:,:,5),nz)
     call readfits(director//'vz_'//tempc//'_full.fits',a(:,:,:,6),nz)

     a(:,:,:,1) = a(:,:,:,1)/diml
     a(:,:,:,2) = a(:,:,:,2)/diml
     a(:,:,:,3) = a(:,:,:,3)/diml
     a(:,:,:,4) = a(:,:,:,4)/dimc
     a(:,:,:,5) = a(:,:,:,5)/dimc
     a(:,:,:,6) = a(:,:,:,6)/dimc

    endif

    if (USE_PML) then
     call readfits(director//'psivz_'//tempc//'_full.fits',pmlvars(:,:,:,1),nzpml)
     call readfits(director//'psip_'//tempc//'_full.fits',pmlvars(:,:,:,2),nzpml)
     if (magnetic) then         
      call readfits(director//'psidzbx_'//tempc//'_full.fits',pmlvars(:,:,:,3),nzpml)
      call readfits(director//'psidzby_'//tempc//'_full.fits',pmlvars(:,:,:,4),nzpml)
      call readfits(director//'psiinductionbx_'//tempc//'_full.fits',pmlvars(:,:,:,5),nzpml)
      call readfits(director//'psiinductionby_'//tempc//'_full.fits',pmlvars(:,:,:,6),nzpml)
     endif

     if (HORIZONTAL_PMLS) then
      write(ci,'(i10)') rank
      call readfits_local(director//'pmlvarsx_1_'//trim(adjustl(trim(ci)))// &
      '_'//tempc//'_full.fits',pmlvarsx(:,:,:,1),2*npmlhor,dim2(rank),nz)

      call readfits_local(director//'pmlvarsx_2_'//trim(adjustl(trim(ci)))// &
      '_'//tempc//'_full.fits',pmlvarsx(:,:,:,2),2*npmlhor,dim2(rank),nz)


      if (PROC_HAS_PML) then
       call readfits_local(director//'pmlvarsy_1_'//trim(adjustl(trim(ci)))// &
       '_'//tempc//'_full.fits',pmlvarsy(:,:,:,1),nx,2*npmlhor,nz)

       call readfits_local(director//'pmlvarsy_2_'//trim(adjustl(trim(ci)))// &
       '_'//tempc//'_full.fits',pmlvarsy(:,:,:,2),nx,2*npmlhor,nz)

      endif

     endif

    endif

  END SUBROUTINE READ_IN_INITIAL_CONDITION


!================================================================================

  SUBROUTINE WRITE_OUT_FULL_STATE(director)

    implicit none
    integer ierr
    character*(*) director
    character*(timestamp_size) tempc
    character*10 ci
   
    call convert_to_string(time,tempc,timestamp_size)

    if (rank==0) call system('find '//director//' -name "*full*" -delete')
    call mpi_barrier(mpi_comm_world, ierr)

    if (.not. displ) then

     call writefits_3d(director//'p_'//tempc//'_full.fits',a(:,:,:,5)*dimrho*dimc**2.,nz)
     call writefits_3d(director//'vy_'//tempc//'_full.fits',a(:,:,:,3)*dimc,nz)
     call writefits_3d(director//'rho_'//tempc//'_full.fits',a(:,:,:,1)*dimrho,nz)
     call writefits_3d(director//'vx_'//tempc//'_full.fits',a(:,:,:,2)*dimc,nz)
     call writefits_3d(director//'vz_'//tempc//'_full.fits',a(:,:,:,4)*dimc,nz)

     if (magnetic) then
      call writefits_3d(director//'bx_'//tempc//'_full.fits',a(:,:,:,6)*dimb,nz)
      call writefits_3d(director//'by_'//tempc//'_full.fits',a(:,:,:,7)*dimb,nz)
      call writefits_3d(director//'bz_'//tempc//'_full.fits',a(:,:,:,8)*dimb,nz)
     endif

    elseif (displ) then

     call writefits_3d(director//'xix_'//tempc//'_full.fits',a(:,:,:,1)*diml,nz)
     call writefits_3d(director//'xiz_'//tempc//'_full.fits',a(:,:,:,3)*diml,nz)
     call writefits_3d(director//'vx_'//tempc//'_full.fits',a(:,:,:,4)*dimc,nz)
     call writefits_3d(director//'vz_'//tempc//'_full.fits',a(:,:,:,6)*dimc,nz)
     if (.not. test_in_2D) then
      call writefits_3d(director//'xiy_'//tempc//'_full.fits',a(:,:,:,2)*diml,nz)
      call writefits_3d(director//'vy_'//tempc//'_full.fits',a(:,:,:,5)*dimc,nz)
     endif
    endif

    if (USE_PML) then
     call writefits_3d(director//'psivz_'//tempc//'_full.fits',pmlvars(:,:,:,1),nzpml)
     call writefits_3d(director//'psip_'//tempc//'_full.fits',pmlvars(:,:,:,2),nzpml)
     if (magnetic) then         
      call writefits_3d(director//'psidzbx_'//tempc//'_full.fits',pmlvars(:,:,:,3),nzpml)
      call writefits_3d(director//'psidzby_'//tempc//'_full.fits',pmlvars(:,:,:,4),nzpml)
      call writefits_3d(director//'psiinductionbx_'//tempc//'_full.fits',pmlvars(:,:,:,5),nzpml)
      call writefits_3d(director//'psiinductionby_'//tempc//'_full.fits',pmlvars(:,:,:,6),nzpml)
     endif

     if (HORIZONTAL_PMLS) then
      write(ci,'(i10)') rank
      call writefits_local(director//'pmlvarsx_1_'//trim(adjustl(trim(ci)))// &
      '_'//tempc//'_full.fits',pmlvarsx(:,:,:,1),2*npmlhor,dim2(rank),nz)

      call writefits_local(director//'pmlvarsx_2_'//trim(adjustl(trim(ci)))// &
      '_'//tempc//'_full.fits',pmlvarsx(:,:,:,2),2*npmlhor,dim2(rank),nz)


      if (PROC_HAS_PML) then
       call writefits_local(director//'pmlvarsy_1_'//trim(adjustl(trim(ci)))// &
       '_'//tempc//'_full.fits',pmlvarsy(:,:,:,1),nx,dim2(rank),nz)

       call writefits_local(director//'pmlvarsy_2_'//trim(adjustl(trim(ci)))// &
       '_'//tempc//'_full.fits',pmlvarsy(:,:,:,2),nx,dim2(rank),nz)

      endif

     endif
    endif


  END SUBROUTINE WRITE_OUT_FULL_STATE


!================================================================================

  SUBROUTINE WRITE_OUT_SLICE(direct)

    implicit none
    integer radial
    real*8 lnorm, cnorm
    character*(timestamp_size) tempc
    character*(*) direct
  
    radial = o_rad
    call convert_to_string(time,tempc,timestamp_size)

    if (.not. displ .and. (.not. TEST_IN_2D)) then

     call writefits_3d(direct//'p_'//tempc//'_slice.fits',a(:,:,radial,5)*dimrho*dimc**2.,1)
     call writefits_3d(direct//'vy_'//tempc//'_slice.fits',a(:,:,radial,3)*dimc,1)
     call writefits_3d(direct//'rho_'//tempc//'_slice.fits',a(:,:,radial,1)*dimrho,1)
     call writefits_3d(direct//'vx_'//tempc//'_slice.fits',a(:,:,radial,2)*dimc,1)
     call writefits_3d(direct//'vz_'//tempc//'_slice.fits',a(:,:,radial,4)*dimc,1)

     if (magnetic) then
      call writefits_3d(direct//'bx_'//tempc//'_slice.fits',a(:,:,radial,6)*dimb,1)
      call writefits_3d(direct//'by_'//tempc//'_slice.fits',a(:,:,radial,7)*dimb,1)
      call writefits_3d(direct//'bz_'//tempc//'_slice.fits',a(:,:,radial,8)*dimb,1)
     endif

    elseif (displ .and. (.not. TEST_IN_2D)) then

     if (kernel_mode) then
      lnorm = 1.
      cnorm = 1.
     else
      lnorm = diml
      cnorm = dimc
     endif

     call writefits_3d(direct//'xix_'//tempc//'_slice.fits',a(:,:,radial,1)*lnorm,1)
     call writefits_3d(direct//'xiy_'//tempc//'_slice.fits',a(:,:,radial,2)*lnorm,1)
     call writefits_3d(direct//'xiz_'//tempc//'_slice.fits',a(:,:,radial,3)*lnorm,1)
     call writefits_3d(direct//'vx_'//tempc//'_slice.fits',a(:,:,radial,4)*cnorm,1)
     call writefits_3d(direct//'vy_'//tempc//'_slice.fits',a(:,:,radial,5)*cnorm,1)
     call writefits_3d(direct//'vz_'//tempc//'_slice.fits',a(:,:,radial,6)*cnorm,1)

    endif


  END SUBROUTINE WRITE_OUT_SLICE


!================================================================================


  SUBROUTINE WRITE_OUT_PARTIAL_STATE(direc)

    implicit none
    character*(timestamp_size) tempc
    character*(timestamp_size-1) tempct
    character*(*) direc

     call convert_to_string(time,tempc,timestamp_size)
     !if (local_time .lt. 0) tempc = '-'//tempct
     !if (local_time .ge. 0) tempc = '+'//tempct

     if (compute_forward) then
      call writefits_3d(direc//'xix_'//tempc//'_partial.fits',a(:,:,st_z:fi_z,1),nz_kern)
      call writefits_3d(direc//'xiz_'//tempc//'_partial.fits',a(:,:,st_z:fi_z,3),nz_kern)


      if (.not. test_in_2d) then
       call writefits_3d(direc//'xiy_'//tempc//'_partial.fits',a(:,:,st_z:fi_z,2),nz_kern)
      endif
     endif 

     call writefits_3d(direc//'acc_x_'//tempc//'_partial.fits',scr(:,:,st_z:fi_z,4),nz_kern)
     call writefits_3d(direc//'acc_z_'//tempc//'_partial.fits',scr(:,:,st_z:fi_z,6),nz_kern)

     call writefits_3d(direc//'vx_'//tempc//'_partial.fits',a(:,:,st_z:fi_z,4),nz_kern)
     call writefits_3d(direc//'vz_'//tempc//'_partial.fits',a(:,:,st_z:fi_z,6),nz_kern)

     if (.not. test_in_2d) then
      call writefits_3d(direc//'acc_y_'//tempc//'_partial.fits',scr(:,:,st_z:fi_z,5),nz_kern)
      call writefits_3d(direc//'vy_'//tempc//'_partial.fits',a(:,:,st_z:fi_z,5),nz_kern)
     endif

     if (rank==0) print *,'Written out partial state; ONLY xi variables for FORWARD'
     if (rank==0) print *,'Written out partial state; ONLY v variables for ADJOINT'
  END SUBROUTINE WRITE_OUT_PARTIAL_STATE


!================================================================================

  SUBROUTINE WRITE_PARTIAL_STATE_BINARIES(direc)

    implicit none
    character*(timestamp_size) tempc
    character*(*) direc
  
    call convert_to_string(time,tempc,timestamp_size)

    call write_binary(direc//'xix_'//tempc//'_partial.dat',a(:,:,st_z:fi_z,1),nz_kern)
    call write_binary(direc//'xiy_'//tempc//'_partial.dat',a(:,:,st_z:fi_z,2),nz_kern)
    call write_binary(direc//'xiz_'//tempc//'_partial.dat',a(:,:,st_z:fi_z,3),nz_kern)
    call write_binary(direc//'vx_'//tempc//'_partial.dat',a(:,:,st_z:fi_z,4),nz_kern)
    call write_binary(direc//'vy_'//tempc//'_partial.dat',a(:,:,st_z:fi_z,5),nz_kern)
    call write_binary(direc//'vz_'//tempc//'_partial.dat',a(:,:,st_z:fi_z,6),nz_kern)


     if (rank==0) print *,'Written out partial binary state'
  END SUBROUTINE WRITE_PARTIAL_STATE_BINARIES


!================================================================================

     SUBROUTINE write_binary(filename,dump_array,dim3)

      implicit none 
      integer k,disp, tag, dim3, reclmax
      integer ierr,stat(MPI_STATUS_SIZE),flag, i, j
      integer temptype, sendtyp, nelements,req
      real*8 dump_array(nx,dim2(rank),dim3)
      real*8, allocatable, dimension(:,:,:) :: temp
      character*(*) filename
      logical simple,extend
      integer (KIND = MPI_ADDRESS_KIND) space, dum, ext
 
      if (rank == 0) then
	 allocate(temp(nx,ny,dim3))
         temp(:, 1:dim2(0), :) = dump_array
         disp = dim2(0) + 1
         print *,'Writing file: ', filename
         call MPI_TYPE_GET_EXTENT(MPI_DOUBLE_PRECISION,dum, ext,ierr)
         reclmax= nx * ny * 2
      endif


       if (rank .ne. 0) then
	 tag = rank
         nelements = nx*dim3*dim2(rank)
         call MPI_ISEND(dump_array, nelements, MPI_DOUBLE_PRECISION, &
			0, tag, MPI_COMM_WORLD, req, ierr)
         !call MPI_TEST(req,flag,stat)
         !print *,flag,rank
       endif

       if (rank == 0) then
          do k=1,numtasks-1
	  space = nx*ny*ext
	  tag = k
	  call MPI_TYPE_VECTOR(nx*dim2(k),1,1,MPI_DOUBLE_PRECISION,temptype,ierr)
	  call MPI_TYPE_CREATE_HVECTOR(dim3,1,space,temptype,sendtyp,ierr)
	  call MPI_TYPE_COMMIT(sendtyp, ierr)
	  call MPI_RECV(temp(1,disp,1), 1, sendtyp, &
			k, tag, MPI_COMM_WORLD, stat, ierr)

	  disp = disp + dim2(k)
         enddo
        endif
 	call MPI_BARRIER(MPI_COMM_WORLD,ierr)

       if (rank == 0) then

        open(44,file=filename,form='unformatted',&
        status='unknown', action='write',&
	access='direct',recl=reclmax)!,RECORDTYPE='fixed')

        do k=1,dim3
 	 write(44, rec=k) ((temp(i,j,k),i=1,nx),j=1,ny)
 	enddo

 	close(44)

        deallocate(temp)

       endif

       call MPI_BARRIER(MPI_COMM_WORLD,ierr)

     end SUBROUTINE write_binary

!================================================================================

   SUBROUTINE read_binary(filename,readarr,dim3)

      implicit none
      integer status,unit,readwrite,blocksize,naxes(3), tag, k, dim3 
      integer group,firstpix, ierr, stat(MPI_STATUS_SIZE), reclmax, i, j
      integer disp, temptyp, sendtyp
      integer nelements
      real*8 nullval,readarr(nx,dim2(rank),dim3)
      real*8, dimension(:,:,:), allocatable :: temp
      logical anynull
      character*(*) filename
      integer (KIND=MPI_ADDRESS_KIND) ext, dum, spacing

      
      if (rank == 0) then

       allocate(temp(nx,ny, dim3))

       reclmax = nx * ny * 2
       print *,'Now reading the file: '//filename
       open(344,file=filename,form='unformatted',status='old', action='read',&
 	access='direct',recl=reclmax)!,recordtype='fixed'

       do k=1,dim3
        read(344, rec=k) ((temp(i,j,k),i=1,nx),j=1,ny)
       enddo

       close(344)


       readarr = temp(:,1:dim2(0),:)  

       call MPI_TYPE_GET_EXTENT(MPI_DOUBLE_PRECISION, dum, ext, ierr)
       spacing = nx*ny*ext
       disp = dim2(0) + 1
      endif

      call MPI_BARRIER(MPI_COMM_WORLD, ierr)       
      do k=1,numtasks-1
                     
        if (rank ==0) then
 	 call MPI_TYPE_VECTOR(nx*dim2(k),1,1,MPI_DOUBLE_PRECISION,temptyp,ierr)
	 call MPI_TYPE_CREATE_HVECTOR(dim3,1,spacing,temptyp,sendtyp,ierr)
	 call MPI_TYPE_COMMIT(sendtyp, ierr)
	 tag = 0
	 call MPI_SEND(temp(1,disp,1), 1, sendtyp, &
			k, tag, MPI_COMM_WORLD, ierr)
	 disp = disp + dim2(k)
        endif
        if (rank  == k) then
  	 nelements = nx * dim2(k) * dim3
	 tag = 0
	 call MPI_RECV(readarr, nelements, MPI_DOUBLE_PRECISION,&
			0, tag, MPI_COMM_WORLD, stat, ierr)
        endif

        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
 
      enddo

      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      if (rank == 0) deallocate(temp)

    end SUBROUTINE read_binary


!================================================================================

     SUBROUTINE write_binary_seq_2D(unitnum,dump_array,dim3, recnum)

      implicit none 
      integer k,disp, tag, dim3, reclmax, unitnum, recnum
      integer ierr,stat(MPI_STATUS_SIZE),flag, i, j
      integer temptype, sendtyp, nelements,req
      real*8 dump_array(nx,dim2(rank),dim3)
      real*8, allocatable, dimension(:,:,:) :: temp
      logical simple,extend
      integer (KIND = MPI_ADDRESS_KIND) space, dum, ext
 
      if (rank == 0) then
	 allocate(temp(nx,ny,dim3))
         temp(:, 1:dim2(0), :) = dump_array
         disp = dim2(0) + 1
         call MPI_TYPE_GET_EXTENT(MPI_DOUBLE_PRECISION,dum, ext,ierr)
         reclmax= nx * ny * 2
      endif


       if (rank .ne. 0) then
	 tag = rank
         nelements = nx*dim3*dim2(rank)
         call MPI_ISEND(dump_array, nelements, MPI_DOUBLE_PRECISION, &
			0, tag, MPI_COMM_WORLD, req, ierr)
         !call MPI_TEST(req,flag,stat)
         !print *,flag,rank
       endif

       if (rank == 0) then
          do k=1,numtasks-1
	  space = nx*ny*ext
	  tag = k
	  call MPI_TYPE_VECTOR(nx*dim2(k),1,1,MPI_DOUBLE_PRECISION,temptype,ierr)
	  call MPI_TYPE_CREATE_HVECTOR(dim3,1,space,temptype,sendtyp,ierr)
	  call MPI_TYPE_COMMIT(sendtyp, ierr)
	  call MPI_RECV(temp(1,disp,1), 1, sendtyp, &
			k, tag, MPI_COMM_WORLD, stat, ierr)

	  disp = disp + dim2(k)
         enddo
        endif
 	call MPI_BARRIER(MPI_COMM_WORLD,ierr)

       if (rank == 0) then


        do k=1,dim3
 	 write(unitnum, rec=(recnum+k-1)) ((temp(i,j,k),i=1,nx),j=1,ny)
 	enddo

! 	close(44)

        deallocate(temp)

       endif

       call MPI_BARRIER(MPI_COMM_WORLD,ierr)

     end SUBROUTINE write_binary_seq_2D

!================================================================================


   SUBROUTINE read_binary_reverse(filename,readarr,dim3)

      implicit none
      integer status,unit,readwrite,blocksize,naxes(3), tag, k, dim3 
      integer group,firstpix, ierr, stat(MPI_STATUS_SIZE), reclmax, i, j
      integer disp, temptyp, sendtyp
      integer nelements
      real*8 nullval,readarr(nx,dim2(rank),dim3)
      real*8, dimension(:,:,:), allocatable :: temp
      logical anynull
      character*(*) filename
      integer (KIND=MPI_ADDRESS_KIND) ext, dum, spacing

      
      if (rank == 0) then

       allocate(temp(nx,ny, dim3))

       reclmax = nx * ny * 2
       print *,'Now reading the file: '//filename
       open(344,file=filename,form='unformatted',status='old', action='read',&
 	access='direct',recl=reclmax)!,recordtype='fixed'

       do k=1,dim3
        read(344, rec=k) ((temp(i,j,dim3-k+1),i=1,nx),j=1,ny)
       enddo

       close(344)


       readarr = temp(:,1:dim2(0),:)  

       call MPI_TYPE_GET_EXTENT(MPI_DOUBLE_PRECISION, dum, ext, ierr)
       spacing = nx*ny*ext
       disp = dim2(0) + 1
      endif

      call MPI_BARRIER(MPI_COMM_WORLD, ierr)       
      do k=1,numtasks-1
                     
        if (rank ==0) then
 	 call MPI_TYPE_VECTOR(nx*dim2(k),1,1,MPI_DOUBLE_PRECISION,temptyp,ierr)
	 call MPI_TYPE_CREATE_HVECTOR(dim3,1,spacing,temptyp,sendtyp,ierr)
	 call MPI_TYPE_COMMIT(sendtyp, ierr)
	 tag = 0
	 call MPI_SEND(temp(1,disp,1), 1, sendtyp, &
			k, tag, MPI_COMM_WORLD, ierr)
	 disp = disp + dim2(k)
        endif
        if (rank  == k) then
  	 nelements = nx * dim2(k) * dim3
	 tag = 0
	 call MPI_RECV(readarr, nelements, MPI_DOUBLE_PRECISION,&
			0, tag, MPI_COMM_WORLD, stat, ierr)
        endif

        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
 
      enddo

      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      if (rank == 0) deallocate(temp)

    end SUBROUTINE read_binary_reverse


!================================================================================



   SUBROUTINE read_binary_writefits(filename, dim3, output_filename)

      implicit none
      integer reclmax, i, j, k, dim3
      real*8 readarr(nx,ny,dim3), tempoa(nx,ny)
      character*(*) filename, output_filename
      inquire(iolength=reclmax) tempoa

       print *,'Now reading the file: ', filename, output_filename
       open(344,file=filename,form='unformatted',status='old', action='read',&
 	access='direct',recl=reclmax)!,recordtype='fixed'

       do k=1,dim3
        read(344, rec=k) ((readarr(i,j,k),i=1,nx),j=1,ny)
       enddo

       close(344)


     if(rank==0) call writefits_local(output_filename, readarr, nx, ny, dim3)
    end SUBROUTINE read_binary_writefits


!================================================================================


end module all_modules

