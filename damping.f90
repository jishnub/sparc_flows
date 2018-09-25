Module Damping

 use initialize
 use all_modules
 implicit none
 
 Contains
!==================================================================================

SUBROUTINE INIT_DAMPING ()

 implicit none
 integer i,j,k
 integer narr(2), inembed(2), onembed(2)
 real(kind=real64) constx, consty, kayref

 allocate(transfermatrix(nx,ny,dimz(rank)), transfertemp(nx/2+1,ny,dimz(rank)),&
    damping_rates(nx/2+1,ny),kay2d(nx/2+1,ny))

 
    constx = 2.*pi/xlength * Rsun
    consty = 2.*pi/ylength * Rsun

     do j=1,ny/2+1
      do i=1,nx/2+1
       kay2d(i,j) = (constx**2.*(i-1.0)**2. + consty**2.*(j-1.0)**2.)**0.5
      enddo
     enddo
   
     do j=ny/2+2,ny
      do i=1,nx/2+1
       kay2d(i,j) = (constx**2.*(i-1.0)**2. + consty**2.*(ny-j+1.0)**2.)**0.5
      enddo
     enddo

     kayref = 902.0
     
     damping_rates = diml/dimc * 10.0**(-6.) * 150.0 * (kay2d/kayref)**2.2 * normx * normy 

      narr(1) = nx
      narr(2) = ny
      inembed(1) = nx
      inembed(2) = ny
      onembed(1) = nx/2+1
      onembed(2) = ny
      call dfftw_plan_many_dft_r2c(fftw_plan_fwd_2D,2,narr,dimz(rank),transfermatrix(1,1,1), &
  & inembed,1,nx*ny,transfertemp(1,1,1),onembed,1,(nx/2+1)*ny,FFTW_MEASURE)


      narr(1) = nx
      narr(2) = ny
      inembed(1) = nx/2+1
      inembed(2) = ny
      onembed(1) = nx
      onembed(2) = ny

      call dfftw_plan_many_dft_c2r(fftw_plan_inv_2D,2,narr,dimz(rank),transfertemp(1,1,1), &
  & inembed,1,(nx/2+1)*ny,transfermatrix(1,1,1),onembed,1,nx*ny,FFTW_MEASURE)

      transfermatrix = 0.
      transfertemp = 0.

      

END SUBROUTINE INIT_DAMPING


!==================================================================================


SUBROUTINE INIT_DAMPING_2D ()

 implicit none
 integer i,j,k
 integer narr(2), inembed(2), onembed(2)
 real(kind=real64) constx, consty, kayref

 allocate(transfermatrix(nx,ny,nz), transfertemp(nx/2+1,ny,nz),&
    damping_rates(nx/2+1,ny),kay2d(nx/2+1,ny))

 
    constx = 2.*pi/xlength * Rsun
    consty = 2.*pi/ylength * Rsun

    j = 1
      do i=1,nx/2+1
       kay2d(i,j) = (constx**2.*(i-1.0)**2. + consty**2.*(j-1.0)**2.)**0.5
      enddo
   
  
     kayref = 902.0
     
     damping_rates = diml/dimc * 10.0**(-6.) * 150.0 * (kay2d/kayref)**2.2 * normx * normy 

      narr(1) = nx
      narr(2) = ny
      inembed(1) = nx
      inembed(2) = ny
      onembed(1) = nx/2+1
      onembed(2) = ny
      call dfftw_plan_many_dft_r2c(fftw_plan_fwd_2D,2,narr,nz,transfermatrix(1,1,1), &
  & inembed,1,nx*ny,transfertemp(1,1,1),onembed,1,(nx/2+1)*ny,FFTW_MEASURE)


      narr(1) = nx
      narr(2) = ny
      inembed(1) = nx/2+1
      inembed(2) = ny
      onembed(1) = nx
      onembed(2) = ny

      call dfftw_plan_many_dft_c2r(fftw_plan_inv_2D,2,narr,nz,transfertemp(1,1,1), &
  & inembed,1,(nx/2+1)*ny,transfermatrix(1,1,1),onembed,1,nx*ny,FFTW_MEASURE)

      transfermatrix = 0.
      transfertemp = 0.

      

END SUBROUTINE INIT_DAMPING_2D


!==================================================================================

SUBROUTINE damp_velocity ()

  implicit none
  integer i,j,k
  real(kind=real64), allocatable, dimension(:,:,:) :: term


  allocate(term(nx,dim2(rank),nz))
  if (.not. test_in_2D) then 

   call TRANSFER_3D_Z(v_x, transfermatrix)
   call dfftw_execute_dft_r2c(fftw_plan_fwd_2D, transfermatrix, transfertemp)
   do k=1,dimz(rank)
    transfertemp(:,:,k) = transfertemp(:,:,k) * damping_rates
   enddo
   call dfftw_execute_dft_c2r(fftw_plan_inv_2D, transfertemp, transfermatrix)
   call INV_TRANSFER_3D_Z(transfermatrix,term)
   RHSv_x = RHSv_x - term

   call TRANSFER_3D_Z(v_y, transfermatrix)
   call dfftw_execute_dft_r2c(fftw_plan_fwd_2D, transfermatrix, transfertemp)
   do k=1,dimz(rank)
    transfertemp(:,:,k) = transfertemp(:,:,k) * damping_rates
   enddo
   call dfftw_execute_dft_c2r(fftw_plan_inv_2D, transfertemp, transfermatrix)
   call INV_TRANSFER_3D_Z(transfermatrix,term)
   RHSv_y = RHSv_y - term

   call TRANSFER_3D_Z(v_z, transfermatrix)
   call dfftw_execute_dft_r2c(fftw_plan_fwd_2D, transfermatrix, transfertemp)
   do k=1,dimz(rank)
    transfertemp(:,:,k) = transfertemp(:,:,k) * damping_rates
   enddo
   call dfftw_execute_dft_c2r(fftw_plan_inv_2D, transfertemp, transfermatrix)
   call INV_TRANSFER_3D_Z(transfermatrix,term)
   RHSv_z = RHSv_z - term


  else

   call dfftw_execute_dft_r2c(fftw_plan_fwd_2D, v_x, transfertemp)
   do k=1,dimz(rank)
    transfertemp(:,:,k) = transfertemp(:,:,k) * damping_rates
   enddo
   call dfftw_execute_dft_c2r(fftw_plan_inv_2D, transfertemp, term)
   RHSv_x = RHSv_x - term

   call dfftw_execute_dft_r2c(fftw_plan_fwd_2D, v_z, transfertemp)
   do k=1,nz
    transfertemp(:,:,k) = transfertemp(:,:,k) * damping_rates
   enddo
   call dfftw_execute_dft_c2r(fftw_plan_inv_2D, transfertemp, term)
   RHSv_z = RHSv_z - term


  endif 


  deallocate(term)

END SUBROUTINE damp_velocity

!==================================================================================

SUBROUTINE TRANSFER_3D_Z_N(input, output)

 ! Splits up the z-direction and accumulates all the y-parts together
 implicit none
 integer i, j, k, sendtag, recvtag, req(2*numtasks-2), ierr, stat(MPI_STATUS_SIZE, 2*numtasks-2)
 integer(kind=8) nelem_send,nelem_recv
 real(kind=real64) input(nx, dim2(rank), nz), output(nx, ny, dimz(rank))
 real(kind=real64), allocatable, dimension(:,:,:) :: temp

 !  Non communicative part:
! output(:,fwd_rdispls(rank):(fwd_rdispls(rank)+dim2(rank)-1),:) = &
 !                        input(:,:,fwd_z_displs(rank):(fwd_z_displs(rank) + dimz(rank)-1))

 do k=1,dimz(rank)
  do j=1,dim2(rank)
   do i=1,nx
    output(i,j+fwd_rdispls(rank)-1,k) = input(i,j,k+fwd_z_displs(rank)-1)
   enddo
  enddo
 enddo

 ! Non-blocking send-recv transpose
 do k=1,numtasks-1
   i = modulo(rank+k,numtasks)
   j = modulo(rank-k,numtasks)

   sendtag = rank
   recvtag = j

   nelem_send = nx*dim2(rank)*dimz(i)
   nelem_recv = nx*dim2(j)*dimz(rank)
!   allocate(temp(nx,dim2(j),dimz(rank)))


   call MPI_ISEND(input(1,1,fwd_z_displs(i)),nelem_send,MPI_DOUBLE_PRECISION,i,sendtag,MPI_COMM_WORLD,req(2*k-1),ierr)
   call MPI_IRECV(output(1,fwd_rdispls(j),1),nelem_recv,MPI_DOUBLE_PRECISION,j,recvtag,MPI_COMM_WORLD,req(2*k),ierr)

 !  output(iin,fwd_rdispls(j):(fwd_rdispls(j)+dim2(j)-1),:) = temp(iin,jin,kin)

 enddo

 call MPI_WAITALL(2*numtasks-2,req,stat,ierr)

END SUBROUTINE TRANSFER_3D_Z_N

!================================================================================================ 

SUBROUTINE INV_TRANSFER_3D_Z_N(input, output)

 ! Splits up the z-direction and accumulates all the y-parts together
 implicit none 
 integer i, j, k, sendtag, recvtag, req(2*numtasks-2), ierr, stat(MPI_STATUS_SIZE, 2*numtasks-2)
 integer iin, jin, kin
 integer(kind=8) nelem_send,nelem_recv
 real(kind=real64) input(nx, ny, dimz(rank)), output(nx, dim2(rank), nz)
 real(kind=real64), dimension(nx, ny, dimz(rank)) :: temp
 
 ! Non-communicative part:

 !    output(:,:,fwd_z_displs(rank):(fwd_z_displs(rank) + dimz(rank)-1)) = &
  !                       input(:,fwd_rdispls(rank):(fwd_rdispls(rank)+dim2(rank)-1),:)

 do k=1,dimz(rank)
  do j=1,dim2(rank)
   do i=1,nx
    output(i,j,fwd_z_displs(rank)+k-1) = &
                         input(i,fwd_rdispls(rank)+j-1,k)
   enddo
  enddo
 enddo
 
!   allocate(temp(nx,dim2(i),dimz(rank)))

 do k=1,numtasks-1

   i = modulo(rank+k,numtasks)
   j = modulo(rank-k,numtasks)

   sendtag = rank
   recvtag = j

   nelem_send = nx*dim2(i)*dimz(rank)
   nelem_recv = nx*dim2(rank)*dimz(j)

!   do kin=1,dimz(rank)
!    do jin=fwd_sdispls(i),dim2(i)+fwd_sdispls(i)
!     do iin=1,nx
!      temp(iin,jin,kin) = input(iin,jin,kin)
!     enddo
!    enddo
!   enddo

   call MPI_ISEND(input(1,fwd_sdispls(i),1),nelem_send,MPI_DOUBLE_PRECISION,i,sendtag,MPI_COMM_WORLD,req(2*k-1),ierr)
   call MPI_IRECV(output(1,1,fwd_z_displs(j)),nelem_recv,MPI_DOUBLE_PRECISION,j,recvtag,MPI_COMM_WORLD,req(2*k),ierr)
 enddo
   
 call MPI_WAITALL(2*numtasks-2,req,stat,ierr)


!   temp = input(:,fwd_sdispls(i):(fwd_sdispls(i)+dim2(i)-1),:)

END SUBROUTINE INV_TRANSFER_3D_Z_N
   
!================================================================================================ 

SUBROUTINE TRANSFER_3D_Z(input, output)

 ! Splits up the z-direction and accumulates all the y-parts together
 implicit none
 integer i, j, k, sendtag, recvtag, req(2), ierr, stat(MPI_STATUS_SIZE, 2)
 integer(kind=8) nelem_send,nelem_recv
 real(kind=real64) input(nx, dim2(rank), nz), output(nx, ny, dimz(rank))
 real(kind=real64), allocatable, dimension(:,:,:) :: temp

 !  Non communicative part:
 output(:,fwd_rdispls(rank):(fwd_rdispls(rank)+dim2(rank)-1),:) = &
                         input(:,:,fwd_z_displs(rank):(fwd_z_displs(rank) + dimz(rank)-1))

 ! Non-blocking send-recv transpose
 do k=1,numtasks-1
   i = modulo(rank+k,numtasks)
   j = modulo(rank-k,numtasks)

   sendtag = 0
   recvtag = 0
   nelem_send = nx*dim2(rank)*dimz(i)
   nelem_recv = nx*dim2(j)*dimz(rank)
   allocate(temp(nx,dim2(j),dimz(rank)))

   call MPI_ISEND(input(1,1,fwd_z_displs(i)),nelem_send,MPI_DOUBLE_PRECISION,i,sendtag,MPI_COMM_WORLD,req(1),ierr)
   call MPI_IRECV(temp,nelem_recv,MPI_DOUBLE_PRECISION,j,recvtag,MPI_COMM_WORLD,req(2),ierr)
   call MPI_WAITALL(2,req,stat,ierr)

   output(:,fwd_rdispls(j):(fwd_rdispls(j)+dim2(j)-1),:) = temp
   deallocate(temp)
 enddo

END SUBROUTINE TRANSFER_3D_Z

!================================================================================================ 

SUBROUTINE INV_TRANSFER_3D_Z(input, output)

 ! Splits up the z-direction and accumulates all the y-parts together
 implicit none 
 integer i, j, k, sendtag, recvtag, req(2), ierr, stat(MPI_STATUS_SIZE, 2)
 integer(kind=8) nelem_send,nelem_recv
 real(kind=real64) input(nx, ny, dimz(rank)), output(nx, dim2(rank), nz)
 real(kind=real64), allocatable,dimension(:,:,:) :: temp
 
 ! Non-communicative part:
 output(:,:,fwd_z_displs(rank):(fwd_z_displs(rank) + dimz(rank)-1)) = &
                         input(:,fwd_rdispls(rank):(fwd_rdispls(rank)+dim2(rank)-1),:)
 
 do k=1,numtasks-1
   i = modulo(rank+k,numtasks)
   j = modulo(rank-k,numtasks)
   sendtag = 0
   recvtag = 0
   nelem_send = nx*dim2(i)*dimz(rank)
   nelem_recv = nx*dim2(rank)*dimz(j)
  
   allocate(temp(nx,dim2(i),dimz(rank)))
   temp = input(:,fwd_sdispls(i):(fwd_sdispls(i)+dim2(i)-1),:)
  
   call MPI_ISEND(temp,nelem_send,MPI_DOUBLE_PRECISION,i,sendtag,MPI_COMM_WORLD,req(1),ierr)
   call MPI_IRECV(output(1,1,fwd_z_displs(j)),nelem_recv,MPI_DOUBLE_PRECISION,j,recvtag,MPI_COMM_WORLD,req(2),ierr)
   call MPI_WAITALL(2,req,stat,ierr)
  
   deallocate(temp)
 enddo
   

END SUBROUTINE INV_TRANSFER_3D_Z
   
!================================================================================================ 


END MODULE DAMPING
