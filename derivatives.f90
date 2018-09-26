MODULE DERIVATIVES

 use initialize
 implicit none

Contains

!================================================================================================ 

SUBROUTINE ddz(var, dvar, bc)
 implicit none
 real(kind=real64) var(nx,dim2(rank),nz), dvar(nx,dim2(rank),nz)
 integer i,j,k,bc

 if (compact_finite_diff) then
  call dbyd2(dvar,var,nx*dim2(rank),nz,bc)
  do k=1,nz
   dvar(:,:,k) = dvar(:,:,k)*stretch(k)
  enddo
 else
  call ddz_ninthorder(var,dvar)
 endif

END SUBROUTINE ddz

!================================================================================================ 

SUBROUTINE ddx(var, dvar,bc)
 implicit none
 real(kind=real64) var(nx,dim2(rank),nz), dvar(nx,dim2(rank),nz)
 integer i,j,k,bc
 complex(kind=real64), allocatable, dimension(:,:,:) :: temp

 if ((PERIODIC) .AND. (USE_FFT)) then
  allocate(temp(nx/2+1,dim2(rank), nz))
  call dfftw_execute_dft_r2c(fftw_plan_fwd_x, var, temp)

  do k=1,nz
   do j=1,dim2(rank)
    do i=1,nx/2+1
     temp(i,j,k) = temp(i,j,k)*eyekx(i)
    enddo
   enddo
  enddo

  call dfftw_execute_dft_c2r(fftw_plan_inv_x, temp, dvar)
  deallocate(temp)

 elseif (PERIODIC .AND. (.NOT. USE_FFT)) then
  call dbyd1(dvar,var,nx,dim2(rank)*nz,5)

 elseif (.NOT. PERIODIC) then
  call dbyd1(dvar,var,nx,dim2(rank)*nz,bc)
  
 endif
 dvar = dvar*stretchx

END SUBROUTINE ddx

!================================================================================================ 

SUBROUTINE ddy(var, dvar,bc)
 implicit none
 integer i,j,k,ierr,sendtag,recvtag,stat(MPI_STATUS_SIZE,2)
 integer n, statt(MPI_STATUS_SIZE), req(2), bc, reqq, tag, nelements
 real(kind=real64) var(nx,dim2(rank),nz), dvar(nx,dim2(rank),nz)
 real(kind=real64), dimension(:,:,:), allocatable :: trans, temp
 real(kind=real64), dimension(:,:), allocatable :: rectemp
 complex(kind=real64), allocatable, dimension(:,:,:) :: tempcomp


 allocate(temp(ny,dim1(rank), nz),trans(ny,dim1(rank),nz),rectemp(nx,nz))

 if ((bc == 4) .AND. (.NOT. periodic)) then
  nelements = nx*nz

  if(ystart(rank) == 1) then
   do k=1,nz
    do i=1,nx
     rectemp(i,k) = dvar(i,1,k) 
    enddo
   enddo
  endif

 
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)

  tag = 0
  do n=1,numtasks-1
   if (rank ==0) then
      call MPI_ISend(rectemp, nelements, MPI_DOUBLE_PRECISION, n, tag, MPI_COMM_WORLD, reqq, ierr)
      call MPI_WAIT(reqq, statt, ierr)
   endif
  enddo
   
  if (rank .NE. 0) &
  call MPI_Recv(rectemp, nelements, MPI_DOUBLE_PRECISION, 0, tag, MPI_COMM_WORLD, statt, ierr)

 ! Non - communicative transpose (local transpose) for the boundary condition
  j = 1
  do k=1,nz
   do i=1,dim1(rank)
     trans(1,i,k) = rectemp(i+fwd_sdispls(rank)-1,k)
   enddo
  enddo

  if(rank == numtasks-1) then
   do k=1,nz
    do i=1,nx
     rectemp(i,k) = dvar(i,dim2(rank),k) 
    enddo
   enddo
  endif

  call MPI_BARRIER(MPI_COMM_WORLD, ierr)

  tag = 0
  do n=1,numtasks-1
   if (rank == numtasks-1) then
      call MPI_ISend(rectemp, nelements, MPI_DOUBLE_PRECISION, n-1, tag, MPI_COMM_WORLD, reqq, ierr)
      call MPI_WAIT(reqq, statt, ierr)
   endif
  enddo
 
  if (rank .NE. numtasks-1) &
  call MPI_Recv(rectemp, nelements, MPI_DOUBLE_PRECISION, numtasks-1, tag, MPI_COMM_WORLD, statt, ierr)

 ! Non - communicative transpose (local transpose) for the boundary condition
  j = ny
  do k=1,nz
   do i=1,dim1(rank)
     trans(ny,i,k) = rectemp(i+fwd_sdispls(rank)-1,k)
   enddo
  enddo

 endif
  deallocate(rectemp)
 

 call transpose_3D_y(var, temp)
 
 if ((PERIODIC) .AND. (USE_FFT)) then
  allocate(tempcomp(ny/2+1,dim1(rank), nz))
  call dfftw_execute_dft_r2c(fftw_plan_fwd_y, temp, tempcomp)

  do k=1,nz
   do j=1,dim2(rank)
    do i=1,ny/2+1
     tempcomp(i,j,k) = tempcomp(i,j,k)*eyeky(i)
    enddo
   enddo
  enddo

  call dfftw_execute_dft_c2r(fftw_plan_inv_y, tempcomp, trans)
  deallocate(tempcomp)
 
 elseif (PERIODIC .AND. (.NOT. USE_FFT)) then
  call dbyd1(trans,temp,ny,dim1(rank)*nz,5)

 elseif (.NOT. PERIODIC) then
  call dbyd1(trans,temp,ny,dim1(rank)*nz,bc)
 endif

 call inv_transpose_3D_y(trans, dvar)
   
 deallocate(trans, temp)

 dvar = dvar*stretchy

END SUBROUTINE ddy


!================================================================================================ 

SUBROUTINE ddxyz(var1, dvar1, var2, dvar2, var3, dvar3, bc)

 ! TO OVERLAP COMMUNICATION AND COMPUTATION
 ! WHILE COMPUTING THE DERIVATIVES, ALSO SIMULATANEOUSLY PERFORM 
 ! PARALLEL ARRAY TRANSFORMS

 implicit none
 real(kind=real64), dimension(nx,dim2(rank),nz) :: var1, dvar1, var2, dvar2, var3, dvar3
 integer i, j, k, sendtag, recvtag, req(2*numtasks-2), req2(2*numtasks-2) , count
 integer ierr, stat(MPI_STATUS_SIZE, 2*numtasks-2), bc, stat2(MPI_STATUS_SIZE, 2*numtasks-2)
 real(kind=real64), dimension(ny, dim1(rank), nz) :: tran, dtran
 complex(kind=real64), allocatable, dimension(:,:,:) :: temp

 count = 2*numtasks-2

 ! TRANSPOSE THE y- variable
 ! Non - communicative transpose

 do j=1,dim1(rank)
  do i=1,dim2(rank)
   tran(j+inv_rdispls(rank)-1,i,:) = var2(i+inv_sdispls(rank)-1,j,:)
  enddo
 enddo

  ! Non-blocking send-recv transpose

 do k=1,numtasks-1
   i = modulo(rank+k,numtasks)
   j = modulo(rank-k,numtasks)
   sendtag = rank 
   recvtag = j 
   call MPI_ISEND(var2(inv_sdispls(i),1,1),1,inv_sendtypes(i),i,sendtag,MPI_COMM_WORLD,req(2*k-1),ierr)
   call MPI_IRECV(tran(inv_rdispls(j),1,1),1,inv_recvtypes(j),j,recvtag,MPI_COMM_WORLD,req(2*k),ierr)
 enddo

 ! WHILE PARALLEL INFORMATION TRANSFER, COMPUTE X-DERIVATIVE

 if ((PERIODIC) .AND. (USE_FFT)) then

  allocate(temp(nx/2+1,dim2(rank), nz))

  call dfftw_execute_dft_r2c(fftw_plan_fwd_x, var1, temp)

  do k=1,nz
   do j=1,dim2(rank)
    do i=1,nx/2+1
     temp(i,j,k) = temp(i,j,k)*eyekx(i)
    enddo
   enddo
  enddo

  call dfftw_execute_dft_c2r(fftw_plan_inv_x, temp, dvar1)

  deallocate(temp)

 elseif (PERIODIC .AND. (.NOT. USE_FFT)) then
  call dbyd1(dvar1,var1,nx,dim2(rank)*nz,5)

 elseif (.NOT. PERIODIC) then
  call dbyd1(dvar1,var1,nx,dim2(rank)*nz,bc)
  
 endif

 dvar1 = dvar1*stretchx


 ! ASSUMING THE PARALLEL TRANSPOSE IS FINISHED
 ! IN PRINCIPLE SHOULD CALL MPI_BARRIER BUT THAT 
 ! IS A WASTE (JUST HOPE THE NETWORK IS SUFFICIENTLY FAST)
 
 call MPI_WAITALL(count,req,stat,ierr)

 if ((PERIODIC) .AND. (USE_FFT)) then
  allocate(temp(ny/2+1,dim1(rank), nz))
  call dfftw_execute_dft_r2c(fftw_plan_fwd_y, tran, temp)

  do k=1,nz
   do j=1,dim2(rank)
    do i=1,ny/2+1
     temp(i,j,k) = temp(i,j,k)*eyeky(i)
    enddo
   enddo
  enddo

  call dfftw_execute_dft_c2r(fftw_plan_inv_y, temp, dtran)
  deallocate(temp)

 elseif (PERIODIC .AND. (.NOT. USE_FFT)) then
  call dbyd1(dtran,tran,ny,dim1(rank)*nz,5)

 elseif (.NOT. PERIODIC) then
  call dbyd1(dtran,tran,ny,dim1(rank)*nz,bc)
 endif

 dtran = dtran * stretchy

 ! INVERSE y- transpose
 ! Non - communicative transpose

 do j=1,dim1(rank)
  do i=1,dim2(rank)
   dvar2(j+inv_rdispls(rank)-1,i,:) = dtran(i+inv_sdispls(rank)-1,j,:)
  enddo
 enddo

  ! Non-blocking send-recv transpose

 do k=1,numtasks-1
   i = modulo(rank+k,numtasks)
   j = modulo(rank-k,numtasks)
   sendtag = rank
   recvtag = j
   call MPI_ISEND(dtran(inv_sdispls(i),1,1),1,inv_sendtypes(i),i,sendtag,MPI_COMM_WORLD,req2(2*k-1),ierr)
   call MPI_IRECV(dvar2(inv_rdispls(j),1,1),1,inv_recvtypes(j),j,recvtag,MPI_COMM_WORLD,req2(2*k),ierr)
 enddo

 ! WHILE INVERSE y-TRANSPOSE PARALLEL TRANSFER
 ! IS OCCURRING, COMPUTE THE Z-derivative

 if (compact_finite_diff) then
  call dbyd2(dvar3,var3,nx*dim2(rank),nz,bc)
  do k=1,nz
   dvar3(:,:,k) = dvar3(:,:,k)*stretch(k)
  enddo
 else
  call ddz_ninthorder(var3,dvar3)
 endif

 call MPI_WAITALL(count,req2,stat2,ierr)

END SUBROUTINE ddxyz


!================================================================================================

SUBROUTINE ddz_ninthorder(input, output)
 implicit none
 integer k
 real(kind=real64), dimension(nx, dim2(rank), nz) :: input, output
 real(kind=real64) consts(0:4)
 
 consts(0) = 0.5*(nz-1.)
 ! Second order differences 
 output(:,:,1) = (-0.5*input(:,:,3) + 2.0*input(:,:,2) - 1.5 * input(:,:,1))*stretch(1)*(nz-1.)
 output(:,:,2) = (- input(:,:,1) + input(:,:,3))*consts(0)*stretch(2)
 output(:,:,3) = (- input(:,:,2) + input(:,:,4))*consts(0)*stretch(3)
 output(:,:,4) = (- input(:,:,3) + input(:,:,5))*consts(0)*stretch(4)

 output(:,:,nz) = (-1.5*input(:,:,nz) +2. * input(:,:,nz-1) - 0.5*input(:,:,nz-2))*stretch(nz)*(nz-1.)
 output(:,:,nz-1) = (- input(:,:,nz-2) + input(:,:,nz))*consts(0)*stretch(nz-1)
 output(:,:,nz-2) = (- input(:,:,nz-3) + input(:,:,nz-1))*consts(0)*stretch(nz-2)
 output(:,:,nz-3) = (- input(:,:,nz-4) + input(:,:,nz-2))*consts(0)*stretch(nz-3)


 ! Implementing the optimized ninth-order centered differences of Bogey, Bailly & (someone else), JCP
 consts(1) = + 0.841570125482*(nz-1.) !* consts(0)
 consts(2) = - 0.244678631765*(nz-1.) !* consts(0)
 consts(3) = + 0.059463584768*(nz-1.) !* consts(0)
 consts(4) = - 0.007650904064*(nz-1.) !* consts(0)
 
 do k=5,nz-4
  output(:,:,k) = ((- input(:,:,k-4) + input(:,:,k+4))*consts(4) +  & 
    (- input(:,:,k-3) + input(:,:,k+3))*consts(3) +  & 
    (- input(:,:,k-2) + input(:,:,k+2))*consts(2) +  &
    (- input(:,:,k-1) + input(:,:,k+1))*consts(1) )*stretch(k)
 enddo

END SUBROUTINE ddz_ninthorder
!================================================================================================



SUBROUTINE d2dx2(var, dvar)
 implicit none
 real(kind=real64) var(nx,dim1(rank),nz), dvar(nx,dim1(rank),nz)
 integer i,j,k
 complex(kind=real64), allocatable, dimension(:,:,:) :: temp

 allocate(temp(nx/2+1,dim2(rank), nz))

 call dfftw_execute_dft_r2c(fftw_plan_fwd_x, var, temp)


 do k=1,nz
 do j=1,dim2(rank)
 do i=1,nx/2+1
  temp(i,j,k) = temp(i,j,k)*eyekx(i)**2.0/normx
 enddo
 enddo
 enddo
 !temp(limitx+1:nx/2+1,:,:) = 0.0

 call dfftw_execute_dft_c2r(fftw_plan_inv_x, temp, dvar)
 deallocate(temp)

END SUBROUTINE d2dx2

!================================================================================================
SUBROUTINE d2dy2(var, dvar)
 implicit none
 integer i,j,k,ierr,sendtag,recvtag,stat(MPI_STATUS_SIZE,2)
 integer n, statt(MPI_STATUS_SIZE), req(2), reqq, tag, nelements
 real(kind=real64) var(nx,dim2(rank),nz), dvar(nx,dim2(rank),nz)
 real(kind=real64), dimension(:,:,:), allocatable :: temp
 complex(kind=real64), dimension(:,:,:), allocatable :: trans

 allocate(temp(ny,dim1(rank), nz),trans(ny/2 + 1,dim1(rank),nz))

 ! Non - communicative transpose (local transpose)
 do k=1,nz
  do j=1,dim2(rank)
   do i=1,dim1(rank)
    temp(j+fwd_rdispls(rank)-1,i,k) = var(i+fwd_sdispls(rank)-1,j,k)
   enddo
  enddo
 enddo

 ! Non-blocking send-recv transpose
 do k=1,numtasks-1
   i = modulo(rank+k,numtasks)
   j = modulo(rank-k,numtasks)
   sendtag = 0
   recvtag = 0
   call MPI_ISEND(var(fwd_sdispls(i),1,1),1,fwd_sendtypes(i),i,sendtag,MPI_COMM_WORLD,req(1),ierr)
   call MPI_IRECV(temp(fwd_rdispls(j),1,1),1,fwd_recvtypes(j),j,recvtag,MPI_COMM_WORLD,req(2),ierr)
   call MPI_WAITALL(2,req,stat,ierr)
 enddo
 call dfftw_execute_dft_r2c(fftw_plan_fwd_y, temp, trans)

 do k=1,nz
  do j=1,dim2(rank)
   do i=1,ny/2+1
    trans(i,j,k) = trans(i,j,k)*eyeky(i)**2.0/normx
   enddo
  enddo
 enddo
 !temp(limity+1:ny/2+1,:,:) = 0.0

 call dfftw_execute_dft_c2r(fftw_plan_inv_y, trans, temp)
 deallocate(trans)


 ! Non - communicative transpose
 do k=1,nz
  do j=1,dim1(rank)
   do i=1,dim2(rank)
    dvar(j+inv_rdispls(rank)-1,i,k) = temp(i+inv_sdispls(rank)-1,j,k)
   enddo
  enddo
 enddo
  ! Non-blocking send-recv transpose

 do k=1,numtasks-1
   i = modulo(rank+k,numtasks)
   j = modulo(rank-k,numtasks)
   sendtag = 0
   recvtag = 0
   call MPI_ISEND(temp(inv_sdispls(i),1,1),1,inv_sendtypes(i),i,sendtag,MPI_COMM_WORLD,req(1),ierr)
   call MPI_IRECV(dvar(inv_rdispls(j),1,1),1,inv_recvtypes(j),j,recvtag,MPI_COMM_WORLD,req(2),ierr)
   call MPI_WAITALL(2,req,stat,ierr)
 enddo
 deallocate(temp)

END SUBROUTINE d2dy2

!================================================================================================


SUBROUTINE CURL(f_x, f_y, f_z, curl_x, curl_y, curl_z)
                                                                                                                                                            
  implicit none
  integer bcx,bcy,bcz
  real(kind=real64), dimension(nx, dim2(rank), nz) :: curl_x, curl_y, curl_z
  real(kind=real64), dimension(nx, dim2(rank), nz) :: f_x, f_y, f_z
  real(kind=real64), dimension(nx, dim2(rank),nz) ::  scr1, scr2, scr3

  call ddxyz(f_y, curl_z, f_z, curl_x, f_x, curl_y,1)
  call ddxyz(f_z, scr1, f_x, scr2, f_y, scr3,1)

  !call ddx(f_y, curl_z, bcx)
  !call ddy(f_x, temp, bcy)
  curl_z = curl_z - scr2
                                                                                                                                                            
  !call ddy(f_z, curl_x, bcy)
  !call ddz(f_y, temp, bcz)
  curl_x = curl_x - scr3
                                                                                                                                                            
  !call ddz(f_x, curl_y,bcz)
  !call ddx(f_z, temp, bcx)
  curl_y = curl_y - scr1
                                                                                                                                                            
                                                                                                                                                            
END SUBROUTINE CURL

!================================================================================================
                                                                                                                                                            
SUBROUTINE CROSS(a_x, a_y, a_z, b_x, b_y, b_z, c_x, c_y, c_z)
                                                                                                                                                            
  implicit none
  real(kind=real64), dimension(nx, dim2(rank), nz) :: a_x, a_y, a_z, b_x, b_y, b_z, c_x, c_y, c_z

  c_x = a_y * b_z - a_z * b_y
  c_y = b_x * a_z - a_x * b_z
  c_z = a_x * b_y - a_y * b_x

END SUBROUTINE CROSS

!================================================================================================


SUBROUTINE TRANSPOSE_3D_Y(input, output)

 implicit none
 integer i, j, k, sendtag, recvtag, req(2*numtasks-2), ierr, stat(MPI_STATUS_SIZE, 2*numtasks-2)
 real(kind=real64) input(nx, dim2(rank), nz), output(ny, dim1(rank), nz)

 ! Non - communicative transpose (local transpose)
 do k=1,nz
  do j=1,dim2(rank)
   do i=1,dim1(rank)
    output(j+fwd_rdispls(rank)-1,i,k) = input(i+fwd_sdispls(rank)-1,j,k)
   enddo
  enddo
 enddo

 ! Non-blocking send-recv transpose
 do k=1,numtasks-1
   i = modulo(rank+k,numtasks)
   j = modulo(rank-k,numtasks)
   sendtag = rank
   recvtag = j
   call MPI_ISEND(input(fwd_sdispls(i),1,1),1,fwd_sendtypes(i),i,sendtag,MPI_COMM_WORLD,req(2*k-1),ierr)
   call MPI_IRECV(output(fwd_rdispls(j),1,1),1,fwd_recvtypes(j),j,recvtag,MPI_COMM_WORLD,req(2*k),ierr)
 enddo
 call MPI_WAITALL(2*numtasks-2,req,stat,ierr)

END SUBROUTINE TRANSPOSE_3D_Y

!================================================================================================ 


SUBROUTINE INV_TRANSPOSE_3D_Y(input, output)

 implicit none
 integer i, j, k, sendtag, recvtag, req(2*numtasks-2), ierr, stat(MPI_STATUS_SIZE, 2*numtasks-2)
 real(kind=real64) input(ny, dim1(rank), nz), output(nx, dim2(rank), nz)

 ! Non - communicative transpose

 do k=1,nz
  do j=1,dim1(rank)
   do i=1,dim2(rank)
    output(j+inv_rdispls(rank)-1,i,k) = input(i+inv_sdispls(rank)-1,j,k)
   enddo
  enddo
 enddo
  ! Non-blocking send-recv transpose

 do k=1,numtasks-1
   i = modulo(rank+k,numtasks)
   j = modulo(rank-k,numtasks)
   sendtag = rank
   recvtag = j
   call MPI_ISEND(input(inv_sdispls(i),1,1),1,inv_sendtypes(i),i,sendtag,MPI_COMM_WORLD,req(2*k-1),ierr)
   call MPI_IRECV(output(inv_rdispls(j),1,1),1,inv_recvtypes(j),j,recvtag,MPI_COMM_WORLD,req(2*k),ierr)
 enddo
 call MPI_WAITALL(2*numtasks-2,req,stat,ierr)

END SUBROUTINE INV_TRANSPOSE_3D_Y


!================================================================================================

END MODULE DERIVATIVES
