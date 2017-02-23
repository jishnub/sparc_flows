
SUBROUTINE COMPUTE_TT_HANASOGE(u0, u, tau, dt, nt)

 implicit none
 integer, intent(in) :: nt
 integer l1,l2, loc, i
 real*8 cc(-nt+1:nt-1), t(-nt+1:nt-1)
 real*8, intent(in) :: u0(nt), u(nt), dt
 real*8 times(3), invmat(3,3), p1, p2, mat(3,3)
 real*8, intent(out) :: tau
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

END SUBROUTINE COMPUTE_TT_HANASOGE

!================================================================================

SUBROUTINE COMPUTE_TT_GIZONBIRCH(u0,u,tau,dt,nt, lef, rig)

    implicit none
    INCLUDE 'fftw3.f'
    integer, parameter :: degree=3
    integer, intent(in) :: nt, lef, rig
    real*8, intent(in) :: u0(nt), u(nt)
    real*8 dt
    real*8, intent(out) :: tau
    real*8 window(nt),functemp(nt)
    integer k,i,num_real_roots
    integer*8 plan
    real*8 polycoeffs(0:degree)
    real*8 :: roots(degree)
    complex*16 wu0(nt/2+1,1:degree+1)
    real*8 u0dot(nt,1:degree+1)
    complex*16, parameter :: eye = (0.0,1.0)
    real*8, parameter :: pi=acos(dble(-1.0))


    polycoeffs = 0

    window = 0.0
    window(lef:rig) = 1.0

    do i=1,degree+1
       call dfftw_plan_dft_r2c_1d(plan,nt,u0,wu0(:,i),&
               FFTW_ESTIMATE)
       call dfftw_execute_dft_r2c(plan, u0, wu0(:,i))
       call dfftw_destroy_plan(plan)

       do k=1,nt/2
           wu0(k,i)=wu0(k,i)*(eye*2*pi*(k-1.0)/(nt*dt))**i
       end do

       call dfftw_plan_dft_c2r_1d(plan,nt,wu0(:,i),u0dot(:,i),&
               FFTW_ESTIMATE)
       call dfftw_execute_dft_c2r(plan, wu0(:,i),u0dot(:,i))
       call dfftw_destroy_plan(plan)
    enddo

    u0dot=u0dot/nt

    ! coefficients of dchi/dt

    do i=0,degree
        functemp = 0
        if (i .eq. 0) then
            functemp = 2*u0dot(:,1)*(u-u0)
        elseif (i .eq. 1) then
            functemp = 2*(u0dot(:,1)**2+(-u+u0)*u0dot(:,2))
        elseif (i .eq. 2) then
            functemp = -3*u0dot(:,1)*u0dot(:,2)+&
                (u-u0)*u0dot(:,3)
        elseif (i .eq. 3) then
            functemp = (3*u0dot(:,2)**2+4*u0dot(:,1)*u0dot(:,3)&
                +(-u+u0)*u0dot(:,4))/3.0

        endif

        call integrate_time(window*functemp,&
                polycoeffs(degree-i),dt,nt)
    end do

    polycoeffs = polycoeffs/polycoeffs(0)

    tau=-polycoeffs(degree)/polycoeffs(degree-1)

    call polyroots(polycoeffs,degree,num_real_roots,roots)

    tau=roots(minloc(abs(roots(1:num_real_roots)-tau),dim=1))


END SUBROUTINE COMPUTE_TT_GIZONBIRCH

!===============================================================================

SUBROUTINE COMPUTE_TT_GIZONBIRCH3(u0,u,tau,dt,nt, lef, rig)

    implicit none
    INCLUDE 'fftw3.f'
    integer, parameter :: degree=3
    integer, intent(in) :: nt, lef, rig
    real*8, intent(in) :: u0(nt), u(nt)
    real*8 dt
    real*8, intent(out) :: tau
    real*8 window(nt),functemp(nt)
    integer k,i,num_real_roots
    integer*8 plan
    real*8 polycoeffs(0:degree)
    real*8 :: roots(degree)
    complex*16 wu0(nt/2+1,1:degree+1)
    real*8 u0dot(nt,1:degree+1)
    complex*16, parameter :: eye = (0.0,1.0)
    real*8, parameter :: pi=acos(dble(-1.0))


    polycoeffs = 0

    window = 0.0
    window(lef:rig) = 1.0

    do i=1,degree+1
       call dfftw_plan_dft_r2c_1d(plan,nt,u0,wu0(:,i),&
               FFTW_ESTIMATE)
       call dfftw_execute_dft_r2c(plan, u0, wu0(:,i))
       call dfftw_destroy_plan(plan)

       do k=1,nt/2
           wu0(k,i)=wu0(k,i)*(eye*2*pi*(k-1.0)/(nt*dt))**i
       end do

       call dfftw_plan_dft_c2r_1d(plan,nt,wu0(:,i),u0dot(:,i),&
               FFTW_ESTIMATE)
       call dfftw_execute_dft_c2r(plan, wu0(:,i),u0dot(:,i))
       call dfftw_destroy_plan(plan)
    enddo

    u0dot=u0dot/nt

    ! coefficients of dchi/dt

    do i=0,degree
        functemp = 0
        if (i .eq. 0) then
            functemp = 2*u0dot(:,1)*(u-u0)
        elseif (i .eq. 1) then
            functemp = 2*(u0dot(:,1)**2+(-u+u0)*u0dot(:,2))
        elseif (i .eq. 2) then
            functemp = -3*u0dot(:,1)*u0dot(:,2)+&
                (u-u0)*u0dot(:,3)
        elseif (i .eq. 3) then
            functemp = (3*u0dot(:,2)**2+4*u0dot(:,1)*u0dot(:,3)&
                +(-u+u0)*u0dot(:,4))/3.0

        endif

        call integrate_time(window*functemp,&
                polycoeffs(degree-i),dt,nt)
    end do

    polycoeffs = polycoeffs/polycoeffs(0)

    tau=-polycoeffs(degree)/polycoeffs(degree-1)

    call polyroots(polycoeffs,degree,num_real_roots,roots)

    tau=roots(minloc(abs(roots(1:num_real_roots)-tau),dim=1))


END SUBROUTINE COMPUTE_TT_GIZONBIRCH3

SUBROUTINE COMPUTE_TT_GIZONBIRCH2(u0,u,tau,dt,nt, lef, rig)

    implicit none
    INCLUDE 'fftw3.f'
    integer, parameter :: degree=2
    integer, intent(in) :: nt, lef, rig
    real*8, intent(in) :: u0(nt), u(nt)
    real*8 dt
    real*8, intent(out) :: tau
    real*8 window(nt),functemp(nt)
    integer k,i,num_real_roots
    integer*8 plan
    real*8 polycoeffs(0:degree)
    real*8 :: roots(degree)
    complex*16 wu0(nt/2+1,1:degree+1)
    real*8 u0dot(nt,1:degree+1)
    complex*16, parameter :: eye = (0.0,1.0)
    real*8, parameter :: pi=acos(dble(-1.0))


    polycoeffs = 0

    window = 0.0
    window(lef:rig) = 1.0

    do i=1,degree+1
       call dfftw_plan_dft_r2c_1d(plan,nt,u0,wu0(:,i),&
               FFTW_ESTIMATE)
       call dfftw_execute_dft_r2c(plan, u0, wu0(:,i))
       call dfftw_destroy_plan(plan)

       do k=1,nt/2
           wu0(k,i)=wu0(k,i)*(eye*2*pi*(k-1.0)/(nt*dt))**i
       end do

       call dfftw_plan_dft_c2r_1d(plan,nt,wu0(:,i),u0dot(:,i),&
               FFTW_ESTIMATE)
       call dfftw_execute_dft_c2r(plan, wu0(:,i),u0dot(:,i))
       call dfftw_destroy_plan(plan)
    enddo

    u0dot=u0dot/nt

    ! coefficients of dchi/dt

    do i=0,degree
        functemp = 0
        if (i .eq. 0) then
            functemp = 2*u0dot(:,1)*(u-u0)
        elseif (i .eq. 1) then
            functemp = 2*(u0dot(:,1)**2+(-u+u0)*u0dot(:,2))
        elseif (i .eq. 2) then
            functemp = -3*u0dot(:,1)*u0dot(:,2)+&
                (u-u0)*u0dot(:,3)
        elseif (i .eq. 3) then
            functemp = (3*u0dot(:,2)**2+4*u0dot(:,1)*u0dot(:,3)&
                +(-u+u0)*u0dot(:,4))/3.0

        endif

        call integrate_time(window*functemp,&
                polycoeffs(degree-i),dt,nt)
    end do

    polycoeffs = polycoeffs/polycoeffs(0)

    tau=-polycoeffs(degree)/polycoeffs(degree-1)

    call polyroots(polycoeffs,degree,num_real_roots,roots)

    tau=roots(minloc(abs(roots(1:num_real_roots)-tau),dim=1))


END SUBROUTINE COMPUTE_TT_GIZONBIRCH2

SUBROUTINE COMPUTE_TT_GIZONBIRCH1(u0,u,tau,dt,nt, lef, rig)

    implicit none
    INCLUDE 'fftw3.f'
    integer, parameter :: degree=1
    integer, intent(in) :: nt, lef, rig
    real*8, intent(in) :: u0(nt), u(nt)
    real*8 dt
    real*8, intent(out) :: tau
    real*8 window(nt),functemp(nt)
    integer k,i,num_real_roots
    integer*8 plan
    real*8 polycoeffs(0:degree)
    real*8 :: roots(degree)
    complex*16 wu0(nt/2+1,1:degree+1)
    real*8 u0dot(nt,1:degree+1)
    complex*16, parameter :: eye = (0.0,1.0)
    real*8, parameter :: pi=acos(dble(-1.0))


    polycoeffs = 0

    window = 0.0
    window(lef:rig) = 1.0

    do i=1,degree+1
       call dfftw_plan_dft_r2c_1d(plan,nt,u0,wu0(:,i),&
               FFTW_ESTIMATE)
       call dfftw_execute_dft_r2c(plan, u0, wu0(:,i))
       call dfftw_destroy_plan(plan)

       do k=1,nt/2
           wu0(k,i)=wu0(k,i)*(eye*2*pi*(k-1.0)/(nt*dt))**i
       end do

       call dfftw_plan_dft_c2r_1d(plan,nt,wu0(:,i),u0dot(:,i),&
               FFTW_ESTIMATE)
       call dfftw_execute_dft_c2r(plan, wu0(:,i),u0dot(:,i))
       call dfftw_destroy_plan(plan)
    enddo

    u0dot=u0dot/nt

    ! coefficients of dchi/dt

    do i=0,degree
        functemp = 0
        if (i .eq. 0) then
            functemp = 2*u0dot(:,1)*(u-u0)
        elseif (i .eq. 1) then
            functemp = 2*(u0dot(:,1)**2+(-u+u0)*u0dot(:,2))
        elseif (i .eq. 2) then
            functemp = -3*u0dot(:,1)*u0dot(:,2)+&
                (u-u0)*u0dot(:,3)
        elseif (i .eq. 3) then
            functemp = (3*u0dot(:,2)**2+4*u0dot(:,1)*u0dot(:,3)&
                +(-u+u0)*u0dot(:,4))/3.0
        endif

        call integrate_time(window*functemp,&
                polycoeffs(degree-i),dt,nt)
    end do

    polycoeffs = polycoeffs/polycoeffs(0)

    tau=-polycoeffs(degree)/polycoeffs(degree-1)

    call polyroots(polycoeffs,degree,num_real_roots,roots)

    tau=roots(minloc(abs(roots(1:num_real_roots)-tau),dim=1))


END SUBROUTINE COMPUTE_TT_GIZONBIRCH1

!================================================================================

SUBROUTINE INTEGRATE_TIME(f,int_f,dt,nt)

    implicit none
    integer, intent(in) :: nt
    real*8, intent(in) :: dt,f(nt)
    real*8, intent(out) ::  int_f

    call simpson_regular(f,dt,nt,int_f)


END SUBROUTINE INTEGRATE_TIME

!================================================================================

SUBROUTINE SIMPSON_REGULAR(f,dvar,nt,int_f)

    implicit none
    integer, intent(in) :: nt
    real*8, intent(in) :: dvar,f(nt)
    integer k
    real*8, intent(out) :: int_f

    do k=1,size(f)

        if ((k == 1) .or. (k == nt)) then
            int_f = int_f + f(k)
            cycle
        end if

        if (mod(k,2) == 0) then
            int_f = int_f + 4*f(k)
        else
            int_f = int_f + 2*f(k)
        endif

    end do

    int_f=int_f*dvar/3.0

END SUBROUTINE SIMPSON_REGULAR
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


subroutine POLYROOTS(polycoeffs,N,num_real_roots,realroots)
 implicit none
 integer, intent(in) :: N
 real*8, intent(in) :: polycoeffs(0:N)
 real*8 companion(N,N),wi(N),wr(N)
 real*8, intent(out) :: realroots(N)
 integer,intent(out) :: num_real_roots
 integer i,info,lwork,lwmax,temp
 parameter(lwmax=100)
 real*8 work(lwmax),vl(N,N),vr(N,N)

 companion=0
 do i=1,N
   companion(1,i) = -polycoeffs(i)
 end do

 do i=2,N
   companion(i,i-1) = 1
 enddo

 LWORK = -1
 CALL DGEEV( 'N', 'N', N, companion, N, WR, WI, VL, N,&
              VR, N, WORK, LWORK, INFO )
 LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )

 call DGEEV('N','N',N,companion,N,wr,wi,vl,N,&
               vr,N,work,lwork,info)

   num_real_roots=0
   do i=1,N
       if (abs(wi(i))<1e-10) num_real_roots=num_real_roots+1
   end do

  temp = 1

  do i=1,N
      if (abs(wi(i))<1e-10) then
          realroots(temp) = wr(i)
          temp = temp+1
      end if
  end do

end subroutine POLYROOTS
