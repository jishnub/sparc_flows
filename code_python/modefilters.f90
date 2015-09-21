SUBROUTINE FMODE_FILTER(nt,nx,fmode)

    implicit none
    integer, intent(in) :: nt,nx
    integer i,j
    real*8, intent(inout) :: fmode(nx, 1, nt)
    real*8 f_low,df,k(nx),dt,f_mode_const
    real*8 Poly(0:2), f0(nx),w(nt),f(nt),f1(nx),d,delta
    real*8 xlength,outputcad,pi
    parameter (xlength = 800.0 * 10**(8))
    parameter (outputcad = 30.0)
    parameter (pi=acos(dble(-1.0)))

    dt = outputcad

    open(123,file='modefilters.data',action='read')
    read(123,*) f_mode_const,Poly(0),Poly(1),Poly(2),f_low,df

    call distmat(nx,1,k) 
    call distmat(nt,1,w) 
    k = abs(k) * 2.*pi/(xlength*10.**(-8.)*nx/(nx-1.))
    w = abs(w) * 2.*pi/(nt*dt)

    f0=f_mode_const*k**0.5
    f1=Poly(0) + Poly(1)*k +Poly(2)*k**2.
    f = w/(2.*pi)*1e3


    fmode = 0.0
    do i=1,nx
    delta = (f1(i) - f0(i))
    do j=1,nt
     d = f(j) - f0(i)
     if ((d .lt. delta) .and. (d>0)) &
        fmode(i,1,j) = 0.5*(1.+cos(pi*(2.0*d/delta-1)))
    enddo
    enddo 

    do j=1,nt
    if (f(j) .lt. f_low+df) &
      fmode(:,1,j) = fmode(:,1,j)*0.5*(1.+cos(pi*(f(j)-(f_low+df))/df) )
    if (f(j) .lt. f_low) fmode(:,1,j) = 0.
    enddo
    
    close(123)

END SUBROUTINE FMODE_FILTER

!=======================================================================

SUBROUTINE PMODE_FILTER(nt,nx, pmode)

    implicit none
    integer, intent(in) :: nt,nx
    integer i,j
    real*8, intent(inout) :: pmode(nx, 1, nt)
    real*8 f_low,df,k(nx),dt,f_mode_const
    real*8 Poly(0:2), f0(nx),w(nt),f(nt),f1(nx),d,delta
    real*8 xlength,outputcad,pi
    parameter (xlength = 800.0 * 10**(8))
    parameter (outputcad = 30.0)
    parameter (pi=acos(dble(-1.0)))

    dt = outputcad
    
    open(123,file='modefilters.data',action='read')
    read(123,*)
    read(123,*) f_mode_const,Poly(0),Poly(1),Poly(2),f_low,df

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
        pmode(i,1,j) = 0.5*(1.+cos(pi*(2.0*d/delta-1)))
     end if   
    enddo
    enddo 

    do j=1,nt
    if (f(j) .lt. f_low) pmode(:,1,j) = 0.
    if (f(j) .lt. f_low+df) &
      pmode(:,1,j) = pmode(:,1,j)*0.5*(1.+cos(pi*(f(j)-(f_low+df))/df) )
    enddo
    
    close(123)

END SUBROUTINE PMODE_FILTER


!=======================================================================

SUBROUTINE P2MODE_FILTER(nt,nx, p2mode)

    implicit none
    integer, intent(in) :: nt,nx
    integer i,j
    real*8, intent(inout) :: p2mode(nx, 1, nt)
    real*8 f_low,df,k(nx),dt,f_mode_const
    real*8 Poly(0:2), f0(nx),w(nt),f(nt),f1(nx),d,delta
    real*8 xlength,outputcad,pi
    parameter (xlength = 800.0 * 10**(8))
    parameter (outputcad = 30.0)
    parameter (pi=acos(dble(-1.0)))

    dt = outputcad
    
    open(123,file='modefilters.data',action='read')
    read(123,*)
    read(123,*)
    read(123,*) f_mode_const,Poly(0),Poly(1),Poly(2),f_low,df

    call distmat(nx,1,k) 
    call distmat(nt,1,w) 
    k = abs(k) * 2.*pi/(xlength*10.**(-8.)*nx/(nx-1.))
    w = abs(w) * 2.*pi/(nt*dt)

!~     f0=P1(0) + P1(1)*k +P1(2)*k**2.
!~     f1=P2(0) + P2(1)*k +P2(2)*k**2.
    f0=f_mode_const*abs(k)**0.5
    f1=Poly(0) + Poly(1)*k +Poly(2)*k**2.
    f = w/(2.*pi)*1e3

    p2mode = 0.0
    do i=1,nx
    delta = (f1(i) - f0(i))
    do j=1,nt
     d = f(j) - f0(i)
     if ((d .lt. delta) .and. (d .gt. 0)) then
        p2mode(i,1,j) = 0.5*(1.+cos(pi*(2.0*d/delta-1)))
     end if   
    enddo
    enddo 

    do j=1,nt
    if (f(j) .lt. f_low) p2mode(:,1,j) = 0.
    if (f(j) .lt. f_low+df) &
      p2mode(:,1,j) = p2mode(:,1,j)*0.5*(1.+cos(pi*(f(j)-(f_low+df))/df) )
    enddo
    
    close(123)

END SUBROUTINE P2MODE_FILTER

!=======================================================================

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
