SUBROUTINE FMODE_FILTER(nt,outputcad,nx,xlength,fmode)
    implicit none
    integer, intent(in) :: nt,nx
    real*8, intent(in) :: xlength,outputcad
    integer i,j
    real*8, intent(out) :: fmode(nx, 1, nt)
    real*8 f_low,df,k(nx),dt!,f_mode_const
    real*8 Poly(0:2),Polylow(0:2), f0(nx),w(nt),f(nt),f1(nx),d,delta
    real*8 pi
    parameter (pi=3.141592654)

    dt = outputcad

!~     f_mode_const=2.056565
  
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
    k = abs(k) * 2.*pi/(xlength *nx/(nx-1.))
    w = abs(w) * 2.*pi/(nt*dt)

!~     f0=f_mode_const*k**0.5
    f0=Polylow(0) + Polylow(1)*k +Polylow(2)*k**2.
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

SUBROUTINE PMODE_FILTER(nt,outputcad,nx,xlength,pmode)
    
    implicit none
    integer, intent(in) :: nt,nx
    real*8, intent(in) :: xlength,outputcad
    integer i,j
    real*8, intent(out) :: pmode(nx, 1, nt)
    real*8 f_low,df,k(nx),dt!,f_mode_const
    real*8 Poly(0:2), f0(nx),w(nt),f(nt),f1(nx),d,delta, Polylow(0:2)
    real*8 pi
    parameter (pi=3.141592654)

    dt = outputcad

!~     f_mode_const=3.15

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
    k = abs(k) * 2.*pi/(xlength*nx/(nx-1.))
    w = abs(w) * 2.*pi/(nt*dt)

!~     f0=f_mode_const*abs(k)**0.5
    f1=Poly(0) + Poly(1)*k +Poly(2)*k**2.
    f0=Polylow(0) + Polylow(1)*k +Polylow(2)*k**2.
    f = w/(2.*pi)*1e3

    pmode = 0.0
    do i=1,nx
    delta = (f1(i) - f0(i))
    do j=1,nt
     d = f(j) - f0(i)
     if ((d .lt. delta) .and. (d .gt. 0)) then
        pmode(i,1,j) = 0.5*(1.+cos(pi*(2.0*d/delta-1)))
!~         pmode(i,1,j) = 1
     end if   
    enddo
    enddo 

!~     do j=1,nt
!~     if (f(j) .lt. f_low) pmode(:,1,j) = 0.
!~     if (f(j) .lt. f_low+df) &
!~       pmode(:,1,j) = pmode(:,1,j)*0.5*(1.+cos(pi*(f(j)-(f_low+df))/df) )
!~     enddo
    
    close(123)

END SUBROUTINE PMODE_FILTER

!=======================================================================

SUBROUTINE P2MODE_FILTER(nt,outputcad,nx,xlength,p2mode)

    implicit none
    integer, intent(in) :: nt,nx
    real*8, intent(in) :: xlength,outputcad
    integer i,j
    real*8, intent(out) :: p2mode(nx, 1, nt)
    real*8 f_low,df,k(nx),dt!,f_mode_const
    real*8 Poly(0:2),Polylow(0:2), f0(nx),w(nt),f(nt),f1(nx),d,delta
    real*8 pi
    parameter (pi=3.141592654)

    dt = outputcad
    
!~     f_mode_const=4.1004

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
    k = abs(k) * 2.*pi/(xlength *nx/(nx-1.))
    w = abs(w) * 2.*pi/(nt*dt)


!~     f0=f_mode_const*abs(k)**0.5
    f0=Polylow(0) + Polylow(1)*k +Polylow(2)*k**2.
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

SUBROUTINE P3MODE_FILTER(nt,outputcad, nx,xlength,pmode)

    implicit none
    integer, intent(in) :: nt,nx
    real*8, intent(in) :: xlength,outputcad
    integer i,j
    real*8, intent(out) :: pmode(nx, 1, nt)
    real*8 f_low,df,k(nx),dt!,f_mode_const
    real*8 Poly(0:2),Polylow(0:2), f0(nx),w(nt),f(nt),f1(nx),d,delta
    real*8 pi
    parameter (pi=3.141592654)

    dt = outputcad

!~     f_mode_const=4.7

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
    k = abs(k) * 2.*pi/(xlength *nx/(nx-1.))
    w = abs(w) * 2.*pi/(nt*dt)

!~     f0=f_mode_const*abs(k)**0.5
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


END SUBROUTINE P3MODE_FILTER

!=======================================================================


SUBROUTINE P4MODE_FILTER(nt,outputcad,nx,xlength,pmode)

    implicit none
    integer, intent(in) :: nt,nx
    real*8, intent(in) :: xlength,outputcad
    integer i,j
    real*8, intent(out) :: pmode(nx, 1, nt)
    real*8 f_low,df,k(nx),dt!,f_mode_const
    real*8 Poly(0:2), f0(nx),w(nt),f(nt),f1(nx),d,delta,Polylow(0:2)
    real*8 pi
    parameter (pi=3.141592654)

    dt = outputcad

!~     f_mode_const=5.4

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
    k = abs(k) * 2.*pi/(xlength *nx/(nx-1.))
    w = abs(w) * 2.*pi/(nt*dt)

!~     f0=f_mode_const*abs(k)**0.5
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
   
END SUBROUTINE P4MODE_FILTER

!=======================================================================


SUBROUTINE P5MODE_FILTER(nt,outputcad,nx,xlength,pmode)

    implicit none
    integer, intent(in) :: nt,nx
    real*8, intent(in) :: xlength,outputcad
    integer i,j
    real*8, intent(out) :: pmode(nx, 1, nt)
    real*8 f_low,df,k(nx),dt,f_mode_const
    real*8 Poly(0:2),Polylow(0:2), f0(nx),w(nt),f(nt),f1(nx),d,delta
    real*8 pi
    parameter (pi=3.141592654)

    dt = outputcad


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
    k = abs(k) * 2.*pi/(xlength *nx/(nx-1.))
    w = abs(w) * 2.*pi/(nt*dt)

!~     f0=f_mode_const*abs(k)**0.5
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
   

END SUBROUTINE P5MODE_FILTER

!=======================================================================


SUBROUTINE P6MODE_FILTER(nt,outputcad,nx,xlength,pmode)

    implicit none
    integer, intent(in) :: nt,nx
    real*8, intent(in) :: xlength,outputcad
    integer i,j
    real*8, intent(out) :: pmode(nx, 1, nt)
    real*8 f_low,df,k(nx),dt,f_mode_const
    real*8 Poly(0:2),Polylow(0:2), f0(nx),w(nt),f(nt),f1(nx),d,delta
    real*8 pi
    parameter (pi=3.141592654)

  dt = outputcad
  

    Poly(0)=2.65
    Poly(1)=5.7
    Poly(2)=-1.2

    Polylow(0)=2.4
    Polylow(1)=5.4
    Polylow(2)=-1.3

    f_low = 1.6
    df = 0.5 
  
  call distmat(nx,1,k) 
  call distmat(nt,1,w) 
  k = abs(k) * 2.*pi/(xlength*nx/(nx-1.))
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
END SUBROUTINE P6MODE_FILTER
!================================================================================

SUBROUTINE P7MODE_FILTER(nt,outputcad,nx,xlength,pmode)

    implicit none
    integer, intent(in) :: nt,nx
    real*8, intent(in) :: xlength,outputcad
    integer i,j
    real*8, intent(out) :: pmode(nx, 1, nt)
    real*8 f_low,df,k(nx),dt,f_mode_const
    real*8 Poly(0:2),Polylow(0:2), f0(nx),w(nt),f(nt),f1(nx),d,delta
    real*8 pi
    parameter (pi=3.141592654)

    dt = outputcad

    Poly(0)=2.9
    Poly(1)=5.9
    Poly(2)=-1.3

    Polylow(0)=2.65
    Polylow(1)=5.7
    Polylow(2)=-1.2

    f_low = 1.6
    df = 0.5  

    call distmat(nx,1,k) 
    call distmat(nt,1,w) 
    k = abs(k) * 2.*pi/(xlength *nx/(nx-1.))
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
   

END SUBROUTINE P7MODE_FILTER

!==========================================================================================

SUBROUTINE LARGE_DIST_PMODE_FILTER(nt,outputcad,nx,xlength,pmode)

  implicit none
    integer, intent(in) :: nt,nx
    real*8, intent(in) :: xlength,outputcad
    integer i,j
    real*8, intent(out) :: pmode(nx, 1, nt)
    real*8 f_low,df,k(nx),dt,f_mode_const
    real*8 Poly(0:2),Polylow(0:2), f0(nx),w(nt),f(nt),f1(nx),d,delta
    real*8 pi
    parameter (pi=3.141592654)

  
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
  k = abs(k) * 2.*pi/(xlength *nx/(nx-1.))
  w = abs(w) * 2.*pi/(nt*dt)
  
  f0=Polylow(0) + Polylow(1)*k +Polylow(2)*k**2.
  f1=Poly(0) + Poly(1)*k +Poly(2)*k**2.
  f = w/(2.*pi)*1e3
  
  pmode = 1.0
!~   do i=1,nx
!~    delta = (f1(i) - f0(i))
!~     do j=1,nt
!~      d = f(j) - f0(i)
!~      if ((d .lt. delta) .and. (d .gt. 0)) then
!~         pmode(i,1,j) = 1
!~      end if   
!~     enddo
!~    enddo 
   
!~    do j=1,nt
!~     if (f(j) .lt. f_low) pmode(:,1,j) = 0.
!~     if (f(j) .lt. f_low+df) &
!~       pmode(:,1,j) = pmode(:,1,j)*0.5*(1.+cos(pi*(f(j)-(f_low+df))/df) )
!~    enddo
   

END SUBROUTINE LARGE_DIST_PMODE_FILTER
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
