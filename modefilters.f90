SUBROUTINE FMODE_FILTER(nt,outputcad,nx,xlength,fmode)
    implicit none
    integer, intent(in) :: nt,nx
    real*8, intent(in) :: xlength,outputcad
    integer i,j
    real*8, intent(out) :: fmode(nx, 1, nt)
    real*8 f_low,df,dt
    real*8 :: f0(nx),f1(nx),k(nx)
    real*8 :: Poly(0:2),Polylow(0:2), w(nt),d,delta,f(nt)
    real*8 pi
    parameter (pi=3.141592654)

    dt = outputcad

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
    k = abs(k) * 2.*pi/xlength
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
        fmode(i,1,j) = 0.5*(1.+cos(pi*(2.0*d/delta-1)))
    enddo
    enddo

    do j=1,nt
    if (f(j) .lt. f_low+df) &
      fmode(:,1,j) = fmode(:,1,j)*0.5*(1.+cos(pi*(f(j)-(f_low+df))/df) )
    if (f(j) .lt. f_low) fmode(:,1,j) = 0.
    enddo
END SUBROUTINE FMODE_FILTER

!=======================================================================

!=======================================================================

SUBROUTINE FMODE_FILTER_ALT(nt,outputcad,nx,xlength,fmode)
    implicit none
    integer, intent(in) :: nt,nx
    real*8, intent(in) :: xlength,outputcad
    integer i,j
    real*8, intent(out) :: fmode(nx, 1, nt)
    real*8 f_low,df,dt
    real*8 :: f0(nx),f1(nx),k(nx)
    real*8 :: Poly(0:2),Polylow(0:2), w(nt),d,delta,f(nt)
    real*8 pi
    parameter (pi=3.141592654)

    dt = outputcad

    Poly(0)=1.1
    Poly(1)=1.9
    Poly(2)=-0.2

    Polylow(0)=0.7
    Polylow(1)=1.7
    Polylow(2)=-0.2

    f_low = 1.1

    call distmat(nx,1,k)
    call distmat(nt,1,w)
    k = abs(k) * 2.*pi/xlength
    w = abs(w) * 2.*pi/(nt*dt)

    f0=Polylow(0) + Polylow(1)*k +Polylow(2)*k**2.
    f1=Poly(0) + Poly(1)*k +Poly(2)*k**2.

    f = w/(2.*pi)*1e3

    df = f(2)-f(1)

    fmode = 0.0
    do i=1,nx
    delta = (f1(i) - f0(i))
    do j=1,nt
     d = f(j) - f0(i)
     if ((d .lt. delta) .and. (d>(delta-2*df))) then
        fmode(i,1,j) = cos(2*pi/(8*df)*(f(j)-(f1(i)-2*df)))
    else if ((d .lt. (2*df)) .and. (d>0)) then
        fmode(i,1,j) = cos(2*pi/(8*df)*(f(j)-(f0(i)+2*df)))
     elseif ((d > 2*df ) .and. (d<(delta-2*df))) then 
       fmode(i,1,j) = 1.0
     endif
    enddo
    enddo

    do j=1,nt
    if (f(j) .lt. f_low+0.5) &
      fmode(:,1,j) = fmode(:,1,j)*0.5*(1.+cos(pi*(f(j)-(f_low+0.5))/0.5) )
    if (f(j) .lt. f_low) fmode(:,1,j) = 0.
    enddo
END SUBROUTINE FMODE_FILTER_ALT

!=======================================================================

SUBROUTINE P1MODE_FILTER(nt,outputcad,nx,xlength,pmode)

    implicit none
    integer, intent(in) :: nt,nx
    real*8, intent(in) :: xlength,outputcad
    integer i,j
    real*8, intent(out) :: pmode(nx, 1, nt)
    real*8 f_low,df,k(nx),dt
    real*8 Poly(0:2), f0(nx),w(nt),f(nt),f1(nx),d,delta, Polylow(0:2)
    real*8 pi
    parameter (pi=3.141592654)

    dt = outputcad

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
     end if
    enddo
    enddo

     do j=1,nt
     if (f(j) .lt. f_low) pmode(:,1,j) = 0.
     if (f(j) .lt. f_low+df) &
       pmode(:,1,j) = pmode(:,1,j)*0.5*(1.+cos(pi*(f(j)-(f_low+df))/df) )
     enddo
END SUBROUTINE P1MODE_FILTER

!=======================================================================

SUBROUTINE P2MODE_FILTER(nt,outputcad,nx,xlength,p2mode)

    implicit none
    integer, intent(in) :: nt,nx
    real*8, intent(in) :: xlength,outputcad
    integer i,j
    real*8, intent(out) :: p2mode(nx, 1, nt)
    real*8 f_low,df,k(nx),dt
    real*8 Poly(0:2),Polylow(0:2), f0(nx),w(nt),f(nt),f1(nx),d,delta
    real*8 pi
    parameter (pi=3.141592654)

    dt = outputcad

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
END SUBROUTINE P2MODE_FILTER

!=======================================================================

SUBROUTINE P3MODE_FILTER(nt,outputcad, nx,xlength,pmode)

    implicit none
    integer, intent(in) :: nt,nx
    real*8, intent(in) :: xlength,outputcad
    integer i,j
    real*8, intent(out) :: pmode(nx, 1, nt)
    real*8 f_low,df,k(nx),dt
    real*8 Poly(0:2),Polylow(0:2), f0(nx),w(nt),f(nt),f1(nx),d,delta
    real*8 pi
    parameter (pi=3.141592654)

    dt = outputcad

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
    real*8 f_low,df,k(nx),dt
    real*8 Poly(0:4), f0(nx),w(nt),f(nt),f1(nx),d,delta,Polylow(0:4)
    real*8 pi
    parameter (pi=3.141592654)

    dt = outputcad

    Poly(0) = 1.41275585491
    Poly(1) = 9.54809436998
    Poly(2) = -10.5814674886
    Poly(3) = 7.99827826844
    Poly(4) = -2.42768573272

    Polylow(0) = 1.25437276419
    Polylow(1) = 8.13839040116
    Polylow(2) = -7.73561854055
    Polylow(3) = 4.96643235694
    Polylow(4) = -1.25914661289

    f_low = 1.6
    df = 0.5

    call distmat(nx,1,k)
    call distmat(nt,1,w)
    k = abs(k) * 2.*pi/(xlength *nx/(nx-1.))
    w = abs(w) * 2.*pi/(nt*dt)

    f0=0.
    f1=0.
    do i=0,size(Polylow)-1
        f0=f0+Polylow(i)*k**i
    end do
    do i=0,size(Poly)-1
        f1=f1+Poly(i)*k**i
    end do

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
    real*8 f_low,df,k(nx),dt
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
    real*8 f_low,df,k(nx),dt
    real*8 Poly(0:2),Polylow(0:2), f0(nx),w(nt),f(nt),f1(nx),d,delta
    real*8 pi
    parameter (pi=3.141592654)

    dt = outputcad

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
    real*8 f_low,df,k(nx),dt
    real*8 Poly(0:4),Polylow(0:4), f0(nx),w(nt),f(nt),f1(nx),d,delta
    real*8 pi
    parameter (pi=3.141592654)

    dt = outputcad

    Poly(0) = 1.6216558409
    Poly(1) = 13.2697412028
    Poly(2) = -14.8197349148
    Poly(3) = 9.44009291408
    Poly(4) = -2.17165011845

    Polylow(0) = 1.69587030179
    Polylow(1) = 10.4058096492
    Polylow(2) = -9.1004586514
    Polylow(3) = 5.14732934967
    Polylow(4) = -1.19380905259

    f_low = 1.6
    df = 0.5

    call distmat(nx,1,k)
    call distmat(nt,1,w)
    k = abs(k) * 2.*pi/(xlength *nx/(nx-1.))
    w = abs(w) * 2.*pi/(nt*dt)

    f0=0.
    f1=0.
    do i=0,size(Polylow)-1
        f0=f0+Polylow(i)*k**i
    end do
    do i=0,size(Poly)-1
        f1=f1+Poly(i)*k**i
    end do
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

SUBROUTINE ALL_PMODE_FILTER(nt,outputcad,nx,xlength,pmode)

  implicit none
    integer, intent(in) :: nt,nx
    real*8, intent(in) :: xlength,outputcad
    integer i,j
    real*8, intent(out) :: pmode(nx, 1, nt)
    real*8 f_low,df,k(nx),dt
    real*8 Poly(0:2),Polylow(0:2), f0(nx),w(nt),f(nt),f1(nx),d,delta
    real*8 pi
    parameter (pi=3.141592654)

    dt=outputcad


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

    pmode = 0.0
    do i=1,nx
    delta = f1(i) - f0(i)

    do j=1,nt
     d = f(j) - f0(i)
     if ((d .lt. delta) .and. (d .gt. 0)) then
        pmode(i,1,j) = 1.
     end if
    enddo
    enddo

    do j=1,nt
    if (f(j) .lt. f_low) pmode(:,1,j) = 0.
    if (f(j) .lt. f_low+df) &
      pmode(:,1,j) = pmode(:,1,j)*0.5*(1.+cos(pi*(f(j)-(f_low+df))/df) )
    enddo
END SUBROUTINE ALL_PMODE_FILTER
!================================================================================

SUBROUTINE FIRST_BOUNCE_FILTER(nt,dt,nx,Lx,srcloc,fb)

    implicit none
    integer, intent(in) :: nt,nx
    real*8, intent(in) :: Lx,dt,srcloc
    real*8, intent(out) :: fb(nx, 1, nt)
    real*8 dist,xcoord
    integer i,timest,timefin

    ! dt in seconds, in driver dt is in minutes

    fb=0.
    do i=1,nx

    xcoord = dble(i)*Lx/nx - Lx/2.
    dist=abs(xcoord - srcloc)


    timest = floor((-0.000417147774671*dist**2+0.313350998096*dist+17.9609631186)* 1./(dt/60))
    timefin = floor((-0.00028361249034*dist**2+0.29270337114*dist+36.4398734578)* 1./(dt/60))

    fb(i,1,timest:timefin) = 1.

    end do
END SUBROUTINE FIRST_BOUNCE_FILTER

SUBROUTINE HIGHPMODE_FILTER(nt,outputcad,nx,xlength,pmode)
    implicit none
    integer, intent(in) :: nt,nx
    real*8, intent(in) :: xlength,outputcad
    integer i,j
    real*8, intent(out) :: pmode(nx, 1, nt)
    real*8 k(nx),dt
    real*8 Polylow(0:2),f0(nx),w(nt),f(nt)
    real*8 pi
    parameter (pi=3.141592654)

    dt = outputcad

    Polylow(0)=2.9
    Polylow(1)=5.9
    Polylow(2)=-1.3

    call distmat(nx,1,k)
    call distmat(nt,1,w)
    k = abs(k) * 2.*pi/(xlength *nx/(nx-1.))
    w = abs(w) * 2.*pi/(nt*dt)

    f0=Polylow(0) + Polylow(1)*k +Polylow(2)*k**2.
    f = w/(2.*pi)*1e3

    pmode = 0.0
    do i=1,nx
    do j=1,nt
    pmode(i,1,j) = 1./(1.+ exp(-(f(j) - f0(i))/0.5))
    enddo
    enddo
END SUBROUTINE HIGHPMODE_FILTER

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
