SUBROUTINE FMODE_FILTER(nt,outputcad,nx,xlength,fmode)
    implicit none
    integer, intent(in) :: nt,nx
    real*8, intent(in) :: xlength,outputcad
    integer i,j
    real*8, intent(out) :: fmode(nx, 1, nt)
    real*8 f_low_cutoff,df,dt
    real*8 :: f_low(nx),f_high(nx),k(nx)
    real*8 :: Poly(0:4),Polylow(0:4), w(nt),d,delta,f(nt)
    real*8 pi
    parameter (pi=3.141592654)

    dt = outputcad

    ! Poly(0)=1.1
    ! Poly(1)=1.9
    ! Poly(2)=-0.2
    !
    ! Polylow(0)=0.7
    ! Polylow(1)=1.7
    ! Polylow(2)=-0.2

    Poly=0
    ! Poly(0)=0.6816995666
    ! Poly(1)=1.97827975135
    ! Poly(2)=-1.12337555632
    ! Poly(3)=0.749996956203
    ! Poly(4)=-0.175140581889

    Poly(0) = 0.783876582013
    Poly(1) = 3.15585935078
    Poly(2) = -1.23549169809
    Poly(3) = 0.320272543559
    Poly(4) = -0.0410749307678

    Polylow=0
    ! Polylow(0)=0.623898151707
    ! Polylow(1)=3.81510329344
    ! Polylow(2)=-2.18073803978
    ! Polylow(3)=0.83870962992
    ! Polylow(4)=-0.126283442819

    Polylow(0) = 0.646451297265
    Polylow(1) = 1.71728869643
    Polylow(2) = -0.0180473695549
    Polylow(3) = -0.241491147354
    Polylow(4) = 0.0745045635964

    call distmat(nx,1,k)
    call distmat(nt,1,w)
    k = abs(k) * 2.*pi/xlength
    w = abs(w) * 2.*pi/(nt*dt)

    f_low=0
    f_high=0

    do i=0,4
    f_low=f_low+Polylow(i)*k**i
    f_high=f_high+Poly(i)*k**i
    enddo

    f = w/(2.*pi)*1e3


    ! fmode = 0.0
    ! do i=1,nx
    ! delta = (f_high(i) - f_low(i))
    ! do j=1,nt
    !  d = f(j) - f_low(i)
    !  if ((d .lt. delta) .and. (d>0)) &
    !     fmode(i,1,j) = 0.5*(1.+cos(pi*(2.0*d/delta-1)))
    ! enddo
    ! enddo

    fmode = 0.0
    df = f(2) - f(1)
    do i=1,nx
    delta = (f_high(i) - f_low(i))
    do j=1,nt
     d = f(j) - f_low(i)
     if ((d < (delta+2*df)) .and. (d>delta)) then
        fmode(i,1,j) = cos(2*pi/(8*df)*(f(j)-f_high(i)))
    else if ((d < 0) .and. (d>(-2*df))) then
        fmode(i,1,j) = cos(2*pi/(8*df)*(f(j)-f_low(i)))
     elseif ((d > 0 ) .and. (d<delta)) then
       fmode(i,1,j) = 1.0
     endif
    enddo
    enddo

    f_low_cutoff = 1.1
    do j=1,nt
    if (f(j) .lt. f_low_cutoff+0.5) &
      fmode(:,1,j) = fmode(:,1,j)*0.5*(1.+cos(pi*(f(j)-(f_low_cutoff+0.5))/0.5) )
    if (f(j) .lt. f_low_cutoff) fmode(:,1,j) = 0.
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
    real*8 f_low_cutoff,df,dt
    real*8 :: f_low(nx),f_high(nx),k(nx)
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

    f_low_cutoff = 1.1

    call distmat(nx,1,k)
    call distmat(nt,1,w)
    k = abs(k) * 2.*pi/xlength
    w = abs(w) * 2.*pi/(nt*dt)

    f_low=Polylow(0) + Polylow(1)*k +Polylow(2)*k**2.
    f_high=Poly(0) + Poly(1)*k +Poly(2)*k**2.

    f = w/(2.*pi)*1e3

    df = f(2)-f(1)

    fmode = 0.0
    do i=1,nx
    delta = (f_high(i) - f_low(i))
    do j=1,nt
     d = f(j) - f_low(i)
     if ((d .lt. delta) .and. (d>(delta-2*df))) then
        fmode(i,1,j) = cos(2*pi/(8*df)*(f(j)-(f_high(i)-2*df)))
    else if ((d .lt. (2*df)) .and. (d>0)) then
        fmode(i,1,j) = cos(2*pi/(8*df)*(f(j)-(f_low(i)+2*df)))
     elseif ((d > 2*df ) .and. (d<(delta-2*df))) then
       fmode(i,1,j) = 1.0
     endif
    enddo
    enddo

    do j=1,nt
    if (f(j) .lt. f_low_cutoff+0.5) &
      fmode(:,1,j) = fmode(:,1,j)*0.5*(1.+cos(pi*(f(j)-(f_low_cutoff+0.5))/0.5) )
    if (f(j) .lt. f_low_cutoff) fmode(:,1,j) = 0.
    enddo
END SUBROUTINE FMODE_FILTER_ALT

!=======================================================================

SUBROUTINE P1MODE_FILTER(nt,outputcad,nx,xlength,pmode)

    implicit none
    integer, intent(in) :: nt,nx
    real*8, intent(in) :: xlength,outputcad
    integer i,j
    real*8, intent(out) :: pmode(nx, 1, nt)
    real*8 f_low_cutoff,df,k(nx),dt
    real*8 Poly(0:4), f_low(nx),w(nt),f(nt),f_high(nx),d,delta, Polylow(0:4)
    real*8 pi
    parameter (pi=3.141592654)

    dt = outputcad

    ! Poly(0)=1.4
    ! Poly(1)=3.0
    ! Poly(2)=-0.5
    !
    ! Polylow(0)=1.1
    ! Polylow(1)=2.4
    ! Polylow(2)=-0.3

    Poly=0
    Poly(0)=1.04946787438
    Poly(1)=4.28500694767
    Poly(2)=-1.79853784214
    Poly(3)=0.523240395536
    Poly(4)=-0.0693528244107

    Polylow = 0
    Polylow(0)=0.88123404174
    Polylow(1)=3.02499866767
    Polylow(2)=-1.36638485558
    Polylow(3)=0.688118939264
    Polylow(4)=-0.16113216593

    call distmat(nx,1,k)
    call distmat(nt,1,w)
    k = abs(k) * 2.*pi/(xlength*nx/(nx-1.))
    w = abs(w) * 2.*pi/(nt*dt)

    f_low=0
    f_high=0

    do i=0,4
    f_low=f_low+Polylow(i)*k**i
    f_high=f_high+Poly(i)*k**i
    enddo

    f = w/(2.*pi)*1e3

    pmode = 0.0
    df = f(2) - f(1)
    do i=1,nx
    delta = (f_high(i) - f_low(i))

    do j=1,nt
     d = f(j) - f_low(i)
     if ((d <(delta+2*df)) .and. (d>delta)) then
        pmode(i,1,j) = cos(2*pi/(8*df)*(f(j)-f_high(i)))
     elseif ((d < 0) .and. (d>(-2*df))) then
        pmode(i,1,j) = cos(2*pi/(8*df)*(f(j)-f_low(i)))
     elseif ((d > 0 ) .and. (d<delta)) then
       pmode(i,1,j) = 1.0
     endif
    enddo
    enddo
    ! pmode = 0.0
    ! do i=1,nx
    ! delta = (f_high(i) - f_low(i))
    ! do j=1,nt
    !  d = f(j) - f_low(i)
    !  if ((d .lt. delta) .and. (d .gt. 0)) then
    !     pmode(i,1,j) = 0.5*(1.+cos(pi*(2.0*d/delta-1)))
    !  end if
    ! enddo
    ! enddo

    f_low_cutoff = 1.1
     do j=1,nt
     if (f(j) .lt. f_low_cutoff) pmode(:,1,j) = 0.
     if (f(j) .lt. f_low_cutoff+0.5) &
       pmode(:,1,j) = pmode(:,1,j)*0.5*(1.+cos(pi*(f(j)-(f_low_cutoff+0.5))/0.5) )
     enddo
END SUBROUTINE P1MODE_FILTER

!=======================================================================

SUBROUTINE P2MODE_FILTER(nt,outputcad,nx,xlength,p2mode)

    implicit none
    integer, intent(in) :: nt,nx
    real*8, intent(in) :: xlength,outputcad
    integer i,j
    real*8, intent(out) :: p2mode(nx, 1, nt)
    real*8 f_low_cutoff,df,k(nx),dt
    real*8 Poly(0:4),Polylow(0:4), f_low(nx),w(nt),f(nt),f_high(nx),d,delta
    real*8 pi
    parameter (pi=3.141592654)

    dt = outputcad

    ! Poly(0)=1.6
    ! Poly(1)=3.8
    ! Poly(2)=-0.65
    !
    ! Polylow(0)=1.4
    ! Polylow(1)=3.3
    ! Polylow(2)=-0.62

    Poly=0
    ! Poly(0)=1.1406600853
    ! Poly(1)=6.17989298938
    ! Poly(2)=-3.90117798051
    ! Poly(3)=1.50707755921
    ! Poly(4)=-0.189573049963

    Poly(0) = 1.10604110752
    Poly(1) = 6.24065357264
    Poly(2) = -4.20067393642
    Poly(3) = 1.99640927807
    Poly(4) = -0.401627562953

    Polylow = 0
    Polylow(0)=0.91217205529
    Polylow(1)=5.16466338364
    Polylow(2)=-3.01937028923
    Polylow(3)=1.11497955424
    Polylow(4)=-0.14466262072

    f_low_cutoff = 1.1
    df = 0.5

    call distmat(nx,1,k)
    call distmat(nt,1,w)
    k = abs(k) * 2.*pi/(xlength *nx/(nx-1.))
    w = abs(w) * 2.*pi/(nt*dt)


    f_low=0
    f_high=0

    do i=0,4
     f_low=f_low+Polylow(i)*k**i
      f_high=f_high+Poly(i)*k**i
    enddo
    f = w/(2.*pi)*1e3

    p2mode = 0.0
    df = f(2) - f(1)
    do i=1,nx
    delta = (f_high(i) - f_low(i))

    do j=1,nt
     d = f(j) - f_low(i)
     if ((d <(delta+2*df)) .and. (d>delta)) then
        p2mode(i,1,j) = cos(2*pi/(8*df)*(f(j)-f_high(i)))
     elseif ((d < 0) .and. (d>(-2*df))) then
        p2mode(i,1,j) = cos(2*pi/(8*df)*(f(j)-f_low(i)))
     elseif ((d > 0 ) .and. (d<delta)) then
       p2mode(i,1,j) = 1.0
     endif
    enddo
    enddo

    ! p2mode = 0.0
    ! do i=1,nx
    ! delta = (f_high(i) - f_low(i))
    ! do j=1,nt
    !  d = f(j) - f_low(i)
    !  if ((d .lt. delta) .and. (d .gt. 0)) then
    !     p2mode(i,1,j) = 0.5*(1.+cos(pi*(2.0*d/delta-1)))
    !  end if
    ! enddo
    ! enddo

    do j=1,nt
    if (f(j) .lt. f_low_cutoff) p2mode(:,1,j) = 0.
    if (f(j) .lt. f_low_cutoff+df) &
      p2mode(:,1,j) = p2mode(:,1,j)*0.5*(1.+cos(pi*(f(j)-(f_low_cutoff+df))/df) )
    enddo
END SUBROUTINE P2MODE_FILTER

!=======================================================================

SUBROUTINE P3MODE_FILTER(nt,outputcad, nx,xlength,pmode)

    implicit none
    integer, intent(in) :: nt,nx
    real*8, intent(in) :: xlength,outputcad
    integer i,j
    real*8, intent(out) :: pmode(nx, 1, nt)
    real*8 f_low_cutoff,df,k(nx),dt
    real*8 Poly(0:4),Polylow(0:4), f_low(nx),w(nt),f(nt),f_high(nx),d,delta
    real*8 pi
    parameter (pi=3.141592654)

    dt = outputcad


    Poly=0
    ! Poly(0)=1.45
    ! Poly(1)=6.87
    ! Poly(2)=-5.325
    ! Poly(3)=3.192
    ! Poly(4)=-0.847

    Poly(0) = 1.32385750705
    Poly(1) = 7.0547233868
    Poly(2) = -4.66879987086
    Poly(3) = 2.10460128455
    Poly(4) = -0.432244695425

    Polylow=0
    ! Polylow(0)=1.214
    ! Polylow(1)=5.906
    ! Polylow(2)=-3.889
    ! Polylow(3)=1.9
    ! Polylow(4)=-0.387
    Polylow(0) = 1.10604110752
    Polylow(1) = 6.24065357264
    Polylow(2) = -4.20067393642
    Polylow(3) = 1.99640927807
    Polylow(4) = -0.401627562953

    ! Poly(0)=2
    ! Poly(1)=4.1
    ! Poly(2)=-0.8
    !
    ! Polylow(0)=2
    ! Polylow(1)=3.55
    ! Polylow(2)=-0.7

    f_low_cutoff = 1.1
    df = 0.5

    call distmat(nx,1,k)
    call distmat(nt,1,w)
    k = abs(k) * 2.*pi/(xlength *nx/(nx-1.))
    w = abs(w) * 2.*pi/(nt*dt)

    f_low=0
    f_high=0

    do i=0,4
     f_low=f_low+Polylow(i)*k**i
      f_high=f_high+Poly(i)*k**i
    enddo
    f = w/(2.*pi)*1e3

    pmode = 0.0

    df = f(2) - f(1)
    do i=1,nx
    delta = (f_high(i) - f_low(i))
    do j=1,nt
     d = f(j) - f_low(i)
     if ((d <(delta+2*df)) .and. (d>delta)) then
        pmode(i,1,j) = cos(2*pi/(8*df)*(f(j)-f_high(i)))
     elseif ((d < 0) .and. (d>(-2*df))) then
        pmode(i,1,j) = cos(2*pi/(8*df)*(f(j)-f_low(i)))
     elseif ((d > 0 ) .and. (d<delta)) then
       pmode(i,1,j) = 1.0
     endif
    enddo
    enddo

    ! pmode = 0.0
    ! do i=1,nx
    !     delta = (f_high(i) - f_low(i))
    !     do j=1,nt
    !         d = f(j) - f_low(i)
    !         if ((d .lt. delta) .and. (d .gt. 0)) then
    !             pmode(i,1,j) = 0.5*(1.+cos(pi*(abs(d)-abs(delta)*0.5)/(abs(delta)*0.5)))
    !         end if
    !     enddo
    ! enddo

    do j=1,nt
        if (f(j) .lt. f_low_cutoff) pmode(:,1,j) = 0.
        if (f(j) .lt. f_low_cutoff+df) &
          pmode(:,1,j) = pmode(:,1,j)*0.5*(1.+cos(pi*(f(j)-(f_low_cutoff+df))/df) )
    enddo
END SUBROUTINE P3MODE_FILTER

!=======================================================================


SUBROUTINE P4MODE_FILTER(nt,outputcad,nx,xlength,pmode)

    implicit none
    integer, intent(in) :: nt,nx
    real*8, intent(in) :: xlength,outputcad
    integer i,j
    real*8, intent(out) :: pmode(nx, 1, nt)
    real*8 f_low_cutoff,df,k(nx),dt
    real*8 Poly(0:4), f_low(nx),w(nt),f(nt),f_high(nx),d,delta,Polylow(0:4)
    real*8 pi
    parameter (pi=3.141592654)

    dt = outputcad

    Poly = 0
    ! Poly(0)=1.31049991687
    ! Poly(1)=9.42269860347
    ! Poly(2)=-8.93928345177
    ! Poly(3)= 5.61901410549
    ! Poly(4)=-1.45560970906

    Poly(0) = 1.42875410041
    Poly(1) = 8.51590951021
    Poly(2) = -6.74658351266
    Poly(3) = 3.37874845461
    Poly(4) = -0.632464159464

    Polylow = 0
    ! Polylow(0)=1.20256877128
    ! Polylow(1)=8.40378321098
    ! Polylow(2)=-8.49203154018
    ! Polylow(3)=6.04350284182
    ! Polylow(4)=-1.75134643804

    Polylow(0) = 1.29998031988
    Polylow(1) = 7.32973515885
    Polylow(2) = -5.1628850738
    Polylow(3) = 2.45226315675
    Polylow(4) = -0.507544900984

    ! Poly(0) = 1.41275585491
    ! Poly(1) = 9.54809436998
    ! Poly(2) = -10.5814674886
    ! Poly(3) = 7.99827826844
    ! Poly(4) = -2.42768573272
    !
    ! Polylow(0) = 1.25437276419
    ! Polylow(1) = 8.13839040116
    ! Polylow(2) = -7.73561854055
    ! Polylow(3) = 4.96643235694
    ! Polylow(4) = -1.25914661289

    f_low_cutoff = 1.1
    df = 0.5

    call distmat(nx,1,k)
    call distmat(nt,1,w)
    k = abs(k) * 2.*pi/(xlength *nx/(nx-1.))
    w = abs(w) * 2.*pi/(nt*dt)

    f_low=0
    f_high=0

    do i=0,4
     f_low=f_low+Polylow(i)*k**i
      f_high=f_high+Poly(i)*k**i
    enddo

    f = w/(2.*pi)*1e3

    pmode = 0.0
    df = f(2) - f(1)

    do i=1,nx
    delta = (f_high(i) - f_low(i))
    do j=1,nt
     d = f(j) - f_low(i)
     if ((d <(delta+2*df)) .and. (d>delta)) then
        pmode(i,1,j) = cos(2*pi/(8*df)*(f(j)-f_high(i)))
     elseif ((d < 0) .and. (d>(-2*df))) then
        pmode(i,1,j) = cos(2*pi/(8*df)*(f(j)-f_low(i)))
     elseif ((d > 0 ) .and. (d<delta)) then
       pmode(i,1,j) = 1.0
     endif
    enddo
    enddo

    ! pmode = 0.0
    ! do i=1,nx
    ! delta = (f_high(i) - f_low(i))
    ! do j=1,nt
    !  d = f(j) - f_low(i)
    !  if ((d .lt. delta) .and. (d .gt. 0)) then
    !     pmode(i,1,j) = 0.5*(1.+cos(pi*(abs(d)-abs(delta)*0.5)/(abs(delta)*0.5)))
    !  end if
    ! enddo
    ! enddo

    do j=1,nt
    if (f(j) .lt. f_low_cutoff) pmode(:,1,j) = 0.
    if (f(j) .lt. f_low_cutoff+df) &
      pmode(:,1,j) = pmode(:,1,j)*0.5*(1.+cos(pi*(f(j)-(f_low_cutoff+df))/df) )
    enddo
END SUBROUTINE P4MODE_FILTER

!=======================================================================


SUBROUTINE P5MODE_FILTER(nt,outputcad,nx,xlength,pmode)

    implicit none
    integer, intent(in) :: nt,nx
    real*8, intent(in) :: xlength,outputcad
    integer i,j
    real*8, intent(out) :: pmode(nx, 1, nt)
    real*8 f_low_cutoff,df,k(nx),dt
    real*8 Poly(0:2),Polylow(0:2), f_low(nx),w(nt),f(nt),f_high(nx),d,delta
    real*8 pi
    parameter (pi=3.141592654)

    dt = outputcad


    Poly(0)=2.35
    Poly(1)=5.6
    Poly(2)=-1.1

    Polylow(0)=2.2
    Polylow(1)=4.7
    Polylow(2)=-1.0

    f_low_cutoff = 1.6
    df = 0.5

    call distmat(nx,1,k)
    call distmat(nt,1,w)
    k = abs(k) * 2.*pi/(xlength *nx/(nx-1.))
    w = abs(w) * 2.*pi/(nt*dt)

    f_low=Polylow(0) + Polylow(1)*k +Polylow(2)*k**2.
    f_high=Poly(0) + Poly(1)*k +Poly(2)*k**2.
    f = w/(2.*pi)*1e3

    pmode = 0.0
    do i=1,nx
    delta = (f_high(i) - f_low(i))
    do j=1,nt
     d = f(j) - f_low(i)
     if ((d .lt. delta) .and. (d .gt. 0)) then
        pmode(i,1,j) = 0.5*(1.+cos(pi*(abs(d)-abs(delta)*0.5)/(abs(delta)*0.5)))
     end if
    enddo
    enddo

    do j=1,nt
    if (f(j) .lt. f_low_cutoff) pmode(:,1,j) = 0.
    if (f(j) .lt. f_low_cutoff+df) &
      pmode(:,1,j) = pmode(:,1,j)*0.5*(1.+cos(pi*(f(j)-(f_low_cutoff+df))/df) )
    enddo
END SUBROUTINE P5MODE_FILTER

!=======================================================================


SUBROUTINE P6MODE_FILTER(nt,outputcad,nx,xlength,pmode)

    implicit none
    integer, intent(in) :: nt,nx
    real*8, intent(in) :: xlength,outputcad
    integer i,j
    real*8, intent(out) :: pmode(nx, 1, nt)
    real*8 f_low_cutoff,df,k(nx),dt
    real*8 Poly(0:2),Polylow(0:2), f_low(nx),w(nt),f(nt),f_high(nx),d,delta
    real*8 pi
    parameter (pi=3.141592654)

    dt = outputcad

    Poly(0)=1.98567115899
    Poly(1)=8.09108986838
    Poly(2)=-3.20316331815

    Polylow(0)=1.80035417224
    Polylow(1)=7.42939105658
    Polylow(2)=-2.84595764385

    f_low_cutoff = 1.6
    df = 0.5

    call distmat(nx,1,k)
    call distmat(nt,1,w)
    k = abs(k) * 2.*pi/(xlength*nx/(nx-1.))
    w = abs(w) * 2.*pi/(nt*dt)

    f_low=Polylow(0) + Polylow(1)*k +Polylow(2)*k**2.
    f_high=Poly(0) + Poly(1)*k +Poly(2)*k**2.
    f = w/(2.*pi)*1e3

    pmode = 0.0
    do i=1,nx
    delta = (f_high(i) - f_low(i))
    do j=1,nt
     d = f(j) - f_low(i)
     if ((d .lt. delta) .and. (d .gt. 0)) then
        pmode(i,1,j) = 0.5*(1.+cos(pi*(abs(d)-abs(delta)*0.5)/(abs(delta)*0.5)))
     end if
    enddo
    enddo

    do j=1,nt
    if (f(j) .lt. f_low_cutoff) pmode(:,1,j) = 0.
    if (f(j) .lt. f_low_cutoff+df) &
      pmode(:,1,j) = pmode(:,1,j)*0.5*(1.+cos(pi*(f(j)-(f_low_cutoff+df))/df) )
    enddo
END SUBROUTINE P6MODE_FILTER
!================================================================================

SUBROUTINE P7MODE_FILTER(nt,outputcad,nx,xlength,pmode)

    implicit none
    integer, intent(in) :: nt,nx
    real*8, intent(in) :: xlength,outputcad
    integer i,j
    real*8, intent(out) :: pmode(nx, 1, nt)
    real*8 f_low_cutoff,df,k(nx),dt
    real*8 Poly(0:4),Polylow(0:4), f_low(nx),w(nt),f(nt),f_high(nx),d,delta
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

    f_low_cutoff = 1.6
    df = 0.5

    call distmat(nx,1,k)
    call distmat(nt,1,w)
    k = abs(k) * 2.*pi/(xlength *nx/(nx-1.))
    w = abs(w) * 2.*pi/(nt*dt)

    f_low=0.
    f_high=0.
    do i=0,size(Polylow)-1
        f_low=f_low+Polylow(i)*k**i
    end do
    do i=0,size(Poly)-1
        f_high=f_high+Poly(i)*k**i
    end do
    f = w/(2.*pi)*1e3

    pmode = 0.0
    do i=1,nx
    delta = (f_high(i) - f_low(i))
    do j=1,nt
     d = f(j) - f_low(i)
     if ((d .lt. delta) .and. (d .gt. 0)) then
        pmode(i,1,j) = 0.5*(1.+cos(pi*(abs(d)-abs(delta)*0.5)/(abs(delta)*0.5)))
     end if
    enddo
    enddo

    do j=1,nt
    if (f(j) .lt. f_low_cutoff) pmode(:,1,j) = 0.
    if (f(j) .lt. f_low_cutoff+df) &
      pmode(:,1,j) = pmode(:,1,j)*0.5*(1.+cos(pi*(f(j)-(f_low_cutoff+df))/df) )
    enddo
END SUBROUTINE P7MODE_FILTER

!==========================================================================================

SUBROUTINE ALL_PMODE_FILTER(nt,outputcad,nx,xlength,pmode)

  implicit none
    integer, intent(in) :: nt,nx
    real*8, intent(in) :: xlength,outputcad
    integer i,j
    real*8, intent(out) :: pmode(nx, 1, nt)
    real*8 f_low_cutoff,df,k(nx),dt
    real*8 Poly(0:2),Polylow(0:2), f_low(nx),w(nt),f(nt),f_high(nx),d,delta
    real*8 pi
    parameter (pi=3.141592654)

    dt=outputcad


    Polylow(0)=1.1
    Polylow(1)=2.4
    Polylow(2)=-0.3

    Poly(0)=3.5
    Poly(1)=6.5
    Poly(2)=-1.3

    f_low_cutoff = 1.6
    df = 0.5

    call distmat(nx,1,k)
    call distmat(nt,1,w)
    k = abs(k) * 2.*pi/(xlength *nx/(nx-1.))
    w = abs(w) * 2.*pi/(nt*dt)

    f_low=Polylow(0) + Polylow(1)*k +Polylow(2)*k**2.
    f_high=Poly(0) + Poly(1)*k +Poly(2)*k**2.
    f = w/(2.*pi)*1e3

    pmode = 0.0
    do i=1,nx
    delta = f_high(i) - f_low(i)

    do j=1,nt
     d = f(j) - f_low(i)
     if ((d .lt. delta) .and. (d .gt. 0)) then
        pmode(i,1,j) = 1.
     end if
    enddo
    enddo

    do j=1,nt
    if (f(j) .lt. f_low_cutoff) pmode(:,1,j) = 0.
    if (f(j) .lt. f_low_cutoff+df) &
      pmode(:,1,j) = pmode(:,1,j)*0.5*(1.+cos(pi*(f(j)-(f_low_cutoff+df))/df) )
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
    real*8 Polylow(0:2),f_low(nx),w(nt),f(nt)
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

    f_low=Polylow(0) + Polylow(1)*k +Polylow(2)*k**2.
    f = w/(2.*pi)*1e3

    pmode = 0.0
    do i=1,nx
    do j=1,nt
    pmode(i,1,j) = 1./(1.+ exp(-(f(j) - f_low(i))/0.5))
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
