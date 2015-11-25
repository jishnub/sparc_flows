program PMODE_FILTER
  implicit none
 
  integer i,j,nt,nx
  parameter (nt=376,nx=512)
  real*8 pmode(nx, 1, nt),f_low,df,k(nx),dt,f_mode_const,pi
  real*8 Poly(0:2), f0(nx),w(nt),f(nt),f1(nx),d,delta,xlength
  
  pi=acos(dble(-1.0))
 
  dt = 30.0
  xlength=800.0 * 10**(8)
  
  open(unit=123,file='pmode_filter_parameters',ACTION='read')

    read(123,*) f_mode_const
    read(123,*) Poly(0)
    read(123,*) Poly(1)
    read(123,*) Poly(2)
    read(123,*) f_low

  df = 0.5
  
  close(123)
  
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
   

    call writefits_3d('pmode_filter.fits',pmode,nx, 1, nt)

  
   
END program PMODE_FILTER

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
