MODULE INTEGRALS

use initialize
!~ use spline
!~ use Gaussm3
contains

REAL*8 FUNCTION INTEGRATE_Z(f)

    implicit none
    real*8, intent(in) :: f(nx,1,nz)
    real*8 :: INTEGRATE_Z(nx,1,nz)
    real*8, dimension(nx,1,nz) :: int_f_spl,int_f_tpz
    real*8 int_spl,int_tpz
    real*8 splncoeff(3,nz)
    integer splerr,xdum,zdum,abcdum,abcdum2
    
    logical use_splines,use_antia
    
    use_splines = .False.
    use_antia = .False.
    
    int_f_spl=0.0
    int_f_tpz=0.0
    
    do xdum=1,nx
        
        if (use_antia) then
            call spline(z,f(xdum,1,:),nz,splncoeff,splerr)
            if (splerr .ne. 0) print *,"Error code in spline",splerr
            do abcdum=1,nz
                do abcdum2=1,3
                    if (isnan(splncoeff(abcdum2,abcdum))) then
                        print *, "nan found in coeff, x=",xdum
                    end if
                end do
            end do
            do zdum=2,240

                call splint(z(1),z(zdum),int_spl,int_tpz,zdum,z(1:zdum),&
                            f(xdum,1,1:zdum),splncoeff(:,1:zdum),splerr)
                            
                int_f_spl(xdum,1,zdum) = int_spl
                int_f_tpz(xdum,1,zdum) = int_tpz
                
                if (splerr .ne. 0) then
                    print *,"Error code in splint",splerr
                    print *,"Integrating from",z(1),"to",z(zdum),"for x=",xdum
                end if
            
            end do
        else
            int_f_tpz(xdum,1,:) = trapz_irregular(f(xdum,1,:),z)
        end if
    end do
    
    if (use_antia .and. use_splines) then
        INTEGRATE_Z = int_f_spl
    else
        INTEGRATE_Z = int_f_tpz
    endif
    

END FUNCTION INTEGRATE_Z


SUBROUTINE INTEGRATE_TIME(f,int_f,dt)

    implicit none
    real*8, intent(in) :: dt,f(:)
    real*8, intent(out) ::  int_f
    
    int_f=simpson_regular(f,dt)
    

END SUBROUTINE INTEGRATE_TIME

REAL*8 FUNCTION TRAPZ_IRREGULAR(f,z)
    implicit none
    real*8, intent(in),dimension(nz) :: f,z
    real*8, dimension(nz) :: int_f
    real*8 TRAPZ_IRREGULAR(nz)
    real*8 temp
    integer k,l
    
    int_f=0.0
    do k=2,240
        do l=1,k-1
            temp = 0.5*(f(l+1)+f(l))*(z(l+1)-z(l))
            int_f(k) = int_f(k) + temp
        end do
    end do
    
    TRAPZ_IRREGULAR = int_f
    
END FUNCTION TRAPZ_IRREGULAR

REAL*8 FUNCTION SIMPSON_REGULAR(f,dvar)

    implicit none
    real*8, intent(in) :: dvar,f(:)
    integer k,nt
    real*8 int_f
    
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
    
    SIMPSON_REGULAR=int_f*dvar/3.0

END FUNCTION SIMPSON_REGULAR

END MODULE INTEGRALS
