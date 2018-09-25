c**********************************************************
c
c       Copyright 1994, Neal Hurlburt, 
c  Lockheed Palo Alto Research Laboratory
c
C***************************************************************************
      SUBROUTINE MTRIDAG(A,B,C,R,U,N1,N2,ISPECIAL)
C
C     This solves a set of tridiagonal systems in parallel
c     Adapted from Press etal. 
c     Revised 3/17/98 for periodic boundaries -- NEH
c     fixed to cope with zero pivots like tridag 6/16/99 -- NEH
C**************************************************************************
      PARAMETER (NMAX=7000)
      integer, parameter :: real64 = kind(1.d0)
      real(kind=real64) A(N2),B(N2),C(N2),R(N1,N2),U(N1,N2)
      real(kind=real64) GAM(NMAX),UD(NMAX),BET
      real(kind=real64) betas(nmax),betau(nmax),betad(nmax)
      real(kind=real64) alpha(nmax),alphab(nmax)
      if (ispecial.eq.1) then
c
c     prepare decomposed matrix
c
         do i=1,n2
            betad(i)=b(i)
            betas(i)=0
            betau(i)=c(i)
            alpha(i)=a(i)
            alphab(i)=0
         enddo
         betas(1)=a(1)
         betas(n2-1)=c(n2-1)
         alphab(1)=c(n2)
         alphab(n2-1)=a(n2)
c
         do i=2,n2-1
            betad(i)=betad(i)-alpha(i-1)*betau(i-1)
            betas(i)=betas(i)-alpha(i-1)*betas(i-1)
            alpha(i)=alpha(i)/betad(i)
            alphab(i)=(alphab(i)-alphab(i-1)*betau(i-1))/betad(i)
         enddo
         do i=1,n2-1
            betad(n2)=betad(n2)-alphab(i)*betas(i)
         enddo
c
c and now do the forward & backward substitution
c
         DO I=1,N1
            ud(1)=r(I,1)
            do j=2,n2-1
               ud(j)=r(i,j)-alpha(j-1)*ud(j-1)
            enddo
            ud(n2)=r(I,n2)
            do j=1,n2-1
               ud(n2)=ud(n2)-alphab(j)*ud(j)
            enddo
            ud(n2)=ud(n2)/betad(n2)
            ud(n2-1)=(ud(n2-1)-betau(n2-1)*ud(n2))/betad(n2-1)
            do j=n2-2,1,-1
               ud(j)=(ud(j)-betau(j)*ud(j+1)-betas(j)*ud(n2))/betad(j)
            enddo
            do j=1,n2
               u(i,j)=ud(j)
            enddo
         enddo
c
      else if (ispecial.eq.2) then
c     Coping with Zero pivots on boundary
c
c do the forward & backward substitution
c
         DO I=1,N1
            ud(1)=r(i,1)/b(1)
            betas(1)=1

            BET=B(2)
            UD(2)=R(i,2)/BET
            gam(2)=c(1)/BET
            betas(2)=a(2)/BET
            DO J=3,N2
               GAM(J)=C(J-1)/BET
               BET=B(J)-A(J)*GAM(J)
               UD(J)=(R(i,J)-A(J)*UD(J-1))/BET
               betas(j)=-betas(j-1)*a(j)/BET
            enddo
            DO  J=N2-1,1,-1
               UD(J)=UD(J)-GAM(J+1)*UD(J+1)
               betas(j)=betas(j)-gam(j+1)*betas(j+1)
            enddo
            U(i,1)=UD(1)/betas(1)
            DO J=2,N2
               U(i,J)=UD(J)-u(i,1)*betas(j)
            enddo
         enddo
      else
         BETAD(1)=B(1)
         DO 11 J=2,N2
            GAM(J)=C(J-1)/BETAD(J-1)
            BETAD(J)=B(J)-A(J)*GAM(J)
 11      CONTINUE
         DO 100 I=1,N1
            UD(1)=R(I,1)/BETAD(1)
            DO 300 J=2,N2
               UD(J)=(R(I,J)-A(J)*UD(J-1))/BETAD(J)
 300        CONTINUE
            DO 12 J=N2-1,1,-1
               UD(J)=UD(J)-GAM(J+1)*UD(J+1)
 12         CONTINUE
            DO 200 J=1,N2
               U(I,J)=UD(J)
 200        CONTINUE
 100     CONTINUE
      endif
      RETURN
      END
