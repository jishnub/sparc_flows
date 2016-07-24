C***************************************************************************
      SUBROUTINE TRIDAG(A,B,C,R,U,N,ISPECIAL)
c     Adapted from Press etal. 
c     Revised 3/17/98 for periodic boundaries -- NEH
C**************************************************************************
      PARAMETER (NMAX=7000)
      REAL*8 A(N),B(N),C(N),R(N),U(N)
      REAL*8 GAM(NMAX),UD(NMAX),BET
      REAL*8 betas(nmax),betau(nmax),betad(nmax)
      REAL*8 alpha(nmax),alphab(nmax)
      if (ispecial.eq.1) then
c
c     Periodic conditions
c
c     prepare decomposed matrix
c
         do i=1,n
            betad(i)=b(i)
            betas(i)=0
            betau(i)=c(i)
            alpha(i)=a(i)
            alphab(i)=0
         enddo
         betas(1)=a(1)
         betas(n-1)=c(n-1)
         alphab(1)=c(n)
         alphab(n-1)=a(n)
c
         do i=2,n-1
            betad(i)=betad(i)-alpha(i-1)*betau(i-1)
            betas(i)=betas(i)-alpha(i-1)*betas(i-1)
            alpha(i)=alpha(i)/betad(i)
            alphab(i)=(alphab(i)-alphab(i-1)*betau(i-1))/betad(i)
         enddo
         do i=1,n-1
            betad(n)=betad(n)-alphab(i)*betas(i)
         enddo
c
c and now do the forward & backward substitution
c
         ud(1)=r(1)
         do i=2,n-1
            ud(i)=r(i)-alpha(i-1)*ud(i-1)
         enddo
         ud(n)=r(n)
         do i=1,n-1
            ud(n)=ud(n)-alphab(i)*ud(i)
         enddo
         ud(n)=ud(n)/betad(n)
         ud(n-1)=(ud(n-1)-betau(n-1)*ud(n))/betad(n-1)
         do i=n-2,1,-1
            ud(i)=(ud(i)-betau(i)*ud(i+1)-betas(i)*ud(n))/betad(i)
         enddo
         do i=1,n
            u(i)=ud(i)
         enddo
c
      else if (ispecial.eq.2) then
c     Coping with Zero pivots on boundary
c
c do the forward & backward substitution
c
         ud(1)=r(1)/b(1)
         betas(1)=1

         BET=B(2)
         UD(2)=R(2)/BET
         gam(2)=c(1)/BET
         betas(2)=a(2)/BET
         DO J=3,N
            GAM(J)=C(J-1)/BET
            BET=B(J)-A(J)*GAM(J)
            UD(J)=(R(J)-A(J)*UD(J-1))/BET
            betas(j)=-betas(j-1)*a(j)/BET
         enddo
         DO  J=N-1,1,-1
            UD(J)=UD(J)-GAM(J+1)*UD(J+1)
            betas(j)=betas(j)-gam(j+1)*betas(j+1)
         enddo
         U(1)=UD(1)/betas(1)
         DO J=2,N
            U(J)=UD(J)-u(1)*betas(j)
         enddo
      else
         BET=B(1)
         UD(1)=R(1)/BET
         DO 11 J=2,N
            GAM(J)=C(J-1)/BET
            BET=B(J)-A(J)*GAM(J)
            UD(J)=(R(J)-A(J)*UD(J-1))/BET
 11      CONTINUE
         DO 12 J=N-1,1,-1
            UD(J)=UD(J)-GAM(J+1)*UD(J+1)
 12      CONTINUE
         DO 20 J=1,N
            U(J)=UD(J)
 20      CONTINUE
      endif
      RETURN
      END
