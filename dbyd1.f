c**********************************************************
c
c       Copyright 1994, Neal Hurlburt, 
c   Lockheed Palo Alto Research Laboratory
c       Revised 3/1998 for periodic conditions
c
C************************************************************
        SUBROUTINE DBYD1(A,B,N1,N2,IBC)
        PARAMETER (NMAX=7000)
        integer, parameter :: real64 = kind(1.d0)
C------------------------------------------------------------
C  This routine evaluates the derivative of array B using 
C  sixth-order compact differences and 
C  stores the result in array A. The derivative is evaluated
C  along the first argument of B.
C
C      IBC<0 returns revision level in a(1,1)
C      IBC=0 uses neumann conditions on boundaries (Dydx=0)
C      IBC=1 calculates value on boundaries
C      IBC=2 neumann on low end, Dirichelet on high
C      IBC=3 vica versa
C      IBC=4 uses neumann conditions on boundaries (Dydx<>0)
C      IBC=5 use periodic conditions
C
C    Note that the neuman BC must be set if either IBC = 2, 3 or 4
c       assuming the total interval is unity
C-------------------------------------------------------------
        real(kind=real64) A(N1,N2),B(N1,N2)
        real(kind=real64) uppr(NMAX),diag(NMAX),lowr(NMAX)
        DATA DIAG/NMAX*1.0/
        data revision/1.0/
C-------------------------------------------------------------
        if (ibc.lt.0) then ! return revision level
           a(1,1)=revision
           return
        endif
        if (ibc.gt.5) then
           print *,' dbyd1: ibc out of range, ibc=',ibc
           stop
        endif
C-------------------------------------------------------------
c
        if (ibc.eq.5) then
           d1i=FLOAT(N1)
        else
           d1i=FLOAT(N1-1)
        endif
C-------------------------------------------------------------

        if (n1.lt.5) then
           if (ibc.eq.0) then
              do j=1,n2
             do i=1,n1
                a(i,j)=0.0
             enddo
              enddo
           else if(ibc.eq.1) then
c       linear fit
              do j=1,n2
             do i=1,n1
                a(i,j)=b(n1,j)-b(1,j)
             enddo
              enddo
           else if (ibc.eq.2) then
              do j=1,n2
             do i=1,n1
                a(i,j)=a(1,j)
             enddo
              enddo
           else if (ibc.eq.3) then
              do j=1,n2
             do i=1,n1
                a(i,j)=a(n1,j)
             enddo
              enddo
           else if (ibc.eq.4) then
              do j=1,n2
             do i=1,n1
                a(i,j)=(float(i-1)*a(1,j)+float(n1-i)*a(n1,j))/d1i
             enddo
              enddo
           endif
           return
        endif
C-------------------------------------------------------------
        C1 = d1i*7.0/9.0
        C2 = d1i/36.0
C-------------------------------------------------------------
C       Generate right hand side and store in A
C       First do interior
C-------------------------------------------------------------
        DO I=1,N2
           DO J=3,N1-2
C       
            A(J,I) = C1*(B(J+1,I)-B(J-1,I))+C2*(B(J+2,I)-B(J-2,I))
C       
           ENDDO
        ENDDO
C--------------------------------------------------------------
C       Now the top and bottom, decreasing to fifth order
C______________________________________________________________
        IF(IBC .LT. 5) THEN
           C1 = -43./96.*d1i
           C2 = -5./6.*d1i
           C3 = 9./8.*d1i
           C4 = 1./6.*d1i
           C5 = -1./96.*d1i
c       
         DO 20 I=1,N2
C       
            A(2,I)    = C1*B(1,I) + C2*B(2,I) + C3*B(3,I) 
     $                  + C4*B(4,I) + C5*B(5,I) 
 20     CONTINUE
           do 25 i=1,n2
              A(N1-1,I) =-(C1*B(N1,I) + C2*B(N1-1,I) + 
     $              C3*B(N1-2,I) + C4*B(N1-3,I) + C5*B(N1-4,I))
 25     continue
C       
       
c       fifth order coeff.
           C1 = -10/3.*d1i
           C2 = -3.0*d1i
           C3 = 6.*d1i
           C4 = 1./3.*d1i
           IF (IBC.EQ.0) THEN
              DO 35 I=1,N2
             A(1,I)=0.0
 35         CONTINUE
c       if ibc = 2,4 then a(*,1) must be set by the calling routine
           ELSE if ((ibc.eq.1).or.(ibc.eq.3)) then
              DO 30 I=1,N2
             A(1,I)=C1*B(1,I)+C2*B(2,I)+C3*B(3,I)+C4*B(4,I)
 30         CONTINUE
           ENDIF
           IF (IBC.EQ.0) THEN
              DO 37 I=1,N2
             A(N1,I)=0.0
 37       CONTINUE
c       if ibc = 3,4 then a(*,n2) must be set by the calling routine
           ELSE if ((ibc.eq.1).or.(ibc.eq.2))  then
              DO 38 I=1,N2
             A(N1,I)=-(C1*B(N1,I)+C2*B(N1-1,I)+
     $             C3*B(N1-2,I)+C4*B(N1-3,I))
 38         CONTINUE
        ENDIF
c       
C       
C--------------------------------------------------------------
C       Now we set up the matrix
C--------------------------------------------------------------
c       here is the sixth order interior value
           alpha=1./3.
           DO 50 J=3,N1-2
            UPPR(J) = alpha
 50     CONTINUE
c       here are the pentadiagonal and fifth order values for the boundary.
           alpha2=3./4.
           gamma2=1./8.
           alpha1=6.0
           beta1=3.
c       precondition the matrix to make it tridiagonal
           const=1./(alpha2-beta1*gamma2)
           up1=(alpha1*alpha2-beta1)*const
           if (mod(ibc,2).ne.0) then
              do 80 i=1,n2
             a(1,i)=(alpha2*a(1,i)-beta1*a(2,i))*const
 80         continue
           endif
           if ((ibc.eq.1).or.(ibc.eq.2)) then   
              do 85 i=1,n2
             a(n1,i)=(alpha2*a(n1,i)-beta1*a(n1-1,i))*const
 85         continue
        endif
c       
           IF (mod(IBC,2).EQ.0) THEN
              UPPR(1) = 0.0
           ELSE
c       fifth order bc.
              uppr(1)=up1
           ENDIF
c       
        uppr(2)=alpha2
c       
        uppr(n1-1)=gamma2
c       
           DO 60 I=1,N1
              LOWR(I)=UPPR(N1-I+1)
 60     CONTINUE
c       
           if (ibc.ge.0) then
              IF ((IBC.NE.1).and.(IBC.NE.2)) THEN
             lowr(n1) = 0.0
              ELSE
             lowr(n1)=up1
              ENDIF
           endif
       
C-------------------------------------------------------------
C       And solve it, storing the results back into A
C-------------------------------------------------------------
C$DOACROSS SHARE(LOWR,UPPR,DIAG,A,N1),LOCAL(I)
!$OMP PARALLEL DO SHARED(LOWR, UPPR, DIAG, A, N1) PRIVATE(I)
           DO 100 I=1,N2
              CALL TRIDAG(LOWR,DIAG,UPPR,A(1,I),A(1,I),N1,0)
 100        CONTINUE
        ELSE
c       periodic conditions
           alpha=1./3.
           C1 = d1i*7.0/9.0
           C2 = d1i/36.0
C-------------------------------------------------------------
C       Generate right hand side and store in A
C       On edges
C-------------------------------------------------------------
        DO I=1,N2
C       
              A(2,I) = C1*(B(3,I)-B(1,I))+C2*(B(4,I)-B(N1,I))
              A(1,I) = C1*(B(2,I)-B(N1,I))+C2*(B(3,I)-B(N1-1,I))
c       
              A(N1-1,I) = C1*(B(N1,I)-B(N1-2,I))+C2*(B(1,I)-B(N1-3,I))
              A(N1,I) = C1*(B(1,I)-B(N1-1,I))+C2*(B(2,I)-B(N1-2,I))
C       
        ENDDO
C-------------------------------------------------------------
C       Now generate matrix elements
C-------------------------------------------------------------
           do j=1,n1
              uppr(j)=alpha
              lowr(j)=alpha
           enddo
C-------------------------------------------------------------
C       and do the job
C-------------------------------------------------------------
C$DOACROSS SHARE(LOWR,UPPR,DIAG,A,N1),LOCAL(I)
!$OMP PARALLEL DO SHARED(LOWR, UPPR, DIAG, A, N1) PRIVATE(I)

           DO  I=1,N2
              CALL TRIDAG(LOWR,DIAG,UPPR,A(1,I),A(1,I),N1,1)
           ENDDO
!$OMP END PARALLEL DO
        ENDIF
        RETURN
        END
