c**********************************************************
c
c       Copyright 1994, Neal Hurlburt, 
c   Lockheed Palo Alto Research Laboratory
c       Revised 3/1998 for periodic conditions
c
C************************************************************
        SUBROUTINE DBYD2(A,B,N1,N2,IBC)
        PARAMETER(NMAX=7000)
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
        data diag/NMAX*1.0/
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
           d2i=FLOAT(N2)
        else
           d2i=FLOAT(N2-1)
        endif
C-------------------------------------------------------------

        if (n2.lt.5) then
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
                a(i,j)=b(i,n2)-b(i,1)
             enddo
              enddo
           else if (ibc.eq.2) then
              do j=1,n2
             do i=1,n1
                a(i,j)=a(i,1)
             enddo
              enddo
           else if (ibc.eq.3) then
              do j=1,n2
             do i=1,n1
                a(i,j)=a(i,n2)
             enddo
              enddo
           else if (ibc.eq.4) then
              do j=1,n2
             do i=1,n1
                a(i,j)=(float(j-1)*a(i,1)+float(n2-j)*a(i,n2))/d2i
             enddo
              enddo
           endif
           return
        endif
        C1 = d2i*7.0/9.0
        C2 = d2i/36.0
C-------------------------------------------------------------
C       Generate right hand side and store in A
C       First do interior
C-------------------------------------------------------------  
        DO I=3,N2-2
            DO J=1,N1
C
                A(J,I) = C1*(B(J,I+1)-B(J,I-1))+C2*(B(J,I+2)-B(J,I-2))
            END DO
        END DO
C
C 10      CONTINUE
C--------------------------------------------------------------
C       Now the top and bottom, decreasing to fifth order
C______________________________________________________________
        IF(IBC .LT. 5) THEN
           C1 = -43./96.*d2i
           C2 = -5./6.*d2i
           C3 = 9./8.*d2i
           C4 = 1./6.*d2i
           C5 = -1./96.*d2i
c       
        DO 20 I=1,N1
C       
            A(I,2)    = C1*B(I,1) + C2*B(I,2) + C3*B(I,3) 
     $          + C4*B(I,4) + C5*B(I,5) 
 20     CONTINUE
           do 25 i=1,n1
              A(I,N2-1) =-(C1*B(I,N2) + C2*B(I,N2-1) + 
     $               C3*B(I,N2-2) + C4*B(I,N2-3) + C5*B(I,N2-4))
 25     continue
C       
       
c       fifth order coeff.
           C1 = -10/3.*d2i
           C2 = -3.0*d2i
           C3 = 6.*d2i
           C4 = 1./3.*d2i
           IF (IBC.EQ.0) THEN
              DO 35 I=1,N1
             A(I,1)=0.0
 35         CONTINUE
c       if ibc = 2,4 then a(*,1) must be set by the calling routine
           ELSE if ((ibc.eq.1).or.(ibc.eq.3)) then
              DO 30 I=1,N1
            A(i,1)=C1*B(i,1)+C2*B(i,2)+C3*B(i,3)+C4*B(i,4)
 30         CONTINUE
           ENDIF
           IF (IBC.EQ.0) THEN
              DO 37 I=1,N1
             A(i,N2)=0.0
 37         CONTINUE
c       if ibc = 3,4 then a(*,n2) must be set by the calling routine
           ELSE if ((ibc.eq.1).or.(ibc.eq.2))  then
              DO 38 I=1,N1
             A(i,N2)=-(C1*B(i,N2)+C2*B(i,N2-1)+
     $               C3*B(i,N2-2)+C4*B(i,N2-3))
 38         CONTINUE
        ENDIF
C
C--------------------------------------------------------------
C       Now we set up the matrix
C--------------------------------------------------------------
c here is the sixth order interior value
           alpha=1./3.
           DO 50 J=3,N2-2
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
              do 80 i=1,n1
             a(i,1)=(alpha2*a(i,1)-beta1*a(i,2))*const
 80         continue
           endif
           IF ((IBC.EQ.1).OR.(IBC.EQ.2)) THEN
              do 85 i=1,n1
             a(i,n2)=(alpha2*a(i,n2)-beta1*a(i,n2-1))*const
 85         continue
        endif
c       
           IF (MOD(IBC,2).EQ.0) THEN
              UPPR(1) = 0.0
           ELSE
c       fifth order bc.
              uppr(1)=up1
           ENDIF
c       
        uppr(2)=alpha2
c       
        uppr(n2-1)=gamma2
c       
           DO 60 I=1,N2
              LOWR(I)=UPPR(N2-I+1)
 60     CONTINUE
c       
           if (ibc.ge.0) then
              IF ((IBC.NE.1).and.(IBC.NE.2)) THEN
             lowr(n2) = 0.0
              ELSE
             lowr(n2)=up1
              ENDIF
           endif
       
C-------------------------------------------------------------
C       And solve it, storing the results back into A
C-------------------------------------------------------------
        CALL MTRIDAG(LOWR,DIAG,UPPR,A,A,N1,N2,0)

C-------------------------------------------------------------
        ELSE
c       periodic conditions
C-------------------------------------------------------------
C       Generate right hand side and store in A
C       On edges
C-------------------------------------------------------------
           C1 = d2i*7.0/9.0
           C2 = d2i/36.0
           DO J=1,N1
C       
              A(J,2) = C1*(B(J,3)-B(J,1))+C2*(B(J,4)-B(J,N2))
              A(J,1) = C1*(B(J,2)-B(J,N2))+C2*(B(J,3)-B(J,N2-1))
c       
              A(J,N2-1) = C1*(B(J,N2)-B(J,N2-2))+C2*(B(J,1)-B(J,N2-3))
              A(J,N2) = C1*(B(J,1)-B(J,N2-1))+C2*(B(J,2)-B(J,N2-2))
C       
        ENDDO
C-------------------------------------------------------------
C       Now generate matrix elements
C-------------------------------------------------------------
c       here is the sixth order interior value
           alpha=1./3.
           DO J=1,N2
              UPPR(J) = alpha
              lowr(j) = alpha
           enddo
C-------------------------------------------------------------
C       and do the job
C-------------------------------------------------------------
           CALL MTRIDAG(LOWR,DIAG,UPPR,A,A,N1,N2,1)
        ENDIF
        RETURN
        END
