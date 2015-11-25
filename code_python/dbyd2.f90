!**********************************************************
!
!       Copyright 1994, Neal Hurlburt,
!   Lockheed Palo Alto Research Laboratory
!       Revised 3/1998 for periodic conditions
!
!************************************************************
        SUBROUTINE DBYD2(A,B,N1,N2,IBC)
        IMPLICIT NONE
        INTEGER N1,N2,IBC,I,J,NMAX
        REAL*8 D2I,C1,C2,C3,C4,C5,alpha,alpha1,alpha2,beta1,gamma2
        REAL*8 const,revision,up1
        PARAMETER(NMAX=700)
!------------------------------------------------------------
!  This routine evaluates the derivative of array B using
!  sixth-order compact differences and
!  stores the result in array A. The derivative is evaluated
!  along the first argument of B.
!
!      IBC<0 returns revision level in a(1,1)
!      IBC=0 uses neumann conditions on boundaries (Dydx=0)
!      IBC=1 calculates value on boundaries
!      IBC=2 neumann on low end, Dirichelet on high
!      IBC=3 vica versa
!      IBC=4 uses neumann conditions on boundaries (Dydx<>0)
!      IBC=5 use periodic conditions
!
!    Note that the neuman BC must be set if either IBC = 2, 3 or 4
!       assuming the total interval is unity
!-------------------------------------------------------------
        REAL*8 A(N1,N2),B(N1,N2)
        REAL*8 uppr(NMAX),diag(NMAX),lowr(NMAX)
        data diag/NMAX*1.0/
        data revision/1.0/
!-------------------------------------------------------------
        if (ibc.lt.0) then ! return revision level
           a(1,1)=revision
           return
        endif
        if (ibc.gt.5) then
           print *,' dbyd2: ibc out of range, ibc=',ibc
           stop
        endif
!-------------------------------------------------------------
!
        if (ibc.eq.5) then
           d2i=FLOAT(N2)
        else
           d2i=FLOAT(N2-1)
        endif
!-------------------------------------------------------------

        if (n2.lt.5) then
           if (ibc.eq.0) then
              do j=1,n2
             do i=1,n1
                a(i,j)=0.0
             enddo
              enddo
           else if(ibc.eq.1) then
!       linear fit
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
!-------------------------------------------------------------
!       Generate right hand side and store in A
!       First do interior
!-------------------------------------------------------------
        DO 10 I=3,N2-2
        DO 10 J=1,N1
!
        A(J,I) = C1*(B(J,I+1)-B(J,I-1))+C2*(B(J,I+2)-B(J,I-2))
!
10      CONTINUE
!--------------------------------------------------------------
!       Now the top and bottom, decreasing to fifth order
!______________________________________________________________
        IF(IBC .LT. 5) THEN
           C1 = -43./96.*d2i
           C2 = -5./6.*d2i
           C3 = 9./8.*d2i
           C4 = 1./6.*d2i
           C5 = -1./96.*d2i
!
           DO 20 I=1,N1
!
              A(I,2)    = C1*B(I,1) + C2*B(I,2) + C3*B(I,3) &
              + C4*B(I,4) + C5*B(I,5)
20      CONTINUE
           do 25 i=1,n1
              A(I,N2-1) =-(C1*B(I,N2) + C2*B(I,N2-1) + &
                  C3*B(I,N2-2) + C4*B(I,N2-3) + C5*B(I,N2-4))
25      continue
!

!       fifth order coeff.
           C1 = -10/3.*d2i
           C2 = -3.0*d2i
           C3 = 6.*d2i
           C4 = 1./3.*d2i
           IF (IBC.EQ.0) THEN
              DO 35 I=1,N1
             A(I,1)=0.0
35      CONTINUE
!       if ibc = 2,4 then a(*,1) must be set by the calling routine
           ELSE if ((ibc.eq.1).or.(ibc.eq.3)) then
              DO 30 I=1,N1
             A(i,1)=C1*B(i,1)+C2*B(i,2)+C3*B(i,3)+C4*B(i,4)
30      CONTINUE
           ENDIF
           IF (IBC.EQ.0) THEN
              DO 37 I=1,N1
             A(i,N2)=0.0
37      CONTINUE
!       if ibc = 3,4 then a(*,n2) must be set by the calling routine
           ELSE if ((ibc.eq.1).or.(ibc.eq.2))  then
              DO 38 I=1,N1
             A(i,N2)=-(C1*B(i,N2)+C2*B(i,N2-1)+ &
                 C3*B(i,N2-2)+C4*B(i,N2-3))
38      CONTINUE
           ENDIF
!
!--------------------------------------------------------------
!       Now we set up the matrix
!--------------------------------------------------------------
! here is the sixth order interior value
           alpha=1./3.
           DO 50 J=3,N2-2
              UPPR(J) = alpha
50      CONTINUE
!       here are the pentadiagonal and fifth order values for the boundary.
           alpha2=3./4.
           gamma2=1./8.
           alpha1=6.0
           beta1=3.
!       precondition the matrix to make it tridiagonal
           const=1./(alpha2-beta1*gamma2)
           up1=(alpha1*alpha2-beta1)*const
           if (mod(ibc,2).ne.0) then
              do 80 i=1,n1
             a(i,1)=(alpha2*a(i,1)-beta1*a(i,2))*const
80       continue
           endif
           IF ((IBC.EQ.1).OR.(IBC.EQ.2)) THEN
              do 85 i=1,n1
             a(i,n2)=(alpha2*a(i,n2)-beta1*a(i,n2-1))*const
85       continue
           endif
!
           IF (MOD(IBC,2).EQ.0) THEN
              UPPR(1) = 0.0
           ELSE
!       fifth order bc.
              uppr(1)=up1
           ENDIF
!
           uppr(2)=alpha2
!
           uppr(n2-1)=gamma2
!
           DO 60 I=1,N2
              LOWR(I)=UPPR(N2-I+1)
60    CONTINUE
!
           if (ibc.ge.0) then
              IF ((IBC.NE.1).and.(IBC.NE.2)) THEN
             lowr(n2) = 0.0
              ELSE
             lowr(n2)=up1
              ENDIF
           endif

!-------------------------------------------------------------
!       And solve it, storing the results back into A
!-------------------------------------------------------------
        CALL MTRIDAG(LOWR,DIAG,UPPR,A,A,N1,N2,0)

!-------------------------------------------------------------
        ELSE
!       periodic conditions
!-------------------------------------------------------------
!       Generate right hand side and store in A
!       On edges
!-------------------------------------------------------------
           C1 = d2i*7.0/9.0
           C2 = d2i/36.0
           DO J=1,N1
!
              A(J,2) = C1*(B(J,3)-B(J,1))+C2*(B(J,4)-B(J,N2))
              A(J,1) = C1*(B(J,2)-B(J,N2))+C2*(B(J,3)-B(J,N2-1))
!
              A(J,N2-1) = C1*(B(J,N2)-B(J,N2-2))+C2*(B(J,1)-B(J,N2-3))
              A(J,N2) = C1*(B(J,1)-B(J,N2-1))+C2*(B(J,2)-B(J,N2-2))
!
           ENDDO
!-------------------------------------------------------------
!       Now generate matrix elements
!-------------------------------------------------------------
!       here is the sixth order interior value
           alpha=1./3.
           DO J=1,N2
              UPPR(J) = alpha
              lowr(j) = alpha
           enddo
!-------------------------------------------------------------
!       and do the job
!-------------------------------------------------------------
           CALL MTRIDAG(LOWR,DIAG,UPPR,A,A,N1,N2,1)
           
        ENDIF
        RETURN
        END
