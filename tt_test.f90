        program tt_test
            implicit none
            integer, parameter :: nt=512
            real*8 u(nt),u0(nt)
            integer i,lef,rig
            real*8 ti,T,w,tau,dt
            real*8,parameter :: pi = acos(dble(-1.0))

            T = 240*60
            dt = T/nt
            w = 2D-3

            do i=1,nt
                ti = i*dt
                u(i) = exp(-(ti-75*60)**2/(2*600**2))
                u0(i) = exp(-(ti-76*60)**2/(2*600**2))
            end do

            lef = int(30*60/T*nt)
            rig = int(100*60/T*nt)

            call compute_tt_gizonbirch(u0,u,tau,dt,nt,lef,rig)
            print*,"Travel time shift",tau
        end program

        SUBROUTINE COMPUTE_TT_GIZONBIRCH(u0,u,tau,dt,nt, lef, rig)
            implicit none
            include 'fftw3.f'
            integer, parameter :: degree=3
            integer, intent(in) :: nt, lef, rig
            real*8, intent(in) :: u0(nt), u(nt)
            real*8 dt
            real*8, intent(out) :: tau
            real*8 window(nt),functemp(nt)
            integer k,i,num_real_roots
            integer*8 plan
            real*8 temp,polycoeffs(0:degree)
            real*8 :: roots(degree)
            complex*16 wu0(nt/2+1,1:degree+1)
            real*8 u0dot(nt,1:degree+1)
            complex*16, parameter :: eye = (0.0,1.0)
            real*8, parameter :: pi=acos(dble(-1.0))

            polycoeffs = 0

            window = 0.0
            window(lef:rig) = 1.0

            do i=1,degree+1
                call dfftw_plan_dft_r2c_1d(plan,nt,u0,wu0(:,i),&
                        FFTW_ESTIMATE)
                call dfftw_execute_dft_r2c(plan, u0, wu0(:,i))
                call dfftw_destroy_plan(plan)

                do k=1,nt/2
                    wu0(k,i)=wu0(k,i)*(eye*2*pi*(k-1.0)/(nt*dt))**i
                end do

                call dfftw_plan_dft_c2r_1d(plan,nt,wu0(:,i),u0dot(:,i),&
                        FFTW_ESTIMATE)
                call dfftw_execute_dft_c2r(plan, wu0(:,i),u0dot(:,i))
                call dfftw_destroy_plan(plan)
            enddo

            u0dot=u0dot/nt

            ! coefficients of dchi/dt

            do i=0,degree
                functemp = 0
                if (i .eq. 0) then
                    functemp = 2*u0dot(:,1)*(u-u0)
                elseif (i .eq. 1) then
                    functemp = 2*(u0dot(:,1)**2+(-u+u0)*u0dot(:,2))
                elseif (i .eq. 2) then
                    functemp = -3*u0dot(:,1)*u0dot(:,2)+&
                        (u-u0)*u0dot(:,3)
                elseif (i .eq. 3) then
                    functemp = (3*u0dot(:,2)**2+4*u0dot(:,1)*u0dot(:,3)&
                        +(-u+u0)*u0dot(:,4))/3.0
                endif

                call integrate_time(window*functemp,&
                        polycoeffs(degree-i),dt,nt)
            end do

            polycoeffs = polycoeffs/polycoeffs(0)

            tau=-polycoeffs(degree)/polycoeffs(degree-1)

            call polyroots(polycoeffs,degree,num_real_roots,roots)

            tau=roots(minloc(abs(roots(1:num_real_roots)-tau),dim=1))

        END SUBROUTINE COMPUTE_TT_GIZONBIRCH

        SUBROUTINE INTEGRATE_TIME(f,int_f,dt,nt)

            implicit none
            integer, intent(in) :: nt
            real*8, intent(in) :: dt,f(nt)
            real*8, intent(out) ::  int_f

            int_f = 0
            call simpson_regular(f,dt,nt,int_f)

        END SUBROUTINE INTEGRATE_TIME

        SUBROUTINE SIMPSON_REGULAR(f,dvar,nt,int_f)

            implicit none
            integer, intent(in) :: nt
            real*8, intent(in) :: dvar,f(nt)
            real*8, intent(out) :: int_f
            integer k

            int_f = 0.0

            do k=1,nt

                if ((k == 1) .or. (k == size(f))) then
                    int_f = int_f + f(k)
                    cycle
                end if

                if (mod(k,2) == 0) then
                    int_f = int_f + 4*f(k)
                else
                    int_f = int_f + 2*f(k)
                endif

            end do

            int_f=int_f*dvar/3.

        END SUBROUTINE SIMPSON_REGULAR


        subroutine POLYROOTS(polycoeffs,N,num_real_roots,realroots)
            implicit none
            integer, intent(in) :: N
            real*8, intent(in) :: polycoeffs(0:N)
            real*8 companion(N,N),wi(N),wr(N)
            real*8, intent(out) :: realroots(N)
            integer,intent(out) :: num_real_roots
            integer i,info,lwork,lwmax,temp
            parameter(lwmax=100)
            real*8 work(lwmax),vl(N,N),vr(N,N)

          companion=0
          do i=1,N
            companion(1,i) = -polycoeffs(i)
          end do

          do i=2,N
            companion(i,i-1) = 1
          enddo

          LWORK = -1
          CALL DGEEV( 'N', 'N', N, companion, N, WR, WI, VL, N,&
                       VR, N, WORK, LWORK, INFO )
          LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )

          call DGEEV('N','N',N,companion,N,wr,wi,vl,N,&
                        vr,N,work,lwork,info)

            num_real_roots=0
            do i=1,N
                if (abs(wi(i))<1e-10) num_real_roots=num_real_roots+1
            end do

           temp = 1

           do i=1,N
               if (abs(wi(i))<1e-10) then
                   realroots(temp) = wr(i)
                   temp = temp+1
               end if
           end do

        end subroutine POLYROOTS
