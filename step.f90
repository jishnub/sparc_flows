Module Time_Step

! --------------------------------------------------------------------------
! MPI Version of the Spherical Acoustic Sun Simulator.
! Copyright 2006, Shravan Hanasoge
                                                                                                                                                         
! Hansen Experimental Physics Laboratory
! 455 Via Palou way, Stanford
! CA 94305, USA
! Email: shravan@stanford.edu
! --------------------------------------------------------------------------

  use initialize
  use mp_physics
  use mp_physics_2d
  use pml
  use displacement
  Implicit None

  real*8, dimension(5) :: betas, optimals

Contains
!=====================================================================================
  SUBROUTINE INITIALIZE_STEP
  
     implicit none

     optimals(1) =  1.0  
     optimals(2) =  0.5
     optimals(3) =  0.1665579725791151184
     optimals(4) =  0.03950410250586729987
     optimals(5) =  0.007810706393337838236

     betas(5) = optimals(2)   
     betas(4) = optimals(3)/betas(5)
     betas(3) = optimals(4)/(betas(5)*betas(4))
     betas(2) = optimals(5)/(betas(5)*betas(4)*betas(3))
     betas(1) = 0.0

     betas = betas * deltat
  end SUBROUTINE INITIALIZE_STEP
!=====================================================================================
  SUBROUTINE STEP

   implicit none
   integer :: i,j,k,l
  ! All new and improved optimized RK 5, 4th order accurate integrator. 50% increase in time step.

   do step_rk = 1,5


     if (step_rk == 1) then
       temp_step = a

       if (USE_PML) psipml = pmlvars

       if (HORIZONTAL_PMLS) then 
	psipmlx = pmlvarsx
	if (PROC_HAS_PML) psipmly = pmlvarsy
       endif
    else  
       temp_step = a + betas(step_rk)*scr
       if (USE_PML) psipml = pmlvars + betas(step_rk) * scrpml

       if (HORIZONTAL_PMLS) then 
	psipmlx = pmlvarsx + betas(step_rk) * scrpmlx
	if (PROC_HAS_PML) psipmly = pmlvarsy + betas(step_rk) * scrpmly
       endif


     endif
 

      if (option==1) then
       call MP_QUIET_SPONGE_3D ()
      elseif (option==2) then
       call MP_QUIET_PML_3D ()
      elseif (option==3) then
       call MP_QUIET_SPONGE_3D_DISPL ()
      elseif (option==4) then
       call MP_QUIET_PML_3D_DISPL ()
      elseif (option==5) then
       call MP_QUIET_SPONGE_2D ()
      elseif (option==6) then
       call MP_QUIET_PML_2D ()
      elseif (option==7) then
       call MP_MHD_SPONGE_3D ()
      elseif (option==8) then
       call MP_MHD_PML_3D ()
      elseif (option==9) then
       call MP_MHD_SPONGE_3D_DISPL ()
      elseif (option==10) then
       call MP_MHD_PML_3D_DISPL ()
      elseif (option==11) then
       call MP_MHD_SPONGE_2D ()
      elseif (option==12) then
       call MP_MHD_PML_2D ()
      elseif (option==13) then
       call MP_FLOWS_SPONGE_3D ()
      elseif (option==14) then
       call MP_FLOWS_PML_3D ()
      elseif (option==15) then
       call MP_FLOWS_SPONGE_3D_DISPL ()
      elseif (option==16) then
       call MP_FLOWS_PML_3D_DISPL ()
      elseif (option==17) then
       call MP_QUIET_PML_2D_DISPL ()
      elseif (option==18) then
       call MP_MHD_PML_2D_DISPL ()
      endif

   enddo

   a = a + deltat*scr

   if (USE_PML) pmlvars = pmlvars + deltat * scrpml
   if ((.not. TEST_IN_2D)  .and. (.not. HORIZONTAL_PMLS)) call filter_vars_3d ()
   if (TEST_IN_2D) call filter_vars_2d ()

   if (HORIZONTAL_PMLS) then 
    pmlvarsx = pmlvarsx + deltat * scrpmlx
    if (PROC_HAS_PML .and. (.not. TEST_IN_2D)) pmlvarsy = pmlvarsy + deltat * scrpmly
    if (mod(time+1,2) ==0 .and. (.not. TEST_IN_2D)) call filter_vars_3d_horiz_pml
   endif



  End SUBROUTINE step
!=====================================================================================
End Module Time_Step
