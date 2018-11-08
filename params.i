! --------------------------------------------------------------------------
! 'params.i' : SIMULATION PARAMETERS
! --------------------------------------------------------------------------
! MPI Version of the Cartesian Solar Wave Simulator.
! Copyright 2008, Shravan Hanasoge
!
! W. W. Hansen Experimental Physics Laboratory
! Stanford University, Stanford
! CA 94305, USA
! Email: shravan@stanford.edu
!
! Also see: http://soi.stanford.edu/
!
! --------------------------------------------------------------------------
!

integer, parameter :: real64 = kind(1.0d0)

! PARAMETERS OF THE SIMULATION
! nx, ny, nz are the numbers of grid points in the x,y,z directions resply.
! nz = 300 generally works well.
! xlength and ylength are the horizontal sides of the box, expressed in cgs
! A timestep of 2 seconds is generally pretty solid. Too large and the sim
! will explode. Too small and the expense is large.
integer nx, ny, nz
parameter (nx = 256, ny = 1, nz = 300)
real*8 xlength, ylength, timestep
parameter (xlength = 300.0 * 10**(8), ylength = xlength, timestep = 2.0)


! DIRECTORY INFORMATION
!
! ENTER THE LOCATION OF THE BACKGROUND MODEL (i.e. the quiet Sun model)
character (LEN = *), PARAMETER :: file_data = 'solar_deep'
!'solar_model'

! IF THE MODEL REQUIRES STABILIZING
logical :: STABILIZE_MODEL = .false.

! DIRECTORY FOR OUTPUT/ SAVED RESTART STATE (ASSUMING THEY ARE THE SAME)
character (LEN = *), PARAMETER :: directory = '/scratch/jishnu/sparc/soundspeed/'

! THE FORCING FUNCTION
character (LEN = *), PARAMETER :: forcingfunc = '/nobackup/shanasog/classic4/ccsource.fits'

! CADENCE OF THE FORCING FUNCTION (i.e. SAMPLING RATE) in SECONDS
real*8 cadforcing
parameter( cadforcing = 30.0)

! NOTE THIS fUNCTION STILL NEEDS TO BE IMPLEMENTED
! I HAVE THE BASIC MACHINERY IN PLACE BUT NEED DAMPING INFO
! DAMPING FUNCTION - IN ORDER TO DAMP WAVES - IN REAL SPACE (x,y,z)
logical, parameter :: DAMP_WAVES = .falsE.


! TYPE OF SIMULATION (DEFAULT IS QUIET OR SOURCE/SOUND-SPEED ANOMALIES):: FLOWS OR MAGNETIC FIELDS
logical, parameter :: sound_speed_perturbation = .True.
logical, parameter :: magnetic = .falsE.
logical, parameter :: FLOWS = .false.
logical :: TEST_IN_2D != .false.


! IF KERNELS
real*8 nupeak, nuwidth
integer st_cc, st_adj, stepskern, nz_kern, nt_kern, st_z, fi_z, totkern, forcing_length
parameter(st_z = 1, fi_z = nz, nz_kern = (fi_z - st_z+1),nupeak=0.0017, nuwidth = 0.0055)
logical BACKGROUND_FLOWS_EXIST, COMPUTE_FORWARD !, COMPUTE_CC_SLICE
logical COMPUTE_ADJOINT, KERNEL_MODE, COMPUTE_POWER_SPEC, CONSTRUCT_KERNELS
parameter(KERNEL_MODE =.true.)
parameter(background_flows_exist = .false.)
logical sound_speed_kernels, flow_kernels, magnetic_kernels, density_kernels


! TYPE OF EQUATION
! IF DISPLACEMENT, SET DISPL = TRUE
logical, parameter :: DISPL = .truE.

! VERTICAL BOUNDARY
! IF USE_PML IS FALSE, A SPONGE WILL BE IMPLEMENTED
logical USE_PML
parameter(USE_PML = .true.)

! HORIZONTAL BOUNDARIES
! IF PML
logical HORIZONTAL_PMLS
parameter(HORIZONTAL_PMLS = .true.)

integer npmlhor
parameter (npmlhor = 9)


! ENTER THE LOCATION OF THE 3D pressure, density, and magnetic background files
! DEFAULT LOCATION IS THE SAME AS OUTPUT DIRECTORY
character (LEN = *), PARAMETER :: dirbackmag = '/scratch/jishnu/flows/start/'
character (LEN = *), PARAMETER :: dirbacktrue = '/scratch/jishnu/flows/true/'

! --------

! STARTING TIME STEP; if starting off, time0 =0
! Otherwise, enter the index of the latest "full" file
integer, parameter :: time0 = 0

! HORIZONTAL BOUNDARY CONDITIONS
! Note since there are two ways to evaluate periodic horizontal
! boundary conditions, these two flags have to be set. If not FFTs,
! compact finite differences will be used to compute the periodic
! derivatives. Of course if periodic is set to false, then normal
! absorbing boundary conditions with derivatives computed using the
! compact finite differences are implemented.
logical, parameter :: PERIODIC = .true.
logical, parameter :: USE_FFT = .truE.

! VERTICAL DERIVATIVE: USING EITHER SIXTH-ORDER ACCURATE COMPACT
! FINITE DIFFERENCES OR NINTH-ORDER ACCURATE OPTIMIZED EXPLICIT
! DIFFERENCES - SET COMPACT_FINITE_DIFF = .FALSE. IF YOU WANT TO
! USE THE LATTER
! RECOMMENDED: USE COMPACT FINITE DIFF.
logical, parameter :: compact_finite_diff = .truE.

! IF USE_PML = .TRUE., THEN SET THE FOLLOWING (DEFAULT VALUES WORK OKAY)
integer npmltop, npmlbot,nzpml
parameter(npmlbot = 10, npmltop = 10, nzpml = (npmlbot + npmltop))


! THE TOTAL WALL-TIME-LENGTH OF THE SIMULATION (IN HOURS)
real *8 wall_time
parameter(wall_time = 96.0)

! SOLAR TIME BEING SIMULATED (in hours)
real*8 solartime
parameter (solartime = 4.0)

! OBSERVATION HEIGHT RELATIVE TO PHOTOSPHERE (in cgs units)
! Generally 200 km above photosphere works well
! The sign indicates whether you wish to extract data above (+)
! or below (-) the photosphere
real*8 obsheight
parameter (obsheight = 20000000.0)

! CADENCE OF SLICE OUTPUT (in seconds; 30 s recommended)
real *8 outputcad
parameter (outputcad = 30.0)

! NUMBER OF ZEROES IN THE OUTPUT FILENAME
integer TIMESTAMP_SIZE
parameter (timestamp_size = 6)

! WAVE EXCIATION DEPTH RELATIVE TO PHOTOSPHERE
! Generally 150 km below the photosphere is a good place to excite the waves
! excitdep is in cgs units
real*8 excitdep
parameter( excitdep = - 15000000.0)

! CONTINUITY
logical enf_cont, psi_cont, vx_cont, vz_cont
parameter ( enf_cont = .TRUE. )
parameter ( psi_cont = .TRUE. )
parameter ( vx_cont = .FALSE. )
parameter ( vz_cont = .FALSE. )
