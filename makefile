## Makefile
## -----------------------------------------------------
## MPI Version of the Spherical Acoustic Sun Simulator.
## Copyright 2006, Shravan Hanasoge
##
## Hansen Experimental Physics Laboratory
## 455 Via Palou way, Stanford
## CA 94305, USA
## Email: shravan@stanford.edu
## -----------------------------------------------------
##

OBJS1=   driver.o        initialize.o    physics.o       dbyd2.o\
        mtridag.o       step.o 	all_modules.o   dbyd1.o tridag.o\
	physics2d.o	derivatives.o	pml.o	displacement.o\
	damping.o	kernels.o	bspline90_22.o integrals.o spline.o \
	splevl.o splint.o
OBJS2= tt_test.o

FC= mpif90
FC77= mpif77

FFLAGS= -DDOUBLE_PRECISION -g ##-p ##-check all ##-fpe0 -traceback -debug #-check bounds
#LIBS1 = -lfftw3 -lcfitsio
LIBS1 = -L/home/apps/fftw-3.2/lib -lfftw3 -lcfitsio -lgfortran -L/home/apps/lapack-3.5 -llapack -lcurl #-lquadmath

COMMAND1= sparc
COMMAND2= tt_test

$(COMMAND1): $(OBJS1) 
	@$(FC) -I $(INCLUDE) $(FFLAGS) -o $(COMMAND1) $(OBJS1) $(LIBS1) 

$(COMMAND2): $(OBJS2)
	@$(FC) $(FFLAGS) -o $(COMMAND2) $(OBJS2) $(LIBS1) 

%.o : %.f
	@$(FC77) $(FFLAGS) -c $< 

%.o : %.f90
	@$(FC) $(FFLAGS) -c $< 

clean:
	@find . -maxdepth 1 -name "*.o" -delete
	@find . -maxdepth 1 -name "*.mod" -delete
	@find . -maxdepth 1 -name "$(COMMAND1)" -delete



initialize.o:	params.i
driver.o:       initialize.o    all_modules.o   physics.o       step.o	\
				kernels.o	physics2d.o	bspline90_22.o derivatives.o \
				integrals.o
physics.o:      initialize.o    all_modules.o	derivatives.o	damping.o
dbyd2.o:        mtridag.o
step.o: 	physics.o	physics2d.o	pml.o	displacement.o	initialize.o
all_modules.o:  initialize.o
pml.o:		initialize.o	derivatives.o	all_modules.o	damping.o
dbyd1.o:        tridag.o
physics2d.o:	derivatives.o	damping.o	initialize.o
derivatives.o:	initialize.o	dbyd1.o	dbyd2.o
displacement.o:	initialize.o	derivatives.o	physics.o	damping.o
damping.o:	initialize.o
process.i:	params.i
kernels.o:	initialize.o	all_modules.o
integrals.o: splint.o
splint.o: splevl.o spline.o
splevl.o: spline.o
