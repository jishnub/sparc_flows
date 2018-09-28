# sparc_flows
Sesimic Propagation through Active Regions and Convection

# Getting started

## Prerequisites

You would need a recent version of gcc (6+ is recommended), LAPACK, CFITSIO, FFTW v3, MPI. 
You might need curl too for some reason with gcc.
It is strongly recommended that you get the Anaconda distribution for python, 
but any distribution of python3 with the appropriate collection of packages would do.

## Installation

Clone the project using 
```
git clone https://github.com/jishnub/sparc_flows.git <optional: your directory name>
```

Edit the makefile to set the path of your libraries (if they are non-standard). 
Edit the file `params.i` and set the parameter `directory` to the location where you want to save your simulation output 
(make sure that there's a trailing slash in the path, and the direcory exists).
Change the other parameters as appropriate. Compile the code using 
```
make
f2py -c mtridag.f dbyd2.f -m dbyd2
f2py -c modefilters.f90 -m modefilters
```
The make command would return warnings about array sizes, these can be ignored.

## Creating directories

You would need the following structure.
```bash
<directory>
├── master.pixels
└── params.0
...
```

The file `master.pixels` has the horizontal coordinate of the sources, each source on a new line. 
The files `params.n` for various values of n contain filtering parameters for different radial orders 
(n=0 being the f-mode, n=1 the p1-mode etc)

```
minimum distance of receiver from source (Mm)
maximum distance of receiver from source (Mm)
temporal window size used to compute travel times (minutes)
```

These parameters can be estimated by using the jupyter notebook `ridgewise_distance_for_tt_calculation.ipynb`, 
but this has to be run after the simulation has been run once to compute the data.

## Generating the true and starting models to use for the simulation
This can be done either in an interactive notebook using `generate_true_and_starting_models.ipynb` or non-interactively 
using `generate_true_and_starting_models.py`. These would create the appropriate fits and other files at the correct locations.

## Setting up locations of executables
Edit the file `data_forward.sh` to set the amount of resources that you would require, 
as well as the location of python and MPI that you would want to use (in case the locations are non-standard). 
Note that we require python v3.x for the code.

Edit the file `data_forward.py` to set the location of libraries that would be searched for by the executable. 
An easy way to check these paths would be by running 
```
ldd sparc
```
and noting down the non-standard locations. A very rudimentary shell is launched by `pbsdsh` that does not have access to the 
usual environment variables such as `LD_LIBRARY_PATH`, so we need to update these variables from within the python script.

Similar modifications need to be carried out for `full.sh/full.py` and `linesearch.sh/linesearch.py`. 
Right now the scripts contain the default settings as on the TIFR cluster, so people running it on the cluster should not need 
to modify these paths at all.

# Running the code
The pbs branch assumes that you are running your code on a cluster that runs on PBS Torque 
(although it might work on other PBS environments since it simply uses `pbsdsh` to spawn multiple jobs across processors.

## Computing data
Submit your code using 
```
qsub data_forward.sh -N <job name>
```
To track the status of the PBS job, one can run 
```
qstat -u $(whoami)
```

To track the progress of the individual simulations, one can run the following:

```
(directory=$(python -c 'import read_params; print(read_params.get_directory())');nsrc=$(cat $directory/master.pixels|wc -l); tail $directory/forward_src0{1..$nsrc}_ls00/out_data_forward) 
```

This will write out the last few lines of output from each simulation that have been logged.

At the end of the simulation, verify that directory/data contains the appropriate fits files. One way to check this is by running
```
(directory=$(python -c 'import read_params; print(read_params.get_directory())'); ls -1 $directory/data)
```
and verifying that the output matches with the output of
```
(directory=$(python -c 'import read_params; print(read_params.get_directory())'); nsrc=$(cat $directory/master.pixels|wc -l); printf '%02d.fits\n' {1..$nsrc} )
```
If this works, the simulation has finished correctly. If it doesn't, check the `error-data_forward` file for error messages.



## Compute filter parameters

Once the data has been computed, run the notebook `ridgewise_distance_for_tt_calculation.ipynb` and get the appropriate 
parameters for the ridge filter. 
This does not write out the parameters immediately, one has to create the appropriate params.n files by hand.
Phase speed filters are also available, modify the `adjont_source_filt` subroutine in `driver.f90` 
appropriately to use them. Phase-speed filter parameters are not implemented in the notebook.


## Compute sensitivity kernels

To compute the data misfits between the true and starting models, and the sensitivity kernels about the starting model, run the 
following:
```
qsub full.sh -N <job name>
```

The output of the forward calculation can be tracked using
```
(directory=$(python -c 'import read_params; print(read_params.get_directory())');nsrc=$(cat $directory/master.pixels|wc -l); tail $directory/forward_src0{1..$nsrc}_ls00/out_forward) 
```
The output of the adjont wavefield calculation can be tracked using 
```
(directory=$(python -c 'import read_params; print(read_params.get_directory())');nsrc=$(cat $directory/master.pixels|wc -l); tail $directory/adjoint_src0{1..$nsrc}/out_adjont) 
```
The output of the kernel computation can be tracked using 
```
(directory=$(python -c 'import read_params; print(read_params.get_directory())');nsrc=$(cat $directory/master.pixels|wc -l); tail $directory/kernel/out_kernel_src0{1..$nsrc}) 
```
Once the computation is over, you should have the file `misfit_00` in <directory>/update. This contains the total misfit summed over radial orders 
for each source. An easy way to check the total misfit summed over all sources is by running 
```
python misfit.py <optional:iteration number, chooses the latest iteration by default>
```

## Compute gradients in parameter space
To compute gradients in parameter space, run 
```
python grad.py algo=bfgs <list of step sizes to use for linesearch, arbitrary but better to use small numbers ~1e-2>
```
Make sure that you enter the correct number of steps, eg for 6 linesearches the command chould be 
```
python grad.py algo=cg 0.01 0.02 0.03 0.04 0.05 0.06
```
The first step would feature a steepest descent, followed by conjugate gradient or some such algorithm at subsequent steps.

## Compute linesearches

Run the following command to compute linesearches:
```
qsub linesearch.sh -N <job name>
```
The status of the individual linesearches can be checked using

```
(directory=$(python -c 'import read_params; print(read_params.get_directory())');nsrc=$(cat $directory/master.pixels|wc -l); numls=$(find $directory -name "forward_src01_ls0[1-9]"|wc -l); tail $directory/forward_src0{1..$nsrc}_ls0{1..$numls}/out_forward)
```

This will however output a wall of text. If you just want to check the progress of one linesearch for each source, you can run 
```
(directory=$(python -c 'import read_params; print(read_params.get_directory())');nsrc=$(cat $directory/master.pixels|wc -l); tail $directory/forward_src0{1..$nsrc}_ls01/out_forward)
```
Once the linesearches finish executing, you can check the misfits for each run using 
```
python lsmisfit.py
```

## Computing the ideal stepsize
If the linesearch misfits contains a minimum, you can skip to the next step. Otherwise we need to estimate the best step size 
(irrespective of whether we are going downwards or upwards on the misfit curve). This is done by carrying out a quadratic fit 
to the misfits. Run 
```
python compute_ideal_stepsize.py
```
This will tell you the step sizes that you should use to obtain the minimum. It is a good idea to re-run the linesearch using 
these new step sizes to verify that there is indeed a minimum. This is especially necessary if the computed step sizes involve 
severe extrapolation. 

To re-run the linesearch, run
``` 
python remove_last_linesearch.py
python grad.py algo=cg <new step sizes>
qsub linesearch.sh -N <job name>
```

## Copy the new model
Once we have obtained a minimum in the linesearch, we need to update our starting model. To do this run 
```
python copy_updated_model.py
```

At this stage, we have successfully completed one iteration, and can start the second iteration by running 

```
qsub full.sh -N <job name>
```


## Iteration status
Because of the large number of steps involved, it is natural to forget which stage we are at. 
There is a script that can remind us what to run next. To check the next step at any stage, run 
```
python iterstatus.py
```

