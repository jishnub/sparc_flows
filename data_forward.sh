#!/bin/bash
#PBS -N  data_f_to_p5_3hr
#PBS -l nodes=1:ppn=24
#PBS -o  output-data_forward
#PBS -e  error-data_forward
#PBS -l walltime=12:00:00
echo $PBS_JOBID
export TERM=xterm
cd $PBS_O_WORKDIR
python $PBS_O_WORKDIR/setup.py
export MPI_TYPE_MAX=1280280
echo "Starting at "`date`
find . -name "linesearch" -delete
touch compute_data
/usr/local/bin/pbsdsh python $PBS_O_WORKDIR/data_forward.py
find . -name "compute_data" -delete
echo "Finished at "`date`
