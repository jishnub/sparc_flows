#!/bin/bash
#PBS -N  data_forward
#PBS -l nodes=1:ppn=24
#PBS -o  output-data_forward
#PBS -e  error-data_forward
#PBS -l walltime=12:00:00
echo $PBS_JOBID
export TERM=xterm
cd $PBS_O_WORKDIR
export MPI_TYPE_MAX=1280280
echo "Starting at "`date`
touch $PBS_O_WORKDIR/compute_data
/usr/local/bin/pbsdsh python $PBS_O_WORKDIR/data_forward.py
find $PBS_O_WORKDIR -name "compute_data"  -exec rm -f {} \;
echo "Finished at "`date`
