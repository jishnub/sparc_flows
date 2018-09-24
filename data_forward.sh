#!/bin/bash
#PBS -l nodes=1:ppn=24
#PBS -o  output-data_forward
#PBS -e  error-data_forward
#PBS -l walltime=12:00:00
echo $PBS_JOBID
export TERM=xterm
cd $PBS_O_WORKDIR
num_src=$(cat master.pixels|wc -l)
num_ls=3
~/anaconda3/bin/python -c "import setup; setup.create_directories($num_src,$num_ls)"
echo "Starting at "`date`
find . -name "linesearch" -delete
touch compute_data
/usr/local/bin/pbsdsh $HOME/anaconda3/bin/python $PBS_O_WORKDIR/data_forward.py
find . -name "compute_data" -delete
echo "Finished at "`date`
