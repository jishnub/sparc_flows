#!/bin/bash
#PBS -l nodes=1:ppn=8
#PBS -o  output-data_forward
#PBS -e  error-data_forward
#PBS -l walltime=12:00:00
cd $PBS_O_WORKDIR
echo $PBS_JOBID
export TERM=xterm

python=$HOME/anaconda3/bin/python

# Setup
directory=$($python -c 'import read_params; print(read_params.get_directory())')
num_src=$(cat $directory/master.pixels|wc -l)
num_ls=5
$python -c "import setup; setup.create_directories($num_src,$num_ls)"

# Generate data
echo "Starting at "`date`
find . -name "linesearch" -delete
touch compute_data
/usr/local/bin/pbsdsh $python $PBS_O_WORKDIR/data_forward.py

# Clean-up
find . -name "compute_data" -delete
echo "Finished at "`date`
