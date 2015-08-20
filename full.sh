#!/bin/bash
#PBS -N  flows_full
#PBS -l nodes=1:ppn=24
#PBS -o  output-full
#PBS -e  error-full
#PBS -l walltime=12:00:00
cd $PBS_O_WORKDIR
echo $PBS_JOBID
export TERM=xterm

echo "Starting at "`date`
source varlist.sh

find . -name "linesearch" -exec rm -f {} \; 
find . -name "compute_data" -exec rm -f {} \; 
find .-name "compute_synth" -exec rm -f {} \; 

iter=`find $directory/update -name 'misfit_[0-9][0-9]'|wc -l`
itername=`printf "%02d" $iter`

/usr/local/bin/pbsdsh python $PBS_O_WORKDIR/full.py

find $directory/status -name "forward*" -exec rm -f {} \;
find $directory/status -name "adjoint*" -exec rm -f {} \;
find $directory/status -name "kernel*" -exec rm -f {} \;

# Concatenate misfit files only after everything is complete
nmasterpixels=`wc -l < $directory/master.pixels`
for src in `seq -f "%02g" 1 $((nmasterpixels))`
do
    cat $directory/kernel/misfit_"$src"_00 >> $directory/update/misfit_$itername
    cat $directory/kernel/misfit_all_"$src"_00 >> $directory/update/misfit_all_$itername
    rm $directory/kernel/misfit_"$src"_00
    rm $directory/kernel/misfit_all_"$src"_00
done

cp $directory/model_psi_ls00.fits $directory/update/model_psi_"$itername".fits
cp $directory/model_c_ls00.fits $directory/update/model_c_"$itername".fits

find . -name "core.*" -exec rm -f {} \; 
find . -name "fort.*" -exec rm -f {} \; 

echo "Finished at "`date`
