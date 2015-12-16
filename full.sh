#!/bin/bash
#PBS -N  full_f_3hr_smooth
#PBS -l nodes=1:ppn=24
#PBS -o  output-full
#PBS -e  error-full
#PBS -l walltime=12:00:00
cd $PBS_O_WORKDIR
echo $PBS_JOBID
export TERM=xterm

[[ -e linesearch ]] && echo "Linesearch running, quitting. Remove the 'linesearch' file if it's not"\
 && exit

[[ -e running_full ]] && echo "Full already running, quitting" && exit

touch running_full

echo "Starting at "`date`

directory=`python -c 'import read_params; print read_params.get_directory()'`


find . -name "compute_data" -delete
find . -name "compute_synth" -delete

iter=`find $directory/update -maxdepth 1 -name 'misfit_[0-9][0-9]'|wc -l`
itername=`printf "%02d" $iter`

/usr/local/bin/pbsdsh python $PBS_O_WORKDIR/full.py

# Concatenate misfit files only after everything is complete
nmasterpixels=`wc -l < $directory/master.pixels`
for src in `seq -f "%02g" 1 $((nmasterpixels))`
do
    [[ -e $directory/kernel/misfit_"$src"_00 ]]  && \
    cat $directory/kernel/misfit_"$src"_00 >> $directory/update/misfit_$itername &&\
    rm $directory/kernel/misfit_"$src"_00
    
    [[ -e $directory/kernel/misfit_all_"$src"_00 ]]  && \
    cat $directory/kernel/misfit_all_"$src"_00 >> $directory/update/misfit_all_$itername &&\    
    rm $directory/kernel/misfit_all_"$src"_00
done

find $directory/status -name "forward*" -delete
find $directory/status -name "adjoint*" -delete
find $directory/status -name "kernel*" -delete

[[ -e $directory/model_psi_ls00.fits ]] && cp $directory/model_psi_ls00.fits $directory/update/model_psi_"$itername".fits
[[ -e $directory/model_c_ls00.fits ]] && cp $directory/model_c_ls00.fits $directory/update/model_c_"$itername".fits
[[ -e vx_00.fits ]] && cp vx_00.fits "$directory"/update/vx_"$itername".fits
[[ -e vz_00.fits ]] && cp vz_00.fits "$directory"/update/vz_"$itername".fits

find . -name "core.*" -delete
find . -name "fort.*" -delete

rm running_full
echo "Finished at "`date`
