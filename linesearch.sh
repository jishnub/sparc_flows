#PBS -l nodes=3:ppn=24
#PBS -o  output-linesearch
#PBS -e  error-linesearch
#PBS -l walltime=12:00:00
cd $PBS_O_WORKDIR
echo $PBS_JOBID
export TERM=xterm

[[ -e running_full ]] && echo "Full running, quitting" && exit
[[ -e linesearch ]] && echo "Linesearch already running, quitting" && exit

directory=$($HOME/anaconda3/bin/python -c 'import read_params; print(read_params.get_directory())')

find $directory -name "compute_data" -exec rm -f {} \; 
find $directory -name "compute_synth" -exec rm -f {} \; 

iter=$(find $directory/update -maxdepth 1 -name 'linesearch_[0-9][0-9]'|wc -l)
itername=$(printf "%02d" $iter)

touch linesearch
echo "Starting iterations at "`date`

[[ numls -eq 0 ]] && numls=$(find $directory/update -maxdepth 1 -regex "$directory/update/test_psi_[0-9]+.fits"|wc -l)
echo "Number of linesearches: "$numls

for lin in $(seq 1 $numls)
do
    linzpd=$(printf "%02d" $lin)
    [[ -e $directory/update/test_psi_"$lin".fits ]] && cp $directory/update/test_psi_"$lin".fits  $directory/model_psi_ls"$linzpd".fits
done

########################################################################
#~ Main computation
########################################################################

/usr/local/bin/pbsdsh $HOME/anaconda3/bin/python $PBS_O_WORKDIR/linesearch.py

########################################################################

nmasterpixels=$(wc -l < $directory/master.pixels)

for lin in $(seq -f "%02g" 1 $numls)
do
    for src in $(seq -f "%02g" 1 $nmasterpixels)
    do
        [[ -e $directory/kernel/misfit_"$src"_"$lin" ]]  && \
        cat $directory/kernel/misfit_"$src"_"$lin" >> $directory/update/linesearch_$itername &&\
        rm $directory/kernel/misfit_"$src"_"$lin"
        
        [[ -e $directory/kernel/misfit_all_"$src"_"$lin" ]] && \
        cat $directory/kernel/misfit_all_"$src"_"$lin" >> $directory/update/linesearch_all_$itername &&\
        rm $directory/kernel/misfit_all_"$src"_"$lin"
    
    done
    
    [[ -e $directory/model_psi_ls"$lin".fits ]] && rm $directory/model_psi_ls"$lin".fits
    
done


find . -name "linesearch" -delete 
find $directory/status -name "forward*" -delete

find . -name "core.*" -delete 
find . -name "fort.*" -delete

echo "Finished at "`date`
