[[ -e running_full ]] && echo "Full running, quitting" && exit
[[ -e linesearch ]] && echo "Linesearch already running, quitting" && exit

directory=`python -c 'import read_params; print(read_params.get_directory())'`
export directory

find $directory -name "compute_data" -exec rm -f {} \; 
find $directory -name "compute_synth" -exec rm -f {} \; 

iter=`find $directory/update -maxdepth 1 -name 'linesearch_[0-9][0-9]'|wc -l`
itername=`printf "%02d" $iter`

touch linesearch
start=$(date +%s)

numls=`find $directory/update -maxdepth 1 -regex "$directory/update/test_psi_[0-9]+.fits"|wc -l`
echo "Number of linesearches: "$numls

for lin in `seq 1 $numls`
do
    linzpd=`printf "%02d" $lin`
    [[ -e $directory/update/test_psi_"$lin".fits ]] && cp $directory/update/test_psi_"$lin".fits  $directory/model_psi_ls"$linzpd".fits
done

nmasterpixels=`wc -l < $directory/master.pixels`

########################################################################
#~ Main computation
########################################################################

compute_forward(){
    cp Spectral Instruction_src"$1"_ls"$2"
    forward=$directory/forward_src"$1"_ls"$2"

    # If you're using anaconda's mpi, make sure to use the mpich version and not openmpi
    # the mpich version distributes correctly across processors, whereas the openmpi version 
    # launches binds rank 0 to processor 0, and launches multiple processes on the same core

    mpiexec -np 1 ./sparc $1 $2 > $forward/out_forward

    find $forward -name "*full*" -delete 
    find $forward -name "*partial*" -delete
    find $directory/status -name forward_src"$1"_ls"$2" -delete
}

export -f compute_forward

parallel compute_forward ::: `seq -f "%02g" $nmasterpixels` ::: `seq -f "%02g" $numls`

########################################################################


for lin in `seq -f "%02g" 1 $numls`
do
    for src in `seq -f "%02g" 1 $((nmasterpixels))`
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


#~ find $directory/update -name "tested*" -exec rm -f {} \; 
#~ find $directory -name "update.fits" -exec rm -f {} \; 
find . -name "linesearch" -delete 
find $directory/status -name "forward*" -delete

find . -name "core.*" -delete 
find . -name "fort.*" -delete

echo "Duration: $((($(date +%s)-$start)/60)) minutes"