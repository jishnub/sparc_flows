[ ! -e sparc ] && echo "The executable sparc doesn't exist, maybe you've not run make?" && exit 1

[[ -e linesearch ]] && echo "Linesearch running, quitting. Remove the 'linesearch' file if it's not" && exit

[[ -e running_full ]] && echo "Full already running, quitting" && exit

touch running_full

start=$(date +%s)

export directory=$(python -c 'import read_params; print(read_params.get_directory())')

find . -name "compute_data" -delete
find . -name "compute_synth" -delete

iter=$(find $directory/update -maxdepth 1 -name 'misfit_[0-9][0-9]'|wc -l)
export itername=$(printf "%02d" $iter)
[ ! -d $directory/tt/iter$itername ] && mkdir -p $directory/tt/iter$itername

export nmasterpixels=$(wc -l < $directory/master.pixels)

###########################################################################
# Main computation
###########################################################################

compute_forward_adjoint_kernel(){
	# forward
	cp Spectral Instruction_src"$1"_ls00
    forward=$directory/forward_src"$1"_ls00
    adjoint=$directory/adjoint_src"$1"

    # If you're using anaconda's mpi, make sure to use the mpich version and not openmpi
    # the mpich version distributes correctly across processors, whereas the openmpi version 
    # launches binds rank 0 to processor 0, and launches multiple processes on the same core
    # mpiexec -np 1 ./sparc $1 00 > $forward/out_forward
    
    for ttfile in $(ls $forward/ttdiff.*); do
    	extension="${ttfile##*.}"
    	filename=$(basename $(echo "${ttfile%.*}"))
    	cp -v $ttfile $directory/tt/iter$itername/"$filename"_src"$1"."$extension"

    done

    #adjoint
    cp Adjoint Instruction_src"$1"_ls00
    # mpiexec -np 1 ./sparc $1 00 > $adjoint/out_adjoint

    #kernel
	# mpiexec -np 1 ./sparc $1 00 > $directory/kernel/out_kernel"$1" 
    
    find $forward -name "*full*" -delete 
    find $forward -name "*partial*" -delete

    find $adjoint -name "*full*" -delete
    find $adjoint -name "*partial*" -delete

    find $directory/status -name forward_src"$1"_ls00 -delete
    find $directory/status -name adjoint_src"$1" -delete
    find $directory/status -name kernel_src"$1" -delete
}

export -f compute_forward_adjoint_kernel

parallel compute_forward_adjoint_kernel ::: $(seq -f "%02g" $nmasterpixels)

# ###########################################################################

# Concatenate misfit files only after everything is complete

for src in $(seq -f "%02g" 1 $((nmasterpixels)))
do
    [[ -e $directory/kernel/misfit_"$src"_00 ]]  && \
    cat $directory/kernel/misfit_"$src"_00 >> $directory/update/misfit_$itername &&\
    rm $directory/kernel/misfit_"$src"_00
    
    [[ -e $directory/kernel/misfit_all_"$src"_00 ]]  && \
    cat $directory/kernel/misfit_all_"$src"_00 >> $directory/update/misfit_all_$itername &&\
    rm $directory/kernel/misfit_all_"$src"_00
done

[[ -e $directory/model_psi_ls00.fits ]] && cp $directory/model_psi_ls00.fits $directory/update/model_psi_"$itername".fits

find . -name "core.*" -delete
find . -name "fort.*" -delete

find . -name running_full -delete

echo "Duration: $((($(date +%s)-$start)/60)) minutes"