export directory=$(python -c 'import read_params; print(read_params.get_directory())')
export nmasterpixels=$(wc -l < $directory/master.pixels)

python -c "import setup; setup.create_directories(\"$directory\",$nmasterpixels,3)"

start=$(date +%s)

find . -name "linesearch" -delete
touch compute_data

compute_forward(){
	forward=$directory/forward_src"$1"_ls00
	
	cp Spectral Instruction_src"$1"_ls00

	# If you're using anaconda's mpi, make sure to use the mpich version and not openmpi
    # the mpich version distributes correctly across processors, whereas the openmpi version 
    # launches binds rank 0 to processor 0, and launches multiple processes on the same core
	mpiexec -np 1 ./sparc $1 00 > $forward/out_data_forward

	[ -e $forward/vz_cc.fits ] && cp $forward/vz_cc.fits $directory/data/"$1".fits

	find $forward -name "*full*" -delete 

    find . -name Instruction_src"$1"_ls00 -delete
    find $directory/status -name forward_src"$1"_ls00 -delete
}

export -f compute_forward

parallel compute_forward ::: $(seq -f "%02g" $nmasterpixels)

find . -name compute_data -delete

echo "Duration: $((($(date +%s)-$start)/60)) minutes"
