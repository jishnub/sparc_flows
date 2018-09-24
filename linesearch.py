import os,shutil,subprocess,read_params
import glob,fnmatch
from datetime import datetime
import multiprocessing
from pathlib import Path
import itertools



# This function carries out one forward calculation for one source.
def compute_forward(src,linesearch_no):

    Instruction = "Instruction_src{:02d}_ls{:02d}".format(src,linesearch_no)
    forward = "forward_src{:02d}_ls{:02d}".format(src,linesearch_no)

    shutil.copyfile("Spectral",Instruction)

    # If you're using anaconda's mpi, make sure to use the mpich version and not openmpi
    # the mpich version distributes correctly across processors, whereas the openmpi version 
    # launches binds rank 0 to processor 0, and launches multiple processes on the same core
    sparccmd="mpiexec -np 1 ./sparc {:02d} {:02d}".format(src,linesearch_no)

    with open(datadir/forward/"out_forward",'w') as outfile:
        fwd=subprocess.call(sparccmd.split(),stdout=outfile,env=env)

    assert fwd==0,"Error in running linesearch for lsno {:02d}".format(linesearch_no)

    ########################################################################################
    # If you've reached this point then the code has finished running correctly. Clean up now
    ########################################################################################

    for f in (datadir/forward).glob("*full*"): os.remove(f)

    file_to_remove=datadir/"status"/"forward_src{:02d}_ls{:02d}".format(src,linesearch_no)
    if file_to_remove.exists(): os.remove(file_to_remove)


if __name__ == "__main__":

    env=dict(os.environ, MPI_TYPE_MAX="1280280")

    datadir=Path(read_params.get_directory())

    #~ Get the total number of source pixels to determine number of linesearches necessary
    with open(datadir/'master.pixels','r') as mpixfile:
        num_src=sum(1 for _ in mpixfile)

    num_ls_per_src = len(fnmatch.filter(os.listdir(datadir),"forward_src01_ls[0-9][1-9]"))

    iterno = len(fnmatch.filter(os.listdir(datadir/"update"),'misfit_[0-9][0-9]'))-1


    # Some status checks

    if Path("linesearch").exists(): 
        print("Linesearch is already running, quitting.")
        quit()

    if Path("running_full").exists(): 
        print("Full is running, quitting. Remove the file 'running_full' if it isn't.")
        quit()

    if Path("compute_data").exists(): os.remove("compute_data")
    if Path("compute_synth").exists(): os.remove("compute_synth")

    # Copy the models to the correct locations
    for linesearch_no in range(num_ls_per_src):
        test_model = datadir/"update"/"test_psi_ls{:d}".format(linesearch_no+1)
        ls_model = datadir/"model_psi_ls{:02d}".format(linesearch_no+1)
        
        if test_model.exists():
            shutil.copyfile(test_model,ls_model)

    #############################################################################################
    # Start the parallel computation
    Path("linesearch").touch()

    total_number_of_processors = multiprocessing.cpu_count()
    number_of_processors_to_use = total_number_of_processors - 1 or 1

    t_start = datetime.now()

    print("Starting {} linesearches each for {} sources on {} processors".format(num_ls_per_src,num_src,number_of_processors_to_use))

    with multiprocessing.Pool(processes=number_of_processors_to_use) as pool:
        pool.starmap(compute_forward,itertools.product(range(1,num_src+1),range(1,num_ls_per_src+1)))

    # pool.close()
    # pool.join()

    delta_t = datetime.now() - t_start

    print("Computation took {}".format(str(delta_t)))

    ##############################################################################################

    # Concatenate misfit files from individual sources
    with open(datadir/"update"/"linesearch_{:02d}".format(iterno),"w") as misfitfile:
        for src in range(num_src):
            for linesearch_no in range(num_ls_per_src):
                
                misfit_src_file_name = datadir/"kernel"/"misfit_{:02d}_{:02d}".format(src+1,linesearch_no+1)
                
                if misfit_src_file_name.exists():
                    with open(misfit_src_file_name,"r") as misfit_src_file:
                        misfitfile.write(misfit_src_file.read())
                    
                    os.remove(misfit_src_file_name)

    # More clean-ups
    for linesearch_no in range(num_ls_per_src):

        ls_model = datadir/"model_psi_ls{:02d}".format(linesearch_no+1)
        if ls_model.exists():   os.remove(ls_model)

    # Some more clean-up
    for f in glob.glob("core.*"): os.remove(f)
    for f in glob.glob("fort.*"): os.remove(f)

    if Path("linesearch").exists(): os.remove("linesearch")