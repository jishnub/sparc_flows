import os,sys,shutil,glob,subprocess,read_params,multiprocessing
from datetime import datetime
from pathlib import Path
import setup


def compute_data(src):

    forward="forward_src{:02d}_ls00".format(src)
    Spectral="Spectral"
    Instruction="Instruction_src{:02d}_ls00".format(src)
    
    shutil.copyfile(Spectral,Instruction)

    # If you're using anaconda's mpi, make sure to use the mpich version and not openmpi
    # the mpich version distributes correctly across processors, whereas the openmpi version 
    # launches binds rank 0 to processor 0, and launches multiple processes on the same core
    sparccmd="mpiexec -np 1 ./sparc {:02d} 00".format(src)
   
    print("Starting computation for src {:02d}".format(src))

    with open(datadir/forward/"out_data_forward",'w') as outfile:
        fwd=subprocess.call(sparccmd.split(),stdout=outfile,env=env)

    print("Finished computation for src {:02d}, cleaning up".format(src))

    ####################################################################
    # if you've reached here then everything has finished correctly. Cleaning up
    ####################################################################
    
    for f in (datadir/forward).glob("*full*"): os.remove(f)
    
    if (datadir/forward/"vz_cc.fits").exists():
        shutil.copyfile(datadir/forward/"vz_cc.fits",
                        datadir/"data"/"{:02d}.fits".format(src))

    if Instruction.exists(): os.remove(Instruction)
    
    file_to_remove=datadir/"status"/"forward_src{:02d}_ls00".format(src)
    if file_to_remove.exists(): os.remove(file_to_remove)

    return True

if __name__ == "__main__":

    env=dict(os.environ, MPI_TYPE_MAX="1280280")

    datadir=Path(read_params.get_directory())

    try:
        with open(datadir/'master.pixels','r') as mp:
            num_src=sum(1 for _ in mp)
    except FileNotFoundError:
        print("No master.pixels file found, create it before running this code")
        quit()

    num_ls_per_src = 3
    setup.create_directories(datadir,num_src,num_ls_per_src)

    # Some status checks
    if Path("linesearch").exists(): os.remove("linesearch")
    if not Path("compute_data").exists(): Path("compute_data").touch()

    #########################################################################################################
    # Start the parallel job
    total_number_of_processors = multiprocessing.cpu_count()
    number_of_processors_to_use = total_number_of_processors - 1 or 1

    print("Using {:d} processors out of {:d}".format(number_of_processors_to_use,total_number_of_processors))

    t_start = datetime.now()

    with multiprocessing.Pool(processes=number_of_processors_to_use) as pool:
        pool.map(compute_data,range(1,num_src+1))

    delta_t = datetime.now() - t_start

    print("Computation took {}".format(str(delta_t)))

    #########################################################################################################

    if Path("compute_data").exists(): os.remove("compute_data")



