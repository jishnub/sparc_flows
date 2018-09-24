import os,shutil,glob,subprocess,read_params,sys,fnmatch
from datetime import datetime
from pathlib import Path
import multiprocessing

def compute_forward_adjoint_kernel(src):

    forward="forward_src{:02d}_ls00".format(src)
    adjoint="adjoint_src{:02d}".format(src)
    kernel="kernel"

    Instruction="Instruction_src{:02d}_ls00".format(src)
    
    # If you're using anaconda's mpi, make sure to use the mpich version and not openmpi
    # the mpich version distributes correctly across processors, whereas the openmpi version 
    # launches binds rank 0 to processor 0, and launches multiple processes on the same core
    sparccmd = "mpiexec -np 1 ./sparc {:02d} 00".format(src)

    print("Starting computation for src {:02d}".format(src))
    
    ####################################################################
    #~ Forward
    ####################################################################
    
    shutil.copyfile("Spectral",Instruction)

    with open(datadir/forward/"out_forward",'w') as outfile:
        fwd=subprocess.call(sparccmd.split(),stdout=outfile,env=env)
    
    assert fwd==0,"Error in running forward"
    
    print("Finished forward computation for src {:02d}".format(src))
    
    for ttfile in fnmatch.filter(os.listdir(os.path.join(datadir,forward)),"ttdiff.[0-9]" ) :
        shutil.copyfile(os.path.join(datadir,forward,ttfile),
                            os.path.join(tt_dir,ttfile))
                        
    
    ####################################################################
    #~ Adjoint
    ####################################################################

    shutil.copyfile("Adjoint",Instruction)

    with open(datadir/adjoint/"out_adjoint",'w') as outfile:
        adj=subprocess.call(sparccmd.split(),stdout=outfile,env=env)

    assert adj==0,"Error in running adjoint"

    print("Finished adjoint computation for src {:02d}".format(src))
            
    ####################################################################
    #~ Kernel
    ####################################################################

    with open(datadir/kernel/"out_kernel{:02d}".format(src),'w') as outfile:
        kern=subprocess.call(sparccmd.split(),stdout=outfile,env=env)

    assert kern==0,"Error in computing kernel"

    print("Finished kernel computation for src {:02d}".format(src))

    ####################################################################
    # if you've reached here then everything has finished correctly. Cleaning up
    ####################################################################
    print("Finished computation for src {:02d}, cleaning up".format(src))

    file_to_remove=datadir/"status"/"forward_src{:02d}_ls00".format(src)
    if file_to_remove.exists(): os.remove(file_to_remove)

    file_to_remove=datadir/"status"/"adjoint_src{:02d}".format(src)
    if file_to_remove.exists(): os.remove(file_to_remove)

    file_to_remove=datadir/"status"/"kernel_src{:02d}".format(src)
    if file_to_remove.exists(): os.remove(file_to_remove)

    for f in (datadir/forward).glob("*full*"): os.remove(f)
    for f in (datadir/forward).glob("*partial*"): os.remove(f)
    for f in (datadir/adjoint).glob("*full*"): os.remove(f)
    for f in (datadir/forward).glob("*partial*"): os.remove(f)


if __name__ == "__main__":

    env=dict(os.environ, MPI_TYPE_MAX="1280280")

    datadir=Path(read_params.get_directory())

    # check the number of iterations that have been done
    # The current iteration is the number already done plus one
    iterno=len(fnmatch.filter(os.listdir(datadir/"update"),'misfit_[0-9][0-9]'))

    # travel time shifts for individual filter gets saved here
    tt_dir = datadir/"tt"/"iter{:02d}".format(iterno)
    tt_dir.mkdir(parents=True,exist_ok=True)

    # Some status checks
    if Path("linesearch").exists(): 
        print("Linesearch is running, quitting. Remove the 'linesearch' file if it's not")
        quit()

    if Path("running_full").exists(): 
        print("Full already running, quitting")
        quit()

    Path("running_full").touch()

    if Path("compute_data").exists(): os.remove("compute_data")
    if Path("compute_synth").exists(): os.remove("compute_synth")

    ##########################################################################################
    # Start the parallel computation
    with open(os.path.join(datadir,'master.pixels'),'r') as mp:
        num_src=sum(1 for _ in mp)

    total_number_of_processors = multiprocessing.cpu_count()
    number_of_processors_to_use = total_number_of_processors - 1 or 1

    print("Using {:d} processors out of {:d}".format(number_of_processors_to_use,total_number_of_processors))

    t_start = datetime.now()

    with multiprocessing.Pool(processes=number_of_processors_to_use) as pool:
        pool.map(compute_forward_adjoint_kernel,range(1,num_src+1))

    delta_t = datetime.now() - t_start

    print("Computation took {}".format(str(delta_t)))

    ##########################################################################################


    # Concatenate misfit files from individual sources
    with open(datadir/"update"/"misfit_{:02d}".format(iterno),"w") as misfitfile:
        for src in range(1,num_src+1):
            misfit_src_file_name = datadir/"kernel"/"misfit_{:02d}_00".format(src)
            if misfit_src_file_name.exists():
                with open(misfit_src_file_name,"r") as misfit_src_file:
                    misfitfile.write(misfit_src_file.read())
                os.remove(misfit_src_file_name)


    with open(datadir/"update"/"misfit_all_{:02d}".format(iterno),"w") as misfitfile:
        for src in range(1,num_src+1):
            misfit_src_file_name = datadir/"kernel"/"misfit_all_{:02d}_00".format(src)
            if misfit_src_file_name.exists():
                with open(misfit_src_file_name,"r") as misfit_src_file:
                    misfitfile.write(misfit_src_file.read())
                os.remove(misfit_src_file_name)


    # Copy the model used to the update directory
    if (datadir/"model_psi_ls00.fits").exists():
        shutil.copyfile(datadir/"model_psi_ls00.fits",
                datadir/"update"/"model_psi_{:02d}.fits".format(iterno))

    # Some more clean-up
    for f in glob.glob("core.*"): os.remove(f)
    for f in glob.glob("fort.*"): os.remove(f)

    if Path("running_full").exists(): os.remove("running_full")