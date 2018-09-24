import os,shutil,glob,subprocess,read_params,sys,fnmatch
from datetime import datetime
from pathlib import Path
import multiprocessing

env=dict(os.environ, MPI_TYPE_MAX="1280280")

codedir=os.path.dirname(os.path.abspath(__file__))
HOME=os.environ["HOME"]

datadir=read_params.get_directory()
updatedir=os.path.join(datadir,"update")

# check the number of iterations that have been done
# The current iteration is the number already done plus one
def get_iter_no():
        # Count the number of misfit_xx files
    return len(fnmatch.filter(os.listdir(updatedir),'misfit_[0-9][0-9]'))
iterno=get_iter_no()

# travel time shifts for individual filter gets saved here
tt_dir = os.path.join(datadir,"tt","iter{:02d}".format(iterno))
if not os.path.exists(tt_dir):
    os.makedirs(tt_dir)


def compute_forward_adjoint_kernel(src):

    forward="forward_src{:02d}_ls00".format(src)
    adjoint="adjoint_src{:02d}".format(src)
    kernel="kernel"

    Spectral=os.path.join(codedir,"Spectral")
    Adjoint=os.path.join(codedir,"Adjoint")
    Instruction=os.path.join(codedir,"Instruction_src{:02d}_ls00".format(src))
    
    # If you're using anaconda's mpi, make sure to use the mpich version and not openmpi
    # the mpich version distributes correctly across processors, whereas the openmpi version 
    # launches binds rank 0 to processor 0, and launches multiple processes on the same core
    sparccmd = "mpiexec -np 1 ./sparc {:02d} 00".format(src)

    print("Starting computation for src {:02d}".format(src))
    
    ####################################################################
    #~ Forward
    ####################################################################
    
    shutil.copyfile(Spectral,Instruction)

    with open(os.path.join(datadir,forward,"out"+forward),'w') as outfile:
        fwd=subprocess.call(sparccmd.split(),stdout=outfile,env=env,cwd=codedir)
    
    assert fwd==0,"Error in running forward"
    
    print("Finished forward computation for src {:02d}".format(src))
    
    for ttfile in fnmatch.filter(os.listdir(os.path.join(datadir,forward)),"ttdiff.[0-9]" ) :
        shutil.copyfile(os.path.join(datadir,forward,ttfile),
                            os.path.join(tt_dir,ttfile))
                        
    
    ####################################################################
    #~ Adjoint
    ####################################################################

    shutil.copyfile(Adjoint,Instruction)

    with open(os.path.join(datadir,adjoint,"out"+adjoint),'w') as outfile:
        adj=subprocess.call(sparccmd.split(),stdout=outfile,env=env,cwd=codedir)

    assert adj==0,"Error in running adjoint"

    print("Finished adjoint computation for src {:02d}".format(src))
            
    ####################################################################
    #~ Kernel
    ####################################################################

    with open(os.path.join(datadir,kernel,"out_kernel{:02d}".format(src)),'w') as outfile:
        kern=subprocess.call(sparccmd.split(),stdout=outfile,env=env,cwd=codedir)

    assert kern==0,"Error in computing kernel"

    print("Finished kernel computation for src {:02d}".format(src))

    ####################################################################
    # if you've reached here then everything has finished correctly. Cleaning up
    ####################################################################
    print("Finished computation for src {:02d}, cleaning up".format(src))

    file_to_remove=os.path.join(datadir,"status","forward_src{:02d}_ls00".format(src))
    if os.path.exists(file_to_remove): os.remove(file_to_remove)

    file_to_remove=os.path.join(datadir,"status","adjoint_src{:02d}".format(src))
    if os.path.exists(file_to_remove): os.remove(file_to_remove)

    file_to_remove=os.path.join(datadir,"status","kernel_src{:02d}".format(src))
    if os.path.exists(file_to_remove): os.remove(file_to_remove)

    for f in glob.glob(os.path.join(datadir,forward,"*full*")): os.remove(f)
    for f in glob.glob(os.path.join(datadir,forward,"*partial*")): os.remove(f)
    for f in glob.glob(os.path.join(datadir,adjoint,"*full*")): os.remove(f)
    for f in glob.glob(os.path.join(datadir,adjoint,"*partial*")): os.remove(f)


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

pool = multiprocessing.Pool(processes=number_of_processors_to_use)
pool.starmap_async(compute_forward_adjoint_kernel,range(1,num_src+1))

pool.close()
pool.join()

t_end = datetime.now()

delta_t = t_end - t_start

print("Computation took {}".format(str(delta_t)))

##########################################################################################


# Concatenate misfit files from individual sources
with open(os.path.join(updatedir,"misfit_{:02d}".format(iterno)),"w") as misfitfile:
    for src in range(1,num_src+1):
        misfit_src_file_name = os.path.join(datadir,"kernel","misfit_{:02d}_00".format(src))
        if Path(misfit_src_file_name).exists():
            with open(misfit_src_file_name,"r") as misfit_src_file:
                misfitfile.write(misfit_src_file.read())
            os.remove(misfit_src_file_name)


with open(os.path.join(updatedir,"misfit_all_{:02d}".format(iterno)),"w") as misfitfile:
    for src in range(1,num_src+1):
        misfit_src_file_name = os.path.join(datadir,"kernel","misfit_all_{:02d}_00".format(src))
        if Path(misfit_src_file_name).exists():
            with open(misfit_src_file_name,"r") as misfit_src_file:
                misfitfile.write(misfit_src_file.read())
            os.remove(misfit_src_file_name)


# Copy the model used to the update directory
if os.path.exists(os.path.join(datadir,"model_psi_ls00.fits")):
    shutil.copyfile(os.path.join(datadir,"model_psi_ls00.fits"),
            os.path.join(datadir,"update","model_psi_{:02d}.fits".format(iterno)))

# Some more clean-up
for f in glob.glob("core.*"): os.remove(f)
for f in glob.glob("fort.*"): os.remove(f)

if Path("running_full").exists(): os.remove("running_full")