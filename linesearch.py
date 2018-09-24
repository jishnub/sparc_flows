import os,shutil,subprocess,read_params
import glob,fnmatch
from datetime import datetime
import multiprocessing
from pathlib import Path
import itertools

env=dict(os.environ, MPI_TYPE_MAX="1280280")

#~ Get the current directory, might need to pass it to the subprocess
codedir=os.path.dirname(os.path.abspath(__file__))
HOME=env["HOME"]

datadir=read_params.get_directory()
updatedir = os.path.join(datadir,"update")


#~ Get the total number of source pixels to determine number of linesearches necessary
with open(os.path.join(datadir,'master.pixels'),'r') as mpixfile:
    num_src=sum(1 for _ in mpixfile)

num_ls_per_src = len(fnmatch.filter(os.listdir(datadir),"forward_src01_ls[0-9][1-9]"))

updatedir=os.path.join(datadir,"update")

iterno = len(fnmatch.filter(os.listdir(updatedir),'misfit_[0-9][0-9]'))-1

# This function carries out one forward calculation for one source.
def compute_forward(src,linesearch_no):

    Spectral = os.path.join(codedir,"Spectral")
    Instruction = os.path.join(codedir,"Instruction_src{:02d}_ls{:02d}".format(src,linesearch_no))
    forward = os.path.join(datadir,"forward_src{:02d}_ls{:02d}".format(src,linesearch_no))

    shutil.copyfile(Spectral,Instruction)

    # If you're using anaconda's mpi, make sure to use the mpich version and not openmpi
    # the mpich version distributes correctly across processors, whereas the openmpi version 
    # launches binds rank 0 to processor 0, and launches multiple processes on the same core
    sparccmd="mpiexec -np 1 ./sparc {:02d} {:02d}".format(src,linesearch_no)

    with open(os.path.join(datadir,forward,"out_forward".format(linesearch_no)),'w') as outfile:
        fwd=subprocess.call(sparccmd.split(),stdout=outfile,env=env,cwd=codedir)

    assert fwd==0,"Error in running linesearch for lsno {:02d}".format(linesearch_no)

    ########################################################################################
    # If you've reached this point then the code has finished running correctly. CLean up now
    ########################################################################################

    partialfiles=glob.glob(os.path.join(forward,"*partial*"))
    for f in partialfiles: os.remove(f)

    fullfiles=glob.glob(os.path.join(forward,"*full*"))
    for f in fullfiles: os.remove(f)

    file_to_remove=os.path.join(datadir,"status","forward_src{:02d}_ls{:02d}".format(src,linesearch_no))
    if os.path.exists(file_to_remove): os.remove(file_to_remove)

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
    test_model = os.path.join(updatedir,"test_psi_ls{:d}".format(linesearch_no+1))
    ls_model = os.path.join(datadir,"model_psi_ls{:02d}".format(linesearch_no+1))
    
    if os.path.exists(test_model):
        shutil.copyfile(test_model,ls_model)

#############################################################################################
# Start the parallel computation
Path("linesearch").touch()

total_number_of_processors = multiprocessing.cpu_count()
number_of_processors_to_use = total_number_of_processors - 1 or 1

t_start = datetime.now()

pool = multiprocessing.Pool(processes=number_of_processors_to_use)

pool.starmap_async(compute_forward,itertools.product(range(1,num_src+1),range(1,num_ls_per_src+1)))

pool.close()
pool.join()

t_end = datetime.now()

delta_t = t_end - t_start

print("Computation took {}".format(str(delta_t)))

##############################################################################################

# Concatenate misfit files from individual sources
with open(os.path.join(updatedir,"linesearch_{:02d}".format(iterno)),"w") as misfitfile:
    for src in range(num_src):
        for linesearch_no in range(num_ls_per_src):
            
            misfit_src_file_name = os.path.join(datadir,"kernel","misfit_{:02d}_{:02d}".format(src+1,linesearch_no+1))
            
            if Path(misfit_src_file_name).exists():
                with open(misfit_src_file_name,"r") as misfit_src_file:
                    misfitfile.write(misfit_src_file.read())
                
                os.remove(misfit_src_file_name)

# with open(os.path.join(updatedir,"linesearch_all_{:02d}".format(iterno)),"w") as misfitfile:
#     for src in range(num_src):
#         for linesearch_no in range(num_ls_per_src):
            
#             misfit_src_file_name = os.path.join(datadir,"kernel","misfit_all_{:02d}_{:02d}".format(src+1,linesearch_no+1))
            
#             if Path(misfit_src_file_name).exists():
#                 with open(misfit_src_file_name,"r") as misfit_src_file:
#                     misfitfile.write(misfit_src_file.read())
                
#                 os.remove(misfit_src_file_name)

# More clean-ups
for linesearch_no in range(num_ls_per_src):

    ls_model = os.path.join(datadir,"model_psi_ls{:02d}".format(linesearch_no+1))
    if os.path.exists(ls_model):   os.remove(ls_model)

# Some more clean-up
for f in glob.glob("core.*"): os.remove(f)
for f in glob.glob("fort.*"): os.remove(f)

if Path("linesearch").exists(): os.remove("linesearch")