import os,sys,shutil,glob,subprocess,read_params,multiprocessing
from datetime import datetime
from pathlib import Path
import setup

env=dict(os.environ, MPI_TYPE_MAX="1280280")

codedir=os.path.dirname(os.path.abspath(__file__))
HOME=os.environ["HOME"]

datadir=read_params.get_directory()

try:
    with open(os.path.join(datadir,'master.pixels'),'r') as mp:
        num_src=sum(1 for _ in mp)
except FileNotFoundError:
    print("No master.pixels file found, create it before running this code")
    quit()

setup.create_directories(datadir,num_src,3)

mpipath=os.path.join(HOME,"anaconda3/bin/mpiexec")
if not Path(mpipath).exists() : 
    print("Could not find mpi, check mpipath specified as {}".format(mpipath))
    quit()

def compute_data(src):

    forward="forward_src{:02d}_ls00".format(src)
    Spectral=os.path.join(codedir,"Spectral")
    Instruction=os.path.join(codedir,"Instruction_src{:02d}_ls00".format(src))
    
    shutil.copyfile(Spectral,Instruction)

    sparccmd=mpipath+" -np 1 ./sparc {:02d} 00".format(src)
   
    print("Starting computation for src {:02d}".format(src))

    with open(os.path.join(datadir,forward,"out_data_forward"),'w') as outfile:
        fwd=subprocess.call(sparccmd.split(),stdout=outfile,env=env,cwd=codedir)

    print("Finished computation for src {:02d}, cleaning up".format(src))

    partialfiles=glob.glob(os.path.join(datadir,forward,"*partial*"))
    for f in partialfiles: os.remove(f)
    
    fullfiles=glob.glob(os.path.join(datadir,forward,"*full*"))
    for f in fullfiles: os.remove(f)
    
    if os.path.exists(os.path.join(datadir,forward,"vz_cc.fits")):
        shutil.copyfile(os.path.join(datadir,forward,"vz_cc.fits"),
                        os.path.join(datadir,"data","{:02d}.fits".format(src)))

    if os.path.exists(Instruction): os.remove(Instruction)
    
    file_to_remove=os.path.join(datadir,"status","forward_src{:02d}_ls00".format(src))
    if os.path.exists(file_to_remove): os.remove(file_to_remove)

# Some status checks
if Path("linesearch").exists(): os.remove("linesearch")
if not Path("compute_data").exists(): Path("compute_data").touch()

#########################################################################################################
# Start the parallel job
total_number_of_processors = multiprocessing.cpu_count()
number_of_processors_to_use = total_number_of_processors - 1 or 1

print("Using {:d} processors out of {:d}".format(number_of_processors_to_use,total_number_of_processors))

t_start = datetime.now()

pool = multiprocessing.Pool(processes=number_of_processors_to_use)
for src in range(1,num_src+1):
    pool.apply_async(compute_data,args=(src,))

pool.close()
pool.join()

t_end = datetime.now()

delta_t = t_end - t_start

print("Computation took {}".format(str(delta_t)))

#########################################################################################################

if Path("compute_data").exists(): os.remove("compute_data")



