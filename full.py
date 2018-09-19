import os,shutil,glob,subprocess,read_params,sys,fnmatch
from datetime import datetime
from pathlib import Path
import multiprocessing

env=dict(os.environ, MPI_TYPE_MAX="1280280")

codedir=os.path.dirname(os.path.abspath(__file__))
HOME=os.environ["HOME"]

datadir=read_params.get_directory()

# check the number of iterations that have been done
# The current iteration is the number already done plus one
def get_iter_no():
    updatedir=os.path.join(datadir,"update")
    # Count the number of misfit_xx files
    return len(fnmatch.filter(os.listdir(updatedir),'misfit_[0-9][0-9]'))
iterno=get_iter_no()



# def safecopy(a,b):
#     try: shutil.copyfile(a,b)
#     except IOError as e:
#         sys.stderr.write("Could not copy "+a+" to "+b+"; "+e.args(1)+"\n")
#         sys.stderr.flush()
        
# def safemkdir(a):
#     if not os.path.exists(a):
#         try: os.makedirs(a)
#         except OSError:
#             if e.errno == 17: pass
#             else: print(e)


mpipath = os.path.join(HOME,"anaconda3/bin/mpiexec")
if not Path(mpipath).exists() : 
    print("Could not find mpi, check mpipath specified as {}".format(mpipath))
    quit()



def compute_forward_adjoint_kernel(src):

    forward="forward_src{:02d}_ls00".format(src)
    adjoint="adjoint_src{:02d}".format(src)
    kernel="kernel"

    Spectral=os.path.join(codedir,"Spectral")
    Adjoint=os.path.join(codedir,"Adjoint")
    Instruction=os.path.join(codedir,"Instruction_src{:02d}_ls00".format(src))
    
    modes={'0':'fmode'}
    for pmodeno in range(1,8): modes.update({str(pmodeno):'p'+str(pmodeno)+'mode'})
    modes['8']='first_bounce_pmode'
    
    ridge_filters = read_params.get_modes_used()
    
    # for ridge_filter in ridge_filters:
    #     if os.path.exists(os.path.join(datadir,forward,"ttdiff."+ridge_filter)):
    #         safecopy(os.path.join(datadir,forward,"ttdiff."+ridge_filter),
    #                 os.path.join(datadir,forward,"ttdiff_prev."+ridge_filter))
    
    sparccmd = mpipath+" -np 1 ./sparc {:02d} 00".format(src)

    print("Starting computation for src {:02d}".format(src))
    
    ####################################################################
    #~ Forward
    ####################################################################
    
    shutil.copyfile(Spectral,Instruction)

    with open(os.path.join(datadir,forward,"out"+forward),'w') as outfile:
        fwd=subprocess.call(sparccmd.split(),stdout=outfile,env=env,cwd=codedir)
    
    assert fwd==0,"Error in running forward"
    print("Finished forward computation for src {:02d}".format(src))
    
    directory = os.path.join(datadir,"tt","iter{:02d}".format(iterno),"windows{:02d}".format(src))
    if not Path(directory).exists(): os.makedirs(directory)

    shutil.copyfile(os.path.join(datadir,forward,"vz_cc.fits"),
                    os.path.join(datadir,"tt","iter{:02d}".format(iterno),"vz_cc_src{:02d}.fits".format(src)))                    

    for ridge_filter in ridge_filters:
        shutil.copyfile(os.path.join(datadir,forward,"ttdiff."+ridge_filter),
                        os.path.join(datadir,"tt","iter{:02d}".format(iterno),
                        "ttdiff_src{:02d}.{}".format(src,modes.get(ridge_filter,ridge_filter)) ))
                        
    
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
for src in range(1,num_src+1):
    pool.apply_async(compute_forward_adjoint_kernel,args=(src,))

pool.close()
pool.join()

t_end = datetime.now()

delta_t = t_end - t_start

print("Computation took {}".format(str(delta_t)))

##########################################################################################


# Concatenate misfit files from individual sources
with open(os.path.join(datadir,"kernel","misfit_{:02d}".format(iterno)),"w") as misfitfile:
    for src in range(1,num_src+1):
        misfit_src_file_name = os.path.join(datadir,"kernel","misfit_{:02d}_00".format(src))
        if Path(misfit_src_file_name).exists():
            with open(misfit_src_file_name,"r") as misfit_src_file:
                misfitfile.write(misfit_src_file.read())
            os.remove(misfit_src_file_name)


with open(os.path.join(datadir,"kernel","misfit_all_{:02d}".format(iterno)),"w") as misfitfile:
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

if os.path.exists(os.path.join(datadir,"vx_00.fits")):
    shutil.copyfile(os.path.join(datadir,"vx_00.fits"),
            os.path.join(datadir,"update","vx_{:02d}.fits".format(iterno)))

if os.path.exists(os.path.join(datadir,"vz_00.fits")):
    shutil.copyfile(os.path.join(datadir,"vz_00.fits"),
            os.path.join(datadir,"update","vz_{:02d}.fits".format(iterno)))

# Some more clean-up
for f in glob.glob("core.*"): os.remove(f)
for f in glob.glob("fort.*"): os.remove(f)

if Path("running_full").exists(): os.remove("running_full")