import os,sys,shutil,glob,time,re,subprocess,read_params,fnmatch

env=dict(os.environ, MPI_TYPE_MAX="1280280")

#~ Get the current directory, might need to pass it to the subprocess
codedir=os.path.dirname(os.path.abspath(__file__))
HOME=env["HOME"]
env["LD_LIBRARY_PATH"]=os.path.join(HOME,"anaconda2/lib")

datadir=read_params.get_directory()

procno=int(env["PBS_VNODENUM"])
nodeno=int(env["PBS_NODENUM"])

#~ Get the total number of source pixels to determine number of linesearches necessary
with open(os.path.join(datadir,'master.pixels'),'r') as mpixfile:
    nsrc=sum(1 for _ in mpixfile)

no_of_ls_per_src=int(sys.argv[1])
total_no_of_jobs=nsrc*no_of_ls_per_src

if procno>=total_no_of_jobs:
    print "Stopping job on node",nodeno,"proc",procno,"at",time.strftime("%H:%M:%S")
    quit()

# if procno<1:
#     print "Processor no "+str(procno)+" temporary stop to check MPI issue"
#     exit()

linesearch_no=procno/nsrc+1
src_no=procno%nsrc+1

print "Running linesearch no",linesearch_no,"for src no",src_no,"on proc",procno+1,"node no",nodeno+1

ls_no=str(linesearch_no).zfill(2)
src=str(src_no).zfill(2)

updatedir=os.path.join(datadir,"update")
iter_no = len(fnmatch.filter(os.listdir(updatedir),'misfit_[0-9][0-9]'))-1

def compute_forward(linesearch_no,src):

    Spectral = os.path.join(codedir,"Spectral")
    Instruction = os.path.join(codedir,"Instruction_src"+src+"_ls"+linesearch_no)
    forward = os.path.join(datadir,"forward_src"+src+"_ls"+linesearch_no)

    shutil.copyfile(Spectral,Instruction)

    #~ Check if process has already run. If process has run correctly, there should be a corresponding vz_cc_iter().fits file
    #~ If it doesn't exist then run the linesearch.
    #~ if not os.path.exists(os.path.join(forward,"vz_cc_iter"+str(iter_no).zfill(2)+".fits")):

    #mpipath="/home/jishnu/lib/openmpi/bin/mpiexec"
    # mpipath="/home/apps/openmpi-1.6.5/bin/mpiexec"
    mpipath=os.path.join(HOME,"anaconda2/bin/mpiexec")
    sparccmd=mpipath+" -np 1 ./sparc "+src+" "+linesearch_no

    t0=time.time()
    with open(os.path.join(datadir,forward,"out_linesearch"),'w') as outfile:
        fwd=subprocess.call(sparccmd.split(),stdout=outfile,env=env,cwd=codedir)

    assert fwd==0,"Error in running linesearch for lsno "+str(linesearch_no)
    t1=time.time()

    # shutil.copyfile(os.path.join(forward,"vz_cc.fits"),os.path.join(forward,"vz_cc_iter"+str(iter_no).zfill(2)+".fits"))

    partialfiles=glob.glob(os.path.join(forward,"*partial*"))
    for f in partialfiles: os.remove(f)

    fullfiles=glob.glob(os.path.join(forward,"*full*"))
    for f in fullfiles: os.remove(f)

    return t1-t0

    #~ else: return 0


evaltime=compute_forward(ls_no,src)
evaltime=divmod(evaltime,60)
print "Finished running linesearch no",linesearch_no,"for src no",src_no,"on proc",procno,"in",evaltime[0],"mins",evaltime[1],"secs"
