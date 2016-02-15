import os,sys,shutil,glob,time,re,subprocess,read_params,fnmatch

env=dict(os.environ, MPI_TYPE_MAX="1280280")

codedir=os.path.dirname(os.path.abspath(__file__))
HOME=env["HOME"]

datadir=read_params.get_directory()

procno=int(env["PBS_VNODENUM"])
nodeno=int(env["PBS_NODENUM"])


with open(os.path.join(datadir,'master.pixels'),'r') as mpixfile:
    nmasterpixels=sum(1 for _ in mpixfile)

total_no_of_linesearches=int(sys.argv[1])
total_no_of_jobs=nmasterpixels*total_no_of_linesearches

if procno>=total_no_of_jobs: 
    print "Stopping job on node",nodeno,"proc",procno,"at",time.strftime("%H:%M:%S")
    quit()

linesearch_no=procno/nmasterpixels+1
src_no=procno%nmasterpixels+1

print "Running linesearch no",linesearch_no,"for src no",src_no,"on proc",procno+1,"node no",nodeno+1

ls_no=str(linesearch_no).zfill(2)
src=str(src_no).zfill(2)

updatedir=os.path.join(datadir,"update")
iter_no = len(fnmatch.filter(os.listdir(updatedir),'misfit_[0-9][0-9]'))-1

def compute_forward(linesearch_no,src):
    
    Spectral=os.path.join(codedir,"Spectral")
    Instruction=os.path.join(codedir,"Instruction_src"+src+"_ls"+linesearch_no)
    forward = os.path.join(datadir,"forward_src"+src+"_ls"+linesearch_no)
    
    shutil.copyfile(Spectral,Instruction)
    
    mpipath=os.path.join(HOME,"anaconda/bin/mpiexec")
    sparccmd=mpipath+" -np 1 ./sparc "+src+" "+linesearch_no
    
    
    t0=time.time()
    with open(os.path.join(datadir,forward,"out_linesearch_"+linesearch_no),'w') as outfile:
        fwd=subprocess.call(sparccmd.split(),stdout=outfile,env=env,cwd=codedir)
        
    assert fwd==0,"Error in running linesearch for lsno "+str(linesearch_no)
    t1=time.time()
    
    shutil.copyfile(os.path.join(forward,"vz_cc.fits"),os.path.join(forward,"vz_cc_iter"+str(iter_no).zfill(2)+".fits"))
    
    partialfiles=glob.glob(os.path.join(forward,"*partial*"))
    for f in partialfiles: os.remove(f)
    
    fullfiles=glob.glob(os.path.join(forward,"*full*"))
    for f in fullfiles: os.remove(f)
    
    return t1-t0
    
evaltime=compute_forward(ls_no,src)
evaltime=divmod(evaltime,60)
print "Finished running linesearch no",linesearch_no,"for src no",src_no,"on proc",procno,"in",evaltime[0],"mins",evaltime[1],"secs"
