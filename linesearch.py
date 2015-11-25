import os,sys,shutil,glob,time,re,subprocess,read_params

env=dict(os.environ, MPI_TYPE_MAX="1280280")

codedir=os.path.dirname(os.path.abspath(__file__))
HOME=os.environ["HOME"]

datadir=read_params.get_directory()

procno=int(env["PBS_VNODENUM"])
nodeno=int(env["PBS_NODENUM"])


with open(os.path.join(datadir,'master.pixels'),'r') as mpixfile:
    nmasterpixels=sum(1 for _ in mpixfile)

total_no_of_linesearches=5
total_no_of_jobs=nmasterpixels*total_no_of_linesearches

if procno>=total_no_of_jobs: 
    print "Stopping job on node",nodeno,"proc",procno,"at",time.strftime("%H:%M:%S")
    quit()

linesearch_no=procno/nmasterpixels+1
src_no=procno%nmasterpixels+1

print "Running linesearch no",linesearch_no,"for src no",src_no,"on proc",procno+1,"node no",nodeno+1

ls_no=str(linesearch_no).zfill(2)
src=str(src_no).zfill(2)

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
    t1=time.time()
    
    
    try:
        shutil.copyfile(os.path.join(forward,"vz_cc.fits"),os.path.join(forward,"vz_cc_"+str(linesearch_no)+".fits"))
    except IOError: pass
    
    partialfiles=glob.glob(os.path.join(forward,"*partial*"))
    for f in partialfiles: os.remove(f)
    
    fullfiles=glob.glob(os.path.join(forward,"*full*"))
    for f in fullfiles: os.remove(f)
    
    return t1-t0
    
evaltime=compute_forward(ls_no,src)
print "Finished running linesearch no",linesearch_no,"for src no",src_no,"on proc",procno+1,"node no",nodeno+1,"in",evaltime,"seconds"
