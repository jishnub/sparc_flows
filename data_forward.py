import os,sys,shutil,glob,subprocess,time,read_params


env=dict(os.environ, MPI_TYPE_MAX="1280280")

codedir=os.path.dirname(os.path.abspath(__file__))
HOME=os.environ["HOME"]

datadir=read_params.get_directory()

procno=int(os.environ["PBS_VNODENUM"])
nodeno=int(os.environ["PBS_NODENUM"])

with open(os.path.join(datadir,'master.pixels'),'r') as mp:
    nmasterpixels=sum(1 for _ in mp)

if procno>=nmasterpixels: 
    print "Stopping job on node",nodeno,"proc",procno,"at",time.strftime("%H:%M:%S")
    quit()

src=str(procno+1).zfill(2)

def compute_data(src):

    forward="forward_src"+src+"_ls00"
    Spectral=os.path.join(codedir,"Spectral")
    Instruction=os.path.join(codedir,"Instruction_src"+src+"_ls00")
    
    shutil.copyfile(Spectral,Instruction)

    mpipath=os.path.join(HOME,"anaconda/bin/mpiexec") 
    sparccmd=mpipath+" -np 1 ./sparc "+src+" 00"
   
    with open(os.path.join(datadir,forward,"out_data_forward"),'w') as outfile:
        fwd=subprocess.call(sparccmd.split(),stdout=outfile,env=env,cwd=codedir)

    partialfiles=glob.glob(os.path.join(datadir,forward,"*partial*"))
    for f in partialfiles: os.remove(f)
    
    fullfiles=glob.glob(os.path.join(datadir,forward,"*full*"))
    for f in fullfiles: os.remove(f)
    
    if not os.path.exists(os.path.join(datadir,"tt","data")):
        os.makedirs(os.path.join(datadir,"tt","data"))
    
    if os.path.exists(os.path.join(datadir,forward,"vz_cc.fits")):
        shutil.move(os.path.join(datadir,forward,"vz_cc.fits"),os.path.join(datadir,forward,"data.fits"))
        shutil.copyfile(os.path.join(datadir,forward,"data.fits"),os.path.join(datadir,"tt","data","data"+src+".fits"))
        shutil.copyfile(os.path.join(datadir,forward,"data.fits"),os.path.join(datadir,"data",src+".fits"))

print "Starting computation on node",nodeno,"proc",procno,"at time",time.strftime("%H:%M:%S")
compute_data(src)

file_to_remove=os.path.join(datadir,"status","forward_src"+src+"_ls00")
if os.path.exists(file_to_remove): os.remove(file_to_remove)

print "Finished computation on node",nodeno,"proc",procno,"at time",time.strftime("%H:%M:%S")

