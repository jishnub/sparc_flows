import os,sys,shutil,glob,subprocess,time,read_params


env=dict(os.environ, MPI_TYPE_MAX="1280280")

codedir=os.path.dirname(os.path.abspath(__file__))
HOME=os.environ["HOME"]

datadir=read_params.get_directory()

procno=int(os.environ["PBS_VNODENUM"])
nodeno=int(os.environ["PBS_NODENUM"])

try:
    with open(os.path.join(datadir,'master.pixels'),'r') as mp:
        nmasterpixels=sum(1 for _ in mp)
except IOError:
    sys.stderr.write("Proc"+str(procno).zfill(2)+
                ": Could not read master.pixels, check if it exists\n") 
    sys.stderr.flush()
    quit()

if procno>=nmasterpixels: 
    print "Stopping job on node",nodeno,"proc",procno,"at",time.strftime("%H:%M:%S")
    quit()

src=str(procno+1).zfill(2)

def safecopy(a,b):
    try: shutil.copyfile(a,b)
    except IOError as e:
        sys.stderr.write("Could not copy "+a+" to "+b+"; "+e.args(1)+"\n")
        sys.stderr.flush()

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
        try: os.makedirs(os.path.join(datadir,"tt","data"))
        except OSError as e:
            if e.errno == 17: pass
            else: print e
    
    if os.path.exists(os.path.join(datadir,forward,"vz_cc.fits")):
        shutil.move(os.path.join(datadir,forward,"vz_cc.fits"),os.path.join(datadir,forward,"data.fits"))
        safecopy(os.path.join(datadir,forward,"data.fits"),os.path.join(datadir,"tt","data","data"+src+".fits"))
        safecopy(os.path.join(datadir,forward,"data.fits"),os.path.join(datadir,"data",src+".fits"))
        
    if os.path.exists(Instruction): os.remove(Instruction)

print "Starting computation on node",nodeno,"proc",procno,"at time",time.strftime("%H:%M:%S")
compute_data(src)

file_to_remove=os.path.join(datadir,"status","forward_src"+src+"_ls00")
if os.path.exists(file_to_remove): os.remove(file_to_remove)

print "Finished computation on node",nodeno,"proc",procno,"at time",time.strftime("%H:%M:%S")

