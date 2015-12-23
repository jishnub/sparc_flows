import os,shutil,glob,subprocess,datetime,time,read_params,sys,fnmatch

env=dict(os.environ, MPI_TYPE_MAX="1280280")

codedir=os.path.dirname(os.path.abspath(__file__))
HOME=os.environ["HOME"]

datadir=read_params.get_directory()

def get_iter_no():
    updatedir=os.path.join(datadir,"update")
    # Count the number of misfit_xx files
    return len(fnmatch.filter(os.listdir(updatedir),'misfit_[0-9][0-9]'))
iterno=get_iter_no()
iterno2dig=str(iterno).zfill(2)

with open(os.path.join(datadir,'master.pixels'),'r') as mp:
    nmasterpixels=sum(1 for _ in mp)

procno=int(os.environ["PBS_VNODENUM"])
nodeno=int(os.environ["PBS_NODENUM"])

if procno>=nmasterpixels: 
    print "Stopping job on node",nodeno,"proc",procno,"at",time.strftime("%H:%M:%S")
    quit()

src=str(procno+1).zfill(2)

def safecopy(a,b):
    try: shutil.copyfile(a,b)
    except IOError as e:
        sys.stderr.write("Could not copy "+a+" to "+b+"; "+e.args(1)+"\n")
        sys.stderr.flush()
        
def safemkdir(a):
    if not os.path.exists(a):
        try: os.makedirs(a)
        except OSError:
            if e.errno == 17: pass
            else: print e

def compute_forward_adjoint_kernel(src):

    forward="forward_src"+src+"_ls00"
    adjoint="adjoint_src"+src
    kernel="kernel"

    Spectral=os.path.join(codedir,"Spectral")
    Adjoint=os.path.join(codedir,"Adjoint")
    Instruction=os.path.join(codedir,"Instruction_src"+src+"_ls00")
    
    modes={'0':'fmode'}
    for pmodeno in xrange(1,8): modes.update({str(pmodeno):'p'+str(pmodeno)+'mode'})
    modes['8']='large_dist_pmode'
    
    ridge_filters = read_params.get_modes_used()
    
    for ridge_filter in ridge_filters:
        if os.path.exists(os.path.join(datadir,forward,"ttdiff."+ridge_filter)):
            safecopy(os.path.join(datadir,forward,"ttdiff."+ridge_filter),
                    os.path.join(datadir,forward,"ttdiff_prev."+ridge_filter))
    
    mpipath=os.path.join(HOME,"anaconda/bin/mpiexec")
    sparccmd=mpipath+" -np 1 ./sparc "+src+" 00"
    
    ####################################################################
    #~ Forward
    ####################################################################
    
    safecopy(Spectral,Instruction)

    with open(os.path.join(datadir,forward,"out"+forward),'w') as outfile:
        fwd=subprocess.call(sparccmd.split(),stdout=outfile,env=env,cwd=codedir)
    
    safemkdir(os.path.join(datadir,"tt","iter"+iterno2dig))
    safemkdir(os.path.join(datadir,"tt","iter"+iterno2dig,"windows"+src))

    safecopy(os.path.join(datadir,forward,"vz_cc.fits"),
                    os.path.join(datadir,"tt","iter"+iterno2dig,"vz_cc_src"+src+".fits"))                    

    for ridge_filter in ridge_filters:
        safecopy(os.path.join(datadir,forward,"ttdiff."+ridge_filter),
                        os.path.join(datadir,"tt","iter"+iterno2dig,
                        "ttdiff_src"+src+"."+modes.get(ridge_filter,ridge_filter)))
                        
    windows_files = glob.glob(os.path.join(datadir,forward,"windows.*"))
    for w in windows_files:
        safecopy(w,os.path.join(datadir,"tt","iter"+iterno2dig,"windows"+src,os.path.basename(w)))

    
    ####################################################################
    #~ Adjoint
    ####################################################################

    safecopy(Adjoint,Instruction)

    with open(os.path.join(datadir,adjoint,"out"+adjoint),'w') as outfile:
        adj=subprocess.call(sparccmd.split(),stdout=outfile,env=env,cwd=codedir)

            
    ####################################################################
    #~ Kernel
    ####################################################################

    with open(os.path.join(datadir,kernel,"out_kernel"+src),'w') as outfile:
        kern=subprocess.call(sparccmd.split(),stdout=outfile,env=env,cwd=codedir)

timestart=datetime.datetime.now()
print "Launching on proc no",procno,"for source",src,"at time",datetime.datetime.strftime(timestart, '%Y-%m-%d %H:%M:%S')
compute_forward_adjoint_kernel(src)
timefin= datetime.datetime.now()
#~ elapsedTime = timefin - timestart
#~ runtime=divmod(elapsedTime.total_seconds(), 60)
print "Finished on proc no",procno,"for source",src,"at",datetime.datetime.strftime(timefin, '%Y-%m-%d %H:%M:%S')
