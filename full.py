import os,shutil,glob,re,subprocess,datetime,time,read_params,sys

env=dict(os.environ, MPI_TYPE_MAX="1280280")

codedir=os.path.dirname(os.path.abspath(__file__))
HOME=os.environ["HOME"]

datadir=read_params.get_directory()

iterno=len([f for f in os.listdir(os.path.join(datadir,'update')) if re.match(r'misfit_[0-9]{2}$',f)])
iterno2dig=str(iterno).zfill(2)

with open(os.path.join(datadir,'master.pixels'),'r') as mp:
    nmasterpixels=sum(1 for _ in mp)

procno=int(os.environ["PBS_VNODENUM"])
nodeno=int(os.environ["PBS_NODENUM"])

if procno>=nmasterpixels: 
    print "Stopping job on node",nodeno,"proc",procno,"at",time.strftime("%H:%M:%S")
    quit()

src=str(procno+1).zfill(2)

def compute_forward_adjoint_kernel(src):

    forward="forward_src"+src+"_ls00"
    adjoint="adjoint_src"+src
    kernel="kernel"

    ttname="vz_cc_src"+src+".fits"
    Spectral=os.path.join(codedir,"Spectral")
    Adjoint=os.path.join(codedir,"Adjoint")
    Instruction=os.path.join(codedir,"Instruction_src"+src+"_ls00")
    
    modes={'0':'fmode','1':'p1mode','2':'p2mode'}
    
    mpipath=os.path.join(HOME,"anaconda/bin/mpiexec")
    sparccmd=mpipath+" -np 1 ./sparc "+src+" 00"
    
    ####################################################################
    #~ Forward
    ####################################################################
    
    if not os.path.exists(os.path.join(datadir,"status",forward)):

        shutil.copyfile(Spectral,Instruction)

        with open(os.path.join(datadir,forward,"out"+forward),'w') as outfile:
            fwd=subprocess.call(sparccmd.split(),stdout=outfile,env=env,cwd=codedir)
        
        if not os.path.exists(os.path.join(datadir,"tt","iter"+iterno2dig)):
            os.makedirs(os.path.join(datadir,"tt","iter"+iterno2dig))
            
        shutil.copyfile(os.path.join(datadir,forward,"vz_cc.fits"),
                        os.path.join(datadir,"tt","iter"+iterno2dig,ttname))
                        
        shutil.copyfile(os.path.join(datadir,forward,"vz_cc.fits"),
                        os.path.join(datadir,forward,"vz_cc_00.fits"))
                        
        shutil.copyfile(os.path.join(codedir,"vz_00.fits"),
                        os.path.join(datadir,"update","vz_"+iterno2dig+".fits"))
                        
        shutil.copyfile(os.path.join(codedir,"vx_00.fits"),
                        os.path.join(datadir,"update","vx_"+iterno2dig+".fits"))
                        
        ridge_filters=read_params.get_ridge_filter()
        
        for ridge_filter in ridge_filters:
            try:
                shutil.copyfile(os.path.join(datadir,forward,"ttdiff."+ridge_filter),
                                os.path.join(datadir,"tt","iter"+iterno2dig,"ttdiff."+modes.get(ridge_filter,ridge_filter)))
            except IOError:
                sys.stderr.write("Could not copy "+os.path.join(datadir,forward,"ttdiff."+ridge_filter))
                sys.stderr.flush()
    
    ####################################################################
    #~ Adjoint
    ####################################################################

    if not os.path.exists(os.path.join(datadir,"status",adjoint)):
        shutil.copyfile(Adjoint,Instruction)

        with open(os.path.join(datadir,adjoint,"out"+adjoint),'w') as outfile:
            adj=subprocess.call(sparccmd.split(),stdout=outfile,env=env,cwd=codedir)
            
    ####################################################################
    #~ Kernel
    ####################################################################

    if not os.path.exists(os.path.join(datadir,"status",kernel+src)):

        with open(os.path.join(datadir,kernel,"out_kernel"+src),'w') as outfile:
            kern=subprocess.call(sparccmd.split(),stdout=outfile,env=env,cwd=codedir)
    
print "Launching on proc no",procno,"for source",src,"at time",datetime.datetime.strftime(datetime.datetime.now(), '%Y-%m-%d %H:%M:%S')

compute_forward_adjoint_kernel(src)
print "Finishing on proc no",procno,"for source",src,"at time",datetime.datetime.strftime(datetime.datetime.now(), '%Y-%m-%d %H:%M:%S')
