import os,shutil,glob,re,subprocess,datetime

env=dict(os.environ, MPI_TYPE_MAX="1280280")

HOME=os.environ['HOME']
codedir=os.path.join(HOME,"sparc")
data="/scratch/shivam/flows/data"

iterno=len([f for f in os.listdir(os.path.join(data,'update')) if re.match(r'misfit_[0-9]{2}$',f)])
itername=str(iterno).zfill(2)

with open(os.path.join(data,'master.pixels'),'r') as mp:
    nmasterpixels=sum(1 for _ in mp)

procno=int(os.environ["PBS_VNODENUM"])

if procno>=nmasterpixels: quit()

src=str(procno+1).zfill(2)

def compute_forward_adjoint_kernel(src):

    forward="forward_src"+src+"_ls00"
    adjoint="adjoint_src"+src
    kernel="kernel"

    ttname="vz_cc_src"+src+".fits"
    tdiff0name="ttdiff_src"+src+".fmode"
    tdiff1name="ttdiff_src"+src+".p1mode"
    tdiff2name="ttdiff_src"+src+".p2mode"
    Spectral=os.path.join(codedir,"Spectral")
    Adjoint=os.path.join(codedir,"Adjoint")
    Instruction=os.path.join(codedir,"Instruction_src"+src+"_ls00")
    
    
    mpipath=os.path.join("/home/apps/openmpi-1.6.5/bin/mpiexec")
    #mpipath=os.path.join("/home/shivam/anaconda/bin/mpiexec")
    sparccmd=mpipath+" -np 1 ./sparc "+src+" 00"
    
    ####################################################################
    #~ Forward
    ####################################################################
    
    if not os.path.exists(os.path.join(data,"status",forward)):

        shutil.copyfile(Spectral,Instruction)

        with open(os.path.join(data,forward,"out"+forward),'w') as outfile:
            fwd=subprocess.call(sparccmd.split(),stdout=outfile,env=env,cwd=codedir)
        
        if not os.path.exists(os.path.join(data,"tt","iter"+itername)):
            os.makedirs(os.path.join(data,"tt","iter"+itername))
            
        shutil.copyfile(os.path.join(data,forward,"vz_cc.fits"),
                        os.path.join(data,"tt","iter"+itername,ttname))
                        
        shutil.copyfile(os.path.join(data,forward,"vz_cc.fits"),
                        os.path.join(data,forward,"vz_cc_00.fits"))
                        
        shutil.copyfile(os.path.join(data,forward,"ttdiff.0"),
                        os.path.join(data,"tt","iter"+itername,tdiff0name))
                        
        shutil.copyfile(os.path.join(data,forward,"ttdiff.1"),
                        os.path.join(data,"tt","iter"+itername,tdiff1name))
                        
        shutil.copyfile(os.path.join(data,forward,"ttdiff.2"),
                        os.path.join(data,"tt","iter"+itername,tdiff2name))
    
    ####################################################################
    #~ Adjoint
    ####################################################################

    if not os.path.exists(os.path.join(data,"status",adjoint)):
        shutil.copyfile(Adjoint,Instruction)

        with open(os.path.join(data,adjoint,"out"+adjoint),'w') as outfile:
            adj=subprocess.call(sparccmd.split(),stdout=outfile,env=env,cwd=codedir)
            
    ####################################################################
    #~ Kernel
    ####################################################################

    if not os.path.exists(os.path.join(data,"status",kernel+src)):

        with open(os.path.join(data,kernel,"out_kernel"+src),'w') as outfile:
            kern=subprocess.call(sparccmd.split(),stdout=outfile,env=env,cwd=codedir)
    
    return 0
    
print "Launching on proc no",procno,"for source",src,"at time",datetime.datetime.strftime(datetime.datetime.now(), '%Y-%m-%d %H:%M:%S')

compute_forward_adjoint_kernel(src)
print "Finishing on proc no",procno,"for source",src,"at time",datetime.datetime.strftime(datetime.datetime.now(), '%Y-%m-%d %H:%M:%S')
