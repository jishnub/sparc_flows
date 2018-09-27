import os,shutil,glob,subprocess,time,read_params,sys,fnmatch
from pathlib import Path
from datetime import datetime

  

def compute_forward_adjoint_kernel(src):

    forward="forward_src{:02d}_ls00".format(src)
    adjoint="adjoint_src{:02d}".format(src)
    kernel="kernel"

    Instruction=codedir/"Instruction_src{:02d}_ls00".format(src)

    ridge_filters = read_params.get_modes_used()
    
    mpipath=HOME/"anaconda3/bin/mpiexec"
    sparccmd="{} -np 1 ./sparc {:02d} 00".format(mpipath,src)
    
    ####################################################################
    #~ Forward
    ####################################################################

    def compute_forward():
    
        shutil.copyfile(codedir/"Spectral",Instruction)

        with open(datadir/forward/"out_forward",'w') as outfile:
            fwd=subprocess.call(sparccmd.split(),stdout=outfile,env=env,cwd=codedir)
        
        assert fwd==0,"Error in running forward"
        
        (datadir/"tt"/"iter{:02d}"/"windows{}".format(src)).mkdir(parents=True,exist_ok=True)

        for ridge_filter in ridge_filters:
            shutil.copyfile(datadir/forward/"ttdiff.{}".format(ridge_filter),
                            datadir/"tt"/"iter{:02d}".format(iterno)/
                            "ttdiff_src{:02d}.{}".format(src,ridge_filter))
    
    ####################################################################
    #~ Adjoint
    ####################################################################

    def compute_adjoint():

        shutil.copyfile(codedir/"Adjoint",Instruction)

        with open(datadir/adjoint/"out_adjoint",'w') as outfile:
            adj=subprocess.call(sparccmd.split(),stdout=outfile,env=env,cwd=codedir)

        assert adj==0,"Error in running adjoint"

    ####################################################################
    #~ Kernel
    ####################################################################

    def compute_kernel():
        with open(datadir/kernel/"out_kernel{:02d}".format(src),'w') as outfile:
            kern=subprocess.call(sparccmd.split(),stdout=outfile,env=env,cwd=codedir)

        assert kern==0,"Error in computing kernel"

    compute_forward()
    compute_adjoint()
    compute_kernel()

if __name__ == "__main__":

    env=dict(os.environ, MPI_TYPE_MAX="1280280")
    env['LD_LIBRARY_PATH'] = ":".join([env.get('LD_LIBRARY_PATH',''),
                                        "/home/apps/gcc-6.1/lib64",
                                        "/home/apps/openmpi-1.6.5/lib",
                                        "/home/apps/lapack-3.5",
                                        "/home/apps/fftw-3.2/lib"])

    codedir=Path(os.path.dirname(os.path.abspath(__file__)))
    HOME=Path(os.environ["HOME"])

    datadir= Path(read_params.get_directory())

    iterno=len(fnmatch.filter(os.listdir(datadir/"update"),'misfit_[0-9][0-9]'))

    with open(datadir/'master.pixels','r') as mp:
        nmasterpixels=sum(1 for _ in mp)

    procno=int(os.environ["PBS_VNODENUM"])

    if procno>=nmasterpixels: quit()

    src = procno + 1

    t_start = datetime.now()

    compute_forward_adjoint_kernel(src)
    delta_t = datetime.now() - t_start

    print(("Finished computing full for source {} in {}".format(src,delta_t)))
