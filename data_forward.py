import os,sys,shutil,glob,subprocess,read_params
from pathlib import Path
from datetime import datetime


def compute_data(src):

    forward="forward_src{:02d}_ls00".format(src)
    Spectral=codedir/"Spectral"
    Instruction=codedir/"Instruction_src{:02d}_ls00".format(src)
    
    shutil.copyfile(Spectral,Instruction)
    
    mpipath=HOME/"anaconda3/bin/mpiexec"
    sparccmd="{} -np 1 ./sparc {:02d} 00".format(mpipath,src)
   
    with open(os.path.join(datadir,forward,"out_data_forward"),'w') as outfile:
        fwd=subprocess.call(sparccmd.split(),stdout=outfile,env=env,cwd=codedir)

    for f in (datadir/forward).glob("*full*"): os.remove(f)
    
    if (datadir/forward/"vz_cc.fits").exists():
        shutil.copyfile(datadir/forward/"vz_cc.fits",datadir/"data"/"{:02d}.fits".format(src))
        
    if Instruction.exists(): os.remove(Instruction)


if __name__ == "__main__":


    env=dict(os.environ, MPI_TYPE_MAX="1280280")
    env['LD_LIBRARY_PATH'] = ":".join([env.get('LD_LIBRARY_PATH',''),
                                        "/home/apps/gcc-6.1/lib64",
                                        "/home/apps/openmpi-1.6.5/lib",
                                        "/home/apps/lapack-3.5"])

    HOME=Path(os.environ["HOME"])
    codedir=Path(os.path.dirname(__file__))


    datadir=Path(read_params.get_directory())

    procno=int(os.environ["PBS_VNODENUM"])

    with open(os.path.join(datadir,'master.pixels'),'r') as mp:
        nsrc=sum(1 for _ in mp)

    if procno>=nsrc: quit()
    src = procno+1

    time_start = datetime.now()

    compute_data(src)

    delta_t = datetime.now() - time_start

    file_to_remove=datadir/"status"/"forward_src{:02d}_ls00".format(src)
    if file_to_remove.exists(): os.remove(file_to_remove)

    print("Finished computing data for source {} in {}".format(src,delta_t))

    

