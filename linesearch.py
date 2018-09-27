import os,sys,shutil,glob,subprocess,read_params,fnmatch
from pathlib import Path
from datetime import datetime



def compute_forward(src,linesearch_no):

    forward="forward_src{:02d}_ls{:02d}".format(src,linesearch_no)

    Instruction=codedir/"Instruction_src{:02d}_ls{:02d}".format(src,linesearch_no)

    shutil.copyfile(codedir/"Spectral",Instruction)

    mpipath=HOME/"anaconda3/bin/mpiexec"
    sparccmd="{} -np 1 ./sparc {:02d} {:02d}".format(mpipath,src,linesearch_no)

    with open(datadir/forward/"out_forward",'w') as outfile:
        fwd=subprocess.call(sparccmd.split(),stdout=outfile,env=env,cwd=codedir)

    assert fwd==0,"Error in running linesearch for lsno {:02d}".format(linesearch_no)
    
    for f in (datadir/forward).glob("*full*"): os.remove(f)

if __name__ == "__main__":

    env=dict(os.environ, MPI_TYPE_MAX="1280280")
    env['LD_LIBRARY_PATH'] = ":".join([env.get('LD_LIBRARY_PATH',''),
                                        "/home/apps/gcc-6.1/lib64",
                                        "/home/apps/openmpi-1.6.5/lib",
                                        "/home/apps/lapack-3.5",
                                        "/home/apps/fftw-3.2/lib"])


    codedir=Path(os.path.dirname(__file__))
    HOME=Path(os.environ["HOME"])

    datadir= Path(read_params.get_directory())

    procno=int(env["PBS_VNODENUM"])

    #~ Get the total number of source pixels to determine number of linesearches necessary
    with open(datadir/'master.pixels','r') as mp:
        nsrc=sum(1 for _ in mp)

    num_ls_per_src = len(fnmatch.filter(os.listdir(datadir),"forward_src01_ls[0-9][1-9]"))

    total_no_of_jobs=nsrc*num_ls_per_src

    if procno>=total_no_of_jobs: quit()

    linesearch_no,src_no = divmod(procno,nsrc)
    linesearch_no += 1
    src_no += 1

    t_start= datetime.now()
    compute_forward(src_no,linesearch_no)
    delta_t = datetime.now() - t_start

    print(("Finished computing full for source {} in {}".format(src,delta_t)))