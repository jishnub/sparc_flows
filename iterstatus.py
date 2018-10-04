import os,read_params,fnmatch,filecmp,numpy as np
from pathlib import Path

datadir=Path(read_params.get_directory())
print(datadir)

if not datadir.exists():
    print(datadir,"doesn't exist. Modify params.i")
    quit()

if not (datadir/"master.pixels").exists():
    print(datadir/"master.pixels","doesn't exist. Please generate this before running the code")
    quit()

with open(datadir/'master.pixels','r') as mp:
        nsrc=sum(1 for _ in mp)

if not (datadir/"data"/"01.fits").exists():
    print("No iterations done")
    print("qsub data_forward.sh")
    quit()

num_ls_per_src = len(fnmatch.filter(os.listdir(datadir),"forward_src01_ls[0-9][1-9]"))
updatedir = datadir/"update"
lsfiles=sorted(fnmatch.filter(os.listdir(updatedir),'linesearch_[0-9][0-9]'))
misfitfiles=sorted(fnmatch.filter(os.listdir(updatedir),'misfit_[0-9][0-9]'))

if len(misfitfiles)==0:
    print("Zero iterations done.")
    print("qsub full.sh")
    quit()

if len(misfitfiles)>len(lsfiles):
    print("Forward computation for iteration",int(misfitfiles[-1][-2:]),"done")
    eps_str = " ".join(["{:.2f}".format(0.01*i) for i in range(1,num_ls_per_src+1)])
    print("python grad.py algo=bfgs {}".format(eps_str))
    print("qsub linesearch.sh")
    quit()

if len(misfitfiles)==len(lsfiles):

    # Read misfit.py and check if model has been updated
    lsdata = np.loadtxt(updatedir/lsfiles[-1],usecols=2)
    misfit=np.array([sum(lsdata[i*nsrc:(i+1)*nsrc]) for i in range(num_ls_per_src)])
    best_model = misfit.argmin()+1
    if filecmp.cmp(datadir/"model_psi_ls00.fits",updatedir/"test_psi_{}.fits".format(best_model)):
        print("iteration",int(lsfiles[-1][-2:])+1,"complete")
        print("qsub full.sh")
    else:
        print("python lsmisfit.py")
        print("If you detect a minimum in misfit, run")
        print("python copy_updated_model.py")
    quit()




