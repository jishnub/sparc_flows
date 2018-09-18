import os,read_params,glob,fnmatch

datadir=read_params.get_directory()
print(datadir)

if not os.path.exists(datadir):
    print(datadir,"doesn't exist. Modify params.i")
    quit()

if not os.path.exists(os.path.join(datadir,"update")):
    print("No iterations done")
    print("qsub data_forward.sh")
    quit()

updatedir = os.path.join(datadir,"update")
lsfiles=sorted([int(f[-2:]) for f in fnmatch.filter(os.listdir(updatedir),'linesearch_[0-9][0-9]')])
misfitfiles=sorted([int(f[-2:]) for f in fnmatch.filter(os.listdir(updatedir),'misfit_[0-9][0-9]')])

if len(misfitfiles)==0:
    print("Zero iterations done.")
    print("qsub full.sh")
    quit()

if len(misfitfiles)>len(lsfiles):
    print("Forward computation for iteration",misfitfiles[-1],"done")
    print("python grad.py algo=cg 0.1 0.2 0.3 0.4 0.5 0.6")
    print("qsub linesearch.sh")
    quit()

if len(misfitfiles)==len(lsfiles):
    print(lsfiles[-1],"iterations complete")
    print("python lsmisfit.py")
    quit()




