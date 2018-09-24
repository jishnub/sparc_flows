import os,read_params,glob,fnmatch

datadir=read_params.get_directory()
print(datadir)

if not os.path.exists(datadir):
    print(datadir,"doesn't exist. Modify params.i")
    quit()

if not os.path.exists(os.path.join(datadir,"data","01.fits")):
    print("No iterations done")
    print("python data_forward.py")
    quit()

updatedir = os.path.join(datadir,"update")
lsfiles=sorted([int(f[-2:]) for f in fnmatch.filter(os.listdir(updatedir),'linesearch_[0-9][0-9]')])
misfitfiles=sorted([int(f[-2:]) for f in fnmatch.filter(os.listdir(updatedir),'misfit_[0-9][0-9]')])

if len(misfitfiles)==0:
    print("Zero iterations done.")
    print("python full.py")
    quit()

if len(misfitfiles)>len(lsfiles):
    print("Forward computation for iteration",misfitfiles[-1],"done")
    num_ls_per_src = len(fnmatch.filter(os.listdir(datadir),"forward_src01_ls[0-9][1-9]"))
    eps_str = " ".join(["{:.2f}".format(0.01*i) for i in range(1,num_ls_per_src+1)])
    print("python grad.py algo=cg {}".format(eps_str))
    print("python linesearch.py")
    quit()

if len(misfitfiles)==len(lsfiles):
    print(lsfiles[-1],"iterations complete")
    print("python lsmisfit.py")
    quit()




