import sys,os,glob,fnmatch
import numpy as np
import read_params

codedir=os.path.dirname(os.path.abspath(__file__))

datadir=read_params.get_directory()

number_args=[x for x in sys.argv if x.isdigit()]
if number_args:
    iterno=number_args[0].zfill(2)
else:
    lsfiles=[os.path.join(datadir,"update",f) for f in 
    fnmatch.filter(os.listdir(os.path.join(datadir,"update")),'linesearch_[0-9][0-9]')]
    nfiles=len(lsfiles)
    if nfiles==0:
        print("No linesearch files found")
        quit()
    else:
        iterno=str(nfiles-1).zfill(2)

num_ls_per_src = len(fnmatch.filter(os.listdir(datadir),"forward_src01_ls[0-9][1-9]"))

lsfile=os.path.join(datadir,"update","linesearch_"+iterno)

if not os.path.exists(lsfile):
    print(lsfile,"doesn't exist")
    quit()

with open(os.path.join(datadir,'master.pixels'),'r') as mpixfile:
    nmasterpixels=sum(1 for _ in mpixfile)

lsdata=np.loadtxt(lsfile,usecols=[2])

misfit=[sum(lsdata[i*nmasterpixels:(i+1)*nmasterpixels]) for i in range(num_ls_per_src)]

print("iteration",int(iterno))

np.set_printoptions(precision=3)
    
for m in misfit: print(m)
