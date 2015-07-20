import sys,os,glob,re
import numpy as np

data="/scratch/shivam/flows/data"

if len(sys.argv)>1:
    iterno=sys.argv[1].zfill(2)
else:
    lsfiles=[f for f in glob.glob(os.path.join(data,"update","linesearch_*")) if "all" not in f]
    nfiles=len(lsfiles)
    if nfiles==0:
        print "No linesearch files found"
        quit()
    else:
        iterno=str(nfiles-1).zfill(2)
no_of_linesearches=5


lsfile=os.path.join(data,"update","linesearch_"+iterno)

with open(os.path.join(data,'master.pixels'),'r') as mpixfile:
    nmasterpixels=sum(1 for _ in mpixfile)

lsdata=np.loadtxt(lsfile,usecols=[2])

misfit=[sum(lsdata[i*nmasterpixels:(i+1)*nmasterpixels]) for i in xrange(no_of_linesearches)]

for m in misfit: print m
