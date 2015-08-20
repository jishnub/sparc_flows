import sys,os,glob,re
import numpy as np

codedir=os.path.dirname(os.path.abspath(__file__))

configvars={}
with open(os.path.join(codedir,"varlist.sh")) as myfile:
    for line in myfile:
        name,var=line.partition("=")[::2]
        configvars[name.strip()]=var.strip().strip('"')

datadir=configvars['directory'].replace('$USER',os.environ['USER'])

if len(sys.argv)>1:
    iterno=next(element for element in sys.argv if element.isdigit()).zfill(2)
else:
    lsfiles=[f for f in glob.glob(os.path.join(datadir,"update","linesearch_*")) if "all" not in f]
    nfiles=len(lsfiles)
    if nfiles==0:
        print "No linesearch files found"
        quit()
    else:
        iterno=str(nfiles-1).zfill(2)
no_of_linesearches=5


lsfile=os.path.join(datadir,"update","linesearch_"+iterno)

if not os.path.exists(lsfile):
    print lsfile,"doesn't exist"
    quit()

with open(os.path.join(datadir,'master.pixels'),'r') as mpixfile:
    nmasterpixels=sum(1 for _ in mpixfile)

lsdata=np.loadtxt(lsfile,usecols=[2])

misfit=[sum(lsdata[i*nmasterpixels:(i+1)*nmasterpixels]) for i in xrange(no_of_linesearches)]

for m in misfit: print m
