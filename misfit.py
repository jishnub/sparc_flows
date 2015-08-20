import sys,os,glob,re
import numpy as np
import pyfits

codedir=os.path.dirname(os.path.abspath(__file__))

configvars={}
with open(os.path.join(codedir,"varlist.sh")) as myfile:
    for line in myfile:
        name,var=line.partition("=")[::2]
        configvars[name.strip()]=var.strip().strip('"')

datadir=configvars['directory'].replace('$USER',os.environ['USER'])

try:
    misfittype=next(element for element in sys.argv if element in ['data','model_psi'])
except StopIteration:
    misfittype="data"

try:
    iterno=next(element for element in sys.argv if element.isdigit()).zfill(2)
except StopIteration:
    lsfiles=[f for f in glob.glob(os.path.join(datadir,"update","misfit_*")) if "all" not in f]
    nfiles=len(lsfiles)
    if nfiles==0:
        print "No misfit files found"
        quit()
    else:
        iterno=str(nfiles-1).zfill(2)

if misfittype=="data":
    
    misfitfile=os.path.join(datadir,"update","misfit_"+iterno)
    
    if not os.path.exists(misfitfile):
        print misfitfile,"doesn't exist"
        quit()

    with open(os.path.join(datadir,'master.pixels'),'r') as mpixfile:
        nmasterpixels=sum(1 for _ in mpixfile)

    misfitdata=np.loadtxt(misfitfile,usecols=[2])

    misfit=sum(misfitdata[:nmasterpixels])

    print misfit

#~ elif misfittype=="model_c":
    #~ 
    #~ truemodel=np.squeeze(pyfits.getdata(os.path.join("/scratch",user,"magnetic/true_c_change_B/soundspeed2D.fits")))
    #~ itermodel=np.squeeze(pyfits.getdata(os.path.join("/scratch",user,"magnetic/data/update","model_c_"+iterno+".fits")))
    #~ 
    #~ datamisfit=np.sqrt(np.sum((truemodel-itermodel)**2))
    #~ 
    #~ print datamisfit
    #~ 
#~ elif misfittype=="model_psi":
    #~ 
    #~ truemodel=np.squeeze(pyfits.getdata(os.path.join("/scratch",user,"magnetic/true_c_change_B/true_psi.fits")))
    #~ itermodel=np.squeeze(pyfits.getdata(os.path.join("/scratch",user,"magnetic/data/update","model_psi_"+iterno+".fits")))
    #~ 
    #~ datamisfit=np.sqrt(np.sum((truemodel-itermodel)**2))
    #~ 
    #~ print datamisfit
    
