import sys,os,glob,re
import numpy as np
import pyfits
import read_params
import warnings

codedir=os.path.dirname(os.path.abspath(__file__))

datadir=read_params.get_directory()

try:
    misfittype=next(element for element in sys.argv if element in ['data','model_psi'])
except StopIteration:
    misfittype="data"
    
try:
    normtype=next(element for element in sys.argv if element in ['--normed','--unnormed'])
except StopIteration:
    normtype="--unnormed"

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

if '--iterno' in sys.argv:
    print "iteration",int(iterno)
warnings.filterwarnings('error')

if misfittype=="data":
    
    misfitfile=os.path.join(datadir,"update","misfit_"+iterno)
    
    if not os.path.exists(misfitfile):
        print misfitfile,"doesn't exist"
        quit()

    try:
        misfitdata=np.loadtxt(misfitfile,usecols=[2],ndmin=1)
    except Warning:
        print "Error reading",misfitfile
    
    misfit_00_found=False
    misfit_00_file=os.path.join(datadir,"update","misfit_00")
    
    if normtype=="--normed":
        try:
            misfit_00=np.loadtxt(misfit_00_file,usecols=[2],ndmin=1)
            misfit_00_found=True
        except IOError:
            print misfit_00_file,"not found, using unnormed"
    
    misfit=sum(misfitdata)
    if normtype=="--normed" and misfit_00_found: misfit/=sum(misfit_00)
    print misfit

#~ elif misfittype=="model_c":
    #~ 
    #~ truemodel=np.squeeze(pyfits.getdata(os.path.join("/scratch",user,"magnetic/true_c_change_B/soundspeed2D.fits")))
    #~ itermodel=np.squeeze(pyfits.getdata(os.path.join("/scratch",user,"magnetic/data/update","model_c_"+iterno+".fits")))
    #~ 
    #~ datamisfit=np.sqrt(np.sum((truemodel-itermodel)**2))
    #~ 
    #~ print datamisfit

elif misfittype=="model_psi":
    
    try:  
        truemodel=np.squeeze(pyfits.getdata("true_psi.fits"))
        truemodel-=truemodel[0,0]
    except IOError:
        print "True model doesn't exist"
        quit()
        
    try:  
        itermodel=np.squeeze(pyfits.getdata(os.path.join(datadir,"update","model_psi_"+iterno+".fits")))
        itermodel-=itermodel[0,0]
    except IOError:
        print "model_psi_"+iterno+".fits doesn't exist"
        quit()
    
    modelmisfit=np.sqrt(np.sum((truemodel-itermodel)**2))
    
    if normtype=="--normed":
        try: 
            model0=np.squeeze(pyfits.getdata(os.path.join(datadir,"update","model_psi_00.fits")))
            model0-=model0[0,0]
            model0_misfit=np.sqrt(np.sum((truemodel-model0)**2))
            modelmisfit/=model0_misfit
        except IOError:
            print "model_psi_00.fits not found, computing unnormed misfits"

    
    print modelmisfit
    
elif misfittype=="model_vx":
    
    try:  
        truemodel=np.squeeze(pyfits.getdata("true_vx.fits"))
    except IOError:
        print "True model doesn't exist"
        quit()
        
    try:  
        itermodel=np.squeeze(pyfits.getdata(os.path.join(datadir,"update","vx_"+iterno+".fits")))
    except IOError:
        print "vx_"+iterno+".fits doesn't exist"
        quit()
    
    modelmisfit=np.sqrt(np.sum((truemodel-itermodel)**2))
    
    if normtype=="--normed":
        try: 
            model0=np.squeeze(pyfits.getdata(os.path.join(datadir,"update","vx_00.fits")))
            model0_misfit=np.sqrt(np.sum((truemodel-model0)**2))
            modelmisfit/=model0_misfit
        except IOError:
            print "vx_00.fits not found, computing unnormed misfits"

    
    print modelmisfit

elif misfittype=="model_vz":
    
    try:  
        truemodel=np.squeeze(pyfits.getdata("true_vz.fits"))
    except IOError:
        print "True model doesn't exist"
        quit()
        
    try:  
        itermodel=np.squeeze(pyfits.getdata(os.path.join(datadir,"update","vz_"+iterno+".fits")))
    except IOError:
        print "vz_"+iterno+".fits doesn't exist"
        quit()
    
    modelmisfit=np.sqrt(np.sum((truemodel-itermodel)**2))
    
    if normtype=="--normed":
        try: 
            model0=np.squeeze(pyfits.getdata(os.path.join(datadir,"update","vz_00.fits")))
            model0_misfit=np.sqrt(np.sum((truemodel-model0)**2))
            modelmisfit/=model0_misfit
        except IOError:
            print "vz_00.fits not found, computing unnormed misfits"

    
    print modelmisfit
