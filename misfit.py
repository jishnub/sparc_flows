import sys,os
import numpy as np
import read_params
from astropy.io import fits
from pathlib import Path
import fnmatch

datadir=Path(read_params.get_directory())

try:
    misfittype=next(element for element in sys.argv if element in ['data','psi','vx','vz'])
except StopIteration:
    misfittype="data"
    
try:
    normtype=next(element for element in sys.argv if element in ['--normed','--unnormed'])
except StopIteration:
    normtype="--unnormed"

try:
    iterno=next(element for element in sys.argv if element.isdigit()).zfill(2)
except StopIteration:
    misfit_files=sorted(fnmatch.filter(os.listdir(datadir/"update"),"misfit_[0-9][0-9]"))
    nfiles=len(misfit_files)
    if nfiles==0:
        print("No misfit files found")
        quit()
    else:
        iterno=str(nfiles-1).zfill(2)

def fitsread(f):
    with fits.open(f) as hdul:
        return hdul[0].data.squeeze()


if misfittype=="data":
    
    misfitfile=datadir/"update"/misfit_files[-1]
    
    if not misfitfile.exists():
        print(misfitfile,"doesn't exist")
        quit()

    
    misfitdata=np.loadtxt(misfitfile,usecols=[2],ndmin=1)
    
    misfit_00_found=False
    misfit_00_file=datadir/"update"/"misfit_00"
    
    if normtype=="--normed":
        try:
            misfit_00=np.loadtxt(misfit_00_file,usecols=[2],ndmin=1)
            misfit_00_found=True
        except IOError:
            print(misfit_00_file,"not found, using unnormed")
    
    misfit=sum(misfitdata)
    if normtype=="--normed" and misfit_00_found: misfit/=sum(misfit_00)

    if '--iterno' in sys.argv:
        print("Iteration",int(iterno),"Misfit",'{:.3g}'.format(misfit))
    else:
        print('{:.3g}'.format(misfit))

elif misfittype=="psi":

    try:  
        truemodel=np.squeeze(pyfits.getdata("true_psi.fits"))

    except IOError:
        print("True model doesn't exist")
        quit()
        
    try:  
        itermodel = fitsread(datadir/"update"/"model_psi_"+iterno+".fits")

    except IOError:
        print("model_psi_"+iterno+".fits doesn't exist")
        quit()
    
    modelmisfit=np.sqrt(np.sum((truemodel-itermodel)**2))
    
    if normtype=="--normed":
        try: 
            model0 = fitsread(datadir/"update"/"model_psi_00.fits")

            model0_misfit=np.sqrt(np.sum((truemodel-model0)**2))
            modelmisfit/=model0_misfit
        except IOError:
            print("model_psi_00.fits not found, computing unnormed misfits")

    
    if '--iterno' in sys.argv:
        print(int(iterno),modelmisfit)
    else:
        print(modelmisfit)
    
elif misfittype=="vx":
    
    try:  
        truemodel=fitsread("true_vx.fits")
    except IOError:
        print("True model doesn't exist")
        quit()
        
    try:  
        itermodel=fitsread(datadir/"update"/"vx_"+iterno+".fits")
    except IOError:
        print("vx_"+iterno+".fits doesn't exist")
        quit()
    
    modelmisfit=np.sqrt(np.sum((truemodel-itermodel)**2))
    
    if normtype=="--normed":
        try: 
            model0=fitsread(datadir/"update"/"vx_00.fits")
            model0_misfit=np.sqrt(np.sum((truemodel-model0)**2))
            modelmisfit/=model0_misfit
        except IOError:
            print("vx_00.fits not found, computing unnormed misfits")

    
    if '--iterno' in sys.argv:
        print(int(iterno),modelmisfit)
    else:
        print(modelmisfit)

elif misfittype=="vz":
    
    try:  
        truemodel=fitsread("true_vz.fits")
    except IOError:
        print("True model doesn't exist")
        quit()
        
    try:  
        itermodel=fitsread(datadir/"update"/"vz_"+iterno+".fits")
    except IOError:
        print("vz_"+iterno+".fits doesn't exist")
        quit()
    
    modelmisfit=np.sqrt(np.sum((truemodel-itermodel)**2))
    
    if normtype=="--normed":
        try: 
            model0=fitsread(datadir/"update"/"vz_00.fits")
            model0_misfit=np.sqrt(np.sum((truemodel-model0)**2))
            modelmisfit/=model0_misfit
        except IOError:
            print("vz_00.fits not found, computing unnormed misfits")

    
    if '--iterno' in sys.argv:
        print(int(iterno),modelmisfit)
    else:
        print(modelmisfit)
    
