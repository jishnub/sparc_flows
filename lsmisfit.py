import sys,os,glob,re,fnmatch
import numpy as np
import read_params
from pathlib import Path

codedir=Path(__file__).parent.absolute()

datadir=Path(read_params.get_directory())

iterno = len(fnmatch.filter(os.listdir(datadir/"update"),'misfit_[0-9][0-9]'))-1
iterno = read_params.parse_cmd_line_params("iterno",default=iterno,mapto=int)

if (datadir/"model_c_ls00.fits").exists():
	var = "c"
elif (datadir/"model_psi_ls00.fits").exists():
	var = "psi"

no_of_linesearches=len(fnmatch.filter(os.listdir(datadir/"update"),'test_{}_[0-9].fits'.format(var)))

lsfile=datadir/"update"/"linesearch_{:02d}".format(iterno)

if not lsfile.exists():
    print(lsfile,("doesn't exist, probably linesearch has not been run yet,"
    " or file nas been renamed/removed"))
    quit()

with open(datadir/'master.pixels','r') as mpixfile:
    nmasterpixels=sum(1 for _ in mpixfile)

lsdata=np.loadtxt(lsfile,usecols=[2])

misfit=[sum(lsdata[i*nmasterpixels:(i+1)*nmasterpixels]) for i in range(no_of_linesearches)]

print("iteration",int(iterno))

for m in misfit:
    print('{:.3g}'.format(m))
