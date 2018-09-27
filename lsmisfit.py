import sys,os,glob,re,fnmatch
import numpy as np
import read_params

codedir=os.path.dirname(os.path.abspath(__file__))

datadir=read_params.get_directory()

iterno = len(fnmatch.filter(os.listdir(os.path.join(datadir,"update")),'misfit_[0-9][0-9]'))-1
iterno = read_params.parse_cmd_line_params("iterno",default=iterno,mapto=int)

no_of_linesearches=len(fnmatch.filter(os.listdir(os.path.join(datadir,"update")),'test_psi_[0-9].fits'))

lsfile=os.path.join(datadir,"update","linesearch_{:02d}".format(iterno))

if not os.path.exists(lsfile):
    print(lsfile,("doesn't exist, probably linesearch has not been run yet,"
    " or file nas been renamed/removed"))
    quit()

with open(os.path.join(datadir,'master.pixels'),'r') as mpixfile:
    nmasterpixels=sum(1 for _ in mpixfile)

lsdata=np.loadtxt(lsfile,usecols=[2])

misfit=[sum(lsdata[i*nmasterpixels:(i+1)*nmasterpixels]) for i in range(no_of_linesearches)]

print("iteration",int(iterno))

for m in misfit:
    print(m)
