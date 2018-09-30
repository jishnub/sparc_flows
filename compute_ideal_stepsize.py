import numpy as np
import os,fnmatch
import read_params
from pathlib import Path

datadir = Path(read_params.get_directory())

iterno = len(fnmatch.filter(os.listdir(os.path.join(datadir,"update")),"misfit_[0-9][0-9]")) - 1

with np.load(datadir/"epslist.npz") as f:
    step_sizes = f[str(iterno)][0]
    lsdata = f[str(iterno)][1]

lsfile=datadir/"update"/"linesearch_{:02d}".format(iterno)
ls_rm_file = datadir/"update"/"ls_{:02d}.rnm".format(iterno)

nmasterpixels = np.loadtxt(datadir/"master.pixels",ndmin=1).size

no_of_linesearches=len(step_sizes)
if os.path.exists(lsfile):
    lsdata=np.loadtxt(lsfile,usecols=[2])
elif os.path.exists(ls_rm_file):
    lsdata = np.loadtxt(ls_rm_file,usecols=[2])
misfit=[sum(lsdata[i*nmasterpixels:(i+1)*nmasterpixels]) for i in range(no_of_linesearches)]

np.set_printoptions(precision=3)

p=np.polyfit(step_sizes,misfit,2)

min_step = -p[1]/(2*p[0])

arbitrary_steps = [min_step*(0.7+0.1*i) for i in range(no_of_linesearches)]
print("Step sizes: ",list(step_sizes))
print("Linesearch misfits: ",misfit)
print("python grad.py algo=bfgs"+("{:10.2E}"*6).format(*arbitrary_steps))
print("Expected misfits: "+("{:10.2E}"*6).format(*np.polyval(p,arbitrary_steps)))
