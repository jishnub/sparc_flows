from __future__ import division
import numpy as np
import plotc
import pyfits
import os,re,sys,glob
from matplotlib.ticker import MaxNLocator
import read_params
plt=plotc.plt

datadir = read_params.get_directory()

def fitsread(fitsfile): return np.squeeze(pyfits.getdata(fitsfile))

try:
    src=next(f for f in sys.argv if (f.startswith("src=") or f.startswith("source=")))
    src=int(src.split("=")[-1])
except StopIteration: src=1

vx_true=fitsread('true_vx.fits')
vz_true=fitsread('true_vz.fits')

vx_00 = fitsread(os.path.join(datadir,'update','vx_00.fits'))
vz_00 = fitsread(os.path.join(datadir,'update','vz_00.fits'))

misfit_vx_00  = np.sum(np.sqrt((vx_00-vx_true)**2),axis=1)
misfit_vz_00  = np.sum(np.sqrt((vz_00-vz_true)**2),axis=1)

iters=[0,2,4,6]

Lx=read_params.get_xlength()

Rsun=695.8
z=np.loadtxt(read_params.get_solarmodel(),usecols=[0])
z=(z-1)*Rsun



for iterno in iters:
    vx_iter = fitsread(os.path.join(datadir,'update','vx_'+str(iterno).zfill(2)+'.fits'))
    vz_iter = fitsread(os.path.join(datadir,'update','vz_'+str(iterno).zfill(2)+'.fits'))
    
    misfit_vx  = np.sum(np.sqrt((vx_iter-vx_true)**2),axis=1)
    misfit_vz  = np.sum(np.sqrt((vz_iter-vz_true)**2),axis=1)
    
    ax1=plt.subplot(121)
    plt.plot(z,misfit_vx,label="iter "+str(iterno))
    plt.legend(loc='best')
    
    ax2=plt.subplot(122)
    plt.plot(z,misfit_vz,label="iter "+str(iterno))
    plt.legend(loc='best')

plt.show()
