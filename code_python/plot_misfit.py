from __future__ import division
import numpy as np
import sys,os,glob,re,fnmatch
import pyfits
import plotc
import warnings

#######################################################################

def fitsread(f): 
    arr=pyfits.getdata(f)
    # If it is 2D, make it 3D. This adds an extra dimension at the end.
    # Bring it to the middle to match the general trend
    # Dimension will be (nz,ny,nx) after this
    if len(arr.shape)==2: arr=np.atleast_3d(arr).transpose(0,2,1)
    # Change dimension to (nx,ny,nz)
    return arr.transpose(2,1,0)

def twodigit(m): return str(m).zfill(2)

def data_misfit(iterno):
    data="/scratch/shivam/flows/data"
    no_of_linesearches=5


    lsfile=os.path.join(data,"update","misfit_"+iterno)

    with open(os.path.join(data,'master.pixels'),'r') as mpixfile:
        nmasterpixels=sum(1 for _ in mpixfile)
    
    lsdata=np.loadtxt(lsfile,usecols=[2],ndmin=1)

    misfit=sum(lsdata[:nmasterpixels])
    return misfit
########################################################################



Rsun=695.9895 # Mm
    
datadir="/scratch/shivam/flows/data/update"

itermax=26

#~ Get shape from a pre-existing file
psi=fitsread(os.path.join(datadir,'model_psi_00.fits'))
nx,ny,nz=psi.shape

array_shape=(nx,ny,nz)

true_psi=np.zeros(array_shape)
model_psi=np.zeros(array_shape)
misfit=np.zeros(array_shape)
tot_misfit=np.zeros(itermax)
dat_misfit=np.zeros(itermax)

true_psi=fitsread(os.path.join(datadir,'true_psi.fits'))

for src in xrange(0,itermax):
    #print src
    psi=fitsread(os.path.join(datadir,'model_psi_'+twodigit(src)+'.fits'))
    model_psi=psi-7e-04
    misfit=(model_psi-true_psi)**2
    #print np.sum(misfit)
    tot_misfit[src]=np.sum(misfit)

mmax=tot_misfit[0]
tot_misfit=tot_misfit/mmax   
x=np.linspace(0,itermax-1,itermax)

plotc.plt.plot(x[15:itermax-1],tot_misfit[15:itermax-1],'o-',color='red')
plotc.plt.plot(x[0:16],tot_misfit[0:16],'o-',color='blue')
plotc.plt.ylim((0.0,1.0))
plotc.plt.ylabel('Normalised Model Misfit Value')
plotc.plt.xlabel('Iteration Number')
plotc.plt.show()

for src in xrange(0,itermax):
    iterno=twodigit(src)
    dat_misfit[src]= data_misfit(iterno)
y=np.linspace(0,itermax-1,itermax)
plotc.plt.plot(y[15:itermax-1],dat_misfit[15:itermax-1],'o-',color='red')
plotc.plt.plot(y[5:16],dat_misfit[5:16],'o-',color='blue')
plotc.plt.ylim((0.0,2.0))    
plotc.plt.ylabel('Data Misfit Value')
plotc.plt.xlabel('Iteration Number')
plotc.plt.show()
