from __future__ import division
import matplotlib.pyplot as plt
import pyfits
import read_params
import os,sys
import numpy as np
import plotc

def fitsread(f): return np.squeeze(pyfits.getdata(f))
datadir = read_params.get_directory()


src=1

data=fitsread(os.path.join(datadir,'forward_src'+str(src).zfill(2)+'_ls00','data.fits'))
vzcc=fitsread(os.path.join(datadir,'forward_src'+str(src).zfill(2)+'_ls00','vz_cc.fits'))

srcloc=np.loadtxt(os.path.join(datadir,'master.pixels'),ndmin=1)[src-1]

nt=data.shape[0]
dt=read_params.get_dt()/60.
t=np.arange(nt)*dt

Lx=read_params.get_xlength()
nx=read_params.get_nx()
x=np.linspace(-Lx/2,Lx/2,num=nx,endpoint=False)

srcpix = abs(x-srcloc).argmin()

print "Source located at",srcloc,"Mm"
print "Source pixel",srcpix

try:
    if filter(lambda x: x.startswith("--dist"),sys.argv): coord='distance'
    elif filter(lambda x: x.startswith("--pix"),sys.argv): coord='pixel'
    if filter(lambda x: x.startswith("--coord"),sys.argv): coord='coord'
except StopIteration:
    coord = 'pixel'
    
if coord=='distance':
    if len(sys.argv)>1:
        dist=float(next(i for i in sys.argv if i.lstrip('-').isdigit()))
        xcoord = (srcloc + dist) % Lx 
        if xcoord > Lx/2: xcoord = xcoord - Lx
        pix = abs(x-xcoord).argmin()
    else:
        pix=130
        
elif coord == 'coord':
    if len(sys.argv)>1:
        xcoord = float(next(i for i in sys.argv if i.lstrip('-').isdigit()))
        pix = abs(x-xcoord).argmin()
    else:
        pix=130
        
elif coord=='pixel':
    if len(sys.argv)>1:
        pix=int(next(i for i in sys.argv if i.isdigit()))
    else:
        pix=130
        print "Using default pixel"


ttfile=np.loadtxt(os.path.join(datadir,'tt','iter00','ttdiff_src'+str(src).zfill(2)+'.p1mode'))
lef=ttfile[ttfile[:,0]==pix][0,2]
rig=ttfile[ttfile[:,0]==pix][0,3]



p1_filter = fitsread('p1mode_filter.fits')
plotc.colorplot(np.fft.ifft2(np.fft.fft2(data)*p1_filter).real,cmap='seismic',centerzero=True,
y=t,x=x-srcloc,title="p1 mode time-distance",colorbar=False)

plotc.draw_hlines(y=[t[lef],t[rig]],linestyle='dashed')
plotc.draw_vlines(x=[x[pix]-srcloc],linestyle='dotted')
plt.ylim(t[10],t[200])
plt.xlim(-70,70)
plt.ylabel("Time (min)",fontsize=16)
plt.xlabel("Horizontal Distance (Mm)",fontsize=16)

plt.figure()
plt.plot(t,data[:,pix],label="observation")
plt.plot(t,vzcc[:,pix],label="model")
plotc.draw_vlines(x=[t[lef],t[rig]],linestyle='dashed')
plt.xlim(t[lef]*0.8,t[rig]*1.2)
plt.xlabel("Time (min)",fontsize=16)
plt.ylabel("Amplitude",fontsize=16)

plt.legend(loc='best')

plt.show()

