from __future__ import division
import matplotlib.pyplot as plt
import pyfits
import read_params
import os,sys,fnmatch
import numpy as np
import plotc


def fitsread(f): return np.squeeze(pyfits.getdata(f))
datadir = read_params.get_directory()

try:
    src=next(f for f in sys.argv if (f.startswith("src=") or f.startswith("source=")))
    src=int(src.split("=")[-1])
except StopIteration: src=1

try:
    ls=next(f for f in sys.argv if (f.startswith("ls=") or f.startswith("linesearch=")))
    ls=int(src.split("=")[-1])
except StopIteration: ls=0

data=fitsread(os.path.join(datadir,'forward_src'+str(src).zfill(2)+'_ls00','data.fits'))
vzcc=fitsread(os.path.join(datadir,'forward_src'+str(src).zfill(2)+'_ls'+str(ls).zfill(2),'vz_cc.fits'))

srcloc=np.loadtxt(os.path.join(datadir,'master.pixels'),ndmin=1)[src-1]

nt=data.shape[0]
dt=read_params.get_dt()/60.
t=np.arange(nt)*dt

Lx=read_params.get_xlength()
nx=read_params.get_nx()
x=np.linspace(-Lx/2,Lx/2,num=nx,endpoint=False)

srcpix = abs(x-srcloc).argmin()

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

modes={'params.0':'fmode'}
for i in xrange(1,6): modes['params.'+str(i)]='p'+str(i)+'mode'
params=fnmatch.filter(os.listdir(datadir),'params.[0-9]')
ridges=[]
for p in params: ridges.append(modes[p])


plotridge=filter(lambda x: x in ['f','p1','p2','p3','p4','p5'],sys.argv)
if plotridge:
    ridges=[r for r in ridges if r[:-4] in plotridge]
else:
    ridges=[ridges[0]]


for ridge in ridges:
    ttfile=np.loadtxt(os.path.join(datadir,'tt','iter00','ttdiff_src'+str(src).zfill(2)+'.'+ridge))
    lef=ttfile[ttfile[:,0]==pix][0,2]
    rig=ttfile[ttfile[:,0]==pix][0,3]
    loc=ttfile[ttfile[:,0]==pix][0,4]
    tmin=ttfile[ttfile[:,0]==pix][0,5]
    tmax=ttfile[ttfile[:,0]==pix][0,6]
    #~ print lef,rig,tmin,tmax
    #~ print t[lef],t[rig],t[tmin],t[tmax]

    modefilter = fitsread(ridge+'_filter.fits')
    data_filtered=np.fft.ifft2(np.fft.fft2(data)*modefilter).real
    vzcc_filtered=np.fft.ifft2(np.fft.fft2(vzcc)*modefilter).real
    
    ax1=plotc.colorplot(data_filtered,cmap='seismic',centerzero=True,
    y=t,x=x-srcloc,colorbar=False,sp=121)[0]
    ax1.set_title(ridge[:-4]+" "+ridge[-4:]+" time-distance",fontsize=20)

    
    plotc.draw_hlines(y=[t[lef],t[rig]],linestyle='dashed')
    plotc.draw_hlines(y=[t[tmin],t[tmax]],linestyle='dashed',color="red")
    plotc.draw_vlines(x=[x[pix]-srcloc],linestyle='dotted')
    plt.ylim(t[10],t[200])
    plt.xlim(-70,70)
    plt.ylabel("Time (min)",fontsize=20)
    plt.xlabel("Horizontal Distance (Mm)",fontsize=20)

    #~ plt.figure()
    plt.subplot(122)
    plt.plot(data_filtered[:,pix],t,color='r',label="observation")
    plt.plot(vzcc_filtered[:,pix],t,color='g',label="model",linestyle='dashed',linewidth=1.5)
    #~ plt.plot(t,data_filtered[:,pix],color='r',label="observation")
    #~ plt.plot(t,vzcc_filtered[:,pix],color='g',label="model",linestyle='dashed',linewidth=1.5)
    plotc.draw_hlines(y=[t[lef],t[rig]],linestyle='dashed')
    plotc.draw_hlines(y=[t[tmin],t[tmax]],linestyle='dashed',color="red")
    #~ plotc.draw_hlines(y=[t[loc]],linestyle='dashed',color="red",linewidth=2)
    plt.ylim(t[lef]*0.5,t[rig]*1.5)
    #~ plt.xlim(t[lef]*0.8,t[rig]*1.2)
    plt.ylabel("Time (min)",fontsize=20)
    #~ plt.xlabel("Time (min)",fontsize=20)
    plt.xlabel("Wave Velocity",fontsize=20)
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    #~ plt.ylabel("Wave Velocity",fontsize=20)
    plt.gca().yaxis.tick_right()
    plt.gca().yaxis.set_label_position("right")

    plt.legend(loc='upper left')
    
    plt.tight_layout()

    plt.show()

