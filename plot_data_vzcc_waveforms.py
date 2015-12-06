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
    ls=int(ls.split("=")[-1])
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

if filter(lambda x: x.startswith("dist"),sys.argv): coord='distance'
elif filter(lambda x: x.startswith("pix"),sys.argv): coord='pixel'
elif filter(lambda x: x.startswith("coord"),sys.argv): coord='coord'
else: coord = 'unspecified'
    
pix=130
if coord=='distance':
    dist=float(filter(lambda x: x.startswith("dist"),sys.argv)[0].split("=")[-1])
    xcoord = (srcloc + dist) % Lx 
    if xcoord > Lx/2: xcoord = xcoord - Lx
    pix = abs(x-xcoord).argmin()
        
elif coord == 'coord':
    xcoord=float(filter(lambda x: x.startswith("coord"),sys.argv)[0].split("=")[-1])
    pix = abs(x-xcoord).argmin()

elif coord=='pixel':
    if len(filter(lambda x: x.startswith("pix"),sys.argv)):
        pix=int(filter(lambda x: x.startswith("pix"),sys.argv)[0].split("=")[-1])
        
else:
    print "Using default pixel",pix
    print "Pass pixel as pix=<pix no>, eg pix=123; or distance as dist=<dist>; or coordinate as coord=<coord>"

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
    
modeind = {'fmode':0}
for i in xrange(1,6): modeind['p'+str(i)+'mode']=i

iterno=0
if filter(lambda x: x.startswith("iter=") or x.startswith("iterno="),sys.argv): 
    iterno=int(filter(lambda x: x.startswith("iter=") or x.startswith("iterno="),sys.argv)[0].split("=")[-1])
    assert iterno>=0,"Iteration number must be greater than zero"

wavespeed = np.loadtxt('wavespeed')[src-1]

for ridge in ridges:
    ttfile=np.loadtxt(os.path.join(datadir,'tt','iter'+str(iterno).zfill(2),'ttdiff_src'+str(src).zfill(2)+'.'+ridge))
    lineno = np.where(ttfile[:,0]==pix)[0]
    if len(lineno)==0: 
        pixold=pix
        pix=int(ttfile[:,0][abs(ttfile[:,0]-pix).argmin()])
        print "Pixel",pixold,"not in file, using pixel",pix,"instead"
    lineno = np.where(ttfile[:,0]==pix)[0][0]
    print "pix",ttfile[lineno,0],"lef",ttfile[lineno,2],\
            "rig",ttfile[lineno,3],"loc",ttfile[lineno,4],"tmin",ttfile[lineno,5],"tmax",ttfile[lineno,6]
    lef=ttfile[lineno,2]
    rig=ttfile[lineno,3]

    modefilter = fitsread(ridge+'_filter.fits')
    data_filtered=np.fft.ifft2(np.fft.fft2(data)*modefilter).real
    vzcc_filtered=np.fft.ifft2(np.fft.fft2(vzcc)*modefilter).real
    
    ax1=plotc.colorplot(data_filtered,cmap='seismic',centerzero=True,
    y=t,x=x-srcloc,colorbar=False,sp=121)[0]
    ax1.set_title(ridge[:-4]+" "+ridge[-4:]+" time-distance",fontsize=20)
    
    mode_speed = wavespeed[modeind[ridge]]
    worldline = abs(x-srcloc)/mode_speed+6
    plt.plot(x-srcloc,worldline,linestyle='dashed',color='red')
    
    plotc.draw_hlines(y=[t[lef],t[rig]],linestyle='dashed')
    plotc.draw_vlines(x=[x[pix]-srcloc],linestyle='dotted')
    plt.ylim(t[10],t[200])
    plt.xlim(-70,70)
    plt.ylabel("Time (min)",fontsize=20)
    plt.xlabel("Horizontal Distance (Mm)",fontsize=20)

    plt.subplot(122)
    
    plt.plot(data_filtered[:,pix],t,color='r',label="observation")
    plt.plot(vzcc_filtered[:,pix],t,color='g',label="model",linestyle='dashed',linewidth=1.5)
    ax=plt.gca()
    xlim_cur=ax.get_xlim()
    plotc.draw_hlines(y=[t[lef],t[rig]],xmin=-1,xmax=1,linestyle='dashed')
    ax.set_xlim(xlim_cur)
    plt.ylim(t[lef]*0.8,t[rig]*1.2)
    plt.ylabel("Time (min)",fontsize=20)
    plt.xlabel("Wave Velocity",fontsize=20)
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    ax.yaxis.tick_right()
    ax.yaxis.set_label_position("right")

    plt.legend(loc='upper left')
    
    plt.tight_layout()

    plt.show()

