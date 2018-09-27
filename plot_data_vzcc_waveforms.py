
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
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

srcloc=np.loadtxt(os.path.join(datadir,'master.pixels'),ndmin=1)[src-1]

nt=data.shape[0]
dt=read_params.get_dt()/60.
t=np.arange(nt)*dt

Lx=read_params.get_xlength()
nx=read_params.get_nx()
x=np.linspace(-Lx/2,Lx/2,num=nx,endpoint=False)

srcpix = abs(x-srcloc).argmin()

if [x for x in sys.argv if x.startswith("dist")]: coord='distance'
elif [x for x in sys.argv if x.startswith("pix")]: coord='pixel'
elif [x for x in sys.argv if x.startswith("coord")]: coord='coord'
else: coord = 'unspecified'
    
pix=130
if coord=='distance':
    dist=float([x for x in sys.argv if x.startswith("dist")][0].split("=")[-1])
    xcoord = (srcloc + dist) % Lx 
    if xcoord > Lx/2: xcoord = xcoord - Lx
    pix = abs(x-xcoord).argmin()
        
elif coord == 'coord':
    xcoord=float([x for x in sys.argv if x.startswith("coord")][0].split("=")[-1])
    pix = abs(x-xcoord).argmin()

elif coord=='pixel':
    if len([x for x in sys.argv if x.startswith("pix")]):
        pix=int([x for x in sys.argv if x.startswith("pix")][0].split("=")[-1])
        
else:
    print("Using default pixel",pix)
    print("Pass pixel as pix=<pix no>, eg pix=123; or distance as dist=<dist>; or coordinate as coord=<coord>")

modes={'0':'f'}
for i in range(1,8): modes[str(i)]='p'+str(i)
modes['8']='first_bounce_p'

modes_used=read_params.get_modes_used()

plotridge=[x for x in sys.argv if x in ['f','p1','p2','p3','p4','p5','p6','p7','first_bounce_p']]
if plotridge:
    ridges=[r+'mode' for r in sys.argv if r in plotridge]
else:
    ridges=[modes[modes_used[0]]+'mode']
    
iterno=0
if [x for x in sys.argv if x.startswith("iter=") or x.startswith("iterno=")]: 
    iterno=int([x for x in sys.argv if x.startswith("iter=") or x.startswith("iterno=")][0].split("=")[-1])
    assert iterno>=0,"Iteration number must be greater than zero"
    
ls=0   
ls_passed = [x for x in sys.argv if x.startswith("ls=")]
if ls_passed:
    ls = int(ls_passed[0].split("=")[-1].strip())

if ls==0:
    vzcc=fitsread(os.path.join(datadir,'tt','iter'+str(iterno).zfill(2),'vz_cc_src'+str(src).zfill(2)+'.fits'))
else:
    forwarddir = 'forward_src'+str(src).zfill(2)+'_ls'+str(ls).zfill(2)
    vzcc=fitsread(os.path.join(datadir,forwarddir,
    'vz_cc_iter'+str(iterno).zfill(2)+'.fits'))

for ridge in ridges:
    ttfile = None
    try:
        ttfile=np.loadtxt(os.path.join(datadir,'tt','iter'+str(iterno).zfill(2),'ttdiff_src'+str(src).zfill(2)+'.'+ridge))
        lineno = np.where(ttfile[:,0]==pix)[0]
        if len(lineno)==0: 
            pixold=pix
            pix=int(ttfile[:,0][abs(ttfile[:,0]-pix).argmin()])
            print("Pixel",pixold,"not in file, using pixel",pix,"instead")
        lineno = np.where(ttfile[:,0]==pix)[0][0]
        print("pix",ttfile[lineno,0],"lef",ttfile[lineno,2],\
                "rig",ttfile[lineno,3],"loc",ttfile[lineno,4],"tmin",ttfile[lineno,5],"tmax",ttfile[lineno,6])
        lef=ttfile[lineno,2]
        rig=ttfile[lineno,3]
    except: pass

    modefilter = fitsread(ridge+'_filter.fits')
    data_filtered=np.fft.ifft2(np.fft.fft2(data)*modefilter).real
    vzcc_filtered=np.fft.ifft2(np.fft.fft2(vzcc)*modefilter).real
    
    ax1=plotc.colorplot(data_filtered,cmap='seismic',centerzero=True,
    y=t,x=x-srcloc,colorbar=False,sp=121)[0]
    ax1.set_title(ridge[:-4]+" "+ridge[-4:],fontsize=20)
    
    if ttfile is not None:
        plt.axhline(t[lef],linestyle='dashed')
        plt.axhline(t[rig],linestyle='dashed')
        
    plt.axvline(x[pix]-srcloc,linestyle='dotted')
    plt.ylim(t[10],t[200])
    plt.xlim(-70,70)
    plt.ylabel("Time (min)",fontsize=20)
    plt.xlabel("Horizontal Distance (Mm)",fontsize=20)
    plt.tick_params(axis='both', which='major', labelsize=14)
    plt.gca().xaxis.set_major_locator(MaxNLocator(4,prune='both'))

    plt.subplot(122)
    
    plt.plot(data_filtered[:,pix],t,color='r',label="observation")
    plt.plot(vzcc_filtered[:,pix],t,color='g',label="model",linestyle='dashed',linewidth=1.5)
    ax=plt.gca()
    xlim_cur=ax.get_xlim()
    if ttfile is not None:
        plt.axhline(t[lef],xmin=-1,xmax=1,linestyle='dashed')
        plt.axhline(t[rig],xmin=-1,xmax=1,linestyle='dashed')
        plt.ylim(t[lef]*0.95,t[rig]*1.05)
    ax.set_xlim(xlim_cur)     
    plt.ylabel("Time (min)",fontsize=20)
    plt.xlabel("Wave Velocity",fontsize=20)
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    ax.yaxis.tick_right()
    ax.yaxis.set_label_position("right")
    plt.tick_params(axis='both', which='major', labelsize=14)
    plt.gca().xaxis.set_major_locator(MaxNLocator(4,prune='both'))

    plt.legend(loc='upper left')
    
    plt.tight_layout()

    plt.show()

