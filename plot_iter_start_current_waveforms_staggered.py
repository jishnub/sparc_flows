from __future__ import division
import numpy as np
import plotc
import pyfits
import os,re,sys,glob
from matplotlib.ticker import MaxNLocator
import read_params
plt=plotc.plt

def fitsread(fitsfile): return np.squeeze(pyfits.getdata(fitsfile))

datadir = read_params.get_directory()

try:
    src=next(f for f in sys.argv if (f.startswith("src=") or f.startswith("source=")))
    src=int(src.split("=")[-1])
except StopIteration: src=1

mp=np.loadtxt(os.path.join(datadir,'master.pixels'),ndmin=1)
srcloc = mp[src]

Lx=read_params.get_xlength()
nx = read_params.get_nx()
x= np.linspace(-Lx/2,Lx/2,nx,endpoint=False)

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

data=fitsread(os.path.join(datadir,'tt','data','data'+str(src).zfill(2)+'.fits'))
nt,nx = data.shape
t=np.arange(nt)*read_params.get_dt()/60

vzcc_start=fitsread(os.path.join(datadir,'tt','iter00','vz_cc_src'+str(src).zfill(2)+'.fits'))

itercutoff = filter(lambda x: x.startswith("iter="),sys.argv)
if len(itercutoff)!=0: itercutoff=int(itercutoff[0].split("=")[-1])
else: itercutoff=np.inf

iters_done = glob.glob(os.path.join(datadir,'tt','iter*'))
niters = len(iters_done)-1
iter_to_plot = min(itercutoff,niters)
iter_dir = os.path.join(datadir,'tt','iter'+str(iter_to_plot).zfill(2))
vzcc_iter = fitsread(os.path.join(iter_dir,'vz_cc_src'+str(src).zfill(2)+'.fits'))

modes={'0':'f'}
for i in xrange(1,8): modes[str(i)]='p'+str(i)
modes['8']='large_dist_p'

modes_used=read_params.get_modes_used()

plotridge=filter(lambda x: x in ['f','p1','p2','p3','p4','p5','p6','p7','large_dist_p'],sys.argv)
if plotridge:
    ridges=[r+'mode' for r in sys.argv if r in plotridge]
else:
    ridges=[modes[modes_used[0]]+'mode']

for ridge in ridges:
    plt.figure()
    filt = fitsread(ridge+'_filter.fits')

    data_filtered = np.fft.ifft2(np.fft.fft2(data)*filt).real
    vzcc_start_filtered = np.fft.ifft2(np.fft.fft2(vzcc_start)*filt).real
    vzcc_iter_filtered = np.fft.ifft2(np.fft.fft2(vzcc_iter)*filt).real

    gl=plotc.gridlist(2,1)
    
    plt.subplot(next(gl))
    plt.plot(t,data_filtered[:,pix],label="True model")
    plt.plot(t,vzcc_start_filtered[:,pix],label="Starting model",linestyle='dashed',linewidth=2,color='red')
    plt.gca().yaxis.set_major_locator(MaxNLocator(4,prune='both'))
    plt.tick_params(axis='x', which='major', label1On=False)
    plt.tick_params(axis='y', which='major', labelsize=15)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    ax1=plt.gca()
    plt.ylabel("wave vz",fontsize=25)
    plt.legend(loc='best')

    plt.subplot(next(gl),sharex=ax1,sharey=ax1)
    plt.plot(t,data_filtered[:,pix],label="True model")
    plt.plot(t,vzcc_iter_filtered[:,pix],label="Iterated model",linestyle='dashed',linewidth=2,color='red')
    plt.gca().yaxis.set_major_locator(MaxNLocator(4,prune='both'))
    plt.tick_params(axis='both', which='major', labelsize=15)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    ax2=plt.gca()
    ax2.yaxis.offsetText.set_visible(False)
    ax2.text(-0.21, 1.01, plt.gca().yaxis.get_offset_text().get_text())
    plt.ylabel("wave vz",fontsize=25)
    plt.legend(loc='best')

    plt.xlabel("Time (min)",fontsize=25)
    
plt.tight_layout()
plt.subplots_adjust(hspace=0)
plt.show()

