from __future__ import division
import numpy as np
import plotc
import pyfits
import os,re,sys,glob,fnmatch
from matplotlib.ticker import MaxNLocator,LogLocator,ScalarFormatter
import read_params
import matplotlib.pyplot as plt
import itertools
from matplotlib import rc

def fitsread(fitsfile): return np.squeeze(pyfits.getdata(fitsfile))
rc("font",family="serif")
rc("text",usetex=True)

datadir = read_params.get_directory()
parentdir = os.path.dirname(datadir)
updatedir = os.path.join(datadir,"update")
strategydir=[os.path.join(parentdir,modeldir) for modeldir in ("f_p1","f_to_p3","f_to_p7_2pixsmooth","f_to_p7_new")]

#~ Get ridges used
ridges=read_params.get_modes_used()
modes={'0':'fmode'}
for i in xrange(1,8): modes[str(i)]='p'+str(i)+'mode'
modes['8']='first_bounce_pmode'

mp=np.loadtxt(os.path.join(datadir,'master.pixels'),ndmin=1)
number_of_sources = len(mp)

#~ First panel, compare different strategies
linestyles = itertools.cycle(('solid','dashed','dotted','dashdot'))
strategylabel = iter(["$\#$"+str(i) for i in xrange(1,5)])
itercutoff=35
ax=plt.subplot(131)
for directory in strategydir:
    misfitfiles = fnmatch.filter(os.listdir(os.path.join(directory,"update")),'misfit_[0-9][0-9]')
    misfitfiles=sorted([os.path.join(directory,"update",f) for f in misfitfiles])[:itercutoff+1]
    total_misfit = np.zeros(len(misfitfiles))
    for fileno,misfitfile in enumerate(misfitfiles):
        total_misfit[fileno] = np.sum(np.loadtxt(misfitfile,usecols=[2]))
        
    plt.semilogy(total_misfit,color='black',ls=next(linestyles),linewidth=2,label=next(strategylabel))

ax.xaxis.set_major_locator(MaxNLocator(5))
plt.tick_params(axis='both', which='major', labelsize=16)
plt.xlabel("Iteration",fontsize=18)
plt.ylabel('Total misfit',fontsize=18)
plt.legend()
plt.tick_params(axis="both",labelsize=16)
ax.grid()

#~ Second panel, different ridges for one strategy
misfitfiles=sorted([os.path.join(updatedir,f) for f in fnmatch.filter(os.listdir(updatedir),'misfit_[0-9][0-9]')])
itercutoff = 35
itercutoff = min(len(misfitfiles),itercutoff)
modemisfit = np.zeros((len(ridges),itercutoff+1))
ax=plt.subplot(132)
linestyles = itertools.cycle(('solid','dashed','dotted','dashdot'))
for ridgeno,ridge in enumerate(ridges[:4]):
    
    nsources_found = 0
    
    for src in xrange(1,number_of_sources+1):
        
        for iterno in xrange(itercutoff+1):
        
            ttfile = os.path.join(strategydir[2],'tt','iter'+str(iterno).zfill(2),
                        'ttdiff_src'+str(src).zfill(2)+'.'+modes[ridge])
            try:
                tt=np.loadtxt(ttfile)
                npix = tt.shape[0]
                modemisfit[ridgeno,iterno] += sum((tt[:,1]/60)**2)/npix
                nsources_found+=1
            except IOError:
                pass
            
    modemisfit[ridgeno] /= nsources_found            

    plt.semilogy(modemisfit[ridgeno],color='black',ls=next(linestyles),label=modes[ridge][:-4],linewidth=2)
   
ax.xaxis.set_major_locator(MaxNLocator(5))
plt.tick_params(axis='both', which='major', labelsize=16)    
plt.xlabel("Iteration",fontsize=18)
plt.ylabel('Mean misfit',fontsize=18)  
plt.legend(ncol=2)
plt.tick_params(axis="both",labelsize=16)
ax.grid()


#~ Third panel, waveforms
ax=plt.subplot(133)
datadir = strategydir[2]
src = read_params.parse_cmd_line_params("src",mapto=int,default=1)
srcloc = mp[src-1]
data=fitsread(os.path.join(datadir,'tt','data','data'+str(src).zfill(2)+'.fits'))
nt,nx = data.shape
filt = fitsread('fmode_filter.fits')
data_filtered = np.fft.ifft2(np.fft.fft2(data)*filt).real

vzcc_start = fitsread(os.path.join(datadir,"tt","iter00","vz_cc_src"+str(src).zfill(2)+".fits"))
vzcc_start_filtered = np.fft.ifft2(np.fft.fft2(vzcc_start)*filt).real

iterdirs = sorted(glob.glob(os.path.join(datadir,'tt','iter*')))
itercutoff = 35
iter_dir = iterdirs[min(itercutoff,len(iterdirs)-1)]
iter_to_plot = iter_dir[-2:]
vzcc_iter = fitsread(os.path.join(iter_dir,'vz_cc_src'+str(src).zfill(2)+'.fits'))
vzcc_iter_filtered = np.fft.ifft2(np.fft.fft2(vzcc_iter)*filt).real

#~ Try to load tt file
ttfile = os.path.join(iter_dir,'ttdiff_src'+str(src).zfill(2)+'.fmode')
ttdiff_array = np.loadtxt(ttfile)
pixel_ttdiff = ttdiff_array[:,0]
lef_ttdiff = ttdiff_array[:,2]
rig_ttdiff = ttdiff_array[:,3]

Lx = read_params.get_xlength()
dt = read_params.get_dt()

t = np.arange(nt)*dt/60.

def dist_to_pix(dist): return int((dist+Lx/2+srcloc)/Lx*nx)
def pix_to_dist(pix): return (pix-nx/2)/nx*Lx-srcloc
pix = dist_to_pix(20.)

lef_ttdiff_cutoff = t[map(int,lef_ttdiff[pixel_ttdiff==pix])]
rig_ttdiff_cutoff = t[map(int,rig_ttdiff[pixel_ttdiff==pix])]

plt.plot(t,data_filtered[:,pix],label="True",color='black')
plt.plot(t,vzcc_start_filtered[:,pix],label="Iter 0",linestyle='dashed',linewidth=2,color='#555555')
plt.plot(t,vzcc_iter_filtered[:,pix],label="Iter "+iter_to_plot,linestyle='dashed',linewidth=2,color='#333333',marker='o')
ax.yaxis.set_major_locator(MaxNLocator(4,prune='both'))
yaxis_formatter = ScalarFormatter(useMathText=False)
#~ yaxis_formatter = ScalarFormatter(useMathText=True)
yaxis_formatter.set_scientific(True)
yaxis_formatter.set_powerlimits((0,0))
ax.yaxis.set_major_formatter(yaxis_formatter)
ax.xaxis.set_major_locator(MaxNLocator(5))

plt.ylabel("Wave velocity\n(arbitrary units)",fontsize=18)
plt.tick_params(axis="both",labelsize=16)
    
xlim_left = lef_ttdiff_cutoff
xlim_right = rig_ttdiff_cutoff
plt.xlim(xlim_left,xlim_right)
plt.legend(loc='best')

plt.xlabel("Time (min)",fontsize=18)

plt.gcf().set_size_inches(12,4)
plt.tight_layout()

save = read_params.parse_cmd_line_params("save")
if save is not None:
    savepath = os.path.join("plots",save)
    print "saving to",savepath
    if not os.path.exists("plots"): os.makedirs("plots")
    plt.savefig(savepath)
else:
    print "Not saving plot to file"

#~ plt.show()
