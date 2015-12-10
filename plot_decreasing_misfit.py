from __future__ import division
import numpy as np
import plotc
import pyfits
import os,re,sys,glob
from matplotlib.ticker import MaxNLocator
import read_params
plt=plotc.plt

datadir = read_params.get_directory()
nfiles=0

lsfiles=[f for f in glob.glob(os.path.join(datadir,"update","misfit_*")) if "all" not in f]
nfiles=len(lsfiles)
if nfiles==0:
    print "No misfit files found"
    quit()

#~ try:
    #~ srcinp=next(f for f in sys.argv if (f.startswith("src=") or f.startswith("source=")))
    #~ plotsrc=int(srcinp.split("=")[-1])
#~ except StopIteration: plotsrc=1

srclocs = np.loadtxt(os.path.join(datadir,'master.pixels'),ndmin=1)
nsrc = len(srclocs)

ridges=read_params.get_modes_used()
modes={'0':'fmode'}
for i in xrange(1,6): modes[str(i)]='p'+str(i)+'mode'

modemisfit = np.zeros((len(ridges),nfiles,nsrc))

subplot_layout = plotc.layout_subplots(len(ridges))[:2]


for plotno,ridge in enumerate(ridges):
    
    for src in xrange(1,nsrc+1):
        for iterno in xrange(nfiles):
        
            ttfile = os.path.join(datadir,'tt','iter'+str(iterno).zfill(2),
                        'ttdiff_src'+str(src).zfill(2)+'.'+modes[ridge])
            tt=np.loadtxt(ttfile)
            npix = tt.shape[0]
            modemisfit[plotno,iterno,src-1] = sum((tt[:,1]/60)**2)/npix
            
    
        subplot_index = (subplot_layout+(plotno+1,))

        plt.subplot(*subplot_index)
        plt.plot(range(nfiles),modemisfit[plotno,:,src-1],'o-',label="x="+str(int(srclocs[src-1]))+" Mm")
        plt.plot(range(nfiles),[0]*nfiles,ls='dashed')

        plt.xlabel("Iteration number",fontsize=20)
        plt.ylabel("Mean misfit per receiver ($s^2$)",fontsize=20)

        plt.gca().yaxis.set_major_locator(MaxNLocator(4,prune='both'))
        plt.tick_params(axis='both', which='major', labelsize=14)
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

        plt.xlim(-0.5,nfiles+0.5)
        plt.legend()
        
    insetax=plt.axes([0.35, 0.5, .3, .3])
    plt.plot(range(nfiles),np.sum(modemisfit[plotno],axis=-1))
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.ylabel("Total Misfit")

plt.show()


