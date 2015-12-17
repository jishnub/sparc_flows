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


typeinp=filter(lambda x: x.startswith("type="),sys.argv)
if len(typeinp)==0: mistype="data"
else: mistype = typeinp[0].split("=")[-1]



srclocs = np.loadtxt(os.path.join(datadir,'master.pixels'),ndmin=1)
nsrc = len(srclocs)

ridges=read_params.get_modes_used()
modes={'0':'fmode'}
for i in xrange(1,10): modes[str(i)]='p'+str(i)+'mode'

if mistype == "data":

    modemisfit = np.zeros((len(ridges),nfiles,nsrc))

    subplot_layout = plotc.layout_subplots(len(ridges))[:2]


    for plotno,ridge in enumerate(ridges):
        
        for src in xrange(1,nsrc+1):
            for iterno in xrange(nfiles):
            
                ttfile = os.path.join(datadir,'tt','iter'+str(iterno).zfill(2),
                            'ttdiff_src'+str(src).zfill(2)+'.'+modes[ridge])
                try:
                    tt=np.loadtxt(ttfile)
                    npix = tt.shape[0]
                    modemisfit[plotno,iterno,src-1] = sum((tt[:,1]/60)**2)/npix
                except IOError:
                    modemisfit[plotno,iterno,src-1] = np.nan
                
        
            subplot_index = (subplot_layout+(plotno+1,))

            plt.subplot(*subplot_index)
            plt.semilogy(range(nfiles),modemisfit[plotno,:,src-1],'o-',label="x="+str(int(srclocs[src-1]))+" Mm")
            plt.semilogy(range(nfiles),[0]*nfiles,ls='dashed')

            plt.title(modes[ridge],fontsize=16,loc='right')
            
            

            #~ plt.gca().yaxis.set_major_locator(MaxNLocator(4,prune='both'))
            plt.tick_params(axis='both', which='major', labelsize=14)
            #~ plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

            plt.xlim(-0.5,nfiles+0.5)
            #~ plt.legend()
            
        #~ insetax=plt.axes([0.35, 0.5, .3, .3])
        #~ plt.plot(range(nfiles),np.sum(modemisfit[plotno],axis=-1))
        #~ plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        #~ plt.ylabel("Total Misfit")

    plt.gcf().text(0.5, 0.02, 'Iteration number', ha='center',fontsize=20)
    plt.gcf().text(0.02, 0.5, 'Mean misfit per receiver ($s^2$)', va='center', rotation='vertical',fontsize=20)
    
    plt.subplots_adjust(hspace=0.3)
    plt.show()

elif mistype == "model":
    try:  
        truemodel=np.squeeze(pyfits.getdata("true_vx.fits"))
    except IOError:
        print "True model doesn't exist"
        quit()
    
    modelmisfit_list = []
    
    for iterno in xrange(nfiles):    
        try:  
            itermodel=np.squeeze(pyfits.getdata(os.path.join(datadir,"update","vx_"+str(iterno).zfill(2)+".fits")))
        except IOError:
            print "vx_"+str(iterno).zfill(2)+".fits doesn't exist"
            quit()
    
        modelmisfit=np.sqrt(np.sum((truemodel-itermodel)**2))

        model0=np.squeeze(pyfits.getdata(os.path.join(datadir,"update","vx_00.fits")))
        model0_misfit=np.sqrt(np.sum((truemodel-model0)**2))
        modelmisfit/=model0_misfit
        
        modelmisfit_list.append(modelmisfit)
    
    plt.plot(range(nfiles),modelmisfit_list,'o-',label="vx")
    plt.xlim(-0.5,nfiles+0.5)
    plt.ylim(0.4,1.05)
    plt.xlabel("Iteration number",fontsize=20)
    plt.ylabel("Normalized model misfit",fontsize=20)
    
    plt.gca().yaxis.set_major_locator(MaxNLocator(5,prune='both'))
    plt.tick_params(axis='both', which='major', labelsize=14)
    

    try:  
        truemodel=np.squeeze(pyfits.getdata("true_vz.fits"))
    except IOError:
        print "True model doesn't exist"
        quit()
    
    modelmisfit_list = []
    
    for iterno in xrange(nfiles):    
        try:  
            itermodel=np.squeeze(pyfits.getdata(os.path.join(datadir,"update","vz_"+str(iterno).zfill(2)+".fits")))
        except IOError:
            print "vz_"+str(iterno).zfill(2)+".fits doesn't exist"
            quit()
    
        modelmisfit=np.sqrt(np.sum((truemodel-itermodel)**2))

        model0=np.squeeze(pyfits.getdata(os.path.join(datadir,"update","vz_00.fits")))
        model0_misfit=np.sqrt(np.sum((truemodel-model0)**2))
        modelmisfit/=model0_misfit
        
        modelmisfit_list.append(modelmisfit)
    
    plt.plot(range(nfiles),modelmisfit_list,'s-',label="vz")
    
    plt.legend(loc="best")
    plt.tight_layout()
    plt.show()

    

