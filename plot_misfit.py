from __future__ import division
import numpy as np
import plotc
import pyfits
import os,re,sys,glob
from matplotlib.ticker import MaxNLocator
import read_params
import matplotlib.pyplot as plt

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

itercutoff = filter(lambda x: x.startswith("iter="),sys.argv)
if len(itercutoff)!=0: itercutoff=int(itercutoff[0].split("=")[-1])
else: itercutoff=np.inf

nfiles = min(itercutoff,nfiles)

srclocs = np.loadtxt(os.path.join(datadir,'master.pixels'),ndmin=1)
nsrc = len(srclocs)

ridges=read_params.get_modes_used()
modes={'0':'fmode'}
for i in xrange(1,8): modes[str(i)]='p'+str(i)+'mode'
modes['8']='large_dist_pmode'

def spaced(mode):
    return mode[:-4]+" "+mode[-4:]

if mistype == "data":

    modemisfit = np.zeros((len(ridges),nfiles,nsrc))

    subplot_layout = plotc.layout_subplots(4)[:2]

    for plotno,ridge in enumerate(ridges[:4]):
        
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
            
        plt.title(spaced(modes[ridge].replace("_"," ")),fontsize=16,loc='right')
        plt.tick_params(axis='both', which='major', labelsize=14)    
        plt.xlim(-0.5,nfiles+0.5)
        plt.grid()
    
    for sp in xrange(subplot_layout[1]):
        sp_ind = (subplot_layout[0]-1)*subplot_layout[1]+sp+1
        subplot_index = (subplot_layout+(sp_ind,))
        plt.subplot(*subplot_index)
        plt.xlabel("Iteration",fontsize=25)
        
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.3,left=0.13)
    plt.gcf().text(0.02, 0.5, 'Mean misfit per receiver ($s^2$)', va='center', rotation='vertical',fontsize=25)

elif mistype == "data_summed":

    modemisfit = np.zeros((len(ridges),nfiles))

    subplot_layout = plotc.layout_subplots(4)[:2]

    for ridgeno,ridge in enumerate(ridges):
        
        nsources_found = 0
        
        for src in xrange(1,nsrc+1):
            
            for iterno in xrange(nfiles):
            
                ttfile = os.path.join(datadir,'tt','iter'+str(iterno).zfill(2),
                            'ttdiff_src'+str(src).zfill(2)+'.'+modes[ridge])
                try:
                    tt=np.loadtxt(ttfile)
                    npix = tt.shape[0]
                    modemisfit[ridgeno,iterno] += sum((tt[:,1]/60)**2)/npix
                    nsources_found+=1
                except IOError:
                    pass
                
        modemisfit[ridgeno] /= nsources_found            

        plt.semilogy(range(nfiles),modemisfit[ridgeno],'o-',label=modes[ridge][:-4])
            
    plt.tick_params(axis='both', which='major', labelsize=14)    
    plt.xlim(-0.5,nfiles+0.5)
    plt.grid()
    
    plt.xlabel("Iteration",fontsize=25)
    plt.ylabel('Mean misfit',fontsize=25)  
    plt.legend(loc='best',ncol=2)  
    plt.tight_layout()
    


elif mistype == "model":
    try:  
        truemodel=np.squeeze(pyfits.getdata("true_vx.fits"))
    except IOError:
        print "True model doesn't exist"
        quit()
    
    modelmisfit_list = []
    
    model0=np.squeeze(pyfits.getdata(os.path.join(datadir,"update","vx_00.fits")))
    model0_misfit=np.sqrt(np.sum((truemodel-model0)**2))
    
    for iterno in xrange(nfiles):
        try:  
            itermodel=np.squeeze(pyfits.getdata(os.path.join(datadir,"update","vx_"+str(iterno).zfill(2)+".fits")))
        except IOError:
            print "vx_"+str(iterno).zfill(2)+".fits doesn't exist"
            continue
    
        modelmisfit=np.sqrt(np.sum((truemodel-itermodel)**2))

        
        modelmisfit/=model0_misfit
        
        modelmisfit_list.append(modelmisfit)
    
    plt.plot(range(len(modelmisfit_list)),modelmisfit_list,'o-',label="vx")
    plt.xlim(-0.5,nfiles+0.5)
    plt.ylim(0.5,1.05)
    plt.xlabel("Iteration number",fontsize=25)
    plt.ylabel("Normalized model misfit",fontsize=25)
    
    plt.gca().yaxis.set_major_locator(MaxNLocator(7,prune='both'))
    plt.tick_params(axis='both', which='major', labelsize=18)
    

    try:  
        truemodel=np.squeeze(pyfits.getdata("true_vz.fits"))
    except IOError:
        print "True model doesn't exist"
        quit()
    
    modelmisfit_list = []
    model0=np.squeeze(pyfits.getdata(os.path.join(datadir,"update","vz_00.fits")))
    model0_misfit=np.sqrt(np.sum((truemodel-model0)**2))
    
    for iterno in xrange(nfiles):    
        try:  
            itermodel=np.squeeze(pyfits.getdata(os.path.join(datadir,"update","vz_"+str(iterno).zfill(2)+".fits")))
        except IOError:
            print "vz_"+str(iterno).zfill(2)+".fits doesn't exist"
            continue
    
        modelmisfit=np.sqrt(np.sum((truemodel-itermodel)**2))

        
        modelmisfit/=model0_misfit
        
        modelmisfit_list.append(modelmisfit)
    
    plt.plot(range(len(modelmisfit_list)),modelmisfit_list,'s-',label="vz")
    plt.grid()
    plt.legend(loc="best",fontsize=20)
    plt.tight_layout()
    
    plt.figure()
    
    try:  
        truemodel=np.squeeze(pyfits.getdata("true_psi.fits"))
    except IOError:
        print "True model doesn't exist"
        quit()
    
    modelmisfit_list = []
    
    truemodel -= truemodel[0,0]
    
    model0=np.squeeze(pyfits.getdata(os.path.join(datadir,"update","model_psi_00.fits")))
    model0-=model0[0,0]
    model0_misfit=np.sqrt(np.sum((truemodel-model0)**2))
    
    for iterno in xrange(nfiles):    
        try:  
            itermodel=np.squeeze(pyfits.getdata(os.path.join(datadir,"update","model_psi_"+str(iterno).zfill(2)+".fits")))
            itermodel -= itermodel[0,0]
        except IOError:
            print "model_psi_"+str(iterno).zfill(2)+".fits doesn't exist"
            continue
    
        modelmisfit=np.sqrt(np.sum((truemodel-itermodel)**2))
        modelmisfit/=model0_misfit
        modelmisfit_list.append(modelmisfit)
    
    plt.plot(range(len(modelmisfit_list)),modelmisfit_list,'s-',label=r"$\psi$")
    
    plt.grid()
    plt.legend(loc="best",fontsize=20)
    
    plt.xlim(-0.5,nfiles+0.5)
    plt.ylim(0.5,1.05)
    plt.xlabel("Iteration number",fontsize=25)
    plt.ylabel("Normalized model misfit",fontsize=25)
    
    plt.gca().yaxis.set_major_locator(MaxNLocator(7,prune='both'))
    plt.tick_params(axis='both', which='major', labelsize=18)
    
    plt.tight_layout()
    

plt.show()

    

