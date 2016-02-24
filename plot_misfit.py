from __future__ import division
import numpy as np
import plotc
import pyfits
import os,re,sys,glob,fnmatch
from matplotlib.ticker import MaxNLocator,LogLocator
import read_params
import matplotlib.pyplot as plt
import itertools
from matplotlib import rc

flows = read_params.if_flows()
sound_speed_perturbed = read_params.if_soundspeed_perturbed()

datadir = read_params.get_directory()
updatedir = os.path.join(datadir,"update")
num_misfit_files=0

#~ Get a list of misfit files in the update directory
misfitfiles=sorted([os.path.join(updatedir,f) for f in fnmatch.filter(os.listdir(updatedir),'misfit_[0-9][0-9]')])
misfit_all_files=sorted([os.path.join(updatedir,f) for f in fnmatch.filter(os.listdir(updatedir),'misfit_all_[0-9][0-9]')])

num_misfit_files=len(misfitfiles)
num_misfit_all_files=len(misfit_all_files)

if num_misfit_files==0:
    print "No misfit files found"
    quit()

#~ If iteration cutoff is specified use it
itercutoff = filter(lambda x: x.startswith("iter="),sys.argv)
if len(itercutoff)!=0: itercutoff=int(itercutoff[0].split("=")[-1])
else: itercutoff=np.inf

num_misfit_files = min(itercutoff,num_misfit_files)
num_misfit_all_files = min(itercutoff,num_misfit_all_files)

#~ What to plot - data or model misfit?
typeinp=filter(lambda x: x.startswith("type="),sys.argv)
if len(typeinp)==0: mistype="data"
else: mistype = typeinp[0].split("=")[-1]

#~ Get source location, useful if plotting sourcewise misfit (data_sourcewise)
srclocs = np.loadtxt(os.path.join(datadir,'master.pixels'),ndmin=1)
num_source = len(srclocs)

#~ Get ridges used
ridges=read_params.get_modes_used()
modes={'0':'fmode'}
for i in xrange(1,8): modes[str(i)]='p'+str(i)+'mode'
modes['8']='first_bounce_pmode'

def spaced(a): 
    b=a[:-4]+" "+a[-4:]
    b=b.replace("_"," ")
    return b

if mistype == "data_sourcewise":

    modemisfit = np.zeros((len(ridges),num_misfit_files,num_source))

    subplot_layout = plotc.layout_subplots(4)[:2]

    for plotno,ridge in enumerate(ridges[:4]):
        
        for src in xrange(1,num_source+1):
            
            for iterno in xrange(num_misfit_files):
            
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
            plt.semilogy(range(num_misfit_files),modemisfit[plotno,:,src-1],'o-',label="x="+str(int(srclocs[src-1]))+" Mm")
            
        plt.title(spaced(modes[ridge].replace("_"," ")),fontsize=16,loc='right')
        plt.tick_params(axis='both', which='major', labelsize=14)    
        plt.xlim(-0.5,num_misfit_files+0.5)
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
    markers = iter(('o', 'v', '8','s','<', 7, '*', 'h', '^', 'D', 'd'))
    linestyles = itertools.cycle(('solid','dashed','dotted'))
    modemisfit = np.zeros((len(ridges),num_misfit_files))

    rc('text', usetex=True)
    rc('font',**{'family':'serif','serif':['Times']})

    plt.subplot(121)

    for ridgeno,ridge in enumerate(ridges[:4]):
        
        nsources_found = 0
        
        for src in xrange(1,num_source+1):
            
            for iterno in xrange(num_misfit_files):
            
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

        plt.semilogy(range(num_misfit_files),modemisfit[ridgeno],marker=next(markers),color='black',
        ls=next(linestyles),label=modes[ridge][:-4])
            
    plt.tick_params(axis='both', which='major', labelsize=14)    
    plt.xlim(-0.5,num_misfit_files+5)
    ax=plt.gca()
    plt.xlabel("Iteration",fontsize=20)
    plt.ylabel('Data misfit',fontsize=20)  
    plt.legend()
    #~ legend=plt.legend(bbox_to_anchor=(1, 0.95),
           #~ bbox_transform=plt.gcf().transFigure)  

    plt.subplot(122)
    for ridgeno,ridge in enumerate(ridges[4:]):
        
        nsources_found = 0
        
        for src in xrange(1,num_source+1):
            
            for iterno in xrange(num_misfit_files):
            
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

        plt.semilogy(range(num_misfit_files),modemisfit[ridgeno],marker=next(markers),color='black',
        ls=next(linestyles),label=modes[ridge][:-4])
            
    plt.tick_params(axis='both', which='major', labelsize=14)    
    plt.xlim(-0.5,num_misfit_files+5)
    ax=plt.gca()
    plt.xlabel("Iteration",fontsize=20)
    plt.ylabel('Data misfit',fontsize=20)  
    plt.legend()
    
    plt.gcf().set_size_inches(8,4)
           
    #~ plotc.apj_2col_format(plt.gcf())
    
    plt.tight_layout()
    if not os.path.exists("plots"): os.makedirs("plots")
    plt.savefig("plots/f6a.eps")
    
elif mistype == "data_freq":
    
    num_lines = np.loadtxt(misfit_all_files[0],comments="#",usecols=map(int,ridges)).shape[0]
    num_freq_bands = num_lines//num_source
    modemisfit = np.zeros((len(ridges),num_misfit_all_files,num_freq_bands))
    
    with open(misfit_all_files[0],'r') as mf:
        all_lines = mf.readlines()
        band_lines = all_lines[:3*num_freq_bands:num_freq_bands]
        for i in xrange(len(band_lines)): band_lines[i]=band_lines[i].strip().split()
        num_rec_line = all_lines[1].strip().split()
    
    for line in band_lines: line.remove("#")    
    num_rec_line.remove("#")

    freqband_list=["{:2.1f} to {:2.1f}".format(float(band_lines[i][1]),float(band_lines[i][2])) for i in xrange(num_freq_bands)]
    
    num_receivers = [int(num_rec_line[i]) for i in map(int,ridges)]
    
    fig,ax = plt.subplots(nrows=3,ncols=3)
    ax= np.array([ax]).flatten()
            
    for iterno,misfit_all_file in enumerate(misfit_all_files):
        misfit_freq=np.loadtxt(misfit_all_file)
        for ridgeno,ridge in enumerate(ridges):
            
            for freqband in xrange(num_freq_bands):
                modemisfit[ridgeno,iterno,freqband] = sum(misfit_freq[freqband::3,int(ridge)])/num_receivers[int(ridge)]


    modemisfit /= num_source
    
    for ridgeno,ridge in enumerate(ridges):
        markers = itertools.cycle(('o', 'v', 's','7','<', '*', 'h', '^', 'D', 'd'))
        linestyles = itertools.cycle(('solid','dashed','dotted'))
        for freqband in xrange(num_freq_bands):
            ax[ridgeno].semilogy(range(num_misfit_all_files),modemisfit[ridgeno,:,freqband],marker=next(markers),color='black',
            ls=next(linestyles),label=freqband_list[freqband])
            #~ ax[ridgeno].legend(loc="best")
            ax[ridgeno].set_xlabel("Iteration")
            ax[ridgeno].set_ylabel("Misfit")
            ax[ridgeno].set_title(spaced(modes[ridge]),loc="right")
    
    handles,labels = ax[0].get_legend_handles_labels()
    ax[-1].legend(handles,labels,loc="center")
    ax[-1].axis("off")
    
    
    plotc.apj_2col_format(plt.gcf())
    plt.tight_layout()
    plt.savefig("plots/f4.eps")


elif mistype == "model":
    
    if flows:
        # vx
        try:  
            truemodel=np.squeeze(pyfits.getdata("true_vx.fits"))
        except IOError:
            print "True vx model doesn't exist"
            quit()
        
        modelmisfit_list = []
        
        model0=np.squeeze(pyfits.getdata(os.path.join(datadir,"update","vx_00.fits")))
        model0_misfit=np.sqrt(np.sum((truemodel-model0)**2))
        
        for iterno in xrange(num_misfit_files):
            try:  
                itermodel=np.squeeze(pyfits.getdata(os.path.join(datadir,"update","vx_"+str(iterno).zfill(2)+".fits")))
            except IOError:
                print "vx_"+str(iterno).zfill(2)+".fits doesn't exist"
                continue
        
            modelmisfit=np.sqrt(np.sum((truemodel-itermodel)**2))

            
            modelmisfit/=model0_misfit
            
            modelmisfit_list.append(modelmisfit)
        
        plt.plot(range(len(modelmisfit_list)),modelmisfit_list,linestyle='solid',marker='^',label="$v_x$",color='black')


        # vz
        try:  
            truemodel=np.squeeze(pyfits.getdata("true_vz.fits"))
        except IOError:
            print "True vz model doesn't exist"
            quit()
        
        modelmisfit_list = []
        model0=np.squeeze(pyfits.getdata(os.path.join(datadir,"update","vz_00.fits")))
        model0_misfit=np.sqrt(np.sum((truemodel-model0)**2))
        
        for iterno in xrange(num_misfit_files):    
            try:  
                itermodel=np.squeeze(pyfits.getdata(os.path.join(datadir,"update","vz_"+str(iterno).zfill(2)+".fits")))
            except IOError:
                print "vz_"+str(iterno).zfill(2)+".fits doesn't exist"
                continue
        
            modelmisfit=np.sqrt(np.sum((truemodel-itermodel)**2))

            
            modelmisfit/=model0_misfit
            
            modelmisfit_list.append(modelmisfit)
        
        plt.plot(range(len(modelmisfit_list)),modelmisfit_list,linestyle='dashed',marker='o',label="$v_z$",color='black')
        
        # Vector potential
        try:  
            truemodel=np.squeeze(pyfits.getdata("true_psi.fits"))
        except IOError:
            print "True psi model doesn't exist"
            quit()
        
        modelmisfit_list = []
        
        truemodel -= truemodel[0,0]
        
        model0=np.squeeze(pyfits.getdata(os.path.join(datadir,"update","model_psi_00.fits")))
        model0-=model0[0,0]
        model0_misfit=np.sqrt(np.sum((truemodel-model0)**2))
        
        for iterno in xrange(num_misfit_files):    
            try:  
                itermodel=np.squeeze(pyfits.getdata(os.path.join(datadir,"update","model_psi_"+str(iterno).zfill(2)+".fits")))
                itermodel -= itermodel[0,0]
            except IOError:
                print "model_psi_"+str(iterno).zfill(2)+".fits doesn't exist"
                continue
        
            modelmisfit=np.sqrt(np.sum((truemodel-itermodel)**2))
            modelmisfit/=model0_misfit
            modelmisfit_list.append(modelmisfit)
        
        plt.plot(range(len(modelmisfit_list)),modelmisfit_list,linestyle='dotted',marker='s',label=r"$\psi$",color='black')
        
        plt.grid()
        plt.legend(loc="best")
        
        plt.xlim(-0.5,num_misfit_files+0.5)
        plt.ylim(0.4,1.1)
        plt.xlabel("Iteration number")
        plt.ylabel("Model misfit")
        
        plt.gca().yaxis.set_major_locator(MaxNLocator(7,prune='both'))

        plotc.apj_1col_format(plt.gcf())
        plt.tight_layout()
        
        if not os.path.exists("plots"): os.makedirs("plots")
        plt.savefig("plots/f6b.eps")
    
    if sound_speed_perturbed:
        try:  
            truemodel=np.squeeze(pyfits.getdata("true_c.fits"))
        except IOError:
            print "True c model doesn't exist"
            quit()
        
        modelmisfit_list = []
        
        model0=np.squeeze(pyfits.getdata(os.path.join(datadir,"update","model_c_00.fits")))
        model0_misfit=np.sqrt(np.sum((truemodel-model0)**2))
        
        for iterno in xrange(num_misfit_files):
            model_c_file = "model_c_"+str(iterno).zfill(2)+".fits"
            try:  
                itermodel=np.squeeze(pyfits.getdata(os.path.join(datadir,"update",model_c_file)))
            except IOError:
                print model_c_file,"doesn't exist"
                continue
        
            modelmisfit=np.sqrt(np.sum((truemodel-itermodel)**2))
            print model_c_file,"misfit",modelmisfit
            
            modelmisfit/=model0_misfit
            
            modelmisfit_list.append(modelmisfit)
        
        plt.plot(range(len(modelmisfit_list)),modelmisfit_list,linestyle='solid',marker='^',label="$c$",color='black')
        
        plt.grid()
        plt.legend(loc="best")
        
        plt.xlim(-0.5,num_misfit_files+0.5)
        plt.ylim(0,1.1)
        plt.xlabel("Iteration number")
        plt.ylabel("Model misfit")
        
        plt.gca().yaxis.set_major_locator(MaxNLocator(7,prune='both'))
        plotc.apj_1col_format(plt.gcf())


plt.show()

    

