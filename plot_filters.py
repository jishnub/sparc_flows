from __future__ import division
import read_params
import numpy as np
import matplotlib.pyplot as plt
import plotc
import modefilters
import os,sys
import pyfits
import warnings

def fitswrite(filename,array):
    warnings.filterwarnings('ignore')
    pyfits.writeto(filename,array,clobber=True)

codedir=os.path.dirname(os.path.abspath(__file__))
datadir=read_params.get_directory()

try:
    src=next(f for f in sys.argv if (f.startswith("src=") or f.startswith("source=")))
    src=int(src.split("=")[-1])
except StopIteration: src=1

srcloc=np.loadtxt(os.path.join(datadir,'master.pixels'),ndmin=1)[src-1]

data=pyfits.getdata(os.path.join(datadir,'forward_src'+str(src).zfill(2)+'_ls00','data.fits'))

nt,_,nx=data.shape
Lx=read_params.get_xlength()
dt=read_params.get_dt()
x=np.linspace(-Lx/2,Lx/2,nx,endpoint=False)
t=np.arange(nt)*dt

nu=np.fft.fftfreq(nt,dt)*1e3
k=np.fft.fftfreq(nx,Lx/nx)*2*np.pi
Rsun=695.8

modes=dict(fmode=False,
p1mode=False,
p2mode=False,
p3mode=False,
p4mode=False,
p5mode=False,
p6mode=False,
p7mode=False,
high_pmodes=False)

if len(sys.argv)>1: 
    for m in sys.argv[1:]:
        for mode in modes:
            if mode.startswith(m): modes[mode]=True
else:
    print "Usage: python plot_filters.py <list of modes>"
    print "eg: python plot_filters.py fmode p3mode p4mode"

plot_time_distance=True

###########################################################################################

for mode,to_plot in modes.items():
    if to_plot:
        mode_filter_function = getattr(modefilters,mode+"_filter")
        mode_filter = mode_filter_function(nt,dt,nx,Lx)
        mode_filter=mode_filter.transpose(2,1,0)

        #~ Load filter parameters form file simply to plot lower and upper boundary lines
        #~ Actual filter computation doesn't depend on what we load here
        Poly=[]
        Polylow=[]
        f_low=0

        with open('modefilters.f90','r') as mf:
            mode_match=False
            for line in mf.readlines():
                if 'subroutine '+mode+'_filter' in line.lower(): mode_match=True
                if not mode_match: continue
                if 'end subroutine '+mode+'_filter' in line.lower(): 
                    mode_match=False
                    break
                if line.strip().startswith('Polylow(') and mode_match:
                    Polylow.append(float(line.strip().split("=")[1]))
                if line.strip().startswith('Poly(') and mode_match:
                    Poly.append(float(line.strip().split("=")[1]))
                if line.lower().strip().startswith('f_low') and mode_match:
                    f_low=float(line.strip().split("=")[1])
           
        Poly=np.array(Poly)
        Polylow=np.array(Polylow)
        
        plt.figure()
                    
        full_powerspectrum=abs(np.fft.fft2(np.squeeze(data)))**2
        powmax = full_powerspectrum.max()
        
        plot_freq_range=[2,6]
                    
        cp=plotc.spectrumplot(full_powerspectrum,x=k*Rsun,y=nu,vmax=powmax*0.1,
        xr=[0,None],yr=plot_freq_range,axes_properties=dict(xscilimits=(-4,4)),sp=121,
        title="Power spectrum",colorbar=True,
        colorbar_properties={'orientation':'horizontal','shrink':0.8})
        
        ax=cp.axis
        
        plt.xlabel("$k R_\odot$",fontsize=14)
        plt.ylabel("Frequency ($mHz$)",fontsize=14)

        f0=sum(p_i*k**i for i,p_i in enumerate(Polylow))
        f1=sum(p_i*k**i for i,p_i in enumerate(Poly))

        plt.plot(k[:len(k)/2]*Rsun,f0[:len(k)/2],color='orangered')
        plt.plot(k[:len(k)/2]*Rsun,f1[:len(k)/2],color='orangered')
        plt.plot(k[:len(k)/2]*Rsun,[f_low]*(len(k)//2),color='orangered',linestyle='dashed')
                
        filtered_powerspectrum=abs(np.fft.fft2(np.squeeze(data))*np.squeeze(mode_filter))**2
        
        plotc.spectrumplot(filtered_powerspectrum,x=k*Rsun,y=nu,
        xr=[0,None],yr=plot_freq_range,sp=122,
        axes_properties=dict(xscilimits=(-4,4),hide_yticklabels=True),
        title=mode,colorbar=True,
        colorbar_properties={'orientation':'horizontal','shrink':0.8},
        subplot_properties={'sharey':ax,'sharex':ax})
        
        plt.xlabel("$k R_\odot$",fontsize=14)

        plt.subplots_adjust(wspace=0)
        
        if plot_time_distance:
            plt.figure()
            td=np.fft.ifft2(np.fft.fft2(np.squeeze(data))*np.squeeze(mode_filter)).real
            plotc.colorplot(td/td.max(),x=x-srcloc,y=t/60,title=mode+" time distance diagram",
                vmax=0.5,centerzero=True,axes_properties=dict(scilimits=(-5,5)),
                xr=[-Lx/2-srcloc,Lx/2-srcloc],yr=[5,t.max()/60*0.85])
            plt.xlabel("Horizontal Distance ($\mathrm{Mm}$)",fontsize=14)
            plt.ylabel("Time (min)",fontsize=14)


#########################################################################################

plotc.show()
