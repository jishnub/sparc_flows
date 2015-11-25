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

src=1
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
high_pmodes=False)

if len(sys.argv)>1: 
    for m in sys.argv[1:]:
        for mode in modes:
            if mode.startswith(m): modes[mode]=True
else:
    print "Usage: python plot_filters.py <list of modes>"
    print "eg: python plot_filters.py fmode p3mode p4mode"

plot_time_distance=True

#########################################################################################


if modes['fmode']:
    #~ f mode

    fmode=modefilters.fmode_filter(nt,dt,nx,Lx)
    fmode=fmode.transpose(2,1,0)
    
    savefilter=False
    if savefilter: fitswrite('fmodefiltertest.fits',fmode)

    #~ f_mode_const=0
    Poly=np.zeros(3)
    Polylow=np.zeros(3)
    f_low=0

    with open('modefilters.f90','r') as mf:
        is_f_mode=False
        for line in mf.readlines():
            if 'subroutine fmode_filter' in line.lower(): is_f_mode=True
            if not is_f_mode: continue
            if 'end subroutine fmode_filter' in line.lower(): 
                is_f_mode=False
                break
            for poly_coeff in xrange(len(Poly)):
                if line.strip().startswith('Poly('+str(poly_coeff)+')') and is_f_mode:
                    Poly[poly_coeff]=float(line.strip().split("=")[1])
            for poly_coeff in xrange(len(Polylow)):
                if line.strip().startswith('Polylow('+str(poly_coeff)+')') and is_f_mode:
                    Polylow[poly_coeff]=float(line.strip().split("=")[1])
            if line.lower().strip().startswith('f_low') and is_f_mode:
                f_low=float(line.strip().split("=")[1])
                
    plt.figure()
    
    full_powerspectrum=abs(np.fft.fft2(np.squeeze(data)))**2
    powmax = full_powerspectrum.max()

    ax1,_=plotc.spectrumplot(full_powerspectrum,x=k*Rsun,y=nu,vmax=powmax*0.5,
    xr=[0,None],yr=[0,8],axes_properties=dict(xscilimits=(-4,4)),sp=121,
    title="Power spectrum",colorbar=True,
    colorbar_properties={'orientation':'horizontal','shrink':0.8})
    
    plt.xlabel("$k R_\odot$",fontsize=14)
    plt.ylabel("Frequency ($mHz$)",fontsize=14)
    
    f0=sum(p_i*k**i for i,p_i in enumerate(Polylow))
    f1=sum(p_i*k**i for i,p_i in enumerate(Poly))
    
    plt.plot(k[:len(k)/2]*Rsun,f0[:len(k)/2],color='orangered')
    plt.plot(k[:len(k)/2]*Rsun,f1[:len(k)/2],color='orangered')
    plt.plot(k[:len(k)/2]*Rsun,[f_low]*(len(k)//2),color='orangered',linestyle='dashed')
    
    filtered_powerspectrum=abs(np.fft.fft2(np.squeeze(data))*np.squeeze(fmode))**2  
    
    if savefilter: fitswrite('fmode_powspec.fits',filtered_powerspectrum)      
    
    ax2,_=plotc.spectrumplot(filtered_powerspectrum,x=k*Rsun,y=nu,
    xr=[0,None],yr=[f_low,8],sp=122,
    axes_properties=dict(xscilimits=(-4,4),hide_yticklabels=True),
    title="f mode",colorbar=True,
    colorbar_properties={'orientation':'horizontal','shrink':0.8},
    subplot_properties={'sharey':ax1,'sharex':ax1})
    
    plt.xlabel("$k R_\odot$",fontsize=14)

    plt.subplots_adjust(wspace=0)
    
    if plot_time_distance:
        plt.figure()
        f_td=np.fft.ifft2(np.fft.fft2(np.squeeze(data))*np.squeeze(fmode)).real
        plotc.colorplot(f_td/f_td.max(),x=x-srcloc,y=t/60,title="f mode time distance diagram",
            vmax=0.5,centerzero=True,axes_properties=dict(scilimits=(-5,5)),
            xr=[-Lx/2-srcloc,Lx/2-srcloc],yr=[5,t.max()/60*0.8])
        plt.xlabel("Horizontal Distance ($\mathrm{Mm}$)",fontsize=14)
        plt.ylabel("Time (min)",fontsize=14)

#############################################################################################

#~ p1 mode

if modes['p1mode']:

    p1mode=modefilters.pmode_filter(nt,dt,nx,Lx)
    p1mode=p1mode.transpose(2,1,0)
    
    savefilter=False
    if savefilter: fitswrite('p1modefiltertest.fits',p1mode)
        

#    f_mode_const=0
    Poly=np.zeros(3)
    Polylow=np.zeros(3)
    f_low=0

    with open('modefilters.f90','r') as mf:
        is_p1_mode=False
        for line in mf.readlines():
            if 'subroutine pmode_filter' in line.lower(): is_p1_mode=True
            if not is_p1_mode: continue
            if 'end subroutine pmode_filter' in line.lower(): 
                is_p1_mode=False
                break
            for poly_coeff in xrange(len(Poly)):
                if line.strip().startswith('Poly('+str(poly_coeff)+')') and is_p1_mode:
                    Poly[poly_coeff]=float(line.strip().split("=")[1])
            for poly_coeff in xrange(len(Polylow)):
                if line.strip().startswith('Polylow('+str(poly_coeff)+')') and is_p1_mode:
                    Polylow[poly_coeff]=float(line.strip().split("=")[1])
            if line.lower().strip().startswith('f_low') and is_p1_mode:
                f_low=float(line.strip().split("=")[1])
           
    
    plt.figure()
                
    full_powerspectrum=abs(np.fft.fft2(np.squeeze(data)))**2
    
    ax1,_=plotc.spectrumplot(full_powerspectrum,x=k*Rsun,y=nu,
    xr=[0,None],yr=[0,8],axes_properties=dict(xscilimits=(-4,4)),sp=121,
    title="Power spectrum",colorbar=True,
    colorbar_properties={'orientation':'horizontal','shrink':0.8})
    
    plt.xlabel("$k R_\odot$",fontsize=14)
    plt.ylabel("Frequency (mHz)",fontsize=14)

    f0=sum(p_i*k**i for i,p_i in enumerate(Polylow))
    f1=sum(p_i*k**i for i,p_i in enumerate(Poly))

    plt.plot(k[:len(k)/2]*Rsun,f0[:len(k)/2],color='orangered')
    plt.plot(k[:len(k)/2]*Rsun,f1[:len(k)/2],color='orangered')
    plt.plot(k[:len(k)/2]*Rsun,[f_low]*(len(k)//2),color='orangered',linestyle='dashed')
            
    filtered_powerspectrum=abs(np.fft.fft2(np.squeeze(data))*np.squeeze(p1mode))**2
    
    if savefilter: fitswrite('p1mode_powspec.fits',filtered_powerspectrum)
    
    ax2,_=plotc.spectrumplot(filtered_powerspectrum,x=k*Rsun,y=nu,
    xr=[0,1450],yr=[f_low,8],sp=122,
    axes_properties=dict(xscilimits=(-4,4),hide_yticklabels=True),
    title="p1 mode",colorbar=True,
    colorbar_properties={'orientation':'horizontal','shrink':0.8},
    subplot_properties={'sharey':ax1,'sharex':ax1})
    
    plt.xlabel("$k R_\odot$",fontsize=14)

    plt.subplots_adjust(wspace=0)
    
    if plot_time_distance:
        plt.figure()
        p1_td=np.fft.ifft2(np.fft.fft2(np.squeeze(data))*np.squeeze(p1mode)).real
        plotc.colorplot(p1_td/p1_td.max(),x=x-srcloc,y=t/60,title="p1 mode time distance diagram",
            vmax=0.5,centerzero=True,axes_properties=dict(scilimits=(-5,5)),
            xr=[-Lx/2-srcloc,Lx/2-srcloc],yr=[5,t.max()/60*0.8],colorbar=False)
        plt.xlabel("Horizontal Distance ($\mathrm{Mm}$)",fontsize=14)
        plt.ylabel("Time (min)",fontsize=14)

#########################################################################################

#~ p2 mode

if modes['p2mode']:

    p2mode=modefilters.p2mode_filter(nt,dt,nx,Lx)
    p2mode=p2mode.transpose(2,1,0)
    
    savefilter=False
    if savefilter: fitswrite('p2modefiltertest.fits',p2mode)

    #~ f_mode_const=0
    Poly=np.zeros(3)
    Polylow=np.zeros(3)
    f_low=0

    with open('modefilters.f90','r') as mf:
        is_p2_mode=False
        for line in mf.readlines():
            if 'subroutine p2mode_filter' in line.lower(): is_p2_mode=True
            if not is_p2_mode: continue
            if 'end subroutine p2mode_filter' in line.lower(): 
                is_p2_mode=False
                break
            for poly_coeff in xrange(len(Poly)):
                if line.strip().startswith('Poly('+str(poly_coeff)+')') and is_p2_mode:
                    Poly[poly_coeff]=float(line.strip().split("=")[1])
            for poly_coeff in xrange(len(Polylow)):
                if line.strip().startswith('Polylow('+str(poly_coeff)+')') and is_p2_mode:
                    Polylow[poly_coeff]=float(line.strip().split("=")[1])
            if line.lower().strip().startswith('f_low') and is_p2_mode:
                f_low=float(line.strip().split("=")[1])
           
    plt.figure()
                
    full_powerspectrum=abs(np.fft.fft2(np.squeeze(data)))**2
    
    ax1,_=plotc.spectrumplot(full_powerspectrum,x=k*Rsun,y=nu,
    xr=[0,None],yr=[0,8],axes_properties=dict(xscilimits=(-4,4)),sp=121,
    title="Power spectrum",colorbar=True,
    colorbar_properties={'orientation':'horizontal','shrink':0.8})
    
    plt.xlabel("$k R_\odot$",fontsize=14)
    plt.ylabel("Frequency ($mHz$)",fontsize=14)

    f0=sum(p_i*k**i for i,p_i in enumerate(Polylow))
    f1=sum(p_i*k**i for i,p_i in enumerate(Poly))

    plt.plot(k[:len(k)/2]*Rsun,f0[:len(k)/2],color='orangered')
    plt.plot(k[:len(k)/2]*Rsun,f1[:len(k)/2],color='orangered')
    plt.plot(k[:len(k)/2]*Rsun,[f_low]*(len(k)//2),color='orangered',linestyle='dashed')
            
    filtered_powerspectrum=abs(np.fft.fft2(np.squeeze(data))*np.squeeze(p2mode))**2
    
    if savefilter: fitswrite('p2mode_powspec.fits',filtered_powerspectrum)
    
    ax2,_=plotc.spectrumplot(filtered_powerspectrum,x=k*Rsun,y=nu,
    xr=[0,None],yr=[0,8],sp=122,
    axes_properties=dict(xscilimits=(-4,4),hide_yticklabels=True),
    title="p2 mode",colorbar=True,
    colorbar_properties={'orientation':'horizontal','shrink':0.8},
    subplot_properties={'sharey':ax1,'sharex':ax1})
    
    plt.xlabel("$k R_\odot$",fontsize=14)

    plt.subplots_adjust(wspace=0)
    
    if plot_time_distance:
        plt.figure()
        p2_td=np.fft.ifft2(np.fft.fft2(np.squeeze(data))*np.squeeze(p2mode)).real
        plotc.colorplot(p2_td/p2_td.max(),x=x-srcloc,y=t/60,title="p2 mode time distance diagram",
            vmax=0.5,centerzero=True,axes_properties=dict(scilimits=(-5,5)),
            xr=[-Lx/2-srcloc,Lx/2-srcloc],yr=[5,t.max()/60*0.8],colorbar=False)
        plt.xlabel("Horizontal Distance ($\mathrm{Mm}$)",fontsize=14)
        plt.ylabel("Time (min)",fontsize=14)

#########################################################################################

#~ p3 mode

if modes['p3mode']:
    p3mode=modefilters.p3mode_filter(nt,dt,nx,Lx)
    p3mode=p3mode.transpose(2,1,0)
    
    savefilter=False
    if savefilter: fitswrite('p3modefiltertest.fits',p3mode)

    f_mode_const=0
    Poly=np.zeros(3)
    Polylow=np.zeros(3)
    f_low=0

    with open('modefilters.f90','r') as mf:
        is_p3_mode=False
        for line in mf.readlines():
            if 'subroutine p3mode_filter' in line.lower(): is_p3_mode=True
            if not is_p3_mode: continue
            if 'end subroutine p3mode_filter' in line.lower(): 
                is_p3_mode=False
                break
            for poly_coeff in xrange(len(Poly)):
                if line.strip().startswith('Poly('+str(poly_coeff)+')') and is_p3_mode:
                    Poly[poly_coeff]=float(line.strip().split("=")[1])
            for poly_coeff in xrange(len(Polylow)):
                if line.strip().startswith('Polylow('+str(poly_coeff)+')') and is_p3_mode:
                    Polylow[poly_coeff]=float(line.strip().split("=")[1])
            if line.lower().strip().startswith('f_low') and is_p3_mode:
                f_low=float(line.strip().split("=")[1])
           
    plt.figure()
    
    full_powerspectrum=abs(np.fft.fft2(np.squeeze(data)))**2
    powmax=full_powerspectrum.max()     
    
    ax1,_=plotc.spectrumplot(full_powerspectrum,x=k*Rsun,y=nu,vmax=powmax*0.5,
    xr=[0,None],yr=[0,8],axes_properties=dict(xscilimits=(-4,4)),sp=121,
    title="Power spectrum",colorbar=True,
    colorbar_properties={'orientation':'horizontal','shrink':0.8})
    
    plt.xlabel("$k R_\odot$",fontsize=14)
    plt.ylabel("Frequency ($mHz$)",fontsize=14)

    f0=sum(p_i*k**i for i,p_i in enumerate(Polylow))
    f1=sum(p_i*k**i for i,p_i in enumerate(Poly))

    plt.plot(k[:len(k)/2]*Rsun,f0[:len(k)/2],color='orangered')
    plt.plot(k[:len(k)/2]*Rsun,f1[:len(k)/2],color='orangered')
    plt.plot(k[:len(k)/2]*Rsun,[f_low]*(len(k)//2),color='orangered',linestyle='dashed')
            
    filtered_powerspectrum=abs(np.fft.fft2(np.squeeze(data))*np.squeeze(p3mode))**2
    
    if savefilter: fitswrite('p3mode_powspec.fits',filtered_powerspectrum)
    
    ax2,_=plotc.spectrumplot(filtered_powerspectrum,x=k*Rsun,y=nu,
    xr=[0,None],yr=[0,8],sp=122,
    axes_properties=dict(xscilimits=(-4,4),hide_yticklabels=True),
    title="p3 mode",colorbar=True,
    colorbar_properties={'orientation':'horizontal','shrink':0.8},
    subplot_properties={'sharex':ax1,'sharey':ax1})
    
    plt.xlabel("$k R_\odot$",fontsize=14)
    
    plt.setp(ax2.get_yticklabels(),visible=False)

    plt.subplots_adjust(wspace=0)
    
    if plot_time_distance:
        plt.figure()
        p3_td=np.fft.ifft2(np.fft.fft2(np.squeeze(data))*np.squeeze(p3mode)).real
        plotc.colorplot(p3_td/p3_td.max(),x=x-srcloc,y=t/60,title="p3 mode time distance diagram",
            vmax=0.5,centerzero=True,axes_properties=dict(scilimits=(-5,5)),
            xr=[-Lx/2-srcloc,Lx/2-srcloc],yr=[5,t.max()/60*0.85])
        plt.xlabel("Horizontal Distance ($\mathrm{Mm}$)",fontsize=14)
        plt.ylabel("Time (min)",fontsize=14)

#########################################################################################

#~ p4 mode

if modes['p4mode']:
    
    p4mode=modefilters.p4mode_filter(nt,dt,nx,Lx)
    p4mode=p4mode.transpose(2,1,0)
    
    savefilter=False
    if savefilter: fitswrite('p4modefiltertest.fits',p4mode)

    #~ f_mode_const=0
    Poly=np.zeros(3)
    Polylow=np.zeros(3)
    f_low=0

    with open('modefilters.f90','r') as mf:
        is_p4_mode=False
        for line in mf.readlines():
            if 'subroutine p4mode_filter' in line.lower(): is_p4_mode=True
            if not is_p4_mode: continue
            if 'end subroutine p4mode_filter' in line.lower(): 
                is_p4_mode=False
                break
            for poly_coeff in xrange(len(Polylow)):
                if line.strip().startswith('Polylow('+str(poly_coeff)+')') and is_p4_mode:
                    Polylow[poly_coeff]=float(line.strip().split("=")[1])
            for poly_coeff in xrange(len(Poly)):
                if line.strip().startswith('Poly('+str(poly_coeff)+')') and is_p4_mode:
                    Poly[poly_coeff]=float(line.strip().split("=")[1])
            if line.lower().strip().startswith('f_low') and is_p4_mode:
                f_low=float(line.strip().split("=")[1])
           
    plt.figure()
                
    full_powerspectrum=abs(np.fft.fft2(np.squeeze(data)))**2
    powmax = full_powerspectrum.max()
    
    plot_freq_range=[2,6]
    
    ax1,_=plotc.spectrumplot(full_powerspectrum,x=k*Rsun,y=nu,vmax=powmax*0.2,
    xr=[0,None],yr=plot_freq_range,axes_properties=dict(xscilimits=(-4,4)),sp=121,
    title="Power spectrum",colorbar=True,
    colorbar_properties={'orientation':'horizontal','shrink':0.8})
    
    plt.xlabel("$k R_\odot$",fontsize=14)
    plt.ylabel("Frequency ($mHz$)",fontsize=14)

    f0=sum(p_i*k**i for i,p_i in enumerate(Polylow))
    f1=sum(p_i*k**i for i,p_i in enumerate(Poly))

    plt.plot(k[:len(k)/2]*Rsun,f0[:len(k)/2],color='orangered')
    plt.plot(k[:len(k)/2]*Rsun,f1[:len(k)/2],color='orangered')
    plt.plot(k[:len(k)/2]*Rsun,[f_low]*(len(k)//2),color='orangered',linestyle='dashed')
            
    filtered_powerspectrum=abs(np.fft.fft2(np.squeeze(data))*np.squeeze(p4mode))**2
    
    
    if savefilter: fitswrite('p4mode_powspec.fits',filtered_powerspectrum)
            
    plotc.spectrumplot(filtered_powerspectrum,x=k*Rsun,y=nu,
    xr=[0,None],yr=plot_freq_range,sp=122,
    axes_properties=dict(xscilimits=(-4,4),hide_yticklabels=True),
    title="p4 mode",colorbar=True,
    colorbar_properties={'orientation':'horizontal','shrink':0.8},
    subplot_properties={'sharey':ax1,'sharex':ax1})
    
    plt.xlabel("$k R_\odot$",fontsize=14)

    plt.subplots_adjust(wspace=0)
    
    if plot_time_distance:
        plt.figure()
        p4_td=np.fft.ifft2(np.fft.fft2(np.squeeze(data))*np.squeeze(p4mode)).real
        plotc.colorplot(p4_td/p4_td.max(),x=x-srcloc,y=t/60,title="p4 mode time distance diagram",
            vmax=0.5,centerzero=True,axes_properties=dict(scilimits=(-5,5)),
            xr=[-Lx/2-srcloc,Lx/2-srcloc],yr=[5,t.max()/60*0.8])
        plt.xlabel("Horizontal Distance ($\mathrm{Mm}$)",fontsize=14)
        plt.ylabel("Time (min)",fontsize=14)

#########################################################################################

#~ p5 mode

if modes['p5mode']:
    
    p5mode=modefilters.p5mode_filter(nt,dt,nx,Lx)
    p5mode=p5mode.transpose(2,1,0)
    
    savefilter=False
    if savefilter: fitswrite('p5modefiltertest.fits',p5mode)

    f_mode_const=0
    Poly=np.zeros(3)
    f_low=0

    with open('modefilters.f90','r') as mf:
        is_p5_mode=False
        for line in mf.readlines():
            if 'subroutine p5mode_filter' in line.lower(): is_p5_mode=True
            if not is_p5_mode: continue
            if 'end subroutine p5mode_filter' in line.lower(): 
                is_p5_mode=False
                break
            if line.lower().strip().startswith('f_mode_const') and is_p5_mode:
                f_mode_const=float(line.strip().split("=")[1])
            for poly_coeff in xrange(len(Poly)):
                if line.strip().startswith('Poly('+str(poly_coeff)+')') and is_p5_mode:
                    Poly[poly_coeff]=float(line.strip().split("=")[1])
            if line.lower().strip().startswith('f_low') and is_p5_mode:
                f_low=float(line.strip().split("=")[1])
           
    plt.figure()
                
    full_powerspectrum=abs(np.fft.fft2(np.squeeze(data)))**2
    powmax = full_powerspectrum.max()
    
    plot_freq_range=[2,6]
                
    ax1,_=plotc.spectrumplot(full_powerspectrum,x=k*Rsun,y=nu,vmax=powmax*0.1,
    xr=[0,None],yr=plot_freq_range,axes_properties=dict(xscilimits=(-4,4)),sp=121,
    title="Power spectrum",colorbar=True,
    colorbar_properties={'orientation':'horizontal','shrink':0.8})
    
    plt.xlabel("$k R_\odot$",fontsize=14)
    plt.ylabel("Frequency ($mHz$)",fontsize=14)

    f0=f_mode_const*abs(k)**0.5
    f1=Poly[0] + Poly[1]*abs(k) +Poly[2]*k**2.

    plt.plot(k[:len(k)/2]*Rsun,f0[:len(k)/2],color='orangered')
    plt.plot(k[:len(k)/2]*Rsun,f1[:len(k)/2],color='orangered')
    plt.plot(k[:len(k)/2]*Rsun,[f_low]*(len(k)//2),color='orangered',linestyle='dashed')
            
    filtered_powerspectrum=abs(np.fft.fft2(np.squeeze(data))*np.squeeze(p5mode))**2
    
    
    if savefilter: fitswrite('p5mode_powspec.fits',filtered_powerspectrum)
    
    plotc.spectrumplot(filtered_powerspectrum,x=k*Rsun,y=nu,
    xr=[0,None],yr=plot_freq_range,sp=122,
    axes_properties=dict(xscilimits=(-4,4),hide_yticklabels=True),
    title="p5 mode",colorbar=True,
    colorbar_properties={'orientation':'horizontal','shrink':0.8},
    subplot_properties={'sharey':ax1,'sharex':ax1})
    
    plt.xlabel("$k R_\odot$",fontsize=14)

    plt.subplots_adjust(wspace=0)
    
    if plot_time_distance:
        plt.figure()
        p5_td=np.fft.ifft2(np.fft.fft2(np.squeeze(data))*np.squeeze(p5mode)).real
        plotc.colorplot(p5_td/p5_td.max(),x=x-srcloc,y=t/60,title="p5 mode time distance diagram",
            vmax=0.5,centerzero=True,axes_properties=dict(scilimits=(-5,5)),
            xr=[-Lx/2-srcloc,Lx/2-srcloc],yr=[5,t.max()/60*0.85])
        plt.xlabel("Horizontal Distance ($\mathrm{Mm}$)",fontsize=14)
        plt.ylabel("Time (min)",fontsize=14)

#########################################################################################

#~ high p modes

if modes['high_pmodes']:

    high_pmodes=modefilters.highpmode_filter(nt,dt,nx,Lx)
    high_pmodes=high_pmodes.transpose(2,1,0)
    
    savefilter=True
    if savefilter: fitswrite('highpmodefiltertest.fits',high_pmodes)

    f_mode_const=0
    Poly=np.zeros(3)
    f_low=0

    with open('modefilters.f90','r') as mf:
        is_hp_mode=False
        for line in mf.readlines():
            if 'subroutine p5mode_filter' in line.lower(): is_hp_mode=True
            if not is_hp_mode: continue
            if 'end subroutine p5mode_filter' in line.lower(): 
                is_hp_mode=False
                break
            if line.lower().strip().startswith('f_mode_const') and is_hp_mode:
                f_mode_const=float(line.strip().split("=")[1])
            for poly_coeff in xrange(len(Poly)):
                if line.strip().startswith('Poly('+str(poly_coeff)+')') and is_hp_mode:
                    Poly[poly_coeff]=float(line.strip().split("=")[1])
            if line.lower().strip().startswith('f_low') and is_hp_mode:
                f_low=float(line.strip().split("=")[1])
           
    plt.figure()
                
    ax1,_=plotc.spectrumplot(abs(np.fft.fft2(np.squeeze(data)))**2,x=k*Rsun,y=nu,
    xr=[0,None],yr=[0,8],axes_properties=dict(xscilimits=(-4,4)),sp=121,
    title="Power spectrum",colorbar=True,
    colorbar_properties={'orientation':'horizontal','shrink':0.8})
    
    plt.xlabel("$k R_\odot$",fontsize=14)
    plt.ylabel("Frequency ($mHz$)",fontsize=14)

    f0=f_mode_const*abs(k)**0.5
    f1=Poly[0] + Poly[1]*abs(k) +Poly[2]*k**2.

    plt.plot(k[:len(k)/2]*Rsun,f0[:len(k)/2],color='orangered')
    plt.plot(k[:len(k)/2]*Rsun,f1[:len(k)/2],color='orangered')
    plt.plot(k[:len(k)/2]*Rsun,[f_low]*(len(k)//2),color='orangered',linestyle='dashed')
    
    filtered_powerspectrum=abs(np.fft.fft2(np.squeeze(data))*np.squeeze(high_pmodes))**2
    if savefilter: fitswrite('p5mode_powspec.fits',filtered_powerspectrum)
            
    ax1,_=plotc.spectrumplot(filtered_powerspectrum,x=k*Rsun,y=nu,
    xr=[0,None],yr=[0,8],sp=122,
    axes_properties=dict(xscilimits=(-4,4),hide_yticklabels=True),
    title="high p modes",colorbar=True,
    colorbar_properties={'orientation':'horizontal','shrink':0.8})
    
    plt.xlabel("$k R_\odot$",fontsize=14)

    plt.subplots_adjust(wspace=0)
    
    if plot_time_distance:
        plt.figure()
        f_td=np.fft.ifft2(np.fft.fft2(np.squeeze(data))*np.squeeze(high_pmodes)).real
        plotc.colorplot(f_td/f_td.max(),x=x,y=t/60,title="f mode time distance diagram",
            vmax=0.5,centerzero=True,axes_properties=dict(scilimits=(-5,5)),
            xr=[-Lx/2-srcloc,Lx/2-srcloc],yr=[5,t.max()/60*0.85])
        plt.xlabel("Horizontal Distance ($\mathrm{Mm}$)",fontsize=14)
        plt.ylabel("Time (min)",fontsize=14)

#########################################################################################

plotc.show()
