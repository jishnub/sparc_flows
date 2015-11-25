from __future__ import division
import numpy as np
import plotc
import pyfits
import os
import modefilters,read_params

Rsun=695.8 # Mm

username=os.environ['USER']
homedir=os.environ['HOME']

codedir=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

datadir=read_params.get_directory()

datafile=os.path.join(datadir,"forward_src01_ls00","data.fits")
data=np.squeeze(pyfits.getdata(datafile))


#~ Read params.i

paramsfile = os.path.join(codedir,"params.i")
with open(paramsfile,'r') as paramsfile:
    for line in paramsfile.readlines():
        if ("parameter" in line.lower()) and ("xlength" in line.lower()):
            xlength=float(line.split()[3])
        if ("parameter" in line.lower()) and ("outputcad" in line.lower()):
            dt=float(line.split()[3].split(")")[0])


nt,nx=data.shape

data_fft=np.fft.fft2(data)
data_ampspec=abs(data_fft)

x=np.linspace(-xlength/2,xlength/2,nx,endpoint=False) # Mm
dx=x[1]-x[0]
k=np.fft.fftfreq(nx,dx)*2*np.pi # Mm^-1

nu=np.fft.fftfreq(nt,dt)*1e3 # mHz

fmode=True
p1mode=True
p2mode=True

if fmode:
    # f-mode

    with open('modefilters.data','r') as f:
        filtercutoffs=map(float,f.readlines()[0].split())
    f_mode_const=filtercutoffs[0]
    Poly=filtercutoffs[1:4]
    f_low = filtercutoffs[-2]

    f0=f_mode_const*np.sqrt(abs(k))
    f1=Poly[0] + Poly[1]*k +Poly[2]*k**2.

    plotc.figure()
    plotc.spectrumplot(data_ampspec,x=k*Rsun,y=nu,sp=121,xr=[0,None],yr=[0,6],
                        axes_properties=dict(title="Full spectrum",xlabel=r"$k R_{\odot}$",ylabel=r"$ \nu $ (mHz)"))

    plotc.plt.plot(k[:len(k)//2]*Rsun,f0[:len(k)//2],color="k")
    plotc.plt.plot(k[:len(k)//2]*Rsun,f1[:len(k)//2],color='k')
    plotc.plt.plot(k[:len(k)//2]*Rsun,[f_low]*len(k[:len(k)//2]))
    
    fmode=np.asfortranarray(np.zeros((nx,1,nt),dtype=float))
    modefilters.fmode_filter(fmode)
    fmode=np.squeeze(fmode).T
    data_f_filtered_ampspec=abs(data_fft*fmode)

    plotc.spectrumplot(data_f_filtered_ampspec,x=k*Rsun,y=nu,sp=122,xr=[0,None],yr=[0,6],
                        axes_properties=dict(title="f-mode",xlabel=r"$k R_{\odot}$",ylabel=r"$ \nu $ (mHz)"))
     
if p1mode:
    # p1 mode

    with open('modefilters.data','r') as f:
        filtercutoffs=map(float,f.readlines()[1].split())
    f_mode_const=filtercutoffs[0]
    Poly=filtercutoffs[1:4]
    f_low = filtercutoffs[-2]

    f0=f_mode_const*np.sqrt(abs(k))
    f1=Poly[0] + Poly[1]*k +Poly[2]*k**2.

    plotc.figure()
    plotc.spectrumplot(data_ampspec,x=k*Rsun,y=nu,sp=121,xr=[0,None],yr=[0,6],
                        axes_properties=dict(title="Full spectrum",xlabel=r"$k R_{\odot}$",ylabel=r"$ \nu $ (mHz)"))


    plotc.plt.plot(k[:len(k)//2]*Rsun,f0[:len(k)//2],color="k")
    plotc.plt.plot(k[:len(k)//2]*Rsun,f1[:len(k)//2],color='k')
    plotc.plt.plot(k[:len(k)//2]*Rsun,[f_low]*len(k[:len(k)//2]))
    
    p1mode=np.asfortranarray(np.zeros((nx,1,nt),dtype=float))
    modefilters.pmode_filter(p1mode)
    p1mode=np.squeeze(p1mode).T
    data_p1_filtered_ampspec=abs(data_fft*p1mode)

    plotc.spectrumplot(data_p1_filtered_ampspec,x=k*Rsun,y=nu,sp=122,xr=[0,None],yr=[0,6],
                        axes_properties=dict(title="p1-mode",xlabel=r"$k R_{\odot}$",ylabel=r"$ \nu $ (mHz)"))
                        
    

if p2mode:
    #p2 mode

    with open('modefilters.data','r') as f:
        filtercutoffs=map(float,f.readlines()[2].split())
    f_mode_const=filtercutoffs[0]
    #~ P1=filtercutoffs[0:3]
    #~ P2=filtercutoffs[3:6]
    Poly=filtercutoffs[1:4]
    f_low = filtercutoffs[-2]
    
    #~ f0=P1[0] + P1[1]*k +P1[2]*k**2.
    #~ f1=P2[0] + P2[1]*k +P2[2]*k**2.
    f0=f_mode_const*np.sqrt(abs(k))
    f1=Poly[0] + Poly[1]*k +Poly[2]*k**2.

    plotc.figure()
    plotc.spectrumplot(data_ampspec,x=k*Rsun,y=nu,sp=121,xr=[0,None],yr=[0,6],
                        axes_properties=dict(title="Full spectrum",xlabel=r"$k R_{\odot}$",ylabel=r"$ \nu $ (mHz)"))
    
    plotc.plt.plot(k[:len(k)//2]*Rsun,f0[:len(k)//2],color="k")
    plotc.plt.plot(k[:len(k)//2]*Rsun,f1[:len(k)//2],color='k')
    plotc.plt.plot(k[:len(k)//2]*Rsun,[f_low]*len(k[:len(k)//2]))
    
    p2mode=np.asfortranarray(np.zeros((nx,1,nt),dtype=float))
    modefilters.p2mode_filter(p2mode)
    p2mode=np.squeeze(p2mode).T
    data_p2_filtered_ampspec=abs(data_fft*p2mode)
    
    plotc.spectrumplot(data_p2_filtered_ampspec,x=k*Rsun,y=nu,sp=122,xr=[0,None],yr=[0,6],
                        axes_properties=dict(title="p2-mode",xlabel=r"$k R_{\odot}$",ylabel=r"$ \nu $ (mHz)"))


plotc.show()
