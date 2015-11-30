from __future__ import division
import plotc
import matplotlib.pyplot as plt
import numpy as np
import read_params
import pyfits
import os,fnmatch

def fitsread(f): return np.squeeze(pyfits.getdata(f))
    
datadir=read_params.get_directory()
codedir=os.path.dirname(os.path.abspath(__file__))

def get_iter_no():
    updatedir=os.path.join(datadir,"update")
    # Count the number of misfit_xx files
    return len(fnmatch.filter(os.listdir(updatedir),'misfit_[0-9][0-9]'))
iterno=get_iter_no()
iterm1=str(iterno-1).zfill(2)

enf_cont = read_params.get_enforced_continuity()
contvar=read_params.get_continuity_variable()

Lx=read_params.get_xlength()
nx=read_params.get_nx()
x=np.linspace(-Lx/2,Lx/2,num=nx,endpoint=False)

Rsun=695.8
z=np.loadtxt(os.path.join(codedir,read_params.get_solarmodel()),usecols=[0])
z=(z-1)*Rsun

if enf_cont and (contvar == 'psi'):
    true_model = fitsread(os.path.join(codedir,'true_psi.fits'))
    current_model = fitsread(os.path.join(datadir,'model_psi_ls00.fits'))
    try:
        update=fitsread(os.path.join(datadir,'update','update_psi_'+iterm1+'.fits'))
    except IOError:
        update=np.zeros_like(current_model)
        print "Could not load update psi"
    _,_,gl=plotc.layout_subplots(3)
    ax1,_=plotc.colorplot(true_model,sp=next(gl),x=x,y=z,title="True psi",
    yr=[-5,None],colorbar_properties={'orientation':'horizontal','shrink':0.8})
    plt.ylabel("Depth (Mm)",fontsize=16)
    plt.xlabel("Horizontal Distance (Mm)",fontsize=16,labelpad=10)
    ax2,_=plotc.colorplot(current_model,sp=next(gl),x=x,y=z,title="Current psi",
    yr=[-5,None],colorbar_properties={'orientation':'horizontal','shrink':0.8},
    axes_properties={'sharey':ax1,'hide_yticklabels':True})
    plt.xlabel("Horizontal Distance (Mm)",fontsize=16,labelpad=10)
    ax3,_=plotc.colorplot(update,sp=next(gl),x=x,y=z,title="Update psi",
    yr=[-5,None],colorbar_properties={'orientation':'horizontal','shrink':0.8},
    axes_properties={'sharey':ax1,'hide_yticklabels':True})
    plt.xlabel("Horizontal Distance (Mm)",fontsize=16,labelpad=10)
    
elif enf_cont and (contvar == 'vx'):
    true_model = fitsread(os.path.join(codedir,'true_vx.fits'))
    current_model = fitsread(os.path.join(datadir,'model_vx_ls00.fits'))
    try:
        update=fitsread(os.path.join(datadir,'update','update_vx_'+iterm1+'.fits'))
    except IOError:
        update=np.zeros_like(current_model)
        print "Could not load update vx"
    _,_,gl=plotc.layout_subplots(3)
    plotc.colorplot(true_model,sp=next(gl),x=x,y=z,title="True vx",
    yr=[-5,None],colorbar_properties={'orientation':'horizontal','shrink':0.8},
    axes_properties={'sharey':ax1,'hide_yticklabels':True})
    plt.ylabel("Depth (Mm)",fontsize=16)
    plt.xlabel("Horizontal Distance (Mm)",fontsize=16,labelpad=10)
    plotc.colorplot(current_model,sp=next(gl),x=x,y=z,title="Current vx",
    yr=[-5,None],colorbar_properties={'orientation':'horizontal','shrink':0.8},
    axes_properties={'sharey':ax1,'hide_yticklabels':True})
    plotc.colorplot(update,sp=next(gl),x=x,y=z,title="Update vx",
    yr=[-5,None],colorbar_properties={'orientation':'horizontal','shrink':0.8},
    axes_properties={'sharey':ax1,'hide_yticklabels':True})
    
elif enf_cont and (contvar == 'vz'):
    true_model = fitsread(os.path.join(codedir,'true_vz.fits'))
    current_model = fitsread(os.path.join(datadir,'model_vz_ls00.fits'))
    try:
        update=fitsread(os.path.join(datadir,'update','update_vz_'+iterm1+'.fits'))
    except IOError:
        update=np.zeros_like(current_model)
        print "Could not load update vz"
    _,_,gl=plotc.layout_subplots(3)
    plotc.colorplot(true_model,sp=next(gl),x=x,y=z,title="True vz",
    yr=[-5,None],colorbar_properties={'orientation':'horizontal','shrink':0.8},
    axes_properties={'sharey':ax1,'hide_yticklabels':True})
    plt.ylabel("Depth (Mm)",fontsize=16)
    plt.xlabel("Horizontal Distance (Mm)",fontsize=16,labelpad=10)
    plotc.colorplot(current_model,sp=next(gl),x=x,y=z,title="Current vz",
    yr=[-5,None],colorbar_properties={'orientation':'horizontal','shrink':0.8},
    axes_properties={'sharey':ax1,'hide_yticklabels':True})
    plotc.colorplot(update,sp=next(gl),x=x,y=z,title="Update vz",
    yr=[-5,None],colorbar_properties={'orientation':'horizontal','shrink':0.8},
    axes_properties={'sharey':ax1,'hide_yticklabels':True})
    
elif not enf_cont:    
    true_model_vx = fitsread(os.path.join(codedir,'true_vx.fits'))
    true_model_vz = fitsread(os.path.join(codedir,'true_vz.fits'))
    
    current_model_vx = fitsread(os.path.join(datadir,'model_vx_ls00.fits'))
    current_model_vz = fitsread(os.path.join(datadir,'model_vz_ls00.fits'))
    
    try:
        update_vx=fitsread(os.path.join(datadir,'update','update_vx_'+iterm1+'.fits'))
    except IOError:
        update_vx=np.zeros_like(current_model)
        print "Could not load update vx"
    try:
        update_vz=fitsread(os.path.join(datadir,'update','update_vz_'+iterm1+'.fits'))
    except IOError:
        update_vz=np.zeros_like(current_model)
        print "Could not load update vz"
        
    _,_,gl=plotc.layout_subplots(6)
    plotc.colorplot(true_model_vx,sp=next(gl),x=x,y=z,title="True vx",
    yr=[-5,None],colorbar_properties={'orientation':'horizontal','shrink':0.8},
    axes_properties={'sharey':ax1,'hide_yticklabels':True})
    plt.ylabel("Depth (Mm)",fontsize=16)
    plotc.colorplot(current_model_vx,sp=next(gl),x=x,y=z,title="Current vx",
    yr=[-5,None],colorbar_properties={'orientation':'horizontal','shrink':0.8},
    axes_properties={'sharey':ax1,'hide_yticklabels':True})
    plotc.colorplot(update_vx,sp=next(gl),x=x,y=z,title="Update vx",
    yr=[-5,None],colorbar_properties={'orientation':'horizontal','shrink':0.8},
    axes_properties={'sharey':ax1,'hide_yticklabels':True})
    
    plotc.colorplot(true_model_vz,sp=next(gl),x=x,y=z,title="True vz",
    yr=[-5,None],colorbar_properties={'orientation':'horizontal','shrink':0.8},
    axes_properties={'sharey':ax1,'hide_yticklabels':True})
    plt.ylabel("Depth (Mm)",fontsize=16)
    plt.xlabel("Horizontal Distance (Mm)",fontsize=16,labelpad=10)
    plotc.colorplot(current_model_vz,sp=next(gl),x=x,y=z,title="Current vz",
    yr=[-5,None],colorbar_properties={'orientation':'horizontal','shrink':0.8},
    axes_properties={'sharey':ax1,'hide_yticklabels':True})
    plt.xlabel("Horizontal Distance (Mm)",fontsize=16,labelpad=10)
    plotc.colorplot(update_vz,sp=next(gl),x=x,y=z,title="Update vz",
    yr=[-5,None],colorbar_properties={'orientation':'horizontal','shrink':0.8},
    axes_properties={'sharey':ax1,'hide_yticklabels':True})
    plt.xlabel("Horizontal Distance (Mm)",fontsize=16,labelpad=10)
    
    plt.subplots_adjust(hspace=0)


plt.subplots_adjust(wspace=0)
plt.show()
    
    


