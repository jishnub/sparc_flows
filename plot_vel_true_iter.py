from __future__ import division
import plotc
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np
import read_params
import pyfits
import os,fnmatch,sys
import warnings

matplotlib.rcParams['axes.unicode_minus'] = False

def fitsread(f): return np.squeeze(pyfits.getdata(f))
    
datadir=read_params.get_directory()
codedir=os.path.dirname(os.path.abspath(__file__))

def get_iter_no():
    updatedir=os.path.join(datadir,"update")
    # Count the number of misfit_xx files
    return len(fnmatch.filter(os.listdir(updatedir),'misfit_[0-9][0-9]'))
iterno=get_iter_no()

itercutoff = filter(lambda x: x.startswith("iter="),sys.argv)
if len(itercutoff)!=0: itercutoff=int(itercutoff[0].split("=")[-1])
else: itercutoff=np.inf

def iterm(i): 
    return str(min(itercutoff,iterno-i)).zfill(2)

enf_cont = read_params.get_enforced_continuity()
contvar=read_params.get_continuity_variable()

Lx=read_params.get_xlength()
nx=read_params.get_nx()
x=np.linspace(-Lx/2,Lx/2,num=nx,endpoint=False)

Rsun=695.8
z=np.loadtxt(os.path.join(codedir,read_params.get_solarmodel()),usecols=[0])
z=(z-1)*Rsun

def forceAspect(ax,aspect=1):
    extent =  ax.get_xlim()+ax.get_ylim()
    ax.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)

if enf_cont and (contvar == 'psi'):
    true_vx = fitsread(os.path.join(codedir,'true_vx.fits'))
    true_vz = fitsread(os.path.join(codedir,'true_vz.fits'))
    current_vx = fitsread(os.path.join(datadir,'update','vx_'+iterm(1)+'.fits'))
    current_vz = fitsread(os.path.join(datadir,'update','vz_'+iterm(1)+'.fits'))
    try:
        prev1_vx = fitsread(os.path.join(datadir,'update','vx_'+iterm(2)+'.fits'))
        prev1_vz = fitsread(os.path.join(datadir,'update','vz_'+iterm(2)+'.fits'))
        prev2_vx = fitsread(os.path.join(datadir,'update','vx_'+iterm(3)+'.fits'))
        prev2_vz = fitsread(os.path.join(datadir,'update','vz_'+iterm(3)+'.fits'))
    except IOError:
        prev1_vx=None
        prev2_vx=None
        prev1_vz=None
        prev2_vz=None
    
    
    warnings.filterwarnings("ignore", message="Unicode equal comparison failed")
    
    yr_plot = [-5,None]
    
    ax1,cb1=plotc.colorplot(true_vx,sp=231,x=x,y=z,yr=yr_plot,
    centerzero=True,colorbar_properties={"orientation":"horizontal","shrink":0.8,"pad":0.25})[0:2]
    
    plt.title(r"True vx",y=1.01)
    
    ax2,cb2,qm=plotc.colorplot(current_vx,sp=232,x=x,y=z,yr=yr_plot,
    centerzero=True,colorbar_properties={"orientation":"horizontal","shrink":0.8,"pad":0.25})
    
    plt.title(r"Iterated vx",y=1.01)
    
    plt.subplot(233)
    model_misfit_initial = np.trapz(true_vx**2,x=x,axis=1)
    max_misfit_init = model_misfit_initial.max()
    model_misfit_initial/= max_misfit_init
    model_misfit = np.trapz((true_vx - current_vx)**2,x=x,axis=1)
    model_misfit/= max_misfit_init
    plt.plot(model_misfit_initial,z,color='black',linestyle='dashed',label="iter 0")
    plt.plot(model_misfit,z,color='black',linestyle='solid',label="iter "+iterm(1))
    plt.legend(loc="lower right",fontsize=14)
    plt.title("Misfit")
    ax=plt.gca()
    ax.xaxis.set_major_locator(MaxNLocator(5,prune="both"))
    ax.yaxis.set_major_locator(MaxNLocator(5,prune="both"))
    plt.ylim(-5,z.max())
    cb3=plt.colorbar(mappable=qm,ax=ax,orientation="horizontal",shrink=0.8,pad=0.25)
    #~ cb3.ax.set_visible(False)
    
    ax4,cb4=plotc.colorplot(true_vz,sp=234,x=x,y=z,yr=yr_plot,
    centerzero=True,colorbar_properties={"orientation":"horizontal","shrink":0.8,"pad":0.25})[0:2]
    
    plt.title(r"True vz",y=1.01)
    
    ax5,cb5=plotc.colorplot(current_vz,sp=235,x=x,y=z,yr=yr_plot,centerzero=True,
    colorbar_properties={"orientation":"horizontal","shrink":0.8,"pad":0.25})[0:2]
    
    plt.title(r"Iterated vz",y=1.01)
    
    plt.subplot(236)
    model_misfit_initial = np.trapz(true_vz**2,x=x,axis=-1)
    max_misfit_init = model_misfit_initial.max()
    model_misfit_initial/= max_misfit_init
    model_misfit = np.trapz((true_vz - current_vz)**2,x=x,axis=-1)
    model_misfit/= max_misfit_init
    plt.plot(model_misfit_initial,z,color='black',linestyle='dashed',label="iter 0")
    plt.plot(model_misfit,z,color='black',linestyle='solid',label="iter "+iterm(1))
    plt.legend(loc="lower right",fontsize=14)
    plt.title("Misfit")
    ax=plt.gca()
    ax.xaxis.set_major_locator(MaxNLocator(5,prune="both"))
    ax.yaxis.set_major_locator(MaxNLocator(5,prune="both"))
    plt.ylim(-5,z.max())
    cb6=plt.colorbar(mappable=qm,ax=ax,orientation="horizontal",shrink=0.8,pad=0.25)
    
    
    for ax in [ax1,ax2,ax4,ax5]:
        ax.set_xlim(-Lx/8,Lx/8)
        ax.set_ylim(-5,z.max())
        ax.set_xlabel("x (Mm)")
    
    ax1.set_ylabel("Depth (Mm)")
    ax4.set_ylabel("Depth (Mm)")
        
    for cb in [cb1,cb2,cb4,cb5]:
        cb.ax.xaxis.set_major_locator(MaxNLocator(4,prune='both'))
        cb.ax.set_ylabel("m/s",rotation=90,fontsize=14)
        
    plotc.apj_2col_format(plt.gcf(),default_fontsize=18)
    plt.tight_layout()
    cb3.ax.set_visible(False)
    cb6.ax.set_visible(False)
    
    if not os.path.exists("plots"): os.makedirs("plots")
    plt.savefig("plots/f7.eps")
    
    #~ plt.figure()
    #~ vx_max_row_index,vx_max_col_index = divmod(true_vx.argmax(),nx)
    #~ vz_max_row_index,vz_max_col_index = divmod(true_vz.argmax(),nx)
    #~ 
    #~ plt.subplot(121)
    #~ plt.plot(z,true_vx[:,vx_max_col_index],label="True")
    #~ plt.plot(z,current_vx[:,vx_max_col_index],label="Iter",
    #~ linestyle='dashed',linewidth=2)
#~ 
    #~ plt.xlabel("Depth (Mm)",fontsize=20)
    #~ plt.ylabel("vx (m/s)",fontsize=20)
    #~ plt.legend(loc='best',fontsize=18)
    #~ plt.xlim(-8,2)
    #~ plt.gca().yaxis.set_major_locator(MaxNLocator(5,prune='both'))
    #~ plt.tick_params(axis='both', which='major', labelsize=14)
#~ 
    #~ plt.subplot(122)
    #~ plt.plot(z,true_vz[:,vz_max_col_index],label="True")
    #~ plt.plot(z,current_vz[:,vz_max_col_index],label="Iter",
    #~ linestyle='dashed',linewidth=2)
#~ 
    #~ plt.xlabel("Depth (Mm)",fontsize=20)
    #~ plt.ylabel("vz (m/s)",fontsize=20)
    #~ plt.legend(loc='best',fontsize=18)
    #~ plt.xlim(-8,2)
    #~ plt.gca().yaxis.set_major_locator(MaxNLocator(5,prune='both'))
    #~ plt.tick_params(axis='both', which='major', labelsize=14)
    #~ 
    #~ plt.tight_layout()
    #~ plt.subplots_adjust(wspace=0.3)
    
    
    
    
    
    
elif enf_cont and (contvar == 'vx'):
    true_model = fitsread(os.path.join(codedir,'true_vx.fits'))
    current_model = fitsread(os.path.join(datadir,'model_vx_ls00.fits'))
    try:
        update=fitsread(os.path.join(datadir,'update','update_vx_'+iterm1+'.fits'))
    except IOError:
        update=np.zeros_like(current_model)
        print "Could not load update vx"
    
    _,_,gl=plotc.layout_subplots(3)
    
    ax1=plotc.colorplot(true_model,sp=next(gl),x=x,y=z,title="True vx",
    yr=[-5,None],colorbar_properties={'orientation':'horizontal','shrink':0.8},
    axes_properties={})[0]
    plt.ylabel("Depth (Mm)",fontsize=16)
    plt.xlabel("Horizontal Distance (Mm)",fontsize=16,labelpad=10)
    
    ax2=plotc.colorplot(current_model,sp=next(gl),x=x,y=z,title="Current vx",
    yr=[-5,None],colorbar_properties={'orientation':'horizontal','shrink':0.8},
    axes_properties={'sharey':ax1,'hide_yticklabels':True})[0]
    plt.xlabel("Horizontal Distance (Mm)",fontsize=16,labelpad=10)
    
    ax3=plotc.colorplot(update,sp=next(gl),x=x,y=z,title="Update vx",
    yr=[-5,None],colorbar_properties={'orientation':'horizontal','shrink':0.8},
    axes_properties={'sharey':ax1,'hide_yticklabels':True})[0]
    plt.xlabel("Horizontal Distance (Mm)",fontsize=16,labelpad=10)
    
    plt.suptitle("Continuity: Compute vz from vx",fontsize=16)
    plt.subplots_adjust(wspace=0)
    
elif enf_cont and (contvar == 'vz'):
    true_model = fitsread(os.path.join(codedir,'true_vz.fits'))
    current_model = fitsread(os.path.join(datadir,'model_vz_ls00.fits'))
    try:
        update=fitsread(os.path.join(datadir,'update','update_vz_'+iterm1+'.fits'))
    except IOError:
        update=np.zeros_like(current_model)
        print "Could not load update vz"
    
    _,_,gl=plotc.layout_subplots(3)
    
    ax1=plotc.colorplot(true_model,sp=next(gl),x=x,y=z,title="True vz",
    yr=[-5,None],colorbar_properties={'orientation':'horizontal','shrink':0.8},
    axes_properties={})[0]
    plt.ylabel("Depth (Mm)",fontsize=16)
    plt.xlabel("Horizontal Distance (Mm)",fontsize=16,labelpad=10)
    
    ax2=plotc.colorplot(current_model,sp=next(gl),x=x,y=z,title="Current vz",
    yr=[-5,None],colorbar_properties={'orientation':'horizontal','shrink':0.8},
    axes_properties={'sharey':ax1,'hide_yticklabels':True})[0]
    plt.xlabel("Horizontal Distance (Mm)",fontsize=16,labelpad=10)
    
    ax3=plotc.colorplot(update,sp=next(gl),x=x,y=z,title="Update vz",
    yr=[-5,None],colorbar_properties={'orientation':'horizontal','shrink':0.8},
    axes_properties={'sharey':ax1,'hide_yticklabels':True})[0]
    plt.xlabel("Horizontal Distance (Mm)",fontsize=16,labelpad=10)
    
    plt.suptitle("Continuity: Compute vx from vz",fontsize=16)
    plt.subplots_adjust(wspace=0)
    
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
    
    ax1=plotc.colorplot(true_model_vx,sp=next(gl),x=x,y=z,title="True vx",
    yr=[-5,None],colorbar_properties={'orientation':'horizontal','shrink':0.8},
    centerzero=True)[0]
    plt.ylabel("Depth (Mm)",fontsize=16)
    
    ax2=plotc.colorplot(current_model_vx,sp=next(gl),x=x,y=z,title="Current vx",
    yr=[-5,None],colorbar_properties={'orientation':'horizontal','shrink':0.8},
    axes_properties={'sharey':ax1,'hide_yticklabels':True},centerzero=True)[0]
    
    ax3=plotc.colorplot(update_vx,sp=next(gl),x=x,y=z,title="Update vx",
    yr=[-5,None],colorbar_properties={'orientation':'horizontal','shrink':0.8},
    axes_properties={'sharey':ax1,'hide_yticklabels':True},centerzero=True)[0]
    
    ax4=plotc.colorplot(true_model_vz,sp=next(gl),x=x,y=z,title="True vz",
    yr=[-5,None],colorbar_properties={'orientation':'horizontal','shrink':0.8},
    axes_properties={'sharey':ax1},centerzero=True)[0]
    plt.ylabel("Depth (Mm)",fontsize=16)
    plt.xlabel("Horizontal Distance (Mm)",fontsize=16,labelpad=5)
    
    ax5=plotc.colorplot(current_model_vz,sp=next(gl),x=x,y=z,title="Current vz",
    yr=[-5,None],colorbar_properties={'orientation':'horizontal','shrink':0.8},
    axes_properties={'sharey':ax1,'hide_yticklabels':True},centerzero=True)[0]
    plt.xlabel("Horizontal Distance (Mm)",fontsize=16,labelpad=5)
    
    ax6=plotc.colorplot(update_vz,sp=next(gl),x=x,y=z,title="Update vz",
    yr=[-5,None],colorbar_properties={'orientation':'horizontal','shrink':0.8},
    axes_properties={'sharey':ax1,'hide_yticklabels':True},centerzero=True)[0]
    plt.xlabel("Horizontal Distance (Mm)",fontsize=16,labelpad=5)
    
    plt.suptitle("No continuity, vx and vz unconstrained",fontsize=16)
    
    plt.tight_layout()
    plt.subplots_adjust(wspace=0)



plt.show()
    
    


