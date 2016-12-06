from __future__ import division
import plotc
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator,ScalarFormatter,NullLocator
import numpy as np
import read_params
import pyfits
import os,fnmatch,sys

def fitsread(f): return np.squeeze(pyfits.getdata(f))

datadir=read_params.get_directory()
codedir=os.path.dirname(os.path.abspath(__file__))

def get_iter_no():
    updatedir=os.path.join(datadir,"update")
    # Count the number of misfit_xx files
    return len(fnmatch.filter(os.listdir(updatedir),'misfit_[0-9][0-9]'))
iterno=get_iter_no()
iterm1=str(iterno-1).zfill(2)

sound_speed_perturbed = read_params.if_soundspeed_perturbed()
flows = read_params.if_flows()
enf_cont = read_params.if_continuity_enforced()
contvar=read_params.get_continuity_variable()

Lx=read_params.get_xlength()
nx=read_params.get_nx()
x=np.linspace(-Lx/2,Lx/2,num=nx,endpoint=False)

Rsun=695.8
z=np.loadtxt(os.path.join(codedir,read_params.get_solarmodel()),usecols=[0])
z=(z-1)*Rsun

#~ Check if specific iter is to be plotted
iter_to_plot=read_params.parse_cmd_line_params("iter",zfill=2,default=iterm1)

if flows:
    if enf_cont and (contvar == 'psi'):
        true_psi = fitsread(os.path.join(codedir,'true_psi.fits'))
        current_psi = fitsread(os.path.join(datadir,'update',
                        'model_psi_'+iter_to_plot+'.fits'))
        current_psi = current_psi-current_psi[0,0]

        true_vx = fitsread('true_vx.fits')
        true_vz = fitsread('true_vz.fits')

        current_vx = fitsread(
        os.path.join(datadir,'update','vx_{}.fits'.format(iter_to_plot)))

        current_vz = fitsread(
        os.path.join(datadir,'update','vz_{}.fits'.format(iter_to_plot)))

        plot_update=True
        gl=plotc.layout_subplots(2)[2]

        cp=plotc.colorplot(true_psi,sp=next(gl),x=x,y=z,
        yr=[-5,None],colorbar_properties={'orientation':'horizontal','shrink':0.8})

        ax1=cp.axis; cb=cp.colorbar
        cb.ax.xaxis.set_major_locator(MaxNLocator(3))

        plt.ylabel("Depth (Mm)",fontsize=20)
        plt.xlabel("x (Mm)",fontsize=20,labelpad=10)
        plt.title(r"True $\psi$",fontsize=20,y=1.01)
        plt.tick_params(axis='both', which='major', labelsize=14)
        ax1.xaxis.set_major_locator(MaxNLocator(4,prune='both'))

        cp=plotc.colorplot(current_psi,sp=next(gl),x=x,y=z,
        yr=[-5,None],colorbar_properties={'orientation':'horizontal','shrink':0.8},
        axes_properties={'sharey':ax1,'hide_yticklabels':True})

        ax2=cp.axis; cb=cp.colorbar
        cb.ax.xaxis.set_major_locator(MaxNLocator(3))

        plt.xlabel("x (Mm)",fontsize=20,labelpad=10)
        plt.title(r"Iterated $\psi$",fontsize=20,y=1.01)
        plt.tick_params(axis='both', which='major', labelsize=14)
        ax2.xaxis.set_major_locator(MaxNLocator(4,prune='both'))

        ax2.yaxis.set_major_locator(NullLocator())

        plt.subplots_adjust(wspace=0)

        plt.figure()

        psi_max_row_index,psi_max_col_index = divmod(true_psi.argmax(),nx)
        curpsi_max_row_index,curpsi_max_col_index = divmod(current_psi.argmax(),nx)

        vx_max_row_index,vx_max_col_index = divmod(true_vx.argmax(),nx)
        curvx_max_row_index,curvx_max_col_index = divmod(current_vx.argmax(),nx)

        vz_max_row_index,vz_max_col_index = divmod(true_vz.argmax(),nx)
        curvz_max_row_index,curvz_max_col_index = divmod(current_vz.argmax(),nx)

        plt.subplot(131)
        plt.plot(z,true_psi[:,curpsi_max_col_index],label="Reference")
        plt.plot(z,current_psi[:,curpsi_max_col_index],label="Iterated",
        linestyle='dashed',linewidth=2)
        plt.xlim(-10,2.5)
        plt.xlabel("Depth (Mm)",fontsize=20)
        plt.title(r"$\psi$",fontsize=20)
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        plt.tick_params(axis='both', which='major', labelsize=14)
        plt.legend(loc='best')

        plt.subplot(132)
        plt.plot(z,true_vx[:,curvx_max_col_index],label="Reference")
        plt.plot(z,current_vx[:,curvx_max_col_index],label="Iterated",
        linestyle='dashed',linewidth=2)
        plt.xlim(-10,2.5)
        plt.xlabel("Depth (Mm)",fontsize=20)
        plt.title("vx",fontsize=20)
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        plt.tick_params(axis='both', which='major', labelsize=14)
        plt.legend(loc='best')

        plt.subplot(133)
        plt.plot(z,true_vz[:,curvz_max_col_index],label="Reference")
        plt.plot(z,current_vz[:,curvz_max_col_index],label="Iterated",
        linestyle='dashed',linewidth=2)
        plt.xlim(-10,2.5)
        plt.xlabel("Depth (Mm)",fontsize=20)
        plt.title("vz",fontsize=20)
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        plt.tick_params(axis='both', which='major', labelsize=14)
        plt.legend(loc='best')

        for ax in plt.gcf().axes:
            ax.xaxis.set_major_locator(MaxNLocator(4))

        plt.tight_layout()


    elif enf_cont and (contvar == 'vx'):
        true_psi = fitsread(os.path.join(codedir,'true_vx.fits'))
        current_psi = fitsread(os.path.join(datadir,'model_vx_ls00.fits'))
        try:
            update=fitsread(os.path.join(datadir,'update','update_vx_'+iterm1+'.fits'))
        except IOError:
            update=np.zeros_like(current_psi)
            print "Could not load update vx"

        _,_,gl=plotc.layout_subplots(3)

        ax1=plotc.colorplot(true_psi,sp=next(gl),x=x,y=z,title="True vx",
        yr=[-5,None],colorbar_properties={'orientation':'horizontal','shrink':0.8},
        axes_properties={})[0]
        plt.ylabel("Depth (Mm)",fontsize=16)
        plt.xlabel("Horizontal Distance (Mm)",fontsize=16,labelpad=10)

        ax2=plotc.colorplot(current_psi,sp=next(gl),x=x,y=z,title="Current vx",
        yr=[-5,None],colorbar_properties={'orientation':'horizontal','shrink':0.8},
        axes_properties={'sharey':ax1,'hide_yticklabels':True})[0]
        plt.xlabel("Horizontal Distance (Mm)",fontsize=16,labelpad=10)

        ax3=plotc.colorplot(update,sp=next(gl),x=x,y=z,title="Update vx",
        yr=[-5,None],colorbar_properties={'orientation':'horizontal','shrink':0.8},
        axes_properties={'sharey':ax1,'hide_yticklabels':True})[0]
        plt.xlabel("Horizontal Distance (Mm)",fontsize=16,labelpad=10)

        plt.suptitle("Continuity: Compute vz from vx",fontsize=16)

        plt.subplots_adjust(wspace=0,left=0.1,right=0.9)

    elif enf_cont and (contvar == 'vz'):
        true_psi = fitsread(os.path.join(codedir,'true_vz.fits'))
        current_psi = fitsread(os.path.join(datadir,'model_vz_ls00.fits'))
        try:
            update=fitsread(os.path.join(datadir,'update','update_vz_'+iterm1+'.fits'))
        except IOError:
            update=np.zeros_like(current_psi)
            print "Could not load update vz"

        _,_,gl=plotc.layout_subplots(3)

        ax1=plotc.colorplot(true_psi,sp=next(gl),x=x,y=z,title="True vz",
        yr=[-5,None],colorbar_properties={'orientation':'horizontal','shrink':0.8},
        axes_properties={})[0]
        plt.ylabel("Depth (Mm)",fontsize=16)
        plt.xlabel("Horizontal Distance (Mm)",fontsize=16,labelpad=10)

        ax2=plotc.colorplot(current_psi,sp=next(gl),x=x,y=z,title="Current vz",
        yr=[-5,None],colorbar_properties={'orientation':'horizontal','shrink':0.8},
        axes_properties={'sharey':ax1,'hide_yticklabels':True})[0]
        plt.xlabel("Horizontal Distance (Mm)",fontsize=16,labelpad=10)

        ax3=plotc.colorplot(update,sp=next(gl),x=x,y=z,title="Update vz",
        yr=[-5,None],colorbar_properties={'orientation':'horizontal','shrink':0.8},
        axes_properties={'sharey':ax1,'hide_yticklabels':True})[0]
        plt.xlabel("Horizontal Distance (Mm)",fontsize=16,labelpad=10)

        plt.suptitle("Continuity: Compute vx from vz",fontsize=16)

        plt.subplots_adjust(wspace=0,left=0.1,right=0.9)

    elif not enf_cont:
        true_model_vx = fitsread(os.path.join(codedir,'true_vx.fits'))
        true_model_vz = fitsread(os.path.join(codedir,'true_vz.fits'))

        current_model_vx = fitsread(os.path.join(datadir,'model_vx_ls00.fits'))
        current_model_vz = fitsread(os.path.join(datadir,'model_vz_ls00.fits'))

        try:
            update_vx=fitsread(os.path.join(datadir,'update','update_vx_'+iterm1+'.fits'))
        except IOError:
            update_vx=np.zeros_like(current_psi)
            print "Could not load update vx"
        try:
            update_vz=fitsread(os.path.join(datadir,'update','update_vz_'+iterm1+'.fits'))
        except IOError:
            update_vz=np.zeros_like(current_psi)
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

if sound_speed_perturbed:
    true_psi = fitsread(os.path.join(codedir,'true_c.fits'))
    current_psi = fitsread(os.path.join(datadir,'update','model_c_'+iter_to_plot+'.fits'))
    nz,nx=true_psi.shape
    oneD_stratification = np.tile(np.loadtxt(read_params.get_solarmodel(),usecols=[1]).reshape(nz,1),(1,nx))

    plot_update=False
    if plot_update:
        try:
            update=fitsread(os.path.join(datadir,'update','update_c_'+iter_to_plot+'.fits'))
        except IOError:
            update=None
            print "Could not load update c"
        gl=plotc.gridlist(1,4)
    else:
        gl=plotc.gridlist(1,3)

    xr_plot = [-Lx/8,Lx/8]
    yr_plot = [-10,z[-1]]

    true_lndeltac = (true_psi-oneD_stratification)/oneD_stratification
    ax1=plt.subplot(next(gl))
    cp=plotc.colorplot(true_lndeltac,ax=ax1,x=x,y=z,
    xr=xr_plot,colorbar_properties={'orientation':'horizontal','shrink':0.8})
    cb=cp.colorbar
    cb.ax.xaxis.set_major_locator(MaxNLocator(3))
    plt.ylabel("Depth (Mm)",fontsize=20)
    plt.xlabel("x (Mm)",fontsize=20,labelpad=10)
    plt.title(r"True $\delta \ln c$",fontsize=20,y=1.01)
    plt.tick_params(axis='both', which='major', labelsize=14)
    ax1.xaxis.set_major_locator(MaxNLocator(4,prune='both'))

    current_lndeltac = (current_psi-oneD_stratification)/oneD_stratification
    cp=plotc.colorplot(current_lndeltac,x=x,y=z,sp=next(gl),centerzero=True,
    xr=xr_plot,colorbar_properties={'orientation':'horizontal','shrink':0.8},
    axes_properties={'hide_yticklabels':True},subplot_properties={'sharey':ax1})
    ax2=cp.axis;cb=cp.colorbar; qm = cp.mappable
    cb.ax.xaxis.set_major_locator(MaxNLocator(3))
    plt.xlabel("x (Mm)",fontsize=20,labelpad=10)
    plt.title(r"Iterated $\delta \ln c$",fontsize=20,y=1.01)
    plt.tick_params(axis='both', which='major', labelsize=14)
    ax2.xaxis.set_major_locator(MaxNLocator(4,prune='both'))

    #~ Misfit
    ax3=plt.subplot(next(gl),sharey=ax1)
    misfit_z_0 = np.trapz(abs(true_lndeltac),x=x)
    misfit_z = np.trapz(abs(true_lndeltac-current_lndeltac),x=x)
    misfit_z /= misfit_z_0.max()
    misfit_z_0/= misfit_z_0.max()
    plt.plot(misfit_z,z,label="iter "+str(iterno-1))
    plt.plot(misfit_z_0,z,label="iter 0")
    plt.legend(loc="lower right")
    plt.ylim(*yr_plot)
    plt.tick_params(axis='x', which='major', labelsize=14)
    plt.tick_params(axis='y', which='major', label1On=False)
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    plt.title("Misfit",fontsize=20,y=1.01)
    #~ Fake hidden colorbar
    cb=plt.colorbar(mappable=qm,orientation="horizontal")
    ax3.xaxis.set_major_locator(MaxNLocator(4,prune='both'))
    ax3.yaxis.set_major_locator(MaxNLocator(5,prune='both'))
    ax3.xaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    cb.ax.set_visible(False)

    if plot_update and update is not None:
        ax4=plt.subplot(next(gl),sharey=ax1)
        cp=plotc.colorplot(update,x=x,y=z,ax=ax4,centerzero=True,
        xr=xr_plot,colorbar_properties={'orientation':'horizontal','shrink':0.8},
        axes_properties={'hide_yticklabels':True})
        ax4.set_ylim(*yr_plot)
        cb=cp.colorbar
        cb.ax.xaxis.set_major_locator(MaxNLocator(3))

        plt.xlabel("x (Mm)",fontsize=20,labelpad=10)
        plt.title(r"Next Update",fontsize=20,y=1.01)
        plt.tick_params(axis='both', which='major', labelsize=14)
        ax4.xaxis.set_major_locator(MaxNLocator(4,prune='both'))

    plt.subplots_adjust(wspace=0)


plt.show()
