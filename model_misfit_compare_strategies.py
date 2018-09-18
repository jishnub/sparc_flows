
import numpy as np
import pyfits
import os,re,sys,glob,fnmatch
import read_params
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator,LogLocator
from matplotlib import rc
import itertools
import fnmatch

rc("font",family="serif")

Rsun = 695.8
Lx = read_params.get_xlength()
nx = read_params.get_nx()
x = np.linspace(-Lx/2,Lx/2,nx,endpoint=False)
z = np.loadtxt(read_params.get_solarmodel(),usecols=[0]); z=(z-1)*Rsun

def model_misfit_depth(directory,iterno=None,var="psi"):
    datadir = os.path.join(os.path.dirname(read_params.get_directory()),directory)
    true_model = np.squeeze(pyfits.getdata("true_"+var+".fits"))
    updatedir = os.path.join(datadir,"update")
    if var=="psi":
        iter_model_list = sorted(fnmatch.filter(os.listdir(updatedir),"model_psi_[0-9][0-9].fits"))
    else:
        iter_model_list = sorted(fnmatch.filter(os.listdir(updatedir),var+"_[0-9][0-9].fits"))
        
    if iterno is None: iter_model_name = iter_model_list[-1]; iterno=len(iter_model_list)
    else: iter_model_name = iter_model_list[iterno]
    iterated_model = np.squeeze(pyfits.getdata(os.path.join(datadir,"update",iter_model_name)))
    iterated_model -= iterated_model[0,0]

    true_misfit = np.trapz(true_model**2,x=x,axis=1)
    iter_misfit = np.trapz((true_model-iterated_model)**2,x=x,axis=1)
    return true_misfit/true_misfit.max(),iter_misfit/true_misfit.max(),iterno

def plot_misfit_z(misfit,**kwargs):
    plt.plot(misfit,z,linewidth=2,color=kwargs.pop("color","black"),**kwargs)
    plt.ylim(-8.,z[-1])
    ax.xaxis.set_major_locator(MaxNLocator(3,prune="both"))
    ax.yaxis.set_major_locator(MaxNLocator(6))
    ax.grid(True)

def model_misfit_strategy(directory,itercutoff=None):
    datadir = os.path.join(os.path.dirname(read_params.get_directory()),directory)
    true_model = np.squeeze(pyfits.getdata("true_psi.fits"))
    updatedir = os.path.join(datadir,"update")
    iter_model_list = sorted(fnmatch.filter(os.listdir(updatedir),"model_psi_[0-9][0-9].fits"))
    
    misfit = []
    true_misfit = np.sum(true_model**2)
    for model_name in iter_model_list:
        iter_model = np.squeeze(pyfits.getdata(os.path.join(updatedir,model_name)))
        iter_model -= iter_model[0,0]
        misfit.append( np.sum((true_model - iter_model)**2)/true_misfit)
        
    return misfit if itercutoff is None else misfit[:itercutoff+1]
    
def model_misfit_psi_v(directory,itercutoff=None,below_surf=False):
    datadir = os.path.join(os.path.dirname(read_params.get_directory()),directory)
    updatedir = os.path.join(datadir,"update")
    true_model = np.squeeze(pyfits.getdata("true_psi.fits"))
    iter_model_list = sorted(fnmatch.filter(os.listdir(updatedir),"model_psi_[0-9][0-9].fits"))
    if itercutoff is not None: iter_model_list=iter_model_list[:itercutoff+1]
    misfit_psi = []
    true_misfit = np.sum(true_model**2)
    for model_name in iter_model_list:
        iter_model = np.squeeze(pyfits.getdata(os.path.join(updatedir,model_name)))
        iter_model -= iter_model[0,0]
        squared_diff = (true_model - iter_model)**2
        if below_surf: squared_diff[z>0]=0
        misfit_psi.append( np.sum(squared_diff)/true_misfit)
        
    true_model = np.squeeze(pyfits.getdata("true_vz.fits"))
    iter_model_list = sorted(fnmatch.filter(os.listdir(updatedir),"vz_[0-9][0-9].fits"))
    if itercutoff is not None: iter_model_list=iter_model_list[:itercutoff+1]
    misfit_vz = []
    true_misfit = np.sum(true_model**2)
    for model_name in iter_model_list:
        iter_model = np.squeeze(pyfits.getdata(os.path.join(updatedir,model_name)))
        squared_diff = (true_model - iter_model)**2
        if below_surf: squared_diff[z>0]=0
        misfit_vz.append( np.sum(squared_diff)/true_misfit)
    
    
    true_model = np.squeeze(pyfits.getdata("true_vx.fits"))
    iter_model_list = sorted(fnmatch.filter(os.listdir(updatedir),"vx_[0-9][0-9].fits"))
    if itercutoff is not None: iter_model_list=iter_model_list[:itercutoff+1]
    misfit_vx = []
    true_misfit = np.sum(true_model**2)
    for model_name in iter_model_list:
        iter_model = np.squeeze(pyfits.getdata(os.path.join(updatedir,model_name)))
        squared_diff = (true_model - iter_model)**2
        if below_surf: squared_diff[z>0]=0
        misfit_vx.append( np.sum(squared_diff)/true_misfit)
        
    return misfit_psi,misfit_vx,misfit_vz
    
ax=plt.subplot2grid((2,3),(0,0),rowspan=2)
plot_misfit_z(model_misfit_depth("f_p1")[1],ls="solid",label=r"#$1$")
plot_misfit_z(model_misfit_depth("f_to_p3",iterno=35)[1],ls="dashed",label=r"#$2$")
plot_misfit_z(model_misfit_depth("f_to_p7_2pixsmooth",iterno=35)[1],ls="dotted",label=r"#$3$")
plot_misfit_z(model_misfit_depth("f_to_p7_new",iterno=35)[1],ls="dotted",label=r"#$4$",marker='o',markersize=3)
plt.legend(loc="lower right",fontsize=20)
plt.ylabel("Depth (Mm)",fontsize=20)
plt.title("Misfit in $\psi$",fontsize=20)
plt.tick_params(axis="both",labelsize=14)
plt.text(0.8,0.92,"(a)",transform=ax.transAxes,fontsize=20,weight="bold")

ax=plt.subplot2grid((2,3),(0,1))
true_misfit,iter_misfit,iterno = model_misfit_depth("f_p1",iterno=15,var="vx")
plot_misfit_z(true_misfit,ls="solid",color="#777777",label="Iter 0")
plot_misfit_z(iter_misfit,ls="solid",label="#$1$, Iter "+str(iterno))
true_misfit,iter_misfit,iterno = model_misfit_depth("f_to_p3",iterno=35,var="vx")
plot_misfit_z(iter_misfit,ls="dotted",label="#$3$, Iter "+str(iterno))
plt.ylabel("Depth (Mm)",fontsize=20)
plt.title("Misfit in $v_x$",fontsize=20)
plt.tick_params(axis="both",labelsize=14)
plt.ylim(-6,z[-1])
plt.legend(loc="lower right")
plt.text(0.8,0.85,"(b)",transform=ax.transAxes,fontsize=20,weight="bold")

ax=plt.subplot2grid((2,3),(0,2))
true_misfit,iter_misfit,iterno = model_misfit_depth("f_p1",iterno=15,var="vz")
plot_misfit_z(true_misfit,ls="solid",color="#777777",label="Iter 0")
plot_misfit_z(iter_misfit,ls="solid",label="#$1$, Iter "+str(iterno))
true_misfit,iter_misfit,iterno = model_misfit_depth("f_to_p3",iterno=35,var="vz")
plot_misfit_z(iter_misfit,ls="dotted",label="#$3$, Iter "+str(iterno))
plt.title("Misfit in $v_z$",fontsize=20)
plt.ylabel("Depth (Mm)",fontsize=20)
plt.ylim(-6,z[-1])
plt.legend(loc="lower right")
plt.tick_params(axis="both",labelsize=14)
plt.text(0.8,0.85,"(c)",transform=ax.transAxes,fontsize=20,weight="bold")

ax=plt.subplot2grid((2,3),(1,1))
plt.plot(model_misfit_strategy("f_p1"),color="black",ls="solid",linewidth=2,label=r"#$1$")
plt.plot(model_misfit_strategy("f_to_p3",itercutoff=35),color="black",ls="dashed",linewidth=2,label=r"#$2$")
plt.plot(model_misfit_strategy("f_to_p7_2pixsmooth",itercutoff=35),color="black",ls="dotted",linewidth=2,label=r"#$3$")
plt.plot(model_misfit_strategy("f_to_p7_new",itercutoff=35),color="black",ls="dotted",linewidth=2,label=r"#$4$",marker='o',markersize=3)
plt.ylim(0,1)
ax.yaxis.set_major_locator(MaxNLocator(5))
plt.legend(loc="upper right",ncol=2,fontsize=16)
plt.title("Total misfit in $\psi$",fontsize=20)
plt.xlabel("Iteration",fontsize=20)
plt.ylabel("$\chi^{\psi}$",fontsize=20)
plt.tick_params(axis="both",labelsize=14)
plt.text(0.1,0.1,"(d)",transform=ax.transAxes,fontsize=20,weight="bold")

ax=plt.subplot2grid((2,3),(1,2))
misfit_psi,misfit_vx,misfit_vz=model_misfit_psi_v("f_to_p3",itercutoff=35)
plt.plot(misfit_psi,color="black",linestyle="solid",linewidth=2,label="$\psi$")
plt.plot(misfit_vx,color="black",linestyle="dashed",linewidth=2,label="$v_x$")
plt.plot(misfit_vz,color="black",linestyle="dashdot",linewidth=2,label="$v_z$")
misfit_psi,misfit_vx,misfit_vz=model_misfit_psi_v("f_to_p3",itercutoff=35,below_surf=True)
plt.plot(misfit_vx,color="black",linestyle="None",marker="o",markersize=3,label=r"$\tilde{v}_x$")
plt.ylim(0.2,1)
ax.yaxis.set_major_locator(MaxNLocator(5))
plt.legend(loc="upper right",ncol=2,fontsize=16)
plt.title("Strategy #$3$",fontsize=20)
plt.xlabel("Iteration",fontsize=20)
plt.ylabel("$\chi^m$",fontsize=20)
plt.tick_params(axis="both",labelsize=14)
plt.text(0.1,0.1,"(e)",transform=ax.transAxes,fontsize=20,weight="bold")

plt.tight_layout()
plt.subplots_adjust(wspace=0.5)
plt.gcf().set_size_inches(15,9)
if not os.path.exists("plots"): os.makedirs("plots")
plt.savefig("plots/strategy_comp.eps")

#~ plt.show()

