
import plotc
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np
import read_params
import pyfits
import os,fnmatch,sys
import warnings

matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['font.family'] = "serif"

def fitsread(f): return np.squeeze(pyfits.getdata(f))
    
datadir=read_params.get_directory()

def get_iter_no():
    updatedir=os.path.join(datadir,"update")
    # Count the number of misfit_xx files
    return len(fnmatch.filter(os.listdir(updatedir),'misfit_[0-9][0-9]'))
iterno=get_iter_no()

itercutoff = read_params.parse_cmd_line_params("iter",default=str(iterno-1).zfill(2))

Lx=read_params.get_xlength()
nx=read_params.get_nx()
x=np.linspace(-Lx/2,Lx/2,num=nx,endpoint=False)
dx = x[1]-x[0]
x_plot = np.append(x-dx/2,x[-1]+dx/2)

Rsun=695.8
z=np.loadtxt(read_params.get_solarmodel(),usecols=[0])
z=(z-1)*Rsun
z_pixels = np.arange(len(z))
z_plot = np.interp(z_pixels[:-1]+0.5,z_pixels,z)
z_plot = np.append(z[0]-(z[1]-z[0])/2,np.append(z_plot,z[-1]+(z[-1]-z[-2])/2))


def forceAspect(ax,aspect=1):
    extent =  ax.get_xlim()+ax.get_ylim()
    ax.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)

true_vx = fitsread('true_vx.fits')
true_vz = fitsread('true_vz.fits')
current_vx = fitsread(os.path.join(datadir,'update','vx_'+itercutoff+'.fits'))
current_vz = fitsread(os.path.join(datadir,'update','vz_'+itercutoff+'.fits'))


warnings.filterwarnings("ignore", message="Unicode equal comparison failed")

yr_plot = [-5,None]

axlist = []
cblist = []
def add_subplot(sp_index,array,title=""):
    ax=plt.subplot(sp_index)
    array_max = abs(array).max()
    qm=ax.pcolorfast(x_plot,z_plot,array,vmax=array_max,vmin=-array_max,cmap="RdBu_r")
    cb=plt.colorbar(mappable=qm,orientation="horizontal",shrink=0.8,pad=0.25,ticks=MaxNLocator(5))
    plt.title(title,y=1.01,fontsize=18)
    axlist.append(ax)
    cblist.append(cb)

add_subplot(221,true_vx,title="True $v_x$")
add_subplot(222,current_vx,title="Iterated $v_x$")
add_subplot(223,true_vz,title="True $v_z$")
add_subplot(224,current_vz,title="Iterated $v_z$")

for ax in axlist:
    ax.set_xlim(-Lx/15,Lx/15)
    ax.set_ylim(-5,z.max())
    ax.set_xlabel("x (Mm)",fontsize=18)
    ax.xaxis.set_major_locator(MaxNLocator(4,prune="both"))
    ax.yaxis.set_major_locator(MaxNLocator(5,prune="both"))
    ax.tick_params(axis='both', labelsize=16)

axlist[0].set_ylabel("Depth (Mm)",fontsize=18)
plt.setp(axlist[1].get_yticklabels(),visible=False)
axlist[2].set_ylabel("Depth (Mm)",fontsize=18)
plt.setp(axlist[3].get_yticklabels(),visible=False)
    
for cb in cblist: 
    cb.ax.set_ylabel("$\mathrm{m}/\mathrm{s}$",rotation=90,fontsize=16)
    cb.ax.tick_params(axis="x",labelsize=16)
    
plt.gcf().set_size_inches(8,6)
plt.tight_layout()
plt.subplots_adjust(wspace=0,hspace=0.3)
    
save = read_params.parse_cmd_line_params("save")
if save is not None:
    savepath = os.path.join("plots",save)
    print("saving to",savepath)
    if not os.path.exists("plots"): os.makedirs("plots")
    plt.savefig(savepath)
else:
    print("Not saving plot to file")

plt.show()
    
    


