from __future__ import division
import plotc
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib import gridspec
import numpy as np
import read_params
import pyfits
import os,fnmatch,sys
import warnings

matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['font.family'] = "serif"

def fitsread(f): return np.squeeze(pyfits.getdata(f))
    
datadir=read_params.get_directory()
parentdir = os.path.dirname(datadir)
strategydir=[os.path.join(parentdir,modeldir) for modeldir in ("f_p1","f_to_p3","f_to_p7_2pixsmooth","f_to_p7_new")]

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

cblist = []
axlist=[]
def add_subplot(sp_index,array,title=""):
    row_ind,col_ind=sp_index
    ax=plt.subplot(2,5,row_ind*5+col_ind+1)
    array_max = abs(array).max()
    qm=ax.pcolorfast(x_plot,z_plot,array,vmax=array_max,vmin=-array_max,cmap="RdBu_r")
    cb=plt.colorbar(mappable=qm,orientation="horizontal",shrink=0.8,pad=0.2,ticks=MaxNLocator(5))
    plt.title(title,y=1.01,fontsize=22)
    cblist.append(cb)
    axlist.append(ax)
    ax.set_xlim(-Lx/15,Lx/15)
    ax.set_ylim(-5,z_plot.max())
    ax.set_xlabel("x (Mm)",fontsize=20)
    ax.xaxis.set_major_locator(MaxNLocator(4,prune="both"))
    ax.yaxis.set_major_locator(MaxNLocator(5,prune="both"))
    ax.tick_params(axis='both', labelsize=18)
    plt.setp(ax.get_yticklabels(),visible=False)

add_subplot((0,0),true_vx,title="True $v_x$")
add_subplot((1,0),true_vz,title="True $v_z$")

#######################################################################################################

itercutoff = [15,35,35,35]


for strategy_index,datadir in enumerate(strategydir):

    current_vx = fitsread(os.path.join(strategydir[strategy_index],'update','vx_'+str(itercutoff[strategy_index])+'.fits'))
    current_vz = fitsread(os.path.join(strategydir[strategy_index],'update','vz_'+str(itercutoff[strategy_index])+'.fits'))

    add_subplot((0,strategy_index+1),current_vx,title="Iterated $v_x$ \#$"+str(strategy_index+1)+"$")
    add_subplot((1,strategy_index+1),current_vz,title="Iterated $v_z$ \#$"+str(strategy_index+1)+"$")


axlist[0].set_ylabel("Depth (Mm)",fontsize=22)
axlist[1].set_ylabel("Depth (Mm)",fontsize=22)
plt.setp(axlist[0].get_yticklabels(),visible=True)
plt.setp(axlist[1].get_yticklabels(),visible=True)
    
for cb in cblist: 
    cb.ax.set_ylabel("$\mathrm{m}/\mathrm{s}$",rotation=90,fontsize=16)
    cb.ax.tick_params(axis="x",labelsize=16)

#################################################################################

fig=plt.gcf()
fig.set_size_inches(14,9)
plt.tight_layout()
plt.subplots_adjust(wspace=0.)

save = read_params.parse_cmd_line_params("save")
if save is not None:
    savepath = os.path.join("plots",save)
    print "saving to",savepath
    if not os.path.exists("plots"): os.makedirs("plots")
    plt.savefig(savepath)
else:
    print "Not saving plot to file"

plt.show()
