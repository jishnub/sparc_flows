
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import rc
import plotc
import read_params
import pyfits
import os


def fitsread(f): return np.squeeze(pyfits.getdata(f))

rc("text",usetex=True)
rc('font',**{'family':'serif'})

datadir=read_params.get_directory()

Rsun=695.8

z=np.loadtxt(read_params.get_solarmodel(),usecols=[0])
z=(z-1)*Rsun

Lx=read_params.get_xlength()
nx = read_params.get_nx()
x=np.linspace(-Lx/2,Lx/2,num=nx,endpoint=False)

vx=fitsread('true_vx.fits')
vz=fitsread('true_vz.fits')

qv=plotc.quiver2D(vx,vz,x=x,y=z,every=[2,2],xr=[-80,80],yr=[-5,None],
scale=8000,key=False,key_properties={'suffix':' m/s','fmt':'{:2.0f}','scale':250},
rasterized=False,color='#777777')

ax=qv.axis

plt.xlabel("Horizontal Distance (Mm)",fontsize=20)
plt.ylabel("Depth (Mm)",fontsize=20)

excite_depth=read_params.get_excitedepth()/1e8
obs_depth = read_params.get_obs_depth()/1e8

src=1
src_locs=np.loadtxt(os.path.join(datadir,'master.pixels'),ndmin=1)
src_x = src_locs[src-1]
plt.plot([src_x],[excite_depth],marker='o',markerfacecolor='#666666',markersize=12)

wavefronts=3
wavefront_r=[1.5+1.*i for i in range(wavefronts)]
pts=50
wfcoords=np.zeros((wavefronts,2,pts))
thetas=np.linspace(np.pi,2*np.pi,pts)
for wfno in range(wavefronts):
    wfcoords[wfno,0]=wavefront_r[wfno]*np.cos(thetas)+src_x
    wfcoords[wfno,1]=wavefront_r[wfno]*np.sin(thetas)+excite_depth
    plt.plot(wfcoords[wfno,0],wfcoords[wfno,1],color='black',linewidth=1.5)


rec_x=18
plt.plot([rec_x],[obs_depth],marker='o',markerfacecolor='#333333',markersize=12)

ax.text(src_x,0.2,"Source",horizontalalignment='center',verticalalignment='bottom',
        fontsize=20)

ax.text(rec_x,0.4,"Receiver",horizontalalignment='center',verticalalignment='bottom',
        fontsize=20)


if not os.path.exists("plots"): os.makedirs("plots")

plt.gcf().set_size_inches(6,4)
plt.tick_params(axis="both",labelsize=18)
plt.tight_layout()

save = read_params.parse_cmd_line_params("save")
if save is not None:
    savepath = os.path.join("plots",save)
    print("saving to",savepath)
    if not os.path.exists("plots"): os.makedirs("plots")
    plt.savefig(savepath)
else:
    print("Not saving plot to file")

plt.show()
