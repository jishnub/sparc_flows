from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import rc
import plotc
from plotc import texfonts
import read_params
import pyfits
import os

def fitsread(f): return np.squeeze(pyfits.getdata(f))

datadir=read_params.get_directory()

Rsun=695.8

z=np.loadtxt(read_params.get_solarmodel(),usecols=[0])
z=(z-1)*Rsun

Lx=read_params.get_xlength()
nx = read_params.get_nx()
x=np.linspace(-Lx/2,Lx/2,num=nx,endpoint=False)

vx=fitsread('true_vx.fits')
vz=fitsread('true_vz.fits')

ax=plotc.quiver2D(vx,vz,x=x,y=z,every=[2,2],xr=[-40,40],yr=[-5,None],
scale=8000,key=True,key_properties={'suffix':' m/s','fmt':'{:2.0f}','scale':250},rasterized=False,color='#777777',
usetex=True)

plt.xlabel("Horizontal Distance (Mm)",fontsize=20)
plt.ylabel("Depth (Mm)",fontsize=20)

excite_depth=read_params.get_excitedepth()/1e8
obs_depth = read_params.get_obs_depth()/1e8

src=1
src_locs=np.loadtxt(os.path.join(datadir,'master.pixels'),ndmin=1)
src_x = src_locs[src-1]
plt.plot([src_x],[excite_depth],'bo',markersize=12)

wavefronts=3
wavefront_r=[1.5+1.*i for i in xrange(wavefronts)]
pts=50
wfcoords=np.zeros((wavefronts,2,pts))
thetas=np.linspace(np.pi,2*np.pi,pts)
col=['red','blue']
for wfno in xrange(wavefronts):
    wfcoords[wfno,0]=wavefront_r[wfno]*np.cos(thetas)+src_x
    wfcoords[wfno,1]=wavefront_r[wfno]*np.sin(thetas)+excite_depth
    plt.plot(wfcoords[wfno,0],wfcoords[wfno,1],color=col[wfno%2],linewidth=1.5)


rec_x=18
plt.plot([rec_x],[obs_depth],'go',markersize=12)

ax.text(src_x,0.2,"Source",horizontalalignment='center',verticalalignment='bottom',
        fontsize=16)
        
ax.text(rec_x,0.35,"Receiver",horizontalalignment='center',verticalalignment='bottom',
        fontsize=16)

plt.tick_params(axis='both', which='major', labelsize=14)

if not os.path.exists("plots"): os.makedirs("plots")
rc('text', usetex=True)
plt.savefig("plots/flow_src_rec.eps")

plt.show()

