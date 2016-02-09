from __future__ import division
import numpy as np
import pyfits
import matplotlib.pyplot as plt
import read_params
import os
from matplotlib.ticker import MaxNLocator,MultipleLocator
from matplotlib import rc,gridspec
import warnings
import plotc

Rsun = 695.8

def fitsread(f): return np.squeeze(pyfits.getdata(f))

datadir = read_params.get_directory()
src=1
srcloc = np.loadtxt(os.path.join(datadir,"master.pixels"),ndmin=1)[src-1]
data = fitsread(os.path.join(datadir,"tt","data","data"+str(src).zfill(2)+".fits"))

data_spec= np.fft.fftshift(abs(np.fft.fft2(data))**2)

nt,nx = data.shape
Lx = read_params.get_xlength()
dt = read_params.get_dt()
dx=Lx/nx

t_edges = (np.arange(nt+1)*dt-dt/2)/3600
x_edges = np.linspace(-Lx/2-dx/2,Lx/2-dx/2,nx+1)-srcloc

k = np.fft.fftshift(np.fft.fftfreq(nx,Lx/nx)*2*np.pi)
dk = k[1]-k[0]
k_edges = np.linspace(k[0]-dk/2,k[-1]+dk/2,nx+1)

nu = np.fft.fftshift(np.fft.fftfreq(nt,dt))*1e3
dnu = nu[1]-nu[0]
nu_edges = np.linspace(nu[0]-dnu/2,nu[-1]+dnu/2,nt+1)

rc("text",usetex=True)
rc('font',**{'family':'serif','serif':['Helvetica']})

if not os.path.exists("plots"): os.makedirs("plots")
###########################################################################################
#~ Time distance
fig=plt.figure()

datamax = abs(data).max()/50
plt.pcolormesh(x_edges,t_edges,data,cmap="Greys",vmax=datamax,vmin=-datamax,rasterized=True)
plt.xlim(x_edges[0],x_edges[-1])
plt.ylim(t_edges[0],t_edges[-1])
plt.xlabel("x (Mm)")
plt.ylabel("Time (hours)")

ax=plt.gca()
ax.xaxis.set_major_locator(MaxNLocator(5))
ax.yaxis.set_major_locator(MultipleLocator(1))

plotc.apj_1col_format(fig)
plt.tight_layout()

plt.savefig("plots/f2.eps")

#############################################################################################
#~ Spectrum

fig=plt.figure()
#~ gs=gridspec.GridSpec(1, 2, width_ratios=[3, 1]) 
#~ gs.update(wspace=0)
#~ ax0 = plt.subplot(gs[0])
plt.pcolormesh(k_edges*Rsun,nu_edges,data_spec,cmap='Greys',vmax=data_spec.max()*0.2,rasterized=True)
ax = plt.gca()
#~ cb=plt.colorbar()

plt.xlim(0,1000)
plt.ylim(0,8)

#~ ax.set_aspect(1000/8)

ax.xaxis.set_major_locator(MaxNLocator(5,prune="lower"))
ax.yaxis.set_major_locator(MaxNLocator(5,prune="both"))
#~ cb.ax.yaxis.set_major_locator(MaxNLocator(5))

plt.xlabel(r"$kR_\odot$")
plt.ylabel(r"Frequency (mHz)")

plt.text(600,2.2,r"f",bbox={"facecolor":"white","edgecolor":"black","pad":10},fontsize=14)
plt.text(620,3.3,r"p$_1$",bbox={"facecolor":"white","edgecolor":"black","pad":6},fontsize=16)
plt.text(700,4.2,r"p$_2$",bbox={"facecolor":"white","edgecolor":"black","pad":6},fontsize=16)
plt.text(620,4.8,r"p$_3$",bbox={"facecolor":"white","edgecolor":"black","pad":6},fontsize=16)
plt.text(547,5.1,r"p$_4$",bbox={"facecolor":"white","edgecolor":"black","pad":6},fontsize=16)
plt.text(480,5.3,r"p$_5$",bbox={"facecolor":"white","edgecolor":"black","pad":6},fontsize=16)
plt.text(400,5.4,r"p$_6$",bbox={"facecolor":"white","edgecolor":"black","pad":6},fontsize=16)

#~ ax1 = plt.subplot(gs[1])
#~ nu_src = 
#~ srcamp = 
#~ sigma = 
#~ soirce_spectrum = srcamp*np.exp(-(nu-nu_src)**2/(2*sigma**2))
#~ ax1.plot(srcamp,nu,color='black')

plotc.apj_1col_format(fig)
plt.tight_layout()
plt.savefig("plots/f3.eps",rasterized=True)
#~ plt.show()




