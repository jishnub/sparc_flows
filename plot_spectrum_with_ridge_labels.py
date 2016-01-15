from __future__ import division
import numpy as np
import pyfits
import matplotlib.pyplot as plt
import read_params
import os
from matplotlib.ticker import MaxNLocator

Rsun = 695.8

def fitsread(f): return np.squeeze(pyfits.getdata(f))

directory = read_params.get_directory()
data = fitsread(os.path.join(directory,"tt","data","data01.fits"))

data_spec= np.fft.fftshift(abs(np.fft.fft2(data))**2)

nt,nx = data.shape
Lx = read_params.get_xlength()
dt = read_params.get_dt()

k = np.fft.fftshift(np.fft.fftfreq(nx,Lx/nx)*2*np.pi)
dk = k[1]-k[0]
k_edges = np.linspace(k[0]-dk/2,k[-1]+dk/2,nx+1)

nu = np.fft.fftshift(np.fft.fftfreq(nt,dt))*1e3
dnu = nu[1]-nu[0]
nu_edges = np.linspace(nu[0]-dnu/2,nu[-1]+dnu/2,nt+1)

plt.pcolormesh(k_edges*Rsun,nu_edges,data_spec,cmap='Oranges',vmax=data_spec.max()*0.2,rasterized=True)
ax = plt.gca()
#~ cb=plt.colorbar()

plt.xlim(0,1000)
plt.ylim(0,8)

#~ ax.set_aspect(1000/8)

plt.tick_params(axis="both",which="major",labelsize=14)

ax.xaxis.set_major_locator(MaxNLocator(5))
ax.yaxis.set_major_locator(MaxNLocator(5,prune="lower"))
#~ cb.ax.yaxis.set_major_locator(MaxNLocator(5))

plt.xlabel(r"$kR_\odot$",fontsize=20)
plt.ylabel(r"$\nu$ (mHz)",fontsize=20)

plt.text(600,2.2,r"f",bbox={"facecolor":"white","edgecolor":"black","pad":10},fontsize=14)
plt.text(605,3.3,r"p$_1$",bbox={"facecolor":"white","edgecolor":"black","pad":6},fontsize=14)
plt.text(615,4.1,r"p$_2$",bbox={"facecolor":"white","edgecolor":"black","pad":6},fontsize=14)
plt.text(620,4.8,r"p$_3$",bbox={"facecolor":"white","edgecolor":"black","pad":6},fontsize=14)
plt.text(547,5.1,r"p$_4$",bbox={"facecolor":"white","edgecolor":"black","pad":6},fontsize=14)
plt.text(480,5.3,r"p$_5$",bbox={"facecolor":"white","edgecolor":"black","pad":6},fontsize=14)
plt.text(414,5.4,r"p$_6$",bbox={"facecolor":"white","edgecolor":"black","pad":6},fontsize=14)

plt.tight_layout()
plt.savefig("power_spectrum.eps",rasterized=True)
plt.show()




