from __future__ import division
import numpy as np
import pyfits
import matplotlib.pyplot as plt
import read_params
import os
from matplotlib.ticker import MaxNLocator,MultipleLocator
from matplotlib import rc
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

#~ rc("text",usetex=True)
rc('font',**{'family':'serif'})


tdfig=None
specfig=None

save = read_params.parse_cmd_line_params("save",return_list=True)
if save is not None: tdfig,specfig = save[0],save[1]

###########################################################################################
#~ Time distance

datamax = abs(data).max()/50
plt.pcolormesh(x_edges,t_edges,data,cmap="Greys",vmax=datamax,vmin=-datamax,rasterized=True)
plt.xlim(x_edges[0],x_edges[-1])
plt.ylim(t_edges[0],t_edges[-1])
plt.xlabel("x (Mm)")
plt.ylabel("Time (hours)")

ax=plt.gca()
ax.xaxis.set_major_locator(MaxNLocator(5))
ax.yaxis.set_major_locator(MultipleLocator(1))

plotc.apj_1col_format(plt.gcf())
plt.tight_layout()

if tdfig is not None:
    savepath = os.path.join("plots",save[0])
    print "saving time distance plot to",savepath
    if not os.path.exists("plots"): os.makedirs("plots")
    plt.savefig(savepath)
else:
    print "Not saving time distance plot to file"


#############################################################################################
#~ Spectrum

fig = plt.figure()
plt.pcolormesh(k_edges*Rsun,nu_edges,data_spec,cmap='Greys',vmax=data_spec.max()*0.2,rasterized=True)
ax = plt.gca()

plt.xlim(0,1000)
plt.ylim(0,8)

ax.xaxis.set_major_locator(MaxNLocator(5,prune="lower"))
ax.yaxis.set_major_locator(MaxNLocator(5,prune="both"))

plt.xlabel(r"$kR_\odot$")
plt.ylabel(r"Frequency (mHz)")

plt.text(943,3.1,r"f",     fontsize=20)
plt.text(943,4.35,r"p$_1$",fontsize=20)
plt.text(860,5.05,r"p$_2$",fontsize=20)
plt.text(750,5.4,r"p$_3$", fontsize=20)
plt.text(667,5.67,r"p$_4$",fontsize=20)
plt.text(587,5.8,r"p$_5$", fontsize=20)
plt.text(521,6.0,r"p$_6$", fontsize=20)
plt.text(460,6.2,r"p$_7$", fontsize=20)

plt.tight_layout()

if specfig is not None:
    savepath = os.path.join("plots",specfig)
    print "saving spectrum plot to",savepath
    if not os.path.exists("plots"): os.makedirs("plots")
    plt.savefig(savepath)
else:
    print "Not saving spectrum plot to file"

plt.show()




