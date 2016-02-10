from __future__ import division
import numpy as np
import read_params
import matplotlib.pyplot as plt
import pyfits,os

Rsun=695.8
z,c=np.loadtxt(read_params.get_solarmodel(),usecols=[0,1],unpack=True); z=(z-1)*Rsun
nz = len(z)
c=c.reshape(nz,1)
sigma_x = 30*np.exp(z/10)+20
sigma_x = sigma_x.reshape(nz,1)

Lx=read_params.get_xlength()
nx=read_params.get_nx()

x=np.linspace(-Lx/2,Lx/2,nx,endpoint=False)
x2D=x.reshape(1,nx)

h = np.exp(-x2D**2/(2*sigma_x**2))

z0 = -1.5 # Mm
sigma_z = 1 # Mm
v = np.exp(-(z-z0)**2/(2*sigma_z)**2)
v = v.reshape(nz,1)

delta_c = 0.5*h*v * 1e5 # cm/s
c_bump = c + delta_c

c_00 = np.tile(c,(1,nx))

#~ plt.pcolormesh(x,z,h*v)
#~ plt.xlim(x[0],x[-1])
#~ plt.ylim(-10,z[-1])

datadir = read_params.get_directory()
pyfits.writeto("true_c.fits",c_bump,clobber=True)
pyfits.writeto(os.path.join(datadir,"model_c_ls00.fits"),c_00,clobber=True)

#~ plt.show()
