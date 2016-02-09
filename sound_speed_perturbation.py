from __future__ import division
import numpy as np
import read_params
import matplotlib.pyplot as plt
import pyfits,os

Rsun=695.8
z=np.loadtxt(read_params.get_solarmodel(),usecols=[0]); z=(z-1)*Rsun
sigma_x = 30*np.exp(z/10)+20
nz = len(z)
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

delta_c = 0.5*h*v # km/s

plt.pcolormesh(x,z,h*v)
plt.xlim(x[0],x[-1])
plt.ylim(-10,z[-1])

datadir = read_params.get_directory()
pyfits.writeto(os.path.join(datadir,"model_c_ls00.fits"),delta_c,clobber=True)
#~ plt.show()
