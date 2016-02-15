from __future__ import division
import numpy as np
import read_params
import matplotlib.pyplot as plt
import pyfits,os,warnings

Rsun=695.8
z,c=np.loadtxt(read_params.get_solarmodel(),usecols=[0,1],unpack=True); z=(z-1)*Rsun
nz = len(z)
c=c.reshape(nz,1)

sigma_x = 15*np.exp(z/10)+10
sigma_x = sigma_x.reshape(nz,1)

Lx=read_params.get_xlength()
nx=read_params.get_nx()
print "nx={:d}, nz={:d}".format(nx,nz)

x=np.linspace(-Lx/2,Lx/2,nx,endpoint=False)
x2D=x.reshape(1,nx)

h = np.exp(-x2D**2/(2*sigma_x**2))

z0 = -2.5 # Mm
sigma_z = 1.4 # Mm
v = np.exp(-(z-z0)**2/(2*sigma_z)**2)
v = v.reshape(nz,1)

c_00 = np.tile(c,(1,nx))

delta_c = 3e-3*h*v *c_00

c_00 = c_00.reshape(nz,1,nx)
delta_c = delta_c.reshape(nz,1,nx)

#~ plt.pcolormesh(x,z,np.squeeze(delta_c/c_00))
#~ plt.title("delta ln c")
#~ plt.colorbar()
#~ plt.xlim(x[0]/2,x[-1]/2)
#~ plt.ylim(-10,z[-1])

datadir = read_params.get_directory()
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    pyfits.writeto(os.path.join(datadir,"model_c_ls00.fits"),c_00,clobber=True)

plt.show()
