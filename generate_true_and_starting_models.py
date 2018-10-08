import numpy as np
import read_params
from dbyd2 import dbyd2
from scipy.special import j1,j0,jn
from scipy import interpolate,fftpack
from astropy.io import fits
from astropy.constants import R_sun
from astropy import units
from pathlib import Path
import shutil

def j2(x): return jn(2,x)
def j1prime(x): return 0.5*(j0(x)-j2(x))


Lx=read_params.get_xlength() * units.Mm
nx=read_params.get_nx()
x=np.linspace(-Lx/2,Lx/2,nx,endpoint=False)

z,c,rho=np.loadtxt(read_params.get_solarmodel(),usecols=[0,1,2],unpack=True)
c *= units.cm/units.s
Rsun = R_sun.to("cm")
z=(z-1)*Rsun

dz = dbyd2(np.atleast_2d(z),1).T

def ddz(arr):
    # always returns 2D
    if len(arr.shape) == 1:
        return dbyd2(np.atleast_2d(arr),1).T / dz
    # assume f(z,x), derivative along first axis
    elif len(arr.shape) == 2:
        return dbyd2(arr,1) / dz
    
dzlnrho = ddz(np.log(rho)).squeeze() / units.Mm


v0 = 240 * units.m/units.s
R = 15 * units.Mm
k = 2*np.pi/(2*R)
z0 = -2.3 * units.Mm
sigmaz = 0.912 * units.Mm

u = v0/k*np.exp(-(z-z0)**2/(2*sigmaz**2))/c
u = u[:,np.newaxis]

hx = np.sign(x)*np.exp(-np.abs(x)/R)*j1(k*np.abs(x))
psi = u*hx

vx = (c*((z-z0)/sigmaz**2 - dzlnrho))[:,None]*psi
vz = c[:,None]*u*k*(j1prime(k*np.abs(x)) - 1/(k*R)*j1(k*np.abs(x)))*np.exp(-np.abs(x)/R)

fits.writeto("true_psi.fits",psi.to("Mm").value,overwrite=True)
fits.writeto("true_vx.fits",vx.to("m/s").value,overwrite=True)
fits.writeto("true_vz.fits",vz.to("m/s").value,overwrite=True)


datadir = Path(read_params.get_directory())
fits.writeto(datadir/"model_psi_ls00.fits",np.zeros_like(psi.value),overwrite=True)

shutil.copyfile(datadir/"model_psi_ls00.fits",datadir/"model_psi_ls00_start.fits")