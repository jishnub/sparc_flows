import numpy as np
from scipy import interpolate,integrate,fftpack
from astropy.io import fits
import read_params
from scipy.special import j1,j0,jn
def j2(z): return jn(2,z)
def j1prime(z): return 0.5*(j0(z)-j2(z))
import os,shutil
import dbyd2
from pathlib import Path

Lx = read_params.get_xlength()
nx = read_params.get_nx()
x = np.linspace(-Lx/2,Lx/2,nx,endpoint=False)
Rsun=6.95989467700E2
z,c_sound,rho = np.loadtxt(read_params.get_solarmodel(),usecols=[0,1,2],unpack=True); z=(z-1)*Rsun;c_sound/=100;


# ### Duvall-Hanasoge model

kDH13 = 2*np.pi/30
RDH13 = 15
sigmazDH13 = 3
z0DH13 = -9
v0DH13 = 1000


# ## Spline only along z


z_cutoff = -25 # Spline lower cutoff
zspline_ind = z>z_cutoff
zspline = z[z>z_cutoff]

psi_z_profile = v0DH13/c_sound/kDH13*np.exp(-(z-z0DH13)**2/(2*sigmazDH13**2))


def coeff_to_model(tck_z):
    f0_x = np.sign(x)*j1(kDH13*abs(x))*np.exp(-abs(x)/RDH13)
    h_z=interpolate.splev(z,tck_z,ext=1)
    return f0_x[None,:]*h_z[:,None]


s_misfit = []
s_list = [1e-1,1e-2,1e-3,1e-4,1e-5,1e-6,1e-7,1e-8]
Nparams = []
for s in s_list:
    tz1D,cz1D,kz=interpolate.splrep(zspline,psi_z_profile[z>z_cutoff],k=2,s=s)
    misfit = (integrate.simps((psi_z_profile - interpolate.splev(z,(tz1D,cz1D,kz),ext=1))**2,x=z)/
            integrate.simps(psi_z_profile**2,x=z))
    s_misfit.append(misfit*100)
    Nparams.append(tz1D.size-kz-1)

f=interpolate.interp1d(s_misfit,s_list)
smoothing_par = f(0.005).item()
np.set_printoptions(precision=3)
print("Smoothing param at 0.1% level {:.1E}".format(smoothing_par))

tz1D,cz1D,kz=interpolate.splrep(zspline,psi_z_profile[z>z_cutoff],k=2,s=smoothing_par)
print("Number of knots:",tz1D.size)
print("Number of coeffs:",tz1D.size-kz-1)

print(', '.join("{:.1f}".format(knot) for knot in tz1D))


# ### Determine the coefficient that corresponds to the surface


b_i_surf = np.zeros_like(cz1D)
deep_z_cutoff = z_cutoff + 2
b_i_deep = np.zeros_like(cz1D)

for i in range(tz1D.size):
    c_i = np.zeros_like(cz1D)
    c_i[i] = 1
    b_i_surf[i] = interpolate.splev(0,(tz1D,c_i,kz))
    b_i_deep[i] = interpolate.splev(deep_z_cutoff,(tz1D,c_i,kz))

c_surf_cutoff = b_i_surf.argmax()
print("c surf cutoff index",c_surf_cutoff)
cz1D_top = np.zeros_like(cz1D)
cz1D_top[c_surf_cutoff:] = cz1D[c_surf_cutoff:]
cz1D_bottom = np.zeros_like(cz1D)
cz1D_bottom[:c_surf_cutoff] = cz1D[:c_surf_cutoff]
print("{:d} parameters to  be fit along z, {:d} parameters clamped, "
      .format(c_surf_cutoff,tz1D.size - c_surf_cutoff))
print("Coeffs below surface",cz1D[:c_surf_cutoff])

deep_z_cutoff_ind = b_i_deep.argmax()
print("Deep z cutoff ind",deep_z_cutoff_ind)


psi_spl_fit_top=interpolate.splev(z,(tz1D,cz1D_top,kz),ext=1);
psi_spl_fit_bottom=interpolate.splev(z,(tz1D,cz1D_bottom,kz),ext=1);



psi_spline = coeff_to_model((tz1D,cz1D_top+cz1D_bottom,kz))
psi_spline_surf = coeff_to_model((tz1D,cz1D_top,kz))

vx_spline = -dbyd2.dbyd2(rho*c_sound*psi_spline.T,1)/dbyd2.dbyd2(np.atleast_2d(z),1)/rho
vx_spline = vx_spline.T
vz_spline = np.zeros_like(psi_spline)
for rowno,row in enumerate(psi_spline):
    vz_spline[rowno] = fftpack.diff(row,period=Lx)*c_sound[rowno]



def peak_and_surf_vel(v_arr,label="v"):
    v_peak_z,v_peak_x = np.unravel_index(v_arr.argmax(),v_arr.shape)
    print("{} peak velocity at z={:.1f} Mm, magnitude {:.1f} m/s".format(label,z[v_peak_z],v_arr[v_peak_z,v_peak_x]))
    print("{} peak surface velocity {:.1f} m/s".format(label,v_arr[abs(z).argmin(),v_peak_x]))

peak_and_surf_vel(vx_spline,label="vx")
peak_and_surf_vel(vz_spline,label="vz")



import read_params

datadir = Path(read_params.get_directory())
    
fits.writeto(datadir/"model_psi_ls00.fits",psi_spline_surf[:,np.newaxis,:],overwrite=True)

shutil.copyfile(datadir/"model_psi_ls00.fits", datadir/"model_psi_ls00_start.fits")

np.savez(datadir/"model_psi_ls00_coeffs.npz",z=np.zeros_like(cz1D))

shutil.copyfile(datadir/"model_psi_ls00_coeffs.npz", datadir/"model_psi_ls00_coeffs_start.npz")

np.savez(datadir/"true_psi_coeffs.npz",tz=tz1D,kz=kz,
    cz_top=cz1D_top,cz_bot=cz1D_bottom,c_surf_cutoff=c_surf_cutoff,z_cutoff = z_cutoff,
        deep_z_cutoff_ind=deep_z_cutoff_ind)

fits.writeto("true_psi.fits",psi_spline,overwrite=True)
fits.writeto("true_vx.fits",vx_spline,overwrite=True)
fits.writeto("true_vz.fits",vz_spline,overwrite=True)

