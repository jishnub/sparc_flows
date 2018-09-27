
import numpy as np
import pyfits
import matplotlib.pyplot as plt
import read_params
import os
from matplotlib.ticker import MaxNLocator,MultipleLocator,ScalarFormatter
from matplotlib import rc
import warnings
import plotc

Rsun = 695.8

def fitsread(f): return np.squeeze(pyfits.getdata(f))

datadir = read_params.get_directory()
src=1
srcloc = np.loadtxt(os.path.join(datadir,"master.pixels"),ndmin=1)[src-1]
data = fitsread(os.path.join(datadir,"data",str(src).zfill(2)+".fits"))

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
rc('font',**{'family':'serif','serif':["Times"]})


tdfig=None
specfig=None

save = read_params.parse_cmd_line_params("save",return_list=True)
if save is not None: tdfig,specfig = save[0],save[1]

###########################################################################################
#~ Time distance
fig = plt.figure()
datamax = abs(data).max()/40
plt.pcolormesh(x_edges,t_edges,data,cmap="Greys",vmax=datamax,vmin=-datamax,rasterized=True)
cb=plt.colorbar(format='%.1e')
#~ cb.ax.yaxis.get_offset_text().set_position((3.5,0))
plt.xlim(x_edges[0],x_edges[-1])
plt.ylim(t_edges[0],t_edges[-1])
plt.xlabel("x (Mm)",fontsize=20)
plt.ylabel("Time (hours)",fontsize=20)

ax=plt.gca()
ax.xaxis.set_major_locator(MaxNLocator(5))
ax.yaxis.set_major_locator(MultipleLocator(1))

plt.tick_params(axis="both",labelsize=18)

fig.set_size_inches(6,4)
plt.tight_layout()

if tdfig is not None:
    savepath = os.path.join("plots",save[0])
    print("saving time distance plot to",savepath)
    if not os.path.exists("plots"): os.makedirs("plots")
    plt.savefig(savepath)
else:
    print("Not saving time distance plot to file")


#############################################################################################
#~ Spectrum
mode = read_params.parse_cmd_line_params(key="mode",default=None)
if mode is not None:
    modefilter = np.fft.fftshift(np.squeeze(pyfits.getdata(os.path.join(
                '{}mode_filter.fits'.format(mode)))))
    data_spec*=modefilter
    # nu_ind_max,k_ind_max=np.unravel_index(data_spec.argmax(),data_spec.shape)

    nu_ind_sigma,k_ind_sigma=np.unravel_index(
            abs(data_spec-data_spec.max()*np.exp(-0.5)).argmin(),
            data_spec.shape)

    power = np.zeros_like(k)
    nu_ridge = np.zeros_like(nu)
    for l_ind,l_row in enumerate(data_spec.T):
        power[l_ind]=l_row.max()
        nu_ridge[l_ind] = nu[l_row.argmax()]

    from scipy.optimize import curve_fit

    def gaus(x,a,x0,sigma):
        return a*np.exp(-(x-x0)**2/(2*sigma**2))

    popt,pcov = curve_fit(gaus,k[(k*Rsun>200) & (k*Rsun<800)]*Rsun,
        power[(k*Rsun>200) & (k*Rsun<800)],
        p0=[2,500,100])


    plt.figure()
    plt.plot(k[k*Rsun>0]*Rsun,power[k*Rsun>0],color="dodgerblue",lw=2)
    plt.plot(k[k*Rsun>0]*Rsun,gaus(k[k*Rsun>0]*Rsun,*popt),color="black",lw=2)
    plt.axhline(popt[0]*np.exp(-0.5),ls='dashed',color='red')
    plt.axhline(popt[0]*0.5,ls='dashed',color='olive')

    print(("max",int(popt[1]),abs(nu_ridge[abs(k*Rsun-popt[1]).argmin()])))
    print(("1 sigma",int(popt[1]-popt[2]),abs(nu_ridge[abs(k*Rsun-popt[1]+popt[2]).argmin()])))
    print(("half power",int(popt[1]-popt[2]*np.sqrt(np.log(4))),
    abs(nu_ridge[abs(k*Rsun-popt[1]+popt[2]*np.sqrt(np.log(4))).argmin()])))

fig = plt.figure()
plt.pcolormesh(k_edges*Rsun,nu_edges,data_spec/data_spec.max(),
cmap='Greys',vmax=1,rasterized=True)
ax = plt.gca()
plt.colorbar()

plt.xlim(0,1000)
plt.ylim(0,8)

ax.xaxis.set_major_locator(MaxNLocator(5,prune="lower"))
ax.yaxis.set_major_locator(MaxNLocator(5,prune="both"))

plt.xlabel(r"$\mathrm{k\,R_\odot}$",fontsize=20,labelpad=5)
plt.ylabel(r"Frequency (mHz)",fontsize=20)

plt.text(943,3.1,r"f",     fontsize=20)
plt.text(943,4.35,r"p$_1$",fontsize=20)
plt.text(860,5.05,r"p$_2$",fontsize=20)
plt.text(750,5.4,r"p$_3$", fontsize=20)
plt.text(667,5.67,r"p$_4$",fontsize=20)
plt.text(587,5.8,r"p$_5$", fontsize=20)
plt.text(521,6.0,r"p$_6$", fontsize=20)
plt.text(460,6.2,r"p$_7$", fontsize=20)

plt.tick_params(axis="both",labelsize=18)
fig.set_size_inches(6,4)
plt.tight_layout()

if specfig is not None:
    savepath = os.path.join("plots",specfig)
    print("saving spectrum plot to",savepath)
    if not os.path.exists("plots"): os.makedirs("plots")
    plt.savefig(savepath)
else:
    print("Not saving spectrum plot to file")

plt.show()
