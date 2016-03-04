from __future__ import division
import matplotlib
from matplotlib import rcParams
import matplotlib.pyplot as plt
import numpy as np
import read_params
import pyfits,os
from matplotlib.ticker import MaxNLocator,NullFormatter,ScalarFormatter
from scipy.ndimage.filters import gaussian_filter1d
import plotc

def fitsread(f): return np.squeeze(pyfits.getdata(f))

Lx = read_params.get_xlength()
nx = read_params.get_nx()
x = np.linspace(-Lx/2,Lx/2,nx,endpoint=False)
z=np.loadtxt(read_params.get_solarmodel(),usecols=[0]);z=(z-1)*695.8

datadir = read_params.get_directory()
parentdir = os.path.dirname(datadir)

arrays_to_plot = []
arrays_to_plot_max_x = []

modes_list = ['f']
for i in xrange(1,8): modes_list.append("p"+str(i))
modes_list.append('first_bounce_p')

modes_found = []

def smooth_z(arr):
    return gaussian_filter1d(arr,sigma=2,axis=0)
    
def smooth_x(arr):
    return gaussian_filter1d(arr,sigma=1,axis=1)

for mode in modes_list: 
    try:
        kernel_psi=fitsread(os.path.join(parentdir,mode,'kernel','kernel_psi_01.fits'))
        kernel_psi = smooth_x(kernel_psi)
        kernel_psi = smooth_z(kernel_psi)
        kernel_psi_mean = np.trapz(kernel_psi,x=x,axis=1)/Lx        
        arrays_to_plot.append(kernel_psi_mean)
        modes_found.append(mode)
    except IOError:
        pass

if not modes_found: 
    print "No modes found"
    exit()

Lregular = 30. # \psi is scaled by Lregular in Mm
true_psi=fitsread('true_psi.fits')*Lregular
psi_max_x_pix=np.unravel_index(true_psi.argmax(),true_psi.shape)[1]
true_psi = true_psi[:,psi_max_x_pix]

iterated_psi = fitsread(os.path.join(datadir,"model_psi_ls00.fits"))*Lregular
iterated_psi-=iterated_psi[0,0]
iterated_psi = iterated_psi[:,psi_max_x_pix]

depth_cutoff = -10

nplots = len(modes_found)+1

for ind,mode in enumerate(modes_found):
    kernel = arrays_to_plot[ind]
    plt.subplot(1,nplots,ind+1)
    plt.plot(kernel,z,color='black')
    plt.plot([0]*len(z),z,linestyle='dotted',color='black')
    plt.ylim(depth_cutoff,z.max())
    ax=plt.gca()
    ax.xaxis.set_major_locator(MaxNLocator(3,prune='upper'))
    #~ ax.xaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    plt.title(mode)
    ax.xaxis.get_offset_text().set_size(11)
    #~ plt.text(0.92, -0.07,sci_power, fontsize=12 , transform = ax.transAxes)
    if ind>0: ax.yaxis.set_major_formatter(NullFormatter())        
    ax.grid()

plt.subplot(1,nplots,1)
plt.ylabel("Depth (Mm)")

plt.subplot(1,nplots,nplots//2+1)
plt.xlabel("Sensitivity Kernel ($s^2\,\mathrm{Mm}^{-3}$)",labelpad=20)

ax=plt.subplot(1,nplots,nplots)
plt.plot(true_psi,z,linewidth=2,label="True",color="black")
plt.xlabel(r"$\psi$ (Mm)",labelpad=20,color="black")
plt.ylim(depth_cutoff,z.max())
ax.xaxis.set_major_locator(MaxNLocator(3,prune='upper'))
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
ax.xaxis.get_offset_text().set_size(11)
ax.yaxis.set_major_formatter(NullFormatter())   

plt.plot(iterated_psi,z,linewidth=2,color="black",linestyle="dashed",label="Iter")
plt.ylim(depth_cutoff,z.max())
ax.xaxis.set_major_locator(MaxNLocator(3,prune='upper'))
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
ax.xaxis.get_offset_text().set_size(11)
ax.yaxis.set_major_formatter(NullFormatter())
ax.grid()

plt.legend(loc="best",fontsize=14)

plotc.apj_2col_format(plt.gcf(),default_fontsize=14)
plt.gcf().set_size_inches(8,5.3)
plt.tight_layout()
plt.subplots_adjust(wspace=0.2)


save = read_params.parse_cmd_line_params("save")
if save is not None:
    savepath = os.path.join("plots",save)
    print "saving to",savepath
    if not os.path.exists("plots"): os.makedirs("plots")
    plt.savefig(savepath)
else:
    print "Not saving plot to file"

plt.show()
