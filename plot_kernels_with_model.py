
import matplotlib
from matplotlib import rcParams
import matplotlib.pyplot as plt
import numpy as np
import read_params
import pyfits,os
from matplotlib import rc
from matplotlib.ticker import MaxNLocator,NullFormatter,ScalarFormatter
from scipy.ndimage.filters import gaussian_filter1d

rc("text",usetex=True)
rc("font",**dict(family="serif",serif="Times"))

def fitsread(f): return np.squeeze(pyfits.getdata(f))

Lx = read_params.get_xlength()
nx = read_params.get_nx()
x = np.linspace(-Lx/2,Lx/2,nx,endpoint=False)
z=np.loadtxt(read_params.get_solarmodel(),usecols=[0]);z=(z-1)*695.8

datadir = read_params.get_directory()
parentdir = os.path.dirname(datadir)

strategydir=[os.path.join(parentdir,modeldir) for modeldir in ("f_p1","f_to_p3","f_to_p7_2pixsmooth","f_to_p7_new")]

arrays_to_plot = []
arrays_to_plot_max_x = []

modes_list = ['f']
for i in range(1,8): modes_list.append("p"+str(i))
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
    print("No modes found")
    exit()

true_psi=fitsread('true_psi.fits')
psi_max_x_pix=np.unravel_index(true_psi.argmax(),true_psi.shape)[1]
true_psi = true_psi[:,psi_max_x_pix]

iterated_psi_1 = fitsread(os.path.join(strategydir[0],"model_psi_ls00.fits"))
iterated_psi_1-=iterated_psi_1[0,0]
iterated_psi_1 = iterated_psi_1[:,psi_max_x_pix]

iterated_psi_2 = fitsread(os.path.join(strategydir[1],"model_psi_ls00.fits"))
iterated_psi_2-=iterated_psi_2[0,0]
iterated_psi_2 = iterated_psi_2[:,psi_max_x_pix]

iterated_psi_3 = fitsread(os.path.join(strategydir[2],"model_psi_ls00.fits"))
iterated_psi_3-=iterated_psi_3[0,0]
iterated_psi_3 = iterated_psi_3[:,psi_max_x_pix]

depth_cutoff = -8

nplots = len(modes_found)+1

rho=np.loadtxt(read_params.get_solarmodel(),usecols=[2])

for ind,mode in enumerate(modes_found):
    kernel = arrays_to_plot[ind]
    plt.subplot(1,nplots,ind+1)
    plt.plot(kernel,z,color='black')
    plt.plot([0]*len(z),z,linestyle='dotted',color='black')
    plt.ylim(depth_cutoff,z.max())
    ax=plt.gca()
    ax.xaxis.set_major_locator(MaxNLocator(3,prune='upper'))
    #~ ax.xaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.xaxis.set_major_formatter(ScalarFormatter())
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    plt.tick_params(axis="both",labelsize=15)
    plt.title(mode)
    ax.xaxis.get_offset_text().set_size(13)
    if ind>0: ax.yaxis.set_major_formatter(NullFormatter())        
    ax.grid()

plt.subplot(1,nplots,1)
plt.ylabel("Depth $(\mathrm{Mm})$",fontsize=18)

plt.subplot(1,nplots,nplots//2+1)
plt.xlabel("Sensitivity Kernel $(s^2\,\mathrm{Mm}^{-3})$",labelpad=20,fontsize=18)

ax=plt.subplot(1,nplots,nplots)
plt.plot(true_psi,z,linewidth=2,label="True",color="#777777")
plt.xlabel(r"$\psi\;(\mathrm{Mm})$",labelpad=20,fontsize=18)
plt.ylim(depth_cutoff,z.max())

plt.tick_params(labelsize=15)
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
ax.xaxis.get_offset_text().set_size(13)
 


plt.plot(iterated_psi_1,z,linewidth=2,color="black",linestyle="dotted",label="$\#1$")
plt.plot(iterated_psi_2,z,linewidth=2,color="black",linestyle="dashdot",label="$\#2$")
plt.plot(iterated_psi_3,z,linewidth=2,color="black",linestyle="dashed",label="$\#3$")
plt.ylim(depth_cutoff,z.max())
ax.xaxis.set_major_locator(MaxNLocator(3,prune='upper'))
xaxis_formatter = ScalarFormatter()
xaxis_formatter.set_scientific(True)
xaxis_formatter.set_powerlimits((0,0))
ax.xaxis.set_major_formatter(xaxis_formatter)
ax.yaxis.set_major_formatter(NullFormatter())  

plt.legend(loc="lower left",fontsize=12)

plt.gcf().set_size_inches(12,6)
plt.tight_layout()
plt.subplots_adjust(wspace=0.2)


save = read_params.parse_cmd_line_params("save")
if save is not None:
    savepath = os.path.join("plots",save)
    print("saving to",savepath)
    if not os.path.exists("plots"): os.makedirs("plots")
    plt.savefig(savepath)
else:
    print("Not saving plot to file")

plt.show()
