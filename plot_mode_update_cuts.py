from __future__ import division
import matplotlib
from matplotlib import rcParams
import matplotlib.pyplot as plt
import numpy as np
import read_params
import pyfits,os
from matplotlib.ticker import MaxNLocator,NullFormatter,ScalarFormatter
import plotc

def fitsread(f): return np.squeeze(pyfits.getdata(f))

Lx = read_params.get_xlength()
nx = read_params.get_nx()
x = np.linspace(-Lx/2,Lx/2,nx,endpoint=False)
z=np.loadtxt(read_params.get_solarmodel(),usecols=[0]);z=(z-1)*695.8

parentdir = os.path.dirname(read_params.get_directory())

arrays_to_plot = []
arrays_to_plot_max_x = []

modes_list = ['f']
for i in xrange(1,8): modes_list.append("p"+str(i))
modes_list.append('first_bounce_p')

modes_found = []

for mode in modes_list: 
    try:
        arrays_to_plot.append(fitsread(os.path.join(parentdir,mode,'update','update_psi_00.fits')))
        modes_found.append(mode)
    except IOError:
        pass

if not modes_found: 
    print "No modes found"
    exit()

arrays_to_plot.append(fitsread('true_psi.fits'))
arrays_to_plot_max_x=divmod(arrays_to_plot[-1].argmax(),arrays_to_plot[-1].shape[1])

modes_found.append(r"$\psi^{ref}$")

nplots=min(9,len(modes_found)+1) # may not want to plot all modes

depth_cutoff = -10

axeslist=[]

for ind,mode in enumerate(modes_found):
    axeslist.append(plt.subplot(1,nplots,ind+1))
    kernel = arrays_to_plot[ind][:,arrays_to_plot_max_x[1]]
    #~ kernel /= kernel.max()
    plt.plot(kernel,z,color='black',linewidth=2 if ind==nplots-1 else 1)
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
plt.xlabel("Sensitivity Kernel ($s^2/$Mm)",labelpad=20)

plt.subplot(1,nplots,nplots)
plt.xlabel(r"$\psi$ (Mm)",labelpad=20)

plotc.apj_2col_format(plt.gcf(),default_fontsize=14)
plt.gcf().set_size_inches(8,5.3)
plt.tight_layout()
plt.subplots_adjust(wspace=0)

ax_index = None

def on_axis_enter(event):
    global ax_index
    ax_index = axeslist.index(event.inaxes)

def on_axis_leave(event):
    global ax_index
    ax_index = None
    
def on_click(event):
    global ax_index
    if ax_index is not None:
        plt.figure()
        plotc.colorplot(arrays_to_plot[ax_index],x=x,y=z,yr=[depth_cutoff,z.max()])
        plt.show(block=False)

plt.connect('axes_enter_event', on_axis_enter)
plt.connect('axes_leave_event', on_axis_leave)
plt.connect('button_press_event', on_click)

if not os.path.exists("plots"): os.makedirs("plots")
plt.savefig("plots/f4.eps")
plt.show()
