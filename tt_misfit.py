from __future__ import division
import numpy as np
import plotc
import pyfits as pf
import os,re,sys
import matplotlib.patches as patches
from matplotlib.ticker import MaxNLocator
import read_params
import matplotlib.pyplot as plt

def fitsread(fitsfile): return np.squeeze(pf.getdata(fitsfile))

codedir=os.path.dirname(os.path.abspath(__file__))
datadir=read_params.get_directory()
masterpixelsfile=os.path.join(datadir,'master.pixels')
masterpixels=np.loadtxt(masterpixelsfile,ndmin=1)

dt=read_params.get_dt()/60.

temp=filter(lambda x: x.startswith("src=") or x.startswith("source="),sys.argv)
if temp: src=int(temp[0].split("=")[-1])
else: src=1

itercutoff = filter(lambda x: x.startswith("iter="),sys.argv)
if len(itercutoff)!=0: itercutoff=int(itercutoff[0].split("=")[-1])
else: itercutoff=np.inf

if src>len(masterpixels): src=len(masterpixels)
src=str(src).zfill(2)


nx=read_params.get_nx()
Lx=read_params.get_xlength()
x=np.linspace(-Lx/2,Lx/2,nx,endpoint=False)

#~ Read in velocity and compute length scales

vxfile=os.path.join(codedir,'true_vx.fits')
vx=fitsread(vxfile)
vx_max_row_index,vx_max_col_index = divmod(vx.argmax(),nx)
vx_max_row = vx[vx_max_row_index]
vx=None
vxmax=vx_max_row.max()
vx_hwhm=np.where(vx_max_row>vxmax/2)[0][-1]-nx//2
vlimleft,vlimright=x[nx//2-vx_hwhm],x[nx//2+vx_hwhm]

datafile=os.path.join(datadir,'forward_src'+src+'_ls00','data.fits')
data=fitsread(datafile)
nt=data.shape[0]
time=np.arange(nt)*dt



modes={'0':'fmode'}
for pmodeno in xrange(1,8): modes.update({str(pmodeno):'p'+str(pmodeno)+'mode'})
modes['8']="first_bounce_pmode"
ridge_filters=sorted(read_params.get_modes_used())

def spaced(a): 
    b=a[:-4]+" "+a[-4:]
    b=b.replace("_"," ")
    return b

vzcc=[]
vzcc_ridges=[]
iters=[]



srcloc=masterpixels[int(src)-1]
srcloc_ind=np.argmin(abs(x-srcloc))
srcloc=round(x[srcloc_ind],1)


ttdiff_ridges={}
for mode in modes: ttdiff_ridges.update({mode:{"iter_no":[],"misfits":[]}})
for iteration_no in xrange(0,100):
    ttpath=os.path.join(datadir,'tt','iter'+str(iteration_no).zfill(2))
    if os.path.exists(ttpath):
        for modeno,mode in enumerate(ridge_filters):
            ttfile=os.path.join(ttpath,'ttdiff_src'+src+'.'+modes[mode])
            if os.path.exists(ttfile):
                ttdiff_ridges[mode]["iter_no"].append(iteration_no)
                ttdiff_ridges[mode]["misfits"].append(np.loadtxt(ttfile))

#~ The actual plotting stuff

modes_to_plot=6

subplot_layout=plotc.layout_subplots(len(ridge_filters[:modes_to_plot]))[:2]
tdfig,tdiffaxes=plotc.plt.subplots(*subplot_layout)
tdiffaxes=np.array(list([tdiffaxes])).flatten()

c=['olivedrab','burlywood','skyblue','dimgray','red','sienna','black','tomato','grey','teal']

for modeno,mode in enumerate(ridge_filters[:modes_to_plot]):
    num_iterations=min(itercutoff,len(ttdiff_ridges[mode]["iter_no"]))
    plot_every=1
    if num_iterations > 4: plot_every = int(np.ceil(np.sqrt(num_iterations)))
    if num_iterations > 10: plot_every = int(num_iterations//2)
    
    iters_to_plot = [ttdiff_ridges[mode]["iter_no"][i] for i in xrange(0,len(ttdiff_ridges[mode]["iter_no"]),plot_every)]
    
    for color_index,iter_index in enumerate(iters_to_plot):
        
        index_of_iteration=ttdiff_ridges[mode]["iter_no"].index(iter_index)
        td=ttdiff_ridges[mode]["misfits"][index_of_iteration]
        td[:,0]-=1 # fortran index starts from 1, change to python zero based index
        
        left_pix=np.where(td[:,0]<srcloc_ind)[0]
        xcoords_left=np.take(x,td[left_pix,0].astype(int))
        right_pix=np.where(td[:,0]>srcloc_ind)[0]
        xcoords_right=np.take(x,td[right_pix,0].astype(int))    
        
        ax=tdiffaxes[modeno]        
        
        #~ Points to the left
        skip_pix = 3
        ax.plot(xcoords_left[::skip_pix],td[left_pix[::skip_pix],1],color=c[color_index],
                    label="iter "+str(iter_index),marker='o',linestyle='-')

        #~ Points to the right
        ax.plot(xcoords_right[::skip_pix],td[right_pix[::skip_pix],1],color=c[color_index],marker='o',linestyle='-')

        #~ Zero misfit line
        ax.axhline(y=[0],ls='--')
        
        ax.set_title(spaced(modes[mode]),fontsize=18,loc='right')
        
        
for ax_no in xrange(0,len(tdiffaxes),subplot_layout[1]):
    tdiffaxes[ax_no].set_ylabel(r"$\Delta \tau$ (sec)",fontsize=30)

for ax_no in xrange(subplot_layout[1]*(subplot_layout[0]-1),len(tdiffaxes)):    
    tdiffaxes[ax_no].set_xlabel("x (Mm)",fontsize=30)
    
   
for ax_ind,ax in enumerate(tdiffaxes):
    if ax_ind>= len(ridge_filters[:modes_to_plot]): 
        ax.axis('off')
        continue
    #~ Source location line
    ax.axvline(x=srcloc,color='black')
    #~ Rectangle with border
    ax.axvline(x=vlimleft,ls='dotted',color='black')
    ax.axvline(x=vlimright,ls='dotted',color='black')
    ax.axvspan(vlimleft,vlimright,color='paleturquoise')
    ax.xaxis.set_major_locator(MaxNLocator(4,prune="both"))
    ax.yaxis.set_major_locator(MaxNLocator(5,prune='both'))
    ax.legend(loc="best")
    ax.xaxis.set_tick_params(which='major', labelsize=16)
    ax.yaxis.set_tick_params(which='major', labelsize=16)

tdfig.tight_layout()
plt.subplots_adjust(wspace=0.3)
#~ plt.suptitle(str(read_params.get_dt_simulation()),fontsize=20)
plt.show()
