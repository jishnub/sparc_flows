from __future__ import division
import numpy as np
import plotc
import pyfits as pf
import os,re,sys
import matplotlib
from matplotlib import rc,gridspec
import matplotlib.patches as patches
from matplotlib.ticker import MaxNLocator
import read_params
import matplotlib.pyplot as plt
import itertools,warnings
from scipy.special import j1

def fitsread(fitsfile): return np.squeeze(pf.getdata(fitsfile))

datadir=read_params.get_directory()
masterpixelsfile=os.path.join(datadir,'master.pixels')
masterpixels=np.loadtxt(masterpixelsfile,ndmin=1)

dt=read_params.get_dt()/60.

src = read_params.parse_cmd_line_params("src",mapto=int,default=1)

if src>len(masterpixels): src=len(masterpixels)
src=str(src).zfill(2)


nx=read_params.get_nx()
Lx=read_params.get_xlength()
x=np.linspace(-Lx/2,Lx/2,nx,endpoint=False)
nz = read_params.get_nz()

#DH13 model parameters
RDH13 = 15
kDH13 = 2*np.pi/(2*RDH13)

flows = read_params.if_flows()
sound_speed_perturbed = read_params.if_soundspeed_perturbed()


datafile=os.path.join(datadir,'data',str(src).zfill(2)+'.fits')
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
            ttfile=os.path.join(ttpath,'ttdiff_src'+src+'.'+modes[mode]+'.0')
            if os.path.exists(ttfile):
                ttdiff_ridges[mode]["iter_no"].append(iteration_no)
                ttdiff_ridges[mode]["misfits"].append(np.loadtxt(ttfile))

fig = plt.figure()

c=['black','#555555','#888888']

# vspans to indicate extent of flow
numspans = 40
spanedges = np.linspace(-3*RDH13,3*RDH13,numspans)
def vx_r_fn(x): return j1(kDH13*abs(x))*np.exp(-abs(x)/RDH13)
vx_x = vx_r_fn(x)*np.sign(x)
norm = matplotlib.colors.Normalize(vx_x.min()*2, vx_x.max()*2)
my_cmap = matplotlib.cm.get_cmap("bwr")

# list of axes
axlist = []

for modeno,mode in enumerate(ridge_filters[:8]):

    iters_to_plot_ridge = ttdiff_ridges[mode]["iter_no"]

    linestyles = itertools.cycle(('solid','dashed','dotted','dashdot'))

    ax = fig.add_subplot(3,3,modeno+1)
    axlist.append(ax)

    if len(iters_to_plot_ridge)>3:
        iters_to_plot_ridge= iters_to_plot_ridge[::3]

    for color_index,iter_no in enumerate(iters_to_plot_ridge):

        index_of_iteration=ttdiff_ridges[mode]["iter_no"].index(iter_no)
        td=ttdiff_ridges[mode]["misfits"][index_of_iteration]
        td[:,0]-=1 # fortran index starts from 1, change to python zero based index

        left_pix=np.where(td[:,0]<srcloc_ind)[0]
        xcoords_left=np.take(x,td[left_pix,0].astype(int))
        right_pix=np.where(td[:,0]>srcloc_ind)[0]
        xcoords_right=np.take(x,td[right_pix,0].astype(int))

        #~ Points to the left
        skip_pix = 1
        ls = next(linestyles)
        color = c[color_index%len(c)]
        p=ax.plot(xcoords_left[::skip_pix],td[left_pix[::skip_pix],1],
                color=color,label=str(iter_no),
                linestyle=ls,linewidth=2,zorder=2)[0]

        #~ Points to the right
        ax.plot(xcoords_right[::skip_pix],td[right_pix[::skip_pix],1],
        color=p.get_color(),linestyle=ls,linewidth=2,zorder=2)

    ax.set_title(spaced(modes[mode]),fontsize=18,loc='right')
    ax.axhline(y=[0],ls='dotted',color='black')
    ax.axvline(x=srcloc,color='#333333',ls='dashed')

    ax.xaxis.set_major_locator(MaxNLocator(4,prune="both"))
    ax.yaxis.set_major_locator(MaxNLocator(5,prune='both'))
    ax.xaxis.set_tick_params(which='major', labelsize=14)
    ax.yaxis.set_tick_params(which='major', labelsize=14)
    ymin,ymax = ax.get_ylim()
    ax.set_ylim(-max(abs(ymin),abs(ymax))*1.5,max(abs(ymin),abs(ymax))*1.5)
    ax.legend(loc="best",fontsize=12,ncol=2)

    for span_left,span_right in zip(spanedges[:-1],spanedges[1:]):
        span_color = my_cmap(norm(vx_r_fn((span_left+span_right)/2)))
        ax.axvspan(span_left,span_right,color=span_color,zorder=0)



for axno,ax in enumerate(axlist):
    ncols = min(3,len(axlist))
    nrows = int(np.ceil(len(axlist)/ncols))
    ax.change_geometry(nrows,ncols,axno+1)

plt.tight_layout()

#~ Get filename to save to
save = read_params.parse_cmd_line_params("save")

if save is not None:
    savepath = os.path.join("plots",save)
    if not os.path.exists("plots"): os.makedirs("plots")
    plt.savefig(savepath)
    print "Saved to",savepath
else:
    print "Not saving plot to file"

plt.show()
