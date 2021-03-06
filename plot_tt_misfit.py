
import numpy as np
import plotc
import pyfits as pf
import os,re,sys
from matplotlib import rc
import matplotlib.patches as patches
from matplotlib.ticker import MaxNLocator
import read_params
import matplotlib.pyplot as plt
import itertools

def fitsread(fitsfile): return np.squeeze(pf.getdata(fitsfile))



datadir=read_params.get_directory()
masterpixelsfile=os.path.join(datadir,'master.pixels')
masterpixels=np.loadtxt(masterpixelsfile,ndmin=1)

dt=read_params.get_dt()/60.

src = read_params.parse_cmd_line_params("src",mapto=int,default=1)
itercutoff = read_params.parse_cmd_line_params("iter",mapto=int,default=np.inf)

if src>len(masterpixels): src=len(masterpixels)
src=str(src).zfill(2)


nx=read_params.get_nx()
Lx=read_params.get_xlength()
x=np.linspace(-Lx/2,Lx/2,nx,endpoint=False)
nz = read_params.get_nz()

flows = read_params.if_flows()
sound_speed_perturbed = read_params.if_soundspeed_perturbed()

#~ Read in velocity and compute length scales
if flows:
    vx=fitsread('true_vx.fits')
    vx_max_row_index,vx_max_col_index = divmod(vx.argmax(),nx)
    vx_max_row = vx[vx_max_row_index]
    vx=None
    vxmax=vx_max_row.max()
    hwhm=np.where(vx_max_row>vxmax/2)[0][-1]-nx//2
    
elif sound_speed_perturbed:
    true_sound_speed = fitsread('true_c.fits')
    one_D_sound_speed = np.tile(np.loadtxt(read_params.get_solarmodel(),usecols=[1],ndmin=2).reshape(nz,1),(1,nx))
    ss_pert = (true_sound_speed - one_D_sound_speed)/one_D_sound_speed
    ss_pert_max_row_index,ss_pert_max_col_index = divmod(ss_pert.argmax(),nx)
    ss_pert_max_row = ss_pert[ss_pert_max_row_index]
    ss_pert_max = ss_pert_max_row.max()
    hwhm = np.where(ss_pert_max_row>ss_pert_max/2)[0][-1]-nx//2

perturbation_left_limit,perturbation_right_limit=x[nx//2-hwhm],x[nx//2+hwhm]

datafile=os.path.join(datadir,'forward_src'+src+'_ls00','data.fits')
data=fitsread(datafile)
nt=data.shape[0]
time=np.arange(nt)*dt



modes={'0':'fmode'}
for pmodeno in range(1,8): modes.update({str(pmodeno):'p'+str(pmodeno)+'mode'})
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
for iteration_no in range(0,100):
    ttpath=os.path.join(datadir,'tt','iter'+str(iteration_no).zfill(2))
    if os.path.exists(ttpath):
        for modeno,mode in enumerate(ridge_filters):
            ttfile=os.path.join(ttpath,'ttdiff_src'+src+'.'+modes[mode])
            if os.path.exists(ttfile):
                ttdiff_ridges[mode]["iter_no"].append(iteration_no)
                ttdiff_ridges[mode]["misfits"].append(np.loadtxt(ttfile))

#~ The actual plotting stuff
rc('text', usetex=True)
rc('font',**{'family':'serif','serif':['Times'],'size':24})

modes_to_plot=min(8,len(ridge_filters))

subplot_layout=plotc.layout_subplots(len(ridge_filters[:modes_to_plot]))[:2]

tdfig,tdiffaxes=plt.subplots(*subplot_layout)
tdiffaxes=np.array(list([tdiffaxes])).flatten()

c=['black','#555555','#888888']

iters_to_plot = sorted(read_params.parse_cmd_line_params("iter",return_list=True,default=[],mapto=int))

for modeno,mode in enumerate(ridge_filters[:modes_to_plot]):
    num_iterations=min(itercutoff,len(ttdiff_ridges[mode]["iter_no"]))
    plot_every=1
    if num_iterations > 4: plot_every = int(np.ceil(np.sqrt(num_iterations)))
    if num_iterations > 10: plot_every = int(num_iterations//2)
    
    iters_to_plot_ridge = [ttdiff_ridges[mode]["iter_no"][i] for i in range(0,len(ttdiff_ridges[mode]["iter_no"]),plot_every)]
    iters_to_plot_temp = iters_to_plot
    if max(iters_to_plot)>iters_to_plot_ridge[-1]:
        for entry in iters_to_plot:
            if entry>iters_to_plot_ridge[-1]: break
        iters_to_plot_temp[iters_to_plot.index(entry):]=[iters_to_plot_ridge[-1]]
    iters_to_plot_ridge = iters_to_plot_temp
     
    
    linestyles = itertools.cycle(('solid','dashed','dotted','dashdot'))
    
    for color_index,iter_index in enumerate(iters_to_plot_ridge):
        
        index_of_iteration=ttdiff_ridges[mode]["iter_no"].index(iter_index)
        td=ttdiff_ridges[mode]["misfits"][index_of_iteration]
        td[:,0]-=1 # fortran index starts from 1, change to python zero based index
        
        left_pix=np.where(td[:,0]<srcloc_ind)[0]
        xcoords_left=np.take(x,td[left_pix,0].astype(int))
        right_pix=np.where(td[:,0]>srcloc_ind)[0]
        xcoords_right=np.take(x,td[right_pix,0].astype(int))    
        
        ax=tdiffaxes[modeno]        
        
        #~ Points to the left
        skip_pix = 1
        ls = next(linestyles)
        ax.plot(xcoords_left[::skip_pix],td[left_pix[::skip_pix],1],color=c[color_index//3],
                    label="iter "+str(iter_index),linestyle=ls,linewidth=2)

        #~ Points to the right
        ax.plot(xcoords_right[::skip_pix],td[right_pix[::skip_pix],1],color=c[color_index//3],linestyle=ls,linewidth=2)

        #~ Zero misfit line
        ax.axhline(y=[0],ls='-',color='black')
        
        ax.set_title(spaced(modes[mode]),fontsize=18,loc='right')
        
        
for ax_no in range(0,len(tdiffaxes),subplot_layout[1]):
    tdiffaxes[ax_no].set_ylabel(r"$\Delta \tau$ (sec)",fontsize=20)

for ax_no in range(subplot_layout[1]*(subplot_layout[0]-1),len(tdiffaxes)):    
    tdiffaxes[ax_no].set_xlabel("x (Mm)",fontsize=20)
    
   
for ax_ind,ax in enumerate(tdiffaxes):
    if ax_ind>= modes_to_plot: 
        ax.axis('off')
        continue
    #~ Source location line
    ax.axvline(x=srcloc,color='black')
    #~ Rectangle with border
    ax.axvline(x=perturbation_left_limit,ls='dotted',color='black')
    ax.axvline(x=perturbation_right_limit,ls='dotted',color='black')
    ax.axvspan(perturbation_left_limit,perturbation_right_limit,color='#CCCCCC')
    ax.xaxis.set_major_locator(MaxNLocator(4,prune="both"))
    ax.yaxis.set_major_locator(MaxNLocator(5,prune='both'))
    ax.xaxis.set_tick_params(which='major', labelsize=16)
    ax.yaxis.set_tick_params(which='major', labelsize=16)

if modes_to_plot<np.prod(subplot_layout):
    # This means that there are empty subplots
    # Place legend in empty spot
    handles,labels = tdiffaxes[0].get_legend_handles_labels()
    tdiffaxes[-1].legend(handles,labels,loc="center",fontsize=18)
else:
    for ax in tdiffaxes[:modes_to_plot]: ax.legend(loc="best")

tdfig.tight_layout()
plt.subplots_adjust(wspace=0.3)

#~ plotc.apj_2col_format(plt.gcf())

plt.gcf().set_size_inches(12,8)
plt.tight_layout()

#~ Get filename to save to
save = read_params.parse_cmd_line_params("save")

if save is not None:
    savepath = os.path.join("plots",save)
    if not os.path.exists("plots"): os.makedirs("plots")
    plt.savefig(savepath)
    print("Saved to",savepath)
else:
    print("Not saving plot to file")

plt.show()
