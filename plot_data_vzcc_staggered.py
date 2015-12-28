from __future__ import division
import numpy as np
import plotc
import pyfits as pf
import os,re,sys
import matplotlib.patches as patches
from matplotlib.ticker import MaxNLocator
import read_params
plt=plotc.plt

def fitsread(fitsfile): return np.squeeze(pf.getdata(fitsfile))

codedir=os.path.dirname(os.path.abspath(__file__))
datadir=read_params.get_directory()
masterpixelsfile=os.path.join(datadir,'master.pixels')
masterpixels=np.loadtxt(masterpixelsfile,ndmin=1)

dt=0.5

temp=filter(lambda x: x.startswith("src=") or x.startswith("source="),sys.argv)
if temp: src=int(temp[0].split("=")[-1])
else: src=1

if src>len(masterpixels): src=len(masterpixels)
src=str(src).zfill(2)


Nx=read_params.get_nx()
Lx=read_params.get_xlength()
x=np.linspace(-Lx/2,Lx/2,Nx,endpoint=False)

#~ Read in velocity and compute length scales

vxfile=os.path.join(codedir,'true_vx.fits')
vx=fitsread(vxfile)
vx_max_row_index,vx_max_col_index = divmod(vx.argmax(),Nx)
vx_max_row = vx[vx_max_row_index]
vx=None
vxmax=vx_max_row.max()
vx_hwhm=np.where(vx_max_row>vxmax/2)[0][-1]-Nx//2
vlimleft,vlimright=x[Nx//2-vx_hwhm],x[Nx//2+vx_hwhm]

datafile=os.path.join(datadir,'forward_src'+src+'_ls00','data.fits')
data=fitsread(datafile)
Nt=data.shape[0]
time=np.arange(Nt)*dt



modes={'0':'fmode'}
for pmodeno in xrange(1,8): modes.update({str(pmodeno):'p'+str(pmodeno)+'mode'})
ridge_filters=sorted(read_params.get_modes_used())

def spaced(a): return a[:-4]+" "+a[-4:]

vzcc=[]
vzcc_ridges=[]
iters=[]



srcloc=masterpixels[int(src)-1]
srcloc_ind=np.argmin(abs(x-srcloc))
srcloc=round(x[srcloc_ind],1)


#~ What to plot

plot_time_distance=False
plot_data_vzcc_staggered=False
plot_travel_time_misfit=True


#~ Time distance plot

if plot_time_distance:
    plotc.colorplot(data_filtered[0],x=x,y=time,vmax=0.2,centerzero=True,sp=121,
                    xr=[x[receivers[0]]-30,x[receivers[-1]]+30],yr=[0,100],
                    axes_properties=dict(locator_properties_x=dict(nbins=5),title="f mode"),
                    colorbar=False)

    plotc.draw_vlines(x=x[receivers],ls='dashed')

    plotc.colorplot(data_filtered[1],x=x,y=time,vmax=0.2,centerzero=True,sp=122,
                    xr=[x[receivers[0]]-30,x[receivers[-1]]+30],yr=[0,100],
                    axes_properties=dict(hide_yticklabels=True,
                    locator_properties_x=dict(nbins=5),title="p mode"),
                    colorbar=False)

    plotc.draw_vlines(x=x[receivers],ls='dashed')
    plotc.plt.subplots_adjust(wspace=0)

#~ Staggered data vzcc plot

if plot_data_vzcc_staggered:
    
    for it in xrange(0,6):
        vzccfile=os.path.join(ttpath,'vz_cc_src'+src+'.fits')
    
    if os.path.exists(vzccfile):
        vzcc.append(fitsread(vzccfile))
        iters.append(it)

    mode_filters=np.zeros((len(ridge_filters),Nt,Nx))

    for modeno,mode in enumerate(ridge_filters):
        filterfile=os.path.join(codedir,modes[mode]+'_filter.fits')
        mode_filters[modeno]=fitsread(filterfile)

    data_filtered=np.zeros((len(ridge_filters),)+data.shape)
    data_ridge_max=np.ones(len(ridge_filters))
    vzcc_filtered=np.zeros((len(ridge_filters),len(vzcc))+vzcc[0].shape)

    for modeno,_ in enumerate(ridge_filters):
        data_filtered[modeno]=np.fft.ifft2(np.fft.fft2(data)*mode_filters[modeno]).real
        data_ridge_max[modeno]=data_filtered[modeno].max()
        data_filtered[modeno]/=data_ridge_max[modeno]
        
        for iterind,vzccarr in enumerate(vzcc):
            vzcc_filtered[modeno,iterind]=np.real(np.fft.ifft2(np.fft.fft2(vzccarr)*mode_filters[modeno]))/data_ridge_max[modeno]

    vzcc=None


    N_receivers=12
    receivers=[Nx//2-20+5*i for i in xrange(N_receivers)]
    
    
    waveforms_f,axesf=plotc.plt.subplots(N_receivers,sharex=True)
    waveforms_p1,axesp=plotc.plt.subplots(N_receivers,sharex=True)

    for receiver_index,receiver in enumerate(receivers):
        
        axf=axesf[receiver_index]
        axp=axesp[receiver_index]
        
        plotc.plot1D(dataf[:,receiver],x=time,ax=axf,xr=[0,100],
                    label="data",color='red',
                    axes_properties=dict(yscilimits=(-1,1),
                    locator_properties_y=dict(nbins=2))
                    )
                    
        plotc.plot1D(datap[:,receiver],x=time,ax=axp,xr=[0,100],
                    label="data",color='red',
                    axes_properties=dict(yscilimits=(-1,1),
                    locator_properties_y=dict(nbins=2))
                    )

        for iterind in xrange(vzcc_filtered.shape[1]):
            plotc.plot1D(vzcc_filtered[0,:,receiver],x=time,ax=axf,xr=[0,100],
                    label="vzcc"+str(iters[iterind]),
                    axes_properties=dict(yscilimits=(-1,1),
                    locator_properties_y=dict(nbins=2))
                    )
                    
            plotc.plot1D(vzcc_filtered[1,:,receiver],x=time,ax=axp,xr=[0,100],
                    label="vzcc"+str(iters[iterind]),
                    axes_properties=dict(yscilimits=(-1,1),
                    locator_properties_y=dict(nbins=2))
                    )
        
        for ax in [axf,axp]:
            ax_ylim=ax.get_ylim()
            yticks=np.linspace(ax_ylim[0],ax_ylim[1],3)
            ax.set_yticks(yticks)
            ax.text(80,ax_ylim[1]*0.5,"x="+str(round(x[receiver],1)))
            
            #~ Exponent
            exp=ax.get_yaxis().get_offset_text()
            exp.set_x(-0.08)
            
            #~ Tick labels
            ytl=ax.get_yticklabels()
            ytl[0].set_verticalalignment('bottom')
            ytl[-1].set_verticalalignment('top')

    data_filtered=None;vzcc_filtered=None

    for axes in [axesf,axesp]:
        axes[0].legend(bbox_to_anchor=(1, 2))
        axes[-1].set_xlabel("Time (mins)")
        axes[0].set_title("source at x="+str(srcloc))

    waveforms_p1.suptitle("p mode")
    waveforms_f.suptitle("f mode")

    waveforms_f.subplots_adjust(hspace=0)
    waveforms_p1.subplots_adjust(hspace=0)


#~ Travel time misfit plot
if plot_travel_time_misfit:
    
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
    
    subplot_layout=plotc.layout_subplots(len(ridge_filters))[:2]
    tdfig,tdiffaxes=plotc.plt.subplots(*subplot_layout)
    tdiffaxes=np.array(list([tdiffaxes])).flatten()
    
    c=['olivedrab','burlywood','skyblue','dimgray','red','sienna','black','tomato','grey','teal']

    for modeno,mode in enumerate(ridge_filters):
        num_iterations=len(ttdiff_ridges[mode]["iter_no"])
        plot_every=1
        if num_iterations > 4: plot_every = int(np.ceil(np.sqrt(num_iterations)))
        if num_iterations > 10: plot_every = int(num_iterations//3)
        for color_index,iter_index in enumerate(ttdiff_ridges[mode]["iter_no"][::plot_every]):
            
            index_of_iteration=ttdiff_ridges[mode]["iter_no"].index(iter_index)
            td=ttdiff_ridges[mode]["misfits"][index_of_iteration]
            td[:,0]-=1 # fortran index starts from 1, change to python zero based index
            
            left_pix=np.where(td[:,0]<srcloc_ind)[0]
            xcoords_left=np.take(x,td[left_pix,0].astype(int))
            right_pix=np.where(td[:,0]>srcloc_ind)[0]
            xcoords_right=np.take(x,td[right_pix,0].astype(int))            
            
            #~ Points to the left
            plotc.plot1D(td[left_pix,1],x=xcoords_left,ax=tdiffaxes[modeno],color=c[color_index],
                        label="iter "+str(iter_index),marker='o',linestyle='-',
                        axes_properties=dict(locator_properties_x=dict(nbins=4)),
                        title=spaced(modes[mode]),title_properties={'fontsize':15,'loc':'right'}
                        )

            #~ Points to the right
            plotc.plot1D(td[right_pix,1],x=xcoords_right,ax=tdiffaxes[modeno],color=c[color_index],
                            marker='o',linestyle='-')

            #~ Zero misfit line
            plotc.draw_hlines(y=[0],ax=tdiffaxes[modeno],ls='--')
            
            
    for ax_no in xrange(0,len(tdiffaxes),subplot_layout[1]):
        tdiffaxes[ax_no].set_ylabel(r"$\Delta \tau$ (sec)",fontsize=20)
    
    for ax_no in xrange(subplot_layout[1]*(subplot_layout[0]-1),len(tdiffaxes)):    
        tdiffaxes[ax_no].set_xlabel("x (Mm)",fontsize=20)
        
       
    for ax_ind,ax in enumerate(tdiffaxes):
        if ax_ind>= len(ridge_filters): 
            ax.axis('off')
            continue
        #~ Source location line
        ax.axvline(x=srcloc,color='black')
        #~ Rectangle with border
        ax.axvline(x=vlimleft,ls='dotted',color='black')
        ax.axvline(x=vlimright,ls='dotted',color='black')
        ax.axvspan(vlimleft,vlimright,color='paleturquoise')
        ax.legend(loc='best')
        ax.yaxis.set_major_locator(MaxNLocator(4,prune='both'))


    tdfig.tight_layout()
plotc.plt.show()
