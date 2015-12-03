from __future__ import division
import numpy as np
import plotc
import pyfits as pf
import os,re,sys
import matplotlib.patches as patches
import read_params
plt=plotc.plt

def fitsread(fitsfile): return np.squeeze(pf.getdata(fitsfile))

codedir=os.path.dirname(os.path.abspath(__file__))
datadir=read_params.get_directory()
masterpixelsfile=os.path.join(datadir,'master.pixels')
masterpixels=np.loadtxt(masterpixelsfile,ndmin=1)

dt=0.5
try:
    src=next(int(f) for f in sys.argv if f.isdigit())
except StopIteration:
    src=1
if src>len(masterpixels): src=len(masterpixels)

src=str(src).zfill(2)


#~ Read in velocity and compute length scales

vxfile=os.path.join(codedir,'true_vx.fits')
vzfile=os.path.join(codedir,'true_vz.fits')

vx=fitsread(vxfile)
vz=fitsread(vzfile)
v_surface=np.sqrt(vx[-1]**2+vz[-1]**2)
vx=None;vz=None

vmax=v_surface.max()
bhigh=np.where(v_surface>vmax/2)[0]
v_high_left,v_high_right=bhigh[[0,-1]]
vfwhm=v_high_right-v_high_left
vfwhm=20
vhwhm=vfwhm//2
v_surface=None


datafile=os.path.join(datadir,'forward_src'+src+'_ls00','data.fits')
data=fitsread(datafile)
Nt,Nx=data.shape

Lx=read_params.get_xlength()

x=np.linspace(-Lx/2,Lx/2,Nx,endpoint=False)
time=np.arange(data.shape[0])*dt

vlimleft,vlimright=x[Nx//2-vhwhm],x[Nx//2+vhwhm]

modes={'0':'fmode'}
for pmodeno in xrange(1,6): modes.update({str(pmodeno):'p'+str(pmodeno)+'mode'})
ridge_filters_driver=read_params.get_ridge_filter()
paramsfiles=[os.path.splitext(f)[1][1:] for f in os.listdir(os.path.join(datadir)) if re.match(r'params.[0-9]$',f)]
ridge_filters=sorted([ridge for ridge in ridge_filters_driver if ridge in paramsfiles])

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
    for iteration_no in xrange(0,20):
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
        for color_index,iter_index in enumerate(ttdiff_ridges[mode]["iter_no"][::3]):
            
            td=ttdiff_ridges[mode]["misfits"][ttdiff_ridges[mode]["iter_no"].index(iter_index)]
            #~ print "iteration",iter_index,modes[mode],"sum diff^2",sum((td[:,1]/60)**2)
            
            left=np.where(td[:,0]<srcloc_ind)[0]
            xcoords_left=np.take(x,td[left,0].astype(int))
            
            right=np.where(td[:,0]>srcloc_ind)[0]
            xcoords_right=np.take(x,td[right,0].astype(int))
            
            #~ Points to the left
            ax=plotc.plot1D(td[left,1],x=xcoords_left,ax=tdiffaxes[modeno],color=c[color_index],
                        label="iter "+str(iter_index),marker='o',linestyle='-',
                        axes_properties=dict(locator_properties_x=dict(nbins=4))
                        )
            plt.title(spaced(modes[mode]),dict(fontsize=20))
                        
            #~ Points to the right
            ax=plotc.plot1D(td[right,1],x=xcoords_right,ax=tdiffaxes[modeno],color=c[color_index],
                            marker='o',linestyle='-')

            #~ Zero misfit line
            ax,_=plotc.draw_hlines(y=[0],ax=ax,ls='--')
            #~ Source location line
            ax,_=plotc.draw_vlines(x=[srcloc],ax=ax)
            
            plt.xlabel("x (Mm)",fontsize=20)
            plt.ylabel("Travel Time Difference (sec)",fontsize=20)
        
       
    for ax in tdiffaxes[:len(ridge_filters)]: 
        plotc.draw_vlines(x=[vlimleft,vlimright],ax=ax,ls='dotted')
        plotc.draw_rectangle(x=[vlimleft,vlimright],ax=ax,color='paleturquoise')
        ax.legend(loc='best')


#~ tdfig.tight_layout()
plotc.plt.show()
