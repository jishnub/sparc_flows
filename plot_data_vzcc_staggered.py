from __future__ import division
import numpy as np
import plotc
import pyfits as pf
import os
import matplotlib.patches as patches

def fitsread(fitsfile): return np.squeeze(pf.getdata(fitsfile))

datadir="/scratch/shivam/flows/data2"
codedir="/home/shivam/sparc2"
masterpixelsfile=os.path.join(datadir,'master.pixels')
masterpixels=np.loadtxt(masterpixelsfile,ndmin=1)

dt=0.5
src=3
if src>len(masterpixels): src=len(masterpixels)

src=str(src).zfill(2)



Bxfile='/home/shivam/sparc2/true_vx.fits'
Bzfile='/home/shivam/sparc2/true_vz.fits'

bx=fitsread(Bxfile)
bz=fitsread(Bzfile)
b_surface=np.sqrt(bx[-1]**2+bz[-1]**2)
bx=None;bz=None
Nx=len(b_surface)

bmax=b_surface.max()
bhigh=np.where(b_surface>bmax/2)[0]
bhigh_left,bhigh_right=bhigh[[0,-1]]
bfwhm=bhigh_right-bhigh_left
bfwhm=20
bhwhm=bfwhm//2
b_surface=None


datafile=os.path.join(datadir,'forward_src'+src+'_ls00','data.fits')
data=fitsread(datafile)

x=np.linspace(-400,400,Nx,endpoint=False)
time=np.arange(data.shape[0])*dt

blimleft,blimright=x[Nx//2-bhwhm],x[Nx//2+bhwhm]

vzcc=[]
vzccf=[]
vzccp=[]
iters=[]
ttdiff_f=[]
ttdiff_p=[]

for it in xrange(14,18):
    ttpath=os.path.join(datadir,'tt','iter'+str(it).zfill(2))
    vzccfile=os.path.join(ttpath,'vz_cc_src'+src+'.fits')
    ttffile=os.path.join(ttpath,'ttdiff_src'+src+'.fmode')
    ttpfile=os.path.join(ttpath,'ttdiff_src'+src+'.p1mode')
    ttp2file=os.path.join(ttpath,'ttdiff_src'+src+'.p2mode')
    if os.path.exists(vzccfile):
        vzcc.append(fitsread(vzccfile))
        iters.append(it)
        ttdiff_f.append(np.loadtxt(ttffile))
        ttdiff_p.append(np.loadtxt(ttpfile))

fmode_filter_file=os.path.join(codedir,'fmode_filter.fits')
fmode=fitsread(fmode_filter_file)

pmode_filter_file=os.path.join(codedir,'pmode_filter.fits')
pmode=fitsread(pmode_filter_file)

srcloc=masterpixels[int(src)-1]
srcloc_ind=np.argmin(abs(x-srcloc))
srcloc=round(x[srcloc_ind],1)

dataf=np.real(np.fft.ifft2(np.fft.fft2(data)*fmode))
datafmax=dataf.max()
dataf/=datafmax

datap=np.real(np.fft.ifft2(np.fft.fft2(data)*pmode))
datapmax=datap.max()
datap/=datapmax

for vzccarr in vzcc:
    vzccfarr=np.real(np.fft.ifft2(np.fft.fft2(vzccarr)*fmode))/datafmax
    vzccf.append(vzccfarr)
    
    vzccparr=np.real(np.fft.ifft2(np.fft.fft2(vzccarr)*pmode))/datapmax
    vzccp.append(vzccparr)

vzcc=None


Nrec=12
recs=[225+5*i for i in xrange(Nrec)]

plotc.colorplot(dataf,x=x,y=time,vmax=0.2,centerzero=True,sp=121,
                xr=[x[recs[0]]-30,x[recs[-1]]+30],yr=[0,100],
                axes_properties=dict(
                locator_properties_x=dict(nbins=5),title="f mode"),
                colorbar=False)

plotc.draw_vlines(x=x[recs],ls='dashed')

plotc.colorplot(datap,x=x,y=time,vmax=0.2,centerzero=True,sp=122,
                xr=[x[recs[0]]-30,x[recs[-1]]+30],yr=[0,100],
                axes_properties=dict(hide_yticklabels=True,
                locator_properties_x=dict(nbins=5),title="p mode"),
                colorbar=False)

plotc.draw_vlines(x=x[recs],ls='dashed')
plotc.plt.subplots_adjust(wspace=0)

waveforms_f,axesf=plotc.plt.subplots(Nrec,sharex=True)
waveforms_p,axesp=plotc.plt.subplots(Nrec,sharex=True)

for recind,rec in enumerate(recs):
    
    axf=axesf[recind]
    axp=axesp[recind]
    
    plotc.plot1D(dataf[:,rec],x=time,ax=axf,xr=[0,100],
                label="data",color='red',
                axes_properties=dict(yscilimits=(-1,1),
                locator_properties_y=dict(nbins=2))
                )
                
    plotc.plot1D(datap[:,rec],x=time,ax=axp,xr=[0,100],
                label="data",color='red',
                axes_properties=dict(yscilimits=(-1,1),
                locator_properties_y=dict(nbins=2))
                )

    for iterind,vzccfarr in enumerate(vzccf):
        plotc.plot1D(vzccfarr[:,rec],x=time,ax=axf,xr=[0,100],
                label="vzcc"+str(iterind),
                axes_properties=dict(yscilimits=(-1,1),
                locator_properties_y=dict(nbins=2))
                )
                
    for iterind,vzccparr in enumerate(vzccp):
        plotc.plot1D(vzccparr[:,rec],x=time,ax=axp,xr=[0,100],
                label="vzcc"+str(iterind),
                axes_properties=dict(yscilimits=(-1,1),
                locator_properties_y=dict(nbins=2))
                )
    
    for ax in [axf,axp]:
        ax_ylim=ax.get_ylim()
        yticks=np.linspace(ax_ylim[0],ax_ylim[1],3)
        ax.set_yticks(yticks)
        ax.text(80,ax_ylim[1]*0.5,"x="+str(round(x[rec],1)))
        
        #~ Exponent
        exp=ax.get_yaxis().get_offset_text()
        exp.set_x(-0.08)
        
        #~ Tick labels
        ytl=ax.get_yticklabels()
        ytl[0].set_verticalalignment('bottom')
        ytl[-1].set_verticalalignment('top')

datap=None;dataf=None;vzccf=None;vzccp=None

for axes in [axesf,axesp]:
    axes[0].legend(bbox_to_anchor=(1, 2))
    axes[-1].set_xlabel("Time (mins)")
    axes[0].set_title("source at x="+str(srcloc))

waveforms_p.suptitle("p mode")
waveforms_f.suptitle("f mode")

waveforms_f.subplots_adjust(hspace=0)
waveforms_p.subplots_adjust(hspace=0)



tdfig,tdiffaxes=plotc.plt.subplots(1,2)
c=['green','blue','cyan','magenta','red','yellow','black','brown','grey','teal']

for i,td in enumerate(ttdiff_f):
    
    print "iteration",i,"f mode sum diff^2",sum((td[:,1]/60)**2)
    
    left=np.where(td[:,0]<srcloc_ind)[0]
    xcoords_left=np.take(x,td[left,0].astype(int))
    
    right=np.where(td[:,0]>srcloc_ind)[0]
    xcoords_right=np.take(x,td[right,0].astype(int))

    ax=plotc.plot1D(td[left,1],x=xcoords_left,ax=tdiffaxes[0],color=c[i],
                label="iter "+str(i),marker='o',linestyle='-',
                axes_properties=dict(locator_properties_x=dict(nbins=4),title="f mode")
                )
    
    ax=plotc.plot1D(td[right,1],x=xcoords_right,ax=ax,
                marker='o',linestyle='-',color=c[i],
                axes_properties=dict(ylabel="$\Delta$t (seconds)",xlabel="Distance (Mm)",
                locator_properties_x=dict(nbins=4))
                )
    
    ax,_=plotc.draw_hlines(y=[0],ax=ax,ls='--')
    ax,_=plotc.draw_vlines(x=[srcloc],ax=ax)
    
        
for i,td in enumerate(ttdiff_p):
    
    print "iteration",i,"p mode sum diff^2",sum((td[:,1]/60)**2)
    
    left=np.where(td[:,0]<srcloc_ind)[0]
    xcoords_left=np.take(x,td[left,0].astype(int))
    
    right=np.where(td[:,0]>srcloc_ind)[0]
    xcoords_right=np.take(x,td[right,0].astype(int))
    
    ax=plotc.plot1D(td[left,1],x=xcoords_left,ax=tdiffaxes[1],
                label="iter "+str(i),marker='o',linestyle='-',color=c[i],
                axes_properties=dict(locator_properties_x=dict(nbins=6),title="p mode")
                )
    
    ax=plotc.plot1D(td[right,1],x=xcoords_right,ax=ax,
                marker='o',linestyle='-',color=c[i],
                axes_properties=dict(ylabel="$\Delta$t (seconds)",xlabel="Distance (Mm)",
                locator_properties_x=dict(nbins=6))                
                )
    
    ax,_=plotc.draw_hlines(y=[0],ax=ax,ls='--')    
    ax,_=plotc.draw_vlines(x=[srcloc],ax=ax)

    
   
for ax in tdiffaxes: 
    plotc.draw_vlines(x=[blimleft,blimright],ax=ax,ls='dotted')
    plotc.draw_rectangle(x=[blimleft,blimright],ax=ax,alpha=0.1)
    ax.legend(loc='best')


#~ tdfig.tight_layout()
plotc.plt.show()
