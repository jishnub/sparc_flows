from __future__ import division
import numpy as np
import plotc
import pyfits
import os,re,sys,glob
from matplotlib.ticker import MaxNLocator,FuncFormatter
import read_params
from matplotlib.widgets import CheckButtons
import matplotlib.pyplot as plt
import matplotlib
import modefilters

matplotlib.rcParams['axes.formatter.useoffset']=False

def fitsread(fitsfile): return np.squeeze(pyfits.getdata(fitsfile))

datadir = read_params.get_directory()

src = read_params.parse_cmd_line_params("src",mapto=int,default=1)

mp=np.loadtxt(os.path.join(datadir,'master.pixels'),ndmin=1)
srcloc = mp[src]

Lx=read_params.get_xlength()
nx = read_params.get_nx()
x= np.linspace(-Lx/2,Lx/2,nx,endpoint=False)

if filter(lambda x: x.startswith("dist"),sys.argv): coord='distance'
elif filter(lambda x: x.startswith("pix"),sys.argv): coord='pixel'
elif filter(lambda x: x.startswith("coord"),sys.argv): coord='coord'
else: coord = 'unspecified'
pix=130

if coord=='distance':
    dist=float(filter(lambda x: x.startswith("dist"),sys.argv)[0].split("=")[-1])
    xcoord = (srcloc + dist) % Lx 
    if xcoord > Lx/2: xcoord = xcoord - Lx
    pix = abs(x-xcoord).argmin()
        
elif coord == 'coord':
    xcoord=float(filter(lambda x: x.startswith("coord"),sys.argv)[0].split("=")[-1])
    pix = abs(x-xcoord).argmin()

elif coord=='pixel':
    if len(filter(lambda x: x.startswith("pix"),sys.argv)):
        pix=int(filter(lambda x: x.startswith("pix"),sys.argv)[0].split("=")[-1])
        
else:
    print "Using default pixel",pix
    print "Pass pixel as pix=<pix no>, eg pix=123; or distance as dist=<dist>; or coordinate as coord=<coord>"

data=fitsread(os.path.join(datadir,'tt','data','data'+str(src).zfill(2)+'.fits'))
nt,nx = data.shape

dt_sec = read_params.get_dt()
dt_min = dt_sec/60
t=np.arange(nt)*dt_min

vzcc_start=fitsread(os.path.join(datadir,'tt','iter00','vz_cc_src'+str(src).zfill(2)+'.fits'))

itercutoff = read_params.parse_cmd_line_params("iter",mapto=int,default=np.inf)

iters_done = glob.glob(os.path.join(datadir,'tt','iter*'))
niters = len(iters_done)-1
iter_to_plot = min(itercutoff,niters)
iter_dir = os.path.join(datadir,'tt','iter'+str(iter_to_plot).zfill(2))
vzcc_iter = fitsread(os.path.join(iter_dir,'vz_cc_src'+str(src).zfill(2)+'.fits'))

modes={'0':'f'}
for i in xrange(1,8): modes[str(i)]='p'+str(i)
modes['8']='first_bounce_p'

modes_used=read_params.get_modes_used()

modes_passed_cmdline = read_params.parse_cmd_line_params("mode",default=["f"],return_list=True)
plotridge=filter(lambda x: x in modes.values(),modes_passed_cmdline)

if plotridge:
    ridges=[r+'mode' for r in modes_passed_cmdline if r in plotridge]
else:
    ridges=[modes[modes_used[0]]+'mode']


for ridge in ridges:
    print ridge
    plt.figure()
    filt = fitsread(ridge+'_filter.fits')

    data_filtered = np.fft.ifft2(np.fft.fft2(data)*filt).real
    vzcc_start_filtered = np.fft.ifft2(np.fft.fft2(vzcc_start)*filt).real
    vzcc_iter_filtered = np.fft.ifft2(np.fft.fft2(vzcc_iter)*filt).real
    
    #~ Apply first bounce filter
    if ridge == "first_bounce_pmode":
        fbfilter = np.squeeze(modefilters.first_bounce_filter(nt,dt_sec,nx,Lx,srcloc)).T
        data_filtered*=fbfilter
        vzcc_start_filtered*=fbfilter
        vzcc_iter_filtered*=fbfilter

    #~ Try to load tt file
    ttfile = os.path.join(iter_dir,'ttdiff_src'+str(src).zfill(2)+'.'+ridge)
    ttdiff_array = np.loadtxt(ttfile)
    pixel_ttdiff = ttdiff_array[:,0]
    lef_ttdiff = ttdiff_array[:,2]
    rig_ttdiff = ttdiff_array[:,3]
    
    lef_ttdiff_cutoff = t[map(int,lef_ttdiff[pixel_ttdiff==pix])]
    rig_ttdiff_cutoff = t[map(int,rig_ttdiff[pixel_ttdiff==pix])]
    
    plt.plot(t,data_filtered[:,pix],label="True",color='black')
    plt.plot(t,vzcc_start_filtered[:,pix],label="Iter 0",linestyle='dashed',linewidth=2,color='#555555')
    plt.plot(t,vzcc_iter_filtered[:,pix],label="Iter "+str(iter_to_plot),linestyle='dashed',linewidth=2,color='#333333',marker='o')
    ax1=plt.gca()
    ax1.yaxis.set_major_locator(MaxNLocator(4,prune='both'))
    ax1.xaxis.set_major_locator(MaxNLocator(6))
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    #~ ax1.get_yaxis().get_offset_text().set_x(-0.07)
    
    #~ plt.legend(loc='best')
    
    #~ ax1_twin=ax1.twiny()
    #~ ax1_twin.plot(t,vzcc_start_filtered[:,pix],label="Starting",linestyle='dashed',linewidth=2,color='red')
    #~ ax1_twin.xaxis.set_major_locator(MaxNLocator(6,prune='both'))
    #~ ax1_twin.xaxis.set_major_formatter(FuncFormatter(lambda x,pos: "{:3.1f}".format(x*60)))
    
    ax1.set_ylabel("Wave velocity\n(arbitrary units)")
    
    xlim_left = lef_ttdiff_cutoff - (rig_ttdiff_cutoff - lef_ttdiff_cutoff)/3
    xlim_right = rig_ttdiff_cutoff + (rig_ttdiff_cutoff - lef_ttdiff_cutoff)/3
    plt.xlim(xlim_left,xlim_right)
    plt.legend(loc='best')
    
    plt.xlabel("Time (min)")
    
    plotc.apj_1col_format(plt.gcf())
  
  
plt.tight_layout()

save = read_params.parse_cmd_line_params("save")
if save is not None:
    savepath = os.path.join("plots",save)
    print "saving to",savepath
    if not os.path.exists("plots"): os.makedirs("plots")
    plt.savefig(savepath)
else:
    print "Not saving plot to file"
    
plt.show()

