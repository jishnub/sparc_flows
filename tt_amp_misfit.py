from __future__ import division
import traveltimes
import numpy as np
import pyfits
import read_params
import os,sys,fnmatch
from matplotlib import pyplot as plt
import plotc

datadir = read_params.get_directory()
updatedir = os.path.join(datadir,'update')

def parse_cmd_line_params(key,mapto=None,default=None):
    cmd_line_param = filter(lambda x: x.startswith(key+"="),sys.argv)
        
    if len(cmd_line_param)!=0: 
        retval=cmd_line_param[0].split("=")[-1]
        retval = retval.strip().lstrip('[').rstrip(']').split(',')
        if mapto is not None: retval=map(mapto,retval)
        if len(retval)==1: retval=retval[0]
            
    else: retval=default
    
    return retval
    

src = parse_cmd_line_params('src',mapto=int,default=1)

misfit_type = parse_cmd_line_params('type',default='tt')
modes = parse_cmd_line_params('mode',default='f')
if type(modes)==str: modes=[modes]

modeindex={'f':'0'}
for i in xrange(1,8): modeindex['p'+str(i)]=str(i)

#~ Get a list of misfit files in the update directory
misfitfiles=sorted([os.path.join(updatedir,f) for f in fnmatch.filter(os.listdir(updatedir),'misfit_[0-9][0-9]')])
num_misfit_files=len(misfitfiles)

#~ If iteration cutoff is specified use it
itercutoff = parse_cmd_line_params('iter',mapto=int,default=np.inf)

iter_no = min(itercutoff,num_misfit_files-1)

datafile = os.path.join(datadir,'tt','data','data'+str(src).zfill(2)+'.fits')
data=pyfits.getdata(datafile)

vzccfile = os.path.join(datadir,'tt','iter'+str(iter_no).zfill(2),'vz_cc_src'+str(src).zfill(2)+'.fits')
vzcc = pyfits.getdata(vzccfile)

nt,_,nx = data.shape
dt = read_params.get_dt()/60 # Minutes
Lx = read_params.get_xlength()

x = np.linspace(-Lx/2,Lx/2,nx,endpoint=False)

print "src",src,"iter",iter_no,"modes",modes

if misfit_type == 'tt':
    
    for mode in modes:
        window_file = os.path.join(datadir,'tt','iter'+str(iter_no).zfill(2),'ttdiff_src'+str(src).zfill(2)+'.'+mode+'mode')
        pixels_list,left_list,rig_list = np.loadtxt(window_file,usecols=[0,2,3],unpack=True)
        pixels_list = pixels_list.astype(int)
        pixels_list-=1 # Fortran index starts from 1
        left_list-=1
        rig_list-=1
        pixels_x = x[pixels_list]
        inner = np.loadtxt(os.path.join(datadir,'params.'+modeindex[mode])).T[0]
        
        modefilter = pyfits.getdata(mode+'mode_filter.fits')

        dat_filtered = np.fft.ifft2(np.fft.fft2(data,axes=(0,2))*modefilter,axes=(0,2)).real
        vzcc_filtered = np.fft.ifft2(np.fft.fft2(vzcc,axes=(0,2))*modefilter,axes=(0,2)).real
        
        #~ plt.figure()
        #~ plotc.colorplot(dat_filtered,sp=121,vmax=1e-4,centerzero=True)
        #~ plotc.colorplot(vzcc_filtered,sp=122,vmax=1e-4,centerzero=True)
        
        tt_misfits = []
        
        for pix,lef,rig in zip(pixels_list,left_list,rig_list):
            
            tt_han = traveltimes.compute_tt_hanasoge(vzcc_filtered[lef:rig+1,0,pix],dat_filtered[lef:rig+1,0,pix],dt)*60
            
            tt_gb = traveltimes.compute_tt_gizonbirch(vzcc_filtered[:,0,pix],dat_filtered[:,0,pix],dt,lef+1,rig+1)*60
            
            tt_misfits.append((tt_han,tt_gb))
        
        tt_misfits = np.array(tt_misfits)
        tt_han = tt_misfits[:,0]
        tt_gb = tt_misfits[:,1]
        
        plt.figure()
        plt.plot(pixels_x[pixels_x<-inner],tt_han[pixels_x<-inner],'bo-',label="Hanasoge")
        plt.plot(pixels_x[pixels_x>inner],tt_han[pixels_x>inner],'bo-')
        
        plt.plot(pixels_x[pixels_x<-inner],tt_gb[pixels_x<-inner],'go-',label="Gizon-Birch")
        plt.plot(pixels_x[pixels_x>inner],tt_gb[pixels_x>inner],'go-')
        
        plt.legend(loc='best')
        plt.title(mode+" mode")
    plt.show()
    
    
elif misfit_type == 'amp':
    pass
