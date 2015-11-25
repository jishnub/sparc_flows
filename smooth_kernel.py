from __future__ import division
import numpy as np
import scipy.interpolate
import pyfits
import os
from scipy.ndimage.filters import gaussian_filter
import matplotlib.pyplot as plt
import plotc
import read_params

codedir=os.path.dirname(os.path.abspath(__file__))
datadir=read_params.get_directory()

Lx=read_params.get_xlength()

kvx=np.squeeze(pyfits.getdata(os.path.join(datadir,'kernel','kernel_vx_01.fits')))
kvz=np.squeeze(pyfits.getdata(os.path.join(datadir,'kernel','kernel_vz_01.fits')))
nz,nx=kvx.shape

x=np.linspace(-Lx/2,Lx/2,num=nx,endpoint=False)

Rsun=695.8
solar_model=read_params.get_solarmodel()
z=np.loadtxt(solar_model,usecols=[0]);z=(z-1)*Rsun

kvx_filtered=gaussian_filter(kvx,(2,2))
kvz_filtered=gaussian_filter(kvz,(2,2))

src_x = np.loadtxt(os.path.join(datadir,'master.pixels'),usecols=[0])
srcpix = abs(x-src_x).argmin()

recpix = 128
with open(os.path.join(codedir,'driver.f90'),'r') as drv:
    adjoint_source_filt=False
    rec_pix_flag=False
    for line in drv:
        if line.strip()=="SUBROUTINE ADJOINT_SOURCE_FILT(nt)": adjoint_source_filt=True
        if line.strip()=="END SUBROUTINE ADJOINT_SOURCE_FILT": 
            adjoint_source_filt=False
            break
        if not adjoint_source_filt: continue
        if "RECEIVER PIXEL FLAG" in line.strip(): 
            rec_pix_flag=True
        if "RECEIVER PIXEL END FLAG" in line.strip():
            rec_pix_flag=False
            break
        if ("do" in line.strip()) and rec_pix_flag:
            recpix=int(line.strip().split(",")[-1])
            break

rec_x = x[recpix]

#~ print src_x,rec_x
#~ print "Source pixel",srcpix,"Receiver pixel",recpix
srcrec_midpoint = (srcpix + recpix)//2
#~ print "midpoint pixel",srcrec_midpoint
#~ print "cell center",nx//2
#~ print srcrec_midpoint,x[srcrec_midpoint]
#~ print "offset",nx//2 - srcrec_midpoint


def asym(arr): return (arr-arr[:,::-1])/2
def sym(arr):
    arrflip = np.fliplr(arr)
    offset = srcrec_midpoint - nx//2
    arrflip = np.roll(arrflip,offset*2,axis=1)
    
    if offset<0: arrflip[:,offset*2:]=0
    elif offset>0: arrflip[:,:2*offset]=0
    
    arrsym = (arr + arrflip)/2
    
    #~ plotc.figure()
    #~ plotc.colorplot(arr,sp=311,x=x,y=z,yr=[-10,None],xr=[src_x-10,rec_x+10],colorbar=False)
    #~ plotc.colorplot(arrflip,sp=312,x=x,y=z,yr=[-10,None],xr=[src_x-10,rec_x+10],colorbar=False)
    #~ plotc.colorplot(arrsym,sp=313,x=x,y=z,yr=[-10,None],xr=[src_x-10,rec_x+10],colorbar=False)
    #~ plotc.show()
    return arrsym
    
    
    
def sym0(arr): return (arr+arr[:,::-1])/2
    

plt.figure()
gl=plotc.gridlist(2,1)

#~ plotc.colorplot(kvx/abs(kvx).max(),x=x,y=z,sp=next(gl),centerzero=True,title="Kvx",
#~ yr=[-15,None],xr=[src_x-10,rec_x+10],
#~ colorbar_properties={"orientation":"horizontal"})
#~ plotc.draw_vlines([src_x,rec_x],linestyle='dashed')
#~ plt.gca().set_aspect(1)

kvx_filtered=sym(kvx_filtered)

ax1,_=plotc.colorplot(kvx_filtered/abs(kvx_filtered).max(),x=x,y=z,sp=next(gl),centerzero=True,
title=r"$\mathrm{K}{vx}$",
yr=[-4,None],xr=[src_x-10,rec_x+10],vmax=0.6,colorbar=False,
axes_properties={'hide_xticklabels':True})

plotc.draw_vlines([src_x,rec_x],linestyle='dashed')
plt.gca().set_aspect(1)



#~ plt.figure()
#~ gl=plotc.gridlist(2,1)
#~ plotc.colorplot(kvz/abs(kvz).max(),sp=next(gl),x=x,y=z,centerzero=True,title="Kvz",
#~ yr=[-15,None],xr=[src_x-10,rec_x+10],
#~ colorbar_properties={"orientation":"horizontal"})
#~ plotc.draw_vlines([src_x,rec_x],linestyle='dashed')
#~ plt.gca().set_aspect(1)

kvz_filtered=sym(kvz_filtered)

ax2,_=plotc.colorplot(kvz_filtered/kvz_filtered.max(),x=x,y=z,sp=next(gl),centerzero=True,
title="$\mathrm{K}{vz}$",yr=[-4,None],xr=[src_x-10,rec_x+10],vmax=0.6,
colorbar_properties={"orientation":"horizontal","pad":0.3})

plotc.draw_vlines([src_x,rec_x],linestyle='dashed')

plt.xlabel("Distance (Mm)",fontsize=14)
plt.gca().set_aspect(1)

plt.subplots_adjust(hspace=0)
plt.tight_layout()

plt.show()
