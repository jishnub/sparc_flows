
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import plotc
import read_params
import os,pyfits

def fitsread(f): return np.squeeze(pyfits.getdata(f))
    
datadir=read_params.get_directory()
codedir=os.path.dirname(os.path.abspath(__file__))

if not os.path.exists(os.path.join(datadir,'plots/vx')):
    os.makedirs(os.path.join(datadir,'plots/vx'))
if not os.path.exists(os.path.join(datadir,'plots/vz')):
    os.makedirs(os.path.join(datadir,'plots/vz'))
if not os.path.exists(os.path.join(datadir,'plots/xiz')):
    os.makedirs(os.path.join(datadir,'plots/xiz'))
if not os.path.exists(os.path.join(datadir,'plots/xix')):
    os.makedirs(os.path.join(datadir,'plots/xix'))
if not os.path.exists(os.path.join(datadir,'plots/accz')):
    os.makedirs(os.path.join(datadir,'plots/accz'))
if not os.path.exists(os.path.join(datadir,'plots/accx')):
    os.makedirs(os.path.join(datadir,'plots/accx'))

Lx=read_params.get_xlength()
nx=read_params.get_nx()
x=np.linspace(-Lx/2,Lx/2,nx,endpoint=False)

z=np.loadtxt(os.path.join(codedir,read_params.get_solarmodel()),usecols=[0])
z=(z-1)*695.8

excite_depth=read_params.get_excitedepth()/1e8
obs_depth = read_params.get_obs_depth()/1e8

src=1
src_locs=np.loadtxt(os.path.join(datadir,'master.pixels'),ndmin=1)
src_x = src_locs[src-1]
rec_x=18

for i in range(280):
    print(i, end=' ')
    tstring = str(i*15).zfill(6)
    try:
        data=fitsread(os.path.join(datadir,'forward_src'+str(src).zfill(2)+'_ls00','vx_'+tstring+'_partial.fits'))
        plotc.colorplot(data/(data.max()*20),centerzero=True,x=x,y=z,yr=[-5,None])
        
        plt.xlabel("Horizontal Distance (Mm)",fontsize=20)
        plt.ylabel("Depth (Mm)",fontsize=20)

        plt.plot(x,[0]*len(x),linestyle='dashed',color='#666666')
        plt.plot([src_x],[excite_depth],'o',color='sandybrown',markersize=10)
        plt.plot([rec_x],[obs_depth],'o',color='sage',markersize=10)
        
        plt.savefig(os.path.join(datadir,'plots/vx/vx_'+str(i).zfill(3)+'.png'))
        plt.clf()
    except IOError: pass

    try:    
        data=fitsread(os.path.join(datadir,'forward_src'+str(src).zfill(2)+'_ls00','vz_'+tstring+'_partial.fits'))
        plotc.colorplot(data/(data.max()*20),centerzero=True,x=x,y=z,yr=[-5,None])
        
        plt.xlabel("Horizontal Distance (Mm)",fontsize=20)
        plt.ylabel("Depth (Mm)",fontsize=20)
        
        plt.plot(x,[0]*len(x),linestyle='dashed',color='#666666')
        plt.plot([src_x],[excite_depth],'o',color='sandybrown',markersize=10)
        plt.plot([rec_x],[obs_depth],'o',color='sage',markersize=10)
        
        plt.savefig(os.path.join(datadir,'plots/vz/vz_'+str(i).zfill(3)+'.png'))
        plt.clf()
    except IOError: pass
#~ 
    #~ try:    
        #~ data=fitsread(os.path.join(datadir,'forward_src'+str(src).zfill(2)+'_ls00','xix_'+tstring+'_partial.fits'))
        #~ plotc.colorplot(data,centerzero=True,x=x,y=z,yr=[-5,None])
        #~ 
        #~ plt.plot(x,[0]*len(x),linestyle='dashed',color='#666666')
        #~ plt.plot([src_x],[excite_depth],'o',color='sandybrown',markersize=10)
        #~ plt.plot([rec_x],[obs_depth],'o',color='sage',markersize=10)
        
        #~ plt.savefig(os.path.join(datadir,'plots/xix/xix_'+str(i).zfill(3)+'.png'))
        #~ plt.clf()
    #~ except IOError: pass
    #~ 
    #~ try:
        #~ data=fitsread(os.path.join(datadir,'forward_src'+str(src).zfill(2)+'_ls00','xiz_'+tstring+'_partial.fits'))
        #~ plotc.colorplot(data,centerzero=True,x=x,y=z,yr=[-5,None])
        #~ 
        #~ plt.plot(x,[0]*len(x),linestyle='dashed',color='#666666')
        #~ plt.plot([src_x],[excite_depth],'o',color='sandybrown',markersize=10)
        #~ plt.plot([rec_x],[obs_depth],'o',color='sage',markersize=10)
        #~ 
        #~ plt.savefig(os.path.join(datadir,'plots/xiz/xiz_'+str(i).zfill(3)+'.png'))
        #~ plt.clf()
    #~ except IOError: pass
#~ 
    #~ try:
        #~ data=fitsread(os.path.join(datadir,'forward_src'+str(src).zfill(2)+'_ls00','acc_x_'+tstring+'_partial.fits'))
        #~ plotc.colorplot(data,centerzero=True,x=x,y=z,yr=[-5,None])
        #~ 
        #~ plt.plot(x,[0]*len(x),linestyle='dashed',color='#666666')
        #~ plt.plot([src_x],[excite_depth],'o',color='sandybrown',markersize=10)
        #~ plt.plot([rec_x],[obs_depth],'o',color='sage',markersize=10)
        
        #~ plt.savefig(os.path.join(datadir,'plots/accx/accx_'+str(i).zfill(3)+'.png'))
        #~ plt.clf()
    #~ except IOError: pass
    #~ 
    #~ try:
        #~ data=fitsread(os.path.join(datadir,'forward_src'+str(src).zfill(2)+'_ls00','acc_z_'+tstring+'_partial.fits'))
        #~ plotc.colorplot(data,centerzero=True,x=x,y=z,yr=[-5,None])
        #~ 
        #~ plt.plot(x,[0]*len(x),linestyle='dashed',color='#666666')
        #~ plt.plot([src_x],[excite_depth],'o',color='sandybrown',markersize=10)
        #~ plt.plot([rec_x],[obs_depth],'o',color='sage',markersize=10)
        #~ 
        #~ plt.savefig(os.path.join(datadir,'plots/accz/accz_'+str(i).zfill(3)+'.png'))
        #~ plt.clf()
    #~ except IOError: pass
    
    
    
