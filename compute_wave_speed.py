from __future__ import division
import numpy as np
import pyfits as pyfits
import plotc
import matplotlib.pyplot as plt
from mpldatacursor import datacursor
from scipy.optimize import curve_fit
from matplotlib.ticker import MaxNLocator
import os,sys,re
import modefilters,read_params

codedir=os.path.dirname(os.path.abspath(__file__))
datadir=read_params.get_directory()

Lx=read_params.get_xlength()
nx=read_params.get_nx()
x=np.linspace(-Lx/2,Lx/2,nx,endpoint=False)
dt_sec=read_params.get_dt()
dt_min=dt_sec/60

nt=np.squeeze(pyfits.getdata(os.path.join(datadir,'forward_src01_ls00','data.fits'))).shape[0]
t=np.arange(nt)*dt_min
Tsol_min = t[-1]

master_pixels=np.loadtxt(os.path.join(datadir,'master.pixels'),ndmin=1)
Nsources=len(master_pixels)

wavespeed=np.zeros((Nsources,10))

modes_to_plot=dict(fmode=False)
for i in xrange(10): modes_to_plot.update({'p'+str(i)+'mode':False})

ridge_filters=read_params.get_modes_used()

modes={'0':'fmode'}
for pmodeno in xrange(1,10): modes.update({str(pmodeno):'p'+str(pmodeno)+'mode'})

#~ for m in ridge_filters: modes_to_plot[modes[m]]=True
modes_to_plot['fmode']=True

try:
    for m in sys.argv[1:]:
        if m.startswith("f"): modes_to_plot['fmode']=True
        for i in xrange(1,10):
            if m.startswith("p"+str(i)): modes_to_plot['p'+str(i)+'mode']=True
except IndexError: pass

#~ Line fit to bottom of wavepacket, coordinates in Mm and min
def bot_time(x,(x1,y1),(x2,y2),sourceloc,source_bot_y):
    return (y2-y1)/(x2-x1)*(x-sourceloc)+source_bot_y
    
#~ Line fit to top of wavepacket, coordinates in Mm and min
def top_time(x,(x1,y1),(x2,y2),sourceloc,source_top_y):
    return (y2-y1)/(x2-x1)*(x-sourceloc)+source_top_y
    


for sourceno in xrange(Nsources):
    sourceloc=master_pixels[sourceno]
    sourcedir='forward_src'+str(sourceno+1).zfill(2)+'_ls00'

    vzcc=np.squeeze(pyfits.getdata(os.path.join(datadir,sourcedir,'data.fits')))
    sourcepix=int(np.floor(sourceloc/Lx*nx)+nx/2)
    sourcepix_loc = x[sourcepix]
    
    #~ print "Source",sourceno+1,"located at",sourceloc,"Mm","pix",sourcepix
    
    colors = ['r','g','b','cyan','brown','skyblue']
    
    if modes_to_plot['fmode']:

        fmode_filter_path=os.path.join(codedir,'fmode_filter.fits')
        if os.path.exists(fmode_filter_path):
            fmode_filter=np.squeeze(pyfits.getdata(fmode_filter_path))
        else:
            fmode_filter=modefilters.fmode_filter(nt,dt_sec,nx,Lx)
            fmode_filter=np.squeeze(fmode_filter).T

        vzcc_f=np.fft.ifft2(np.fft.fft2(vzcc)*fmode_filter).real
        #~ plotc.colorplot(vzcc_f,x=x,y=t,centerzero=True,colorbar=False,sp=121,xr=[-Lx/4,0])
        #~ datacursor(display='multiple',draggable='True')
        #~ plt.show()
        #~ exit()

        Npoints=4
        x_pix=np.array([-30-ind*6 for ind in xrange(Npoints)])+sourcepix
        x_coord=(x_pix-nx//2+1)/nx*Lx

        #~ plotc.draw_vlines(x_coord)

        t_cutoff=np.zeros((Npoints,3))

        #~ Line fit to bottom of wavepacket, x in Mm, y in min
        #~ Just select points from the plot, absolute shifts by source loc doesn't matter in slope
        x1=-124;y1=228;
        x2=-35.3;y2=42.4;
        source_bot=-20
        t_cutoff[:,0]=bot_time(x_coord,(x1,y1),(x2,y2),sourceloc,source_bot)
        
        
        #~ Line fit to top of wavepacket, x in Mm, y in min
        x1=-65.5;y1=225;
        x2=-25.2;y2=93.9;
        source_top=-10
        t_cutoff[:,1]=top_time(x_coord,(x1,y1),(x2,y2),sourceloc,source_top)
        

        #~ plt.plot(x_coord,t_cutoff[:,0],'b-')
        #~ plt.plot(x_coord,t_cutoff[:,1],'g-')
        #~ plt.show()
        #~ exit()
        
        def local_max_index(arr):
            return (np.diff(np.sign(np.diff(arr))) < 0).nonzero()[0] + 1

        peakcenter=[]
        peakamplitude=[]
        peakwidth=[]

        def gaus(x,a,x0,sigma):
            return a*np.exp(-(x-x0)**2/(2*sigma**2))

        for xind,pix in enumerate(x_pix):
            tseries=vzcc_f[:,pix]
            t_pix=local_max_index(tseries)
            t_coord = t_pix*dt_min
            local_max_in_range=np.where((t_coord>t_cutoff[xind,0]) & (t_coord<t_cutoff[xind,1]))[0]
            t_pix=t_pix[local_max_in_range]
            t_coord=t_coord[local_max_in_range]
            wave_amp=tseries[t_pix]
            
            #~ plt.subplot(121)
            #~ plt.plot([x_coord[xind]]*len(t_coord),t_coord,'o',color=colors[xind])
            #~ plt.subplot(122)
            #~ plt.plot(t_coord,wave_amp,'o',color=colors[xind])
            #~ plt.gca().yaxis.set_major_locator(MaxNLocator(5,prune='both'))
            #~ plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
            #~ plt.show()
            #~ exit()
            
            popt,_ = curve_fit(gaus,t_coord,wave_amp,p0=[wave_amp[len(wave_amp)//2],t_coord[len(t_coord)//2],(t_coord[-1]-t_coord[0])*0.4])
            
            #~ gfit = gaus(t_coord,*popt)
            #~ plt.plot(t_coord,gfit,color=colors[xind])
            
            peakamplitude.append(popt[0])
            peakcenter.append(popt[1])
            peakwidth.append(popt[2])

        peakamplitude=np.array(peakamplitude)
        peakcenter=np.array(peakcenter)
        peakwidth=np.array(peakwidth)
        v,_=np.polyfit(peakcenter,x_coord,1)
        
        #~ print peakcenter,x_coord

        print "Source",sourceno+1,"velocity of f mode",abs(v),"Mm/min"
            
        #~ plt.show()
        #~ exit()
        #~ 
        wavespeed[sourceno,0]=abs(v)
    
    if modes_to_plot['p1mode']:
        
        pmode_filter_path=os.path.join(codedir,'p1mode_filter.fits')
        if os.path.exists(pmode_filter_path):
            pmode_filter=np.squeeze(pyfits.getdata(pmode_filter_path))
        else:
            pmode_filter=modefilters.pmode_filter(nt,dt_sec,nx,Lx)
            pmode_filter=np.squeeze(pmode_filter).T

        
        vzcc_p=np.fft.ifft2(np.fft.fft2(vzcc)*pmode_filter).real
        #~ plt.figure(1)
        #~ plotc.colorplot(vzcc_p,centerzero=True,colorbar=False)
        #~ plt.figure(2)
        #~ plotc.colorplot(vzcc_p,x=x,y=t,sp=121,xr=[-Lx/4,0],centerzero=True,colorbar=False)
        #~ datacursor(display='multiple',draggable='True')
        #~ plt.show()
        #~ exit()

        Npoints=4
        x_pix=np.array([-30-ind*6 for ind in xrange(Npoints)])+sourcepix
        x_coord=x[x_pix]
        #~ np.set_printoptions(threshold=np.nan)
        #~ print vzcc_p[140:200,214]
        
        #~ plt.figure(2)
        #~ plotc.draw_vlines(x_coord)
        #~ plt.figure(1)
        #~ plotc.draw_vlines(x_pix)

        t_cutoff=np.zeros((Npoints,2))

        #~ Line fit to bottom of wavepacket, x in Mm, y in min
        x1=-119;y1=174;
        x2=-40.6;y2=52.2;
        source_bot=-8
        t_cutoff[:,0]=bot_time(x_coord,(x1,y1),(x2,y2),sourcepix_loc,source_bot)
        
        #~ Line fit to top of wavepacket, x in Mm, y in min
        x1=-93.1;y1=186;
        x2=-34.6;y2=76.8;
        source_top=20
        t_cutoff[:,1]=top_time(x_coord,(x1,y1),(x2,y2),sourcepix_loc,source_top)
        
        #~ plt.figure(1)
        #~ plt.plot(x_pix,t_cutoff[:,0]/dt_min,'b-')
        #~ plt.plot(x_pix,t_cutoff[:,1]/dt_min,'g-')
        #~ plt.figure(2)
        #~ plt.plot(x_coord,t_cutoff[:,0],'b-')
        #~ plt.plot(x_coord,t_cutoff[:,1],'g-')
        #~ plt.show()
        #~ quit()

        def local_max_index(arr):
            return (np.diff(np.sign(np.diff(arr))) < 0).nonzero()[0] + 1

        peakcenter=[]
        peakamplitude=[]
        peakwidth=[]

        def gaus(x,a,x0,sigma):
            return a*np.exp(-(x-x0)**2/(2*sigma**2))

        for xind,pix in enumerate(x_pix):
            tseries = vzcc_p[:,pix]
            t_pix=local_max_index(tseries)
            t_coord = t_pix*dt_min
            local_max_in_range=np.where((t_coord>t_cutoff[xind,0]) & (t_coord<t_cutoff[xind,1]))[0]
            t_pix=t_pix[local_max_in_range]
            t_coord=t_coord[local_max_in_range]
            wave_amp=tseries[t_pix]
            
            #~ plt.subplot(121)
            #~ plt.plot([x_coord[xind]]*len(t_coord),t_coord,'o',color=colors[xind])
            #~ plt.subplot(122)
            #~ plt.plot(t_coord,wave_amp,'o',color=colors[xind])
            #~ plt.gca().yaxis.set_major_locator(MaxNLocator(5,prune='both'))
            #~ plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
            #~ plt.show()
            #~ exit()
            
            popt,_ = curve_fit(gaus,t_coord,wave_amp,p0=[wave_amp[len(wave_amp)//2],t_coord[len(t_coord)//2],(t_coord[-1]-t_coord[0])*0.4])
            
            gfit = gaus(t_coord,*popt)
            #~ plt.plot(t_coord,gfit,color=colors[xind])
            
            
            peakamplitude.append(popt[0])
            peakcenter.append(popt[1])
            peakwidth.append(popt[2])
            
        peakamplitude=np.array(peakamplitude)
        peakcenter=np.array(peakcenter)
        peakwidth=np.array(peakwidth)
            
        
            
        v,_=np.polyfit(peakcenter,x_coord,1)
        
        print "Source",sourceno+1,"velocity of p mode",abs(v),"Mm/min"
        #~ plt.show()
        #~ exit()
        wavespeed[sourceno,1]=abs(v)
    
    if modes_to_plot['p2mode']:
        
        p2mode_filter_path=os.path.join(codedir,'p2mode_filter.fits')
        if os.path.exists(p2mode_filter_path):
            p2mode_filter=np.squeeze(pyfits.getdata(p2mode_filter_path))
        else:
            p2mode_filter=modefilters.p2mode_filter(nt,dt_sec,nx,Lx)
            p2mode_filter=np.squeeze(p2mode_filter).T

        vzcc_p=np.real(np.fft.ifft2(np.fft.fft2(vzcc)*p2mode_filter))
        #~ plotc.colorplot(vzcc_p,x=x,y=t,centerzero=True,xr=[-Lx/5,0],sp=121,colorbar=False)
        #~ datacursor(display='multiple',draggable='True')
        #~ plt.show()
        #~ quit()


        Npoints=4
        x_pix=np.array([-30-ind*6 for ind in xrange(Npoints)])+sourcepix
        x_coord=(x_pix-nx//2+1)/nx*Lx

        #~ plotc.draw_vlines(x_coord)

        t_cutoff=np.zeros((Npoints,2))

        #~ Line fit to bottom of wavepacket, x in Mm, y in min
        x1=-53.7;y1=33.4;
        x2=-129;y2=111;
        source_bot=0
        t_cutoff[:,0]=bot_time(x_coord,(x1,y1),(x2,y2),sourceloc,source_bot)
        
        #~ Line fit to top of wavepacket, x in Mm, y in min
        x1=-102;y1=135;
        x2=-46.1;y2=56.8;
        source_top=16
        t_cutoff[:,1]=top_time(x_coord,(x1,y1),(x2,y2),sourceloc,source_top)
        
        #~ plt.plot(x_coord,t_cutoff[:,0],'b-')
        #~ plt.plot(x_coord,t_cutoff[:,1],'g-')
        #~ plt.show()
        #~ quit()

        xdat=vzcc_p[:,x_pix].T

        def local_max_index(arr):
            return (np.diff(np.sign(np.diff(arr))) < 0).nonzero()[0] + 1

        peakcenter=[]
        peakamplitude=[]
        peakwidth=[]

        def gaus(x,a,x0,sigma):
            return a*np.exp(-(x-x0)**2/(2*sigma**2))

        
        for xind,pix in enumerate(x_pix):
            tseries=vzcc_f[:,pix]
            t_pix=local_max_index(tseries)
            t_coord = t_pix*dt_min
            local_max_in_range=np.where((t_coord>t_cutoff[xind,0]) & (t_coord<t_cutoff[xind,1]))[0]
            t_pix=t_pix[local_max_in_range]
            t_coord=t_coord[local_max_in_range]
            wave_amp=tseries[t_pix]    

            #~ plt.subplot(121)
            #~ plt.plot([x_coord[xind]]*len(t_coord),t_coord,'o',color=colors[xind])
            #~ plt.subplot(122)
            #~ plt.plot(t_coord,wave_amp,'o',color=colors[xind])
            #~ plt.gca().yaxis.set_major_locator(MaxNLocator(5,prune='both'))
            #~ plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
            
            popt,_ = curve_fit(gaus,t_coord,wave_amp,p0=[wave_amp[len(wave_amp)//2],t_coord[len(t_coord)//2],(t_coord[-1]-t_coord[0])*0.4])
            
            #~ gfit = gaus(t_coord,*popt)            
            #~ plt.plot(t_coord,gfit,color=colors[xind])
            
            peakamplitude.append(popt[0])
            peakcenter.append(popt[1])
            peakwidth.append(popt[2])
            
        peakamplitude=np.array(peakamplitude)
        peakcenter=np.array(peakcenter)
        peakwidth=np.array(peakwidth)
        
        v,_=np.polyfit(peakcenter,x_coord,1)
        print "Source",sourceno+1,"velocity of p2 mode",abs(v),"Mm/min"
        #~ plt.show()
        #~ exit()
        wavespeed[sourceno,2]=abs(v)

    if modes_to_plot['p3mode']:
        
        p3mode_filter_path=os.path.join(codedir,'p3mode_filter.fits')
        if os.path.exists(p3mode_filter_path):
            p3mode_filter=np.squeeze(pyfits.getdata(p3mode_filter_path))
        else:
            p3mode_filter=modefilters.p3mode_filter(nt,dt_sec,nx,Lx)
            p3mode_filter=np.squeeze(p3mode_filter).T

        vzcc_p=np.real(np.fft.ifft2(np.fft.fft2(vzcc)*p3mode_filter))
        #~ plotc.colorplot(vzcc_p,x=x,y=t,centerzero=True,sp=121,colorbar=False,xr=[-Lx/5,0])
        #~ datacursor(display='multiple',draggable='True')
        #~ plt.show()
        #~ quit()


        Npoints=4
        x_pix=np.array([-30-ind*6 for ind in xrange(Npoints)])+sourcepix
        x_coord=(x_pix-nx//2+1)/nx*Lx

        #~ plotc.draw_vlines(x_coord)

        t_cutoff=np.zeros((Npoints,2))

        #~ Line fit to bottom of wavepacket, x in Mm, y in min
        x1=-147;y1=112;
        x2=-60;y2=35;
        source_bot=0
        t_cutoff[:,0]=bot_time(x_coord,(x1,y1),(x2,y2),sourceloc,source_bot)
        
        #~ Line fit to top of wavepacket, x in Mm, y in min
        x1=-99.7;y1=126;
        x2=-41.5;y2=50.5;
        source_top=20
        t_cutoff[:,1]=top_time(x_coord,(x1,y1),(x2,y2),sourceloc,source_top)

        #~ plt.plot(x_coord,t_cutoff[:,0],'b-')
        #~ plt.plot(x_coord,t_cutoff[:,1],'g-')
        #~ plt.show()
        #~ quit()

        xdat=vzcc_p[:,x_pix].T

        def local_max_index(arr):
            return (np.diff(np.sign(np.diff(arr))) < 0).nonzero()[0] + 1

        peakcenter=[]
        peakamplitude=[]
        peakwidth=[]

        def gaus(x,a,x0,sigma):
            return a*np.exp(-(x-x0)**2/(2*sigma**2))

        for xind,pix in enumerate(x_pix):
            tseries=vzcc_f[:,pix]
            t_pix=local_max_index(tseries)
            t_coord = t_pix*dt_min
            local_max_in_range=np.where((t_coord>t_cutoff[xind,0]) & (t_coord<t_cutoff[xind,1]))[0]
            t_pix=t_pix[local_max_in_range]
            t_coord=t_coord[local_max_in_range]
            wave_amp=tseries[t_pix]    

            #~ plt.subplot(121)
            #~ plt.plot([x_coord[xind]]*len(t_coord),t_coord,'o',color=colors[xind])
            #~ plt.subplot(122)
            #~ plt.plot(t_coord,wave_amp,'o',color=colors[xind])
            #~ plt.gca().yaxis.set_major_locator(MaxNLocator(5,prune='both'))
            #~ plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
            
            popt,_ = curve_fit(gaus,t_coord,wave_amp,p0=[wave_amp[len(wave_amp)//2],t_coord[len(t_coord)//2],(t_coord[-1]-t_coord[0])*0.4])
            
            #~ gfit = gaus(t_coord,*popt)            
            #~ plt.plot(t_coord,gfit,color=colors[xind])
            
            peakamplitude.append(popt[0])
            peakcenter.append(popt[1])
            peakwidth.append(popt[2])
            
        peakamplitude=np.array(peakamplitude)
        peakcenter=np.array(peakcenter)
        peakwidth=np.array(peakwidth)
        
        v,_=np.polyfit(peakcenter,x_coord,1)
        print "Source",sourceno+1,"velocity of p3 mode",abs(v),"Mm/min"
        #~ plt.show()
        #~ exit()
        wavespeed[sourceno,3]=abs(v)

    if modes_to_plot['p4mode']:
        
        p4mode_filter_path=os.path.join(codedir,'p4mode_filter.fits')
        if os.path.exists(p4mode_filter_path):
            p4mode_filter=np.squeeze(pyfits.getdata(p4mode_filter_path))
        else:
            p4mode_filter=modefilters.p4mode_filter(nt,dt_sec,nx,Lx)
            p4mode_filter=np.squeeze(p4mode_filter).T

        vzcc_p=np.real(np.fft.ifft2(np.fft.fft2(vzcc)*p4mode_filter))
        plotc.colorplot(vzcc_p,x=x,y=t,centerzero=True,colorbar=False,sp=121,xr=[-Lx/5,0])
        datacursor(display='multiple',draggable='True')
        #~ plt.show()
        #~ quit()


        Npoints=4
        x_pix=np.array([-30-ind*6 for ind in xrange(Npoints)])+sourcepix
        x_coord=(x_pix-nx//2+1)/nx*Lx

        plotc.draw_vlines(x_coord)

        t_cutoff=np.zeros((Npoints,2))

        #~ Line fit to bottom of wavepacket, x in Mm, y in min
        x1=-60.6;y1=28.2;
        x2=-152;y2=85.7;
        source_bot=3
        t_cutoff[:,0]=bot_time(x_coord,(x1,y1),(x2,y2),sourceloc,source_bot)
        
        #~ Line fit to top of wavepacket, x in Mm, y in min
        x1=-132;y1=127;
        x2=-56.4;y2=56.6;
        source_top=20
        t_cutoff[:,1]=top_time(x_coord,(x1,y1),(x2,y2),sourceloc,source_top)

        plt.plot(x_coord,t_cutoff[:,0],'b-')
        plt.plot(x_coord,t_cutoff[:,1],'g-')
        #~ plt.show()
        #~ quit()

        #~ xdat=vzcc_p[:,x_pix].T

        def local_max_index(arr):
            return (np.diff(np.sign(np.diff(arr))) < 0).nonzero()[0] + 1

        peakcenter=[]
        peakamplitude=[]
        peakwidth=[]

        def gaus(x,a,x0,sigma):
            return a*np.exp(-(x-x0)**2/(2*sigma**2))


        for xind,pix in enumerate(x_pix):
            tseries = vzcc_p[:,pix]
            t_pix=local_max_index(tseries)
            t_coord = t_pix*dt_min
            local_max_in_range=np.where((t_coord>t_cutoff[xind,0]) & (t_coord<t_cutoff[xind,1]))[0]
            t_pix=t_pix[local_max_in_range]
            t_coord=t_coord[local_max_in_range]
            wave_amp=tseries[t_pix]    

            plt.subplot(121)
            plt.plot([x_coord[xind]]*len(t_coord),t_coord,'o',color=colors[xind])
            plt.subplot(122)
            plt.plot(t_coord,wave_amp,'o',color=colors[xind])
            plt.gca().yaxis.set_major_locator(MaxNLocator(5,prune='both'))
            plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
            
            popt,_ = curve_fit(gaus,t_coord,wave_amp,p0=[wave_amp[len(wave_amp)//2],t_coord[len(t_coord)//2],(t_coord[-1]-t_coord[0])*0.4])
            
            gfit = gaus(t_coord,*popt)            
            plt.plot(t_coord,gfit,color=colors[xind])
            
            peakamplitude.append(popt[0])
            peakcenter.append(popt[1])
            peakwidth.append(popt[2])
            
        peakamplitude=np.array(peakamplitude)
        peakcenter=np.array(peakcenter)
        peakwidth=np.array(peakwidth)
        

            
        v,_=np.polyfit(peakcenter,x_coord,1)
        print "Source",sourceno+1,"velocity of p4 mode",abs(v),"Mm/min"
        #~ plt.show()
        #~ exit()
        wavespeed[sourceno,4]=abs(v)

    if modes_to_plot['p5mode']:
        
        p5mode_filter_path=os.path.join(codedir,'p5mode_filter.fits')
        if os.path.exists(p5mode_filter_path):
            p5mode_filter=np.squeeze(pyfits.getdata(p5mode_filter_path))
        else:
            p5mode_filter=modefilters.p5mode_filter(nt,dt_sec,nx,Lx)
            p5mode_filter=np.squeeze(p5mode_filter).T

        vzcc_p=np.real(np.fft.ifft2(np.fft.fft2(vzcc)*p5mode_filter))
        #~ plotc.colorplot(vzcc_p,x=x-sourcepix_loc,y=t,centerzero=True,sp=121,colorbar=False,xr=[-Lx/3,0])
        #~ datacursor(display='multiple',draggable='True')
        #~ plt.show()
        #~ quit()


        Npoints=6
        x_pix=np.array([-80-ind*12 for ind in xrange(Npoints)])+sourcepix
        x_coord=(x_pix-nx//2+1)/nx*Lx

        #~ plotc.draw_vlines(x_coord)

        t_cutoff=np.zeros((Npoints,2))

        #~ Line fit to bottom of wavepacket, x in Mm, y in min
        x1=-269;y1=155;
        x2=-61.5;y2=35;
        source_bot=4
        t_cutoff[:,0]=bot_time(x_coord,(x1,y1),(x2,y2),sourceloc,source_bot)
        
        #~ Line fit to top of wavepacket, x in Mm, y in min
        x1=-188;y1=171;
        x2=-29.2;y2=41.7;
        source_top=28
        t_cutoff[:,1]=top_time(x_coord,(x1,y1),(x2,y2),sourceloc,source_top)

        #~ plt.plot(x_coord,t_cutoff[:,0],'b-')
        #~ plt.plot(x_coord,t_cutoff[:,1],'g-')
        #~ plt.show()
        #~ quit()

        def local_max_index(arr):
            return (np.diff(np.sign(np.diff(arr))) < 0).nonzero()[0] + 1

        peakcenter=[]
        peakamplitude=[]
        peakwidth=[]

        def gaus(x,a,x0,sigma):
            return a*np.exp(-(x-x0)**2/(2*sigma**2))


        for xind,pix in enumerate(x_pix):
            tseries = vzcc_p[:,pix]
            t_pix=local_max_index(tseries)
            t_coord = t_pix*dt_min
            local_max_in_range=np.where((t_coord>t_cutoff[xind,0]) & (t_coord<t_cutoff[xind,1]))[0]
            t_pix=t_pix[local_max_in_range]
            t_coord=t_coord[local_max_in_range]
            wave_amp=tseries[t_pix]    

            #~ plt.subplot(121)
            #~ plt.plot([x_coord[xind]]*len(t_coord),t_coord,'o',color=colors[xind])
            #~ plt.subplot(122)
            #~ plt.plot(t_coord,wave_amp,'o',color=colors[xind])
            #~ plt.gca().yaxis.set_major_locator(MaxNLocator(5,prune='both'))
            #~ plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
            
            popt,_ = curve_fit(gaus,t_coord,wave_amp,p0=[wave_amp[len(wave_amp)//2],t_coord[len(t_coord)//2],(t_coord[-1]-t_coord[0])*0.4])
            
            #~ gfit = gaus(t_coord,*popt)            
            #~ plt.plot(t_coord,gfit,color=colors[xind])
            
            peakamplitude.append(popt[0])
            peakcenter.append(popt[1])
            peakwidth.append(popt[2])
            
        peakamplitude=np.array(peakamplitude)
        peakcenter=np.array(peakcenter)
        peakwidth=np.array(peakwidth)
        
            
        v,_=np.polyfit(peakcenter,x_coord,1)
        print "Source",sourceno+1,"velocity of p5 mode",abs(v),"Mm/min"
        #~ plt.show()
        #~ exit()
        wavespeed[sourceno,5]=abs(v)

    if modes_to_plot['p6mode']:
        
        p6mode_filter_path=os.path.join(codedir,'p6mode_filter.fits')
        if os.path.exists(p6mode_filter_path):
            p6mode_filter=np.squeeze(pyfits.getdata(p6mode_filter_path))
        else:
            p6mode_filter=modefilters.p6mode_filter(nt,dt_sec,nx,Lx)
            p6mode_filter=np.squeeze(p6mode_filter).T

        vzcc_p=np.real(np.fft.ifft2(np.fft.fft2(vzcc)*p6mode_filter))
        #~ plotc.colorplot(vzcc_p,x=x,y=t,centerzero=True,sp=121,colorbar=False,xr=[-Lx/2,0])
        #~ datacursor(display='multiple',draggable='True')
        #~ plt.show()
        #~ quit()


        Npoints=6
        x_pix=np.array([-80-ind*12 for ind in xrange(Npoints)])+sourcepix
        x_coord=x[x_pix]

        #~ plotc.draw_vlines(x_coord)

        t_cutoff=np.zeros((Npoints,2))

        #~ Line fit to bottom of wavepacket, x in Mm, y in min
        x1=-71.8;y1=28.2;
        x2=-324;y2=151;
        source_bot=4
        t_cutoff[:,0]=bot_time(x_coord,(x1,y1),(x2,y2),sourceloc,source_bot)
        
        #~ Line fit to top of wavepacket, x in Mm, y in min
        x1=-226;y1=171;
        x2=-59.4;y2=47.1;
        source_top=30
        t_cutoff[:,1]=top_time(x_coord,(x1,y1),(x2,y2),sourceloc,source_top)

        #~ plt.plot(x_coord,t_cutoff[:,0],'b-')
        #~ plt.plot(x_coord,t_cutoff[:,1],'g-')
        #~ plt.show()
        #~ quit()

        xdat=vzcc_p[:,x_pix].T

        def local_max_index(arr):
            return (np.diff(np.sign(np.diff(arr))) < 0).nonzero()[0] + 1

        peakcenter=[]
        peakamplitude=[]
        peakwidth=[]

        def gaus(x,a,x0,sigma):
            return a*np.exp(-(x-x0)**2/(2*sigma**2))


        for xind,tseries in enumerate(xdat):
            t_pix=local_max_index(tseries)
            t_coord = t_pix*dt_min
            local_max_in_range=np.where((t_coord>t_cutoff[xind,0]) & (t_coord<t_cutoff[xind,1]))[0]
            t_pix=t_pix[local_max_in_range]
            t_coord=t_coord[local_max_in_range]
            wave_amp=tseries[t_pix]    

            #~ plt.subplot(121)
            #~ plt.plot([x_coord[xind]]*len(t_coord),t_coord,'o',color=colors[xind])
            #~ plt.subplot(122)
            #~ plt.plot(t_coord,wave_amp,'o',color=colors[xind])
            #~ plt.gca().yaxis.set_major_locator(MaxNLocator(5,prune='both'))
            #~ plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
            
            popt,_ = curve_fit(gaus,t_coord,wave_amp,p0=[wave_amp[len(wave_amp)//2],t_coord[len(t_coord)//2],(t_coord[-1]-t_coord[0])*0.4])
            
            #~ gfit = gaus(t_coord,*popt)            
            #~ plt.plot(t_coord,gfit,color=colors[xind])
            
            peakamplitude.append(popt[0])
            peakcenter.append(popt[1])
            peakwidth.append(popt[2])
            
        peakamplitude=np.array(peakamplitude)
        peakcenter=np.array(peakcenter)
        peakwidth=np.array(peakwidth)
        
        v,_=np.polyfit(peakcenter,x_coord,1)
        print "Source",sourceno+1,"velocity of p6 mode",abs(v),"Mm/min"
        #~ plt.show()
        #~ exit()
        wavespeed[sourceno,6]=abs(v)

    if modes_to_plot['p7mode']:
        
        p7mode_filter_path=os.path.join(codedir,'p7mode_filter.fits')
        if os.path.exists(p7mode_filter_path):
            p7mode_filter=np.squeeze(pyfits.getdata(p7mode_filter_path))
        else:
            p7mode_filter=modefilters.p7mode_filter(nt,dt_sec,nx,Lx)
            p7mode_filter=np.squeeze(p7mode_filter).T

        vzcc_p=np.real(np.fft.ifft2(np.fft.fft2(vzcc)*p7mode_filter))
        plotc.colorplot(vzcc_p,x=x,y=t,centerzero=True,sp=121,colorbar=False,xr=[-Lx/2,0])
        datacursor(display='multiple',draggable='True')
        #~ plt.show()
        #~ quit()


        Npoints=6
        x_pix=np.array([-80-ind*12 for ind in xrange(Npoints)])+sourcepix
        x_coord=(x_pix-nx//2+1)/nx*Lx

        plotc.draw_vlines(x_coord)

        t_cutoff=np.zeros((Npoints,2))

        #~ Line fit to bottom of wavepacket, x in Mm, y in min
        x1=-324;y1=141;
        x2=-109;y2=42.4;
        source_bot=0
        t_cutoff[:,0]=bot_time(x_coord,(x1,y1),(x2,y2),sourceloc,source_bot)
        
        #~ Line fit to top of wavepacket, x in Mm, y in min
        x1=-226;y1=156;
        x2=-68.2;y2=53.2;
        source_top=25
        t_cutoff[:,1]=top_time(x_coord,(x1,y1),(x2,y2),sourceloc,source_top)

        plt.plot(x_coord,t_cutoff[:,0],'b-')
        plt.plot(x_coord,t_cutoff[:,1],'g-')
        #~ plt.show()
        #~ quit()

        xdat=vzcc_p[:,x_pix].T

        def local_max_index(arr):
            return (np.diff(np.sign(np.diff(arr))) < 0).nonzero()[0] + 1

        peakcenter=[]
        peakamplitude=[]
        peakwidth=[]

        def gaus(x,a,x0,sigma):
            return a*np.exp(-(x-x0)**2/(2*sigma**2))


        for xind,tseries in enumerate(xdat):
            t_pix=local_max_index(tseries)
            t_coord = t_pix*dt_min
            local_max_in_range=np.where((t_coord>t_cutoff[xind,0]) & (t_coord<t_cutoff[xind,1]))[0]
            t_pix=t_pix[local_max_in_range]
            t_coord=t_coord[local_max_in_range]
            wave_amp=tseries[t_pix]    

            plt.subplot(121)
            plt.plot([x_coord[xind]]*len(t_coord),t_coord,'o',color=colors[xind])
            plt.subplot(122)
            plt.plot(t_coord,wave_amp,'o',color=colors[xind])
            plt.gca().yaxis.set_major_locator(MaxNLocator(5,prune='both'))
            plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
            
            popt,_ = curve_fit(gaus,t_coord,wave_amp,p0=[wave_amp[len(wave_amp)//2],t_coord[len(t_coord)//2],(t_coord[-1]-t_coord[0])*0.4])
            
            gfit = gaus(t_coord,*popt)            
            plt.plot(t_coord,gfit,color=colors[xind])
            
            peakamplitude.append(popt[0])
            peakcenter.append(popt[1])
            peakwidth.append(popt[2])
            
        peakamplitude=np.array(peakamplitude)
        peakcenter=np.array(peakcenter)
        peakwidth=np.array(peakwidth)
            
        v,_=np.polyfit(peakcenter,x_coord,1)
        print "Source",sourceno+1,"velocity of p7 mode",abs(v),"Mm/min"
        plt.show()
        exit()
        wavespeed[sourceno,7]=abs(v)


np.savetxt('wavespeed',wavespeed,fmt="%14.8f")
