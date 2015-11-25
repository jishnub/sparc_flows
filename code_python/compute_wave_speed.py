from __future__ import division
import numpy as np
import pyfits as pf
import plotc
from mpldatacursor import datacursor
from scipy.optimize import curve_fit
import os
import modefilters,read_params

codedir=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

Lx=read_params.get_xlength()
dt=read_params.get_dt()/60
datadir=read_params.get_directory()

master_pixels=np.loadtxt(os.path.join(datadir,'master.pixels'),ndmin=1)
Nsources=len(master_pixels)

wavespeed=np.zeros((Nsources,3))

fmode=False
p1mode=False
p2mode=True

for sourceno in xrange(Nsources):
    sourceloc=master_pixels[sourceno]
    sourcedir='forward_src'+str(sourceno+1).zfill(2)+'_ls00'

    vzcc=np.squeeze(pf.getdata(os.path.join(datadir,sourcedir,'vz_cc.fits')))
    Nt,Nx=vzcc.shape
    sourcepix=int(np.floor(sourceloc/Lx*Nx)+Nx/2-1)
    
    
    if fmode:
        ####################################################################
        #~ f-mode

        fmode_filter_path=os.path.join(codedir,'fmode_filter.fits')
        if os.path.exists(fmode_filter_path):
            fmode_filter=np.squeeze(pf.getdata(fmode_filter_path))
        else:
            fmode_filter=np.asfortranarray(np.zeros((Nx,1,Nt),dtype=float))
            modefilters.fmode_filter(fmode_filter)
            fmode_filter=np.squeeze(fmode_filter).T

        vzcc_f=np.real(np.fft.ifft2(np.fft.fft2(vzcc)*fmode_filter))
        vzcc_f/=abs(vzcc_f).max()
        #~ plotc.colorplot(vzcc_f)
        #~ datacursor(display='multiple',draggable='True')


        Npoints=4
        xpix=np.array([-30-ind*6 for ind in xrange(Npoints)])+sourcepix
        xpoints=(xpix-Nx//2+1)/Nx*Lx

        #~ plotc.draw_vlines(xpix)

        t_cutoff=np.zeros((Npoints,3))

        #~ Line fit to bottom of wavepacket
        def bot_time(x):
            return (332-65.7)/(210-258)*(x-sourcepix)+20

        #~ Line fit to top of wavepacket
        def top_time(x):
            return (362-268)/(228-241)*(x-sourcepix)+60


        t_cutoff[:,0]=bot_time(xpix)
        t_cutoff[:,1]=top_time(xpix)
        
        y= (-0.465*(xpix-sourcepix)/dt+240)
        #y=-5*(t-240)+230
        

        #~ plotc.plt.plot(xpix,t_cutoff[:,0],'b-')
        #~ plotc.plt.plot(xpix,t_cutoff[:,1],'g-')
        #~ plotc.plt.plot(xpix,y,'r-') 
        #~ plotc.plt.show()
        
        xdat=vzcc_f[:,xpix].T

        def local_max_index(arr):
            return (np.diff(np.sign(np.diff(arr))) < 0).nonzero()[0] + 1

        peakcenter=[]
        peakamplitude=[]
        peakwidth=[]

        def gaus(x,a,x0,sigma):
            return a*np.exp(-(x-x0)**2/(2*sigma**2))


        for xind,tseries in enumerate(xdat):
            tp=local_max_index(tseries)
            tp=tp[np.where(np.logical_and(tp<t_cutoff[xind,1],tp>t_cutoff[xind,0]))[0]]
            yp=tseries[tp]    
            
            #plotc.plt.plot([xpix[xind]]*len(tp),tp,'o',color='orange')

            popt,_ = curve_fit(gaus,tp,yp,p0=[yp[len(yp)//2],tp[len(tp)//2],(tp[-1]-tp[0])*0.4])
            
            peakamplitude.append(popt[0])
            peakcenter.append(popt[1])
            peakwidth.append(popt[2])

        peakamplitude=np.array(peakamplitude)
        peakcenter=np.array(peakcenter)*dt
        peakwidth=np.array(peakwidth)*dt
        v,_=np.polyfit(peakcenter,xpoints,1)
        
        #~ print peakcenter,xpoints

        print "Source",sourceno+1,"velocity of f mode",abs(v),"Mm/min"
        
        #plotc.plt.figure()
        for xind,tseries in enumerate(xdat):
            tp=local_max_index(tseries)
            tp=tp[np.where(np.logical_and(tp<t_cutoff[xind,1],tp>t_cutoff[xind,0]))[0]]
            yp=tseries[tp]    
            #~ plotc.plt.plot(tp,yp,'o')
            gfit=gaus(tp,peakamplitude[xind],peakcenter[xind]/dt,peakwidth[xind]/dt)
            #print gfit
            #~ plotc.plt.plot(tp,gfit)
            
        #~ plotc.plt.show()
        
        
        wavespeed[sourceno,0]=abs(v)
    
    if p1mode:
        ####################################################################
        #~ p1-mode
        
        pmode_filter_path=os.path.join(codedir,'pmode_filter.fits')
        if os.path.exists(pmode_filter_path):
            pmode_filter=np.squeeze(pf.getdata(pmode_filter_path))
        else:
            pmode_filter=np.asfortranarray(np.zeros((Nx,1,Nt),dtype=float))
            modefilters.pmode_filter(pmode_filter)
            pmode_filter=np.squeeze(pmode_filter).T

        vzcc_p=np.real(np.fft.ifft2(np.fft.fft2(vzcc)*pmode_filter))
        vzcc_p/=abs(vzcc_p).max()
        #~ plotc.colorplot(vzcc_p,vmax=0.5,centerzero=True)
        

        #~ datacursor(display='multiple',draggable='True')
        #plotc.plt.show()
        #quit()


        Npoints=4
        xpix=np.array([-30-ind*6 for ind in xrange(Npoints)])+sourcepix
        xpoints=(xpix-Nx//2+1)/Nx*Lx

        #~ plotc.draw_vlines(xpix)

        t_cutoff=np.zeros((Npoints,2))

        #~ Line fit to bottom of wavepacket
        def bot_time(x):
            return (202-80.8)/(199-228)*(x-sourcepix)

        #~ Line fit to top of wavepacket
        def top_time(x):
            return (340-235)/(215-232)*(x-sourcepix)+58

        t_cutoff[:,0]=bot_time(xpix)
        t_cutoff[:,1]=top_time(xpix)
        
        #~ plotc.plt.plot(xpix,t_cutoff[:,0],'b-')
        #~ plotc.plt.plot(xpix,t_cutoff[:,1],'g-')
        #plotc.plt.show()
        #quit()

        xdat=vzcc_p[:,xpix].T

        def local_max_index(arr):
            return (np.diff(np.sign(np.diff(arr))) < 0).nonzero()[0] + 1

        peakcenter=[]
        peakamplitude=[]
        peakwidth=[]

        def gaus(x,a,x0,sigma):
            return a*np.exp(-(x-x0)**2/(2*sigma**2))


        for xind,tseries in enumerate(xdat):
            tp=local_max_index(tseries)
            tp=tp[np.where(np.logical_and(tp<t_cutoff[xind,1],tp>t_cutoff[xind,0]))[0]]
            yp=tseries[tp]    

            #~ plotc.plot([xpix[xind]]*len(tp),tp,'o',color='orange')
            #~ plotc.show()
            
            popt,_ = curve_fit(gaus,tp,yp,p0=[yp[len(yp)//2],tp[len(tp)//2],(tp[-1]-tp[0])*0.4])
            
            peakamplitude.append(popt[0])
            peakcenter.append(popt[1])
            peakwidth.append(popt[2])
            
        peakamplitude=np.array(peakamplitude)
        peakcenter=np.array(peakcenter)*dt
        peakwidth=np.array(peakwidth)*dt
            
        v,_=np.polyfit(peakcenter,xpoints,1)
        #~ plotc.plt.plot(peakcenter,xpoints,'r-')
        #~ plotc.plt.show()
        print "Source",sourceno+1,"velocity of p mode",abs(v),"Mm/min"

        wavespeed[sourceno,1]=abs(v)
    
    if p2mode:
        ####################################################################
        #~ p2-mode
        
        p2mode_filter_path=os.path.join(codedir,'p2mode_filter.fits')
        if os.path.exists(p2mode_filter_path):
            p2mode_filter=np.squeeze(pf.getdata(p2mode_filter_path))
        else:
            p2mode_filter=np.asfortranarray(np.zeros((Nx,1,Nt),dtype=float))
            modefilters.p2mode_filter(p2mode_filter)
            p2mode_filter=np.squeeze(p2mode_filter).T

        vzcc_p=np.real(np.fft.ifft2(np.fft.fft2(vzcc)*p2mode_filter))
        vzcc_p/=abs(vzcc_p).max()
        #~ plotc.colorplot(vzcc_p,vmax=0.5,centerzero=True)
        

        #~ datacursor(display='multiple',draggable='True')
        #plotc.plt.show()
        #quit()


        Npoints=4
        xpix=np.array([-30-ind*6 for ind in xrange(Npoints)])+sourcepix
        xpoints=(xpix-Nx//2+1)/Nx*Lx

        #~ plotc.draw_vlines(xpix)

        t_cutoff=np.zeros((Npoints,2))

        #~ Line fit to bottom of wavepacket
        def bot_time(x):
            return (209-57.9)/(172-226)*(x-sourcepix)+6

        #~ Line fit to top of wavepacket
        def top_time(x):
            return (89.1-275)/(229-185)*(x-sourcepix)+31

        t_cutoff[:,0]=bot_time(xpix)
        t_cutoff[:,1]=top_time(xpix)
        
        #~ plotc.plt.plot(xpix,t_cutoff[:,0],'b-')
        #~ plotc.plt.plot(xpix,t_cutoff[:,1],'g-')
        #~ plotc.plt.show()
        #~ quit()

        xdat=vzcc_p[:,xpix].T

        def local_max_index(arr):
            return (np.diff(np.sign(np.diff(arr))) < 0).nonzero()[0] + 1

        peakcenter=[]
        peakamplitude=[]
        peakwidth=[]

        def gaus(x,a,x0,sigma):
            return a*np.exp(-(x-x0)**2/(2*sigma**2))


        for xind,tseries in enumerate(xdat):
            tp=local_max_index(tseries)
            tp=tp[np.where(np.logical_and(tp<t_cutoff[xind,1],tp>t_cutoff[xind,0]))[0]]
            yp=tseries[tp]    

            #~ plotc.plot([xpix[xind]]*len(tp),tp,'o',color='orange')
            #~ plotc.show()
            
            popt,_ = curve_fit(gaus,tp,yp,p0=[yp[len(yp)//2],tp[len(tp)//2],(tp[-1]-tp[0])*0.4])
            
            peakamplitude.append(popt[0])
            peakcenter.append(popt[1])
            peakwidth.append(popt[2])
            
        peakamplitude=np.array(peakamplitude)
        peakcenter=np.array(peakcenter)*dt
        peakwidth=np.array(peakwidth)*dt
        
        ##~ plotc.plt.figure()
        #for xind,tseries in enumerate(xdat):
            #tp=local_max_index(tseries)
            #tp=tp[np.where(np.logical_and(tp<t_cutoff[xind,1],tp>t_cutoff[xind,0]))[0]]
            #yp=tseries[tp]    
            ##~ plotc.plot(tp,yp,'o')
            #gfit=gaus(tp,peakamplitude[xind],peakcenter[xind]/dt,peakwidth[xind]/dt)
            ##~ plotc.plot(tp,gfit)
        
        ##~ plotc.show()
            
        v,_=np.polyfit(peakcenter,xpoints,1)
        print "Source",sourceno+1,"velocity of p2 mode",abs(v),"Mm/min"

        wavespeed[sourceno,1]=abs(v)

#~ np.savetxt(os.path.join(datadir,'wavespeed'),wavespeed,fmt="%10.8f")
