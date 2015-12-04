from __future__ import division
import numpy as np
import pyfits as pyfits
import plotc
from mpldatacursor import datacursor
from scipy.optimize import curve_fit
import os,sys,re
import modefilters,read_params

codedir=os.path.dirname(os.path.abspath(__file__))
datadir=read_params.get_directory()

Lx=read_params.get_xlength()
dt_sec=read_params.get_dt()
dt_min=dt_sec/60

master_pixels=np.loadtxt(os.path.join(datadir,'master.pixels'),ndmin=1)
Nsources=len(master_pixels)

wavespeed=np.zeros((Nsources,6))

modes_to_plot=dict(
fmode=False,
p1mode=False,
p2mode=False,
p3mode=False,
p4mode=False,
p5mode=False)

ridge_filters_driver=read_params.get_ridge_filter()
paramsfiles=[os.path.splitext(f)[1][1:] for f in os.listdir(os.path.join(datadir)) if re.match(r'params.[0-9]$',f)]
ridge_filters=[ridge for ridge in ridge_filters_driver if ridge in paramsfiles]

modes={'0':'fmode'}
for pmodeno in xrange(1,6): modes.update({str(pmodeno):'p'+str(pmodeno)+'mode'})

for m in ridge_filters: modes_to_plot[modes[m]]=True

try:
    for m in sys.argv[1:]:
        if m.startswith("f"): modes_to_plot['fmode']=True
        elif m.startswith("p1"): modes_to_plot['p1mode']=True
        elif m.startswith("p2"): modes_to_plot['p2mode']=True
        elif m.startswith("p3"): modes_to_plot['p3mode']=True
        elif m.startswith("p4"): modes_to_plot['p4mode']=True
        elif m.startswith("p5"): modes_to_plot['p5mode']=True
except IndexError: pass

#~ Line fit to bottom of wavepacket, coordinates in pixels
def bot_time(x,(x1,y1),(x2,y2),source_bot):
    return (y2-y1)/(x2-x1)*(x-sourcepix)+source_bot
    
#~ Line fit to top of wavepacket, coordinates in pixels
def top_time(x,(x1,y1),(x2,y2),source_top):
    return (y2-y1)/(x2-x1)*(x-sourcepix)+source_top
    


for sourceno in xrange(Nsources):
    sourceloc=master_pixels[sourceno]
    sourcedir='forward_src'+str(sourceno+1).zfill(2)+'_ls00'

    vzcc=np.squeeze(pyfits.getdata(os.path.join(datadir,sourcedir,'data.fits')))
    Nt,Nx=vzcc.shape
    sourcepix=int(np.floor(sourceloc/Lx*Nx)+Nx/2-1)
    
    print "Source",sourceno,"located at",sourceloc,"Mm"
    
    if modes_to_plot['fmode']:
        ####################################################################
        #~ f-mode

        fmode_filter_path=os.path.join(codedir,'fmode_filter.fits')
        if os.path.exists(fmode_filter_path):
            fmode_filter=np.squeeze(pyfits.getdata(fmode_filter_path))
        else:
            fmode_filter=modefilters.fmode_filter(Nt,dt_sec,Nx,Lx)
            fmode_filter=np.squeeze(fmode_filter).T

        vzcc_f=np.fft.ifft2(np.fft.fft2(vzcc)*fmode_filter).real
        #~ plotc.colorplot(vzcc_f)
        #~ datacursor(display='multiple',draggable='True')


        Npoints=4
        xpix=np.array([-30-ind*6 for ind in xrange(Npoints)])+sourcepix
        xpoints=(xpix-Nx//2+1)/Nx*Lx

        #~ plotc.draw_vlines(xpix)

        t_cutoff=np.zeros((Npoints,3))

        #~ Line fit to bottom of wavepacket
        x1=64;y1=209;
        x2=108;y2=94.3;
        source_bot=10
        t_cutoff[:,0]=bot_time(xpix,(x1,y1),(x2,y2),source_bot)
        
        #~ Line fit to top of wavepacket
        x1=48.8;y1=266;
        x2=105;y2=66.6;
        source_top=66
        t_cutoff[:,1]=top_time(xpix,(x1,y1),(x2,y2),source_top)
        

        #~ plotc.plt.plot(xpix,t_cutoff[:,0],'b-')
        #~ plotc.plt.plot(xpix,t_cutoff[:,1],'g-')

        
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
            
            #~ plotc.plt.plot([xpix[xind]]*len(tp),tp,'o',color='orange')

            popt,_ = curve_fit(gaus,tp,yp,p0=[yp[len(yp)//2],tp[len(tp)//2],(tp[-1]-tp[0])*0.4])
            
            peakamplitude.append(popt[0])
            peakcenter.append(popt[1])
            peakwidth.append(popt[2])

        peakamplitude=np.array(peakamplitude)
        peakcenter=np.array(peakcenter)*dt_min
        peakwidth=np.array(peakwidth)*dt_min
        v,_=np.polyfit(peakcenter,xpoints,1)
        
        #~ print peakcenter,xpoints

        print "Source",sourceno+1,"velocity of f mode",abs(v),"Mm/min"
        
        #plotc.plt.figure()
        for xind,tseries in enumerate(xdat):
            tp=local_max_index(tseries)
            tp=tp[np.where(np.logical_and(tp<t_cutoff[xind,1],tp>t_cutoff[xind,0]))[0]]
            yp=tseries[tp]    
            #~ plotc.plt.plot(tp,yp,'o')
            gfit=gaus(tp,peakamplitude[xind],peakcenter[xind]/dt_min,peakwidth[xind]/dt_min)
            #print gfit
            #~ plotc.plt.plot(tp,gfit)
            
        #~ plotc.plt.show()
        
        
        wavespeed[sourceno,0]=abs(v)
    
    if modes_to_plot['p1mode']:
        ####################################################################
        #~ p1-mode
        
        pmode_filter_path=os.path.join(codedir,'p1mode_filter.fits')
        if os.path.exists(pmode_filter_path):
            pmode_filter=np.squeeze(pyfits.getdata(pmode_filter_path))
        else:
            pmode_filter=modefilters.pmode_filter(Nt,dt_sec,Nx,Lx)
            pmode_filter=np.squeeze(pmode_filter).T

        vzcc_p=np.real(np.fft.ifft2(np.fft.fft2(vzcc)*pmode_filter))
        #~ plotc.colorplot(vzcc_p,centerzero=True)

        #~ datacursor(display='multiple',draggable='True')

        Npoints=4
        xpix=np.array([-30-ind*6 for ind in xrange(Npoints)])+sourcepix
        xpoints=(xpix-Nx//2+1)/Nx*Lx

        #~ plotc.draw_vlines(xpix)

        t_cutoff=np.zeros((Npoints,2))

        #~ Line fit to bottom of wavepacket
        x1=76.9;y1=52.7;
        x2=20.8;y2=153;
        source_bot=6
        t_cutoff[:,0]=bot_time(xpix,(x1,y1),(x2,y2),source_bot)
        
        #~ Line fit to top of wavepacket
        x1=31.8;y1=254;
        x2=86;y2=99.4;
        source_top=55
        t_cutoff[:,1]=top_time(xpix,(x1,y1),(x2,y2),source_top)
        
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
        peakcenter=np.array(peakcenter)*dt_min
        peakwidth=np.array(peakwidth)*dt_min
            
        v,_=np.polyfit(peakcenter,xpoints,1)
        #~ plotc.plt.plot(peakcenter,xpoints,'r-')
        #~ plotc.plt.show()
        print "Source",sourceno+1,"velocity of p mode",abs(v),"Mm/min"

        wavespeed[sourceno,1]=abs(v)
    
    if modes_to_plot['p2mode']:
        ####################################################################
        #~ p2-mode
        
        p2mode_filter_path=os.path.join(codedir,'p2mode_filter.fits')
        if os.path.exists(p2mode_filter_path):
            p2mode_filter=np.squeeze(pyfits.getdata(p2mode_filter_path))
        else:
            p2mode_filter=modefilters.p2mode_filter(Nt,dt_sec,Nx,Lx)
            p2mode_filter=np.squeeze(p2mode_filter).T

        vzcc_p=np.real(np.fft.ifft2(np.fft.fft2(vzcc)*p2mode_filter))
        #~ plotc.colorplot(vzcc_p,centerzero=True)
        

        #~ datacursor(display='multiple',draggable='True')
        #plotc.plt.show()
        #quit()


        Npoints=4
        xpix=np.array([-30-ind*6 for ind in xrange(Npoints)])+sourcepix
        xpoints=(xpix-Nx//2+1)/Nx*Lx

        #~ plotc.draw_vlines(xpix)

        t_cutoff=np.zeros((Npoints,2))

        #~ Line fit to bottom of wavepacket
        x1=70.5;y1=49.1;
        x2=8.5;y2=140;
        source_bot=6
        t_cutoff[:,0]=bot_time(xpix,(x1,y1),(x2,y2),source_bot)
        
        #~ Line fit to top of wavepacket
        x1=39.5;y1=194;
        x2=86;y2=76.8;
        source_top=35
        t_cutoff[:,1]=top_time(xpix,(x1,y1),(x2,y2),source_top)
        
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
        peakcenter=np.array(peakcenter)*dt_min
        peakwidth=np.array(peakwidth)*dt_min
        
        ##~ plotc.plt.figure()
        #for xind,tseries in enumerate(xdat):
            #tp=local_max_index(tseries)
            #tp=tp[np.where(np.logical_and(tp<t_cutoff[xind,1],tp>t_cutoff[xind,0]))[0]]
            #yp=tseries[tp]    
            ##~ plotc.plot(tp,yp,'o')
            #gfit=gaus(tp,peakamplitude[xind],peakcenter[xind]/dt_min,peakwidth[xind]/dt_min)
            ##~ plotc.plot(tp,gfit)
        
        #~ plotc.show()
            
        v,_=np.polyfit(peakcenter,xpoints,1)
        print "Source",sourceno+1,"velocity of p2 mode",abs(v),"Mm/min"

        wavespeed[sourceno,2]=abs(v)

    if modes_to_plot['p3mode']:
        ####################################################################
        #~ p2-mode
        
        p3mode_filter_path=os.path.join(codedir,'p3mode_filter.fits')
        if os.path.exists(p3mode_filter_path):
            p3mode_filter=np.squeeze(pyfits.getdata(p3mode_filter_path))
        else:
            p3mode_filter=modefilters.p3mode_filter(Nt,dt_sec,Nx,Lx)
            p3mode_filter=np.squeeze(p3mode_filter).T

        vzcc_p=np.real(np.fft.ifft2(np.fft.fft2(vzcc)*p3mode_filter))
        #~ vzcc_p/=abs(vzcc_p).max()
        #~ plotc.colorplot(vzcc_p,centerzero=True)
#~ 
#~ 
        #~ datacursor(display='multiple',draggable='True')
        #~ plotc.plt.show()
        #~ quit()


        Npoints=4
        xpix=np.array([-30-ind*6 for ind in xrange(Npoints)])+sourcepix
        xpoints=(xpix-Nx//2+1)/Nx*Lx

        #~ plotc.draw_vlines(xpix)

        t_cutoff=np.zeros((Npoints,2))

        #~ Line fit to bottom of wavepacket
        x1=78.2;y1=32.3;
        x2=10.5;y2=117;
        source_bot=6
        t_cutoff[:,0]=bot_time(xpix,(x1,y1),(x2,y2),source_bot)
        
        #~ Line fit to top of wavepacket
        x1=39.5;y1=194;
        x2=86;y2=76.8;
        source_top=28
        t_cutoff[:,1]=top_time(xpix,(x1,y1),(x2,y2),source_top)

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
        peakcenter=np.array(peakcenter)*dt_min
        peakwidth=np.array(peakwidth)*dt_min
        
        ##~ plotc.plt.figure()
        #for xind,tseries in enumerate(xdat):
            #tp=local_max_index(tseries)
            #tp=tp[np.where(np.logical_and(tp<t_cutoff[xind,1],tp>t_cutoff[xind,0]))[0]]
            #yp=tseries[tp]    
            ##~ plotc.plot(tp,yp,'o')
            #gfit=gaus(tp,peakamplitude[xind],peakcenter[xind]/dt_min,peakwidth[xind]/dt_min)
            ##~ plotc.plot(tp,gfit)
        
        #~ plotc.show()
            
        v,_=np.polyfit(peakcenter,xpoints,1)
        print "Source",sourceno+1,"velocity of p3 mode",abs(v),"Mm/min"

        wavespeed[sourceno,3]=abs(v)

    if modes_to_plot['p4mode']:
        ####################################################################
        #~ p2-mode
        
        p4mode_filter_path=os.path.join(codedir,'p4mode_filter.fits')
        if os.path.exists(p4mode_filter_path):
            p4mode_filter=np.squeeze(pyfits.getdata(p4mode_filter_path))
        else:
            p4mode_filter=modefilters.p4mode_filter(Nt,dt_sec,Nx,Lx)
            p4mode_filter=np.squeeze(p4mode_filter).T

        vzcc_p=np.real(np.fft.ifft2(np.fft.fft2(vzcc)*p4mode_filter))
        vzcc_p/=abs(vzcc_p).max()
        #~ plotc.colorplot(vzcc_p,vmax=0.5,centerzero=True)
#~ 
        #~ datacursor(display='multiple',draggable='True')
        #~ plotc.plt.show()
        #~ quit()


        Npoints=4
        xpix=np.array([-30-ind*6 for ind in xrange(Npoints)])+sourcepix
        xpoints=(xpix-Nx//2+1)/Nx*Lx

        #~ plotc.draw_vlines(xpix)

        t_cutoff=np.zeros((Npoints,2))

        #~ Line fit to bottom of wavepacket
        x1=69.8;y1=36;
        x2=11.8;y2=95.8;
        source_bot=4
        t_cutoff[:,0]=bot_time(xpix,(x1,y1),(x2,y2),source_bot)
        
        #~ Line fit to top of wavepacket
        x1=102;y1=56.4;
        x2=46.6;y2=140;
        source_top=53
        t_cutoff[:,1]=top_time(xpix,(x1,y1),(x2,y2),source_top)

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
        peakcenter=np.array(peakcenter)*dt_min
        peakwidth=np.array(peakwidth)*dt_min
        
        ##~ plotc.plt.figure()
        #for xind,tseries in enumerate(xdat):
            #tp=local_max_index(tseries)
            #tp=tp[np.where(np.logical_and(tp<t_cutoff[xind,1],tp>t_cutoff[xind,0]))[0]]
            #yp=tseries[tp]    
            ##~ plotc.plot(tp,yp,'o')
            #gfit=gaus(tp,peakamplitude[xind],peakcenter[xind]/dt_min,peakwidth[xind]/dt_min)
            ##~ plotc.plot(tp,gfit)
        
        #~ plotc.show()
            
        v,_=np.polyfit(peakcenter,xpoints,1)
        print "Source",sourceno+1,"velocity of p4 mode",abs(v),"Mm/min"

        wavespeed[sourceno,4]=abs(v)

    if modes_to_plot['p5mode']:
        ####################################################################
        #~ p2-mode
        
        p5mode_filter_path=os.path.join(codedir,'p5mode_filter.fits')
        if os.path.exists(p5mode_filter_path):
            p5mode_filter=np.squeeze(pyfits.getdata(p5mode_filter_path))
        else:
            p5mode_filter=modefilters.p5mode_filter(Nt,dt_sec,Nx,Lx)
            p5mode_filter=np.squeeze(p5mode_filter).T

        vzcc_p=np.real(np.fft.ifft2(np.fft.fft2(vzcc)*p5mode_filter))
        #~ vzcc_p/=abs(vzcc_p).max()
        #~ plotc.colorplot(vzcc_p,centerzero=True)
        

        #~ datacursor(display='multiple',draggable='True')
        #~ plotc.plt.show()
        #~ quit()


        Npoints=4
        xpix=np.array([-30-ind*6 for ind in xrange(Npoints)])+sourcepix
        xpoints=(xpix-Nx//2+1)/Nx*Lx

        #~ plotc.draw_vlines(xpix)

        t_cutoff=np.zeros((Npoints,2))

        #~ Line fit to bottom of wavepacket
        x1=70.6;y1=28.8;
        x2=7.9;y2=86.8;
        source_bot=0
        t_cutoff[:,0]=bot_time(xpix,(x1,y1),(x2,y2),source_bot)
        
        #~ Line fit to top of wavepacket
        x1=86;y1=54.9;
        x2=33;y2=123;
        source_top=35
        t_cutoff[:,1]=top_time(xpix,(x1,y1),(x2,y2),source_top)

        #~ plotc.plot(xpix,t_cutoff[:,0],'b-')
        #~ plotc.plot(xpix,t_cutoff[:,1],'g-')
        #~ plotc.show()
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
        peakcenter=np.array(peakcenter)*dt_min
        peakwidth=np.array(peakwidth)*dt_min
        
        ##~ plotc.plt.figure()
        #for xind,tseries in enumerate(xdat):
            #tp=local_max_index(tseries)
            #tp=tp[np.where(np.logical_and(tp<t_cutoff[xind,1],tp>t_cutoff[xind,0]))[0]]
            #yp=tseries[tp]    
            ##~ plotc.plot(tp,yp,'o')
            #gfit=gaus(tp,peakamplitude[xind],peakcenter[xind]/dt_min,peakwidth[xind]/dt_min)
            ##~ plotc.plot(tp,gfit)
        
        #~ plotc.show()
            
        v,_=np.polyfit(peakcenter,xpoints,1)
        print "Source",sourceno+1,"velocity of p5 mode",abs(v),"Mm/min"

        wavespeed[sourceno,5]=abs(v)


np.savetxt('wavespeed',wavespeed,fmt="%14.8f")
