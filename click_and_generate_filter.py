from __future__ import division
import read_params
import numpy as np
import matplotlib.pyplot as plt
import plotc
import os,sys
import pyfits
import warnings

poly=[]
plotted_pts=[]
plotted_lines=[]
shift_pressed = False
lines_fitted = 0

np.set_printoptions(precision=3)

def line_picker(line, mouseevent):
    """
    find the points within a certain distance from the mouseclick in
    data coords and attach some extra attributes, pickx and picky
    which are the data points that were picked
    """
    if mouseevent.xdata is None:
        return False, dict()
    xdata = line.get_xdata()
    ydata = line.get_ydata()
    maxd = 0.05
    d = np.sqrt((xdata - mouseevent.xdata)**2. + (ydata - mouseevent.ydata)**2.)

    ind = np.nonzero(np.less_equal(d, maxd))
    if len(ind):
        pickx = np.take(xdata, ind)[0]
        picky = np.take(ydata, ind)[0]
        props = dict(line=line,ind=ind[0], pickx=pickx, picky=picky)
        return True, props
    else:
        return False, dict()

def fit_curve():
    
    global poly
    global plotted_pts
    global plotted_lines
    
    if not len(poly): return
    xcoords,ycoords = zip(*poly)
    if len(xcoords)<3: 
        if len(plotted_lines)>lines_fitted:
            plotted_lines[-1].pop(0).remove()
            del plotted_lines[-1]
            
        return
    
    if len(plotted_lines)>lines_fitted:
        plotted_lines[-1].pop(0).remove()
        del plotted_lines[-1]

    pfit = np.polyfit(xcoords,ycoords,2)
    
    global k
    fit_nu = np.polyval(pfit,abs(k))
    l=plt.plot(k,fit_nu,'g')
    
    plotted_lines.append(l)
    return "{0}*k**2 + {1}*k + {2}".format(*pfit)

def onpick(event):
    if len(event.ind) and len(poly):
        #~ Remove corresponding point from poly array
        xdata,ydata=map(np.array,zip(*poly))
        match=np.where(np.isclose(xdata,event.pickx) & np.isclose(ydata,event.picky))[0]
        if len(match):
            match=match[0]
        else: return
        event.line.remove()
        del poly[match]
        del plotted_pts[match]
        fit_curve()
        plt.draw()

def onclick(event):
    global poly
    global shift_pressed
    if event.button == 1 and shift_pressed:
        eventx,eventy = event.xdata, event.ydata
        pt=plt.plot([eventx],[eventy],marker='o',color='b',linestyle='none',picker=line_picker)
        plotted_pts.append(pt)
        poly.append((eventx,eventy))
        fit_curve()
        plt.draw()
        
def onpress(event):
    
    global poly
    global plotted_pts
    global plotted_lines
    global lines_fitted
    
    if event.key == 'c':
        polystr = fit_curve()
        print polystr
        poly=[]
        plotted_lines[-1][0].set_color('black')
        lines_fitted+=1
        
    elif event.key=="d":
        
        for obj in plotted_pts:
            line=obj.pop(0)
            line.remove()
        
        for obj in plotted_lines:
            line=obj.pop(0)
            line.remove()
        
        plotted_pts=[]
        plotted_lines=[]
        poly=[]
        lines_fitted = 0
      
    elif event.key=="shift":
        global shift_pressed
        shift_pressed = True
    
    plt.draw()

def onrelease(event):
    global shift_pressed
    if event.key=='shift':
        shift_pressed = False

def fitswrite(filename,array):
    warnings.filterwarnings('ignore')
    pyfits.writeto(filename,array,clobber=True)

codedir=os.path.dirname(os.path.abspath(__file__))
datadir=read_params.get_directory()

try:
    src=next(f for f in sys.argv if (f.startswith("src=") or f.startswith("source=")))
    src=int(src.split("=")[-1])
except StopIteration: src=1

srcloc=np.loadtxt(os.path.join(datadir,'master.pixels'),ndmin=1)[src-1]

data=np.squeeze(pyfits.getdata(os.path.join(datadir,'forward_src'+str(src).zfill(2)+'_ls00','data.fits')))

nt,nx=data.shape
Lx=read_params.get_xlength()
dt=read_params.get_dt()
x=np.linspace(-Lx/2,Lx/2,nx,endpoint=False)
t=np.arange(nt)*dt

nu=np.fft.fftshift(np.fft.fftfreq(nt,dt))*1e3
dnu=nu[1]-nu[0]
nu_edges = np.linspace(nu[0]-dnu/2,nu[-1]+dnu/2,num=nt+1)
k=np.fft.fftshift(np.fft.fftfreq(nx,Lx/nx))*2*np.pi
dk = k[1]-k[0]
k_edges = np.linspace(k[0]-dk/2,k[-1]+dk/2,num=nx+1)

spectrum = np.fft.fftshift(abs(np.fft.fft2(data)))

#########################################################################################


plt.pcolormesh(k_edges,nu_edges,spectrum/spectrum.max(),cmap='Oranges',vmax=0.3)

plt.connect('button_press_event', onclick)
plt.connect('key_press_event', onpress)
plt.connect('key_release_event', onrelease)
plt.connect('pick_event', onpick)

plt.xlim(0,1.5)
plt.ylim(1.2,6)

plt.show()

