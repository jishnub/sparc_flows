
import read_params
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import plotc
import os,sys
import pyfits
import warnings

class Point():
    def __init__(self,Line2D):
        self.x  = Line2D.get_xdata()[0]
        self.y  = Line2D.get_ydata()[0]
        self.artist = Line2D
        
    def remove(self):
        self.artist.remove()
        self.artist=None

class Poly():
    def __init__(self):
        self.points=[]
        self.fitted_curve=None
        
    def fit(self,order=4):
    
        #~ If the number of points are not sufficient, remove fitted polynomial if it exists
        if len(self.points)<order+1: 
            if self.fitted_curve is not None: self.fitted_curve.remove()
            self.fitted_curve=None
            return (None,None)

        #~ If there are sifficient points, fit the points with a polynomial function        
        xcoords = [pt.x for pt in self.points]
        ycoords = [pt.y for pt in self.points]
        pfit = np.polyfit(xcoords,ycoords,order)
        
        #~ Generate points on a fine grid along the fitted polynomial
        #~ This is to plot a continuous line
        fit_nu = np.polyval(pfit,abs(k))
        
        #~ Update fitted curve
        if self.fitted_curve is not None: self.fitted_curve.remove()
        self.fitted_curve,=plt.plot(k,fit_nu,'g')
        
        #~ String to format fitted polynomial like p_2*k**2 + p_1*k +  p_0
        fmtstr=" + ".join(["{}*k**"+str(i) if i>1 else "{}*k" if i==1 else "{}"  for i in range(len(pfit))[::-1]])
        polystr=fmtstr.format(*pfit).replace("+ -","- ")
        
        fitstr=[]
        for index,p_i in enumerate(pfit[::-1]):
            fitstr.append("index("+str(index)+") = "+str(p_i))
        fitstr="\n".join(fitstr)
        
        return polystr,fitstr

    def remove_match(self,artist):
        ''' Find and remove the point that was clicked on '''
        for point in self.points:
            if point.artist == artist:
                point.remove()
                self.poly.points.remove(point)
                break
    
    def clear(self):
        ''' Refresh the working slate by deleting plots and lines.
        Points and lines already finalized are untouched. '''
        for point in self.points:
            point.remove()
        if self.fitted_curve is not None:
            self.fitted_curve.remove()
            self.fitted_curve=None
        self.points=[]

class Track_interactions():
    def __init__(self,figure):
        self.button_press_event=figure.canvas.mpl_connect('button_press_event', self.onclick)
        self.key_press_event=figure.canvas.mpl_connect('key_press_event', self.onpress)
        self.key_release_event=figure.canvas.mpl_connect('key_release_event', self.onrelease)
        self.pick_event=figure.canvas.mpl_connect('pick_event', self.onpick)
        self.shift_pressed=False
        self.poly=Poly()
        
    def onpick(self,event):
        ''' Remove a point when it is clicked on '''
        self.poly.remove_match(event.artist)
        self.poly.fit()
        plt.draw()
            
    def onclick(self,event):
        ''' Add a point at the (x,y) coordinates of click '''
        if event.button == 1 and self.shift_pressed:
            pt_artist,=plt.plot([event.xdata],[event.ydata],marker='o',color='b',linestyle='none',picker=5)
            self.poly.points.append(Point(pt_artist))
            self.poly.fit()
            plt.draw()
        
    def onpress(self,event):
    
        if event.key == 'c':
            polystr,fitstr = self.poly.fit()
            if polystr is None: return
            print(polystr,"\n",fitstr,"\n")
            self.poly.fitted_curve.set_color("black")
            plt.draw()
            self.poly=Poly()
        
        elif event.key=="d":
            self.poly.clear()
            
        elif event.key=="shift":
            self.shift_pressed = True
        
        plt.draw()

    def onrelease(self,event):
        if event.key=='shift':
            self.shift_pressed = False


datadir=read_params.get_directory()

src=read_params.parse_cmd_line_params("src",mapto=int,default=1)

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

figure=plt.figure()

plt.pcolormesh(k_edges,nu_edges,spectrum/spectrum.max(),cmap='Oranges',vmax=0.6)
plt.text(1.4,2.9,"f",fontsize=20)
plt.text(1.42,4.28,"p1",fontsize=20)
plt.text(1.4,5.1,"p2",fontsize=20)
plt.text(1.33,5.8,"p3",fontsize=20)
plt.text(1.04,5.78,"p4",fontsize=20)
plt.text(0.83,5.8,"p5",fontsize=20)
plt.text(0.7,5.8,"p6",fontsize=20)
plt.text(0.58,5.8,"p7",fontsize=20)

plt.xlim(0,1.5)
plt.ylim(1.2,6)

_=Track_interactions(figure)

plt.show()

