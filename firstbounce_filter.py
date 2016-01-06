from __future__ import division
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import read_params
import os
import pyfits

def fitsread(f): return np.squeeze(pyfits.getdata(f))

poly=[]
plotted_pts=[]
plotted_lines=[]

def onclick(event):
    global poly
    if event.button == 1:
        eventx,eventy = event.xdata, event.ydata
        pt=plt.plot([eventx],[eventy],marker='o',color='b',linestyle='none')
        plotted_pts.append(pt)
        plt.draw()
        poly.append((eventx,eventy))
    
def onpress(event):
    
    global poly
    global plotted_pts
    global plotted_lines
    
    if event.key=="c":
        
        if not poly: return
        xcoords,ycoords = zip(*poly)
        if len(xcoords)<3: return
        pfit = np.polyfit(xcoords,ycoords,2)

        polystr=""
        for i,p_i in enumerate(pfit):
            polystr+=str(p_i)+"*x**"+str(len(pfit)-i-1)+"+"
        
        polystr=polystr.replace("*x**0","")
        polystr=polystr.replace("x**1","x")
        polystr=polystr.rstrip("+")
        polystr=polystr.rstrip("*")

        print polystr
        
        global x
        global srcloc
        fit_t = np.polyval(pfit,abs(x))
        l=plt.plot(x,fit_t,'g')
        plotted_lines.append(l)
        poly=[]
        
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
        
    elif event.key=="r":
        
        if plotted_pts:
            plotted_pts[-1].pop(0).remove()
            plotted_pts.pop()
            if poly: poly.pop()
        elif plotted_lines:
            plotted_lines[-1].pop(0).remove()
            plotted_lines.pop()
            
    
    plt.draw()
    
datadir=read_params.get_directory()

src = 4
srcloc = np.loadtxt(os.path.join(datadir,"master.pixels"))[src+1]

data = fitsread(os.path.join(datadir,'tt','data','data'+str(src).zfill(2)+'.fits'))
data/=data.max()

nx = read_params.get_nx()
Lx = read_params.get_xlength()
dx = Lx/nx
nt = data.shape[0]
dt_sec = read_params.get_dt()
dt_minutes = dt_sec/60

x = np.linspace(-Lx/2-dx/2,Lx/2+dx/2,num=nx+1,endpoint=False)
t= np.linspace(-dt_minutes/2,(nt+0.5)*dt_minutes,num=nt+1)

cmap_cutoff = 1e-1
plt.pcolormesh(x,t,data,cmap='RdBu_r',vmax=cmap_cutoff,vmin=-cmap_cutoff)

plt.xlim(x[0],x[-1])
plt.ylim(t[0],t[-1])

plt.draw()
leftclick = plt.connect('button_press_event', onclick)
shiftpress = plt.connect('key_press_event', onpress)

plt.show()



