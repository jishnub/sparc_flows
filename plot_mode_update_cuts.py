from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
import read_params
import pyfits
from matplotlib.ticker import MaxNLocator

def fitsread(f): return np.squeeze(pyfits.getdata(f))

z=np.loadtxt(read_params.get_solarmodel(),usecols=[0]);z=(z-1)*695.8

update_f = fitsread('/scratch/jishnu/flows/f_3hr/update/update_psi_00.fits')
update_p1 = fitsread('/scratch/jishnu/flows/p1_3hr/update/update_psi_00.fits')
update_p2 = fitsread('/scratch/jishnu/flows/p2_3hr/update/update_psi_00.fits')

true_psi = fitsread('true_psi.fits')

update_f_argmax=divmod(update_f.argmax(),update_f.shape[1])
update_p1_argmax=divmod(update_p1.argmax(),update_p1.shape[1])
update_p2_argmax=divmod(update_p2.argmax(),update_p2.shape[1])
true_psi_argmax=divmod(true_psi.argmax(),true_psi.shape[1])

nplots=4
gl=iter(int("1"+str(nplots)+str(i)) for i in xrange(1,nplots+1))

plt.subplot(next(gl))
plt.plot(update_f[:,update_f_argmax[1]],z)
plt.plot([0]*len(z),z,linestyle='dotted')
plt.ylim(-10,z.max())
plt.gca().xaxis.set_major_locator(MaxNLocator(4,prune='upper'))
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.title("f mode",fontsize=20)
plt.tick_params(axis='both', which='major', labelsize=14)
plt.ylabel("Depth (Mm)",fontsize=20)

plt.subplot(next(gl))
plt.plot(update_p1[:,update_p1_argmax[1]],z)
plt.plot([0]*len(z),z,linestyle='dotted')
plt.ylim(-10,z.max())
plt.gca().xaxis.set_major_locator(MaxNLocator(4,prune='upper'))
plt.setp(plt.gca().get_yticklabels(),visible=False)
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.title("p1 mode",fontsize=20)
plt.tick_params(axis='both', which='major', labelsize=14)

plt.subplot(next(gl))
plt.plot(update_p2[:,update_p2_argmax[1]],z)
plt.plot([0]*len(z),z,linestyle='dotted')
plt.ylim(-10,z.max())
plt.gca().xaxis.set_major_locator(MaxNLocator(4,prune='upper'))
plt.setp(plt.gca().get_yticklabels(),visible=False)
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.title("p2 mode",fontsize=20)
plt.tick_params(axis='both', which='major', labelsize=14)

plt.subplot(next(gl))
plt.plot(true_psi[:,true_psi_argmax[1]],z,color='red')
plt.ylim(-10,z.max())
plt.gca().xaxis.set_major_locator(MaxNLocator(4,prune='both'))
plt.setp(plt.gca().get_yticklabels(),visible=False)
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.title(r"True $\psi$",fontsize=20)
plt.tick_params(axis='both', which='major', labelsize=14)

plt.tight_layout()
plt.subplots_adjust(wspace=0)


plt.show()
