from __future__ import division
import numpy as np
import pyfits
import os,re,sys,glob,fnmatch
import matplotlib
from matplotlib import ticker,pyplot as plt,rc
import read_params
import itertools
import dbyd2
from scipy import fftpack,interpolate,integrate
from scipy.special import j1

flows = read_params.if_flows()

datadir = read_params.get_directory()
updatedir = os.path.join(datadir,"update")
num_misfit_files=0

#~ Get a list of misfit files in the update directory
misfitfiles=sorted([os.path.join(updatedir,f) for f in fnmatch.filter(os.listdir(updatedir),'misfit_[0-9][0-9]')])
misfit_all_files=sorted([os.path.join(updatedir,f) for f in fnmatch.filter(os.listdir(updatedir),'misfit_all_[0-9][0-9]')])

num_misfit_files=len(misfitfiles)
num_misfit_all_files=len(misfit_all_files)

if num_misfit_files==0:
    print "No misfit files found"
    quit()

#~ If iteration cutoff is specified use it
itercutoff = read_params.parse_cmd_line_params('iter',mapto=int,default=np.inf)

num_misfit_files = min(itercutoff,num_misfit_files)
num_misfit_all_files = min(itercutoff,num_misfit_all_files)

#~ What to plot - data or model misfit?
mistype = read_params.parse_cmd_line_params("type",default="data")

#~ Get source location, useful if plotting sourcewise misfit (data_sourcewise)
srclocs = np.loadtxt(os.path.join(datadir,'master.pixels'),ndmin=1)
num_source = len(srclocs)

#~ Get ridges used
ridges=read_params.get_modes_used()
modes={'0':'f'}
for i in xrange(1,8): modes[str(i)]='p'+str(i)
modes['8']='first_bounce_p'

np.set_printoptions(linewidth=200,precision=4)

#~ Get filename to save to
save = read_params.parse_cmd_line_params("save")

if mistype == "data":
    markers = iter(('o', 'v', '8','s','<', 7, '*', 'h', '^', 'D', 'd'))
    linestyles = itertools.cycle(('solid','dashed','dotted'))
    modemisfit = np.zeros((len(ridges),num_misfit_files))

    rc('text', usetex=True)
    rc('font',**{'family':'serif','serif':['Times']})

    plt.subplot(121)

    for ridgeno,ridge in enumerate(ridges[:4]):

        nsources_found = 0

        for src in xrange(1,num_source+1):

            for iterno in xrange(num_misfit_files):

                ttfile = os.path.join(datadir,'tt','iter'+str(iterno).zfill(2),
                            'ttdiff_src'+str(src).zfill(2)+'.'+modes[ridge]+'mode.0')
                try:
                    tt=np.loadtxt(ttfile)
                    npix = tt.shape[0]
                    modemisfit[ridgeno,iterno] += sum((tt[:,1]/60)**2)/npix
                    nsources_found+=1
                except IOError:
                    pass

        if nsources_found>0:
            modemisfit[ridgeno] /= nsources_found
        else:
            print ridge,nsources_found

        try:
            plt.semilogy(range(num_misfit_files),modemisfit[ridgeno],marker=next(markers),color='black',
            ls=next(linestyles),label=modes[ridge])
        except: pass

    plt.tick_params(axis='both', which='major', labelsize=14)
    plt.ylabel('Mean misfit',fontsize=20)
    plt.legend(ncol=2)


    total_misfit = np.zeros(num_misfit_files)
    for fileno,misfitfile in enumerate(misfitfiles):
        total_misfit[fileno] = np.sum(np.loadtxt(misfitfile,usecols=[2]))

    plt.subplot(122)

    plt.semilogy(range(num_misfit_files),total_misfit,color='black',marker='o',ls='solid',zorder=1)
    plt.tick_params(axis='both', which='major', labelsize=14)
    plt.ylabel('Total misfit',fontsize=20)

    try:
        with open(os.path.join(datadir,'iter_modes')) as f:
            all_lines = f.readlines()
            for l in all_lines:
                plt.axvline(float(l.split()[0]),color='grey',ls='dashed',zorder=0)
    except IOError:
        pass

    for ax in plt.gcf().axes:
        ax.set_xlim(-0.5,num_misfit_files+2)
        ax.set_ylim(ax.dataLim.extents[1]/2,ax.dataLim.extents[3]*1.5)
        ax.set_xlabel("Iteration",fontsize=20)

    plt.gcf().set_size_inches(8,4)

    plt.tight_layout()

elif mistype == "model":

    z,c,rho = np.loadtxt(read_params.get_solarmodel(),usecols=[0,1,2],unpack=True)
    z = (z-1)*695.8
    c /= 1e2

    # true model
    kDH13 = 2*np.pi/30
    RDH13 = 15
    sigmazDH13 = 0.912
    z0DH13 = -2.3
    v0DH13 = 240

    def ddz(z_profile):
        darr = dbyd2.dbyd2(np.asfortranarray(np.atleast_2d(z_profile)),1)
        dz2d = dbyd2.dbyd2(np.asfortranarray(np.atleast_2d(z)),1)
        return np.squeeze(darr/dz2d)

    with np.load(os.path.join(datadir,'true_psi_coeffs.npz')) as f:
        tz = f['tz']
        cz_top = f['cz_top']
        cz_true = f['cz_top']+f['cz_bot']
        kz = f['kz']

    psi_true = interpolate.splev(z,(tz,cz_true,kz),ext=1)

    psi_norm = integrate.simps(psi_true**2,x=z)
    vz_norm = integrate.simps(c**2*psi_true**2,x=z)
    vx_norm = integrate.simps(1/rho**2*ddz(rho*c*psi_true)**2,x=z)

    misfit_psi = []
    misfit_vx = []
    misfit_vz = []

    for iterno in xrange(num_misfit_files):
        with np.load(os.path.join(datadir,'update',
            'model_psi_{:02d}_coeffs.npz'.format(iterno))) as f:
            cz_model = f['z'] + cz_top
            psi_model = interpolate.splev(z,(tz,cz_model,kz),ext=1)

        misfit_psi_i = integrate.simps((psi_true-psi_model)**2,x=z)/psi_norm
        misfit_psi.append(misfit_psi_i)

        misfit_vz_i = integrate.simps(c**2*(psi_true-psi_model)**2,x=z)/vz_norm
        misfit_vz.append(misfit_vz_i)

        misfit_vx_i = integrate.simps(1/rho**2*ddz(rho*c*(psi_true-psi_model))**2,x=z)/vx_norm
        misfit_vx.append(misfit_vx_i)

    plt.semilogy(misfit_psi,linestyle='solid',marker='o',label="$\psi$",color='black',zorder=1)
    plt.semilogy(misfit_vx,linestyle='dashed',marker='^',label="$v_x$",color='black',zorder=1)
    plt.semilogy(misfit_vz,linestyle='dotted',marker='s',label="$v_z$",color='black',zorder=1)

    plt.gca().yaxis.set_major_locator(ticker.LogLocator(subs=[1,3,5,7]))
    def percent(y,position):
        s = str(100 * y)

        # The percent symbol needs escaping in latex
        if matplotlib.rcParams['text.usetex'] is True:
            return s + r'$\%$'
        else:
            return s + '%'

    plt.gca().yaxis.set_major_formatter(ticker.FuncFormatter(percent))

    plt.grid(color='dimgrey')
    plt.legend(loc="best")


    try:
        with open(os.path.join(datadir,'iter_modes')) as f:
            mode_iter = f.readlines()
        for line in mode_iter:
            words=line.split()
            plt.axvline(float(words[0]),ls='dashed',color='grey',lw=1.5,zorder=0)
            modestr = ' '.join([modes[i] for i in words[1:]])
            import textwrap
            modestr = textwrap.fill(modestr,5)
            modestr = modestr.replace(' ',',').replace('\n',',\n')
            plt.text(int(words[0])+1,1,modestr,fontdict={'size':16},
            verticalalignment='top',bbox={'facecolor':'white'})

    except IOError:
        pass

    plt.xlim(-1,num_misfit_files+0.5)
    plt.ylim(plt.gca().dataLim.extents[1]/2,1.5)
    plt.xlabel("Iteration number")
    plt.ylabel("Model misfit")

    plt.tight_layout()


if save is not None:
    savepath = os.path.join("plots",save)
    if not os.path.exists("plots"): os.makedirs("plots")
    plt.savefig(savepath)
    print "Saved to",savepath
else:
    print "Not saving plot to file"

plt.show()
