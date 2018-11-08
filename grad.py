import numpy as np
from scipy import interpolate,integrate
from scipy.special import j1,j0,jn
def j2(z): return jn(2,z)
def j1prime(z): return 0.5*(j0(z)-j2(z))
import os,fnmatch,sys
from astropy.io import fits
import read_params
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt,ticker
from pathlib import Path

datadir = Path(read_params.get_directory())
updatedir = datadir/"update"

########################################################################

def fitsread(f):
    
    with fits.open(f) as hdul:
        return hdul[0].data

def fitswrite(f,arr):
    fits.writeto(f,arr,overwrite=True)


def antisymmetrize(arr):   arr[:]=0.5*(arr[:]-arr[::-1])
def symmetrize(arr):
    return 0.5*(arr + np.flip(arr,axis=2))


def rms(arr): return np.sqrt(np.sum(arr**2)/np.prod(arr.shape))


########################################################################

def main():

    iterno=len(fnmatch.filter(os.listdir(datadir/"update"),'misfit_[0-9][0-9]'))-1

    args=sys.argv[1:]
    optimization_algo=read_params.parse_cmd_line_params(key="algo",default="bfgs")

    steepest_descent = optimization_algo.lower()=='sd'
    conjugate_gradient = optimization_algo.lower()=='cg'
    BFGS = optimization_algo.lower()=='bfgs'

    if not (BFGS or steepest_descent or conjugate_gradient):
        print("No matching optimization algorithm, quitting")
        quit()

    def isfloat(value):
        try:
            float(value)
            return True
        except ValueError:
            return False

    eps = list(map(float,list(filter(isfloat,args))))
    num_ls_per_src = len(fnmatch.filter(os.listdir(datadir),"forward_src01_ls[0-9][1-9]"))

    if len(eps) > num_ls_per_src:
        print("{} step sizes specified by {} linesearches to be done. Ignoring the last {} values".format(
                                                        len(eps),num_ls_per_src,len(eps)-num_ls_per_src))
        eps = eps[:num_ls_per_src]

    elif len(eps)==0:
        eps=[0.01*i for i in range(1,num_ls_per_src+1)]

    Rsun=6.95989467700E2 # Mm
    z = np.loadtxt(read_params.get_solarmodel(),usecols=[0]); z=(z-1)*Rsun

    #~ Get shape
    nx = read_params.get_nx()
    ny = 1
    nz = read_params.get_nz()

    Lx = read_params.get_xlength()
    dx = Lx/nx

    x = np.linspace(-Lx/2,Lx/2,nx,endpoint=False)

    def integrate_2D(arr):
        # Assume (nz,nx) format
        return integrate.simps(integrate.simps(arr,dx=dx,axis=1),x=z,axis=0)

    solarmodel=np.loadtxt(read_params.get_solarmodel())
    c = solarmodel[:,1]
    c_surface = c[abs(z).argmin()]
    
    c_quiet_3D = np.tile(c[:,None,None],(1,1,nx))

    num_src=len(np.loadtxt(datadir/'master.pixels',ndmin=1))

    ############################################################################

    def read_model(var='psi',iterno=iterno):
        modelfile = updatedir/'model_{}_{:02d}.fits'.format(var,iterno)
        return fitsread(modelfile)
        # modelfile = updatedir/('model_'+var+'_'+str(iterno).zfill(2)+'_coeffs.npz')
        # with np.load(modelfile) as f:
        #     return dict(list(f.items()))

    def read_grad(var='psi',iterno=iterno):
        modelfile = updatedir/'gradient_{}_{:02d}.fits'.format(var,iterno)
        return fitsread(modelfile)

        # modelfile = updatedir/('gradient_'+var+'_'+str(iterno).zfill(2)+'.npz')
        # with np.load(modelfile) as f:
        #     return dict(list(f.items()))

    def read_update(var='psi',iterno=iterno):
        modelfile = updatedir/'update_{}_{:02d}.fits'.format(var,iterno)
        return fitsread(modelfile)

        # modelfile = updatedir/('update_'+var+'_'+str(iterno).zfill(2)+'.npz')
        # with np.load(modelfile) as f:
        #     return dict(list(f.items()))

    # def read_BFGS_hessian(var='psi',iterno=iterno):
        
    #     modelfile = updatedir/('BFGS_hessian_'+var+'_'+str(iterno).zfill(2)+'.npz')
    #     with np.load(modelfile) as f:
    #         return dict(list(f.items()))

    def read_kern(var='psi',src=1):
        return fitsread(datadir/'kernel'/'kernel_{}_{:02d}.fits'.format(var,src))

    current_model = read_model(var="c",iterno=iterno)
    
    true_dc = fitsread("dc_true.fits")

    fig = plt.figure()

    ############################################################################
    # Gradient computation
    ############################################################################

    for var in ["c"]:

        kernel=np.zeros_like(current_model)
        hess=np.zeros_like(kernel)

        for src in range(1,num_src+1):

            kernel_i = read_kern(var=var,src=src)

            # fig.clf()
            # plt.pcolormesh(x,z,kernel_i.squeeze(),vmax=abs(kernel_i).max(),vmin=-abs(kernel_i).max(),cmap="RdBu_r")
            # plt.savefig(str(updatedir/"kernel_{}_src{:d}.png".format(var,src)))
            kernel += kernel_i

            hessian_i = abs(fitsread(datadir/'kernel'/'hessian_{:02d}.fits'.format(src)))
            hess += hessian_i

        hess=hess * np.reshape(solarmodel[:,2],(nz,1,1))
        hess = hess/abs(hess).max()
        hess[hess<5e-3]=5e-3

        kernel /= hess
        kernel = symmetrize(kernel)

        fitswrite(updatedir/"gradient_{}_{:02d}.fits".format(var,iterno),kernel)

        fig.clf()
        plt.pcolormesh(x,z,kernel.squeeze(),cmap="RdBu_r",vmax=abs(kernel).max(),vmin=-abs(kernel).max())
        plt.savefig(str(updatedir/"gradient_{}_{:02d}.png".format(var,iterno)))


    # fitswrite(updatedir/"grad_psi_{:02d}.fits".format(iterno),kernel)

    # psi_true = fitsread("true_psi.fits").squeeze().T # shape would be (nz,nx)

    # Plot gradient (kernel)
    # f=plt.figure()
    # ax1=plt.subplot(121)
    # plt.pcolormesh(x,z,psi_true,cmap="RdBu_r")
    # plt.title("True psi",fontsize=16)
    # plt.xlim(-100,100)
    # plt.ylim(-20,z[-1])
    # plt.xlabel("x (Mm)",fontsize=16)
    # plt.ylabel("z (Mm)",fontsize=16)

    # ax2=plt.subplot(122)
    # plt.pcolormesh(x,z,kernel/abs(kernel).max(),cmap="RdBu_r",vmax=1,vmin=-1)
    # plt.title("Gradient",fontsize=16)
    # ax2.set_xlim(ax1.get_xlim())
    # ax2.set_ylim(ax1.get_ylim())
    # plt.xlabel("x (Mm)",fontsize=16)
    # plt.ylabel("z (Mm)",fontsize=16)

    # f.set_size_inches(8,3.5)
    # plt.tight_layout()
    # plt.savefig(str(updatedir/"grad_{:02d}.png".format(iterno)))

    # compute basis coefficient gradients
    # def compute_grad_basis(coeffs):

    #     hs = interpolate.splev(z,(tz,coeffs["z"].iterated+cz_ref_top,kz),ext=1)
    #     g_c = f0_x[None,:]

    #     grad = dict.fromkeys(list(coeffs.keys()))
    #     for key in grad:
    #         grad[key] = np.zeros(coeffs[key].size)

    #     for j in coeffs["z"].get_range():
    #         bj = np.zeros_like(coeffs["z"].iterated)
    #         bj[j] = 1
    #         bz = interpolate.splev(z,(tz,bj,kz),ext=1)
    #         grad["z"][j] = integrate_2D(kernel*bz[:,None]*g_c)

    #     return grad

    #~ Write out gradients for this iteration
    
    # np.savez(updatedir/'gradient_psi_{:02d}.npz'.format(iterno),**compute_grad_basis(coeffs))

    ############################################################################
    # Optimization
    ############################################################################

    def normalize(d):
        if d.max()!=0:
            d /= d.max()
        # for key,value in list(d.items()):
        #     if value.max()!=0:
        #         d[key] = value/value.max()

    # update={}

    #~ Get update direction based on algorithm of choice
    if (iterno==0) or steepest_descent:

        def sd_update(var='psi'):
            grad = read_grad(var=var,iterno=iterno)
            update = - grad

            print("Steepest descent")
            return update

        update = sd_update(var='c')

        if iterno > 0: print('Forcing steepest descent')

    elif conjugate_gradient:

        def get_beta(grad,lastgrad,lastupdate):
            #~ Choose algo
            polak_ribiere = True
            hestenes_stiefel = True and (not polak_ribiere)

            np.seterr(divide="raise")
            try:
                beta = np.dot(grad.flatten(),(grad - lastgrad).flatten())
                if polak_ribiere:
                    beta/= np.sum(lastgrad**2.)
                elif hestenes_stiefel:
                    beta/= np.dot((grad - lastgrad).flatten(),lastupdate.flatten())
            except FloatingPointError:
                beta=0

            return beta

        def cg_update(var='psi'):
            grad_k = read_grad(var=var,iterno=iterno)
            grad_km1 = read_grad(var=var,iterno=iterno-1)
            p_km1 = read_update(var=var,iterno=iterno-1)

            beta_k = get_beta(grad_k,grad_km1,p_km1)
            update = -grad_k +  beta_k*p_km1

            np.set_printoptions(linewidth=200,precision=4)

            # beta_k = {}
            # update = {}
            # for param in list(grad_k.keys()):
            #     beta_k[param] = get_beta(grad_k[param],grad_km1[param],p_km1[param])
            #     update[param] = -grad_k[param] +  beta_k[param]*p_km1[param]
            
            print("beta",beta_k)

            return update

        update = cg_update(var='c')

    # elif BFGS:
    #     def BFGS_update(var='psi'):
    #         print("Using BFGS")
    #         grad_k = read_grad(var=var,iterno=iterno)
    #         grad_km1 = read_grad(var=var,iterno=iterno-1)

    #         model_k  = read_model(var=var,iterno=iterno)
    #         model_km1  = read_model(var=var,iterno=iterno-1)

    #         try:
    #             Hkm1 = read_BFGS_hessian(var=var,iterno=iterno)
    #         except IOError:
    #             Hkm1 = {}
    #             for key in list(grad_k.keys()):
    #                 Hkm1[key] = np.identity(grad_k[key].size)


    #         Hk = {}
    #         for key in list(grad_k.keys()):
    #             y_km1 = grad_k[key] - grad_km1[key]
    #             s_km1 = model_k[key] - model_km1[key]

    #             if (not y_km1.any()) or (not s_km1.any()):
    #                 Hk[key] = Hkm1[key]
    #             else:

    #                 rho_km1 = 1/y_km1.dot(s_km1)

    #                 left = np.identity(s_km1.size) - rho_km1*np.outer(s_km1,y_km1)
    #                 right = np.identity(s_km1.size) - rho_km1*np.outer(y_km1,s_km1)

    #                 Hk[key] = left.dot(Hkm1[key]).dot(right) + rho_km1*np.outer(s_km1,s_km1)

    #             update[key] = -Hk[key].dot(grad_k[key])

    #         hessfile = updatedir/('BFGS_hessian_{}_{:02d}.npz'.format(var,iterno))
    #         np.savez(hessfile,**Hk)

    #         return update

    #     update = BFGS_update(var='psi')

    
    ############################################################################

    # Deep z cutoff
    # update["z"] *= 1/(1+np.exp(-(np.arange(cz_ref_bot.size)-c_deep_z_cutoff_index)/0.2))

    normalize(update)
    update *= c_surface
    fitswrite(updatedir/"update_c_{:02d}.fits".format(iterno),update)
    # np.savez(updatedir/'update_psi_{:02d}.npz'.format(iterno),**update)

    fig.clf()
    plt.subplot(121)
    plt.pcolormesh(x,z,true_dc.squeeze(),cmap="RdBu_r",vmax=abs(true_dc).max(),vmin=-abs(true_dc).max())
    plt.xlim(-50,50)
    plt.ylim(-10,z[-1])
    plt.colorbar()
    plt.title("true dc")
    plt.subplot(122)
    plt.pcolormesh(x,z,update.squeeze(),cmap="RdBu_r",vmax=abs(update).max(),vmin=-abs(update).max())
    plt.xlim(-50,50)
    plt.ylim(-10,z[-1])
    plt.colorbar()
    plt.title("Update")
    plt.tight_layout()
    plt.savefig(str(updatedir/"update_c_{:02d}.png".format(iterno)))

    # plt.figure()
    # plt.subplot(211)
    # plt.bar(range(len(cz_ref_top)),cz_ref_top,facecolor="lightgreen",edgecolor="green",label="above the surface")
    # plt.bar(range(len(cz_ref_bot)),cz_ref_bot,facecolor="skyblue",edgecolor="blue",label="below the surface")
    # plt.ylabel("Coefficient value")
    # plt.title("Spline coefficients of true model")
    # plt.legend(loc="best")
    # plt.subplot(212)
    # plt.bar(range(len(update["z"])),update["z"])
    # plt.title("Normalized update")
    # plt.xlabel("Coefficient index")
    # plt.ylabel("Coefficient value")
    # plt.tight_layout()
    # plt.savefig(str(updatedir/"update.png"))
    # plt.clf()

    ############################################################################

    

    def create_ls_model(var='psi',eps=0,kind='linear'):

        # ls_cz = coeffs.get("z").iterated.copy()

        if kind=='linear':
            
            lsmodel = current_model + eps*update

        # elif kind=='exp':
        #     ls_cz *= 1+eps*update


        # print("eps {} coeffs".format(eps),ls_cz[:coeff_surf_cutoff_ind])

        # h_z=interpolate.splev(z,(tz,ls_cz+cz_ref_top,kz),ext=1)
        # lsmodel = f0_x[None,:]*h_z[:,None] # shape would be (nz,nx)

        # lsmodel = lsmodel[:,np.newaxis,:]

        # coeffdict = {"z":ls_cz}

        return lsmodel



    #~ Create models for linesearch
    for var in ["c"]:
        for i,eps_i in enumerate(eps):

            lsmodel = create_ls_model(var=var,eps=eps_i,kind='linear')

            fig.clf()
            plt.subplot(121)
            plt.pcolormesh(x,z,true_dc.squeeze(),cmap="RdBu_r",vmax=abs(true_dc).max(),vmin=-abs(true_dc).max())
            plt.xlim(-50,50)
            plt.ylim(-10,z[-1])
            plt.colorbar()
            plt.subplot(122)
            test_dc = (lsmodel- c_quiet_3D ).squeeze()
            plt.pcolormesh(x,z,test_dc, cmap="RdBu_r",vmax = abs(test_dc).max(),vmin = -abs(test_dc).max())
            plt.xlim(-50,50)
            plt.ylim(-10,z[-1])
            plt.colorbar()
            

            plt.tight_layout()
            plt.savefig(str(updatedir/"test_{}_{:d}.png".format(var,i+1)))
            
            # test_model_vertical_slice_ax.plot(z,lsmodel[:,0,peak_value_x],'o',ms=3,label="model {:d}".format(i+1),zorder=3)

            fitswrite(updatedir/'test_{}_{:d}.fits'.format(var,i+1), lsmodel)
        # np.savez(updatedir/'test_psi_{:d}_coeffs.npz'.format(i+1),**model_coeffs)

    # plt.legend(loc="best")
    # test_model_vertical_slice_fig.tight_layout()
    # test_model_vertical_slice_fig.savefig(str(updatedir/'test_psi_vertical.png'.format(i+1)))

    #~ Update epslist
    epslist_path = datadir/"epslist.npz"
    try:
        epslist=np.load(epslist_path)
        iter_done_list = epslist.files
        if str(iterno) in iter_done_list:
            ls_iter = epslist[str(iterno)]
            ls_iter[2:4] = ls_iter[0:2]
            ls_iter[0] = eps
            ls_iter[1] = 0
            epslist = {i:epslist[i] for i in iter_done_list}
            epslist.update({str(iterno):ls_iter})
        else:
            arr = np.zeros((4,num_ls_per_src))
            arr[0] = eps
            epslist = {i:epslist[i] for i in iter_done_list}
            epslist[str(iterno)] = arr
    except IOError:
        arr = np.zeros((4,num_ls_per_src))
        arr[0] = eps
        epslist = {str(iterno):arr}

    np.savez(epslist_path,**epslist)


########################################################################

if __name__ == "__main__":
    main()
