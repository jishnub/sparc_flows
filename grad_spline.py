from __future__ import division
import numpy as np
from scipy.interpolate import UnivariateSpline
from scipy.ndimage.filters import gaussian_filter1d
from scipy import interpolate,integrate
import os,fnmatch,sys
import pyfits
import warnings
import read_params
from matplotlib import pyplot as plt

#######################################################################

def fitsread(f):
    try:
        arr=pyfits.getdata(f)
    except IOError:
        raise IOError
    # If it is 2D, make it 3D.
    # Dimension will be (nz,1,nx) after this
    if len(arr.shape)==2: arr=arr[:,np.newaxis,:]
    # Change dimension to (nx,ny,nz)
    return arr.transpose(2,1,0)

def fitswrite(f,arr):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        # Change dimensions from (nx,ny,nz) to (nz,ny,nx)
        arr=arr.transpose(2,1,0)
        # clobber=True rewrites a pre-existing file
        pyfits.writeto(f,arr,clobber=True)

def get_iter_no():
    datadir = read_params.get_directory()
    updatedir=os.path.join(datadir,"update")
    # Count the number of misfit_xx files
    return len(fnmatch.filter(os.listdir(updatedir),'misfit_[0-9][0-9]'))-1

def get_number_of_sources():
    # Count the number of lines in master.pixels
    # Safest to use loadtxt, since it removes newline correctly
    datadir = read_params.get_directory()
    return len(np.loadtxt(os.path.join(datadir,'master.pixels'),ndmin=1))

def filterx(kern,nk):

    Lx = read_params.get_xlength()
    nx = read_params.get_nx()
    k = np.fft.rfftfreq(nx,1./nx)
    sigmak = nk # retain an arbitrary number of wavenumbers
    filt_k = np.exp(-k**2/(2*sigmak**2))

    filt_k = filt_k[:,np.newaxis,np.newaxis]

    return np.fft.irfft(np.fft.rfft(kern,axis=0)*filt_k,axis=0).real

def filterz(arr,algo='spline',sp=1.0):
    nx,ny,nz=arr.shape
    z_ind=np.arange(nz)
    b=arr.copy()

    if algo=='smooth':
        coeffs=np.zeros(6)
        coeffs[0]=0.75390625000000
        coeffs[1]=0.41015625000000
        coeffs[2]=-0.23437500000000
        coeffs[3]=0.08789062500000
        coeffs[4]=-0.01953125000000
        coeffs[5]=0.00195312500000
        sigmaz=4

        temp=np.zeros_like(arr)

        kst=1

        for k in xrange(5,nz-5):
            temp[:,:,k] =( coeffs[0]*arr[:,:,k] + 0.5*coeffs[1]*(arr[:,:,k-1] + arr[:,:,k+1])
                        +  0.5*coeffs[2]*(arr[:,:,k-2] + arr[:,:,k+2])
                        +  0.5*coeffs[3]*(arr[:,:,k-3] + arr[:,:,k+3])
                        +  0.5*coeffs[4]*(arr[:,:,k-4] + arr[:,:,k+4])
                        +  0.5*coeffs[5]*(arr[:,:,k-5] + arr[:,:,k+5])  )

        if (kst==1):

            temp[:,:,nz-5:nz-1] = 0.5 * (arr[:,:,nz-4:] + arr[:,:,nz-6:nz-2])
            temp[:,:,1:4] = 0.5 * (arr[:,:,:3] + arr[:,:,2:5])

        arr[:,:,kst:nz-kst+1]=temp[:,:,kst:nz-kst+1]

    elif algo=='spline':
        for x_ind in xrange(nx):
            for y_ind in xrange(ny):
                arrz=arr[x_ind,y_ind]
                arrzmax=arrz.max()
                arrz=arrz/arrzmax

                s=UnivariateSpline(z_ind,arrz,s=sp)
                arr[x_ind,y_ind]=s(z_ind)*arrzmax

    elif algo=='gaussian':
        arr[:]=gaussian_filter1d(arr,sigma=sp,axis=-1)

def antisymmetrize(arr):   arr[:]=0.5*(arr[:]-arr[::-1])
def symmetrize(arr):   arr[:]=0.5*(arr[:]+arr[::-1])

def updatedir(filename):
    datadir = read_params.get_directory()
    return os.path.join(datadir,'update',filename)

def rms(arr): return np.sqrt(np.sum(arr**2)/np.prod(arr.shape))

def filter_and_symmetrize(totkern,hess,sym=None,z_filt_algo='gaussian',z_filt_pix=0.3):
    kern = totkern/hess
    kern = filterx(kern,30)
    if sym=='sym':
        symmetrize(kern)
    elif sym=='asym':
        antisymmetrize(kern)
    filterz(kern,algo=z_filt_algo,sp=z_filt_pix)
    return kern



########################################################################

def main():

    datadir=read_params.get_directory()
    iterno=get_iter_no()

    args=sys.argv[1:]
    optimization_algo=filter(lambda temp: temp.startswith('algo='),args)
    if len(optimization_algo)>0: optimization_algo = optimization_algo[0].split("=")[-1]
    else: optimization_algo="conjugate gradient"

    steepest_descent = optimization_algo=='sd' or optimization_algo=='steepest descent'
    conjugate_gradient = optimization_algo=='cg' or optimization_algo=='conjugate gradient'
    LBFGS = optimization_algo=='LBFGS'

    def isfloat(value):
        try:
            float(value)
            return True
        except ValueError:
            return False

    eps = map(float,filter(isfloat,args))
    if eps==[]: eps=[0.1*i for i in xrange(1,7)]

    Rsun=695.9895 # Mm
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

    back=np.loadtxt(read_params.get_solarmodel())

    num_src=get_number_of_sources()

    ############################################################################

    psi_true = np.squeeze(pyfits.getdata(read_params.get_true_psi_filename()))

    ############################################################################
    # Spline
    ############################################################################

    def coeff_to_model(coeffs):
        model = interpolate.bisplev(xspline,zspline,(tx_ref,tz_ref,coeffs,kx_ref,kz_ref))
        model_fullsize = np.zeros_like(psi_true)
        model_fullsize.put(spline_ind_1D,model.T.flatten())
        return model_fullsize

    with np.load(os.path.join(datadir,"true_psi_coeffs.npz")) as f:

        x_cutoff = f["x_cutoff"]; z_cutoff = f["z_cutoff"] # Mm
        xspline_index = abs(x)<x_cutoff
        xspline_index_int = np.where(xspline_index)[0]
        xspline = x[xspline_index]
        zspline_index = z>z_cutoff
        zspline_index_int = np.where(zspline_index)[0]
        zspline = z[zspline_index]

        xspline_index_mgrid = np.array([xspline_index_int for _ in zspline_index_int])
        zspline_index_mgrid = np.array([[zj]*len(xspline_index_int) for zj in zspline_index_int])
        spline_ind_1D = np.ravel_multi_index([zspline_index_mgrid.flatten(),
                                    xspline_index_mgrid.flatten()],(nz,nx))

        xspline_coord_mgrid = np.array([xspline for _ in zspline_index_int])
        zspline_coord_mgrid = np.array([[zj]*len(xspline) for zj in zspline])

        coeff_surf_cutoff_ind = f["coeff_surf_cutoff_ind"]

        tx_ref = f["tx"]
        tz_ref = f["tz"]
        c_ref_above_surface = f["c_upper"]
        c_ref_below_surface = f["c_lower"]
        kx_ref = f["kx"]
        kz_ref = f["kz"]
        spl_c_shape_xz_2D = (len(tx_ref)-kx_ref-1,len(tz_ref)-kz_ref-1)

    ############################################################################
    # Gradient computation
    ############################################################################

    array_shape=(nx,ny,nz)
    totkern_psi=np.zeros(array_shape)
    hess=np.zeros(array_shape)

    def read_model(var='psi',iterno=iterno):
        modelfile = updatedir('model_'+var+'_'+str(iterno).zfill(2)+'_coeffs.npz')
        with np.load(modelfile) as f:
            return dict(f.items())

    def read_grad(var='psi',iterno=iterno):
        return pyfits.getdata(updatedir('gradient_'+var+'_'+str(iterno).zfill(2)+'.fits'))

    def read_update(var='psi',iterno=iterno):
        return pyfits.getdata(updatedir('update_'+var+'_'+str(iterno).zfill(2)+'.fits'))

    def read_kern(var='psi',src=1):
        return fitsread(os.path.join(datadir,'kernel','kernel_'+var+'_'+str(src).zfill(2)+'.fits'))

    for src in xrange(1,num_src+1):

        totkern_psi += read_kern(var='psi',src=src)

        hess+=abs(fitsread(os.path.join(datadir,'kernel','hessian_'+str(src).zfill(2)+'.fits')))

    hess=hess*np.atleast_3d(back[:,2]).transpose(0,2,1)
    hess = hess/abs(hess).max()
    hess[hess<5e-3]=5e-3


    # Filter and smooth
    totkern_psi = filter_and_symmetrize(totkern_psi,hess,z_filt_algo='gaussian',z_filt_pix=2.,sym='asym')

    # Basis kernels
    grad_spline = np.zeros_like(c_ref_below_surface)
    grad = np.squeeze(totkern_psi).T
    # fig=plt.figure()
    # plotdir = os.path.join(datadir,"update","plots_grad_iter"+str(iterno).zfill(2))
    # if not os.path.exists(plotdir): os.makedirs(plotdir)
    for j,_ in enumerate(c_ref_below_surface):
        _,zind = np.unravel_index(j,spl_c_shape_xz_2D)
        if zind>=coeff_surf_cutoff_ind: continue
        c_only_j = np.zeros_like(c_ref_below_surface)
        c_only_j[j] = 1
        bspline_j = coeff_to_model(c_only_j)
        # plt.subplot(131)
        # plt.pcolormesh(x,z,grad,cmap="RdBu")
        # plt.xlim(-50,50)
        # plt.ylim(-8,z[-1])
        # plt.subplot(132)
        # vmax = abs(bspline_j).max()
        # plt.pcolormesh(x,z,bspline_j,cmap="RdBu",vmax=vmax,vmin=-vmax)
        # plt.title("B-spline {:d}".format(j))
        # plt.xlim(-50,50)
        # plt.ylim(-8,z[-1])
        grad_spline[j] = integrate_2D(grad*bspline_j)
        # plt.subplot(133)
        # vmax = abs(grad*bspline_j).max()
        # plt.pcolormesh(x,z,grad*bspline_j,cmap="RdBu",vmax=vmax,vmin=-vmax)
        # plt.xlim(-50,50)
        # plt.ylim(-8,z[-1])
        # plt.title("{:.1E}".format(grad_spline[j]))
        # fig.set_size_inches(11,4)
        # plt.tight_layout()

        # plt.savefig(os.path.join(plotdir,str(j).zfill(4)+".png"))
        # plt.clf()


    plt.figure()
    for i in [(len(tz_ref)-kz_ref-1)*j for j in xrange(0,len(tx_ref)-kx_ref)]:
        plt.axvspan(i+coeff_surf_cutoff_ind,i+(len(tz_ref)-kz_ref-1),
        facecolor="honeydew",edgecolor="lightsage")
        plt.axvline(i,ls="dotted",color="black")
    plt.plot(grad_spline,'o-',markersize=3)
    plt.xlim(0,(len(tz_ref)-kz_ref-1)*(len(tx_ref)-kx_ref-1))
    plt.savefig(os.path.join(datadir,"update","grad_"+str(iterno).zfill(2)+".png"))

    #~ Write out gradients for this iteration
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pyfits.writeto(updatedir('gradient_psi_'+str(iterno).zfill(2)+'.fits'),grad_spline,clobber=True)

    ############################################################################
    # Optimization
    ############################################################################

    #~ Get update direction based on algorithm of choice
    if (iterno==0) or steepest_descent:

        def sd_update(var='psi'):
            grad = read_grad(var=var,iterno=iterno)
            update = -grad
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                pyfits.writeto(updatedir('update_'+var+'_'+str(iterno).zfill(2)+'.fits'),update,clobber=True)
                # fitswrite(updatedir('update_'+var+'_'+str(iterno).zfill(2)+'.fits'),update)

            print "Steepest descent"

        sd_update(var='psi')

        if iterno > 0: print 'Forcing steepest descent'

    elif conjugate_gradient:

        def get_beta(grad,lastgrad,lastupdate):
            #~ Choose algo
            polak_ribiere = True
            hestenes_stiefel = True and (not polak_ribiere)

            if polak_ribiere:
                print "Conjugate gradient, Polak Ribiere method"
                beta = np.sum(grad*(grad - lastgrad))
                beta/= np.sum(lastgrad**2.)

            elif hestenes_stiefel:
                print "Conjugate gradient, Hestenes Stiefel method"
                beta = np.sum(grad*(grad - lastgrad))
                beta/= np.sum((grad - lastgrad)*lastupdate)

            print "beta",beta
            if beta==0: print "Conjugate gradient reduces to steepest descent"
            elif beta<0:
                print "Stepping away from previous update direction"

            return beta

        def cg_update(var='psi'):
            grad_k = read_grad(var=var,iterno=iterno)
            grad_km1 = read_grad(var=var,iterno=iterno-1)
            p_km1 = read_update(var=var,iterno=iterno-1)

            beta_k = get_beta(grad_k,grad_km1,p_km1)
            update=-grad_k +  beta_k*p_km1

            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                pyfits.writeto(updatedir('update_'+var+'_'+str(iterno).zfill(2)+'.fits'),update,clobber=True)

        cg_update(var='psi')

    elif LBFGS:

        m=4
        k = iterno

        def LBFGS_update(var='psi'):
            model_k= read_model(var=var,iterno=iterno)
            grad_k = read_grad(var=var,iterno=iterno)
            model_iplus1 = model_k
            grad_iplus1 = grad_k
            q = grad_k

            for i in xrange(k-1,k-m-1,-1):

                if var+'_y_'+str(i) in LBFGS_data.keys():
                    y_i=LBFGS_data[var+'_y_'+str(i)]
                else:
                    grad_i= read_grad(var=var,iterno=i)
                    y_i = grad_iplus1 - grad_i
                    LBFGS_data[var+'_y_'+str(i)] = y_i
                    grad_iplus1 = grad_i

                if var+'_s_'+str(i) in LBFGS_data.keys():
                    s_i=LBFGS_data[var+'_s_'+str(i)]
                else:
                    model_i= read_grad(var=var,iterno=i)
                    s_i = model_iplus1 - model_i
                    LBFGS_data[var+'_s_'+str(i)] = s_i
                    model_iplus1 = model_i

                if var+'_rho_'+str(i) in LBFGS_data.keys():
                    rho_i = LBFGS_data[var+'_rho_'+str(i)]
                else:
                    rho_i = 1/np.sum(y_i*s_i)
                    LBFGS_data[var+'_rho_'+str(i)] = rho_i

                if var+'_alpha_'+str(i) in LBFGS_data.keys():
                    alpha_i = LBFGS_data[var+'_alpha_'+str(i)]
                else:
                    alpha_i = rho_i*np.sum(s_i*q)
                    LBFGS_data[var+'_alpha_'+str(i)] = rho_i

                q = q - alpha_i * y_i

            model_iplus1=None; grad_iplus1=None; model_i=None; grad_i=None

            s_kminus1 = LBFGS_data[var+'_s_'+str(i)]
            y_kminus1 = LBFGS_data[var+'_y_'+str(i)]

            H0_k = np.dot(s_kminus1*y_kminus1)/np.sum(y_kminus1**2)
            r = np.dot(H0_k,q)

            for i in xrange(k-m,k):

                rho_i = LBFGS_data[var+'_rho_'+str(i)]
                y_i = LBFGS_data[var+'_y_'+str(i)]
                s_i = LBFGS_data[var+'_s_'+str(i)]
                alpha_i = LBFGS_data[var+'_alpha_'+str(i)]

                beta = rho_i* np.sum(y_i*r)
                r = r + s_i*(alpha_i - beta)

            for i in xrange(0,k-m):
                if var+'_s_'+str(i) in LBFGS_data: del LBFGS_data[var+'_s_'+str(i)]
                if var+'_y_'+str(i) in LBFGS_data: del LBFGS_data[var+'_y_'+str(i)]
                if var+'_rho_'+str(i) in LBFGS_data: del LBFGS_data[var+'_rho_'+str(i)]
                if var+'_alpha_'+str(i) in LBFGS_data: del LBFGS_data[var+'_alpha_'+str(i)]


            fitswrite(updatedir('update_'+var+'_'+str(iterno).zfill(2)+'.fits'),r)

        LBFGS_data_file = updatedir('LBFGS_data.npz')
        try:
            LBFGS_data={}
            data=np.load(LBFGS_data_file)
            for key,val in data.items(): LBFGS_data[key]=val
            data=None
        except IOError: LBFGS_data={}

        LBFGS_update(var='psi')

        np.savez(LBFGS_data_file,**LBFGS_data)

    ############################################################################

    large_x_cutoff = 40
    deep_z_cutoff = -5

    def create_ls_model(i=0,var='psi',eps=0,kind='linear'):
        model = read_model(var=var,iterno=iterno)
        model_back = model["back"]
        c_model_below_surface = model["c_lower"]
        if i==0:
            fig=plt.figure()
            plt.plot(c_ref_above_surface+c_ref_below_surface,'o-',markersize=3)
            plt.plot(c_ref_above_surface+c_model_below_surface,'o-',markersize=3)
            for i in [(len(tz_ref)-kz_ref-1)*j for j in xrange(0,len(tx_ref)-kx_ref)]:
                plt.axvspan(i+coeff_surf_cutoff_ind,i+(len(tz_ref)-kz_ref-1),
                facecolor="honeydew",edgecolor="lightsage")
                plt.axvline(i,ls="dotted",color="black")
            plt.xlim(0,(len(tz_ref)-kz_ref-1)*(len(tx_ref)-kx_ref-1))
            plt.savefig(os.path.join(datadir,"update","coeffs_1D_"+str(iterno).zfill(2)+".png"))

            plt.clf()
            plt.subplot(131)
            vmax = abs(c_ref_above_surface+c_ref_below_surface).max()
            plt.pcolormesh((c_ref_above_surface+c_ref_below_surface).reshape(spl_c_shape_xz_2D).T,
            cmap="RdBu",vmax=vmax,vmin=-vmax)
            plt.xlim(0,spl_c_shape_xz_2D[0])
            plt.ylim(0,spl_c_shape_xz_2D[1])
            plt.axhline(10,ls="dotted",color="black")
            plt.title("True")

            plt.subplot(132)

            plt.pcolormesh((c_ref_above_surface+c_model_below_surface).reshape(spl_c_shape_xz_2D).T,
            cmap="RdBu",vmax=vmax,vmin=-vmax)
            plt.axhline(10,ls="dotted",color="black")
            plt.xlim(0,spl_c_shape_xz_2D[0])
            plt.ylim(0,spl_c_shape_xz_2D[1])
            plt.title("Iterated")

            plt.subplot(133)

            plt.pcolormesh(
            (c_ref_below_surface-c_model_below_surface).reshape(spl_c_shape_xz_2D).T,
            cmap="RdBu")

            plt.axhline(10,ls="dotted",color="black")
            plt.xlim(0,spl_c_shape_xz_2D[0])
            plt.ylim(0,spl_c_shape_xz_2D[1])
            plt.title("Difference")

            fig.set_size_inches(11,4)
            plt.tight_layout()
            plt.savefig(os.path.join(datadir,"update","coeffs_2D_"+str(iterno).zfill(2)+".png"))


        update = read_update(var=var,iterno=iterno)

        updatemax=update.max()
        if updatemax!=0: update/=updatemax

        if kind=='linear':
            model_scale = rms(c_model_below_surface)
            if model_scale ==0: model_scale=c_ref_above_surface.max()
            c_model_below_surface = c_model_below_surface+eps*update*model_scale

        elif kind=='exp':
            c_model_below_surface= c_model_below_surface*(1+eps*update)


        model_below_surface = coeff_to_model(c_model_below_surface)
        cutoff_x = 1/(1+np.exp((abs(x)-large_x_cutoff)/3))
        cutoff_z = 1/(1+np.exp(-(z - deep_z_cutoff)/1))
        model_below_surface *= cutoff_x[None,:]*cutoff_z[:,None]
        model_above_surface = coeff_to_model(c_ref_above_surface)
        model = model_below_surface + model_above_surface + model_back
        model = model[:,np.newaxis,:]
        model = np.transpose(model,(2,1,0))

        return model,c_model_below_surface,model_back

    #~ Create models for linesearch
    for i,eps_i in enumerate(eps):
        lsmodel,ls_coeffs_below_surface,back = create_ls_model(i=i,var='psi',eps=eps_i,kind='linear')
        fitswrite(updatedir('test_psi_'+str(i+1)+'.fits'), lsmodel)

        plt.subplot(121)
        plt.pcolormesh(x,z,psi_true,cmap="RdBu_r",vmax=abs(psi_true).max(),vmin=-abs(psi_true).max())
        plt.colorbar()
        plt.xlim(-70,70)
        plt.ylim(-7,z[-1])
        plt.subplot(122)
        arr_to_plot = np.squeeze(lsmodel).T-back
        plt.pcolormesh(x,z,arr_to_plot,cmap="RdBu_r",
        vmax=abs(arr_to_plot).max(),vmin=-abs(arr_to_plot).max())
        plt.colorbar()
        plt.axvline(-x_cutoff,color="black",ls="dotted")
        plt.axvline(x_cutoff,color="black",ls="dotted")
        plt.axvline(-large_x_cutoff,color="brown",ls="dotted")
        plt.axvline(large_x_cutoff,color="brown",ls="dotted")
        plt.axhline(z_cutoff,color="black",ls="dotted")
        plt.axhline(deep_z_cutoff,color="brown",ls="dotted")
        plt.xlim(-70,70)
        plt.ylim(-7,z[-1])

        plt.gcf().set_size_inches(8,3)
        plt.tight_layout()
        plt.savefig(updatedir("test_psi_"+str(i)+".png"))
        plt.clf()


        np.savez(updatedir('test_psi_'+str(i+1)+'_coeffs.npz'),c_lower=ls_coeffs_below_surface,back=back)

    #~ Update epslist
    epslist_path = os.path.join(datadir,"epslist.npz")
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
            arr = np.zeros((4,6))
            arr[0] = eps
            epslist = {i:epslist[i] for i in iter_done_list}
            epslist[str(iterno)] = arr
    except IOError:
        arr = np.zeros((4,6))
        arr[0] = eps
        epslist = {str(iterno):arr}

    #~ np.set_printoptions(precision=5)
    #~ print "Updated epslist in grad"
    #~ print epslist[str(iterno)]
    np.savez(epslist_path,**epslist)


########################################################################

if __name__ == "__main__":
    main()
