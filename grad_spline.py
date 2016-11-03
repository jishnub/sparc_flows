from __future__ import division
import numpy as np
from scipy.ndimage.filters import gaussian_filter1d
from scipy import interpolate,integrate
from scipy.special import j1,j0,jn
def j2(z): return jn(2,z)
def j1prime(z): return 0.5*(j0(z)-j2(z))
import os,fnmatch,sys
import pyfits
import warnings
import read_params
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt,ticker

########################################################################

class supergranule():
    def __init__(self,**kwargs):
        self.__dict__.update(kwargs)

DH13 = supergranule(R = 15,
                    k = 2*np.pi/30,
                    sigmaz = 0.912,
                    z0 = -2.3,
                    v0 = 240)

########################################################################

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

                s=interpolate.UnivariateSpline(z_ind,arrz,s=sp)
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
    if len(optimization_algo)>0:
        optimization_algo = optimization_algo[0].split("=")[-1]
    else: optimization_algo="cg"

    steepest_descent = optimization_algo=='sd'
    conjugate_gradient = optimization_algo=='cg'
    LBFGS = optimization_algo=='LBFGS'
    BFGS = optimization_algo=='BFGS'

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

    def read_model(var='psi',iterno=iterno):
        modelfile = updatedir('model_'+var+'_'+str(iterno).zfill(2)+'_coeffs.npz')
        with np.load(modelfile) as f:
            return dict(f.items())

    def read_grad(var='psi',iterno=iterno):
        # return pyfits.getdata(updatedir('gradient_'+var+'_'+str(iterno).zfill(2)+'.fits'))
        modelfile = updatedir('gradient_'+var+'_'+str(iterno).zfill(2)+'.npz')
        with np.load(modelfile) as f:
            return dict(f.items())

    def read_update(var='psi',iterno=iterno):
        # return pyfits.getdata(updatedir('update_'+var+'_'+str(iterno).zfill(2)+'.fits'))
        modelfile = updatedir('update_'+var+'_'+str(iterno).zfill(2)+'.npz')
        with np.load(modelfile) as f:
            return dict(f.items())

    def read_BFGS_hessian(var='psi',iterno=iterno):
        # return pyfits.getdata(updatedir('update_'+var+'_'+str(iterno).zfill(2)+'.fits'))
        modelfile = updatedir('BFGS_hessian_'+var+'_'+str(iterno).zfill(2)+'.npz')
        with np.load(modelfile) as f:
            return dict(f.items())

    def read_kern(var='psi',src=1):
        return fitsread(os.path.join(datadir,'kernel','kernel_'+var+'_'+str(src).zfill(2)+'.fits'))

    psi_true = np.squeeze(pyfits.getdata(read_params.get_true_psi_filename()))

    ############################################################################
    # Spline
    ############################################################################

    def coeff_to_model(tck_z,tck_R):
        f0_x = np.sign(x)*j1(DH13.k*abs(x))*np.exp(-abs(x)/DH13.R)
        f0_x_max = f0_x.max()
        f0_x/=f0_x_max

        h_z=interpolate.splev(z,tck_z,ext=1)

        # f1_x is the derivative of f0_x wrt R
        f1_x = x*np.exp(-abs(x)/DH13.R)/DH13.R**2*(j1(DH13.k*abs(x))-np.pi*j1prime(DH13.k*abs(x)))
        f1_x/=f0_x_max
        if tck_R[1] is not None:
            R1_z = interpolate.splev(z,tck_R,ext=1)
        else:
            R1_z = np.zeros_like(h_z)


        return (f0_x[None,:]+f1_x[None,:]*R1_z[:,None])*h_z[:,None]


    f= dict(np.load(os.path.join(datadir,"true_psi_coeffs.npz")))
    coeff_surf_cutoff_ind = f.get("c_surf_cutoff").item()
    tR = f.get("tR")
    tz = f.get("tz")
    cz_ref_top = f.get("cz_top")
    cz_ref_bot = f.get("cz_bot")
    kR = f.get("kR")
    kz = f.get("kz")

    z_spl_cutoff = f.get("z_spline_cutoff")
    R_surf_cutoff = f.get("R_surf_cutoff")

    class spline_basis_coeffs():
        def __init__(self,true_coeffs,iter_coeffs,low_ind=0,high_ind=None):
            self.true = true_coeffs
            self.iterated = iter_coeffs
            self.low_ind = low_ind
            if high_ind is not None:
                self.high_ind = high_ind if high_ind>0 else self.true.size+high_ind
            else:
                self.high_ind = None
            self.grad = None
            self.update = None
            self.size = self.iterated.size if self.iterated is not None else None


        def get_true(self,low_ind=None,high_ind=None):
            low_ind = self.low_ind if low_ind is None else low_ind
            high_ind = self.high_ind if high_ind is None else high_ind
            return self.true[low_ind:high_ind]

        def get_iterated(self,low_ind=None,high_ind=None):
            low_ind = self.low_ind if low_ind is None else low_ind
            high_ind = self.high_ind if high_ind is None else high_ind
            return self.iterated[low_ind:high_ind]

        def get_update(self,low_ind=None,high_ind=None):
            low_ind = self.low_ind if low_ind is None else low_ind
            high_ind = self.high_ind if high_ind is None else high_ind
            return self.update[low_ind:high_ind]

        def get_range(self):
            return np.array(range(self.low_ind,self.high_ind))

    iter_model = read_model()

    coeffs = {
    "z":spline_basis_coeffs(true_coeffs=cz_ref_bot,iter_coeffs=iter_model["z"],
    low_ind=kz+1,high_ind=coeff_surf_cutoff_ind)}

    if iter_model.get("cR") is not None:
        coeffs["R"]=spline_basis_coeffs(true_coeffs=np.zeros_like(iter_model.get("R")),
        iter_coeffs=iter_model.get("R"),low_ind=kR+1,
        high_ind=R_surf_cutoff)

    ############################################################################
    # Gradient computation
    ############################################################################

    array_shape=(nx,ny,nz)
    totkern_psi=np.zeros(array_shape)
    hess=np.zeros(array_shape)

    for src in xrange(1,num_src+1):

        totkern_psi += read_kern(var='psi',src=src)

        hess+=abs(fitsread(os.path.join(datadir,'kernel','hessian_'+str(src).zfill(2)+'.fits')))

    hess=hess*np.atleast_3d(back[:,2]).transpose(0,2,1)
    hess = hess/abs(hess).max()
    hess[hess<5e-3]=5e-3

    # Filter and smooth
    totkern_psi = filter_and_symmetrize(totkern_psi,hess,
                    z_filt_algo='gaussian',z_filt_pix=2.,sym='asym')

    kernel = np.squeeze(totkern_psi).T

    large_x_cutoff = 40
    cutoff_x = 1/(1+np.exp((abs(x)-large_x_cutoff)/5))
    kernel = kernel*cutoff_x[None,:]

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pyfits.writeto(updatedir("grad_psi_"+str(iterno).zfill(2)+".fits"),
                        kernel,clobber=True)

    # Plot gradient (kernel)
    f=plt.figure()
    plt.subplot(121)
    plt.pcolormesh(x,z,psi_true,cmap="RdBu_r")
    plt.title("True psi",fontsize=16)
    plt.xlim(-50,50)
    plt.ylim(-6,z[-1])
    plt.xlabel("x (Mm)",fontsize=16)
    plt.ylabel("z (Mm)",fontsize=16)

    plt.subplot(122)
    plt.pcolormesh(x,z,kernel/abs(kernel).max(),cmap="RdBu_r",vmax=1,vmin=-1)
    plt.title("Gradient",fontsize=16)
    plt.xlim(-50,50)
    plt.ylim(-6,z[-1])
    plt.xlabel("x (Mm)",fontsize=16)
    plt.ylabel("z (Mm)",fontsize=16)

    f.set_size_inches(8,3.5)
    plt.tight_layout()
    plt.savefig(os.path.join(datadir,"update","grad_"+str(iterno).zfill(2)+".png"))

    # compute basis coefficient gradients
    def compute_grad_basis(coeffs):

        f0_x = np.sign(x)*j1(DH13.k*abs(x))*np.exp(-abs(x)/DH13.R)
        f0_x_max = f0_x.max()
        f0_x/=f0_x_max
        hs = interpolate.splev(z,(tz,coeffs["z"].iterated+cz_ref_top,kz),ext=1)
        g_c = f0_x[None,:]

        grad = dict.fromkeys(coeffs.keys())
        for key in grad:
            grad[key] = np.zeros(coeffs[key].size)

        if "R" in coeffs.keys() and coeffs["R"].iterated is not None:
            # f1_x is the derivative of f0_x wrt R
            f1_x = x*np.exp(-abs(x)/DH13.R)/DH13.R**2*(j1(DH13.k*abs(x))-np.pi*j1prime(DH13.k*abs(x)))
            f1_x/=f0_x_max
            R1 = interpolate.splev(z,(tR,coeffs["R"].iterated,kR),ext=1)
            g_c = g_c + f1_x[None,:]*R1[:,None]
            g_R = f1_x[None,:]*hs[:,None]

            for j in coeffs["R"].get_range():
                qj = np.zeros_like(coeffs["R"].iterated)
                qj[j] = 1
                qz = interpolate.splev(z,(tR,qj,kR),ext=1)
                grad["R"][j] = integrate_2D(kernel*qz[:,None]*g_R)

        if not os.path.exists(os.path.join(datadir,"update","basis_kernels")):
            os.makedirs(os.path.join(datadir,"update","basis_kernels"))
        for j in coeffs["z"].get_range():
            bj = np.zeros_like(coeffs["z"].iterated)
            bj[j] = 1
            bz = interpolate.splev(z,(tz,bj,kz),ext=1)
            grad["z"][j] = integrate_2D(kernel*bz[:,None]*g_c)

            plt.figure()
            plt.subplot(1,2,1)
            plt.pcolormesh(x,z,kernel/abs(kernel).max(),cmap="RdBu_r",vmax=1,vmin=-1)
            plt.title("Gradient",fontsize=16)
            plt.xlim(-50,50)
            plt.ylim(-6,z[-1])
            plt.xlabel("x (Mm)",fontsize=16)
            plt.ylabel("z (Mm)",fontsize=16)

            plt.subplot(1,2,2)
            plt.pcolormesh(x,z,kernel*bz[:,None]*g_c/abs(kernel*bz[:,None]*g_c).max(),
                            cmap="RdBu_r",vmax=1,vmin=-1)
            plt.title("Integral: {:.1e}".format(grad["z"][j]),fontsize=16)
            plt.xlim(-50,50)
            plt.ylim(-6,z[-1])
            plt.xlabel("x (Mm)",fontsize=16)
            plt.ylabel("z (Mm)",fontsize=16)

            plt.gcf().set_size_inches(8,3.5)
            plt.tight_layout()
            plt.savefig(os.path.join(datadir,"update",
                        "basis_kernels","{:02d}.png".format(j)))
            plt.clf()

        return grad

    #~ Write out gradients for this iteration
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        np.savez(updatedir('gradient_psi_'+str(iterno).zfill(2)+'.npz'),
        **compute_grad_basis(coeffs))

    ############################################################################
    # Optimization
    ############################################################################

    #~ Get update direction based on algorithm of choice
    if (iterno==0) or steepest_descent:

        def sd_update(var='psi'):
            grad = read_grad(var='psi',iterno=iterno)
            update = {p:-grad[p] for p in grad.keys()}
            updatefile = updatedir('update_'+var+'_'+str(iterno).zfill(2)+'.npz')
            np.savez(updatefile,**update)
            print "Steepest descent"

        sd_update(var='psi')

        if iterno > 0: print 'Forcing steepest descent'

    elif conjugate_gradient:

        def get_beta(grad,lastgrad,lastupdate):
            #~ Choose algo
            polak_ribiere = True
            hestenes_stiefel = True and (not polak_ribiere)

            np.seterr(divide="raise")
            try:
                beta = np.sum(grad*(grad - lastgrad))
                if polak_ribiere:
                    beta/= np.sum(lastgrad**2.)
                elif hestenes_stiefel:
                    beta/= np.sum((grad - lastgrad)*lastupdate)
            except FloatingPointError:
                beta=0

            return beta

        def cg_update(var='psi'):
            grad_k = read_grad(var=var,iterno=iterno)
            grad_km1 = read_grad(var=var,iterno=iterno-1)
            p_km1 = read_update(var=var,iterno=iterno-1)

            beta_k = {}
            update = {}
            for param in grad_k.keys():
                beta_k[param] = get_beta(grad_k[param],grad_km1[param],p_km1[param])
                update[param] = -grad_k[param] +  beta_k[param]*p_km1[param]

            print "beta",dict((k,round(v,2)) for k,v in beta_k.iteritems())

            updatefile = updatedir('update_'+var+'_'+str(iterno).zfill(2)+'.npz')
            np.savez(updatefile,**update)

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

    elif BFGS:
        def BFGS_update(var='psi'):
            print "Using BFGS"
            grad_k = read_grad(var=var,iterno=iterno)
            grad_km1 = read_grad(var=var,iterno=iterno-1)

            model_k  = read_model(var=var,iterno=iterno)
            model_km1  = read_model(var=var,iterno=iterno-1)

            try:
                Hkm1 = read_BFGS_hessian(var=var,iterno=iterno)
            except IOError:
                Hkm1 = {}
                for key in grad_k.keys():
                    Hkm1[key] = np.identity(grad_k[key].size)


            Hk = {}
            update = {}
            for key in grad_k.keys():
                y_km1 = grad_k[key] - grad_km1[key]
                s_km1 = model_k[key] - model_km1[key]

                if (not y_km1.any()) or (not s_km1.any()):
                    Hk[key] = Hkm1[key]
                else:

                    rho_km1 = 1/y_km1.dot(s_km1)

                    left = np.identity(s_km1.size) - rho_km1*np.outer(s_km1,y_km1)
                    right = np.identity(s_km1.size) - rho_km1*np.outer(y_km1,s_km1)

                    Hk[key] = left.dot(Hkm1[key]).dot(right) + rho_km1*np.outer(s_km1,s_km1)

                update[key] = -Hk[key].dot(grad_k[key])

            updatefile = updatedir('update_'+var+'_'+str(iterno).zfill(2)+'.npz')
            np.savez(updatefile,**update)

            hessfile = updatedir('BFGS_hessian_'+var+'_'+str(iterno).zfill(2)+'.npz')
            np.savez(hessfile,**Hk)

        BFGS_update(var='psi')


    ############################################################################

    # Add update to coeff dictionary
    update = read_update(var='psi',iterno=iterno)
    for param,coeff in coeffs.items():
        coeff.update = update[param]

    ############################################################################

    deep_z_cutoff = -4

    b_i_surf = np.zeros_like(coeffs.get("z").true)
    for i in xrange(b_i_surf.size):
        c_i = np.zeros_like(b_i_surf)
        c_i[i] = 1
        b_i_surf[i] = interpolate.splev(deep_z_cutoff,(tz,c_i,kz))

    c_deep_z_cutoff_index = b_i_surf.argmax()

    coeffs["z"].update *= 1/(1+np.exp(-(np.arange(b_i_surf.size)-c_deep_z_cutoff_index)/1))

    ############################################################################

    # Plot update coefficients

    model_grad_coeffs_fig=plt.figure()
    ax = [plt.subplot(len(coeffs.keys()),1,i) for i in xrange(1,len(coeffs.keys())+1)]

    for ind,(param,coeff) in enumerate(coeffs.items()):

        # True model coefficients
        ax[ind].plot(coeff.get_range(),coeff.get_true(),
        'o-',markersize=4,color="teal",label="True")
        # Iterated model coefficients
        ax[ind].plot(coeff.get_range(),coeff.get_iterated(),
        'o-',markersize=4,color="brown",label="Iter")

        # Read grad and plot in twin axis

        ax2 = ax[ind].twinx()
        ax2.bar(coeff.get_range()-0.3,coeff.get_update(),width=0.6,bottom=0,
        color="goldenrod",label="Update",edgecolor="peru",alpha=0.4)

        s = ticker.ScalarFormatter()
        s.set_scientific(True)
        s.set_powerlimits((0,0))
        ax2.yaxis.set_major_formatter(s)
        ax[ind].yaxis.set_major_formatter(s)



        ax[ind].set_title(param,fontsize=16)

        handles1,labels1 = ax[ind].get_legend_handles_labels()
        handles2,labels2 = ax2.get_legend_handles_labels()
        ax[ind].legend(handles1+handles2,labels1+labels2,loc="upper right")

    sp_ind_z = map(lambda x: x.get_title(),ax).index("z")

    ax[sp_ind_z].plot(range(coeff_surf_cutoff_ind-1,cz_ref_top.size-kz-1),
    (cz_ref_top + cz_ref_bot)[coeff_surf_cutoff_ind-1:cz_ref_top.size-kz-1],
    'o--',markersize=4,color="teal")

    ax[sp_ind_z].plot(range(coeff_surf_cutoff_ind-1,cz_ref_top.size-kz-1),
    (cz_ref_top + coeffs["z"].iterated)[coeff_surf_cutoff_ind-1:cz_ref_top.size-kz-1],
    '--',color="brown")

    ax[sp_ind_z].axvspan(0,kz+0.5,color="thistle",zorder=0)
    ax[sp_ind_z].axvspan(cz_ref_top.size-kz-1.5,cz_ref_top.size-1,color="thistle",
    zorder=0,label="set to zero")
    ax[sp_ind_z].axvspan(coeff_surf_cutoff_ind-0.5,cz_ref_top.size-kz-1.5,color="lightgrey",
    zorder=0,label="clamped")

    ax[sp_ind_z].plot(range(kz+2),coeffs["z"].iterated[:kz+2],'--',color="brown",zorder=1)
    ax[sp_ind_z].plot(range(kz+2),coeffs["z"].true[:kz+2],'--',color="teal",zorder=1)
    ax[sp_ind_z].plot(range(kz+1),coeffs["z"].true[:kz+1],marker="o",mfc="thistle",
    ls="None",zorder=2)

    ax[sp_ind_z].plot(range(coeffs["z"].iterated.size-(kz+2),coeffs["z"].iterated.size),
    (cz_ref_top + coeffs["z"].true)[-(kz+2):],'--',color="brown")
    ax[sp_ind_z].plot(range(coeffs["z"].true.size-(kz+2),coeffs["z"].true.size),
    (cz_ref_top + coeffs["z"].true)[-(kz+2):],'--',color="teal")
    ax[sp_ind_z].plot(range(coeffs["z"].true.size-(kz+1),coeffs["z"].true.size),
    (cz_ref_top + coeffs["z"].true)[-(kz+1):],marker="o",mfc="thistle",
    ls="None",zorder=2)

    for ax_i in ax:
        ax_i.margins(x=0.1)
        ax_i.legend(loc="best")

    plt.tight_layout()
    plt.savefig(os.path.join(datadir,"update","coeffs_1D_"+str(iterno).zfill(2)+".png"))

    ############################################################################



    def create_ls_model(var='psi',eps=0,kind='linear'):

        for param in update.keys():
            if abs(update[param]).max()!=0:
                update[param]/=abs(update[param]).max()

        ls_cz = coeffs.get("z").iterated
        ls_cR = coeffs.get("R")
        if ls_cR is not None: ls_cR = ls_cR.iterated

        if kind=='linear':
            cz_scale = rms(coeffs["z"].get_iterated())
            if cz_scale==0: cz_scale=cz_ref_top.max()
            ls_cz += eps*update["z"]*cz_scale

            if ls_cR is not None:
                cR_scale = rms(coeffs.get("R").get_iterated())
                if cR_scale==0: cR_scale=0.1
                ls_cR += eps*update.get("R")*cR_scale

        elif kind=='exp':
            ls_cz *= 1+eps*update
            if ls_cR is not None:
                ls_cR *= 1+eps*update

        lsmodel = coeff_to_model((tz,ls_cz+cz_ref_top,kz),(tR,ls_cR,kR))

        # lsmodel *= cutoff_x[None,:]
        lsmodel += iter_model["back"]

        lsmodel = lsmodel[:,np.newaxis,:]
        lsmodel = np.transpose(lsmodel,(2,1,0))

        coeffdict = {"z":ls_cz,"back":iter_model["back"]}
        if ls_cR is not None:
            coeffdict["R"] = ls_cR
        return lsmodel,coeffdict

    #~ Create models for linesearch
    for i,eps_i in enumerate(eps):
        lsmodel,model_coeffs = create_ls_model(var='psi',eps=eps_i,kind='linear')
        fitswrite(updatedir('test_psi_'+str(i+1)+'.fits'), lsmodel)
        np.savez(updatedir('test_psi_'+str(i+1)+'_coeffs.npz'),**model_coeffs)


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
