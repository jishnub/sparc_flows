from __future__ import division
import numpy as np
from scipy.interpolate import UnivariateSpline
from scipy.ndimage.filters import gaussian_filter1d
import os,fnmatch,sys
import pyfits
import warnings
import read_params

#######################################################################

def fitsread(f): 
    try:
        arr=pyfits.getdata(f)
    except IOError:
        raise IOError
    # If it is 2D, make it 3D. This adds an extra dimension at the end.
    # Bring it to the middle to match the general trend
    # Dimension will be (nz,ny,nx) after this
    if len(arr.shape)==2: arr=np.atleast_3d(arr).transpose(0,2,1)
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
    updatedir=os.path.join(datadir,"update")
    # Count the number of misfit_xx files
    return len(fnmatch.filter(os.listdir(updatedir),'misfit_[0-9][0-9]'))-1
    
def get_number_of_sources():
    # Count the number of lines in master.pixels
    # Safest to use loadtxt, since it removes newline correctly
    return len(np.loadtxt(os.path.join(datadir,'master.pixels'),ndmin=1))

def filterx(kern):
    temp=np.zeros_like(kern)
    nx,ny,nz=kern.shape
    
    x=np.fft.fftfreq(nx)*nx
    smooth_x = 4
    x = np.exp(-x**2./(2.*smooth_x**2.))
    
    filtx=np.fft.rfft(x)
    filtx=np.atleast_3d(filtx).transpose(1,0,2)
    
    temp=np.fft.rfft(kern,axis=0)
    temp*=filtx
    temp=np.fft.irfft(temp,axis=0).real
    
    kern[:]=temp[:]

def filterz(arr,algo='spline',sp=1.0):
    nx,ny,nz=arr.shape
    z=np.arange(nz)
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
        for x in xrange(nx):
            for y in xrange(ny):
                arrz=arr[x,y]
                arrzmax=arrz.max()
                arrz=arrz/arrzmax
                
                s=UnivariateSpline(z,arrz,s=sp)
                arr[x,y]=s(z)*arrzmax
                
    elif algo=='gaussian':
        arr[:]=gaussian_filter1d(arr,sigma=sp,axis=-1)
    
def antisymmetrize(arr):   arr[:]=0.5*(arr[:]-arr[::-1])
def symmetrize(arr):   arr[:]=0.5*(arr[:]+arr[::-1])

def updatedir(filename): return os.path.join(datadir,'update',filename)

def rms(arr): return np.sqrt(np.sum(arr**2)/np.prod(arr.shape)) 

def filter_and_symmetrize(totkern,hess,sym=None,z_filt_algo='gaussian',z_filt_pix=0.3):
    kern = totkern/hess
    filterx(kern)
    if sym=='sym':
        symmetrize(kern)
    elif sym=='asym':
        antisymmetrize(kern)
    filterz(kern,algo=z_filt_algo,sp=z_filt_pix)
    return kern
    

########################################################################

datadir=read_params.get_directory()
iterno=get_iter_no()

def main():

    args=sys.argv[1:]
    algo=filter(lambda x: x.startswith('algo='),args)[0].split("=")[-1]
    
    steepest_descent=False
    conjugate_gradient = False
    LBFGS = False
    if algo=='sd' or algo=='steepest descent': steepest_descent=True
    elif algo=='cg' or algo=='conjugate gradient': conjugate_gradient=True
    elif algo=='LBFGS': LBFGS=True
    
    def isfloat(value):
        try:
            float(value)
            return True
        except ValueError:
            return False
    
    eps = map(float,filter(isfloat,args))
    

    Rsun=695.9895 # Mm

    #~ Get shape
    nx = read_params.get_nx()
    ny = 1
    nz = read_params.get_nz()

    back=np.loadtxt(read_params.get_solarmodel())

    num_src=get_number_of_sources()

    array_shape=(nx,ny,nz)
    totkern_c=np.zeros(array_shape)
    totkern_psi=np.zeros(array_shape)
    totkern_vx=np.zeros(array_shape)
    totkern_vz=np.zeros(array_shape)
    hess=np.zeros(array_shape)
    
    
    enf_cont=read_params.get_enforced_continuity()
    cont_var=read_params.get_continuity_variable()
    psi_cont = False
    vx_cont = False
    vz_cont = False
    if enf_cont:
        if cont_var == 'psi': psi_cont = True 
        elif cont_var == 'vx': vx_cont = True 
        elif cont_var == 'vz': vz_cont = True 
        else: 
            print "Continuity variable unknown, check params.i"
            quit()

    for src in xrange(1,num_src+1):
        
        kern=fitsread(os.path.join(datadir,'kernel','kernel_c_'+str(src).zfill(2)+'.fits'))
        totkern_c+=kern
        
        if enf_cont and psi_cont:
            kern=fitsread(os.path.join(datadir,'kernel','kernel_psi_'+str(src).zfill(2)+'.fits'))
            totkern_psi+=kern
        
        elif enf_cont and vx_cont:
            kern=fitsread(os.path.join(datadir,'kernel','kernel_vx_'+str(src).zfill(2)+'.fits'))
            totkern_vx+=kern
            
        elif enf_cont and vz_cont:
            kern=fitsread(os.path.join(datadir,'kernel','kernel_vz_'+str(src).zfill(2)+'.fits'))
            totkern_vz+=kern
        else:
            kern=fitsread(os.path.join(datadir,'kernel','kernel_vx_'+str(src).zfill(2)+'.fits'))
            totkern_vx+=kern
            kern=fitsread(os.path.join(datadir,'kernel','kernel_vz_'+str(src).zfill(2)+'.fits'))
            totkern_vz+=kern
        
        kern=fitsread(os.path.join(datadir,'kernel','hessian_'+str(src).zfill(2)+'.fits'))
        hess+=abs(kern)

    hess=hess*np.atleast_3d(back[:,2]).transpose(0,2,1)
    hess = hess/abs(hess).max()
    hess[hess<5e-3]=5e-3

    #~ Sound speed kernel    
    totkern_c = filter_and_symmetrize(totkern_c,hess,z_filt_algo='smooth',z_filt_pix=0.3,sym=None)

    #~ Vector potential/stream function kernel    
    totkern_psi = filter_and_symmetrize(totkern_psi,hess,z_filt_algo='gaussian',z_filt_pix=8.,sym='asym')
    
    #~ Velocity kernels    
    totkern_vx = filter_and_symmetrize(totkern_vx,hess,z_filt_algo='gaussian',z_filt_pix=2.,sym='asym')
    totkern_vz = filter_and_symmetrize(totkern_vz,hess,z_filt_algo='gaussian',z_filt_pix=2.,sym='sym')
    
    fitswrite(updatedir('gradient_c_'+str(iterno).zfill(2)+'.fits'),totkern_c)
    
    if enf_cont and psi_cont:
        fitswrite(updatedir('gradient_psi_'+str(iterno).zfill(2)+'.fits'),totkern_psi)
    elif enf_cont and vx_cont:
        fitswrite(updatedir('gradient_vx_'+str(iterno).zfill(2)+'.fits'),totkern_vx)
    elif enf_cont and vz_cont:
        fitswrite(updatedir('gradient_vz_'+str(iterno).zfill(2)+'.fits'),totkern_vx)
    elif not enf_cont:
        fitswrite(updatedir('gradient_vz_'+str(iterno).zfill(2)+'.fits'),totkern_vz)
        fitswrite(updatedir('gradient_vx_'+str(iterno).zfill(2)+'.fits'),totkern_vx)

    #~ Get update direction based on algorithm of choice
    if (iterno==0) or steepest_descent:
        
        def sd_update(var='psi'):
            update = fitsread(updatedir('gradient_'+var+'_'+str(iterno).zfill(2)+'.fits'))
            updatemax = update.max()
            if updatemax!=0: update/=updatemax
            fitswrite(updatedir('update_'+var+'_'+str(iterno).zfill(2)+'.fits'),update)
        
        if os.path.exists(updatedir('model_c_'+str(iterno).zfill(2)+'.fits')):
            sd_update(var='c')
        
        if enf_cont and psi_cont: sd_update(var='psi')
                
        elif enf_cont and vx_cont: sd_update(var='vx')
                
        elif enf_cont and vz_cont: sd_update(var='vz')

        elif not enf_cont:
            sd_update(var='vz')
            sd_update(var='vx')
        
        if iterno > 0: print 'Forcing steepest descent'

    elif (conjugate_gradient):
        
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
            
            #~ beta=max(beta,0)
                
            print "beta",beta
            if beta==0: print "Conjugate gradient reduces to steepest descent"
                
            return beta
        
        def cg_update(var='psi'):
            grad=fitsread(updatedir('gradient_'+var+'_'+str(iterno).zfill(2)+'.fits'))
            lastgrad=fitsread(updatedir('gradient_'+var+'_'+str(iterno-1).zfill(2)+'.fits'))
            lastupdate=fitsread(updatedir('update_'+var+'_'+str(iterno-1).zfill(2)+'.fits'))
            
            beta = get_beta(grad,lastgrad,lastupdate)
            update=grad +  beta*lastupdate
            updatemax = update.max()
            if updatemax!=0: update/=updatemax
            fitswrite(updatedir('update_'+var+'_'+str(iterno).zfill(2)+'.fits'),update)
            
        
        if enf_cont and psi_cont: cg_update(var='psi')
            
        elif enf_cont and vx_cont: cg_update(var='vx')
            
        elif enf_cont and vz_cont: cg_update(var='vz')
            
        elif not enf_cont:
            cg_update(var='vz')
            cg_update(var='vx')
        
        try: cg_update(var='c')
        except IOError: pass
            
            
        
    elif LBFGS:
        #needs to be corrected
        m=4
        k = iterno
        
            
        def LBFGS_update(var='psi'):
            model_k=fitsread(updatedir('model_'+var+'_'+str(k).zfill(2)+'.fits'))
            grad_k =fitsread(updatedir('gradient_'+var+'_'+str(k).zfill(2)+'.fits'))
            model_iplus1 = model_k
            grad_iplus1 = grad_k
            q = grad_k
            
            for i in xrange(k-1,k-m-1,-1):
                
                if var+'_y_'+str(i) in LBFGS_data.keys():
                    y_i=LBFGS_data[var+'_y_'+str(i)]
                else:
                    grad_i=fitsread(updatedir('gradient_'+var+'_'+str(i).zfill(2)+'.fits'))
                    y_i = grad_iplus1 - grad_i
                    LBFGS_data[var+'_y_'+str(i)] = y_i
                    grad_iplus1 = grad_i
                    
                if var+'_s_'+str(i) in LBFGS_data.keys():
                    s_i=LBFGS_data[var+'_s_'+str(i)]
                else:
                    model_i=fitsread(updatedir('model_'+var+'_'+str(i).zfill(2)+'.fits'))
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
            
            
            updatemax = r.max()
            if updatemax!=0: r/=updatemax
            
            fitswrite(updatedir('update_'+var+'_'+str(iterno).zfill(2)+'.fits'),r)

        LBFGS_data_file = updatedir('LBFGS_data.npz')
        try:
            LBFGS_data={} 
            data=np.load(LBFGS_data_file)
            for key,val in data.items(): LBFGS_data[key]=val
            data=None
        except IOError: LBFGS_data={}
        
        if enf_cont and psi_cont: LBFGS_update(var='psi')
            
        elif enf_cont and vx_cont: LBFGS_update(var='vx')
            
        elif enf_cont and vz_cont: LBFGS_update(var='vz')
            
        elif not enf_cont:
            LBFGS_update(var='vx')
            LBFGS_update(var='vz')
            
        try: LBFGS_update(var='c')
        except IOError: pass
        
        np.savez(LBFGS_data_file,**LBFGS_data)
            
            
            
    #~ Create new models to be used for linesearch

    if enf_cont and psi_cont:
        model_psi = fitsread(updatedir('model_psi_'+str(iterno).zfill(2)+'.fits'))
        update_psi = fitsread(updatedir('update_psi_'+str(iterno).zfill(2)+'.fits'))
        psi_scale=rms(model_psi)
        
    elif enf_cont and vx_cont:
        model_vx=fitsread(updatedir('model_vx_'+str(iterno).zfill(2)+'.fits'))
        update_vx = fitsread(updatedir('update_vx_'+str(iterno).zfill(2)+'.fits'))
        vx_scale = 2000
        
    elif enf_cont and vz_cont:
        model_vz=fitsread(updatedir('model_vz_'+str(iterno).zfill(2)+'.fits'))
        update_vz = fitsread(updatedir('update_vz_'+str(iterno).zfill(2)+'.fits'))
        vz_scale=100
        
    elif not enf_cont:
        update_vx = fitsread(updatedir('update_vx_'+str(iterno).zfill(2)+'.fits'))
        vx_scale = 200
        
        update_vz = fitsread(updatedir('update_vz_'+str(iterno).zfill(2)+'.fits'))
        vz_scale = 100
    
    if os.path.exists(updatedir('model_c_'+str(iterno).zfill(2)+'.fits')):
        model_c = fitsread(updatedir('model_c_'+str(iterno).zfill(2)+'.fits'))
        update_c = fitsread(updatedir('update_c_'+str(iterno).zfill(2)+'.fits'))
        
    for i,eps_i in enumerate(eps):
        
        if os.path.exists(updatedir('model_c_'+str(iterno).zfill(2)+'.fits')):
            lsmodel = model_c + eps_i * update_c
            fitswrite(updatedir('test_c_'+str(i+1)+'.fits'), lsmodel)
            
        if enf_cont and psi_cont:
            lsmodel = model_psi + eps_i * psi_scale * update_psi
            fitswrite(updatedir('test_psi_'+str(i+1)+'.fits'), lsmodel)
            
        elif enf_cont and vx_cont:
            lsmodel = lastmodel_vx - vx_scale * eps_i * update_vx
            fitswrite(updatedir('test_vx_'+str(i+1)+'.fits'), lsmodel)

        elif enf_cont and vz_cont:
            lsmodel = model_vz - vz_scale * eps_i * update_vz
            fitswrite(updatedir('test_vz_'+str(i+1)+'.fits'), lsmodel)
            
        elif not enf_cont:
            lsmodel = model_vx - vx_scale*eps_i * update_vx
            fitswrite(updatedir('test_vx_'+str(i+1)+'.fits'), lsmodel)
            
            lsmodel = model_vz - vz_scale*eps_i * update_vz
            fitswrite(updatedir('test_vz_'+str(i+1)+'.fits'), lsmodel)


    try:
        epslist=np.load('epslist.npz')
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
    np.savez('epslist.npz',**epslist)


########################################################################

if __name__ == "__main__":
    main()
