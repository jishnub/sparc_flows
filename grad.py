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
    optimization_algo=filter(lambda x: x.startswith('algo='),args)
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

    #~ Get shape
    nx = read_params.get_nx()
    ny = 1
    nz = read_params.get_nz()
    
    Lx = read_params.get_xlength()

    back=np.loadtxt(read_params.get_solarmodel())

    num_src=get_number_of_sources()

    array_shape=(nx,ny,nz)
    totkern_c=np.zeros(array_shape)
    totkern_psi=np.zeros(array_shape)
    totkern_vx=np.zeros(array_shape)
    totkern_vz=np.zeros(array_shape)
    hess=np.zeros(array_shape)
    
    sound_speed_perturbed = read_params.if_soundspeed_perturbed()
    flows = read_params.if_flows()
    continuity_enforced=read_params.if_continuity_enforced()
    cont_var=read_params.get_continuity_variable()
    psi_cont = False
    vx_cont = False
    vz_cont = False
    if continuity_enforced:
        if cont_var == 'psi': psi_cont = True 
        elif cont_var == 'vx': vx_cont = True 
        elif cont_var == 'vz': vz_cont = True 
        else: 
            print "Continuity variable unknown, check params.i"
            quit()

    def read_model(var='psi',iterno=iterno):
        return fitsread(updatedir('model_'+var+'_'+str(iterno).zfill(2)+'.fits'))
        
    def read_grad(var='psi',iterno=iterno):
        return fitsread(updatedir('gradient_'+var+'_'+str(iterno).zfill(2)+'.fits'))
    
    def read_update(var='psi',iterno=iterno):
        return fitsread(updatedir('update_'+var+'_'+str(iterno).zfill(2)+'.fits'))
        
    def read_kern(var='psi',src=1):
        return fitsread(os.path.join(datadir,'kernel','kernel_'+var+'_'+str(src).zfill(2)+'.fits'))

    for src in xrange(1,num_src+1):
        
        if sound_speed_perturbed:
            totkern_c += read_kern(var='c',src=src)
        
        if flows:
            if continuity_enforced and psi_cont:
                totkern_psi += read_kern(var='psi',src=src)
            
            elif continuity_enforced and vx_cont:
                totkern_vx += read_kern(var='vx',src=src)
                
            elif continuity_enforced and vz_cont:
                totkern_vz += read_kern(var='vz',src=src)
                
            elif not continuity_enforced:
                totkern_vx += read_kern(var='vx',src=src)
                totkern_vz += read_kern(var='vz',src=src)
            
        hess+=abs(fitsread(os.path.join(datadir,'kernel','hessian_'+str(src).zfill(2)+'.fits')))

    hess=hess*np.atleast_3d(back[:,2]).transpose(0,2,1)
    hess = hess/abs(hess).max()
    hess[hess<5e-3]=5e-3

    #~ Smoothing and symmetrization 
    if sound_speed_perturbed: 
        totkern_c = filter_and_symmetrize(totkern_c,hess,z_filt_algo='gaussian',z_filt_pix=5,sym='sym')
    if flows:
        if continuity_enforced and psi_cont:
            totkern_psi = filter_and_symmetrize(totkern_psi,hess,z_filt_algo='gaussian',z_filt_pix=2.,sym='asym')
        elif (continuity_enforced and vx_cont) or (not continuity_enforced):
            totkern_vx = filter_and_symmetrize(totkern_vx,hess,z_filt_algo='gaussian',z_filt_pix=2.,sym='asym')
        elif (continuity_enforced and vz_cont) or (not continuity_enforced):
            totkern_vz = filter_and_symmetrize(totkern_vz,hess,z_filt_algo='gaussian',z_filt_pix=2.,sym='sym')
    
    #~ Write out gradients for this iteration
    if sound_speed_perturbed:
        fitswrite(updatedir('gradient_c_'+str(iterno).zfill(2)+'.fits'),-totkern_c)
    if flows:
        if continuity_enforced and psi_cont:
            fitswrite(updatedir('gradient_psi_'+str(iterno).zfill(2)+'.fits'),totkern_psi)
        elif continuity_enforced and vx_cont:
            fitswrite(updatedir('gradient_vx_'+str(iterno).zfill(2)+'.fits'),totkern_vx)
        elif continuity_enforced and vz_cont:
            fitswrite(updatedir('gradient_vz_'+str(iterno).zfill(2)+'.fits'),totkern_vx)
        elif not continuity_enforced:
            fitswrite(updatedir('gradient_vz_'+str(iterno).zfill(2)+'.fits'),totkern_vz)
            fitswrite(updatedir('gradient_vx_'+str(iterno).zfill(2)+'.fits'),totkern_vx)

    #~ Get update direction based on algorithm of choice
    if (iterno==0) or steepest_descent:
        
        def sd_update(var='psi'):
            grad = read_grad(var=var,iterno=iterno)
            update = -grad
            fitswrite(updatedir('update_'+var+'_'+str(iterno).zfill(2)+'.fits'),update)
            print "Steepest descent"
        
        if sound_speed_perturbed: sd_update(var='c')
        
        if flows:
            if continuity_enforced and psi_cont: sd_update(var='psi')
                    
            elif continuity_enforced and vx_cont: sd_update(var='vx')
                    
            elif continuity_enforced and vz_cont: sd_update(var='vz')

            elif not continuity_enforced:
                sd_update(var='vz')
                sd_update(var='vx')
            
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
            
            fitswrite(updatedir('update_'+var+'_'+str(iterno).zfill(2)+'.fits'),update)
        
        if sound_speed_perturbed: cg_update(var='c')    
        
        if flows:
            if continuity_enforced and psi_cont: cg_update(var='psi')
                
            elif continuity_enforced and vx_cont: cg_update(var='vx')
                
            elif continuity_enforced and vz_cont: cg_update(var='vz')
                
            elif not continuity_enforced:
                cg_update(var='vz')
                cg_update(var='vx')
        
            
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

        if sound_speed_perturbed: LBFGS_update(var='c')
        
        if flows:
            if continuity_enforced and psi_cont: LBFGS_update(var='psi')
                
            elif continuity_enforced and vx_cont: LBFGS_update(var='vx')
                
            elif continuity_enforced and vz_cont: LBFGS_update(var='vz')
                
            elif not continuity_enforced:
                LBFGS_update(var='vx')
                LBFGS_update(var='vz')
        
        np.savez(LBFGS_data_file,**LBFGS_data)
            

    def create_ls_model(var='psi',eps=0,kind='linear'):
        model = read_model(var=var,iterno=iterno)
        update = read_update(var=var,iterno=iterno)
        
        updatemax=update.max()
        if updatemax!=0: update/=updatemax
        
        cutoff_switch,cutoff_dist=read_params.get_cutoff_dist()
        cutoff = 1.0
        
        if cutoff_switch:
            pix = np.arange(nx)
            xcutoffpix = cutoff_dist/Lx*nx
            cutoff =  1./(1+np.exp((pix-(nx/2+xcutoffpix))/2.))+1./(1+np.exp(-(pix-(nx/2-xcutoffpix))/2.))-1.
            cutoff = cutoff.reshape(nx,1,1)
        
        if kind=='linear':
            model_scale = rms(model)
            if model_scale == 0: model_scale = 100
            test= model+eps*update*model_scale*cutoff
        elif kind=='exp':
            test= model*(1+eps*update*cutoff)
        
        return test
    
    #~ Create models for linesearch
    for i,eps_i in enumerate(eps):
        
        if sound_speed_perturbed:
            lsmodel = create_ls_model(var='c',eps=eps_i,kind='exp')
            fitswrite(updatedir('test_c_'+str(i+1)+'.fits'), lsmodel)
        
        if flows:    
            if continuity_enforced and psi_cont:
                lsmodel = create_ls_model(var='psi',eps=eps_i,kind='linear')
                fitswrite(updatedir('test_psi_'+str(i+1)+'.fits'), lsmodel)
                
            elif continuity_enforced and vx_cont:
                lsmodel = create_ls_model(var='vx',eps=eps_i,kind='linear')
                fitswrite(updatedir('test_vx_'+str(i+1)+'.fits'), lsmodel)

            elif continuity_enforced and vz_cont:
                lsmodel = create_ls_model(var='vz',eps=eps_i,kind='linear')
                fitswrite(updatedir('test_vz_'+str(i+1)+'.fits'), lsmodel)
                
            elif not continuity_enforced:
                lsmodel = create_ls_model(var='vx',eps=eps_i,kind='linear')
                fitswrite(updatedir('test_vx_'+str(i+1)+'.fits'), lsmodel)
                
                lsmodel = create_ls_model(var='vz',eps=eps_i,kind='linear')
                fitswrite(updatedir('test_vz_'+str(i+1)+'.fits'), lsmodel)


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
