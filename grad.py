from __future__ import division
import numpy as np
from scipy.interpolate import UnivariateSpline
from scipy.ndimage.filters import gaussian_filter1d
import os,fnmatch
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
    return len(fnmatch.filter(os.listdir(updatedir),'misfit_[0-9][0-9]'))
    
def check_BFGS_stencil(iterno,stencilBFGS):
    try: assert (iterno - 1 - stencilBFGS)>1
    except AssertionError:    
        print "BFGS stencil too wide, quitting."
        quit()

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

def twodigit(m): return str(m).zfill(2)
    
def iterminus(m): return twodigit(iterno-m)

def updatedir(filename): return os.path.join(datadir,'update',filename)

def rms(arr): return np.sqrt(np.sum(arr**2)/np.prod(arr.shape)) 

########################################################################

eps = 1e-7

codedir=os.path.dirname(os.path.abspath(__file__))
datadir=read_params.get_directory()
iterno=get_iter_no()

def main(eps):


    Rsun=695.9895 # Mm

    #~ Decide on algorithm
    steepest_descent = False
    conjugate_gradient=True and not steepest_descent
    LBFGS = True and not (steepest_descent or conjugate_gradient)

    #~ Arbitrarily chosen maximum number of iterations. Hope we get to 5
    itermax=25

    #~ Get shape from a pre-existing file
    kern=fitsread(os.path.join(datadir,'kernel','kernel_c_01.fits'))
    nx,ny,nz=kern.shape

    

    back=np.loadtxt('polytrope')

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
        
        kern=fitsread(os.path.join(datadir,'kernel','kernel_c_'+twodigit(src)+'.fits'))
        totkern_c+=kern
        
        if enf_cont and psi_cont:
            kern=fitsread(os.path.join(datadir,'kernel','kernel_psi_'+twodigit(src)+'.fits'))
            totkern_psi+=kern
        
        elif enf_cont and vx_cont:
            kern=fitsread(os.path.join(datadir,'kernel','kernel_vx_'+twodigit(src)+'.fits'))
            totkern_vx+=kern
            
        elif enf_cont and vz_cont:
            kern=fitsread(os.path.join(datadir,'kernel','kernel_vz_'+twodigit(src)+'.fits'))
            totkern_vz+=kern
        else:
            kern=fitsread(os.path.join(datadir,'kernel','kernel_vx_'+twodigit(src)+'.fits'))
            totkern_vx+=kern
            kern=fitsread(os.path.join(datadir,'kernel','kernel_vz_'+twodigit(src)+'.fits'))
            totkern_vz+=kern
            
        
        kern=fitsread(os.path.join(datadir,'kernel','hessian_'+twodigit(src)+'.fits'))
        hess+=abs(kern)

    hess=hess*np.atleast_3d(back[:,2]).transpose(0,2,1)
    hess = hess/abs(hess).max()
    hess[hess<5e-3]=5e-3


    #~ Sound speed kernel
    kern = totkern_c/hess
    filterx(kern)
    filterz(kern,algo='smooth',sp=0.3)
    totkern_c=kern

    #~ Vector potential/stream function kernel
    kern = totkern_psi/hess
    filterx(kern)
    antisymmetrize(kern)
    filterz(kern,algo='gaussian',sp=1.0)
    totkern_psi=kern
    
    #~ Velocity kernels
    kern = totkern_vx/hess
    #~ filterx(kern)
    #~ antisymmetrize(kern)
    #~ filterz(kern,algo='gaussian',sp=2.0)
    totkern_vx = kern
    
    kern = totkern_vz/hess
    filterx(kern)
    symmetrize(kern)
    filterz(kern,algo='gaussian',sp=2.0)
    totkern_vz = kern
    
    fitswrite(updatedir('gradient_c_'+iterminus(1)+'.fits'),totkern_c)
    
    if enf_cont and psi_cont:
        fitswrite(updatedir('gradient_psi_'+iterminus(1)+'.fits'),totkern_psi)
    elif enf_cont and vx_cont:
        fitswrite(updatedir('gradient_vx_'+iterminus(1)+'.fits'),totkern_vx)
    elif enf_cont and vz_cont:
        fitswrite(updatedir('gradient_vz_'+iterminus(1)+'.fits'),totkern_vx)
    elif not enf_cont:
        fitswrite(updatedir('gradient_vz_'+iterminus(1)+'.fits'),totkern_vz)
        fitswrite(updatedir('gradient_vx_'+iterminus(1)+'.fits'),totkern_vx)

    model_c_exists=True
    model_psi_exists=True
    model_vx_exists=True
    model_vz_exists=True

    #~ Get update direction based on algorithm of choice
    if (iterno==1) or steepest_descent:
        try:
            lastmodel_c=fitsread(updatedir('model_c_'+iterminus(1)+'.fits'))
        except IOError: model_c_exists=False
            
        try:
            lastmodel_psi=fitsread(updatedir('model_psi_'+iterminus(1)+'.fits'))
        except IOError: model_psi_exists=False
        
        try:
            lastmodel_vx=fitsread(updatedir('model_vx_'+iterminus(1)+'.fits'))
        except IOError: model_vx_exists=False
        
        try:
            lastmodel_vz=fitsread(updatedir('model_vz_'+iterminus(1)+'.fits'))
        except IOError: model_vz_exists=False
        
        if model_c_exists:
            fitswrite(updatedir('update_c_'+iterminus(1)+'.fits'),totkern_c)
        
        if enf_cont and psi_cont:
            fitswrite(updatedir('update_psi_'+iterminus(1)+'.fits'),totkern_psi)
        elif enf_cont and vx_cont:
            fitswrite(updatedir('update_vx_'+iterminus(1)+'.fits'),totkern_vx)
        elif enf_cont and vz_cont:
            fitswrite(updatedir('update_vz_'+iterminus(1)+'.fits'),totkern_vz)
        
        elif not enf_cont:
            fitswrite(updatedir('update_vz_'+iterminus(1)+'.fits'),totkern_vz)
            fitswrite(updatedir('update_vx_'+iterminus(1)+'.fits'),totkern_vx)
        
        update_c = totkern_c 
        update_psi = totkern_psi
        update_vx = totkern_vx
        update_vz = totkern_vz
        
        if (iterno > 1): print 'Forcing steepest descent'

    elif (iterno>1 and conjugate_gradient):
        print 'Conjugate Gradient'
        try:
            lastmodel_c=fitsread(updatedir('model_c_'+iterminus(1)+'.fits'))
        except IOError: model_c_exists=False
        
        try:
            lastmodel_psi=fitsread(updatedir('model_psi_'+iterminus(1)+'.fits'))
        except IOError: model_psi_exists=False
        
        try:
            lastmodel_vx=fitsread(updatedir('model_vx_'+iterminus(1)+'.fits'))
        except IOError: model_vx_exists=False
        
        try:
            lastmodel_vz=fitsread(updatedir('model_vz_'+iterminus(1)+'.fits'))
        except IOError: model_vz_exists=False
        
        if enf_cont and psi_cont:
            grad_psi=fitsread(updatedir('gradient_psi_'+iterminus(1)+'.fits'))
            lastgrad_psi=fitsread(updatedir('gradient_psi_'+iterminus(2)+'.fits'))
            lastupdate_psi=fitsread(updatedir('update_psi_'+iterminus(2)+'.fits'))
        
            con = np.sum(grad_psi*(grad_psi - lastgrad_psi))
            den = np.sum(lastgrad_psi**2.)
            
        elif enf_cont and vx_cont:
            grad_vx=fitsread(updatedir('gradient_vx_'+iterminus(1)+'.fits'))
            lastgrad_vx=fitsread(updatedir('gradient_vx_'+iterminus(2)+'.fits'))
            lastupdate_vx=fitsread(updatedir('update_vx_'+iterminus(2)+'.fits'))
            
            con = np.sum(grad_vx*(grad_vx - lastgrad_vx))
            den = np.sum(lastgrad_vx**2.)
            
        elif enf_cont and vz_cont:
            grad_vz=fitsread(updatedir('gradient_vz_'+iterminus(1)+'.fits'))
            lastgrad_vz=fitsread(updatedir('gradient_vz_'+iterminus(2)+'.fits'))
            lastupdate_vz=fitsread(updatedir('update_vz_'+iterminus(2)+'.fits'))
            
            con = np.sum(grad_vz*(grad_vz - lastgrad_vz))
            den = np.sum(lastgrad_vz**2.)
        
        elif not enf_cont:
            grad_vz=fitsread(updatedir('gradient_vz_'+iterminus(1)+'.fits'))
            lastgrad_vz=fitsread(updatedir('gradient_vz_'+iterminus(2)+'.fits'))
            lastupdate_vz=fitsread(updatedir('update_vz_'+iterminus(2)+'.fits'))
            
            grad_vx=fitsread(updatedir('gradient_vx_'+iterminus(1)+'.fits'))
            lastgrad_vx=fitsread(updatedir('gradient_vx_'+iterminus(2)+'.fits'))
            lastupdate_vx=fitsread(updatedir('update_vx_'+iterminus(2)+'.fits'))
            
            con = np.sum(grad_vx*(grad_vx - lastgrad_vx))
            den = np.sum(lastgrad_vx**2.)
        
        if model_c_exists:
            grad_c=fitsread(updatedir('gradient_c_'+iterminus(1)+'.fits'))
            lastgrad_c=fitsread(updatedir('gradient_c_'+iterminus(2)+'.fits'))
            lastupdate_c=fitsread(updatedir('update_c_'+iterminus(2)+'.fits'))
        
        
        con=con/den
        den=con
        if con<0: con=0
        print con,den
        
        if model_c_exists:
            update_c = totkern_c +  con* lastupdate_c
            fitswrite(updatedir('update_c_'+iterminus(1)+'.fits'),update_c)
            
        if enf_cont and psi_cont:
            update_psi = totkern_psi +  con* lastupdate_psi
            fitswrite(updatedir('update_psi_'+iterminus(1)+'.fits'),update_psi)
        elif enf_cont and vx_cont:
            update_vx = totkern_vx +  con* lastupdate_vx
            fitswrite(updatedir('update_vx_'+iterminus(1)+'.fits'),update_vx)
        elif enf_cont and vz_cont:
            update_vz = totkern_vz +  con* lastupdate_vz
            fitswrite(updatedir('update_vz_'+iterminus(1)+'.fits'),update_vz)
        elif not enf_cont:
            update_vx = totkern_vx +  con* lastupdate_vx
            fitswrite(updatedir('update_vx_'+iterminus(1)+'.fits'),update_vx)
            
            update_vz = totkern_vz +  con* lastupdate_vz
            fitswrite(updatedir('update_vz_'+iterminus(1)+'.fits'),update_vz)
        
        
    elif (iterno>1 and LBFGS):
        #needs to be corrected
        stencilBFGS=4
        check_BFGS_stencil(iterno,stencilBFGS)
        minBFGS=iterno-1-stencilBFGS
        
        alph=np.zeros(iteration)
        
        try:
            model_c=fitsread(updatedir('model_c_'+iterminus(1)+'.fits'))
        except: model_c_exists=False
        
        grad_c=fitsread(updatedir('gradient_c_'+iterminus(1)+'.fits'))
        update_c=grad_c
        
        model_psi=fitsread(updatedir('model_psi_'+iterminus(1)+'.fits'))
        grad_psi=fitsread(updatedir('gradient_psi_'+iterminus(1)+'.fits'))
        update_psi=grad_psi
     ## for c  
        if model_c_exists:
            for i in xrange(iterno-2,minBFGS-1,-1):
                
                previterind=twodigit(i)
                try:
                    lastmodel=fitsread(updatedir('model_c_'+previterind+'.fits'))
                except IOError: model_c_exists=False
                
                lastgrad_c=fitsread(updatedir('gradient_c_'+previterind+'.fits'))
                
                alph[i] = np.sum((model-lastmodel)*update_c) /np.sum((model-lastmodel)*(grad_c-lastgrad_c))
                update_c = update_c - alph[i] * (grad_c - lastgrad_c) 
                
                model = lastmodel
                grad_c = lastgrad_c
                
        lastmodel=fitsread(updatedir('model_c_'+twodigit(iterno-2)+'.fits'))
        lastgrad_psi=fitsread(updatedir('gradient_c_'+twodigit(iterno-2)+'.fits'))   
        model=fitsread(updatedir('model_c_'+twodigit(iterno-1)+'.fits'))
        grad_psi=fitsread(updatedir('gradient_c_'+twodigit(iterno-1)+'.fits'))     
        
        hessinv_c=((grad_c-lastgrad_c)*(model-lastmodel))/((grad_c-lastgrad_c)*(grad_c-lastgrad_c))
        update_c = update_c * hessinv_c
         
        for i in xrange(minBFGS,iterno-1):
                model=fitsread(updatedir('model_c_'+twodigit(i)+'.fits'))
                grad_c=fitsread(updatedir('gradient_c_'+twodigit(i)+'.fits'))
                
                alph[i] = alph[i] - np.sum((grad_c-lastgrad_c)*update_c) /np.sum((model-lastmodel)*(grad_c-lastgrad_c))
                update_c = update_c + alph[i] * (model - lastmodel) 


                lastmodel = model
                lastgrad_c = grad_c
     ## for psi
        for i in xrange(iterno-2,minBFGS-1,-1):
        
            lastmodel=fitsread(updatedir('model_psi_'+previterind+'.fits'))
            lastgrad_psi=fitsread(updatedir('gradient_psi_'+previterind+'.fits'))
            
            alph[i] = np.sum((model-lastmodel)*update_psi) /np.sum((model-lastmodel)*(grad_psi-lastgrad_psi))    
            update_psi = update_psi - alph[i] * (model - lastmodel)  
            
            model = lastmodel
            grad_psi = lastgrad_psi 
                
                
        lastmodel=fitsread(updatedir('model_psi_'+twodigit(iterno-2)+'.fits'))
        lastgrad_psi=fitsread(updatedir('gradient_psi_'+twodigit(iterno-2)+'.fits'))   
        model=fitsread(updatedir('model_psi_'+twodigit(iterno-1)+'.fits'))
        grad_psi=fitsread(updatedir('gradient_psi_'+twodigit(iterno-1)+'.fits'))     
        
        hessinv_psi=((grad_psi-lastgrad_psi)*(model-lastmodel))/((grad_psi-lastgrad_psi)*(grad_psi-lastgrad_psi))
        update_psi = update_psi * hessinv_psi
        
            
        #~ for i in xrange(minBFGS,iterno-1): 
    #~ Create new models to be used for linesearch

    if enf_cont and psi_cont:
        psimax = abs(update_psi).max()
        update_psi = update_psi/psimax 
        if model_c_exists:
            update_c = update_c/psimax
        psi_scale=rms(lastmodel_psi)
        
    elif enf_cont and vx_cont:
        vx_max = abs(update_vx).max()
        update_vx/=vx_max
        if model_c_exists:
            update_c = update_c/vx_max
        vx_scale = 2000
        
    elif enf_cont and vz_cont:
        vz_max = abs(update_vz).max()
        update_vz/=vz_max
        if model_c_exists:
            update_c = update_c/vz_max
        vz_scale = 1e-8
        
    elif not enf_cont:
        vx_max = abs(update_vx).max()
        update_vx/=vx_max
        if model_c_exists:
            update_c = update_c/vx_max
        vx_scale = 200
        
        update_vz/=vx_max
        vz_scale = 100
        if model_c_exists:
            update_c = update_c/vx_max
    
    for i in xrange(1,6):
        if model_c_exists:
            update = lastmodel_c
            fitswrite(updatedir('test_c_'+str(i)+'.fits'), update)
        if enf_cont and psi_cont:
            update = lastmodel_psi + eps * i * psi_scale * update_psi
            fitswrite(updatedir('test_psi_'+str(i)+'.fits'), update)
        elif enf_cont and vx_cont:
            if model_vx_exists:
                update = lastmodel_vx - vx_scale * eps * i * update_vx
                fitswrite(updatedir('test_vx_'+str(i)+'.fits'), update)
            else:
                print "model vx doesn't exist"
        elif enf_cont and vz_cont:
            if model_vz_exists:
                update = lastmodel_vz - vz_scale * eps * i * update_vz
                fitswrite(updatedir('test_vz_'+str(i)+'.fits'), update)
            else:
                print "model vz doesn't exist"
        elif not enf_cont:
            if model_vx_exists:
                
                update = lastmodel_vx - vx_scale*eps * i * update_vx
                fitswrite(updatedir('test_vx_'+str(i)+'.fits'), update)
            else: print "model vx doesn't exist"
            if model_vz_exists:
                update = lastmodel_vz - vz_scale*eps * i * update_vz
                fitswrite(updatedir('test_vz_'+str(i)+'.fits'), update)
            else: print "model vz doesn't exist"


    try:
        epslist=np.load('epslist.npz')['epslist']
        iterind=np.where(epslist[:,0]==iterno)[0]
        if len(iterind)==0:
            np.append(epslist,[[iterno,eps]],axis=0)
        else:
            epslist[iterind,1]=eps
    except IOError:
        epslist=[[iterno,eps]]

    np.savez('epslist.npz',epslist=epslist)


########################################################################

if __name__ == "__main__":
    main(eps)
