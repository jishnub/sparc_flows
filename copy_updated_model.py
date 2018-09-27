import sys,os,shutil,read_params,fnmatch,numpy as np

datadir=read_params.get_directory()

lsfiles = sorted(fnmatch.filter(os.listdir(os.path.join(datadir,"update")),
                "linesearch_[0-9][0-9]"))
ls_last = np.loadtxt(os.path.join(datadir,"update",lsfiles[-1]))

ls_last = np.array([sum(ls_last[i*8:(i+1)*8,2])
        for i in range(ls_last.shape[0]//8)])

np.set_printoptions(precision=1)
print("Linesearch misfits",ls_last)
num = read_params.parse_cmd_line_params(key="num",
        default=ls_last.argmin()+1,mapto=str)
if os.path.exists(os.path.join(datadir,"model_psi_ls00_coeffs.npz")):
    basis="spline"
basis = read_params.parse_cmd_line_params(key="basis",default=basis)

print(("Copying model {}, basis {}".format(num,basis)))
updated_model_c=os.path.join(datadir,'update','test_c_'+num+'.fits')
updated_model_vectorpsi=os.path.join(datadir,'update','test_psi_'+num+'.fits')
updated_model_vectorpsi_coeff = os.path.join(datadir,'update','test_psi_'+num+'_coeffs.npz')

model_c=os.path.join(datadir,'model_c_ls00.fits')
model_vectorpsi=os.path.join(datadir,'model_psi_ls00.fits')
model_vectorpsi_coeffs=os.path.join(datadir,'model_psi_ls00_coeffs.npz')

#~ shutil.copyfile(updated_model_c,model_c)
shutil.copyfile(updated_model_vectorpsi,model_vectorpsi)
if basis=="spline":
    shutil.copyfile(updated_model_vectorpsi_coeff,model_vectorpsi_coeffs)
