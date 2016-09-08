import sys,os,shutil,read_params

num = read_params.parse_cmd_line_params(key="num")
basis = read_params.parse_cmd_line_params(key="basis")

datadir=read_params.get_directory()

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
