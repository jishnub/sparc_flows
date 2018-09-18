import sys,os,shutil,read_params
   
try:
    num=sys.argv[1]
    int(num)
except:
    print("Usage: python copy_updated_model.py <model number>")
    print("Example: python copy_updated_model.py 3")
    quit()

datadir=read_params.get_directory()

updated_model_c=os.path.join(datadir,'update','test_c_'+num+'.fits')
updated_model_vectorpsi=os.path.join(datadir,'update','test_psi_'+num+'.fits')

model_c=os.path.join(datadir,'model_c_ls00.fits')
model_vectorpsi=os.path.join(datadir,'model_psi_ls00.fits')

#~ shutil.copyfile(updated_model_c,model_c)
shutil.copyfile(updated_model_vectorpsi,model_vectorpsi)


