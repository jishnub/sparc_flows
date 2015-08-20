import sys,os,shutil

if len(sys.argv)==1:
    print "Usage: python copy_updated_model.py <model number>"
    print "Example: python copy_updated_model.py 3"
    quit()

num=sys.argv[1]

codedir=os.path.dirname(os.path.abspath(__file__))
configvars={}
with open(os.path.join(codedir,"varlist.sh")) as myfile:
    for line in myfile:
        name,var=line.partition("=")[::2]
        configvars[name.strip()]=var.strip().strip('"')

datadir=configvars['directory']

updated_model_c=os.path.join(datadir,'update','test_c_'+num+'.fits')
updated_model_vectorpsi=os.path.join(datadir,'update','test_psi_'+num+'.fits')

model_c=os.path.join(datadir,'model_c_ls00.fits')
model_vectorpsi=os.path.join(datadir,'model_psi_ls00.fits')

#~ shutil.copyfile(updated_model_c,model_c)
shutil.copyfile(updated_model_vectorpsi,model_vectorpsi)


