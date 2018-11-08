import os,shutil,read_params,fnmatch,numpy as np
from pathlib import Path

datadir= Path(read_params.get_directory())

lsfiles = sorted(fnmatch.filter(os.listdir(datadir/"update"),"linesearch_[0-9][0-9]"))
ls_last = np.loadtxt(datadir/"update"/lsfiles[-1])

ls_last = np.array([sum(ls_last[i*8:(i+1)*8,2])
        for i in range(ls_last.shape[0]//8)])

print("Linesearch misfits: "+("{:.1e} "*len(ls_last)).format(*ls_last))
num = read_params.parse_cmd_line_params(key="num",default=ls_last.argmin()+1,mapto=int)

for var in ["c","psi"]:

	basis = None
	if (datadir/"model_{}_ls00_coeffs.npz".format(var)).exists():  basis="spline"

	basis = read_params.parse_cmd_line_params(key="basis",default=basis)

	

	test_model = datadir/'update'/'test_{}_{:d}.fits'.format(var,num)
	if test_model.exists():
		print(("Copying {} model {}, basis {}".format(var,num,basis)))
		shutil.copyfile(test_model,datadir/'model_{}_ls00.fits'.format(var))

	if basis=="spline":
	    shutil.copyfile(datadir/'update'/'test_{}_{:d}_coeffs.npz'.format(var,num),
	    				datadir/'model_{}_ls00_coeffs.npz'.format(var))
