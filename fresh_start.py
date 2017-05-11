import os,shutil,glob
import read_params

datadir = read_params.get_directory()

for uf in glob.glob(os.path.join(datadir,"update","*")):
    try:
        os.remove(uf)
    except:
        shutil.rmtree(uf)

for uf in glob.glob(os.path.join(datadir,"kernel","*")):
    try:
        os.remove(uf)
    except:
        shutil.rmtree(uf)

for uf in glob.glob(os.path.join(datadir,"tt","*")):
    try:
        os.remove(uf)
    except:
        shutil.rmtree(uf)

os.remove(os.path.join(datadir,"model_psi_ls00.png"))
os.remove(os.path.join(datadir,"model_psi_ls00.fits"))
os.remove(os.path.join(datadir,"model_psi_ls00_coeffs.npz"))

shutil.copyfile(os.path.join(datadir,"model_psi_ls00_coeffs_start.npz"),
os.path.join(datadir,"model_psi_ls00_coeffs.npz"))

shutil.copyfile(os.path.join(datadir,"model_psi_ls00_start.fits"),
os.path.join(datadir,"model_psi_ls00.fits"))
