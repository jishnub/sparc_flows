
import matplotlib
from matplotlib import pyplot as plt
import numpy as np
import read_params
import os,fnmatch,itertools

datadir = read_params.get_directory()

coeff_files = sorted(fnmatch.filter(os.listdir(os.path.join(datadir,"update")),
                "model_psi_[0-9][0-9]_coeffs.npz"))

spline_basis_coeffs = [np.load(os.path.join(datadir,"update",f)) for f in coeff_files]
spline_basis_z = [s["z"] for s in spline_basis_coeffs]

true_coeffs = np.load(os.path.join(datadir,"true_psi_coeffs.npz"))
true_coeffs_z = true_coeffs["cz_bot"]
true_coeffs_top = true_coeffs["cz_top"]
c_surf_cutoff = true_coeffs['c_surf_cutoff']

modenames = {'0':'f'}
for i in range(1,8): modenames[str(i)] = 'p{:d}'.format(i)

iterno = read_params.parse_cmd_line_params(key="iter",
        default=len(spline_basis_z)-1,mapto=int)


plt.plot(spline_basis_z[iterno][:c_surf_cutoff],ls="solid",
marker="o",markersize=6,label="Iterated\ncoefficients",
color="saddlebrown",lw=1.5,zorder=4)

bridge = [spline_basis_z[iterno][c_surf_cutoff-1],true_coeffs_top[c_surf_cutoff]]

plt.plot(np.arange(c_surf_cutoff-1,c_surf_cutoff+1),
bridge,ls="dashed",color='darkmagenta',lw=1,zorder=3)

plt.plot(np.arange(c_surf_cutoff,true_coeffs_top.size),
true_coeffs_top[c_surf_cutoff:],ls="dashed",marker='o',ms=8,color='darkmagenta',
mec='darkmagenta',mfc="lightcoral",lw=1,zorder=3,
label="Coefficients\nabove surface,\nclamped")

plt.gca().margins(y=0.1)

plt.bar(np.arange(true_coeffs_z.size)-0.4,true_coeffs_z+true_coeffs_top,width=0.8,
        facecolor="burlywood",edgecolor="peru",label="Coefficients\nfor reference\nmodel"
        ,zorder=1)

plt.xlabel("Coefficient number",fontsize=20)
plt.ylabel("Coefficient value",fontsize=20)

plt.tick_params(labelsize=13)

plt.legend(loc="best")

if not os.path.exists("plots"): os.makedirs("plots")

plt.savefig("plots/iterated_coefficients.eps")

plt.show()
