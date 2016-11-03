from __future__ import division,print_function
import matplotlib
from matplotlib import pyplot as plt
import numpy as np
import read_params
import os,fnmatch

datadir = read_params.get_directory()

coeff_files = sorted(fnmatch.filter(os.listdir(os.path.join(datadir,"update")),
                "model_psi_[0-9][0-9]_coeffs.npz"))

spline_basis_coeffs = [np.load(os.path.join(datadir,"update",f)) for f in coeff_files]
spline_basis_z = [s["z"] for s in spline_basis_coeffs]

true_coeffs = np.load(os.path.join(datadir,"true_psi_coeffs.npz"))
true_coeffs_z = true_coeffs["cz_bot"]

plt.bar(np.arange(true_coeffs_z.size)-0.4,true_coeffs_z,width=0.8,
        facecolor="papayawhip",edgecolor="peru",label="true")

colors = ["chocolate","olive","cadetblue","seagreen","darkgoldenrod","slateblue"]

for iterno,sz in enumerate(spline_basis_z):
    plt.plot(sz,ls="-",marker="o",markersize=3,label="iter "+str(iterno),
    color=colors[iterno])

plt.legend(loc="best")
plt.show()
