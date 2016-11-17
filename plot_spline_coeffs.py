from __future__ import division,print_function
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



colors = itertools.cycle(["mediumseagreen","darkgoldenrod","midnightblue","chocolate"])

for iterno in xrange(0,len(spline_basis_z)-1,max(1,len(spline_basis_z)//3)):
    plt.plot(spline_basis_z[iterno],ls="dotted",marker="o",
    markersize=3,label="iter "+str(iterno),
    color=next(colors),zorder=1)

iterno = len(spline_basis_z)-1
plt.plot(spline_basis_z[iterno],ls="solid",marker="o",markersize=3,label="iter "+str(iterno),
color=next(colors),lw=1.5,zorder=2)

plt.gca().margins(y=0.1)

plt.bar(np.arange(true_coeffs_z.size)-0.4,true_coeffs_z,width=0.8,
        facecolor="papayawhip",edgecolor="peru",label="true",zorder=0)

plt.legend(loc="best")
plt.show()
