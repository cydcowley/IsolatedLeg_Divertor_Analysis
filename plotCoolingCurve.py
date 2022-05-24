from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np

from AnalyticCoolingCurves import LfuncN
Exhemptfiles = ["baserun","notes.txt","fi75E-3","fi280E-3norecomb","fi280E-3"]

plt.rcParams["font.family"] = "serif"
params = {'legend.fontsize': 'medium',
         'axes.labelsize': 'medium',
         'axes.titlesize':'medium',
         'xtick.labelsize':'medium',
         'ytick.labelsize':'medium',
        #  'figure.figsize': (4,3.2),
         }
plt.rcParams.update(params)
T = np.linspace(1,80,10000)
L = []
for t in T:
    L.append(LfuncN(t))
plt.plot(T,L,linewidth=2)
plt.xlabel("T (eV)")
ylabel = "L"+r"$_{N}$"+" (Wm"+r"$^{3}$"+")"
plt.ylabel(ylabel)
plt.savefig("Figures/coolingCurve.png",dpi=400,bbox_inches='tight')
plt.show()


