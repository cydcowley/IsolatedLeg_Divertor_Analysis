from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import os
import matplotlib.image as m
import imageio as io
from scipy import interpolate
from PIL import Image
import cv2
import glob
from scipy.integrate import trapz
import mpltools
from mpltools import special
import sys
sys.path.append('d:\\my stuff\\PhD\\DLS-model\\')
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
ylabel = "Q"+r"$_{N}$"+" (Wm"+r"$^{3}$"+")"
plt.ylabel(ylabel)
plt.savefig("Figures/coolingCurve.png",dpi=400,bbox_inches='tight')
plt.show()


rootgrp = Dataset(str("D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L1_Angle0\ImpurityScanqpll5E7ne1E19_GaussCooling_NoFluxLim\\fi5000E-3\\balance.nc"), "r", format="NETCDF4")
imprad = np.array(rootgrp["b2stel_she_bal"][1]/(np.array(rootgrp['ne'])**2))
imprad = imprad/np.array(rootgrp['vol'])
te = np.array(rootgrp["te"])
te = te/(1.60*10**(-19))
plt.plot(te.flatten(),imprad.flatten(),marker="o",linestyle="")
plt.show()