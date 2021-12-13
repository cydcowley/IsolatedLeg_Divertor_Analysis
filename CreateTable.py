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
from UnpackSOLPS import unpackSOLPS
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

fileNames = [
    "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L0.75Angle0\ImpurityScanqpll5E7ne1E19_Ncooling_Recomb\\fi160E-3\\balance.nc",
    "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L1_Angle0\ImpurityScanqpll5E7ne1E19_Ncooling_Recomb\\fi110E-3Backward\\balance.nc",
"D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L1_Angle90\ImpurityScanqpll5E7ne1E19_Ncooling_Recomb_v0BC\\fi20E-3Backwards\\balance.nc",
]

for fileName in fileNames:
    quantities2d,SOLring1 = unpackSOLPS(fileName, -1,14)

    indexfront = np.argmax(np.gradient(SOLring1.cond/SOLring1.B))
    PE = SOLring1.ne*SOLring1.te
    pressureDiff =PE[indexfront]/PE[-1]
    qdiff = trapz(SOLring1.qf,SOLring1.Spar)/trapz(SOLring1.B,SOLring1.cond/SOLring1.B)
    Ttest = np.linspace(0,100,10000)
    Q = []
    for t in Ttest:
        Q.append(LfuncN(t))
    tint = trapz(Ttest**(0.5)*Q,Ttest)
    plt.plot(SOLring1.te,PE/PE[-1])
    # plt.plot(SOLring1.Spar[indexfront],PE[indexfront]/PE[-1],marker="o")
    # plt.plot(SOLring1.Spar,SOLring1.cond)
    print(pressureDiff)
    print(qdiff)
    print("t int is",tint)
plt.show()