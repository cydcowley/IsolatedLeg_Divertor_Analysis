from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import os
from natsort import natsorted
import matplotlib.image as m
import imageio as io
from scipy import interpolate
from scipy.integrate import quad,trapz, cumtrapz, odeint, solve_ivp
from PIL import Image
import cv2
import glob
import mpltools
from mpltools import special
import sys
from UnpackSOLPS import unpackSOLPS
sys.path.append('D:\\my stuff\PhD\\recreation of detachment sensitivity paper\\')
sys.path.append('D:\\my stuff\PhD\\ThermalFrontFigures\\')
sys.path.append('d:\\my stuff\\PhD\\DLS-model\\')


Exhemptfiles = ["baserun","notes.txt","ref","INPUT.m",
    # "fi435E-3",
    # "fi405E-3",
    # "fi375E-3",
    # "fi345E-3",
    # "fi315E-3",
    # "fi290E-3",
    # "fi155E-3",
    # "fi175E-3",
    # "fi135E-3",
    # "fi115E-3",
    "fi100E-3",
    # "fi85E-3",
    # "fi70E-3",
    # "fi60E-3",
    "fi17E-3",
    "fi19E-3",
    # "fi22E-3Backwards",
    "fi30E-3Backwards",
    "fi33E-3",
    "fi38E-3",
    "fi43E-3",
    "fi50E-3",
    "fi55E-3Backwards",
    "fi60E-3",
    "fi70E-3",
    # "fi80E-3",
    # "fi90E-3",
    # "fi140E-3Backward",
    ]

params = {'legend.fontsize': 'medium',
         'axes.labelsize': 'medium',
         'axes.titlesize':'medium',
         'xtick.labelsize':'medium',
         'ytick.labelsize':'medium',
        #  'figure.figsize': (4,3.2),
         }
plt.rcParams.update(params)

def makePlot(folder):
    Rrange = [0.5,1.7]
    Zrange = [-2.1,-1.3]
    # Rrange = [0.2,0.7]
    # Zrange = [-1.5,-1.3]
    reverse = -1
    Files = os.listdir(str(folder))
    Files = natsorted(Files)
    if "powerScan" in folder:
        Files = Files[::-1]
    # define frame containers for making movies/gifs
    image_listte = []
    image_listte2d = []
    image_listn0 = []
    image_listheat = []
    image_listne = []
    image_listSource = []
    image_listheat2d = []
    counter = 0
    
    for File in Files:
        print(File)
        if File in Exhemptfiles:
            print("hello")
            continue
        print(File)
        fileName = str(folder)+"/"+str(File)+"/" + str("balance.nc")
        Xpoint = -1
        quantities2d,SOLring = unpackSOLPS(fileName, reverse,14)

        plotquantities = ["te2d"]
        if "n02d" in plotquantities:
            plt.pcolor(quantities2d["r"],quantities2d["z"],quantities2d["n0"])
            plt.plot( quantities2d["r"][quantities2d["ring"]][::reverse][0:Xpoint], quantities2d["z"][quantities2d["ring"]][::reverse][0:Xpoint],color="C1")
            plt.axis('equal')
            plt.clim([10**18,10**19])
            plt.colorbar()
            plt.xlim(Rrange)
            plt.ylim(Zrange)
            plt.title("Neutral Density, "+str(File))
            plt.xlabel("R")
            plt.ylabel("Z")
            imgname = "ne2d.png"
            plt.savefig(imgname,dpi=400)
            imgt2 = cv2.imread(imgname)
            image_listn0.append(imgt2)
            # plt.clim([0,1.75*10**20])
            plt.show()
        
        if "te2d" in plotquantities:
            plt.set_cmap("viridis")
            plt.pcolor(quantities2d["r"],quantities2d["z"],
            quantities2d["te"])
            # plt.plot( quantities2d["r"][quantities2d["ring"]][::reverse][0:Xpoint], quantities2d["z"][quantities2d["ring"]][::reverse][0:Xpoint],color="C1")
            plt.clim([2,50])
            plt.colorbar()
            
            plt.xlim([0.90,1.2])
            plt.ylim([-1,0])
            plt.axis('scaled')
            # plt.tight_layout()
            # plt.title("Electron Temperature, "+str(File))
            plt.title("Electron Temperature")
            imgname = "te2d"+str(counter)+".png"
            plt.savefig(imgname,dpi=400)

            imgt2 = cv2.imread(imgname)
            image_listte2d.append(imgt2)
            plt.show()

        if "heats2d" in plotquantities:
            plt.pcolor(quantities2d["r"],quantities2d["z"],-1*quantities2d["imprad"])
            plt.colorbar()
            plt.set_cmap("plasma")
            # plt.clim([200,500])
            # plt.xlim(Rrange)
            # plt.ylim(Zrange)
            plt.title("Impurity Radiation")
            plt.xlabel("R (m)")
            plt.ylabel("Z (m)")
            plt.axis('scaled')
            plt.savefig("heat2d.png",dpi=400)
            
            imgheat2d = cv2.imread("heat2d.png")
            image_listheat2d.append(cv2.cvtColor(imgheat2d, cv2.COLOR_RGB2BGR))
            plt.show()
        counter = counter+1
    if "te1d" in plotquantities:
        frame0 = np.array(image_listte[0])
        height, width, layers = frame0.shape
        size = (width,height)
        out = cv2.VideoWriter('te.avi',cv2.VideoWriter_fourcc(*'DIVX'), 15, size)
        for img in image_listte:
            out.write(img)
        out.release()
        cv2.destroyAllWindows()
    if "te2d" in plotquantities:
        giflist = []
        
        frame0 = np.array(image_listte2d[0])
        height, width, layers = frame0.shape
        size = (width,height)
        out = cv2.VideoWriter('te.avi',cv2.VideoWriter_fourcc(*'DIVX'), 2, size)
        for img in image_listte2d:
            giflist.append(cv2.cvtColor(img, cv2.COLOR_BGR2RGB))
            out.write(img)
        out.release()
        io.mimsave('heat.gif',giflist,duration=1)
        cv2.destroyAllWindows()
    if "n02d" in plotquantities:
        frame0 = np.array(image_listn0[0])
        height, width, layers = frame0.shape
        size = (width,height)
        out = cv2.VideoWriter('n0.avi',cv2.VideoWriter_fourcc(*'DIVX'), 15, size)
        for img in image_listn0:
            out.write(img)
        out.release()
        cv2.destroyAllWindows()
    if "partSource" in plotquantities:
        frame0 = np.array(image_listSource[0])
        height, width, layers = frame0.shape
        size = (width,height)
        out = cv2.VideoWriter('source.avi',cv2.VideoWriter_fourcc(*'DIVX'), 15, size)
        for img in image_listSource:
            out.write(img)
        out.release()
        cv2.destroyAllWindows()
    if "heats" in plotquantities:
        frame0 = np.array(image_listheat[0])
        height, width, layers = frame0.shape
        size = (width,height)
        out = cv2.VideoWriter('heat.avi',cv2.VideoWriter_fourcc(*'DIVX'), 15, size)
        for img in image_listheat:
            out.write(img)
        out.release()
        cv2.destroyAllWindows()
    if "ne1d" in plotquantities:
        io.mimsave('ne.gif', image_listne, duration=0.5)
    if "heats2d" in plotquantities:
        io.mimsave('heat.gif',image_listheat2d,duration=0.2)


# makePlot(folder="D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L1_Angle0\ImpurityScanqpll5E7ne1E19_Ncooling_Recomb")
makePlot(folder="D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L1_Angle90\ImpurityScanqpll5E7ne1E19_Ncooling_Recomb_v0BC")
# makePlot(folder="D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L1_Angle0BulgeExperimentNoBulge_Lpar20\ImpurityScanqpll5E7ne1E19_Ncooling_Recomb")
# makePlot(folder="D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L1_Angle0BulgeExperimentBulge\ImpurityScanqpll5E7ne1E19_Ncooling_Recomb")
