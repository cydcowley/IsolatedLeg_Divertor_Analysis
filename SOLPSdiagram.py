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
from mpl_toolkits.axes_grid1 import make_axes_locatable
import sys
from UnpackSOLPS import unpackSOLPS
from unpackConfigurations import unpackConfiguration
from matplotlib.collections import LineCollection
import re
plt.rcParams["font.family"] = "serif"
plt.rcParams['font.size'] = 9
params = {'legend.fontsize': 'large',
         'axes.labelsize': 'large',
         'axes.titlesize':'large',
         'xtick.labelsize':'large',
         'ytick.labelsize':'large',
        #  'figure.figsize': (8,6),
         }
params = {'legend.fontsize': 'small',
         'axes.labelsize': 'small',
         'axes.titlesize':'small',
         'xtick.labelsize':'small',
         'ytick.labelsize':'small',

        #  'figure.figsize': (4,3.2),
         }
plt.rcParams.update(params)


def plotWALL(filewalls,axis):
    datawall = open(filewalls)
    tline = datawall.readlines(1)
    linenum = 0
    while True:
        tline = datawall.readlines(1)
        if  '*** 3b. Data for additional surfaces' in tline[0]:
            break
    ns = int(datawall.readlines(1)[0])
    coords = []
    ic = 1
    for i in range(ns):
        displayname = datawall.readlines(1)
        rlbnd = int(datawall.readlines(1)[0][1])
        tmp = datawall.readlines(1)
        iliin = int(tmp[0][1:6])

        if iliin <0:
            tmp = datawall.readlines(1)
            continue
        if rlbnd == 2:
            check = datawall.readlines(1)[0]
            co = re.findall("\D+\d+\.\d+\D+\d+",check)
            for i in range(len(co)):
                co[i] = float(co[i])/100
            coords.append(co)
            stype = datawall.readlines(1)
            axis.plot([co[0],co[3]],[co[1],co[4]],color="C1")



def plotGrids(routine):
    gridcolors = ["#53ba83","#059b9a","#095169","#0c0636","000000"]
    Files = []
    Files = ["D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L1_Angle30\ImpurityScanqpll5E7ne1E19_Ncooling_Recomb\\fi80E-3"]

    C = []
    Bt = []
    tu = []




    # fig, axs = plt.subplots(1, len(Files),sharey='row')
    fig = plt.figure()
    axs1 = fig.add_subplot(121)
    axs2 = fig.add_subplot(222)
    axs3 = fig.add_subplot(224)
    axs1.set_aspect('equal')
    colornum = 0
    for File in Files:
        #unpack position data from SOLPS balance file (SOL ring chosen was to best match starting coordinates):

        rootgrp = Dataset(str(File)+"/balance.nc", "r", format="NETCDF4")

        # unpack dimensions of super-x
        r = rootgrp['crx']
        z = rootgrp['cry']
        R = (rootgrp['crx'][0]+rootgrp['crx'][1]+rootgrp['crx'][2]+rootgrp['crx'][3])/4
        Z = (rootgrp['cry'][0]+rootgrp['cry'][1]+rootgrp['cry'][2]+rootgrp['cry'][3])/4
        T = rootgrp["te"]
        ne = np.array(rootgrp["ne"])*10**(-19)
        bb = rootgrp['bb']
        dv = np.array(rootgrp['vol'])
        quantities2d = {}
        
        #grid coordinates:
        quantities2d["r"] = (rootgrp['crx'][0]+rootgrp['crx'][1]+rootgrp['crx'][2]+rootgrp['crx'][3])/4
        quantities2d["z"] = (rootgrp['cry'][0]+rootgrp['cry'][1]+rootgrp['cry'][2]+rootgrp['cry'][3])/4
        hx = rootgrp['hx']

        #total magnetic field:
        quantities2d["TotalField"] = np.array(bb[3])
        quantities2d["Bpol"] = np.array(bb[0])
        #length of grid
        s = np.array(hx)*np.abs(np.array(bb[3])/np.array(bb[0]))
        quantities2d["sdiff"] = s
        quantities2d["Area"] = np.array(dv)/s
        cond = rootgrp['fhe_cond'][0]/quantities2d["Area"]
        plt.set_cmap("viridis")
        print(R[:,0][40])
        print(Z[:,0][40])
        axs3.plot(np.sqrt((R[:,0]-R[:,0][0])**2+(Z[:,0]-Z[:,0][0])**2),cond[:,1]/(10**6))
        axs2.plot(np.sqrt((R[:,0]-R[:,0][0])**2+(Z[:,0]-Z[:,0][0])**2),ne[:,1])
        xpoint = [np.sqrt((R[:,0][10]-R[:,0][0])**2+(Z[:,0][10]-Z[:,0][0])**2),np.sqrt((R[:,0][10]-R[:,0][0])**2+(Z[:,0][10]-Z[:,0][0])**2)]
        axs2.plot(xpoint,[-1,2],linestyle="--",color="r")
        axs3.plot(xpoint,[-10,70],linestyle="--",color="r")
        
        mesh = axs1.pcolormesh(R,Z,np.array(T)/(1.60*10**(-19)))

        cbar = fig.colorbar(mesh, ax=axs1,shrink=0.7,orientation="horizontal")
        cbar.set_label('T (eV)')
        # plt.colorbar(mesh,cax=axs1)

        #unpack the grids for straight box
        rbl = r[0]
        rbr = r[1]
        rtl = r[2]
        rtr = r[3]
        zbl = z[0]
        zbr = z[1]
        ztl = z[2]
        ztr = z[3]
        istring = "I-"
        #create grid lines
        segs1 = np.stack((rbr,zbr), axis=2)
        segs2 = segs1.transpose(1,0,2)
        #plot the grid lines

        istring = "IV-"

        colornum +=1
        

    #plotting and layout specifics for flaring case
    if routine=="Flaring":
        filewalls = "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L1_Angle30\ImpurityScanqpll5E7ne1E19_Ncooling_Recomb\\fi20E-3\\input.dat"
        plotWALL(filewalls,axs1)
        axs1.plot(np.transpose(R)[-1],np.transpose(Z)[-1],color="C1")
        axs1.scatter(R[:,0][10],Z[:,0][10],s=60,facecolors='none', edgecolors='r',linestyle="--")

        axs1.set_xlim([0.84,1.88])
        axs1.set_ylim([-1.03,0.28])  
        axs2.set_xlim([-0.01,0.26])
        axs3.set_xlim([-0.01,0.26])
        axs2.set_ylim([0.1,1.4])  
        axs3.set_ylim([-1,50])
        
        #     axs.xaxis.set_ticks(np.arange(1, 1.2, 0.249))
        # p0,=axs.plot(-10,-10,color=gridcolors[0])
        # p1,=axs.plot(-10,-10,color=gridcolors[3])
        # p2,=axs.plot(-10,-10,color="C1")
        axs1.annotate("pump", xy=(1.65, -0.73), xytext=(1.6, -0.55),arrowprops=dict(arrowstyle="->"))
        axs1.annotate("x-point", xy=(1, 0.02), xytext=(0.9, 0.16),arrowprops=dict(arrowstyle="->"))
        axs1.annotate("west wall", xy=(1, 0), xytext=(0.88, -0.6))
        axs1.annotate("east wall", xy=(1, 0), xytext=(1.45, -0.3))
        axs1.annotate("target", xy=(1.5, -0.85), xytext=(1.55, -0.9))
        axs1.annotate("a)", xy=(1.5, -0.85), xytext=(1.65, 0.15))
  
        axs1.set_ylabel("Z (m)")
        axs1.set_xlabel("R (m)")  
        axs3.set_xlabel("s"+r"$_{⟂}$" +" (m)")  
        axs3.set_ylabel("q"+r"$_{||}$" +" (MWm"+r"$^{-2}$" +")") 
        axs2.set_ylabel("n (10"+r"$^{19}$" +"m"+r"$^{-3}$" +")")   
        axs2.set_xticks([])
        axs2.annotate("b)", xy=(0.235, 1.2), xytext=(0.235, 1.2))
        axs2.annotate("x-point", xy=(0.05, 0.5), xytext=(0.07, 0.3),arrowprops=dict(arrowstyle="->"))
        
        axs3.annotate("c)", xy=(0.1, 30), xytext=(0.235, 43))
        axs3.annotate("x-point", xy=(0.05, 12), xytext=(0.07, 4),arrowprops=dict(arrowstyle="->"))
        
        plt.tight_layout(pad=0.2)
        # plt.legend([p0,p1,p2],["Straightdown","Flared","Wall Surface"],bbox_to_anchor=(0.8, 0.9))
        plt.savefig("Figures/SOLPSdiagram.png",dpi=400,bbox_inches='tight')



routine = "Flaring"
# routine = "Flux Expansion"
# routine = "Connection Length"
# routine = "Kink"
# routine = "Stability"
plotGrids(routine)