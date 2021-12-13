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
from unpackConfigurations import unpackConfiguration
from matplotlib.collections import LineCollection
import re
plt.rcParams["font.family"] = "serif"
params = {'legend.fontsize': 'large',
         'axes.labelsize': 'large',
         'axes.titlesize':'large',
         'xtick.labelsize':'large',
         'ytick.labelsize':'large',
        #  'figure.figsize': (8,6),
         }
params = {'legend.fontsize': 'medium',
         'axes.labelsize': 'medium',
         'axes.titlesize':'medium',
         'xtick.labelsize':'medium',
         'ytick.labelsize':'medium',
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

    if routine =="Flux Expansion":
        Files = ["D:\\my stuff\\PhD\\SOLPSRUNS\\isolatedBox\\L1_Angle90\\ImpurityScanqpll5E7ne1E19_Ncooling_Recomb_v0BC\\fi17E-3",
        "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L1_Angle30\ImpurityScanqpll5E7ne1E19_Ncooling_Recomb\\fi20E-3",
        "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L1_Angle11\ImpurityScanqpll5E7ne1E19_Ncooling_Recomb\\fi65E-3",
        "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L1_Angle0\ImpurityScanqpll5E7ne1E19_Ncooling_Recomb\\fi70E-3"
        ]
    if routine =="Connection Length":
        Files = ["D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L0.75Angle0\\fi100E-3",
        "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L1_Angle0\ImpurityScanqpll5E7ne1E19_Ncooling_Recomb\\fi70E-3",
        "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L1.25Angle0\\fi75E-3",
        "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L1.50Angle0\\fi55E-3",
        ]
    if routine =="Flaring":
        Files = [
        # "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L1_Angle0BulgeExperimentNoBulge_Lpar20\ImpurityScanqpll5E7ne1E19_Ncooling_Recomb\\fi135E-3",
        # "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L1_Angle0BulgeExperimentBulge_Lpar20\ImpurityScanqpll5E7ne1E19_Ncooling_Recomb\\fi60E-3",
        "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L1_Angle0BulgeExperimentNoBulge_Lpar20_TightGrid\ImpurityScanqpll5E7ne1.5E19_Ncooling_Recomb\\fi25E-3",
        "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L1_Angle0BulgeExperimentBulge_Lpar20_TightGrid\ImpurityScanqpll5E7ne1.5E19_Ncooling_Recomb\\fi30E-3",

        ]
    if routine =="Kink":
        Files = [
        "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\KinkAngle90Angle0\ImpurityScanqpll5E7ne1E19_Ncooling_Recomb\\fi9E-3",
        "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\KinkAngle67Angle23\ImpurityScanqpll5E7ne1E19_Ncooling_Recomb\\fi10E-3",
        "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\KinkAngle23Angle67\ImpurityScanqpll5E7ne1E19_Ncooling_Recomb\\fi10E-3",
        "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\KinkAngle0Angle90\ImpurityScanqpll5E7ne1E19_Ncooling_Recomb\\fi9E-3",
        ]   
    if routine =="Stability":
        Files = ["D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L1_Angle-10\ImpurityScanqpll5E7ne1E19_Ncooling_Recomb\\fi70E-3"]   
    if routine =="PSI":
        Files = [
        # "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L1_Angle0BulgeExperimentNoBulge_Lpar20\ImpurityScanqpll5E7ne1E19_Ncooling_Recomb\\fi135E-3",
        "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L1.50Angle0\ImpurityScanqpll5E7ne1E19_Ncooling\\fi55E-3",
        "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L1.25Angle0\ImpurityScanqpll5E7ne1E19_Ncooling\\fi75E-3",
        "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L1_Angle0\ImpurityScanqpll5E7ne1E19_Ncooling_Recomb\\fi70E-3",
        "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L0.75Angle0\ImpurityScanqpll5E7ne1E19_Ncooling_Recomb\\fi110E-3",
        "D:\\my stuff\\PhD\\SOLPSRUNS\\isolatedBox\\L1_Angle90\\ImpurityScanqpll5E7ne1E19_Ncooling_Recomb_v0BC\\fi17E-3",
        "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L1_Angle30\ImpurityScanqpll5E7ne1E19_Ncooling_Recomb\\fi20E-3",
        "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L1_Angle11\ImpurityScanqpll5E7ne1E19_Ncooling_Recomb\\fi65E-3",
        "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L1_Angle0\ImpurityScanqpll5E7ne1E19_Ncooling_Recomb\\fi70E-3",
        "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\KinkAngle90Angle0\ImpurityScanqpll5E7ne1E19_Ncooling_Recomb\\fi9E-3",
        "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\KinkAngle67Angle23\ImpurityScanqpll5E7ne1E19_Ncooling_Recomb\\fi10E-3",
        "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\KinkAngle23Angle67\ImpurityScanqpll5E7ne1E19_Ncooling_Recomb\\fi10E-3",
        "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\KinkAngle0Angle90\ImpurityScanqpll5E7ne1E19_Ncooling_Recomb\\fi9E-3",
        "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L1_Angle0BulgeExperimentNoBulge_Lpar20\ImpurityScanqpll5E7ne1E19_Ncooling_Recomb\\fi215E-3",
        "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L1_Angle0BulgeExperimentBulge_Lpar20_TightGrid\ImpurityScanqpll5E7ne1.5E19_Ncooling_Recomb\\fi30E-3",

        ]
    C = []
    Bt = []
    tu = []

    if routine in ["Flux Expansion","Kink","Stability"]:
        fig, axs = plt.subplots(1, 1)
        axs.set_aspect('equal')
    elif routine == "Connection Length":
        fig = plt.figure(figsize=(3.5,5))
        
        gs = fig.add_gridspec(1, len(Files), hspace=0, wspace=0)
        axs = gs.subplots(sharey='row')
        # fig, axs = plt.subplots(1, len(Files),sharey='row')
        for i in range(len(axs)):
            axs[i].set_aspect('equal')
    elif routine == "Flaring":
        fig = plt.figure(figsize=(5,7))
        
        gs = fig.add_gridspec(1, len(Files), hspace=0, wspace=0)
        axs = gs.subplots(sharey='row')
        # fig, axs = plt.subplots(1, len(Files),sharey='row')
        for i in range(len(axs)):
            axs[i].set_aspect('equal')
    elif routine == "PSI":
        fig = plt.figure(figsize=(15,25))
        
        gs = fig.add_gridspec(1, 4, hspace=0, wspace=0)
        axs = gs.subplots(sharey='row')
        # fig, axs = plt.subplots(1, len(Files),sharey='row')
        for i in range(len(axs)):
            axs[i].set_aspect('equal')
    PSIAXIS = 0
    colornum = 0
    for File in Files:
        #unpack position data from SOLPS balance file (SOL ring chosen was to best match starting coordinates):

        rootgrp = Dataset(str(File)+"/balance.nc", "r", format="NETCDF4")

        # unpack dimensions of super-x
        r = rootgrp['crx']
        z = rootgrp['cry']

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
        if routine=="Flaring":
            istring = "IV-"
            axs[colornum].add_collection(LineCollection(segs1,linewidth=0.5,color =gridcolors[2*colornum],label=istring+str(colornum+1)))
            axs[colornum].add_collection(LineCollection(segs2,linewidth=0.5,color =gridcolors[2*colornum]))
        if routine=="PSI":
            istring = "IV-"
            axs[PSIAXIS].add_collection(LineCollection(segs1,linewidth=0.5,color =gridcolors[colornum],label=istring+str(colornum+1)))
            axs[PSIAXIS].add_collection(LineCollection(segs2,linewidth=0.5,color =gridcolors[colornum]))
        if routine=="Connection Length":
            istring = "IV-"
            axs[colornum].add_collection(LineCollection(segs1,linewidth=0.5,color =gridcolors[colornum],label=istring+str(colornum+1)))
            axs[colornum].add_collection(LineCollection(segs2,linewidth=0.5,color =gridcolors[colornum]))
        if routine == "Flux Expansion":
            istring = "II-"
            axs.add_collection(LineCollection(segs1,linewidth=0.5,color =gridcolors[colornum],label=istring+str(colornum+1)))
            axs.add_collection(LineCollection(segs2,linewidth=0.5,color =gridcolors[colornum]))
        if routine == "Kink":
            istring = "III-"
            # axs.plot(rbr[14],zbr[14])
            axs.add_collection(LineCollection(segs1,linewidth=0.5,color =gridcolors[colornum],label=istring+str(colornum+1)))
            axs.add_collection(LineCollection(segs2,linewidth=0.5,color =gridcolors[colornum]))
        if routine == "Stability":
            axs.add_collection(LineCollection(segs1,linewidth=0.5,color =gridcolors[2],label=istring+str(colornum+1)))
            axs.add_collection(LineCollection(segs2,linewidth=0.5,color =gridcolors[2]))
        colornum +=1
        if colornum==4:
            colornum=0
            PSIAXIS = PSIAXIS+1
        
    #plotting and layout specifics for expansion case
    if routine=="Flux Expansion":
        plt.xlabel("R (m)")
        plt.xlim([0.9,2.03])
        plt.ylim([-1.05,0.25])  
        plt.ylabel("Z (m)")
        plt.tight_layout()
        plt.legend(bbox_to_anchor=(0.98, 0.73))
        plt.savefig("Figures/gridExpansion.png",dpi=400,bbox_inches='tight')
    
    if routine=="Kink":
        plt.xlabel("R (m)")
        plt.xlim([0.9,2.22])
        plt.ylim([-1.15,0.25])  
        plt.ylabel("Z (m)")
        plt.tight_layout()
        plt.legend(bbox_to_anchor=(0.98, 0.99))
        plt.savefig("Figures/gridavB.png",dpi=400,bbox_inches='tight')
    if routine=="Stability":
        plt.xlabel("R (m)")
        plt.xlim([0.58,1.2])
        plt.ylim([-1.1,0.1])  
        plt.ylabel("Z (m)")
        plt.tight_layout()
        plt.savefig("Figures/gridStability.png",dpi=400,bbox_inches='tight')

    #plotting and layout specifics for flaring case
    if routine=="Flaring":
        filewalls = "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L1_Angle0BulgeExperimentNoBulge_Lpar20\ImpurityScanqpll5E7ne1E19_Ncooling_Recomb\\fi135E-3\input.dat"
        filewalls = "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L1_Angle0BulgeExperimentNoBulge_Lpar20_TightGrid\ImpurityScanqpll5E7ne1.5E19_Ncooling_Recomb\\fi25E-3\input.dat"
        plotWALL(filewalls,axs[0])
        filewalls = "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L1_Angle0BulgeExperimentBulge_Lpar20\ImpurityScanqpll5E7ne1E19_Ncooling_Recomb\\fi60E-3\input.dat"
        filewalls = "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L1_Angle0BulgeExperimentBulge_Lpar20_TightGrid\ImpurityScanqpll5E7ne1.5E19_Ncooling_Recomb\\fi20E-3\input.dat"
        plotWALL(filewalls,axs[1])
        for i in range(len(Files)):
            axs[i].set_xlim([0.5,1.55])
            axs[i].set_ylim([-1.03,0.05])  
            axs[i].xaxis.set_ticks(np.arange(1, 1.2, 0.249))
        p0,=axs[0].plot(-10,-10,color=gridcolors[0])
        p1,=axs[0].plot(-10,-10,color=gridcolors[3])
        p2,=axs[0].plot(-10,-10,color="C1")
        axs[0].set_ylabel("Z (m)")
        axs[0].set_xlabel("R (m)",x=1)  
        fig.tight_layout(pad=.0)
        plt.legend([p0,p1,p2],["Straightdown","Flared","Wall Surface"],bbox_to_anchor=(0.8, 0.9))
        plt.savefig("Figures/gridBulge.png",dpi=400,bbox_inches='tight')
    if routine=="PSI":
        for i in range(4):
            axs[i].set_xlim([0.5,2.3])
            axs[i].set_ylim([-1.55,0.25])  
            axs[i].xaxis.set_ticks(np.arange(1, 1.2, 0.249))
        p0,=axs[0].plot(-10,-10,color=gridcolors[0])
        p1,=axs[0].plot(-10,-10,color=gridcolors[3])
        p2,=axs[0].plot(-10,-10,color="C1")
        axs[0].set_ylabel("Z (m)")
        axs[0].set_xlabel("R (m)",x=2)  
        fig.tight_layout(pad=.0)
        plt.savefig("Figures/gridPSI.png",dpi=400,bbox_inches='tight')

    #plotting and layout specifics for length case
    if routine=="Connection Length":
        for i in range(len(Files)):
            axs[i].set_xlim([0.92,1.23])
            axs[i].set_ylim([-1.55,0.05])  
            axs[i].xaxis.set_ticks(np.arange(1, 1.5, 5))
        p0,=axs[0].plot(-10,-10,color=gridcolors[0])
        p1,=axs[0].plot(-10,-10,color=gridcolors[1])
        p2,=axs[0].plot(-10,-10,color=gridcolors[2])
        p3,=axs[0].plot(-10,-10,color=gridcolors[3])
        axs[0].set_ylabel("Z (m)")
        axs[0].set_xlabel("R (m)",x=2)  
        fig.tight_layout(pad=.0)
        plt.legend([p0,p1, p2,p3],["I-1", "I-2","I-3","I-4"],bbox_to_anchor=(-1.9, 0.3))
        plt.savefig("Figures/gridLength.png",dpi=400,bbox_inches='tight')

    plt.show()


routine = "Flaring"
# routine = "Flux Expansion"
# routine = "Connection Length"
# routine = "Kink"
# routine = "Stability"
routine = "PSI"
plotGrids(routine)