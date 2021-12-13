from matplotlib import lines
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import os
import json
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
from parse import parse
sys.path.append('D:\\my stuff\\PhD\\Theoretical_Detachment_Control_Scripts')
from LipschultzDLS import ChInt,ChIntThermalForce,averageB
from UnpackSOLPS import unpackSOLPS,SOLring
from AnalyticCoolingCurves import LfuncN
import re
from scipy.optimize import curve_fit
plt.rcParams["font.family"] = "serif"
params = {'legend.fontsize': 'medium',
         'axes.labelsize': 'medium',
         'axes.titlesize':'medium',
         'xtick.labelsize':'medium',
         'ytick.labelsize':'medium',
        #  'figure.figsize': (4,3.2),
         }
plt.rcParams.update(params)

customOrder = 0

colors = ["#356288","#fe1100","#aacfdd","#fe875d"]
gridcolors = ["#53ba83","#059b9a","#095169","#0c0636","000000"]


def determineC0(Spar,C,L):

    for k in range(len(C)):

        if Spar[k] >= L/100:
            return k-1



def scanThresholds(routine):
    folderList = []
    #define folder list for each routine
    if routine=="Flux Expansion":
        folderList = [
        # "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L1_Angle90\ImpurityScanqpll5E7ne1E19_Ncooling_Recomb_v0BC",
        # "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L1_Angle30\ImpurityScanqpll5E7ne1E19_Ncooling_Recomb",
        # "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L1_Angle11\ImpurityScanqpll5E7ne1E19_Ncooling_Recomb",
        # "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L1_Angle0\ImpurityScanqpll5E7ne1E19_Ncooling_Recomb",
        "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L1_Angle90\DensityScanqpll5E7fi70E-3_Ncooling",
        "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L1_Angle30\DensityScanqpll5E7fi70E-3_Ncooling",
        "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L1_Angle11\DensityScanqpll5E7fi70E-3_Ncooling",
        "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L1_Angle0\DensityScanqpll5E7fi70E-3_Ncooling",  
        ]
    if routine=="Connection Length" or routine=="Length Difference":
        folderList = [
            "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L0.75Angle0\ImpurityScanqpll5E7ne1E19_Ncooling_Recomb",
            "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L1_Angle0\ImpurityScanqpll5E7ne1E19_Ncooling_Recomb",
            "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L1.25Angle0\ImpurityScanqpll5E7ne1E19_Ncooling",
            "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L1.50Angle0\ImpurityScanqpll5E7ne1E19_Ncooling",
        # "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L0.75Angle0\DensityScanqpll5E7fi70E-3_Ncooling",
        # "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L1_Angle0\DensityScanqpll5E7fi70E-3_Ncooling",
        # "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L1.25Angle0\DensityScanqpll5E7fi70E-3_Ncooling",
        # "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L1.50Angle0\DensityScanqpll5E7fi70E-3_Ncooling",
        
        ]
    if routine=="Averaged Field":
        folderList = [
            "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\KinkAngle90Angle0\ImpurityScanqpll5E7ne1E19_Ncooling_Recomb",
            "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\KinkAngle67Angle23\ImpurityScanqpll5E7ne1E19_Ncooling_Recomb",
            "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\KinkAngle23Angle67\ImpurityScanqpll5E7ne1E19_Ncooling_Recomb",
            "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\KinkAngle0Angle90\ImpurityScanqpll5E7ne1E19_Ncooling_Recomb"
        ]


    reverse = -1
    counter = 0
    #create empty arrays to store information about each box, such as detachment threshold
    Thresholds = []
    ThresholdError = []
    averageDensity = []
    averageimpurityPower = []
    CSimple = []
    avB = []
    Fr = []
    Lt = []
    percImpurity = []

    Exhemptfiles = ["baserun","ref","notes.txt","INPUT.m",]#files to be skipped

    #iterate through each geometry
    for folder in folderList:
        print(folder)
        Files = os.listdir(str(folder))
        Files = natsorted(Files)

        if "powerScan" in folder:
            Files = Files[::-1]

        #create empty arrays to store information about each simulation, such as front position
        C = []
        Sh = []
        Sherr = []
        Bt = []
        Bx = []
        Heats = []
        nav = []
        qftot = []

        #iterate through each simulation in the scan
        for File in Files:

            if File in Exhemptfiles:
                continue

            fileName = str(folder)+"/"+str(File)+"/" + str("balance.nc")
            rootgrp = Dataset(str(fileName), "r", format="NETCDF4")
            # print(rootgrp)
            Xpoint = -1
            SOLring1 = 0

            #unpack data from SOLPS
            quantities2d,SOLring1 = unpackSOLPS(fileName, reverse,15)

            #determine front position and errors
            pos,upper,lower = SOLring1.calcFrontqpll()

            Sh.append(pos)
            Sherr.append([pos-lower,upper-pos])
            # calculate the control parameter for this SOL ring
            C.append(SOLring1.determineC())
            Bt.append(SOLring1.B[0])
            Bx.append(SOLring1.B[-1])
            Heats.append(SOLring1.performHeatAnalysis())
            nav.append( trapz(SOLring1.qf,SOLring1.Spar)/trapz(SOLring1.qf/SOLring1.ne,SOLring1.Spar))
            qftot.append(trapz(SOLring1.qf,SOLring1.Spar))
            if lower >20/100:
                print(File)
    
        Sh = np.array(Sh)
        C = np.array(C)
        QF = []
        QI = []
        for i in range(len(Heats)):
            QF.append(Heats[i]["impurity"])
            QI.append(Heats[i]["conducted"])

        #determine the detachment threshold of the data

        L = SOLring1.Spar[-1]
        index0 = determineC0(Sh-np.array(Sherr)[:,0],C,L)

        averageDensity.append(nav[index0])
        averageimpurityPower.append(qftot[index0])
        # print(Sh-np.array(Sherr)[:,0])
        Thresholds.append((C[index0]))

        ThresholdError.append(C[index0+1] - C[index0])
        
        Fr.append(np.mean(Bt)/np.mean(Bx))
        Lt.append(SOLring1.Spar[-1])
        percImpurity.append(-1*QF[index0]*np.mean(Bt)/(50E6*0.5))
        CSimple.append(ChInt(SOLring1.Spar, SOLring1.B, SOLring1.Spar[-1],SOLring1.Spar[-1] ,0))
        avB.append(averageB(SOLring1.Spar, SOLring1.B, SOLring1.Spar[-1],SOLring1.Spar[-1] ,0))
        counter = counter+1

        

    avB = np.array(avB)
    Lt = np.array(Lt)
    #plotting specifics for box rotation routine
    if routine=="Flux Expansion":
        label0 =  "C"+r"$_{t,SOLPS}$"
        for i in range(len(Thresholds)):
            if i==0:
                plt.annotate("I-"+str(i+1),(Fr[i],Thresholds[i]/Thresholds[0]+0.08))
                plt.errorbar((Fr[i]),(Thresholds[i]/Thresholds[0]),[[0],[ThresholdError[i]/Thresholds[0]]],label=label0,marker="o",color = gridcolors[i])
            else:
                plt.annotate("I-"+str(i+1),(Fr[i]-0.04,Thresholds[i]/Thresholds[0]))
                plt.errorbar((Fr[i]),(Thresholds[i]/Thresholds[0]),[[Thresholds[i]/Thresholds[0]-Thresholds[i]/(Thresholds[0]+ThresholdError[0])],[ThresholdError[i]/Thresholds[0]]],marker="o",color = gridcolors[i])
        label1 = "F"+r"$_{R}^{-1}$"
        label2 = "C"+r"$_{t, DLS}$"
        label3 = r"$\left< B \right>_{t}^{-2/7}$"
        plt.plot((Fr),(Fr/Fr[0]),color="C0",linestyle="--",label=label1)
        plt.plot((Fr),(CSimple/CSimple[0]),color="C1",label=label2)
        plt.plot((Fr),(avB/avB[0])**(-2/7),color="C4",linestyle="-.",label=label3)
        plt.xlabel(label1)
        plt.legend()
        plt.tight_layout()
        plt.savefig("Figures/thresholdExpansion.png",dpi=400,bbox_inches='tight')
        plt.show()

    #plotting specifics for connection length variation routine
    if routine=="Connection Length":
        label0 =  "C"+r"$_{t,SOLPS}$"
        for i in range(len(Thresholds)):
            if i==0:
                plt.annotate("I-"+str(i+1),(Lt[i],Thresholds[i]/Thresholds[0]-0.05))
                plt.errorbar(Lt[i],(Thresholds[i]/Thresholds[0]),[[0],[ThresholdError[i]/Thresholds[0]]],label=label0,marker="o",color = gridcolors[i])
            else:
                plt.annotate("I-"+str(i+1),(Lt[i]-1,Thresholds[i]/Thresholds[0]))
                plt.errorbar(Lt[i],(Thresholds[i]/Thresholds[0]),[[Thresholds[i]/Thresholds[0]-Thresholds[i]/(Thresholds[0]+ThresholdError[0])],[ThresholdError[i]/Thresholds[0]]],marker="o",color = gridcolors[i])
        label1 = "L"+r"$_{||}^{-2/7}$"
        label2 = "L"+r"$_{||}$"+" (m)"
        
        Lplot = np.linspace(np.amin(Lt),np.amax(Lt),100)
        plt.plot(Lplot,(Lplot/Lplot[0])**(-2/7),color="C3",label=label1)
        def testfunc(x,a):
            return x**a

        popt, pcov = curve_fit(testfunc, Lt/Lt[0], Thresholds/Thresholds[0],p0=[-5/7])
        label3 = "L"+r"$_{||}^{"+str(int(popt[0]*70)/10)+r"/7}$"
        print(popt)
        plt.plot(Lplot, testfunc(Lplot/Lplot[0], *popt),label=label3,color="C6",linestyle="--")
        plt.xlabel(label2)
        plt.legend()
        plt.tight_layout()
        plt.savefig("Figures/thresholdLength.png",dpi=400,bbox_inches='tight')
        plt.show()
    # if routine=="Length Difference":
    #     label0 =  "C"+r"$_{t,SOLPS}$"
    #     for i in range(len(Thresholds)):
    #         if i==0:
    #             plt.annotate("I-"+str(i+1),(Lt[i],Thresholds[i]/Thresholds[0]-0.05))
    #             plt.errorbar(Lt[i],(Thresholds[i]/Thresholds[0]),ThresholdError[i]/Thresholds[0],label=label0,marker="o",color = gridcolors[i])
    #         else:
    #             plt.annotate("I-"+str(i+1),(Lt[i]-1,Thresholds[i]/Thresholds[0]))
    #             plt.errorbar(Lt[i],(Thresholds[i]/Thresholds[0]),ThresholdError[i]/Thresholds[0],marker="o",color = gridcolors[i])
    #     label1 = "L"+r"$_{t}^{-2/7}$"
    #     Lplot = np.linspace(np.amin(Lt),np.amax(Lt),100)
    #     plt.plot(Lplot,(Lplot/Lplot[0])**(-2/7),color="C3",label=label1)
    #     plt.plot(Lt,averageDensity[0]/averageDensity,label="average ne")
    #     plt.plot(Lt,averageimpurityPower/averageimpurityPower[0],label="QF")
        
    #     plt.xlabel(label1)
    #     plt.legend()
    #     plt.savefig("Figures/thresholdLength.png",dpi=400,bbox_inches='tight')
    #     plt.show()

    #plotting specifics for kink routine
    if routine=="Averaged Field":
        label0 =  "C"+r"$_{t,SOLPS}$"
        for i in range(len(Thresholds)):
            if i==0:
                plt.annotate("III-"+str(i+1),(avB[i],Thresholds[i]/Thresholds[0]-0.05))
                plt.errorbar(avB[i],(Thresholds[i]/Thresholds[0]),[[0],[ThresholdError[i]/Thresholds[0]]],label=label0,marker="o",color = gridcolors[i])
            else:
                plt.annotate("III-"+str(i+1),(avB[i]-0.02,Thresholds[i]/Thresholds[0]))
                plt.errorbar(avB[i],(Thresholds[i]/Thresholds[0]),[[0],[ThresholdError[i]/Thresholds[0]]],marker="o",color = gridcolors[i])
        label0 = r"$\left<" +"B"+r" \right>_{"+"above t"+r"}^{-2/7}$"
        label1 = "L"+r"$_{t}^{-2/7}$"
        label2 = "C"+r"$_{t, DLS}$"
        label3 = r"$\left< B \right>_{"+"above t"+r"}^{-2/7}$"
        plt.plot(avB,(Lt/Lt[0])**(-2/7),color="C3",label=label1)
        plt.plot(avB,(avB/avB[0])**(-2/7),color="C4",linestyle="-.",label=label3)
        plt.plot(avB,(CSimple/CSimple[0]),color="C1",label=label2)
        plt.xlabel(label0)
        plt.legend()
        plt.savefig("Figures/thresholdavB.png",dpi=400,bbox_inches='tight')
        plt.show()


#pick the appropriate routine

routine = "Flux Expansion"
# routine = "Connection Length"
# routine = "Averaged Field"
# routine = "Length Difference"

scanThresholds(routine)