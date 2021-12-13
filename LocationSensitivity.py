from operator import index
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
from scipy.optimize import curve_fit
import re
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
# customOrder = ["fi200E-3","fi250E-3","fi260E-3","fi262E-3","fi264E-3","fi270E-3","fi262E-3Backward","fi250E-3Backward","fi245E-3Backward","fi240E-3Backward","fi235E-3Backward","fi230E-3Backward","fi225E-3Backward","fi220E-3Backward"]

def determineC0(Spar,C):
    for k in range(len(C)):
        if Spar[k] > Spar[-1]/100:

            return k-1


def scanThresholds(routine,plotspace):
    folderList = []
    if routine=="Flaring":
        folderList = [
        "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L1_Angle0BulgeExperimentNoBulge_Lpar20\ImpurityScanqpll5E7ne1E19_Ncooling_Recomb",
        "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L1_Angle0BulgeExperimentBulge_Lpar20\ImpurityScanqpll5E7ne1E19_Ncooling_Recomb",
        # "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L1_Angle0BulgeExperimentNoBulge_Lpar20_TightGrid\ImpurityScanqpll5E7ne1.5E19_Ncooling_Recomb",
        # "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L1_Angle0BulgeExperimentBulge_Lpar20_TightGrid\ImpurityScanqpll5E7ne1.5E19_Ncooling_Recomb",
        ]
    if routine=="FlaringAll":
        folderList = [
        "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L1_Angle0BulgeExperimentNoBulge_Lpar20\ImpurityScanqpll5E7ne1E19_Ncooling_Recomb",
        "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L1_Angle0BulgeExperimentBulge_Lpar20\ImpurityScanqpll5E7ne1E19_Ncooling_Recomb",
        # "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L1_Angle0BulgeExperimentNoBulge_Lpar20_TightGrid\ImpurityScanqpll5E7ne1.5E19_Ncooling_Recomb",
        # "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L1_Angle0BulgeExperimentBulge_Lpar20_TightGrid\ImpurityScanqpll5E7ne1.5E19_Ncooling_Recomb",
        ]
    if routine == "GradB" or routine == "GradBSensitivity" or routine == "GradBAbsolute":
        folderList = [
        "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L1_Angle90\ImpurityScanqpll5E7ne1E19_Ncooling_Recomb_v0BC",
        # "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L1_Angle0\ImpurityScanqpll5E7ne1E19_NoFluxLim_noDiff",
        "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L1_Angle0\ImpurityScanqpll5E7ne1E19_Ncooling_Recomb"
        ]
    if routine == "Different C":
        folderList = [
        "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L1_Angle90\ImpurityScanqpll5E7ne1E19_Ncooling_Recomb_v0BC",
        "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L1_Angle90\powerScanne1E19fi17E-3_Ncooling",
        "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L1_Angle90\DensityScanqpll5E7fi17E-3_Ncooling"
        ]
        
    reverse = -1
    Thresholds = []
    ThresholdError = []
    Fr = []
    percImpurity = []
    Exhemptfiles = ["baserun","ref","notes.txt","INPUT.m"]
    counter = 0
    FlaringPitch = []
    FlaringDistance = []
    for folder in folderList:
        Files = os.listdir(str(folder))
        Files = natsorted(Files)

        if "powerScan" in folder:
            Files = Files[::-1]
        C = []
        Sh = []
        Sherr = []
        Bt = []
        Heats = []
        N0RAD = []
        threshFile = 0
        for File in Files:
            if File in Exhemptfiles:
                continue

            fileName = str(folder)+"/"+str(File)+"/" + str("balance.nc")
            rootgrp = Dataset(str(fileName), "r", format="NETCDF4")
            # print(rootgrp)
            Xpoint = -1
            SOLring1 = 0
            RING = 14
            if routine=="Flaring" or  routine=="FlaringAll":
                RING = 25

            quantities2d,SOLring1 = unpackSOLPS(fileName, reverse,RING)
            pos,upper,lower = SOLring1.calcFrontqpll()
            partoPol = interpolate.interp1d(SOLring1.Spar,SOLring1.Spol,kind='cubic',fill_value=0,bounds_error=False)

            if plotspace=="pol":
                Sh.append(partoPol(pos))
                Sherr.append([partoPol(pos)-partoPol(lower),partoPol(upper)-partoPol(pos)])
            else:
                Sh.append(pos)
                Sherr.append([pos-lower,upper-pos])
            # calculate the control parameter for this SOL ring
            C.append(SOLring1.determineC())
            N0RAD.append(np.sqrt(trapz(SOLring1.qf,SOLring1.Spar)/trapz(SOLring1.qf/SOLring1.ne**2,SOLring1.Spar)))


            # N0RAD.append(np.mean(SOLring1.ne[0:20]))
            Bt.append(SOLring1.B[0])
            Heats.append(SOLring1.performHeatAnalysis())
            if lower > 0 and threshFile==0:
                threshFile = prevFile
            prevFile = File
                # plt.plot(SOLring1.Spar,SOLring1.ne)
                # plt.show()
        print(threshFile)
        quantities2d,SOLring1 = unpackSOLPS(str(folder)+"/"+str(threshFile)+"/" + str("balance.nc"), reverse,RING)
        FlaringDistance.append(SOLring1.Spol)
        FlaringPitch.append(SOLring1.Bpol/SOLring1.B)


        Sh = np.array(Sh)
        C = np.array(C)
        QF = []
        QI = []
        for i in range(len(Heats)):
            QF.append(Heats[i]["impurity"])
            QI.append(Heats[i]["conducted"])

        #determine the detachment threshold of the data

        index0 = determineC0(np.array(Sh)-np.array(Sherr)[:,0],C)
        print(index0)
        Cplot = C[index0:]
        print(len(Cplot))
        Splot = Sh[index0:]
        SplotErr =Sherr[index0:]
        # print(index0)
        # for i in range(len(C)):
        #     if Sh[i] >= 1.001*Sh[index0]:
        #         Cplot.append(C[i])
        #         Splot.append(Sh[i])
        #         SplotErr.append(Sh[:,1][i])
        Cplot = np.array(Cplot)
        Splot =np.array(Splot)
        SplotErr = np.array(SplotErr)

        Thresholds.append(C[index0])
        ThresholdError.append(C[index0+1] - C[index0])
        Fr.append(np.mean(Bt))

        percImpurity.append(QF[index0]/QI[index0])
        
        if counter ==0:
            label = r"$C_{f,SOLPS}$"+" horizontal"
        else:
            label = r"$C_{f,SOLPS}$"+" vertical"
        Ct = (C[index0])

        Ctplusone = C[index0+1]
        xerr = C/Ct-C/Ctplusone
        for i in range(0,index0+1):
            Sherr[i] = [0,0]

        SHSIMPLE = 0
        if routine=="Flaring":
            SHSIMPLE =np.linspace(0.1,12,20)
            if counter==0:
                plt.errorbar(C/Ct,Sh,xerr=[xerr,xerr*0],yerr=np.transpose(np.array(Sherr)),linestyle = "",marker = "o",label="Straightdown SOLPS",color ="#53ba83")
            else:
                plt.errorbar(C/Ct,Sh,xerr=[xerr,xerr*0],yerr=np.transpose(np.array(Sherr)),linestyle = "",marker = "s",label="Flared SOLPS",color="#095169")
        if routine=="FlaringAll":
            print(Ct)
            SHSIMPLE =np.linspace(0.1,12,20)
            label = ""
            if counter ==0:
                label = "straightdown fixed wall"
            if counter ==1:
                label = "bulged fixed wall"
            if counter ==2:
                label = "straightdown tight wall"
            if counter ==3:
                label = "bulged tight wall"
            plt.plot(C/C[index0],Sh,linestyle = "",marker = "o",label=label,color =colors[counter])
            # plt.plot(N0RAD[index0+1:],Sh[index0+1:],linestyle = "",marker = "o",label=label,color =colors[counter])
            
            # plt.errorbar(C/C[index0],Sh,np.transpose(np.array(Sherr)),linestyle = "",marker = "o",label=label,color =colors[counter])

        if routine == "GradB" :
            SHSIMPLE =np.linspace(0.1,8,20)
            if counter==0:
                plt.errorbar(C/Ct,Sh,xerr=[xerr,xerr*0],yerr=np.transpose(np.array(Sherr)),linestyle = "",marker = "o",label="Horizontal SOLPS",color ="#53ba83")
            else:
                plt.errorbar(C/Ct,Sh,xerr=[xerr,xerr*0],yerr=np.transpose(np.array(Sherr)),linestyle = "",marker = "s",label="Vertical SOLPS",color="#0c0636")
        if routine == "GradBAbsolute" :
            SHSIMPLE =np.linspace(0.1,9,20)
            if counter==0:
                plt.errorbar(C,Sh,xerr=[xerr*0,xerr*0],yerr=np.transpose(np.array(Sherr)),linestyle = "",marker = "o",label="Horizontal SOLPS",color ="#53ba83")
            else:
                plt.errorbar(C,Sh,xerr=[xerr*0,xerr*0],yerr=np.transpose(np.array(Sherr)),linestyle = "",marker = "s",label="Vertical SOLPS",color="#0c0636")
       

        if routine == "Different C":
            SHSIMPLE =np.linspace(0.1,10,20)
            if counter==0:
                plt.errorbar(C/Ct,Sh,xerr=[xerr,xerr*0],yerr=np.transpose(np.array(Sherr)),linestyle = "",marker = "o",label="Impurity scan",color ="#53ba83")
            elif counter==1:
                plt.errorbar(C/Ct,Sh,xerr=[xerr,xerr*0],yerr=np.transpose(np.array(Sherr)),linestyle = "",marker = "s",label="Power scan",color="#048BA8")
            else:
                plt.errorbar(C/Ct,Sh,xerr=[xerr,xerr*0],yerr=np.transpose(np.array(Sherr)),linestyle = "",marker = "^",label="Density scan",color="#4C2C69")

        CSimple = []
        CSimple0 = []
        TotalField = SOLring1.B
        avB = []
        L = []
        FLUX = []
        for sh in SHSIMPLE:
            CSimple.append(ChIntThermalForce(SOLring1.Spar, TotalField, SOLring1.Spar[-1],SOLring1.Spar[-1] ,sh))
            CSimple0.append(ChInt(SOLring1.Spar, TotalField, SOLring1.Spar[-1],SOLring1.Spar[-1] ,sh))
            fieldinterp = interpolate.interp1d(SOLring1.Spar,TotalField)
            FLUX.append(fieldinterp(sh)/SOLring1.B[-1])
            L.append(SOLring1.Spar[-1]-sh)
            avB.append(averageB(SOLring1.Spar, TotalField, SOLring1.Spar[-1],SOLring1.Spar[-1] ,sh)**(2/7))
        FLUX = np.array(FLUX)
        L = np.array(L)
        avB = np.array(avB)
        C0 = 7**(-2/7)*(2500)**(-3/14)*(1.26E-29)**(-1/2)
        def fluxlimfunc(X,kappa0):
            [ni,T,qcond] = X
            mi = 2*1.67*10**(-27)
            return qcond/(T**(5/2)*kappa0-qcond*kappa0*T**(1/2)/(0.3*ni*np.sqrt(mi)))
        popt, pcov = curve_fit(fluxlimfunc, [SOLring1.ne,SOLring1.te,SOLring1.cond], np.gradient(SOLring1.te)/np.gradient(SOLring1.Spar),p0=[1000])
        print(popt)
        ylabel = "s"+r'$_{f,||}$'+" (m)"
        # ylabel = "ne"
        if plotspace =="pol":
            SHSIMPLE = partoPol(SHSIMPLE)
            ylabel = "s"+r'$_{f,pol}$'+" (m)"
        if routine=="Flaring":
            if counter==0:
                plt.plot((CSimple0/CSimple0[0]),(SHSIMPLE),label="Straightdown DLS",linestyle="--",color = "#F96E46",linewidth = 2)
            else:
                plt.plot((CSimple0/CSimple0[0]),(SHSIMPLE),label="Flared DLS",linestyle="-",color = "#F9C846",linewidth = 2)

        if routine == "GradB":
            if counter==0:
                plt.plot((CSimple0/CSimple0[0]),(SHSIMPLE),label="Horizontal DLS",linestyle="--",color = "#F96E46",linewidth = 2)
            else:
                plt.plot((CSimple0/CSimple0[0]),(SHSIMPLE),label="Vertical DLS",linestyle="-",color = "#F9C846",linewidth = 2)
        if routine == "GradBAbsolute":
            
            if counter==0:
                plt.plot(C0*FLUX*L**(-2/7)*avB,(SHSIMPLE),label="Horizontal DLS",linestyle="--",color = "#F96E46",linewidth = 2)
            else:

                plt.plot(C0*FLUX*L**(-2/7)*avB,(SHSIMPLE),label="Vertical DLS",linestyle="-",color = "#F9C846",linewidth = 2)

        if routine == "Different C":
            if counter==0:
                plt.plot((CSimple0/CSimple0[0]),(SHSIMPLE),label="Horizontal DLS",linestyle="--",color = "#F96E46",linewidth = 2)


        counter = counter+1
    if routine=="Flaring":
        plt.legend()
        plt.ylabel(ylabel)
        plt.xlabel(r"$C/C_{t}$")
        plt.xlim([0.98,1.73])
        plt.tight_layout()
        if plotspace =="pol":
            plt.savefig("Figures/CprofileBulge.png",dpi=400,bbox_inches='tight')
        else:
            plt.savefig("Figures/CprofileBulgePar.png",dpi=400,bbox_inches='tight')
    if routine=="FlaringAll":
        plt.legend()
        plt.ylabel(ylabel)
        plt.xlabel(r"$C$")
        # plt.xlim([0.9,1.68])
        plt.tight_layout()
        if plotspace =="pol":
            plt.savefig("Figures/CprofileBulge.png",dpi=400,bbox_inches='tight')
        else:
            plt.savefig("Figures/CprofileBulgePar.png",dpi=400,bbox_inches='tight')
    
    if routine == "GradB":
        print(percImpurity)
        ylabel = "s"+r'$_{f,pol}$'+" (m)"
        plt.ylabel(ylabel)
        plt.xlabel(r"$C/C_{t}$")
        plt.xlim([0.95,1.94])
        plt.ylim([-0.01,0.43])
        plt.legend(loc="lower right")
        plt.tight_layout()
        plt.savefig("Figures/CprofileGradB.png",dpi=400,bbox_inches='tight')
    if routine == "GradBSensitivity":
        ylabel = "s"+r'$_{f,pol}$'+" (m)"
        plt.ylabel(ylabel)
        plt.xlabel("sensitivity")
        # plt.xlim([0.95,1.7])
        # plt.ylim([-0.01,0.37])
        plt.legend(loc="upper right")
        plt.tight_layout()
        plt.savefig("Figures/CprofileGradB.png",dpi=400,bbox_inches='tight')
    if routine == "GradBAbsolute":
        ylabel = "s"+r'$_{f,pol}$'+" (m)"
        plt.ylabel(ylabel)
        plt.xlabel(r"$C$"+ " (W"+r"$^{-5/7}$"+"m"+r"$^{-11/7}$"+")")
        # plt.xlim([0.95,1.7])
        # plt.ylim([-0.01,0.37])
        plt.legend(loc="upper left")
        plt.tight_layout()
        plt.savefig("Figures/CprofileGradBAbsolute.png",dpi=400,bbox_inches='tight')
    if routine == "Different C":
        plt.legend()
        ylabel = "s"+r'$_{f,pol}$'+" (m)"
        plt.ylabel(ylabel)
        plt.xlabel(r"$C/C_{t}$")
        plt.xlim([0.95,1.94])
        plt.legend(loc="upper left")
        plt.tight_layout()
        plt.savefig("Figures/DifferentC.png",dpi=400,bbox_inches='tight')
    plt.show()
    if routine=="Flaring":
        fig = plt.figure(figsize=(3.7,2.8))
        plt.plot(FlaringPitch[0],FlaringDistance[0],color="#53ba83",linewidth=2.3,label="Straightdown")
        plt.plot(FlaringPitch[1],FlaringDistance[1],color="#095169",linewidth=2.3,label="Flared")
        # plt.plot(FlaringDistance[1],np.sqrt(SOLring1.B**2),color="#53ba83",linewidth=2.3,label="Straightdown")
        # plt.ylim([0.3,0.6])
        plt.tight_layout(pad=.0)
    
        plt.legend()
        plt.ylabel("s"+r'$_{pol}$'+" (m)")
        plt.xlabel("B"+r'$_{pol}$'+"/B")
        
        plt.savefig("Figures/FieldBulge.png",dpi=400,bbox_inches='tight')
        plt.show()

# plt.plot(np.transpose(coords))
# plt.show()
# print(datawall.read())
# while '*** 3b. Data for additional surfaces' not in tline:
#         tline = datawall.readlines(1)
#         print(tline)

routine = "Flaring"
# routine = "FlaringAll"
plotspace = "par"
# plotspace = "pol"
# routine = "GradB"
# routine = "GradBAbsolute"
# routine = "Different C"

scanThresholds(routine,plotspace)