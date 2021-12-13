from operator import index
# from matplotlib.lines import _LineStyle
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
from scipy.signal import find_peaks
sys.path.append('D:\\my stuff\\PhD\\Theoretical_Detachment_Control_Scripts')
from LipschultzDLS import ChInt,ChIntThermalForce
from UnpackSOLPS import unpackSOLPS,SOLring
from AnalyticCoolingCurves import LfuncN
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
        if Spar[k] > Spar[-1]/15:
            print(k)
            return k-1

def tempMethod(Te,Spar,tempthresh):
    if np.amin(Te) >tempthresh:
        return 0
    tefunc = interpolate.interp1d(Te,Spar,kind='cubic',fill_value=0,bounds_error=False)
    # plt.plot(Spar,Te)
    # plt.plot(tefunc(5),5,marker="o",label="5ev")
    # plt.plot(tefunc(6),6,marker="o",label="6ev")
    # plt.plot(tefunc(4),4,marker="o",label="4ev")
    # plt.ylabel("Te (eV)")
    # plt.xlabel("Spar (m)")
    # plt.legend()
    # plt.tight_layout()
    # plt.savefig("tedef.png",dpi=400)
    # plt.show()
    return tefunc(tempthresh)


def condvectionMethod(cond,conv,Spar):
    ratio = cond/conv

    if np.amin(ratio) >1:
        return 0
    ratioFunc = interpolate.interp1d(ratio,Spar,kind='cubic',fill_value=0,bounds_error=False)
    convFunc = interpolate.interp1d(Spar,conv,kind='cubic',fill_value=0,bounds_error=False)
    
    # plt.plot(Spar,cond,label="conducted")
    # plt.plot(Spar,conv,label="convected")
    # plt.plot(ratioFunc(1),convFunc(ratioFunc(1)),marker="o",label="front")
    # plt.ylabel("Heat flux density (Wm^-2)")
    # plt.xlabel("Spar (m)")
    # plt.legend()
    # plt.tight_layout()
    # plt.savefig("convdef.png",dpi=400)
    # plt.show()
    return ratioFunc(1)

def impurityPeakMethod(qf,Spar):
    qf = np.abs(qf)
    
    ind = np.argmax(qf)
    plt.plot(Spar,qf)
    plt.plot(Spar[ind],qf[ind],marker="o")
    plt.ylabel("Heat flux density (Wm^-2)")
    plt.xlabel("Spar (m)")
    plt.legend()
    plt.tight_layout()
    plt.savefig("convdef.png",dpi=400)
    plt.show()
    return Spar[ind]

def PeakLossPress(press,Spar):
    # cond = interpolate.interp1d(Spar,cond,kind='cubic',fill_value=0)
    # Spar = np.linspace(np.amin(Spar),np.amax(Spar),10000)
    # cond = cond(Spar)
    loss = np.gradient(press)/np.gradient(Spar)
    ind = np.argmax(loss)
    lossFuncUpper = interpolate.interp1d(loss[ind:],Spar[ind:])
    thresh = 0.9
    peaks, _ = find_peaks(loss, width=15)
    maxLoss = np.amax(loss)

    posupper = lossFuncUpper(thresh*maxLoss)

    poslower = 0
    if loss[0] < thresh*maxLoss:
        lossFuncLower = interpolate.interp1d(loss[0:ind],Spar[0:ind])
    
        poslower = lossFuncLower(thresh*maxLoss)
    if not peaks:
        peaks = [0]
    # plt.plot(Spar,loss)
    # plt.plot(Spar[peaks],loss[peaks],marker="o",linestyle="")
    # plt.plot(Spar[ind],loss[ind],marker="o",label="peak")
    # plt.plot(posupper,thresh*maxLoss,marker="o",color="C3",label="90% peak")
    # plt.plot(poslower,thresh*maxLoss,marker="o",color="C3",)
    # plt.ylabel("Heat flux density loss (Wm^-3)")
    # plt.xlabel("Spar (m)")
    # plt.legend()
    # plt.tight_layout()
    # plt.savefig("peakLoss.png",dpi=400)
    # plt.show()
    return Spar[peaks[0]],Spar[peaks[0]],Spar[peaks[0]]

def PeakLossMethod(cond,Spar):
    # cond = interpolate.interp1d(Spar,cond,kind='cubic',fill_value=0)
    # Spar = np.linspace(np.amin(Spar),np.amax(Spar),10000)
    # cond = cond(Spar)
    loss = np.gradient(cond)/np.gradient(Spar)
    ind = np.argmax(loss)
    lossFuncUpper = interpolate.interp1d(loss[ind:],Spar[ind:])
    thresh = 0.9

    maxLoss = np.amax(loss)

    posupper = lossFuncUpper(thresh*maxLoss)

    poslower = 0
    if loss[0] < thresh*maxLoss:
        lossFuncLower = interpolate.interp1d(loss[0:ind],Spar[0:ind])
    
        poslower = lossFuncLower(thresh*maxLoss)

    # plt.plot(Spar[ind],loss[ind],marker="o",label="peak")
    # plt.plot(posupper,thresh*maxLoss,marker="o",color="C3",label="90% peak")
    # plt.plot(poslower,thresh*maxLoss,marker="o",color="C3",)
    # plt.ylabel("Heat flux density loss (Wm^-3)")
    # plt.xlabel("Spar (m)")
    # plt.legend()
    # plt.tight_layout()
    # plt.savefig("peakLoss.png",dpi=400)
    # plt.show()
    return Spar[ind],posupper,poslower


def WindowLossMethod(heat,Spar):

    for windowlength in range(1,len(Spar)):
        heatchange = []
        for i in range(len(Spar)-windowlength):
            heatchange.append(heat[i+windowlength]-heat[i])
        
        if np.amax(heatchange) >= heat[-1]/2:
            windowindex = np.argmax(heatchange)
            # plt.plot(Spar,heat,linewidth=2)
            # frontindex = np.argmax(np.gradient(heat)/np.gradient(Spar))
            # # plt.plot(Spar[windowindex+windowlength],heat[windowindex+windowlength],marker="o",label="upper lim")
            # plt.plot(Spar[frontindex],heat[frontindex],marker="o",label="peak in conducted heat loss",color="C3",linestyle="")
            # plt.fill_betweenx([-6E7,9E7],Spar[windowindex+windowlength],Spar[windowindex],alpha=0.4,color="C1",label="50% conducted heat loss region")
            # plt.ylabel(r"$\frac{q_{||}B_{X}}{B}$ (Wm$^{-2}$)")
            # # plt.ylabel("neTe")
            # plt.ylim([-0.1E7,7.5E7])
            # plt.xlim([-0.5,20.5])
            # plt.xlabel(r"s$_{||}$ (m)")
            # plt.legend(loc="lower right")
            # plt.tight_layout()
            # plt.savefig("windowLoss.png",dpi=400)
            # plt.show()
            return (Spar[windowindex]+Spar[windowindex+windowlength])/2,Spar[windowindex+windowlength],Spar[windowindex]


def pressureMethod(press,Spar):
    pressLoss = np.gradient(press)/np.gradient(Spar)

    i = np.argmax(pressLoss[5:])
    plt.plot(Spar,pressLoss)
    # plt.plot(Spar[i],pressLoss[i],marker="o")
    plt.show()
    return Spar[i],Spar[i],Spar[i]
    



def scanThresholds(routine,plotspace):

    folder = "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L1_Angle90\ImpurityScanqpll5E7ne1E19_Ncooling_Recomb_v0BC"
    # folder = "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L1_Angle0\ImpurityScanqpll5E7ne1E19_Ncooling_Recomb"
    # folder = "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L1_Angle90\DensityScanqpll5E7fi17E-3_Ncooling"
        
    reverse = -1
    Thresholds = []
    ThresholdError = []
    Fr = []
    percImpurity = []
    Exhemptfiles = ["baserun","ref","notes.txt","INPUT.m","fi115E-3"]
    counter = 0
    FlaringPitch = []
    FlaringDistance = []
    for counter in [0,1,3,2,4]:
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
            prevFile = File
            fileName = str(folder)+"/"+str(File)+"/" + str("balance.nc")
            rootgrp = Dataset(str(fileName), "r", format="NETCDF4")

            SOLring1 = 0

            quantities2d,SOLring1 = unpackSOLPS(fileName, reverse,14)
            # plt.plot(SOLring1.Spar,SOLring1.ne)
            # plt.ylim([0,0.5E20])
            # plt.show()
            if counter ==0:
                pos,upper,lower = tempMethod(SOLring1.te,SOLring1.Spar,5),tempMethod(SOLring1.te,SOLring1.Spar,6),tempMethod(SOLring1.te,SOLring1.Spar,4)
            elif counter ==1:
                pos = condvectionMethod(SOLring1.cond,SOLring1.conv,SOLring1.Spar)
                upper = pos
                lower = pos
            elif counter ==2:
                pos,upper,lower = PeakLossMethod((SOLring1.cond)/SOLring1.B,SOLring1.Spar)
            elif counter ==3:
                pos,upper,lower = PeakLossPress(SOLring1.ne*SOLring1.te,SOLring1.Spar)

            elif counter ==4:
                pos,upper,lower = WindowLossMethod((SOLring1.cond)/SOLring1.B,SOLring1.Spar)

            else:
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


            if lower > 0 and threshFile==0:
                threshFile = prevFile
                # plt.plot(SOLring1.Spar,SOLring1.ne)
                # plt.show()


        print(threshFile)
        Sh = np.array(Sh)
        C = np.array(C)
        QF = []
        QI = []
        for i in range(len(Heats)):
            QF.append(Heats[i]["impurity"])
            QI.append(Heats[i]["conducted"])

        #determine the detachment threshold of the data

        index0 = determineC0(np.array(Sh)-np.array(Sherr)[:,0],C)
        Cplot = C[index0:]
        Splot = Sh[index0:]
        SplotErr =Sherr[index0:]

        Cplot = np.array(Cplot)
        Splot =np.array(Splot)
        SplotErr = np.array(SplotErr)

        Thresholds.append(C[index0])
        ThresholdError.append(C[index0+1] - C[index0])


        if counter ==0:
            label = r"$C_{f,SOLPS}$"+" horizontal"
        else:
            label = r"$C_{f,SOLPS}$"+" vertical"
        Ct = (C[index0])
        SHSIMPLE = 0

        if routine == "GradB" :
            SHSIMPLE =np.linspace(0.1,7,20)
            if counter==0:
                plt.errorbar(C/Ct,Sh,np.transpose(np.array(Sherr)),linestyle = "",marker = "o",label="5ev point",color ="C0")
            elif counter ==1:
                plt.errorbar(C/Ct,Sh,np.transpose(np.array(Sherr)),linestyle = "",marker = "v",label="convection = conduction",color="C5")
            elif counter ==2:
                plt.errorbar(C/Ct,Sh,np.transpose(np.array(Sherr)),linestyle = "",marker = "s",label="peak in conducted heat loss",color="C2")
            elif counter ==3:
                plt.errorbar(C/Ct,Sh,np.transpose(np.array(Sherr)),linestyle = "",marker = "*",label="peak in electon static pressure loss",color="C3")

            else:
                plt.errorbar(C/Ct,Sh,np.transpose(np.array(Sherr)),linestyle = "",marker = "",label="power loss window 50%",capsize=6,color="C4")

        CSimple = []
        CSimple0 = []
        TotalField = SOLring1.B
        for sh in SHSIMPLE:
            CSimple.append(ChIntThermalForce(SOLring1.Spar, TotalField, SOLring1.Spar[-1],SOLring1.Spar[-1] ,sh))
            CSimple0.append(ChInt(SOLring1.Spar, TotalField, SOLring1.Spar[-1],SOLring1.Spar[-1] ,sh))
        ylabel = "s"+r'$_{f,||}$'+" (m)"
        
        if plotspace =="pol":
            SHSIMPLE = partoPol(SHSIMPLE)
            ylabel = "s"+r'$_{f,pol}$'+" (m)"

        if routine == "GradB":
            if counter==0:
                plt.plot((CSimple0/CSimple0[0]),(SHSIMPLE),label="Horizontal DLS",linestyle="--",color = "#F96E46",linewidth = 2)

        counter = counter+1

    plt.legend(loc="upper left")

    plt.ylabel(ylabel)
    plt.xlabel(r"$C/C_{t}$")
    plt.xlim([0.95,1.7])
    plt.tight_layout()
    plt.savefig("Figures/FrontMethods.png",dpi=400,bbox_inches='tight')


# plt.plot(np.transpose(coords))
# plt.show()
# print(datawall.read())
# while '*** 3b. Data for additional surfaces' not in tline:
#         tline = datawall.readlines(1)
#         print(tline)

# routine = "Flaring"
# routine = "FlaringAll"
# plotspace = "par"
plotspace = "pol"
routine = "GradB"
# routine = "Different C"
scanThresholds(routine,plotspace)
