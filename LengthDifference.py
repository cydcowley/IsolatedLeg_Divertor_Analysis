from math import e
from matplotlib import lines
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import os
import json
from natsort import natsorted
import matplotlib.image as m
import imageio as io
from numpy.core.fromnumeric import shape
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
sys.path.append('D:\\my stuff\\PhD\\DLS-model')
from LipschultzDLS import ChInt,ChIntThermalForce,averageB
from UnpackSOLPS import unpackSOLPS,SOLring
from scipy.optimize import curve_fit
from AnalyticCoolingCurves import LfuncN
import re
from AnalyticCoolingCurves import LfunLengFunccGauss, LfuncN,ratesAmjul,ratesAmjulCX
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
def determineC0(Spar,C):
    for k in range(len(C)):

        if Spar[k] > 0:
            return k-1

def fitfunc(L,A):
    return L**A
exponents = []

def scanThresholds(routine):
    folderList = []
    if routine=="Length Difference":
        folderList = [
            "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L0.75Angle0\ImpurityScanqpll5E7ne1E19_Ncooling_Recomb",
            "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L1_Angle0\ImpurityScanqpll5E7ne1E19_Ncooling_Recomb",
            "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L1.25Angle0",
            "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L1.50Angle0\ImpurityScanqpll5E7ne1E19_Ncooling"
            ]
    if routine=="Different Equations":
        folderList = [
            # "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L0.75Angle0\ImpurityScanqpll5E7ne1E19_Ncooling_Recomb",
            # "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L1_Angle0\ImpurityScanqpll5E7ne1E19_Ncooling_Recomb",
            # "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L1_Angle0\ImpurityScanqpll5E7ne1E19_NoFluxLim_noDiff",
            # "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L1.25Angle0\ImpurityScanqpll5E7ne1E19_Ncooling",
            # "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L1.50Angle0\ImpurityScanqpll5E7ne1E19_Ncooling",
            "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L0.75Angle0\DensityScanqpll5E7fi70E-3_Ncooling",
            "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L1_Angle0\DensityScanqpll5E7fi70E-3_Ncooling",
            "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L1.25Angle0\DensityScanqpll5E7fi70E-3_Ncooling",
            "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L1.50Angle0\DensityScanqpll5E7fi70E-3_Ncooling",
            ]
    if routine=="Length Fit":
        folderList = [
            "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L1_Angle0BulgeExperimentNoBulge_Lpar20\ImpurityScanqpll5E7ne1E19_Ncooling_Recomb",
            "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L0.75Angle0\ImpurityScanqpll5E7ne1E19_Ncooling_Recomb",
            "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L1_Angle0\ImpurityScanqpll5E7ne1E19_Ncooling_Recomb",
            "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L1.25Angle0",
            "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L1.50Angle0\ImpurityScanqpll5E7ne1E19_Ncooling"]
    if routine=="Kink":
        folderList = [
            "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\KinkAngle90Angle0\ImpurityScanqpll5E7ne1E19_Ncooling_Recomb",
            "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\KinkAngle67Angle23\ImpurityScanqpll5E7ne1E19_Ncooling_Recomb",
            "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\KinkAngle23Angle67\ImpurityScanqpll5E7ne1E19_Ncooling_Recomb",
            "D:\my stuff\PhD\SOLPSRUNS\isolatedBox\KinkAngle0Angle90\ImpurityScanqpll5E7ne1E19_Ncooling_Recomb"
]
        
    reverse = -1
    Thresholds = []
    ThresholdError = []
    averageDensity = []
    nu = []
    qpll = []
    averageimpurityPower = []
    CSimple = []
    avB = []
    Fr = []
    percImpurity = []
    Exhemptfiles = ["baserun","ref","notes.txt","INPUT.m","fi115E-3"]
    counter = 0
    ficalcs = []

    Lt = []
    avB = []
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
        nav = []
        qftot = []
        threshFile = 0
        prevFile = 0
        ring = 14
        for File in Files:

            if File in Exhemptfiles:
                continue

            fileName = str(folder)+"/"+str(File)+"/" + str("balance.nc")
            rootgrp = Dataset(str(fileName), "r", format="NETCDF4")
            # print(rootgrp)
            Xpoint = -1
            if folder =="D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L1_Angle0BulgeExperimentNoBulge_Lpar20\ImpurityScanqpll5E7ne1E19_Ncooling_Recomb":
                SOLring1 =25
            quantities2d,SOLring1 = unpackSOLPS(fileName, reverse,ring)
            pos,upper,lower = SOLring1.calcFrontqpll()

            Sh.append(pos)
            Sherr.append([pos-lower,upper-pos])
            # calculate the control parameter for this SOL ring
            C.append(SOLring1.determineC())
            Bt.append(SOLring1.B[0])
            Heats.append(SOLring1.performHeatAnalysis())
            nav.append( trapz(SOLring1.qf,SOLring1.Spar)/trapz(SOLring1.qf/SOLring1.ne,SOLring1.Spar))
            qftot.append(trapz(SOLring1.qf,SOLring1.Spar))


            if lower > 0 and threshFile==0:
                threshFile = prevFile
                print(threshFile)
                print(File)
            if lower>0:
                plt.plot(SOLring1.te,SOLring1.ne*SOLring1.te/(SOLring1.ne*SOLring1.te)[-1])

            prevFile = File            

        # print(threshFile)

        plt.show()
        Sh = np.array(Sh)
        C = np.array(C)

        QF = []
        QI = []
        for i in range(len(Heats)):
            QF.append(Heats[i]["impurity"])
            QI.append(Heats[i]["conducted"])

        #determine the detachment threshold of the data 
        
        index0 = determineC0(Sh-np.array(Sherr)[:,0],C)

        label1 = "L"+r"$_{t}^{-2/7}$"
        label2 = "L"+r"$_{t}^{-4.5/7}$"
        Leff = (SOLring1.Spar[-1]-Sh[index0:])
    
        if routine=="Length Fit":
            plt.errorbar(C[index0:]/C[index0],Sh[index0:],np.transpose(Sherr[index0:]),marker="o",linestyle="",label="SOLPS")

            plt.plot((Leff/Leff[0])**(-4.5/7),Sh[index0:],label=label2)
            plt.plot((Leff/Leff[0])**(-2/7),Sh[index0:],color="C3",label=label1)
            plt.xlabel("C/Ct")
            plt.ylabel("Sf (m)")
            plt.legend()
            plt.tight_layout()
            plt.savefig(str(counter)+".png",dpi=400,bbox_inches='tight')
            plt.show()
            popt, pcov = curve_fit(fitfunc, Leff/Leff[0], C[index0:]/C[index0])
            exponents.append(popt[0]) 
        averageDensity.append(nav[index0])
        averageimpurityPower.append(qftot[index0])
        # print(Sh-np.array(Sherr)[:,0])
        Thresholds.append(C[index0])
        ThresholdError.append((C[index0+1] - C[index0]))
        Fr.append(np.mean(Bt))
        Lt.append(SOLring1.Spar[-1])
        percImpurity.append(-1*QF[index0]*np.mean(Bt)/(50E6*0.5))
        CSimple.append(ChInt(SOLring1.Spar, SOLring1.B, SOLring1.Spar[-1],SOLring1.Spar[-1] ,0))
        avB.append(averageB(SOLring1.Spar, SOLring1.B, SOLring1.Spar[-1],SOLring1.Spar[-1] ,0))
        counter = counter+1

        fileName = str(folder)+"/"+str(threshFile)+"/" + str("balance.nc")
        Xpoint = -1

        print(fileName)
        if folder =="D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L1_Angle0BulgeExperimentNoBulge_Lpar20\ImpurityScanqpll5E7ne1E19_Ncooling_Recomb":
            ring =25
        quantities2d,SOLring1 = unpackSOLPS(fileName, reverse,ring)
        # quantities2d,SOLring2 = unpackSOLPS(fileName, reverse,ring+1)
        rootgrp = Dataset(str(fileName), "r", format="NETCDF4")
        dv = np.array(rootgrp['vol'])
        J = rootgrp["fna_tot"][1][0][ring][::reverse][1:Xpoint]
        # J = rootgrp["fne"][ring][::reverse][1:Xpoint]
        convectFlux = rootgrp["fmo_flua"][1][0][ring][::reverse][1:Xpoint]#rootgrp["fmo_flua"][0][0][ring][::reverse][1:Xpoint]#/rootgrp["vol"][ring][::reverse][1:Xpoint]
        flowvel = rootgrp["ua"][1][ring][::reverse][1:Xpoint]
        # J = np.sum(rootgrp["fna_tot"][0][1][ring][::reverse][1:Xpoint]
        # print(rootgrp["ua"])

        EIRENEPARTLOSS = (np.sum(rootgrp["eirene_mc_papl_sna_bal"],axis=0)[1]/dv)[ring][1:Xpoint][::reverse]
        EIRENEMOMLOSS = (np.sum(rootgrp["eirene_mc_mapl_smo_bal"],axis=0)[1]/dv)[ring][1:Xpoint][::reverse]

        ficalcs.append(SOLring1.returnCalculationsfi(0))
        mi = 2*1.67E-27
        e = 1.60E-19
        nu.append(SOLring1.ne[-1])
        qpll.append(SOLring1.cond[-1])
        electronPress = SOLring1.ne*SOLring1.te*e
        rampress = SOLring1.ne*mi*SOLring1.FlowVelocity**2
        totpress = SOLring1.ne*(SOLring1.te*e+SOLring1.ti*e)#+rampress
        pressurecalc = np.multiply(SOLring1.ne,0)+electronPress[-1]
        dtds = np.gradient(SOLring1.te)/np.gradient(SOLring1.Spar)
        pressurecalc = pressurecalc-0.3*cumtrapz(SOLring1.ne[::-1]*dtds[::-1],SOLring1.Spar[::-1],initial=0)
        dtotds = np.gradient(totpress[::-1])/np.gradient(SOLring1.Spar[::-1])
        convected = -1*np.gradient(1.67*10**(-27)*20*SOLring1.ne*(flowvel)**2)/np.gradient(SOLring1.Spol)
        pressurecalc = pressurecalc+0.3*cumtrapz(dtotds,SOLring1.Spar[::-1],initial=0)
        fna_nanom = rootgrp["fna_nanom"]
        ionStaticPressureGrad = np.gradient(SOLring1.ni*SOLring1.te*e)/SOLring1.Sdiff
        fnbx_nanom = np.sum(fna_nanom,axis=0)[0][ring][1:Xpoint]
        fnby_nanom = np.sum(fna_nanom,axis=0)[1][ring][1:Xpoint]
        radTrans = fnby_nanom-np.sum(fna_nanom,axis=0)[1][rootgrp["topiy"][ring][1]][rootgrp["topix"][ring][1:Xpoint]]
        v0 = 0.25*(8/3.1415)**0.5
        v0 = v0*(SOLring1.te[0]*1.60E-19)**0.5
        v0 = 40*v0*(2*1.67E-27)**(-0.5)
        v0 = SOLring1.ni[0]*SOLring1.FlowVelocity[0]/SOLring1.n0[0]

        n0calc = SOLring1.ne[0]*SOLring1.FlowVelocity[0]/v0

        dn0calc = cumtrapz(-1*SOLring1.ionisSource+SOLring1.recombSource,SOLring1.Spar,initial=0)
        J0 = cumtrapz(-0.8*SOLring1.ionisSource+SOLring1.recombSource,SOLring1.Spar,initial=0)+SOLring1.n0[0]*v0
        # plt.plot(SOLring1.te,J0)
        # plt.plot(SOLring1.te,SOLring1.ni*SOLring1.FlowVelocity)
        # plt.plot(SOLring1.te,electronPress/totpress)

        # plt.plot(SOLring1.te,SOLring1.ionisSource)
        v0calc = (dn0calc/(SOLring1.n0-SOLring1.n0[0]))
        # plt.plot(SOLring1.te[:],J)
        # plt.plot(SOLring1.te[:],electronPress)
        # print(SOLring1.FlowVelocity[0]/np.sqrt(1.60E-19*(SOLring1.ti[0])/(2*1.67E-27)))
        # plt.plot(SOLring1.te,-1*cumtrapz(EIRENEMOMLOSS[::-1],initial=0)[::-1])
        # plt.plot(SOLring1.te,cumtrapz(ionStaticPressureGrad[::-1],SOLring1.Spar[::-1],initial=0)[::-1])
        # plt.plot(SOLring1.te[:],cumtrapz((radTrans[::-1]),initial=0)[::-1])
        
        # plt.plot(SOLring1.te[100:],(SOLring1.ne*SOLring1.n0*SOLring1.FlowVelocity*ratesAmjul("D:\\my stuff\\PhD\\DLS-model\\rates/AMJUEL_H4-2.1.5.txt",SOLring1.te,SOLring1.ne)/(SOLring1.te[-1]*SOLring1.ne[-1]))[100:])
        momentumLoss = 2*1.67E-27*cumtrapz((SOLring1.FlowVelocity[::-1])*SOLring1.ne[::-1]*SOLring1.n0[::-1]*ratesAmjulCX("D:\my stuff\PhD\DLS-model\\rates/AMJUEL_H3-3.1.8.txt",SOLring1.te,SOLring1.te)[::-1],SOLring1.Spar[::-1],initial=0)[::-1]
        # momentumLoss = momentumLoss+2*1.67E-27*cumtrapz(SOLring1.FlowVelocity[::-1]*SOLring1.ne[::-1]*SOLring1.n0[::-1]*ratesAmjul("D:\my stuff\PhD\DLS-model\\rates//AMJUEL_H4-2.1.8.txt",SOLring1.te,SOLring1.ne)[::-1],SOLring1.Spar[::-1],initial=0)[::-1]
        # momentumLoss = momentumLoss-2*1.67E-27*cumtrapz(np.sqrt(3E-19/(3.E-27))*SOLring1.ne[::-1]*SOLring1.n0[::-1]*ratesAmjul("D:\my stuff\PhD\DLS-model\\rates//AMJUEL_H4-2.1.5.txt",SOLring1.te,SOLring1.ne)[::-1],SOLring1.Spar[::-1],initial=0)[::-1]
        # plt.plot(SOLring1.te,momentumLoss)
        # plt.plot(SOLring1.te,-1*cumtrapz(EIRENEMOMLOSS[::-1],SOLring1.Spar[::-1],initial=0)[::-1])

        # plt.plot(SOLring1.te[100:],(cumtrapz(SOLring1.CXmomLoss,initial=0))[100:])
        # plt.plot(SOLring1.te[100:],(n0calc+dn0calc/(v0))[100:])
        # print(SOLring1.ne[0]*SOLring1.FlowVelocity[0]/SOLring1.n0[0])
        # plt.plot(SOLring1.te,(SOLring1.ne*SOLring1.n0*ratesAmjul("D:\\my stuff\\PhD\\DLS-model\\rates/AMJUEL_H4-2.1.5.txt",SOLring1.te,SOLring1.ne)))
        
        # plt.plot(SOLring1.te,cumtrapz(SOLring1.ionisSource[::-1],-1*SOLring1.Spar[::-1],initial=0)[::-1])
        # plt.plot(SOLring1.te,SOLring1.recombSource)
        # plt.plot(SOLring1.te,(SOLring1.ne*(SOLring1.te+SOLring1.ti))/(SOLring1.ne[-1]*(SOLring1.te[-1]+SOLring1.ti[-1])))
        # plt.plot(SOLring1.te,(SOLring1.ne*SOLring1.te+1.67*10**(-27)*SOLring1.ne*SOLring1.FlowVelocity**2)/(SOLring1.ne[-1]*SOLring1.te[-1]+1.67*10**(-27)*SOLring1.ne[-1]*SOLring1.FlowVelocity[-1]**2))
        # plt.plot(SOLring1.te[130:-50],SOLring1.ne[130:-50]-SOLring2.ne[130:-50])

        # plt.plot(SOLring1.Spar,SOLring1.te/(SOLring1.te[-1]))
        # plt.plot(SOLring1.Spar,SOLring1.cond/(SOLring1.cond[-1]))
    # plt.xlim([20,80])
    # plt.ylim([0,2E23])
    plt.xlabel("Te (eV)")
    # plt.ylabel("q/qu")
    plt.ylabel("J")
    plt.savefig("pressurevar.png",dpi=400)
    plt.show()
        

    avB = np.array(avB)
    Lt = np.array(Lt)
    averageimpurityPower = np.array(averageimpurityPower)
    averageDensity = np.array(averageDensity)
    ficalcs = np.array(ficalcs)
    nu = np.array(nu)
    qpll = np.array(qpll)
    if routine=="Length Difference":
        label0 =  "C"+r"$_{t,SOLPS}$"
        for i in range(len(Thresholds)):
            if i==0:
                plt.annotate("I-"+str(i+1),(Lt[i],Thresholds[i]/Thresholds[0]-0.05))
                plt.errorbar(Lt[i],(Thresholds[i]/Thresholds[0]),ThresholdError[i]/Thresholds[0],label=label0,marker="o",color = gridcolors[i])
            else:
                plt.annotate("I-"+str(i+1),(Lt[i]-1,Thresholds[i]/Thresholds[0]))
                plt.errorbar(Lt[i],(Thresholds[i]/Thresholds[0]),ThresholdError[i]/Thresholds[0],marker="o",color = gridcolors[i])

        Lplot = np.linspace(np.amin(Lt),np.amax(Lt),100)
        plt.plot(Lplot,(Lplot/Lplot[0])**(-2/7),color="C3",label=label1)
        plt.plot(Lplot,(Lplot/Lplot[0])**(-4.5/7),color="C1",label=label2)
        # plt.plot(Lt,averageDensity[0]/averageDensity,label="average ne")
        # plt.plot(Lt,averageimpurityPower/averageimpurityPower[0],label="QF")
        plt.xlabel(label1)
        plt.legend()
        plt.savefig("LengthFit.png",dpi=400,bbox_inches='tight')
        plt.show()
        popt, pcov = curve_fit(fitfunc, Lt/Lt[0],Thresholds/Thresholds[0])
        exponents.append(popt[0]) 
    if routine=="Different Equations":
        label0 =  "C"+r"$_{t,SOLPS}$"
        for i in range(len(Thresholds)):
            if i==0:
                plt.annotate("I-"+str(i+1),(Lt[i],Thresholds[i]/Thresholds[0]-0.05))
                plt.errorbar(Lt[i],(Thresholds[i]/Thresholds[0]),[[0],[ThresholdError[i]/Thresholds[0]]],label=label0,marker="o",color = gridcolors[i])
            else:
                plt.annotate("I-"+str(i+1),(Lt[i]-1,Thresholds[i]/Thresholds[0]))
                plt.errorbar(Lt[i],(Thresholds[i]/Thresholds[0]),[[0],[ThresholdError[i]/Thresholds[0]]],marker="o",color = gridcolors[i])

        Lplot = np.linspace(np.amin(Lt),np.amax(Lt),100)
        plt.plot(Lplot,(Lplot/Lplot[0])**(-2/7),color="C3",label=label1)
        # plt.plot(Lt,averageDensity[0]/averageDensity,label="average ne")
        # plt.plot(Lt,averageimpurityPower/averageimpurityPower[0],label="QF")
        plt.plot(Lt,ficalcs[:,0]/ficalcs[:,0][0],label="non constant lengyel integral")
        # plt.plot(Lt,ficalcs[:,1]/ficalcs[:,1][0],label="FI2")
        plt.plot(Lt,ficalcs[:,2]/ficalcs[:,2][0],label="real Tu used")
        plt.plot(Lt,ficalcs[:,3]/ficalcs[:,3][0],label="integrate over distance")
        plt.plot(Lt,ficalcs[:,4]/ficalcs[:,4][0],label="non-constant pressure")
        plt.plot(Lt,ficalcs[:,5]/ficalcs[:,5][0],label="all power sinks")
        plt.xlabel("Lt (m)")
        plt.legend(bbox_to_anchor=(1, 1))
        plt.savefig("DifferentEquations.png",dpi=400,bbox_inches='tight')
        plt.show()

    if routine=="Kink":
        label0 =  "C"+r"$_{t,SOLPS}$"
        for i in range(len(Thresholds)):
            if i==0:
                plt.annotate("III-"+str(i+1),(Lt[i],Thresholds[i]/Thresholds[0]-0.05))
                plt.errorbar(Lt[i],(Thresholds[i]/Thresholds[0]),ThresholdError[i]/Thresholds[0],label=label0,marker="o",color = gridcolors[i])
            else:
                plt.annotate("III-"+str(i+1),(Lt[i]-1,Thresholds[i]/Thresholds[0]))
                plt.errorbar(Lt[i],(Thresholds[i]/Thresholds[0]),ThresholdError[i]/Thresholds[0],marker="o",color = gridcolors[i])
        label1 = "L"+r"$_{t}^{-2/7}$"
        label1 = "L"+r"$_{t}^{-4.5/7}$"
        Lplot = np.linspace(Lt[0],Lt[1],100)

        plt.plot(Lplot,(Lplot/Lplot[0])**(-2/7),color="C3",label=label1)
        plt.plot(Lplot,(Lplot/Lplot[0])**(-4.5/7),color="C1",label=label1)

        plt.xlabel(label1)
        plt.legend()
        plt.savefig("Kinkfit.png",dpi=400,bbox_inches='tight')
        plt.show()
        popt, pcov = curve_fit(fitfunc, Lt/Lt[0],Thresholds/Thresholds[0])
        exponents.append(popt[0]) 






# routines = ["Length Difference","Length Fit","Kink"]
routine = "Different Equations"
# for routine in routines:
#     scanThresholds(routine)

# exponents = np.array(exponents)

# exponents = exponents.flatten()

# print(np.mean(exponents)*7)
# print(np.std(exponents)*7)

scanThresholds(routine)