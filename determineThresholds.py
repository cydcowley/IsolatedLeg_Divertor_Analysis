
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import os
from natsort import natsorted
from scipy.integrate import trapz
from LipschultzDLS import ChInt,averageB
from UnpackSOLPS import unpackSOLPS
from scipy.optimize import curve_fit

# set plt figure parameters
plt.rcParams["font.family"] = "serif"
params = {'legend.fontsize': 'medium',
         'axes.labelsize': 'medium',
         'axes.titlesize':'medium',
         'xtick.labelsize':'medium',
         'ytick.labelsize':'medium',

         }
plt.rcParams.update(params)


colors = ["#356288","#fe1100","#aacfdd","#fe875d"]
gridcolors = ["#53ba83","#059b9a","#095169","#0c0636","000000"]


def determineC0(Spar,C,L):
    for k in range(len(C)):
        if Spar[k] > 0:
            return k-1



def scanThresholds(routine,DLSAnalysis):
    folderList = []
    #define folder list for each routine
    if routine=="Flux Expansion":
        folderList = [
        "balFiles\\L1_Angle90",
        "balFiles\L1_Angle30",
        "balFiles\L1_Angle11",
        "balFiles\L1_Angle0",
        ]
    if routine=="Connection Length" or routine=="Length Difference":
        folderList = [

        "balFiles\L0.75Angle0",
        "balFiles\L1_Angle0",
        "balFiles\L1.25Angle0",
        "balFiles\L1.50Angle0",   
        ]
    if routine=="Averaged Field":
        folderList = [
            "balFiles\KinkAngle90Angle0",
            "balFiles\KinkAngle67Angle23",
            "balFiles\KinkAngle23Angle67",
            "balFiles\KinkAngle0Angle90"
        ]


    reverse = -1
    counter = 0
    #create empty arrays to store information about each box, such as detachment threshold
    Thresholds = []
    ThresholdError = []
    averageDensity = []
    ThresholdsCalc = []
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
        Ccalcs = []
        Sh = []
        Sherr = []
        Bt = []
        Bx = []


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

            #calculate control parameter using additional models
            Ccalcs.append(SOLring1.returnCalculationsC(0))
            #calculate Bx and Bt
            Bt.append(SOLring1.B[0])
            Bx.append(SOLring1.B[-1])

        Sh = np.array(Sh)
        C = np.array(C)

        #determine the detachment threshold of the data

        L = SOLring1.Spar[-1]
        index0 = determineC0(Sh-np.array(Sherr)[:,0],C,L)

        # print(Sh-np.array(Sherr)[:,0])
        Thresholds.append((C[index0]))
        ThresholdsCalc.append(Ccalcs[index0])

        ThresholdError.append(C[index0+1] - C[index0])
        
        Fr.append(np.mean(Bt)/np.mean(Bx))
        Lt.append(SOLring1.Spar[-1])
        CSimple.append(ChInt(SOLring1.Spar, SOLring1.B, SOLring1.Spar[-1],SOLring1.Spar[-1] ,0))
        avB.append(averageB(SOLring1.Spar, SOLring1.B, SOLring1.Spar[-1],SOLring1.Spar[-1] ,0))
        counter = counter+1

        

    avB = np.array(avB)
    Lt = np.array(Lt)
    ThresholdsCalc = np.array(ThresholdsCalc)
    #plotting specifics for box rotation routine
    if routine=="Flux Expansion":
        label0 =  "C"+r"$_{t,SOLPS}$"
        for i in range(len(Thresholds)):
            if i==0:
                plt.annotate("II-"+str(i+1),(Fr[i],Thresholds[i]/Thresholds[0]+0.08))
                plt.errorbar((Fr[i]),(Thresholds[i]/Thresholds[0]),[[0],[ThresholdError[i]/Thresholds[0]]],label=label0,marker="o",color = gridcolors[i])
            else:
                plt.annotate("II-"+str(i+1),(Fr[i]-0.04,Thresholds[i]/Thresholds[0]))
                plt.errorbar((Fr[i]),(Thresholds[i]/Thresholds[0]),[[Thresholds[i]/Thresholds[0]-Thresholds[i]/(Thresholds[0]+ThresholdError[0])],[ThresholdError[i]/Thresholds[0]]],marker="o",color = gridcolors[i])
        label1 = "F"+r"$_{R}^{-1}$"
        label2 = "C"+r"$_{t, DLS}$"
        label3 = r"$\left< B \right>_{t}^{-2/7}$"
        
        if DLSAnalysis:
            plt.plot((Fr),(CSimple/CSimple[0]),color="C1",label=label2)
            plt.plot(Fr,ThresholdsCalc[:,1]/ThresholdsCalc[:,1][0],label="true front width",color="#ED6868",linestyle="dotted")
            plt.plot(Fr,ThresholdsCalc[:,2]/ThresholdsCalc[:,2][0],label="true pressure variation",color="#C72810",linestyle="dashed")
            plt.plot(Fr,Thresholds/Thresholds[0],label="true heat sources",color="#661515",linestyle="dashdot")
        else:
            plt.plot((Fr),(Fr/Fr[0]),color="C0",linestyle="--",label=label1)
            plt.plot((Fr),(CSimple/CSimple[0]),color="C1",label=label2)
            plt.plot((Fr),(avB/avB[0])**(-2/7),color="C4",linestyle="-.",label=label3)
        

        plt.xlabel(label1)
        plt.legend()
        plt.tight_layout()
        if DLSAnalysis:
            plt.savefig("Figures/thresholdExpansionFurther.png",dpi=400,bbox_inches='tight')

        else:
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
        Lplot = np.linspace(np.amin(Lt),np.amax(Lt),100)
        if DLSAnalysis:
            plt.plot(Lplot,(Lplot/Lplot[0])**(-2/7),color="C1",label="C"+r"$_{t, DLS}$")
        else:

            plt.plot(Lplot,(Lplot/Lplot[0])**(-2/7),color="C3",label=label1)
        if DLSAnalysis:
            plt.plot(Lt,ThresholdsCalc[:,1]/ThresholdsCalc[:,1][0],label="SOLPS T"+r'$_{e}$'+" profile",color="#ED6868",linestyle="dotted")
            plt.plot(Lt,ThresholdsCalc[:,2]/ThresholdsCalc[:,2][0],label="SOLPS pressure variation",color="#C72810",linestyle="dashed")
            plt.plot(Lt,Thresholds/Thresholds[0],label="SOLPS heat sources/sinks",color="#661515",linestyle="dashdot")
            
        label1 = "L"+r"$_{||}^{-2/7}$"
        label2 = "L"+r"$_{||}$"+" (m)"
        
        def testfunc(x,a):
            return x**a

        popt, pcov = curve_fit(testfunc, Lt/Lt[0], Thresholds/Thresholds[0],p0=[-5/7])
        label3 = "L"+r"$_{||}^{"+str(int(popt[0]*70)/10)+r"/7}$"
        print(popt)
        if not DLSAnalysis:
            plt.plot(Lplot, testfunc(Lplot/Lplot[0], *popt),label=label3,color="C6",linestyle="--")
        plt.xlabel(label2)
        plt.legend()
        plt.tight_layout()
        if DLSAnalysis:
            plt.savefig("Figures/thresholdLengthFurther.png",dpi=400,bbox_inches='tight')
        else:
            plt.savefig("Figures/thresholdLength.jpg",dpi=400,bbox_inches='tight')
        plt.show()

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
        plt.plot(avB,(CSimple/CSimple[0]),color="C1",label=label2)
        if DLSAnalysis:
            plt.plot(avB,ThresholdsCalc[:,1]/ThresholdsCalc[:,1][0],label="true front width",color="#ED6868",linestyle="dotted")
            plt.plot(avB,ThresholdsCalc[:,2]/ThresholdsCalc[:,2][0],label="true pressure variation",color="#C72810",linestyle="dashed")
            plt.plot(avB,Thresholds/Thresholds[0],label="true heat sources",color="#661515",linestyle="dashdot")

        else:
            plt.plot(avB,(Lt/Lt[0])**(-2/7),color="C3",label=label1)
            plt.plot(avB,(avB/avB[0])**(-2/7),color="C4",linestyle="-.",label=label3)


        plt.xlabel(label0)
        plt.legend()
        if DLSAnalysis:
            plt.legend(loc=1)
            plt.savefig("Figures/thresholdavBFurther.png",dpi=400,bbox_inches='tight')
        else:
            plt.savefig("Figures/thresholdavB.png",dpi=400,bbox_inches='tight')
        plt.show()


#pick the appropriate routine

# routine = "Flux Expansion"
routine = "Connection Length"
# routine = "Averaged Field"
# routine = "Length Difference"
DLSAnalysis = 1
scanThresholds(routine,DLSAnalysis)