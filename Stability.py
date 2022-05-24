import matplotlib.pyplot as plt
import numpy as np
import os
from natsort import natsorted
from scipy import interpolate
from LipschultzDLS import ChInt,ChIntThermalForce
from UnpackSOLPS import unpackSOLPS,SOLring


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
        if Spar[k] > 0:
            print(k)
            return k-1


def scanThresholds(folderList):
    Rrange = [0.5,1.7]
    Zrange = [-2.1,-1.3]
    # Rrange = [0.2,0.7]
    # Zrange = [-1.5,-1.3]
    reverse = -1
    Thresholds = []
    ThresholdError = []
    Fr = []

    Exhemptfiles = ["baserun","ref","notes.txt","INPUT.m"]
    counter = 0
    for folder in folderList:

        Files = os.listdir(str(folder))
        Files = natsorted(Files)#[::-1]
        FilesForward = [ file for file in Files if "Back" not in file ]
        FilesBackward = [ file for file in Files if "Back" in file ]
        if "L1_Angle11" in folder:
            Files = [ file for file in Files if "115" not in file ]
        if "powerScan" in folder:
            Files = Files[::-1]
        # define frame containers for making movies/gifs
        if customOrder != 0:
            Files = customOrder
        C = []
        Sh = []
        Sherr = []
        Bt = []
        Heats = []
        for File in FilesForward[4:]:
            if File in Exhemptfiles:
                continue
            fileName = str(folder)+"/"+str(File)+"/" + str("balance.nc")
            Xpoint = -1
            quantities2d,SOLring1 = unpackSOLPS(fileName, reverse,14)
            pos,upper,lower = SOLring1.calcFrontqpll()
            partoPol = interpolate.interp1d(SOLring1.Spar,SOLring1.Spol,kind='cubic',fill_value=0,bounds_error=False)
            plotspace = "pol"
            if plotspace=="pol":
                Sh.append(partoPol(pos))
                Sherr.append([partoPol(pos)-partoPol(lower),partoPol(upper)-partoPol(pos)])
            else:
                Sh.append(pos)
                Sherr.append([pos-lower,upper-pos])
            # calculate the control parameter for this SOL ring
            C.append(SOLring1.determineC())
            Bt.append(SOLring1.B[0])

        Cback = []
        Shback = []
        Sherrback = []
        for File in FilesBackward:
            if File in Exhemptfiles:
                continue
            fileName = str(folder)+"/"+str(File)+"/" + str("balance.nc")
            Xpoint = -1
            quantities2d,SOLring1 = unpackSOLPS(fileName, reverse,14)
            pos,upper,lower = SOLring1.calcFrontqpll()
            partoPol = interpolate.interp1d(SOLring1.Spar,SOLring1.Spol,kind='cubic',fill_value=0,bounds_error=False)
            plotspace = "pol"
            if plotspace=="pol":
                Shback.append(partoPol(pos))
                Sherrback.append([partoPol(pos)-partoPol(lower),partoPol(upper)-partoPol(pos)])
            else:
                Shback.append(pos)
                Sherrback.append([pos-lower,upper-pos])
            # calculate the control parameter for this SOL ring
            Cback.append(SOLring1.determineC())
            Bt.append(SOLring1.B[0])




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

        

        Thresholds.append(C[index0])
        ThresholdError.append(C[index0+1] - C[index0])
        Fr.append(np.mean(Bt))
        if counter ==0:
            label = r"$C_{f,SOLPS}$"+" horizontal"
        else:
            label = r"$C_{f,SOLPS}$"+" vertical"
        Ct = (C[index0])
        Ctplusone = C[index0+1]
        xerr = C/Ct-C/Ctplusone
        xerrback = Cback/Ct-Cback/Ctplusone
        for i in range(len(Sh)):
            if Sh[i]<0.01:
                Sherr[i] = [0,0]
        for i in range(len(Shback)):
            if Shback[i]<0.01:
                Sherrback[i] = [0,0]
        plt.errorbar(C/Ct,Sh,xerr=[xerr,xerr*0],yerr=np.transpose(np.array(Sherr)),linestyle = "",marker = "o",label="SOLPS forward",color="#53ba83")
        plt.errorbar(Cback/Ct,Shback,xerr=[xerrback,xerrback*0],yerr=np.transpose(np.array(Sherrback)),linestyle = "",marker = "^",label="SOLPS backward",color="#4C2C69")
        ylabel = "s"+r'$_{f,pol}$'+" (m)"
        plt.ylabel(ylabel)
        plt.xlabel(r"$C/C_{t}$")
        plt.legend(loc="upper left")
        plt.savefig("Figures/StabilitySOLPS.png",dpi=400)
        plt.show()
    # plt.ylabel("Ct/Ct0")
        CSimple = []
        CSimple0 = []
        SHSIMPLE =np.linspace(0.1,11,50)
        TotalField = SOLring1.B
        for sh in SHSIMPLE:
            CSimple.append(ChIntThermalForce(SOLring1.Spar, TotalField, SOLring1.Spar[-1],SOLring1.Spar[-1] ,sh))
            CSimple0.append(ChInt(SOLring1.Spar, TotalField, SOLring1.Spar[-1],SOLring1.Spar[-1] ,sh))

        # plt.plot((CSimple/CSimple[0]),SHSIMPLE)
        label = r"$C_{f,theory}$"+" horizontal"
        label = "DLS theory"

        plt.plot((CSimple0/CSimple0[0]),partoPol(SHSIMPLE),label=label,color = "#095169",linewidth = 2.5)
        plt.plot([0.972,1],[0.,0.],color =  "#095169",linewidth = 2.5,linestyle="-")
        # plt.plot([np.amin(CSimple0/CSimple0[0]),1],[0.25,0.25],color = colors[counter+2],linewidth = 2.5,linestyle="--")

        plt.text(0.99, 0.1, "unstable region",
             ha="center",
             size="medium",rotation=-30,
             bbox=dict(boxstyle="darrow", fc="w", ec="k"))

        plt.text(0.99, 0.32, "detached",
             ha="center",
             size="medium",rotation=20,
             bbox=dict(boxstyle="larrow", fc="w", ec="C1"))
        plt.text(0.975, 0.05, "attached",
             ha="center",
             size="medium",rotation=0,
             bbox=dict(boxstyle="darrow", fc="w", ec="C2"))
        plt.text(0.986, 0.05, "attached",
             ha="center",
             size="medium",rotation=0,
             bbox=dict(boxstyle="rarrow", fc="w", ec="C0"))
        plt.text(1.01, 0.45, "detached",
             ha="center",
             size="medium",rotation=15,
             bbox=dict(boxstyle="darrow", fc="w", ec="C2"))
        counter = counter+1
    plt.legend()
    ylabel = "s"+r'$_{f,pol}$'+" (m)"
    plt.ylabel(ylabel)
    plt.xlabel(r"$C/C_{t}$")
    plt.tight_layout()
    # plt.savefig("CprofileGradB.png",dpi=400)
    plt.savefig("Figures/CprofileStability.png",dpi=400,bbox_inches='tight')
    plt.show()


folderList = []

folderList.append("balFiles\\L1_Angle-10\\")
scanThresholds(folderList)