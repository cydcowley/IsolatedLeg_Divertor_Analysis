from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import os
from natsort import natsorted
import matplotlib.image as m
import imageio as io
from scipy import interpolate
from scipy.interpolate import splrep, splev
from scipy.integrate import quad,trapz, cumtrapz, odeint, solve_ivp
from PIL import Image
import cv2
import glob
import mpltools
from mpltools import special
import sys
from scipy.interpolate.interpolate import interp1d

from scipy.linalg.decomp_svd import null_space
sys.path.append('D:\\my stuff\\PhD\\DLS-model\\')
from UnpackSOLPS import unpackSOLPS,SOLring
from AnalyticCoolingCurves import ratesAmjul,ratesAmjulCX

tdata = np.linspace(.1,100,100000)
recombData = ratesAmjul("D:\\my stuff\\PhD\\DLS-model\\rates/AMJUEL_H4-2.1.8.txt",tdata,1E19)
Rrecomb = interpolate.interp1d(tdata,recombData,fill_value=0,bounds_error=False)
tdata = np.linspace(0.1,100,10000)
ionisdata = ratesAmjul("D:\\my stuff\\PhD\\DLS-model\\rates/AMJUEL_H4-2.1.5.txt",tdata,1E19)
Rion = interpolate.interp1d(tdata,ionisdata,fill_value=0,bounds_error=False)
tdata = np.linspace(0.1,100,10000)
CXdata = ratesAmjulCX("D:\\my stuff\\PhD\\DLS-model\\rates/AMJUEL_H3-3.1.8.txt",tdata,tdata)
cxfunc = interpolate.interp1d(tdata,CXdata,fill_value=0,bounds_error=False)

plt.rcParams["font.family"] = "serif"
params = {'legend.fontsize': 'large',
         'axes.labelsize': 'large',
         'axes.titlesize':'large',
         'xtick.labelsize':'large',
         'ytick.labelsize':'large',
        #  'figure.figsize': (4,3.2),
         }
plt.rcParams.update(params)


v0 = 10*0.25*((8/np.pi)*1.60E-19*5/(2*1.67E-27))**(0.5)
        
def determineC0(Spar,C):
    index0 = 0
    gradientC = np.diff(Spar)/np.diff(C)
    gradientC2 = np.gradient(gradientC)/np.diff(C)
    index0 = np.argmax(gradientC)
    # plt.plot(C,Spar)
    # plt.plot(C[index0],Spar[index0],marker="o")
    # plt.show()
    # plt.plot(gradientC)
    # plt.plot(index0,gradientC[index0],marker="o")
    # plt.show()
    # plt.plot(C,gradientC2)
    # plt.show()
    return index0

def scanThresholds(folderList):
    Rrange = [0.5,1.7]
    Zrange = [-2.1,-1.3]
    # Rrange = [0.2,0.7]
    # Zrange = [-1.5,-1.3]
    reverse = -1
    Thresholds = []
    ThresholdError = []
    Fr = []
    percImpurity = []
    Exhemptfiles = ["baserun","notes.txt"]
    for folder in folderList:

        Files = os.listdir(str(folder))
        Files = natsorted(Files)#[::-1]
        if "L1_Angle-10" in folder:
            Files = [ file for file in Files if "Back" not in file ]
        if "L1_Angle11" in folder:
            Files = [ file for file in Files if "130" not in file ]
        if "powerScan" in folder:
            Files = Files[::-1]
        # define frame containers for making movies/gifs
        C = []
        Sh = []
        J = []
        Tt = []
        TU = []
        NU = []
        NF = []
        G = []
        P = []
        N0 = []
        flux = []
        Cd = []
        Vt = []
        NT = []
        N0UCALC = []
        QT = []
        e = 1.60E-19
        mi = 2*1.67E-27
        for File in Files:
            if File in Exhemptfiles:
                continue
            fileName = str(folder)+"/"+str(File)+"/" + str("balance.nc")
            Xpoint = -1
            quantities2d,SOLring1 = unpackSOLPS(fileName, reverse)
            # calculate front position
            # SOLring1.calcFrontTemp(5)
            SOLring1.calcFrontTemp(1)
            # SOLring1.calcFrontTemp(6)
            Sh.append(SOLring1.returnFrontPosition())
            print(SOLring1.returnFrontPosition())
            # calculate the control parameter for this SOL ring
            C.append(SOLring1.determineC())
            PU = SOLring1.te[-1]*SOLring1.ne[-1]
            J.append(SOLring1.FlowVelocity[0]*SOLring1.ne[0])
            # J.append(SOLring1.n0[0])
            Tt.append(SOLring1.te[0])
            NT.append(SOLring1.ne[0])
            G.append(SOLring1.cond[0])
            P.append(SOLring1.te[-1]*SOLring1.ne[-1])
            N0.append(SOLring1.n0[0])
            flux.append(SOLring1.ne[0]*SOLring1.FlowVelocity[0])
            Vt.append(SOLring1.FlowVelocity[0])
            QT.append(SOLring1.cond[0])#+SOLring1.FlowVelocity[0]*SOLring1.te[0]*5*e*SOLring1.ne[0])
            Cd.append(np.sqrt(2*SOLring1.te[0]*e/mi))
            N0UCALC.append((Rrecomb(SOLring1.te[-1])*SOLring1.ne[-1])/(Rion(SOLring1.te[-1])))
            NU.append(SOLring1.ne[-1])
            TU.append(SOLring1.te[-1])
            ninterp = interp1d(SOLring1.te,SOLring1.ne,kind='cubic')
            # NF.append(ninterp(10))
            # plt.plot(np.diff(2*e*SOLring1.ne*SOLring1.te+SOLring1.ne*mi*SOLring1.FlowVelocity**2)/np.diff(SOLring1.Spar))
            # print(np.amax(SOLring1.recombSource))
            # RHS = mi*(SOLring1.FlowVelocity-v0)*SOLring1.ne*SOLring1.n0*cxfunc(SOLring1.te)
            VN = SOLring1.FlowVelocity*SOLring1.ne/SOLring1.n0
            VN = VN/20
            RHS = mi*(SOLring1.FlowVelocity)*SOLring1.ne*SOLring1.n0*ratesAmjulCX("D:\\my stuff\\PhD\\DLS-model\\rates/AMJUEL_H3-3.1.8.txt",SOLring1.te,SOLring1.te)
            # RHS = RHS +mi*SOLring1.FlowVelocity*ratesAmjul("D:\\my stuff\\PhD\\DLS-model\\rates/AMJUEL_H4-2.1.8.txt",SOLring1.te,SOLring1.ne)*SOLring1.ne**2
            # RHS = RHS - Rion(SOLring1.te)*SOLring1.ne*SOLring1.n0*mi*VN
            # # plt.plot(RHS)
            # plt.show()
            print(RHS[10])

            # plt.plot(SOLring1.Spar, p(SOLring1.Spar), label="fitted")
            # plt.plot(SOLring1.Spar, SOLring1.n0/SOLring1.n0[0], label="fitted")
            # plt.plot(SOLring1.Spar, SOLring1.FlowVelocity/SOLring1.FlowVelocity[0], label="fitted")
            # plt.plot(SOLring1.Spar, SOLring1.ne*SOLring1.FlowVelocity, label="fitted")


            PKE = 0.5*SOLring1.ne*mi*SOLring1.FlowVelocity**2#+0.5*(1/SOLring1.n0)*mi*(J0)**2
            # plt.plot(SOLring1.Spar,(e*SOLring1.ne*SOLring1.te +e*SOLring1.ne*SOLring1.ti +PKE),label="SOLPS")
            # plt.plot(SOLring1.Spar,(2*e*SOLring1.ne[-1]*SOLring1.te[-1]*SOLring1.B/SOLring1.B[-1]),label="SOLPS")
            
            # plt.plot(SOLring1.Spar,(2*e*SOLring1.ne*SOLring1.te ),label="2neT")
            # plt.plot(SOLring1.Spar,np.add(-1*e*SOLring1.n0*SOLring1.te,(2*e*SOLring1.ne[-1]*SOLring1.te[-1]+e*SOLring1.n0[-1]*SOLring1.te[-1])),label="Krash calculated")
            # print("pressure loss is",SOLring1.ne[0]*SOLring1.te[0]/(SOLring1.ne[-1]*SOLring1.te[-1]))
            # plt.plot(SOLring1.Spar,e*SOLring1.ne[0]*SOLring1.te[0]+e*SOLring1.ne[0]*SOLring1.ti[0]+PKE[0]+cumtrapz(RHS,SOLring1.Spar,initial=0),label="kallenbach calculated")
            # plt.plot(SOLring1.Spar,SOLring1.CXmomLoss)
            # plt.plot(SOLring1.Spar,(SOLring1.te+SOLring1.ti)*SOLring1.ne+PKE)
            plt.plot(SOLring1.Spar,SOLring1.FlowVelocity)
            # plt.plot(SOLring1.Spar,)
            # plt.plot(SOLring1.Spar[::-1],SOLring1.ne[-1]*SOLring1.te[-1]-cumtrapz(0.71*SOLring1.ne[::-1],SOLring1.te[::-1],initial=0))
            # plt.plot(SOLring1.Spar[200:],SOLring1.ne[-1]*SOLring1.te[-1]/SOLring1.te[200:])
            # plt.plot(SOLring1.Spar[200:],SOLring1.ne[-1]*SOLring1.te[-1]/SOLring1.te[200:])
            
            # plt.plot(2*SOLring1.ne[-1]*SOLring1.te[-1]/SOLr ing1.te -SOLring1.n0*5/SOLring1.te )
            # plt.plot(2*SOLring1.ne) 
            # plt.plot(2*e*SOLring1.ne*SOLring1.te[-1],label="thermal")
            # plt.plot(SOLring1.ne*mi*SOLring1.FlowVelocity**2)
            # v0 = v0/10
            # P0 = 0.5*mi*SOLring1.n0[0]*v0**2 
            # P0 = P0+cumtrapz(mi*SOLring1.FlowVelocity*SOLring1.ne*SOLring1.n0*cxfunc(SOLring1.te),initial=0)
            # plt.plot(SOLring1.Spar,SOLring1.n0,label=str(File))
            # dn0 = SOLring1.ne**2*ratesAmjul("D:\\my stuff\\PhD\\DLS-model\\rates/AMJUEL_H4-2.1.8.txt",SOLring1.te,SOLring1.ne)/(np.sqrt((mi*0.5*SOLring1.n0[0]*v0**2-e*SOLring1.n0*SOLring1.te)/(mi*0.5*SOLring1.n0)))
            # dn0 = dn0-SOLring1.ne*SOLring1.n0*ratesAmjul("D:\\my stuff\\PhD\\DLS-model\\rates/AMJUEL_H4-2.1.5.txt",SOLring1.te,SOLring1.ne)/(np.sqrt((mi*0.5*SOLring1.n0[0]*v0**2-e*SOLring1.n0*SOLring1.te)/(mi*0.5*SOLring1.n0)))
            # plt.plot(SOLring1.te,SOLring1.n0[0]+cumtrapz(dn0,SOLring1.Spar,initial=0))
            # plt.plot(SOLring1.te,SOLring1.n0)
            # plt.plot(SOLring1.te,np.sqrt((mi*0.5*SOLring1.n0[0]*v0**2-e*SOLring1.n0*SOLring1.te)/(mi*0.5*SOLring1.n0)))
            # plt.legend()
            plt.xlabel("s (m)")
            plt.ylabel("2neT")
            plt.show()
            # plt.plot(SOLring1.Spar,SOLring1.qf/np.min(SOLring1.qf),label="impurity radiation")
            # ploss = np.gradient(SOLring1.ne*SOLring1.te)/np.gradient(SOLring1.Spar)
            # plt.plot(SOLring1.Spar,ploss/np.amax(ploss),label="thermal pressure loss")
            # plt.plot(SOLring1.Spar,SOLring1.ionisSource)
            # plt.plot(SOLring1.Spar,-1*SOLring1.recombSource)
            # plt.plot(SOLring1.Spar,SOLring1.ne*SOLring1.FlowVelocity)
            # plt.plot(SOLring1.Spar,SOLring1.te*5E21)
            # plt.plot(SOLring1.Spar,SOLring1.ne*np.sqrt(e*SOLring1.te/mi))
            # # plt.plot(SOLring1.FlowVelocity)
            # plt.legend()
            # plt.xlabel("s (m)")
            # plt.ylabel("Thermal plasma pressure")
            # plt.savefig("pressureloss.png",dpi=400)
            # plt.show()
            end = 200
            SOLring1.cond = SOLring1.cond[:end]
            SOLring1.conv = SOLring1.conv[:end]
            SOLring1.te = SOLring1.te[:end]
            SOLring1.Spar = SOLring1.Spar[:end]
            SOLring1.FlowVelocity = SOLring1.FlowVelocity[:end]
            SOLring1.ne = SOLring1.ne[:end]
            SOLring1.qf = SOLring1.qf[:end]
            SOLring1.ionisLoss = SOLring1.ionisLoss[:end]
            SOLring1.recombLoss = SOLring1.recombLoss[:end]
            HEAT = np.abs(SOLring1.qf)+np.abs(SOLring1.ionisLoss)+np.abs(SOLring1.recombLoss)
            # plt.plot(SOLring1.Spar,SOLring1.te*1E6)
            # plt.plot(SOLring1.Spar,SOLring1.cond+SOLring1.conv)
            # plt.plot(SOLring1.Spar,SOLring1.cond[0]+SOLring1.conv[0]+cumtrapz(HEAT,SOLring1.Spar,initial=0))
            # plt.plot(SOLring1.Spar,2000*SOLring1.te**(5/2)*np.gradient(SOLring1.te)/np.gradient(SOLring1.Spar))
            # plt.plot(SOLring1.Spar,SOLring1.FlowVelocity*SOLring1.te*5*e*SOLring1.ne)
            # plt.plot(SOLring1.Spar,SOLring1.cond)
            # plt.plot(SOLring1.Spar,SOLring1.te*100)
        plt.legend()
        plt.show()

        Sh = np.array(Sh)
        C = np.array(C)
        index0 = determineC0(Sh[:,0],C)
        print(index0)
        index0 = -1
        J = np.array(J[index0+1:])
        flux = flux[index0+1:]
        Sh = Sh[:,0][index0+1:]
        Tt = np.array(Tt[index0+1:])
        NT = np.array(NT[index0+1:])
        P = P[index0+1:]
        N0 = N0[index0+1:]
        N0UCALC = N0UCALC[index0+1:]
        G = G[index0+1:]
        QT = QT[index0+1:]
        NU = np.array(NU)
        TU = np.array(TU)
        # plt.plot(N0)
        # plt.plot(Tt)
        # plt.errorbar(J,Sh[:,0],Sh[:,1])
        # plt.plot(Sh,N0,marker="o",linestyle = "",color="C0",)
        # plt.plot(Sh,np.array(flux)/(v0),color="C1")

        # plt.plot(Sh,np.array(P)/np.array(Tt)-2*np.array(NT),color="C2")
        # plt.plot(Sh,(1**(7/2)-Tt**(7/2))*2000*(2/7)/QT,marker="o",linestyle = "",color="C0")
        
        plt.plot(Sh,NU*TU/10/(NU[0]*TU[0]/10))
        plt.plot(Sh,NU*TU**(1.7)/10**(1.7)/(NU[0]*TU[0]**(1.7)/10**(1.7)))
        plt.plot(Sh,NF/NF[0])
        # plt.plot(Sh,5.7*J*Tt*e,color="C1")
        # plt.plot(Sh,J,marker="o",linestyle = "")
        # plt.plot(Sh,2*NT*np.sqrt(Tt*e/mi))
        # plt.plot(Sh[:,0],Tt,marker="o",linestyle = "")
    # plt.plot(Sh,N0,marker="o",linestyle = "",color="C0",label="target n0")
    # plt.plot(Sh,np.array(flux)/(v0),color="C1",label="fast kallenbach neutral prediction")
    # plt.plot(Sh,np.array(P)/np.array(Tt)-2*np.array(NT),color="C2",label="krash neutral prediction")
    # plt.plot(Sh,QT,marker="o",label = "SOLPS conv+cond",linestyle = "",color="C0")
    # plt.plot(Sh,5.7*J*Tt*e,label="sheath transmittion gamma = 5.7",color="C1")
    plt.legend()
    plt.xlabel("Front position (m)")
    plt.ylabel("n0 (m^-3)")
    plt.ylabel("qt")
    plt.tight_layout()
    plt.savefig("neutralpred.png",dpi=500)
    plt.show()


folderList = []
# folderList.append("D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L1_Angle90\ImpurityScanqpll5E7ne1E19_Ncooling_Recomb_v0BC")
# folderList.append("D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L1_Angle30\ImpurityScanqpll5E7ne1E19_Ncooling_Recomb")
# folderList.append("D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L1_Angle11\ImpurityScanqpll5E7ne1E19_Ncooling_Recomb")
# folderList.append("D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L1_Angle90\powerScanne1E19fi17E-3_Ncooling")
folderList.append("D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L1_Angle0\ImpurityScanqpll5E7ne1E19_Ncooling_Recomb")
# folderList.append("D:\my stuff\PhD\SOLPSRUNS\isolatedBox\L1_Angle-10\ImpurityScanqpll5E7ne1E19_Ncooling_Recomb")
scanThresholds(folderList)