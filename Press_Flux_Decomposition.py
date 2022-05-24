
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import  cumtrapz
from UnpackSOLPS import unpackSOLPS

plt.rcParams["font.family"] = "serif"
params = {'legend.fontsize': 'medium',
         'axes.labelsize': 'medium',
         'axes.titlesize':'medium',
         'xtick.labelsize':'medium',
         'ytick.labelsize':'medium',
        #  'figure.figsize': (4,3.2),
         }
plt.rcParams.update(params)


colors = ["#356288","#fe1100","#aacfdd","#fe875d"]

def determineC0(Spar,C):
    for k in range(len(C)):
        if Spar[k] > Spar[-1]/100:
            return k-1


fileNames = ["balFiles\\L1_Angle90\\fi17E-3\\balance.nc","balFiles\\L1_Angle90\\fi50E-3\\balance.nc"]
for fileName in fileNames:
    if fileName == "balFiles\\L1_Angle90\\fi17E-3\\balance.nc":
        figname = ""
    else:
        figname = "Detached"
    rootgrp = Dataset(str(fileName), "r", format="NETCDF4")

    Xpoint = -1
    SOLring1 = 0
    RING = 14

    quantities2d,SOLring1 = unpackSOLPS(fileName, -1,RING)
    pos,upper,lower = SOLring1.calcFrontqpll()
    plt.plot(SOLring1.Spar,1.60*10**(-19)*SOLring1.ne*SOLring1.te,label="electron static",color="#3296FA",linewidth=2)
    plt.plot(SOLring1.Spar,1.60*10**(-19)*SOLring1.ne*SOLring1.ti,label="ion static",color="#53A600",linestyle="--",linewidth=2)
    plt.plot(SOLring1.Spar,(1*1.67*10**(-27))*SOLring1.ne*np.sign(SOLring1.FlowVelocity)*SOLring1.FlowVelocity**2,label="ram",color="#1C61A6",linestyle="-.",linewidth=2)
    if fileName != "balFiles\\L1_Angle90\\fi17E-3\\balance.nc":
        plt.fill_betweenx([-6E7,9E7],lower,upper,alpha=0.4,color="C1",label="detachment front region")    
        plt.plot([pos,pos],[-6E7,9E7],color="C3",label="effective front position")
    plt.ylim(0,np.amax(1.1*1.60*10**(-19)*SOLring1.ne*SOLring1.ti))
    plt.xlabel("s"+r'$_{||}$'+" (m)")
    plt.ylabel("pressure (Pa)")
    plt.legend()
    plt.savefig("Figures/otherPressures"+str(figname)+".png",dpi=400,bbox_inches='tight')
    plt.show()


    plt.plot(SOLring1.Spar,SOLring1.cond,label="electron conduction",color="#3296FA",linewidth=2)
    plt.plot(SOLring1.Spar,SOLring1.condI,label="ion conduction",color="#53A600",linestyle="--",linewidth=2)
    plt.plot(SOLring1.Spar,SOLring1.conv,label="convection",color="#1C61A6",linestyle="-.",linewidth=2)
    if fileName != "balFiles\\L1_Angle90\\fi17E-3\\balance.nc":
        plt.fill_betweenx([-6E7,9E7],lower,upper,alpha=0.4,color="C1",label="detachment front region")
        plt.plot([pos,pos],[-6E7,9E7],color="C3",label="effective front position")
    plt.ylim(0,np.amax(1.1*SOLring1.cond))
    plt.legend()
    plt.xlabel("s"+r'$_{||}$'+" (m)")
    plt.ylabel("heat flux density (Wm"+r"$^{-2}$"+")")
    plt.savefig("Figures/otherFluxes"+str(figname)+".png",dpi=400,bbox_inches='tight')
    plt.show()

    plt.plot(SOLring1.Spar,-1*cumtrapz(SOLring1.qf,SOLring1.Spar,initial=0),label="cumulative impurity losses",color="#3296FA",linewidth=2)
    plt.plot(SOLring1.Spar,-1*cumtrapz(SOLring1.ionisLoss+SOLring1.recombLoss,SOLring1.Spar,initial=0),label="cumulative deuterium losses",color="#53A600",linestyle="--",linewidth=2)
    plt.plot(SOLring1.Spar,cumtrapz(SOLring1.radTrans,SOLring1.Spar,initial=0),label="cumulative radial transport losses",color="#1C61A6",linestyle="-.",linewidth=2)
    if fileName != "balFiles\\L1_Angle90\\fi17E-3\\balance.nc":
        plt.fill_betweenx([-6E7,9E7],lower,upper,alpha=0.4,color="#FF7F00",label="detachment front region")
        plt.plot([pos,pos],[-6E7,9E7],color="C3",label="effective front position")
    plt.xlabel("s"+r'$_{||}$'+" (m)")
    plt.ylabel("line-integrated loss (Wm"+r"$^{-2}$"+")")
    plt.ylim(1.5*np.amin(cumtrapz(SOLring1.radTrans,SOLring1.Spar,initial=0)),np.amax(-1.1*cumtrapz(SOLring1.qf,SOLring1.Spar,initial=0)))
    plt.legend()
    plt.savefig("Figures/otherSinks"+str(figname)+".png",dpi=400,bbox_inches='tight')
    plt.show()
