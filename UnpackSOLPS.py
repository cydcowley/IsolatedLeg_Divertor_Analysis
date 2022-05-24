from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import os
import matplotlib.image as m
import imageio as io
from scipy import interpolate
from PIL import Image
import glob
from scipy.integrate import trapz

import sys
sys.path.append('d:\\my stuff\\PhD\\DLS-model\\')
from AnalyticCoolingCurves import LfuncN
Exhemptfiles = ["baserun","notes.txt","fi75E-3","fi280E-3norecomb","fi280E-3"]

params = {'legend.fontsize': 'medium',
         'axes.labelsize': 'medium',
         'axes.titlesize':'medium',
         'xtick.labelsize':'medium',
         'ytick.labelsize':'medium',
        #  'figure.figsize': (4,3.2),
         }
plt.rcParams.update(params)

def twoPointDensity(fr):
    return (fr**2)*(np.log(fr)/(fr-1))**(6/7)

def twoPointTemperature(fr):
    return (fr**(-2))*(np.log(fr)/(fr-1))**(-4/7)   


def determineC0(Spar,C):
    index0 = 0
    gradientC = np.gradient(Spar)/np.gradient(C)
    index0 = np.argmax(gradientC)-1
    return index0

class SOLring:
    def __init__(self,Spar,Spol,Sdiff,ParArea,te,ti,cond,condI,ne,ni,n0,B,Bpol,qf,n0rad,ionisLoss,recombLoss,molecularloss,divuloss,conv,radTrans,FlowVelocity,recombSource,ionisSource,CXmomLoss):
        self.Spar = Spar # parallel distance, [m]
        self.Spol = Spol # poloidal distance, [m]
        self.Sdiff = Sdiff # change in parallel distance betwen adjcent grids, [m]
        self.ParArea = ParArea # parallel area
        self.te = te # electron temperature, [eV]
        self.ti = ti #ion temperature [eV]
        self.cond = cond # conducted electron heat flux density, [Wm^{-2}]
        self.condI = condI # conducted ion heat flux density, [Wm^{-2}]
        self.ne = ne # electron density, [m^{-3}]
        self.ni = ni # ion density, [m^{-3}]
        self.n0 = n0
        self.B = B # total field strength, [T]
        self.Bpol= Bpol #poloidal field strength [T]
        self.qf = qf # volumetric impurity radiation, [Wm^{-3}]
        self.n0rad = n0rad # neutral radiation, [Wm^{-3}]
        self.ionisLoss = ionisLoss # ionisation volumetric heat flux loss [Wm^{-3}]
        self.recombLoss = recombLoss # recombination heat flux loss
        self.molecularloss = molecularloss # molecular heat flux loss
        self.divuloss = divuloss # heat flux loss due to grad v
        self.conv = conv # convective heat flux, wm^-2
        self.radTrans = radTrans # heat flux loss due to radial transport
        self.FlowVelocity = FlowVelocity 
        self.recombSource = recombSource
        self.ionisSource = ionisSource
        self.frontPositions = []
        self.heatsTotal = 0
        self.coolingCurve = interpolate.interp1d(self.te,self.qf/self.ne**2, kind='cubic')
        self.volumetricSinks = -1*self.ionisLoss-self.recombLoss-self.molecularloss-self.radTrans#+self.divuloss
        self.CXmomLoss = CXmomLoss

    def determineC(self):
        """determines the detachment control parameter of a SOL ring"""
        # first determine fi by dividing the integral of the real cooling curve by that of the analytical cooling curve
        fi = []
        Lcurve0 = []
        Lcurve1 = []
        for t in range(len(self.te)):
            Lcurve0.append(LfuncN(self.te[t]))
            Lcurve1.append(self.qf[t]/self.ne[t]**2)
        # Lcurve0 = interpolate.interp1d(self.te, Lcurve0, kind='cubic')
        # Lcurve1 = interpolate.interp1d(self.te, Lcurve1, kind='cubic')
        int0 = trapz(Lcurve0,self.te)
        int1 = trapz(self.qf/self.ne**2,self.te)
        # int0 = trapz(LfuncN,self.te)
        # int1 = trapz(self.qf/self.ne**2,self.te)
        fi = np.abs(int1/int0)
        C = np.sqrt(np.abs(fi))
        # add the effects of nu and qpllu into C
        # print("fi is",fi)
        C = C*self.ne[-1]
        C = C/(self.cond[-1])**(5/7)
        return C

    def plotConv(self):
        term1 = -1*5*self.ne*self.te*1.60*10**(-19)*self.FlowVelocity
        term2 = -1*0.5*1.67*10**(-27)*self.ne*self.FlowVelocity**3
        Static = (2*self.ne*self.te*1.60*10**(-19))
        Dynamic = (self.ne*2*1.67*10**(-27)*self.FlowVelocity**2)
        # plt.legend()
        
        # plt.show()
    
    def performHeatAnalysis(self):
        Conducted = trapz(self.B,self.cond/self.B)
        Impurity = trapz(self.qf,self.Spar)
        Ionisation = trapz(self.ionisLoss,self.Spar)
        Recombination = trapz(self.recombLoss,self.Spar)
        Molecular = trapz(self.molecularloss,self.Spar)
        Convected = trapz(self.B,self.conv/self.B)
        # plt.plot(self.Spar,self.ne[-1]*self.te[-1]-self.ne*self.te)


        self.heatsTotal = {
            "conducted": Conducted,
            "convected": Convected,
            "impurity": Impurity,
            "ionisation": Ionisation,
            "recombination": Recombination,
            "molecular": Molecular,
            "J": self.FlowVelocity[0]*self.ne[0]

        }
        return self.heatsTotal

    def calcFrontTemp(self,TempPoint):
        if np.amin(self.te) < 5:
            sfunc = interpolate.interp1d(self.te, self.Spar, kind='cubic',bounds_error=False,fill_value=0)
            frontPos = sfunc(TempPoint)
        else:
            frontPos = 0
        self.frontPositions.append(frontPos)
    
    def calcFrontqpll(self):
        stest = self.Spar

        perc = 0.9
        gradq =(np.gradient(self.cond/self.B)/np.gradient(self.Spar))
        frontIndex = np.argmax(gradq)
        # plt.plot(self.Spar,self.te)
        # plt.show()
        frontPos = stest[frontIndex]
        upperfunc = interpolate.interp1d(gradq[frontIndex:],stest[frontIndex:])
        upperLim = upperfunc(np.amax(gradq)*perc)
        lowerLim = 0

        if frontIndex >0:
                if np.amin(gradq[:frontIndex+1]) < np.amax(gradq)*perc:
                    lowerFunc = interpolate.interp1d(gradq[:frontIndex+1],stest[:frontIndex+1])
                    lowerLim = lowerFunc(np.amax(gradq)*perc)

        return frontPos, upperLim, lowerLim
    def calcFrontqpll(self):
        stest = self.Spar
        gradqoverB = np.gradient(self.cond/self.B)/np.gradient(stest)
        qoverB =interpolate.interp1d(stest,self.cond/self.B,fill_value="extrapolate",bounds_error=False)
        diff = 0.01

        for windowlength in np.linspace(diff,np.amax(stest),int((np.amax(stest)-diff)/diff)):
            heatchange = []
            s0s = np.linspace(0,np.amax(stest)-windowlength)
            for s0 in s0s:
                heatchange.append(qoverB(s0+windowlength)-qoverB(s0))

            if np.amax(heatchange) >= qoverB(np.amax(stest))/2:
                windowindex = np.argmax(heatchange)
                return stest[np.argmax(gradqoverB)],s0s[windowindex]+windowlength,s0s[windowindex]
                    


    def calcFrontqf(self,frac):
        qfint = cumtrapz(self.qf,self.Spar,initial=0)
        # plt.plot(self.Spar,self.cond)
        qfintInterp = interpolate.interp1d(qfint,self.Spar)
        # plt.plot(qfintInterp(qfint),qfint)
        frontPos = qfintInterp(qfint[-1]*frac)
        # plt.plot(frontPos,qfint[-1]*frac,marker="o")
        # plt.show()
        self.frontPositions.append(frontPos)
    
    def returnFrontPosition(self):
        return np.mean(self.frontPositions), np.amax(self.frontPositions)-np.amin(self.frontPositions)
    def returnCalculationsfi(self,frontpos):
        #0 is most theoretical
        Tinterp = interpolate.interp1d(self.te,self.Spar,kind='cubic')
        if np.amin(self.te) >3:
            coldpos = 0
        else: 
            coldpos = Tinterp(3)
        Coldindex = (np.abs(self.Spar-coldpos)).argmin()

        tb = trapz(self.B[Coldindex:]*self.cond[Coldindex:],self.cond[Coldindex:]/self.B[Coldindex:])
        tb = tb - trapz(self.cond[Coldindex:],self.cond[Coldindex:])
        Lfunc = []
        for t in self.te:
            Lfunc.append(LfuncN(t))

        kappa = np.polyfit(self.te[1:]**(5/2)*np.diff(self.te)/np.diff(self.Spar), self.cond[1:], 1)[0]
        kappa0 = 2600
        # plt.plot(kappa0*self.te[1:]**(5/2)*np.diff(self.te)/np.diff(self.Spar),self.cond[1:],color="C0")
        kfluxlim = kappa0*(self.te[1:]**(5/2))
        kfluxlim = kfluxlim/(1+(kappa0*self.te[1:]**(1/2)*np.diff(self.te)/np.diff(self.Spar))/(6*self.ne[1:]*np.sqrt(9E-31)))
        
        # plt.plot(kfluxlim*np.diff(self.te)/np.diff(self.Spar),self.cond[1:],color="C1")
        # plt.show()

        # kappa = 1000
        tother =  trapz(kappa*self.te[Coldindex:]**(5/2)*self.volumetricSinks[Coldindex:],self.te[Coldindex:])+trapz((self.conv[Coldindex:])*self.B[Coldindex:],(self.conv[Coldindex:])/self.B[Coldindex:]) 
        FieldInterp = interpolate.interp1d(self.Spar, self.B,kind='cubic',bounds_error=False, fill_value=(self.B[0],self.B[-1]))
        Sfront = np.linspace(frontpos,self.Spar[-1],1000)
        Tucalc0 = self.cond[-1]*7/(2*kappa0)
        Tucalc0 = Tucalc0*trapz(FieldInterp(Sfront)/self.B[-1],Sfront)
        Tucalc0 = Tucalc0**(2/7)
        ionisSource = interpolate.interp1d(self.Spar, self.ionisSource,kind='cubic',bounds_error=False, fill_value=(self.B[0],self.B[-1]))
        # plt.plot(self.Spar,self.cond)
        kappatest = self.cond/(np.gradient(self.te)/self.Sdiff)
        e = 1.60E-19
        mi = 2*1.67E-27

        # plt.ylim([0,100000])
        # plt.plot(self.Spar,self.conv+self.cond)
        vt = self.FlowVelocity[0]
        qtpred = 6*1.60*10**(-19)*self.ne[0]*self.te[0]*vt*0.5-1.60*10**(-19)*self.te[0]*self.n0[0]*vt/0.1
        NTfunc = np.poly1d(np.polyfit(self.Spar,self.n0*self.te,deg=5))
        # plt.plot(self.Spar,(self.n0*self.te))
        # plt.plot(self.Spar,np.gradient(NTfunc(self.Spar))/self.Sdiff)
        # plt.plot(self.Spar,np.gradient(NTfunc(self.Spar)))
        # plt.show()


        # calculation 0 is the most simple
        Calc0 = (self.cond[-1]**2/2)*(FieldInterp(frontpos)/self.B[-1])**2
        
        int0 = kappa0*Tucalc0**2*self.ne[-1]**2*trapz(Lfunc[Coldindex:]*self.te[Coldindex:]**(1/2),self.te[Coldindex:])
        
        Calc0 = Calc0/int0
        Calc0 = np.sqrt(Calc0)
        print(Calc0)
        Calc0 = Calc0*self.ne[-1]
        print(Calc0)
        Calc0 = Calc0/(self.cond[-1])**(5/7)
        #calculation 1 adds real tb
        Calc1 = (self.cond[-1]**2/2)+tb
        int1 =kappa0*Tucalc0**2*self.ne[-1]**2*trapz(Lfunc[Coldindex:]*self.te[Coldindex:]**(1/2),self.te[Coldindex:])
        Calc1 = Calc1/int0
        Calc1 = np.sqrt(Calc1)
        Calc1 = Calc1*self.ne[-1]
        Calc1 = Calc1/(self.cond[-1])**(5/7)
        #Calculation 2 uses the real Tu
        Calc2 = (self.cond[-1]**2/2)+tb
        int2 = kappa0*self.te[-1]**2*self.ne[-1]**2*trapz(Lfunc[Coldindex:]*self.te[Coldindex:]**(1/2),self.te[Coldindex:])
        Calc2 = Calc2/int2
        Calc2 = np.sqrt(Calc2)
        Calc2 = Calc2*self.ne[-1]
        Calc2 = Calc2/(self.cond[-1])**(5/7)
        # Calculation 3 uses correct integral and target heat
        Calc3 = (self.cond[-1]**2/2)+tb
        int3 = self.ne[-1]**2*self.te[-1]**2*trapz(self.cond[Coldindex:]*Lfunc[Coldindex:]*self.te[Coldindex:]**(-2),self.Spar[Coldindex:])
        Calc3 = Calc3/int3
        Calc3 = np.sqrt(Calc3)
        Calc3 = Calc3*self.ne[-1]
        Calc3 = Calc3/(self.cond[-1])**(5/7)
        # Calculation 4 adds q at the virtual target
        Calc4 = (self.cond[-1]**2/2)+tb
        int4 = trapz(self.cond[Coldindex:]*Lfunc[Coldindex:]*self.ne[Coldindex:]**2,self.Spar[Coldindex:])
        Calc4 = Calc4/int4
        Calc4 = np.sqrt(Calc4)
        Calc4 = Calc4*self.ne[-1]
        Calc4 = Calc4/(self.cond[-1])**(5/7)

        Calc5 = trapz(self.qf/self.ne**2,self.te)#(self.cond[-1]**2/2-self.cond[Coldindex]**2/2)+tb-tother
        int5 = trapz(Lfunc,self.te)
        Calc5 = np.abs(Calc5/int5)
        Calc5 = np.sqrt(Calc5)
        Calc5 = Calc5*self.ne[-1]
        Calc5 = Calc5/(self.cond[-1])**(5/7)
        
        return([Calc0,Calc1,Calc2,Calc3,Calc4,Calc5])

    def returnCalculationsC(self,frontpos):
        #0 is most theoretical
        Tinterp = interpolate.interp1d(self.te,self.Spar,kind='cubic')
        Lfunc = []
        for t in self.te:
            Lfunc.append(LfuncN(t))

        kappa = np.polyfit(self.te[1:]**(5/2)*np.diff(self.te)/np.diff(self.Spar), self.cond[1:], 1)[0]
        kappa0 = 2600
        # plt.plot(kappa0*self.te[1:]**(5/2)*np.diff(self.te)/np.diff(self.Spar),self.cond[1:],color="C0")

        # kappa = 1000
        FieldInterp = interpolate.interp1d(self.Spar, self.B,kind='cubic',bounds_error=False, fill_value=(self.B[0],self.B[-1]))

        Coldindex = 0
        Calc0 = (self.cond[-1]**2/2)*(FieldInterp(frontpos)/self.B[-1])**2
        
        int0 = kappa0*self.te[-1]**2*self.ne[-1]**2*trapz(Lfunc[Coldindex:]*self.te[Coldindex:]**(1/2),self.te[Coldindex:])
        
        Calc0 = Calc0/int0
        fi = Calc0
        sqrtfi = np.sqrt(Calc0)
        Calc0 = np.sqrt(Calc0)
        
        Calc0 = Calc0*self.ne[-1]
        Calc0 = Calc0/(self.cond[-1])**(5/7)
        #calculation 1 adds real tb
        densityTF = self.te[-1]**2
        densityTF = densityTF*self.ne[-1]**2
        densityTF = np.sqrt(densityTF/(self.te**2))

        int1 =self.B[-1]**2*trapz(2*self.cond[Coldindex:]*densityTF[Coldindex:]**2*Lfunc[Coldindex:]*fi/self.B[Coldindex:]**2,self.Spar[Coldindex:])
        q1 = np.sqrt(int1)

        Calc1 = sqrtfi*self.ne[-1]*(q1)**(-5/7)
        #Calculation 2 uses the real Tu

        int2 =self.B[-1]**2*trapz(2*self.cond[Coldindex:]*self.ne[Coldindex:]**2*Lfunc[Coldindex:]*fi/self.B[Coldindex:]**2,self.Spar[Coldindex:])
        q2 = np.sqrt(int2)
        Calc2 = sqrtfi*self.ne[-1]*(int2)**(-5/7)
        # Calculation 3 uses correct integral and target heat


        return([Calc0,Calc1,Calc2])

def unpack2d(rootgrp):
    """this function unpacks the raw balance file into a dictionary of 2d arrays containing important quantities"""
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

    quantities2d["sdiffpol"] = np.array(hx)
    # Parallel area:
    quantities2d["Area"] = np.array(dv)/s
    #Grid volume
    quantities2d["V"] = dv
    #specific flux ring to focus on
    # print(rootgrp['jsep'][0])
    sep = rootgrp['jsep'][0]
    ring = sep+5
    quantities2d["ring"] = ring
    #electron density (m^{-3})
    quantities2d["ne"] = rootgrp["ne"]
    #ion density (m^{-3})
    quantities2d["ni"] = rootgrp["na"][1]
    #Conductive electron heat flux (Wm^{-2}):
    fhe_cond = rootgrp['fhe_cond'][0]/quantities2d["Area"]
    quantities2d["cond"] = fhe_cond
    #Conductive ion heat flux (Wm^{-2}):
    fhi_cond = rootgrp['fhi_cond'][0]/quantities2d["Area"]
    quantities2d["condI"] = fhi_cond
    #electron temperature (eV):
    te = np.array(rootgrp["te"])
    quantities2d["te"] = te/(1.60*10**(-19))
    #ion temperature (eV)
    ti = np.array(rootgrp["ti"])
    quantities2d["ti"] = ti/(1.60*10**(-19))
    #artificial impurity radiation (W):
    #imprad = np.sum(rootgrp["b2stel_she_bal"][2:],axis=0)
    imprad = rootgrp["b2stel_she_bal"][1]
    quantities2d["imprad"] = imprad
    #flow velocity
    vfluid = rootgrp["ua"][1]
    quantities2d["vfluid"] = vfluid
    #recombination particle source 
    quantities2d["recombSource"] = np.sum(rootgrp['eirene_mc_pppl_sna_bal'],axis=0)[1]/dv
    # RecombRate = quantities2d["recombSource"]/(np.array(quantities2d["ne"])**2)
    # trueRate = ratesAmjul("D:\my stuff\PhD\DLS-model\\rates/AMJUEL_H4-2.1.8.txt",quantities2d["te"].flatten(),np.array(quantities2d["ne"]).flatten())
    # plt.plot(quantities2d["te"].flatten(),RecombRate.flatten(),linestyle="",marker = "o")
    # plt.show()
    #ionisation particle source
    quantities2d["ionisSource"] = np.sum(rootgrp['eirene_mc_papl_sna_bal'],axis=0)[1]/dv
    #convected electron heat flux (Wm^{-2})
    fhe_conv = (te*np.array(rootgrp["ne"])*5 +0.5*1.67*10**(-27)*np.array(rootgrp["ne"])*vfluid**2)*vfluid
    quantities2d["conv"] = fhe_conv#fhe_conv
    #neutral density
    nDatom = rootgrp['dab2'][0][:,:]
    quantities2d["n0"] = nDatom
    #neutral radiation heat source/sink
    quantities2d["neutralrad"] = np.sum(rootgrp["b2stel_she_bal"][0:2],axis=0)
    #ionisation source/sink
    quantities2d["ionisationloss"] = np.sum(rootgrp['eirene_mc_eael_she_bal'],axis=0)#13.6*1.602E-19*rootgrp["b2stel_sna_ion_bal"][0]
    #molecular eriene sink
    quantities2d["molecularloss"] = np.sum(rootgrp['eirene_mc_emel_she_bal'],axis=0)
    #recombination heat source/sink
    quantities2d["recombloss"] = np.sum(rootgrp['eirene_mc_epel_she_bal'],axis=0)#np.sum(rootgrp['eirene_mc_eiel_she_bal'],axis=0)
    # plt.plot(quantities2d["te"].flatten(),quantities2d["recombloss"].flatten(),linestyle="",marker = "o")
    # plt.show()
    #heat loss due to div of atom velocity
    quantities2d["Divuloss"] = np.array(rootgrp["b2sihs_divua_bal"])+np.array(rootgrp["b2sihs_divue_bal"])

    quantities2d["radialHeatFlux"] = (np.array(rootgrp['fhe_52'])+np.array(rootgrp['fhe_32'])+np.array(rootgrp['fhe_thermj']) + np.array(rootgrp['fhe_cond']) + np.array(rootgrp['fhe_dia']) + np.array(rootgrp['fhe_ecrb']) + 
    np.array(rootgrp['fhe_strange'])+np.array(rootgrp['fhe_pschused']))[1]#+np.array(rootgrp['fhi_32'])+np.array(rootgrp['fhi_52'])+np.array(rootgrp['fhi_cond'])+np.array(rootgrp['fhi_dia'])+np.array(rootgrp['fhi_ecrb'])+np.array(rootgrp['fhi_strange'])+
    #np.array(rootgrp['fhi_pschused'])+np.array(rootgrp['fhi_inert'])+np.array(rootgrp['fhi_vispar'])+np.array(rootgrp['fhi_anml'])+np.array(rootgrp['fhi_kevis']))[1]

    #charge exchange momentum loss
    # print(rootgrp)
    quantities2d["CXmomLoss"] = np.array(rootgrp['b2stcx_smq_bal'][0])

    return quantities2d

def unpackSOLPS(fileName,reverse,ring):
    rootgrp = Dataset(str(fileName), "r", format="NETCDF4")
    quantities2d = unpack2d(rootgrp)
    Xpoint = -1
    # plt.plot(rootgrp['bb'][2][ring],quantities2d["z"][ring])
    # plt.show()
    # plt.plot(quantities2d['r'][ring],quantities2d["z"][ring])
    # plt.show()
    SOLring1 = SOLring(Spar= np.cumsum(quantities2d["sdiff"][ring][::reverse][1:Xpoint]),
    Spol= np.cumsum(quantities2d["sdiffpol"][ring][::reverse][1:Xpoint]),
    Sdiff = quantities2d["sdiff"][ring][::reverse][1:Xpoint],
    ParArea = quantities2d["Area"][ring][::reverse][1:Xpoint],
    te = quantities2d["te"][ring][::reverse][1:Xpoint],
    ti = quantities2d["ti"][ring][::reverse][1:Xpoint],
    cond = -1*reverse*quantities2d["cond"][ring][::reverse][1:Xpoint],
    condI = -1*reverse*quantities2d["condI"][ring][::reverse][1:Xpoint],
    ne = quantities2d["ne"][ring][::reverse][1:Xpoint],
    ni = quantities2d["ni"][ring][::reverse][1:Xpoint],
    n0 = quantities2d["n0"][ring][::reverse][1:Xpoint],
    B= quantities2d["TotalField"][ring][::reverse][1:Xpoint],
    Bpol= quantities2d["Bpol"][ring][::reverse][1:Xpoint],
    qf = (quantities2d["imprad"]/quantities2d["V"])[ring][::reverse][1:Xpoint],
    n0rad = (quantities2d["neutralrad"]/quantities2d["V"])[ring][::reverse][1:Xpoint],
    ionisLoss = (quantities2d["ionisationloss"]/quantities2d["V"])[ring][::reverse][1:Xpoint],
    recombLoss = (quantities2d["recombloss"]/quantities2d["V"])[ring][::reverse][1:Xpoint],
    molecularloss = (quantities2d["molecularloss"]/quantities2d["V"])[ring][::reverse][1:Xpoint],
    divuloss = (quantities2d["Divuloss"]/quantities2d["V"])[ring][::reverse][1:Xpoint],
    conv = -1*reverse*quantities2d["conv"][ring][::reverse][1:Xpoint],
    radTrans=(quantities2d["radialHeatFlux"]/quantities2d["V"])[ring][::reverse][1:Xpoint],
    FlowVelocity = (quantities2d["vfluid"])[ring][::reverse][1:Xpoint],
    recombSource = (quantities2d["recombSource"])[ring][::reverse][1:Xpoint],
    ionisSource = (quantities2d["ionisSource"])[ring][::reverse][1:Xpoint],
    CXmomLoss = (quantities2d["CXmomLoss"])[ring][::reverse][1:Xpoint])
    return quantities2d,SOLring1
