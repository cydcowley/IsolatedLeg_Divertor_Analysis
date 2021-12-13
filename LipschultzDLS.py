import scipy as sp
from scipy.optimize import fsolve
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy import interpolate
import winsound


def ChInt(s, B_field, sx, L,sh = 0):
    B_field = interpolate.interp1d(s,B_field,kind='cubic',fill_value="extrapolate")
    answer = quad(integrand,sh,L,args = (sx, L, B_field), epsabs = 0.0000000000000000001)[0]
    answer = (answer)**(-2/7)
    answer = answer*B_field(sh)/B_field(sx)
    return answer

def ChIntThermalForce(s, B_field, sx, L,sh = 0):
    B_field = interpolate.interp1d(s,B_field,kind='cubic')
    answer = quad(integrand,sh,L,args = (sx, L, B_field), epsabs = 0.0000000000000000001)[0]
    answer = ((answer)**(-2/7))**(1.71)
    answer = answer*B_field(sh)/B_field(sx)
    return answer

def Tu(s, B_field, sx, L,Beta = 1, sh = 0):
    B_field = interpolate.interp1d(s,B_field,kind='cubic')
    intTop = quad(integrand,sh,L,args = (sx, L, B_field), epsabs = 0.0000000000000000001)[0]
    answer = (intTop)**(2/7)
    return answer

def Leff(s,B_field,sx,L,sh=0):
    B_field = interpolate.interp1d(s,np.add(np.zeros(len(B_field)),1),kind='cubic')
    answer = quad(integrand,sh,L,args = (sx, L, B_field), epsabs = 0.0000000000000000001)[0]
    return answer   

def averageB(s, B_field, sx, L, sh = 0):
    B_field = interpolate.interp1d(s,B_field,kind='cubic',fill_value="extrapolate")
    answer = quad(integrand,sh,L,args = (sx, L, B_field), epsabs = 0.0000000000000000001)[0]
    answer = answer/(L-sh)
    return answer

    return answer

def secondTerm(z, B_field, zx, L,Beta = 1, zh = 0):
    gradB = interpolate.interp1d(z,np.gradient(B_field)/np.gradient(z),kind='cubic')
    B_field = interpolate.interp1d(z,B_field,kind='cubic')
    intTop = quad(integrand,zh,L,args = (zx, L, B_field), epsabs = 0.0000000000000000001)[0]
    answer = (2/7)*L*(B_field(zh)**2)*qi(zh,zx,L)/intTop
    return L*(gradB(zh))/(B_field(zh)),-1*answer

def stability(z, B_field, zx, L,Beta = 1, zh = 0):
    B_field = interpolate.interp1d(z,B_field)
    zdh = zh +0.001
    int0 = quad(integrand,zh,L,args = (zx, L, B_field), epsabs = 0.0000000000000000001)[0]
    int1 = quad(integrand,zdh,L,args = (zx, L, B_field), epsabs = 0.0000000000000000001)[0]
    answer0 = (int0**(2/7))/B_field(zh)
    answer1 = (int1**(2/7))/B_field(zdh)
    return (answer1-answer0)

def qi(z,zx,L):
    if z<zx:
        return -(L-zx)
    elif z>=zx:
        return -(L-z)

def qf(z,zx,L,B_field):
    qf = quad(integrand,z,L,args = (zx, L, B_field), epsabs = 0.0000000000000000001)[0]
    qf = qf**(2/7)/B_field(z)
    return -qf


def plotqfqi(z, B_field, zx, L,Beta = 1, zh = 0):
    # B_field = np.multiply(B_field,0)
    # B_field = np.add(B_field,1)
    B_field = interpolate.interp1d(z,B_field)
    qI = []
    qF = []

    
    testqfs = []
    A = np.linspace(0,50,200)
    for a in A:
        testqfs.append(qf(zh,zx,L,B_field)*a)
    f = interpolate.interp1d(testqfs,A)
    Amplitude = f(qi(zh,zx,L))
    for i in range(len(z)):
        qI.append(qi(z[i],zx,L))
        qF.append(Amplitude*qf(z[i],zx,L,B_field))
    plt.plot(z/L,qI,label="qi")
    plt.plot(z/L,qF, label = "qf")
    plt.xlabel("zh/L")
    plt.ylabel("Power Flux")
    title = "Zx/L = " + str(zx/L) + ", Bx/Bt = " + str(B_field(zx)/B_field(zh))+ "zh = " + str(zh)
    plt.legend()
    plt.title(title)
    plt.show()



def integrand(s,sx, L, B_field):
    if s<sx:
        return (B_field(s)/B_field(sx))
    else:
        return  (L-s)*(B_field(s)/B_field(sx))/(L-sx)
    


