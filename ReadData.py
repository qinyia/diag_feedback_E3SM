
import cdms2 as cdms
import cdutil
import MV2 as MV
import numpy as np
import pylab as pl
import matplotlib as mpl
import sys
import matplotlib.pyplot as plt
import os
import pandas as pd
import cdtime
import time
import ReadData as RD
from genutil import statistics
import genutil
from scipy import stats
import numpy.ma as ma
import scipy as sp
import scipy.interpolate

# Definition of thermodynamic constant:
Cpd = 1005.7  # Heat capacity at constant pressure for dry air [J kg^-1 K^-1]
Cpv = 1870.0  # Heat capacity at constant pressure of water vapor [J kg^-1 K^-1]
Cl = 4190.0   # Heat capacity of liquid water [J kg^-1 K^-1]
Rv = 461.5    # Gas constant of water vapor [J kg^-1 K^-1]
Rd = 287.04   # Gas constant of dry air [J kg^-1 K^-1]
Lv0 = 2.501e6 # Latent heat of vaporization at 0 deg C [J kg^-1]
g = 9.80616   # Accelleration of gravity [m s^-2]
rho_l = 1000.0 # Density of water [kg m^-3]
p_0 = 100000.  # Reference pressure [Pa]
epsilon = Rd/Rv

### March 12: get pressure thickness
def get_dp(psurf,ptop,level):

    # ---- we need the level variable is from bottom to top like: (1000, 925, 850, ...)
    nlev = len(level)

    # define new array to save the 'fake' interface level
    p1 = np.zeros((nlev+1))
    p1[0] = psurf
    p1[1:] = level

    p2 = np.zeros((nlev+1))
    p2[:-1] = level
    p2[-1] = ptop

    # get the interface level here
    ilev1 = (p1+p2)/2.

    # the final pressure thickness: 
    dp = ilev1[:-1] - ilev1[1:]

    # quick check the impact of dp
#    dp = [37.5, 37.5, 75, 112.5, 125, 100, 100, 100, 75, 50, 50, 50, 40, 25, 20, 15, 10, 15,10]
#    print(dp)

    return dp

#   Saturation vapor pressure over water (Pa) (Emanuel) from temperature T (K)
def e_star(T):
#    xx = 100.*np.exp(53.67957 - 6743.769/T - 4.8451*np.log(T))
    print(type(T))
    esw = 6.1078*np.exp(17.2693882*(T-273.16)/(T-35.86)) # unit: hPa
    esi = 6.1078*np.exp(21.8745584*(T-276.16)/(T-7.66)) # unit: hPa
    es = np.where(T<273.16,esi,esw)
    es = MV.masked_where(T.mask == True,es*100.)  # convert to Pa
    print(genutil.minmax(es))
    return es

#   saturation mixing ratio rs (Kg/Kg) from temperature t (K) and pressure P (Pa)
def r_star(p, T):
    rr = epsilon*e_star(T)/(p - e_star(T))
    hus = rr/(1.+rr)
    return hus

def e_star_GG(T):
   # use the Goff-Gratch equation: unit of T is K, unit of e_star is hPa. 
   # to liquid (-49.9 ~ 100 degC)
   T00 = 273.16
   t1 = 10.79574*(1.-T00/T)
   t2 = 5.02800*MV.log10(T/T00)
   t3 = 1.50475*1e-4*(1.-10**(-8.2969*(T/T00-1)))
   t4 = 0.42873*1e-3*(10**(4.76955*(1.-T/T00))-1.)
   t5 = 0.78614
   esw = 10**(t1-t2+t3+t4+t5)

   # to ice (-100 ~ 0 degC)
   t1 = -9.09685*(T00/T-1.)
   t2 = -3.56654*MV.log10(T00/T)
   t3 = 0.87682*(1-T/T00)
   t4 = 0.78614
   esi = 10**(t1+t2+t3+t4)

   es = np.where(T<T00-49.9,esi,esw)
   es = MV.masked_where(T.mask == True, es*100.) # convert to Pa
   print(genutil.minmax(es))
   return es

def r_star_GG(p, T):
    rr = epsilon*e_star_GG(T)/(p - e_star_GG(T))
    hus = rr/(1.+rr)
    return hus

# April 18 
#   potential temperature [K] from pressure P [Pa] and temperature T [K]
def theta(p, T): return T*(p_0/p)**(Rd/Cpd)

#   temperature [K] from pressure P [Pa] and potential temperature theta [K]
def theta_to_T(theta, p): return theta*(p/p_0)**(0.2854)


def get_tropopause_pressure(Ta0):

    # April 27: steal from Mark's codes
    # Perform tropopause pressure calculation
    # generated via: f2py -c -m tropo tropo.f90
    # tropo.f90 came from http://www.inscc.utah.edu/~reichler/research/projects/TROPO/code.txt

    import tropo

    print("Get_tropopause_pressure: Ta0's shape is", Ta0.shape)

    PP = Ta0.getLevel()
    if PP.units == 'Pa':
        PP=Ta0.getLevel()[:]/100. #kern_plev/100.
    else:
        PP=Ta0.getLevel()[:]
    print('PP is',PP)

    kern_lons = Ta0.getLongitude()[:]
    plimu=45000
    pliml=7500
    plimlex=7500
    dofill=0
    # fortran wants this in lon,lat,ctp but lets swap time with lon
    tropp=MV.zeros((Ta0[:,0,:].shape))
    tropp=MV.masked_where(tropp==0,tropp)
    for LO in range(len(kern_lons)): 
        temp=Ta0[:,:,:,LO](order=[0,2,1])
        tp,tperr=tropo.tropo(temp,PP,plimu,pliml,plimlex,dofill)

        tropp[:,:,LO]=tp/100. # convert to hPa
#        print('tropp minmax is',genutil.minmax(tropp[:,:,LO]))

    tropp=MV.masked_where(tropp>plimu/100.,tropp)
    tropp=MV.masked_where(tropp<0.,tropp)

    return tropp


