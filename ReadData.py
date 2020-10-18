
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

# =====================================================

def get_weights_SPC(ps, trop, field):
    """
    from Stephen Po-Chedley
	copy from Mark's CMIP6_utils.py

    wts = get_weights(ps, trop, field)

    Function returns the layer weights (layer thickness) given the upper
    integration boundary (e.g., the tropopause, trop) and the lower bound
    (e.g., surface pressure, ps). Uses the input field to get pressure level
    information.
                    ps - response field [time, lat, lon]
                    trop - radiative kernel [time, lat, lon]
                    field - field that will be integrated [time, plev, lat, lon]

    """

    PAL = MV.zeros(ps.shape)
    PAT = MV.ones(ps.shape)
    trop_wts = MV.zeros(field.shape - MV.array([0,1,0,0]))
    atm_wts = MV.zeros(field.shape - MV.array([0,1,0,0]))
    plev = field.getLevel()[:]
    
    # make sure pressures are in Pa.
    if MV.maximum(plev)<=2000:
        plev=100*plev
    if MV.maximum(ps)<=2000:
        ps=100*ps
    if MV.maximum(trop)<=2000:
        trop=100*trop
    
    if plev[0] < plev[1]:
        raise ValueError('This script assumes that plev[0] > plev[1].')
    for i in range(len(plev)-1):
        # allocate slice
        sfield = MV.zeros(ps.shape)
        tfield = MV.zeros(ps.shape)
        # get first level above surface
        p = plev[i]
        pp = plev[i+1]
        ISURF = MV.greater(ps, p) # 1 where ps>p, 0 where ps<p
        II = MV.equal(PAL, 0)
        IA = MV.logical_and(ISURF, II)
        PAL = MV.where(IA, pp, PAL)
        # get first layer below tropopause
        ITROP = MV.less(trop, p)
        II = MV.greater_equal(trop, pp)
        IB = MV.logical_and(ITROP, II)
        PAT = MV.where(IB, p, PAT)
        # layers above tropopause or below surface (i.e., zero weight)
        above_trop = MV.logical_not(ITROP)
        below_surf = MV.logical_not(ISURF)
        IC = MV.logical_or(below_surf,above_trop)
        # layers adjacent to both tropopause and surface (i.e., layer is between tropopause and surface)
        ID = MV.logical_and(IA, IB)
        # layers not adjacent to tropopause or surface (i.e., full layers)
        IE = MV.logical_or(IA, IB)
        IE = MV.logical_and(MV.logical_not(IC), MV.logical_not(IE))
        # layers not adjacent to surface (i.e., full layers)
        IF = MV.logical_and(MV.logical_not(below_surf), MV.logical_not(IA))
        # TROPOSPHERIC VALUES
        sfield = MV.where(IA, ps-PAL, sfield)
        sfield = MV.where(IB, PAT-trop, sfield)
        sfield = MV.where(IC, 0., sfield)
        sfield = MV.where(ID, ps-trop, sfield)
        sfield = MV.where(IE, p - pp, sfield)
        # store field and weight by per 100 hPa (1/100 for Pa to hPa and 1/100 for per *100* hPa)
        trop_wts[:, i, :, :] = sfield / 10000.
        
        # ATMOSPHERIC VALUES
        tfield = MV.where(IA, ps-PAL, tfield)
        tfield = MV.where(below_surf, 0., tfield)
        tfield = MV.where(IF, p - pp, tfield)
        # store field and weight by per 100 hPa (1/100 for Pa to hPa and 1/100 for per *100* hPa)
        atm_wts[:, i, :, :] = tfield / 10000.

    return trop_wts,atm_wts

