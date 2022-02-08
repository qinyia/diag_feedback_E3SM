#****************************************************************
#
#    Filename: cal_EIS.py
#
#    Author: Yi Qin - qin4@llnl.gov
#    Description: calculate EIS
#    Input: post-processed output from raw model data
#    Output: global_EIS_xxxx.nc
#    Create: 2022-02-01
#    Last Modified: 2022-02-01 12:57:44
#****************************************************************

#IMPORT STUFF:
#=====================
import cdms2 as cdms
import cdutil
import MV2 as MV
import numpy as np
import pylab as pl
import matplotlib as mpl
import sys

## qinyi 
import matplotlib.pyplot as plt
import os
import pandas as pd
import cdtime
import time
import ReadData as RD
import genutil
from genutil import statistics
import numpy.ma as ma
from global_land_mask import globe 
 
###########################################################################
# HELPFUL FUNCTIONS FOLLOW
###########################################################################
def calc_EIS(Ts, SLP, T700, z700):
    '''
    Input:
           Ts: surface temperature (K)
           SLP: sea level pressure (Pa)
           T700: temperature at 700 hPa (K)
           z700: height at 700 hPa (m)
    Output: 
           EIS: estimated inversion strengh (K)
           LTS: lower troposheric strength (K)

    # Notes:
    #  - surface relative humidity is assumed to be constant
    #  - potential temperature is calcualted with respect to the 1000 hPa
    #    reference level
    #  - traditionally, surface air temperature is used in place of surface
    #    temperature.
    #  - physical constants are from Curry and Webster, Appendix B
    #
    # Physical Constants:
    #  - Gas constant for water vapor = 461 J/K/kg
    #  - Latent heat of vaporization at 273.15 K = 2.5*10^6 J/kg
    #  - Gas constant for dry air = 287 J/K/kg
    #  - Specific heat at constant pressure = 1004 J/K/kg
    #
    #References:
    #
    #Curry, J. A., & Webster, P. J. (1998). Thermodynamics of Atmospheres and 
    #Oceans. Elsevier.
    #
    #Georgakakos, K. P., & Bras, R. L. (1984). A hydrologically useful station
    #precipitation model: 1. Formulation. Water Resources Research, 20(11), 
    #1585-1596, https://doi.org/10.1029/WR020i011p01585
    #
    #Wood, R., & Bretherton, C. S. (2006). On the relationship between 
    #stratiform low cloud cover and lower-tropospheric stability. Journal of 
    #Climate, 19(24), 6425-6432, https://doi.org/10.1175/JCLI3988.1

    '''

    #calculate the LCL height
    es=6.11*np.exp(((2.5*10**6)/461)*((1/273.15)-(1/Ts))) #Clausius-Clapeyron Equation (Equation 4.23 from Curry and Webster)
    e=es*80/100 #assume RH is 80%, as in WB06
    Td=((1/273.15)-(461/(2.5*10**6))*np.log(e/6.11))**(-1) #form of Clausius-Clapeyron Equation
    p_LCL=SLP*(((Ts-Td)/223.15)+1)**(-3.5) #Equation (12) from Georgakakos and Bras (1984)
    T_LCL=Ts*(((Ts-Td)/223.15)+1)**(-1) #Equation (13) from Georgakakos and Bras (1984)
    LCL=((287*(Ts+T_LCL)/2)/9.8)*np.log(SLP/p_LCL) #Hypsometric equation

    #calculate LTS
    theta700=T700*(1000/700)**(287/1004)
    theta_slp=Ts*(1000./(SLP*0.01))**(287/1004)
    LTS=theta700-theta_slp

    #calculate the moist adaiabtic lapse rate at 850 hPa
    T_bar=(Ts+T700)/2 #approximation of the temperature at 850 hPa, following WB06
    es_bar=6.11*np.exp(((2.5*10**6)/461)*((1/273.15)-(1/T_bar))) #Clausius-Clapeyron Equation (Equation 4.23 from Curry and Webster)
    qs_bar=0.622*es_bar/(850-es_bar) #Equation 4.37 of Curry and Webster
    gamma_m=(9.8/1004)*(1-((1+(((2.5*10**6)*qs_bar)/(287*T_bar)))/(1+((((2.5*10**6)**2)*qs_bar)/(1004*461*(T_bar**2)))))) #Equation (5) of WB06
    
    EIS=LTS-gamma_m*(z700-LCL); #Equation (4) of WB06

    return EIS, LTS

#################################################################################
def cal_EIS(direc_data,case_stamp,yearS,yearE,fname1,fname2,outdir,figdir,exp1,exp2):

    outfile = outdir+'global_EIS_'+case_stamp+'.nc'

    if os.path.isfile(outfile):
        print('cal_EIS', case_stamp, 'output is ready. Please continue. ')
    else:
        yearS_4d = "{:04d}".format(yearS)
        yearE_4d = "{:04d}".format(yearE)
        nyears = yearE - yearS + 1
    
        direc_data1 = direc_data+'/'+fname1+'/'
        direc_data2 = direc_data+'/'+fname2+'/'
    
        used_models = 'E3SM-1-0'
        
        yrS=yearS
        yrE=yearE
        monS=1
        monE=12
        
        yrS_4d='{:04d}'.format(yrS)
        yrE_4d='{:04d}'.format(yrE)
        monS_2d='{:02d}'.format(monS)
        monE_2d='{:02d}'.format(monE)
    
        lats = np.arange(-90,92.5,2.5)
        lons = np.arange(1.25,360,2.5)
        grid = cdms.createGenericGrid(lats,lons)

        var3d = ['Z3','T','OMEGA']
        var2d = ['psl','ts','tas']

        var = var2d + var3d 

        # =============================================
        # read 2D variables
        # =============================================
        dic_all = {}
        for svar in var:
            if svar in var:
                print(svar)
                f1 = cdms.open(direc_data+fname1+'/'+svar+'_'+exp1+'_'+yearS_4d+'01-'+yearE_4d+'12.nc')
                pi_raw = f1(svar)
                f1.close()

                f2 = cdms.open(direc_data+fname2+'/'+svar+'_'+exp2+'_'+yearS_4d+'01-'+yearE_4d+'12.nc')
                ab_raw = f2(svar)
                f2.close()

                # reverse lev direction
                if svar in var3d:
                    # by default, the level for raw model ouput is in unit of hPa. Be cautious.
                    if svar in ['Z3','T']:
                        spec_lev = 700
                    elif svar in ['OMEGA']:
                        spec_lev = 500

                    pi_raw = pi_raw.pressureRegrid(cdms.createAxis([spec_lev]))[:,0,:,:]
                    ab_raw = ab_raw.pressureRegrid(cdms.createAxis([spec_lev]))[:,0,:,:]

                #----------------------------------------------------------
                # regrid data                 
                #----------------------------------------------------------
                pi_raw_grd = pi_raw.regrid(grid,regridTool='esmf',regridMethod='linear')
                ab_raw_grd = ab_raw.regrid(grid,regridTool='esmf',regridMethod='linear')

                print('pi_raw_grd.shape = ',pi_raw_grd.shape,genutil.minmax(pi_raw_grd))
                print('ab_raw_grd.shape = ',ab_raw_grd.shape,genutil.minmax(ab_raw_grd))

                dic_all[svar+'_ano'] = ab_raw_grd - pi_raw_grd
                dic_all[svar+'_pi'] = pi_raw_grd
                dic_all[svar+'_ab'] = ab_raw_grd

                dic_all[svar+'_ano'].setAxisList(pi_raw_grd.getAxisList())
                dic_all[svar+'_pi'].setAxisList(pi_raw_grd.getAxisList())
                dic_all[svar+'_ab'].setAxisList(pi_raw_grd.getAxisList())

                del(pi_raw, ab_raw, pi_raw_grd, ab_raw_grd)
            else:
                print('we dont find this variable:',svar,' in your file. Please check it!')
 
        #----------------------------------------------------------
        # calculate global mean surface air temperature anomaly  
        #----------------------------------------------------------
        anomtas = cdutil.ANNUALCYCLE.climatology(dic_all['tas_ano']) #(12, 90, 144)
        avgdtas = cdutil.averager(MV.average(anomtas,axis=0), axis='xy', weights='weighted') # (scalar)

        print('avgdtas = ',avgdtas)

        #----------------------------------------------------------
        # calculate EIS
        #----------------------------------------------------------
        EIS_pi,LTS_pi = calc_EIS(dic_all['ts_pi'], dic_all['psl_pi'], dic_all['T_pi'], dic_all['Z3_pi'])
        EIS_ab,LTS_ab = calc_EIS(dic_all['ts_ab'], dic_all['psl_ab'], dic_all['T_ab'], dic_all['Z3_ab'])

        dic_all['EIS_pi'] = EIS_pi
        dic_all['EIS_ab'] = EIS_ab
        dic_all['EIS_ano'] = EIS_ab - EIS_pi

        AXL = dic_all['ts_pi'].getAxisList()
        dic_all['EIS_ano'].setAxisList(AXL)
        dic_all['EIS_pi'].setAxisList(AXL)
        dic_all['EIS_ab'].setAxisList(AXL)

        #----------------------------------------------------------
        # get time-series of annual mean for regression if exp is coupled
        #----------------------------------------------------------
        dic_all2 = {}
        for svar in dic_all.keys():
            if '_ano' in svar:
                cdutil.setTimeBoundsMonthly(dic_all[svar])
                dic_all2[svar+'_ann'] = cdutil.YEAR(dic_all[svar])

        tas_ano_ann_gm = cdutil.averager(dic_all2['tas_ano_ann'], axis='xy', weights='weighted')

        #----------------------------------------------------------
        # calculate climatological control state 
        #----------------------------------------------------------
        dic_out = {}
        for svar in dic_all.keys():
            if 'EIS' not in svar and 'OMEGA' not in svar: # only output EIS
                continue

            if '_pi' in svar:
                dic_out[svar+'_clim'] = MV.average(dic_all[svar],axis=0)
            if '_ab' in svar:
                dic_out[svar+'_clim'] = MV.average(dic_all[svar],axis=0)
            elif '_ano' in svar:
                if 'coupled' not in case_stamp:
                    dic_out[svar+'_clim'] = MV.average(dic_all[svar],axis=0)/avgdtas
                else:
                    print('doing regression...')
                    dic_out[svar+'_clim'],intercept = statistics.linearregression(dic_all2[svar+'_ann'], x=tas_ano_ann_gm)

        print('dic_out.keys() = ',dic_out.keys())

        #----------------------------------------------------------
        # save data into file     
        #----------------------------------------------------------
        value = 0
        cdms.setNetcdfShuffleFlag(value) ## where value is either 0 or 1
        cdms.setNetcdfDeflateFlag(value) ## where value is either 0 or 1
        cdms.setNetcdfDeflateLevelFlag(value) ## where value is a integer between 0 and 9 included

        fout = cdms.open(outfile,'w')

        for svar in dic_out.keys():
            print('svar = ', svar)
            tmp = dic_out[svar] 
            fout.write(tmp, id = svar)
            fout.comment = '_pi_clim is the climatological control state, _ab_clim is the climatological warming state and xxx_ano_clim is the anomaly normalized by global mean tas anomaly'

        fout.close()
        
