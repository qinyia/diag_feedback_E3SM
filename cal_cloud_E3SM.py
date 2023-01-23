#****************************************************************
#
#    Filename: cal_cloud_E3SM.py
#
#    Author: Yi Qin - qin4@llnl.gov
#    Description: calculate climatology of control simulation; calculate response per global-mean surface temperature change
#    Input: post-processed output from raw model data
#    Output: global_cloud_xxxx.nc
#    Create: 2021-07-18
#    Last Modified: 2021-07-26 18:08:44
#                   Jul 26, 2021: add low cloud fraction, liquid water path and ice water path
#                   Feb 01, 2022: add output of '_ab'
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
from joblib import Parallel,delayed
import multiprocessing
 
###########################################################################
# HELPFUL FUNCTIONS FOLLOW
###########################################################################
def cal_cloud(direc_data,case_stamp,yearS,yearE,fname1,fname2,outdir,figdir,exp1,exp2):

    outfile = outdir+'global_cloud_'+case_stamp+'.nc'

    if os.path.isfile(outfile):
        print('cal_cloud', case_stamp, 'output is ready. Please continue. ')
        return 

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

    var3d = ['CLDLIQ', 'CLDICE', 'CLOUD']
    var2d = ['tas','TGCLDIWP','TGCLDLWP','CLDLOW','CLDMED','CLDHGH','CLDTOT','INCLDLWP','INCLDIWP']

    var = var2d + var3d 
    # =============================================
    # read 2D variables
    # =============================================
    dic_all = {}
    for svar in var:
        if svar == 'INCLDLWP':
            print(svar)
            svar_a = 'TGCLDLWP'
            f1 = cdms.open(direc_data+fname1+'/'+svar_a+'_'+exp1+'_'+yearS_4d+'01-'+yearE_4d+'12.nc')
            pi_raw_a = f1(svar_a)
            f1.close()
            f2 = cdms.open(direc_data+fname2+'/'+svar_a+'_'+exp2+'_'+yearS_4d+'01-'+yearE_4d+'12.nc')
            ab_raw_a = f2(svar_a)
            f2.close()
            svar_b = 'CLDTOT'
            f1 = cdms.open(direc_data+fname1+'/'+svar_b+'_'+exp1+'_'+yearS_4d+'01-'+yearE_4d+'12.nc')
            pi_raw_b = f1(svar_b)
            f1.close()
            f2 = cdms.open(direc_data+fname2+'/'+svar_b+'_'+exp2+'_'+yearS_4d+'01-'+yearE_4d+'12.nc')
            ab_raw_b = f2(svar_b)
            f2.close()

            pi_raw = pi_raw_a/pi_raw_b
            ab_raw = ab_raw_a/ab_raw_b
        elif svar == 'INCLDIWP':
            print(svar)
            svar_a = 'TGCLDIWP'
            f1 = cdms.open(direc_data+fname1+'/'+svar_a+'_'+exp1+'_'+yearS_4d+'01-'+yearE_4d+'12.nc')
            pi_raw_a = f1(svar_a)
            f1.close()
            f2 = cdms.open(direc_data+fname2+'/'+svar_a+'_'+exp2+'_'+yearS_4d+'01-'+yearE_4d+'12.nc')
            ab_raw_a = f2(svar_a)
            f2.close()
            svar_b = 'CLDTOT'
            f1 = cdms.open(direc_data+fname1+'/'+svar_b+'_'+exp1+'_'+yearS_4d+'01-'+yearE_4d+'12.nc')
            pi_raw_b = f1(svar_b)
            f1.close()
            f2 = cdms.open(direc_data+fname2+'/'+svar_b+'_'+exp2+'_'+yearS_4d+'01-'+yearE_4d+'12.nc')
            ab_raw_b = f2(svar_b)
            f2.close()

            pi_raw = pi_raw_a/pi_raw_b
            ab_raw = ab_raw_a/ab_raw_b

        else:
            print(svar)
            f1 = cdms.open(direc_data+fname1+'/'+svar+'_'+exp1+'_'+yearS_4d+'01-'+yearE_4d+'12.nc')
            pi_raw = f1(svar)
            f1.close()

            f2 = cdms.open(direc_data+fname2+'/'+svar+'_'+exp2+'_'+yearS_4d+'01-'+yearE_4d+'12.nc')
            ab_raw = f2(svar)
            f2.close()

            # reverse lev direction
            if svar in var3d:
                pi_raw = pi_raw[:,::-1,:,:]
                ab_raw = ab_raw[:,::-1,:,:]

        #----------------------------------------------------------
        # regrid data                 
        #----------------------------------------------------------
        pi_raw_grd = pi_raw.regrid(grid,regridTool='esmf',regridMethod='linear')
        ab_raw_grd = ab_raw.regrid(grid,regridTool='esmf',regridMethod='linear')

        print('pi_raw_grd.shape = ',pi_raw_grd.shape)
        print('ab_raw_grd.shape = ',ab_raw_grd.shape)

        dic_all[svar+'_ano'] = ab_raw_grd - pi_raw_grd
        dic_all[svar+'_pi'] = pi_raw_grd
        dic_all[svar+'_ab'] = ab_raw_grd

        dic_all[svar+'_ano'].setAxisList(pi_raw_grd.getAxisList())
        dic_all[svar+'_pi'].setAxisList(pi_raw_grd.getAxisList())
        dic_all[svar+'_ab'].setAxisList(pi_raw_grd.getAxisList())

        del(pi_raw, ab_raw, pi_raw_grd, ab_raw_grd)

    #----------------------------------------------------------
    # calculate global mean surface air temperature anomaly  
    #----------------------------------------------------------
    anomtas = cdutil.ANNUALCYCLE.climatology(dic_all['tas_ano']) #(12, 90, 144)
    avgdtas = cdutil.averager(MV.average(anomtas,axis=0), axis='xy', weights='weighted') # (scalar)

    print('avgdtas = ',avgdtas)

    # get time-series of annual mean
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
        
