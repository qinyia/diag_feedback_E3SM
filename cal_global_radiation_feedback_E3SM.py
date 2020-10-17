#!/usr/bin/env python
# coding: utf-8

# This script is used to calculate radiation response from AMIP and AMIP-p4K for E3SM.
# June 5, 2020 --- created by Yi Qin
# June 30, 2020 --- modified into a function to be called by main.py 

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
import cartopy.crs as ccrs
import cartopy
import matplotlib.pyplot as plt
import os
import pandas as pd
import cdtime
import time
import genutil
import numpy.ma as ma
from genutil import statistics

########## MAIN SUBROUTINE STARTS HERE ....
def Global_RadFeedback(direc_kernel,direc_data,case_stamp,yearS,yearE,fname1,fname2,outdir):

    yearS_4d = "{:04d}".format(yearS)
    yearE_4d = "{:04d}".format(yearE)

    nyears = yearE - yearS + 1
    
    #phases = ['CMIP5','CMIP6']
    phases = ['CMIP6']
    
    for phase in phases:
    #    print('----------Hi, here is processing ',phase,' Data----------------')
    
    #    direc_kernel='/work/qin4/script/cloud-radiative-kernels/'
    #    direc_data = '/work/qin4/Data/E3SM_data/'
    
        All_project_cntl = ['CMIP','CMIP']
        All_exp_cntl = ['piControl','amip']
        All_project_new = ['CMIP','CFMIP']
        if phase == 'CMIP5':
            All_exp_new = ['abrupt4xCO2','amip4K']
        else:
            All_exp_new = ['abrupt-4xCO2','amip-p4K']
    
        if phase == 'CMIP5':
            All_nyears=[150,27]
        else:
    #        All_nyears=[150,36]
            All_nyears=[150,nyears]
    
        All_startyear = [1,1]
        
        comp='atmos'
        freq='mon'
        var = ["rlut","rsnt","rlutcs","rsntcs","ts"]
        var3d = []
    
        used_models = ['E3SM-1-0']
        used_ripfs = ['r1i1p1f1']
    
        # --------------------------- March 3: modified --- based on needed variables, identify models which provide these: End --------
        for icase in range(1,2):
        
            project_cntl = All_project_cntl[icase]
            exp_cntl = All_exp_cntl[icase]
            project_new = All_project_new[icase]
            exp_new = All_exp_new[icase]
            nyears = All_nyears[icase]
            startyear = All_startyear[icase]
    
            for imod in range(len(used_models)): 
       
                # read 2D variables

                dic_all = {}
                for svar in var:
        
                    if svar in var:
                        print(svar)
                        f1 = cdms.open(direc_data+fname1+'/'+svar+'_FC5_'+yearS_4d+'01-'+yearE_4d+'12.nc')
                        pi_raw = f1(svar)
                        f1.close()
    
                        f2 = cdms.open(direc_data+fname2+'/'+svar+'_FC5_4K_'+yearS_4d+'01-'+yearE_4d+'12.nc')
                        ab_raw = f2(svar)
                        f2.close()
    
                        # reverse lev direction
                        if svar in var3d:
                            pi_raw = pi_raw[:,::-1,:,:]
                            ab_raw = ab_raw[:,::-1,:,:]

                        print(pi_raw.shape)
                        print(ab_raw.shape)
   
                        dic_all[svar+'_ano'] = ab_raw - pi_raw
                        dic_all[svar+'_pi'] = pi_raw
                        dic_all[svar+'_ab'] = ab_raw

                        dic_all[svar+'_ano'].setAxisList(pi_raw.getAxisList())
                        dic_all[svar+'_pi'].setAxisList(pi_raw.getAxisList())
                        dic_all[svar+'_ab'].setAxisList(pi_raw.getAxisList())

                        del(pi_raw, ab_raw)
                    else:
                        print('we dont find this variable:',svar,' in your file. Please check it!')
                
                ### get SWCF, LWCF and their anomalies
                dic_all['SWCRE_ano'] = dic_all['rsnt_ano'] - dic_all['rsntcs_ano']
                dic_all['LWCRE_ano'] = dic_all['rlutcs_ano'] - dic_all['rlut_ano']
                dic_all['netCRE_ano'] = dic_all['SWCRE_ano'] + dic_all['LWCRE_ano']
                dic_all['SWCLR_ano'] = dic_all['rsntcs_ano']
                dic_all['LWCLR_ano'] = -1.*dic_all['rlutcs_ano']
                dic_all['FTOA_ano'] = dic_all['rsnt_ano']- dic_all['rlut_ano']
                dic_all['FSNT_ano'] = dic_all['rsnt_ano']
                dic_all['FLNT_ano'] = -1.*dic_all['rlut_ano']
                dic_all['FTOACLR_ano'] = dic_all['rsntcs_ano'] - dic_all['rlutcs_ano']
                dic_all['FSNTCLR_ano'] = dic_all['rsntcs_ano']
                dic_all['FLNTCLR_ano'] = -1.*dic_all['rlutcs_ano']
            
                dic_all['SWCRE_pi'] = dic_all['rsnt_pi'] - dic_all['rsntcs_pi']
                dic_all['LWCRE_pi'] = dic_all['rlutcs_pi'] - dic_all['rlut_pi']
                dic_all['netCRE_pi'] = dic_all['SWCRE_pi'] + dic_all['LWCRE_pi']
                dic_all['SWCLR_pi'] = dic_all['rsntcs_pi']
                dic_all['LWCLR_pi'] = -1.*dic_all['rlutcs_pi']
                dic_all['FTOA_pi'] = dic_all['rsnt_pi'] - dic_all['rlut_pi']
                dic_all['FSNT_pi'] = dic_all['rsnt_pi']
                dic_all['FLNT_pi'] = -1.*dic_all['rlut_pi']
                dic_all['FTOACLR_pi'] = dic_all['rsntcs_pi'] - dic_all['rlutcs_pi']
                dic_all['FSNTCLR_pi'] = dic_all['rsntcs_pi']
                dic_all['FLNTCLR_pi'] = -1.*dic_all['rlutcs_pi']
                
                dic_all['SWCRE_ab'] = dic_all['rsnt_ab'] - dic_all['rsntcs_ab']
                dic_all['LWCRE_ab'] = dic_all['rlutcs_ab'] - dic_all['rlut_ab']
                dic_all['netCRE_ab'] = dic_all['SWCRE_ab'] + dic_all['LWCRE_ab']
                dic_all['SWCLR_ab'] = dic_all['rsntcs_ab']
                dic_all['LWCLR_ab'] = -1.*dic_all['rlutcs_ab']
                dic_all['FTOA_ab'] = dic_all['rsnt_ab'] - dic_all['rlut_ab']
                dic_all['FSNT_ab'] = dic_all['rsnt_ab']
                dic_all['FLNT_ab'] = -1.*dic_all['rlut_ab']
                dic_all['FTOACLR_ab'] = dic_all['rsntcs_ab'] - dic_all['rlutcs_ab']
                dic_all['FSNTCLR_ab'] = dic_all['rsntcs_ab']
                dic_all['FLNTCLR_ab'] = -1.*dic_all['rlutcs_ab']
                
                AXL = dic_all['rlutcs_pi'].getAxisList()
        
                thesevars = []
                thesevars = np.append(thesevars, ['ts','SWCRE','LWCRE','netCRE','SWCLR','LWCLR','FTOA','FSNT','FLNT','FTOACLR','FSNTCLR','FLNTCLR'])  
                df_gavg = pd.DataFrame()

                dic_new = {}
                ### get monthly climatology
                for ivar,svar in enumerate(thesevars):
                    print('here we are processing ',svar)
  
                    # Mar 11: add setting time bounds to avoid the 'Nontype error from time.py' while get annual mean using cdutil.YEAR
                    cdutil.setTimeBoundsMonthly(dic_all[svar+'_ano'])
    
                    # get regression onto global mean surface air temperature anomaly
                    # 0. get annual mean of each variable
                    dic_new[svar+'_ano_annual'] = cdutil.YEAR(dic_all[svar+'_ano'])
                
                    # 1. get global mean surface air temperature series 
                    if svar=='ts':
                        dic_new[svar+'_ano_annual_gavg'] = cdutil.averager(dic_new[svar+'_ano_annual'],axis='xy',weights='weighted')

                    # 2. regression on ts anomaly loop over each variable
                    dic_new[svar+'_ano_slope'],dic_new[svar+'_ano_intercept'] = statistics.linearregression(dic_new[svar+'_ano_annual'],x=dic_new['ts_ano_annual_gavg'])
                    print(dic_new['ts_ano_slope'].shape, dic_new['ts_ano_intercept'].shape)

                    # 3. get global mean of slope
                    dic_new[svar+'_ano_slope_gavg'] = cdutil.averager(dic_new[svar+'_ano_slope'],axis='xy',weights='weighted')

                    # 4. print and test 
                    if ivar == len(thesevars)-1:
                        print(dic_new['ts_ano_slope_gavg'],dic_new['SWCRE_ano_slope_gavg'], dic_new['LWCRE_ano_slope_gavg'])
                    
                    ### get climatology 
                    dic_new[svar+'_ano_aavg'] = MV.average(dic_all[svar+'_ano'], axis=0)
                    dic_new[svar+'_pi_aavg'] = MV.average(dic_all[svar+'_pi'], axis=0)
                    dic_new[svar+'_ab_aavg'] = MV.average(dic_all[svar+'_ab'], axis=0)

                    ### get global mean climatology
                    dic_new[svar+'_ano_gaavg'] = cdutil.averager(dic_new[svar+'_ano_aavg'], axis='xy',weights='weighted')
                    dic_new[svar+'_pi_gaavg'] = cdutil.averager(dic_new[svar+'_pi_aavg'], axis='xy',weights='weighted')
                    dic_new[svar+'_ab_gaavg'] = cdutil.averager(dic_new[svar+'_ab_aavg'], axis='xy',weights='weighted')
                    
                    ### get annual mean anomaly per K
                    dic_new[svar+'_ano_aavg_perK'] = dic_new[svar+'_ano_aavg']/dic_new['ts_ano_gaavg']
                    
                    ### print them
                    ano_gaavg = dic_new[svar+'_ano_gaavg']
                    pi_gaavg = dic_new[svar+'_pi_gaavg']
                    ab_gaavg = dic_new[svar+'_ab_gaavg']
                    
                    print(svar+'_ano_gaavg is: '+str(np.round(ano_gaavg,3)))
                    print(svar+'_pi_gaavg is: '+str(np.round(pi_gaavg,3)))
                    print(svar+'_ab_gaavg is: '+str(np.round(ab_gaavg,3)))
                    
                    ## output to CSV file for further check
                    tmp = pd.DataFrame([[svar,str(np.round(ano_gaavg,3)),str(np.round(ano_gaavg/dic_new['ts_ano_gaavg'],3))]],columns=['varname','anomaly','anomaly_perK'])
                    print(tmp)
                    df_gavg = pd.concat([df_gavg,tmp])
                    print(df_gavg)
                    
                print(df_gavg)
    
                df_gavg.to_csv(outdir+'global_mean_features_'+exp_cntl+'_'+used_models[imod]+'_'+case_stamp+'.csv')


