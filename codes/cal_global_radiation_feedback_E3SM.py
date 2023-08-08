#!/usr/bin/env python
# coding: utf-8

# This script is used to calculate radiation response from AMIP and AMIP-p4K for E3SM.
# June 5, 2020 --- created by Yi Qin
# June 30, 2020 --- modified into a function to be called by main.py 

# Aug 24, 2021: 
# -------------------------------------------------------------------------------
# https://cdat.llnl.gov/documentation/utilities/utilities-1.html
# "The default action of the setTimeBoundsMonthly function assumes that the
# time point is at the beginning of the month. To compute bounds assuming that
# the time point at the end of the month."
# Line 170 -- the time point show at the end of the month. so I need to add '1' 
# there to keep the result correctly.
# -------------------------------------------------------------------------------



#IMPORT STUFF:
#=====================
import xarray as xr
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
import numpy.ma as ma
sys.path.append('../')
import cases_lookup as CL
import PlotDefinedFunction as PDF


########## MAIN SUBROUTINE STARTS HERE ....
def Global_RadFeedback(direc_data,case_stamp,yearS,yearE,fname1,fname2,outdir,exp1,exp2):

    if os.path.isfile(outdir+'global_mean_features_'+case_stamp+'.csv'):
        print('Global_RadFeedback is done.')
        return 

    yearS_4d = "{:04d}".format(yearS)
    yearE_4d = "{:04d}".format(yearE)

    nyears = yearE - yearS + 1
    
    #phases = ['CMIP5','CMIP6']
    phases = ['CMIP6']
    
    for phase in phases:
    #    print('----------Hi, here is processing ',phase,' Data----------------')
    
    
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
                        f1 = xr.open_dataset(direc_data+fname1+'/'+svar+'_'+exp1+'_'+yearS_4d+'01-'+yearE_4d+'12.nc')
                        pi_raw = f1[svar]
                        f1.close()
    
                        f2 = xr.open_dataset(direc_data+fname2+'/'+svar+'_'+exp2+'_'+yearS_4d+'01-'+yearE_4d+'12.nc')
                        ab_raw = f2[svar]
                        f2.close()
    
                        # reverse lev direction
                        if svar in var3d:
                            pi_raw = pi_raw[:,::-1,:,:]
                            ab_raw = ab_raw[:,::-1,:,:]

                        print(pi_raw.shape)
                        print(ab_raw.shape)
   
                        dic_all[svar+'_ano'] = xr.DataArray(ab_raw - pi_raw, coords=pi_raw.coords)
                        dic_all[svar+'_pi'] = pi_raw
                        dic_all[svar+'_ab'] = ab_raw

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
                
                # Define new time coordinate
                newtime = pd.date_range(start='1850-01-01', periods=dic_all['rsnt_pi'].shape[0], freq='MS')
        
                thesevars = []
                thesevars = np.append(thesevars, ['ts','SWCRE','LWCRE','netCRE','SWCLR','LWCLR','FTOA','FSNT','FLNT','FTOACLR','FSNTCLR','FLNTCLR'])  

                df_gavg = pd.DataFrame()

                dic_new = {}
                ### get monthly climatology
                for ivar,svar in enumerate(thesevars):
                    print('here we are processing ',svar)

                    dic_all[svar+'_ano'] = dic_all[svar+'_ano'].assign_coords({'time':("time",newtime)})
  
                    # get regression onto global mean surface air temperature anomaly
                    # 0. get annual mean of each variable
                    dic_new[svar+'_ano_ann'] = PDF.weighted_annual_mean(dic_all[svar+'_ano'].time,dic_all[svar+'_ano'])
                
                    # 1. get global mean surface air temperature series 
                    if svar=='ts':
                        dic_new[svar+'_ano_ann_gavg'] = PDF.area_averager(dic_new[svar+'_ano_ann'])

                    # get global and climatological mean - Note: average of annual mean and average of all monthly data result in slightly different results.
                    #dic_new[svar+'_ano_gaavg'] = PDF.area_averager((dic_new[svar+'_ano_ann'].mean(axis=0)))
                    dic_new[svar+'_ano_gaavg'] = PDF.area_averager((dic_all[svar+'_ano'].mean(axis=0)))

                    # 2. regression on ts anomaly loop over each variable
                    dic_new[svar+'_ano_slope'],dic_new[svar+'_ano_intercept'] = PDF.linearregression_nd(dic_new[svar+'_ano_ann'],x=np.reshape(dic_new['ts_ano_ann_gavg'].values, (dic_new[svar+'_ano_ann'].shape[0],1,1)))

                    # 3. get global mean of slope
                    dic_new[svar+'_ano_slope'] = xr.DataArray(dic_new[svar+'_ano_slope'], coords=dic_all[svar+'_ano'][0,:,:].coords)
                    dic_new[svar+'_ano_intercept'] = xr.DataArray(dic_new[svar+'_ano_intercept'], coords=dic_all[svar+'_ano'][0,:,:].coords)

                    dic_new[svar+'_ano_slope_gavg'] = PDF.area_averager(dic_new[svar+'_ano_slope'])

                    ### get global mean climatology
                    dic_new[svar+'_ano_gaavg_perK'] = dic_new[svar+'_ano_gaavg']/dic_new['ts_ano_gaavg']
                    

                    ## output to CSV file for further check
                    if 'coupled' in case_stamp:
                        ano_out = dic_new[svar+'_ano_gaavg']
                        ano_perK_out = dic_new[svar+'_ano_slope_gavg']
                    else:
                        ano_out = dic_new[svar+'_ano_gaavg']
                        ano_perK_out = dic_new[svar+'_ano_gaavg_perK']

                    tmp = pd.DataFrame([[svar,str(np.round(ano_out.values,3)),str(np.round(ano_perK_out.values,3))]],columns=['varname','anomaly','anomaly_perK'])

                    print(tmp)
                    df_gavg = pd.concat([df_gavg,tmp])
                    print(df_gavg)
                    
                print(df_gavg)
    
                df_gavg.to_csv(outdir+'global_mean_features_'+case_stamp+'.csv')


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if __name__ == "__main__":

    #direc_data = '/compyfs/qiny108/diag_feedback_E3SM_postdata/'
    direc_data = '/compyfs/qiny108/colla/diag_feedback_E3SM_postdata/'

    case_stamp = 'v2test'
    yearS = 2
    yearE = 3
    fname1,_,_ = CL.get_lutable(case_stamp,'amip')
    fname2,_,_ = CL.get_lutable(case_stamp,'amip4K')
    outdir = './'
    figdir = './'
    exp1 = 'FC5'
    exp2 = 'FC5_4K'
    
    Global_RadFeedback(direc_data,case_stamp,yearS,yearE,fname1,fname2,outdir,exp1,exp2)


