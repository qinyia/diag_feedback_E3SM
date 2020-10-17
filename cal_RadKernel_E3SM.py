
## this script is used to calculate radiative kernel following Soden et al. and Shell et al.
## some trival but important things need to be taken care:
## 1. albedo should be in percent (%) unit;
## 2. pay attention to the masking, which should be level < p_tropopause;

## created: March 12, 2020 by Yi Qin
## modified on June 30, 2020 --- change into a function used by main.py

import cdms2 as cdms
import cdutil
import MV2 as MV
import numpy as np
import pylab as pl
import matplotlib as mpl
mpl.use('Agg')
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
from scipy.interpolate import interp1d
import ReadData as RD


cdms.axis.latitude_aliases.append("Y")
cdms.axis.longitude_aliases.append("X")

########## MAIN SUBROUTINE STARTS HERE ....

def get_fdbk(DATA,AXL2d,ts_ano_grd_avg,ts_ano_grd_ann,name,exp_cntl,dic_gfdbk,dic_feedback):
    if exp_cntl == 'amip':
        cdutil.setTimeBoundsMonthly(DATA)
        DATA_ann = cdutil.YEAR(DATA)
        newdata1 = MV.average(DATA_ann,axis=0)/MV.average(ts_ano_grd_avg)
        print(newdata1.shape)
        newdata1.setAxisList(AXL2d)
        dic_gfdbk[name+'_gfdbk'] = newdata1
        newdata2 = cdutil.averager(newdata1,axis='xy',weights='weighted')
        dic_feedback[name+'_feedback'] = newdata2
    elif exp_cntl == 'piControl':
        cdutil.setTimeBoundsMonthly(DATA)
        DATA_ann = cdutil.YEAR(DATA)
        slope,intercept = genutil.statistics.linearregression(DATA_ann,x = ts_ano_grd_ann)
        dic_gfdbk[name+'_gfdbk'] = slope
        newdata = cdutil.averager(slope,axis='xy',weights='weighted')
        dic_feedback[name+'_feedback'] = newdata
    return dic_gfdbk,dic_feedback
 

def RadKernel(kernel_dir,direc_data,case_stamp,yearS,yearE,fname1,fname2,outdir,figdir):

    yearS_4d = "{:04d}".format(yearS)
    yearE_4d = "{:04d}".format(yearE)
    nyears = yearE - yearS + 1
 
    get_tropo_Soden_method = True
    get_tropo_Mark_method = False
    
    phases = ['CMIP6']

    # figure font
    fh = 17
    
    for phase in phases:
        print('----------Hi, here is processing ',phase,' Data----------------')
    
        ### Jan 23: merge amip and cmip experiments in the same script. here.
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
            All_nyears=[150,nyears]
        
        comp='atmos'
        freq='mon'
    #    var2d = ["ts","rlut","rsut","rlutcs","rsutcs","rsdt","rsus","rsds","rsuscs","rsdscs"]
        var2d = ['ts','rsnt','rsntcs','rlut','rlutcs','rsdscs','rsuscs','rsds','rsus']
        var3d = ['ta','hus']
        var = var2d + var3d
    
        usedin = 'RadKernel'
    
        used_models = ['E3SM-1-0']
    
        #for icase in range(len(All_exp_cntl)):
    #    for icase in range(0,1): # piControl
        for icase in range(1,2): # amip
        
            project_cntl = All_project_cntl[icase]
            exp_cntl = All_exp_cntl[icase]
            project_new = All_project_new[icase]
            exp_new = All_exp_new[icase]
        
            for imod in range(len(used_models)): 
       
                dic_all = {}

                # read 2D variables
                for svar in var:
       
                    if svar in var:
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
                dic_all['SWCRE_ano']  = dic_all['rsnt_ano'] - dic_all['rsntcs_ano']
                dic_all['LWCRE_ano']  = dic_all['rlutcs_ano'] - dic_all['rlut_ano']
                dic_all['netCRE_ano'] = dic_all['SWCRE_ano'] + dic_all['LWCRE_ano']
           
                dic_all['SWCRE_pi'] = dic_all['rsnt_pi'] - dic_all['rsntcs_pi']
                dic_all['LWCRE_pi'] = dic_all['rlutcs_pi'] - dic_all['rlut_pi']
                dic_all['netCRE_pi'] = dic_all['SWCRE_pi'] + dic_all['LWCRE_pi']
               
                dic_all['SWCRE_ab'] = dic_all['rsnt_ab'] - dic_all['rsntcs_ab']
                dic_all['LWCRE_ab'] = dic_all['rlutcs_ab'] - dic_all['rlut_ab']
                dic_all['netCRE_ab'] = dic_all['SWCRE_ab'] + dic_all['LWCRE_ab']
    
                ## get albedo
                dic_all['alb_pi'] = dic_all['rsus_pi']/dic_all['rsds_pi'] * 100.
                dic_all['alb_ab'] = dic_all['rsus_ab']/dic_all['rsds_ab'] *100.
                dic_all['alb_pi'] = MV.masked_outside(dic_all['alb_pi'], 0.0, 100.)
                dic_all['alb_ab'] = MV.masked_outside(dic_all['alb_ab'], 0.0, 100.)
                dic_all['alb_ano'] = dic_all['alb_ab'] - dic_all['alb_pi'] 
                dic_all['alb_pi'].setAxisList(dic_all['rsus_pi'].getAxisList())
                dic_all['alb_ab'].setAxisList(dic_all['rsus_pi'].getAxisList())
                dic_all['alb_ano'].setAxisList(dic_all['rsus_pi'].getAxisList())
                
                dic_all['alb_clr_pi'] = dic_all['rsuscs_pi']/dic_all['rsdscs_pi'] * 100.
                dic_all['alb_clr_ab'] = dic_all['rsuscs_ab']/dic_all['rsdscs_ab'] *100.
                dic_all['alb_clr_pi'] = MV.masked_outside(dic_all['alb_clr_pi'], 0.0, 100.)
                dic_all['alb_clr_ab'] = MV.masked_outside(dic_all['alb_clr_ab'], 0.0, 100.)
                dic_all['alb_clr_ano'] = dic_all['alb_clr_ab'] - dic_all['alb_clr_pi']
                dic_all['alb_clr_pi'].setAxisList(dic_all['rsus_pi'].getAxisList())
                dic_all['alb_clr_ab'].setAxisList(dic_all['rsus_pi'].getAxisList())
                dic_all['alb_clr_ano'].setAxisList(dic_all['rsus_pi'].getAxisList())
     
    ##================================================================================###
                ###-------------------------__Read kernel data _-----------------------###
                
                ### attention: the kernel data's latitude is from north to south [90....-90]
                ### we need to reverse it.

                f1 = cdms.open(kernel_dir+"RRTMG_ts_toa_cld_highR.nc")
                dic_all['ts_kernel_cld'] = f1('lwkernel')[:,::-1,:]
                print(dic_all['ts_kernel_cld'].shape)
                f1.close()
                f1 = cdms.open(kernel_dir+"RRTMG_ts_toa_clr_highR.nc")
                dic_all['ts_kernel_clr'] = f1('lwkernel')[:,::-1,:]
                f1.close()
    
                f1 = cdms.open(kernel_dir+"RRTMG_t_toa_cld_highR.nc")
                dic_all['t_kernel_cld'] = f1('lwkernel')[:,:,::-1,:]
                f1.close()
                f1 = cdms.open(kernel_dir+"RRTMG_t_toa_clr_highR.nc")
                dic_all['t_kernel_clr'] = f1('lwkernel')[:,:,::-1,:]
                f1.close()
    
                f1 = cdms.open(kernel_dir+"RRTMG_wv_sw_toa_cld_highR.nc")
                dic_all['wv_sw_kernel_cld'] = f1('swkernel')[:,:,::-1,:]
                f1.close()
                f1 = cdms.open(kernel_dir+"RRTMG_wv_sw_toa_clr_highR.nc")
                dic_all['wv_sw_kernel_clr'] = f1('swkernel')[:,:,::-1,:]
                f1.close()
    
                f1 = cdms.open(kernel_dir+"RRTMG_wv_lw_toa_cld_highR.nc")
                dic_all['wv_lw_kernel_cld'] = f1('lwkernel')[:,:,::-1,:]
                print(dic_all['wv_lw_kernel_cld'].shape)
                f1.close()
                f1 = cdms.open(kernel_dir+"RRTMG_wv_lw_toa_clr_highR.nc")
                dic_all['wv_lw_kernel_clr'] = f1('lwkernel')[:,:,::-1,:]
                f1.close()
    
                f1 = cdms.open(kernel_dir+"RRTMG_alb_toa_cld_highR.nc")
                dic_all['alb_kernel_cld'] = f1('swkernel')[:,::-1,:]
                f1.close()
                f1 = cdms.open(kernel_dir+"RRTMG_alb_toa_clr_highR.nc")
                dic_all['alb_kernel_clr'] = f1('swkernel')[:,::-1,:]
                f1.close()
    
                if len(dic_all['t_kernel_cld'].getLevel()[:]) > len(dic_all['ta_ano'].getLevel()[:]):
                    HAX = dic_all['t_kernel_cld'].getAxisList()
                    stdlevs=dic_all['ta_ano'].getLevel()[:]
                    if dic_all['ta_ano'].getLevel().units == 'Pa':
                        stdlevs=stdlevs/100.
                    kern_plev = cdms.createAxis(stdlevs)
                    kern_plev.designateLevel()
                    kern_plev.units = "hPa"
                    kern_plev.id = "plev"
                    print(stdlevs)
                    ## vertically regrid kernel pressure to model pressure, because model pressure grid is less than kernel one
                    ## model: 19 levels; kernel: 24 levels
                    these_vars = ['t_kernel_cld','t_kernel_clr','wv_sw_kernel_cld','wv_sw_kernel_clr','wv_lw_kernel_cld','wv_lw_kernel_clr']
                else:
                    HAX = dic_all['ta_ano'].getAxisList()
                    stdlevs=dic_all['t_kernel_cld'].getLevel()[:]
                    if dic_all['t_kernel_cld'].getLevel().units == 'Pa':
                        stdlevs = stdlevs/100.
                    kern_plev = cdms.createAxis(stdlevs)
                    kern_plev.designateLevel()
                    kern_plev.units = "hPa"
                    kern_plev.id = "plev"
                    print(stdlevs)
                    ## vertically regrid model pressure to kernel pressure
                    these_vars = ['ta_ano','hus_ano','ta_pi','ta_ab','hus_pi','hus_ab']
    
                for ivar in these_vars:
                    print(ivar)
                    DATA = dic_all[ivar]
                    rawlevs = DATA.getLevel()[:]
                    if DATA.getLevel().units == 'Pa':
                        rawlevs = rawlevs/100.
                    #print(rawlevs)
                    #print(stdlevs)
                    # "steal" from Mark's script: load_TOA_kernels.py -- March 12
                    # I use pressureRegrid function before, but it always informs the missing value after regrid.
                    # just use Mark's method here, and later check what's wrong with the pressureRegrid. 
                    DATA2=MV.where(DATA.mask,np.nan,DATA)
                    f = interp1d(rawlevs,DATA2,axis=1,fill_value="extrapolate")
                    newdata = f(stdlevs.data)
                    newdata=MV.masked_where(np.isnan(newdata),newdata)
                    newdata.setAxis(0,HAX[0])
                    newdata.setAxis(1,kern_plev)
                    newdata.setAxis(2,HAX[2])
                    newdata.setAxis(3,HAX[3])
                    dic_all[ivar] = newdata
    
                # --------------------------------------------------------------------
                # expand kernel data into several years * 12 month
                # --------------------------------------------------------------------
                dic_mon = {}
                stdtime=dic_all['ta_ano'].getTime()[:]
                kern_time = cdms.createAxis(stdtime)
                kern_time.designateTime()
                kern_time.units = dic_all['ta_ano'].getTime().units
                kern_time.id = "time"
    
                these_vars = ['ts_kernel_cld','ts_kernel_clr','alb_kernel_cld','alb_kernel_clr','t_kernel_cld','t_kernel_clr','wv_sw_kernel_cld','wv_sw_kernel_clr','wv_lw_kernel_cld','wv_lw_kernel_clr']
    
                HAX = dic_all['t_kernel_cld'].getAxisList()
    
                for ivar in these_vars:
                    DATA = dic_all[ivar]
                    if len(DATA.shape) == 3: # (time,lat,lon)
                        newdata2 = np.tile(DATA,(nyears,1,1))
                        newdata2.setAxis(0,kern_time)
                        newdata2.setAxis(1,HAX[2])
                        newdata2.setAxis(2,HAX[3])
                    elif len(DATA.shape) == 4: # (time,lev,lat,lon)
                        newdata2 = np.tile(DATA,(nyears,1,1,1))
                        newdata2.setAxis(0,kern_time)
                        newdata2.setAxis(1,HAX[1])
                        newdata2.setAxis(2,HAX[2])
                        newdata2.setAxis(3,HAX[3])
                    dic_mon[ivar+'_monthly'] = newdata2

                # --------------------------------------------------------------------
                # regrid the input data to kernel grid
                # --------------------------------------------------------------------
                dic_grd = {}
                if get_tropo_Soden_method:
                    these_vars = ['ts_ano','alb_ano','alb_clr_ano','ta_ano','hus_ano','ta_pi','hus_pi','hus_ab','rsntcs_ano','rlutcs_ano','SWCRE_ano','LWCRE_ano','netCRE_ano','rlut_ano','rsnt_ano']
                elif get_tropo_Mark_method:
                    these_vars = ['ts_ano','alb_ano','alb_clr_ano','ta_pi','ta_ab','hus_pi','hus_ab','rsntcs_ano','rlutcs_ano','SWCRE_ano','LWCRE_ano','netCRE_ano','rlut_ano','rsnt_ano']
                
                kernel_grid = dic_all['ts_kernel_cld'].getGrid() # get new grid info
                for ivar in these_vars:
                    DATA = dic_all[ivar]
                    newdata = DATA.regrid(kernel_grid,regridTool='esmf',regridMethod='linear')
                    print(newdata.shape) #(12*nyears, 64, 128)
                    dic_grd[ivar+'_grd'] = newdata
                    print(ivar,newdata.shape)

                del(dic_all)
    
    
                #################################################################################
                #### Temperature feedback calculation
                #################################################################################
                dic_kern = {}
    
                lat = dic_mon['t_kernel_cld_monthly'].getLatitude()[:]
                lon = dic_mon['t_kernel_cld_monthly'].getLongitude()[:]
                time = dic_mon['t_kernel_cld_monthly'].getTime()[:]
                lev = dic_mon['t_kernel_cld_monthly'].getLevel()[:]
                nlat = len(lat)
                nlon = len(lon)
                nlev = len(lev)
                ntime = len(time)
    
                AXL4d = dic_mon['t_kernel_cld_monthly'].getAxisList()
                AXL3d = dic_mon['ts_kernel_cld_monthly'].getAxisList()
                AXL2d = dic_mon['ts_kernel_cld_monthly'][0,:,:].getAxisList()
               
                # multiply monthly mean TS change by the TS kernels (function of lat, lon, month) (units W/m2)
    
                dic_grd['ts_ano_grd'][dic_grd['ts_ano_grd'].mask] = 0
                dic_kern['dLW_ts'] = dic_mon['ts_kernel_cld_monthly'] * dic_grd['ts_ano_grd']
                dic_kern['dLW_ts'].setAxisList(AXL3d)
    
                # for clear-sky
                dic_kern['dLW_ts_clr'] = dic_mon['ts_kernel_clr_monthly'] * dic_grd['ts_ano_grd']
                dic_kern['dLW_ts_clr'].setAxisList(AXL3d)
    
                # =========================================================================================================#
                # get tropopause height 
                # =========================================================================================================#
    
                # April 27, 2020: update the calculation of tropopause pressure by using the time-varying one. 
                # refer to Mark et al. (2020)
    
                if get_tropo_Soden_method:
                    # crude tropopause estimate: 100 hPa in the tropics, lowering with 
                    # cosine to 300 hPa at the poles. !!!!!!!!!--------- need to be updated into varying tropopause calculation.
                    pi = 3.1415926
                    x = np.cos(lat*pi/100.)
                    p_tropopause_zonalmean = 300-200*x
                    dic_kern['p_tropopause'] = cdms.asVariable(np.transpose(np.tile(p_tropopause_zonalmean,(ntime,nlev,nlon,1)),(0,1,3,2)))
                    dic_kern['p_tropopause'].setAxisList(AXL4d)
                    dic_kern['lev_4d'] = np.transpose(np.tile(lev,(ntime,nlat,nlon,1)),(0,3,1,2))
                    print(dic_kern['lev_4d'].shape)
    
                elif get_tropo_Mark_method:
                    p_tropopause_3d = RD.get_tropopause_pressure(ta_pi_grd)
                    # expand to 4D variable including pressure level
                    dic_kern['p_tropopause'] = np.transpose(np.tile(p_tropopause_3d,(nlev,1,1,1)),(1,0,2,3))
                    print(dic_kern['p_tropopause'].shape)
    
                    dic_kern['lev_4d'] = np.transpose(np.tile(lev,(ntime,nlat,nlon,1)),(0,3,1,2))
                    print(dic_kern['lev_4d'].shape)
    
                    dic_grd['ta_pi_grd'] = MV.masked_where(dic_kern['lev_4d']<=dic_kern['p_tropopause'],dic_grd['ta_pi_grd'])
                    dic_grd['hus_pi_grd'] = MV.masked_where(dic_kern['lev_4d']<=dic_kern['p_tropopause'],dic_grd['hus_pi_grd'])
    
                    p_tropopause_3d = RD.get_tropopause_pressure(dic_grd['ta_ab_grd'])
                    # expand to 4D variable including pressure level
                    dic_kern['p_tropopause'] = np.transpose(np.tile(p_tropopause_3d,(nlev,1,1,1)),(1,0,2,3))
                    print(dic_kern['p_tropopause'].shape)
    
                    dic_grd['ta_ab_grd'] = MV.masked_where(dic_kern['lev_4d']<=dic_kern['p_tropopause'],dic_grd['ta_ab_grd'])
                    dic_grd['hus_ab_grd'] = MV.masked_where(dic_kern['lev_4d']<=dic_kern['p_tropopause'],dic_grd['hus_ab_grd'])
    
                    # get ta and hus anomaly again
                    dic_grd['ta_ano_grd'] = dic_grd['ta_ab_grd'] - dic_grd['ta_pi_grd']
                    dic_grd['hus_ano_grd'] = dic_grd['hus_ab_grd'] - dic_grd['hus_pi_grd']
                    dic_grd['ta_ano_grd'].setAxisList(AXL4d)
                    dic_grd['hus_ano_grd'].setAxisList(AXL4d)
    
                # set the temperature change to zero in the stratosphere (mask out stratosphere)
                if get_tropo_Soden_method:
                    dic_grd['ta_ano_grd'] = MV.masked_where(dic_kern['lev_4d']<=dic_kern['p_tropopause'],dic_grd['ta_ano_grd']) 
                    dic_mon['t_kernel_cld_monthly'] = MV.masked_where(dic_kern['lev_4d']<=dic_kern['p_tropopause'],dic_mon['t_kernel_cld_monthly']) 
                    dic_mon['t_kernel_clr_monthly'] = MV.masked_where(dic_kern['lev_4d']<=dic_kern['p_tropopause'],dic_mon['t_kernel_clr_monthly']) 
    
                # get air temperature feedback
                dic_grd['ta_ano_grd'][dic_grd['ta_ano_grd'].mask]=0
                dic_kern['dLW_ta'] = dic_mon['t_kernel_cld_monthly'] * dic_grd['ta_ano_grd']
                dic_kern['dLW_ta'].setAxisList(AXL4d)
    
                # for clear-sky
                dic_kern['dLW_ta_clr'] = dic_mon['t_kernel_clr_monthly'] * dic_grd['ta_ano_grd']
                dic_kern['dLW_ta_clr'].setAxisList(AXL4d)
    
           
                # get pressure thickness
                psurf = 1000.
                ptop = min(lev)
                dp = RD.get_dp(psurf, ptop, lev)
                dic_kern['dp_4d'] = np.transpose(np.tile(np.tile(dp,[1,1,1,1]),[ntime,nlat,nlon,1]),(0,3,1,2))
    
                # vertical integral (sum) -- first get weights: dp
                dic_kern['dLW_ta_psum'] = MV.sum(dic_kern['dLW_ta']*dic_kern['dp_4d']/100.,axis=1) 
                print(dic_kern['dLW_ta_psum'].shape)
                # for clear-sky
                dic_kern['dLW_ta_clr_psum'] = MV.sum(dic_kern['dLW_ta_clr']*dic_kern['dp_4d']/100.,axis=1) 
    
                # get annual global-mean surface temperature anomaly time series
                ts_ano_grd_avg = cdutil.averager(dic_grd['ts_ano_grd'],axis='xy',weights='weighted')
                cdutil.setTimeBoundsMonthly(ts_ano_grd_avg)
                ts_ano_grd_ann = cdutil.YEAR(ts_ano_grd_avg)

                # --------------------- kernel feedback ---------------------------------------------
                dic_gfdbk = {}
                dic_feedback = {}

                names = ['dLW_ts','dLW_ta_psum','dLW_ts_clr','dLW_ta_clr_psum']
                for name in names:
                    DATA = dic_kern[name]
                    result = get_fdbk(DATA,AXL2d,ts_ano_grd_avg,ts_ano_grd_ann,name,exp_cntl,dic_gfdbk,dic_feedback)
   
                dic_kern['dLW_t_clr_psum'] = dic_kern['dLW_ta_clr_psum'] + dic_kern['dLW_ts_clr']
                dic_kern['dLW_t_psum'] = dic_kern['dLW_ta_psum'] + dic_kern['dLW_ts']
                dic_kern['dLW_t_clr_psum'].setAxisList(AXL3d)
                dic_kern['dLW_t_psum'].setAxisList(AXL3d)
    
                dic_gfdbk['dLW_t_gfdbk'] = dic_gfdbk['dLW_ts_gfdbk'] + dic_gfdbk['dLW_ta_psum_gfdbk']
                dic_gfdbk['dLW_t_clr_gfdbk'] = dic_gfdbk['dLW_ts_clr_gfdbk'] + dic_gfdbk['dLW_ta_clr_psum_gfdbk']
                dic_gfdbk['dLW_t_gfdbk'].setAxisList(AXL2d)
                dic_gfdbk['dLW_t_clr_gfdbk'].setAxisList(AXL2d)
    
                dic_gfdbk['t_gfdbk'] = dic_gfdbk['dLW_t_gfdbk']
                dic_gfdbk['t_clr_gfdbk'] = dic_gfdbk['dLW_t_clr_gfdbk']
    
                dic_feedback['t_feedback'] = dic_feedback['dLW_ts_feedback'] + dic_feedback['dLW_ta_psum_feedback']
                dic_feedback['t_clr_feedback'] = dic_feedback['dLW_ts_clr_feedback'] + dic_feedback['dLW_ta_clr_psum_feedback']
    
                print("Temperature feedback: ",dic_feedback['t_feedback'], "W/m2/K")
                print("clr Temperature feedback: ",dic_feedback['t_clr_feedback'], "W/m2/K")
    
                #################################################################################
                #### Planck feedback calculation
                #################################################################################
    
                # extend to 4-d
                dic_grd['ts_ano_grd_4d'] = cdms.asVariable(np.transpose(np.tile(dic_grd['ts_ano_grd'],(nlev,1,1,1)),(1,0,2,3)))
                print(dic_grd['ts_ano_grd_4d'].shape)
    
                # mask stratosphere
                if get_tropo_Soden_method:
                    dic_grd['ts_ano_grd_4d'] = MV.masked_where(dic_kern['lev_4d'] <=dic_kern['p_tropopause'],dic_grd['ts_ano_grd_4d'])
                dic_grd['ts_ano_grd_4d']   = MV.masked_where(dic_grd['ta_ano_grd'].mask==True,dic_grd['ts_ano_grd_4d'])
            
                dic_grd['ts_ano_grd_4d'][dic_grd['ts_ano_grd_4d'].mask]=0
                dic_kern['dLW_planck'] = dic_mon['t_kernel_cld_monthly'] * dic_grd['ts_ano_grd_4d']
                dic_kern['dLW_planck'].setAxisList(AXL4d)
    
                dic_kern['dLW_planck_clr'] = dic_mon['t_kernel_clr_monthly'] * dic_grd['ts_ano_grd_4d']
                dic_kern['dLW_planck_clr'].setAxisList(AXL4d)
    
                dic_kern['dLW_planck_psum'] = MV.sum(dic_kern['dLW_planck']*dic_kern['dp_4d']/100.,axis=1) + dic_kern['dLW_ts']
                dic_kern['dLW_planck_clr_psum'] = MV.sum(dic_kern['dLW_planck_clr']*dic_kern['dp_4d']/100.,axis=1) + dic_kern['dLW_ts_clr']
                dic_kern['dLW_planck_psum'].setAxisList(AXL3d)
                dic_kern['dLW_planck_clr_psum'].setAxisList(AXL3d)
    
                names = ['dLW_planck_psum','dLW_planck_clr_psum']
                for name in names:
                    DATA = dic_kern[name]
                    result = get_fdbk(DATA,AXL2d,ts_ano_grd_avg,ts_ano_grd_ann,name,exp_cntl,dic_gfdbk,dic_feedback)
 
                dic_gfdbk['planck_gfdbk'] = dic_gfdbk['dLW_planck_psum_gfdbk']
                dic_gfdbk['planck_clr_gfdbk'] = dic_gfdbk['dLW_planck_clr_psum_gfdbk']
                
                dic_feedback['planck_feedback'] = dic_feedback['dLW_planck_psum_feedback']
                dic_feedback['planck_clr_feedback'] = dic_feedback['dLW_planck_clr_psum_feedback']
    
                print('Planck feedback: ',dic_feedback['planck_feedback'],'W/m2/K')
                print('clr Planck feedback: ',dic_feedback['planck_clr_feedback'],'W/m2/K')
    
                #################################################################################
                #### Lapse rate feedback calculation
                #################################################################################
    
                # difference between ta and ts
                dic_grd['dt_ano'] = dic_grd['ta_ano_grd'] - dic_grd['ts_ano_grd_4d']
                dic_grd['dt_ano'].setAxisList(AXL4d)
    
                # mask stratosphere
                if get_tropo_Soden_method:
                    dic_grd['dt_ano'] = MV.masked_where(dic_kern['lev_4d'] <=dic_kern['p_tropopause'],dic_grd['dt_ano'])
            
                dic_grd['dt_ano'][dic_grd['dt_ano'].mask]=0
                dic_kern['dLW_lapserate'] = dic_mon['t_kernel_cld_monthly'] * dic_grd['dt_ano']
                dic_kern['dLW_lapserate_clr'] = dic_mon['t_kernel_clr_monthly'] * dic_grd['dt_ano']
                dic_kern['dLW_lapserate'].setAxisList(AXL4d)
                dic_kern['dLW_lapserate_clr'].setAxisList(AXL4d)
    
                dic_kern['dLW_lapserate_psum']     = MV.sum(dic_kern['dLW_lapserate']*dic_kern['dp_4d']/100.,axis=1)
                dic_kern['dLW_lapserate_clr_psum'] = MV.sum(dic_kern['dLW_lapserate_clr']*dic_kern['dp_4d']/100.,axis=1)
                dic_kern['dLW_lapserate_psum'].setAxisList(AXL3d)
                dic_kern['dLW_lapserate_clr_psum'].setAxisList(AXL3d)
    
                names = ['dLW_lapserate_psum','dLW_lapserate_clr_psum']
                for name in names:
                    DATA = dic_kern[name]
                    result = get_fdbk(DATA,AXL2d,ts_ano_grd_avg,ts_ano_grd_ann,name,exp_cntl,dic_gfdbk,dic_feedback)
 
                dic_gfdbk['lapserate_gfdbk']     = dic_gfdbk['dLW_lapserate_psum_gfdbk']
                dic_gfdbk['lapserate_clr_gfdbk'] = dic_gfdbk['dLW_lapserate_clr_psum_gfdbk']
    
                dic_feedback['lapserate_feedback']     = dic_feedback['dLW_lapserate_psum_feedback']
                dic_feedback['lapserate_clr_feedback'] = dic_feedback['dLW_lapserate_clr_psum_feedback']
                
                print('Lapse rate feedback: ',dic_feedback['lapserate_feedback'],'W/m2/K')
                print('clr Lapse rate feedback: ',dic_feedback['lapserate_clr_feedback'],'W/m2/K')
    
                print('------Summary of temperature feedback -------------------------')
                print('Lapse rate feedback: ',dic_feedback['lapserate_feedback'],'W/m2/K')
                print('Planck feedback: ',dic_feedback['planck_feedback'],'W/m2/K')
                print('Lapse rate + Planck feedback: ',dic_feedback['lapserate_feedback']+dic_feedback['planck_feedback'],'W/m2/K')
                print('Temperature feedback: ',dic_feedback['t_feedback'],'W/m2/K')
                print('clr Lapse rate feedback: ',dic_feedback['lapserate_clr_feedback'],'W/m2/K')
                print('clr Planck feedback: ',dic_feedback['planck_clr_feedback'],'W/m2/K')
                print('clr Lapse rate + clr Planck feedback: ',dic_feedback['lapserate_clr_feedback']+dic_feedback['planck_clr_feedback'],'W/m2/K')
                print('clr Temperature feedback: ',dic_feedback['t_clr_feedback'],'W/m2/K')
    
                #################################################################################
                #### Albedo feedback calculation
                #################################################################################
    
                dic_grd['alb_ano_grd'][dic_grd['alb_ano_grd'].mask]=0
                dic_kern['dSW_alb'] = dic_grd['alb_ano_grd'] * dic_mon['alb_kernel_cld_monthly']
                dic_kern['dSW_alb'].setAxisList(AXL3d)
    
                dic_grd['alb_clr_ano_grd'][dic_grd['alb_clr_ano_grd'].mask]=0
                dic_kern['dSW_alb_clr'] = dic_grd['alb_clr_ano_grd'] * dic_mon['alb_kernel_clr_monthly']
                dic_kern['dSW_alb_clr'].setAxisList(AXL3d)
    
                names = ['dSW_alb','dSW_alb_clr']
                for name in names:
                    DATA = dic_kern[name]
                    result = get_fdbk(DATA,AXL2d,ts_ano_grd_avg,ts_ano_grd_ann,name,exp_cntl,dic_gfdbk,dic_feedback)
   
                dic_gfdbk['alb_gfdbk']     = dic_gfdbk['dSW_alb_gfdbk']
                dic_gfdbk['alb_clr_gfdbk'] = dic_gfdbk['dSW_alb_clr_gfdbk']
    
                dic_feedback['alb_feedback'] = dic_feedback['dSW_alb_feedback']
                dic_feedback['alb_clr_feedback'] = dic_feedback['dSW_alb_clr_feedback']
    
                print("Surface albedo feedback: ",dic_feedback['alb_feedback'], "W/m2/K")
                print("clr Surface albedo feedback: ",dic_feedback['alb_clr_feedback'], "W/m2/K")
    
                #################################################################################
                #### Water vapor feedback calculation
                #################################################################################
    
                ## calculate the change in mixing ratio for 1K warming at constant RH
                # 'logq2-logq1' -- get logq by direct method: dlogq = log(q2) - log(q1)
                # 'angie_github' -- this method produces the water vapor feedback so large. It cannot pass the clear-sky linear test
    
                dlogq_method = 'logq2-logq1' # 'angie_github'
    
                qs0 = RD.r_star_GG(dic_kern['lev_4d']*100.,dic_grd['ta_pi_grd'])
                rh0 = dic_grd['hus_pi_grd']/qs0
    
                ta1k = dic_grd['ta_pi_grd']+1.0
                qs1 = RD.r_star_GG(dic_kern['lev_4d']*100.,ta1k)
                dic_kern['q1k'] = rh0*qs1
                del(qs0,rh0,ta1k,qs1)
    
                ## prevent negative values of moisture and its change
                q1k = MV.where(dic_kern['q1k'] < dic_grd['hus_pi_grd'], dic_grd['hus_pi_grd'], dic_kern['q1k'])
    
                if dlogq_method == 'logq2-logq1':
                    dic_kern['dlogq1k'] = np.log(dic_kern['q1k']) - np.log(dic_grd['hus_pi_grd'])
                elif dlogq_method == 'angie_github':
                    dic_kern['dq1k'] = dic_kern['q1k'] - dic_grd['hus_pi_grd']
                    dic_kern['dlogq1k'] = dic_kern['dq1k']/dic_grd['hus_pi_grd']

                dic_kern['dlogq1k'].setAxisList(AXL4d)
    
                if dlogq_method == 'logq2-logq1':
                    dic_kern['dlogq'] = np.log(dic_grd['hus_ab_grd']) - np.log(dic_grd['hus_pi_grd'])
                elif dlogq_method == 'angie_github':
                    dic_kern['dlogq'] = dic_grd['hus_ano_grd']/dic_grd['hus_pi_grd']
    
                dic_kern['dlogq2'] = dic_kern['dlogq']/dic_kern['dlogq1k']
    
                ## mask the statrosphere
                if get_tropo_Soden_method:
                    dic_kern['dlogq1k'] = MV.masked_where(dic_kern['lev_4d']<=dic_kern['p_tropopause'],dic_kern['dlogq1k'])
                    dic_kern['dlogq'] = MV.masked_where(dic_kern['lev_4d']<=dic_kern['p_tropopause'],dic_kern['dlogq']) 
                    dic_kern['dlogq2'] = MV.masked_where(dic_kern['lev_4d']<=dic_kern['p_tropopause'],dic_kern['dlogq2']) 
        
                    dic_mon['wv_lw_kernel_cld_monthly'] = MV.masked_where(dic_kern['lev_4d']<=dic_kern['p_tropopause'],dic_mon['wv_lw_kernel_cld_monthly'])
                    dic_mon['wv_sw_kernel_cld_monthly'] = MV.masked_where(dic_kern['lev_4d']<=dic_kern['p_tropopause'],dic_mon['wv_sw_kernel_cld_monthly'])
                    dic_mon['wv_lw_kernel_clr_monthly'] = MV.masked_where(dic_kern['lev_4d']<=dic_kern['p_tropopause'],dic_mon['wv_lw_kernel_clr_monthly'])
                    dic_mon['wv_sw_kernel_clr_monthly'] = MV.masked_where(dic_kern['lev_4d']<=dic_kern['p_tropopause'],dic_mon['wv_sw_kernel_clr_monthly'])
    
    
                dic_kern['dlogq'][dic_kern['dlogq'].mask]=0
                dic_kern['dLW_q'] = dic_mon['wv_lw_kernel_cld_monthly'] / dic_kern['dlogq1k'] * dic_kern['dlogq'] #hus_ano_grd
                dic_kern['dSW_q'] = dic_mon['wv_sw_kernel_cld_monthly'] / dic_kern['dlogq1k'] * dic_kern['dlogq'] # hus_ano_grd
                dic_kern['dLW_q'].setAxisList(AXL4d)
                dic_kern['dSW_q'].setAxisList(AXL4d)
    
                dic_kern['dLW_q_clr'] = dic_mon['wv_lw_kernel_clr_monthly'] / dic_kern['dlogq1k'] * dic_kern['dlogq'] #hus_ano_grd
                dic_kern['dSW_q_clr'] = dic_mon['wv_sw_kernel_clr_monthly'] / dic_kern['dlogq1k'] * dic_kern['dlogq'] #hus_ano_grd
                dic_kern['dLW_q_clr'].setAxisList(AXL4d)
                dic_kern['dSW_q_clr'].setAxisList(AXL4d)
    
   
                ## vertical integral (sum)
                dic_kern['dLW_q_psum'] = MV.sum(dic_kern['dLW_q'] * dic_kern['dp_4d']/100., axis=1)
                dic_kern['dSW_q_psum'] = MV.sum(dic_kern['dSW_q'] * dic_kern['dp_4d']/100., axis=1)
    
                dic_kern['dLW_q_clr_psum'] = MV.sum(dic_kern['dLW_q_clr'] * dic_kern['dp_4d']/100., axis=1)
                dic_kern['dSW_q_clr_psum'] = MV.sum(dic_kern['dSW_q_clr'] * dic_kern['dp_4d']/100., axis=1)
    
                dic_kern['dLW_q_psum'].setAxisList(AXL3d)
                dic_kern['dSW_q_psum'].setAxisList(AXL3d)
                dic_kern['dLW_q_clr_psum'].setAxisList(AXL3d)
                dic_kern['dSW_q_clr_psum'].setAxisList(AXL3d)
    
                names = ['dSW_q_psum','dSW_q_clr_psum','dLW_q_psum','dLW_q_clr_psum']
                for name in names:
                    DATA = dic_kern[name]
                    result = get_fdbk(DATA,AXL2d,ts_ano_grd_avg,ts_ano_grd_ann,name,exp_cntl,dic_gfdbk,dic_feedback)
 
                dic_gfdbk['q_gfdbk'] = dic_gfdbk['dSW_q_psum_gfdbk'] + dic_gfdbk['dLW_q_psum_gfdbk']
                dic_gfdbk['q_clr_gfdbk'] = dic_gfdbk['dSW_q_clr_psum_gfdbk'] + dic_gfdbk['dLW_q_clr_psum_gfdbk']
    
                dic_feedback['q_feedback'] = dic_feedback['dSW_q_psum_feedback'] + dic_feedback['dLW_q_psum_feedback']
                dic_feedback['q_clr_feedback'] = dic_feedback['dSW_q_clr_psum_feedback'] + dic_feedback['dLW_q_clr_psum_feedback']

                # June 30, 2020 add: q_sw and q_lw separately
                dic_feedback['q_sw_feedback'] = dic_feedback['dSW_q_psum_feedback']
                dic_feedback['q_lw_feedback'] = dic_feedback['dLW_q_psum_feedback']
                dic_feedback['q_clr_sw_feedback'] = dic_feedback['dSW_q_clr_psum_feedback']
                dic_feedback['q_clr_lw_feedback'] = dic_feedback['dLW_q_clr_psum_feedback']

                dic_gfdbk['q_sw_gfdbk'] = dic_gfdbk['dSW_q_psum_gfdbk']
                dic_gfdbk['q_lw_gfdbk'] = dic_gfdbk['dLW_q_psum_gfdbk']
                dic_gfdbk['q_clr_sw_gfdbk'] = dic_gfdbk['dSW_q_clr_psum_gfdbk']
                dic_gfdbk['q_clr_lw_gfdbk'] = dic_gfdbk['dLW_q_clr_psum_gfdbk']
               
                print("Water vapor feedback: ",dic_feedback['q_feedback'],"W/m2/K")
                print("clr Water vapor feedback: ",dic_feedback['q_clr_feedback'],"W/m2/K")
    
                print('------ Hi, summary all feedbacks except for cloud feedback------------')
                print('Planck feedback: ',dic_feedback['planck_feedback'],'W/m2/K')
                print('Lapse rate feedback: ',dic_feedback['lapserate_feedback'],'W/m2/K')
                print('Lapse rate + Planck feedback: ',dic_feedback['lapserate_feedback']+dic_feedback['planck_feedback'],'W/m2/K')
                print('Temperature feedback: ',dic_feedback['t_feedback'],'W/m2/K')
                print('Water vapor feedback: ',dic_feedback['q_feedback'],'W/m2/K')
                print("Surface albedo feedback: ",dic_feedback['alb_feedback'], "W/m2/K")
    
                print('--------clear-sky component----------------------------------------------')
                print('clr Planck feedback: ',dic_feedback['planck_clr_feedback'],'W/m2/K')
                print('clr Lapse rate feedback: ',dic_feedback['lapserate_clr_feedback'],'W/m2/K')
                print('clr Lapse rate + clr Planck feedback: ',dic_feedback['lapserate_clr_feedback']+dic_feedback['planck_clr_feedback'],'W/m2/K')
                print('clr Temperature feedback: ',dic_feedback['t_clr_feedback'],'W/m2/K')
                print('clr Water vapor feedback: ',dic_feedback['q_clr_feedback'],'W/m2/K')
                print("clr Surface albedo feedback: ",dic_feedback['alb_clr_feedback'], "W/m2/K")
    
                sums_sw_feedback = dic_feedback['dSW_q_clr_psum_feedback'] + dic_feedback['dSW_alb_clr_feedback']
                sums_lw_feedback = dic_feedback['dLW_q_clr_psum_feedback'] + dic_feedback['dLW_ts_clr_feedback'] + dic_feedback['dLW_ta_clr_psum_feedback']
    
                print('sum of all clear-sky sw feedbacks',sums_sw_feedback,'W/m2/K')
                print('sum of all clear-sky lw feedbacks',sums_lw_feedback,'W/m2/K')
    
                # get the TOA anomalies from direct model data
                dic_kern['anos_sw'] = dic_grd['rsntcs_ano_grd']
                dic_kern['anos_lw'] = -1.*dic_grd['rlutcs_ano_grd']
                dic_kern['anos_sw'].setAxisList(AXL3d)
                dic_kern['anos_lw'].setAxisList(AXL3d)
    
                names = ['anos_sw','anos_lw']
                for name in names:
                   DATA = dic_kern[name]
                   result = get_fdbk(DATA,AXL2d,ts_ano_grd_avg,ts_ano_grd_ann,name,exp_cntl,dic_gfdbk,dic_feedback)
                    
                print('TOA direct clear-sky sw radiation feedback',dic_feedback['anos_sw_feedback'],'W/m2/K')
                print('TOA direct clear-sky lw radiation feedback',dic_feedback['anos_lw_feedback'],'W/m2/K')
    
                ## plot lat-lon figure as Shell et al. (2006) Figure2
                # get lat-lon data
                dic_kern['dLW_clr_sum'] = dic_kern['dLW_t_clr_psum'] + dic_kern['dLW_q_clr_psum']
                dic_kern['dSW_clr_sum'] = dic_kern['dSW_alb_clr'] + dic_kern['dSW_q_clr_psum']
                dic_kern['dLW_clr_sum'].setAxisList(AXL3d)
                dic_kern['dSW_clr_sum'].setAxisList(AXL3d)
    
                dic_kern['dLW_clr_dir'] = -1. * dic_grd['rlutcs_ano_grd']
                dic_kern['dSW_clr_dir'] = dic_grd['rsntcs_ano_grd']
                dic_kern['dLW_clr_dir'].setAxisList(AXL3d)
                dic_kern['dSW_clr_dir'].setAxisList(AXL3d)
    
                # difference between direct and kernel-calculation
                dic_kern['dLW_clr_dir_sum'] = dic_kern['dLW_clr_dir'] - dic_kern['dLW_clr_sum']
                dic_kern['dSW_clr_dir_sum'] = dic_kern['dSW_clr_dir'] - dic_kern['dSW_clr_sum']
                dic_kern['dLW_clr_dir_sum'].setAxisList(AXL3d)
                dic_kern['dSW_clr_dir_sum'].setAxisList(AXL3d)
    
                names = ['dLW_clr_dir','dSW_clr_dir','dLW_clr_sum','dSW_clr_sum','dLW_clr_dir_sum','dSW_clr_dir_sum']
                for name in names:
                   DATA = dic_kern[name]
                   result = get_fdbk(DATA,AXL2d,ts_ano_grd_avg,ts_ano_grd_ann,name,exp_cntl,dic_gfdbk,dic_feedback)
   
                # get the vertically-averaged kernel for later plotting
                dic_mon['t_kernel_cld_psum'] = MV.sum(dic_mon['t_kernel_cld_monthly']*dic_kern['dp_4d']/100,axis=1)
                dic_mon['t_kernel_cld_psum'].setAxisList(AXL3d)
                dic_mon['wv_lw_kernel_cld_psum'] = MV.sum(dic_mon['wv_lw_kernel_cld_monthly']*dic_kern['dp_4d']/100,axis=1)
                dic_mon['wv_sw_kernel_cld_psum'] = MV.sum(dic_mon['wv_sw_kernel_cld_monthly']*dic_kern['dp_4d']/100,axis=1)
                dic_mon['wv_kernel_cld_psum'] = dic_mon['wv_lw_kernel_cld_psum'] + dic_mon['wv_sw_kernel_cld_psum']
                dic_mon['wv_kernel_cld_psum'].setAxisList(AXL3d)
                dic_mon['ts_kernel_cld_psum'] = dic_mon['ts_kernel_cld_monthly']
    
                # test -------------------------------------------------------------------------
                ## zonal-mean 
                fig = plt.figure(figsize=(18,9))
                plt.suptitle(case_stamp,fontsize=fh,y=0.97)
    
                names1 = ['dLW_clr_sum_gfdbk','dSW_clr_sum_gfdbk']
                names2 = ['dLW_clr_dir_gfdbk','dSW_clr_dir_gfdbk']

                names1_out = ['LW_clrsky_KernelSum','SW_clrsky_KernelSum']
                names2_out = ['LW_clrsky_Direct','SW_clrsky_Direct']
                
                for n,name in enumerate(names1):
                    DATA1 = MV.average(dic_gfdbk[names1[n]],axis=1)
                    DATA2 = MV.average(dic_gfdbk[names2[n]],axis=1)
                    ax2 = fig.add_subplot(1,2,n+1)
                    plt.plot(lat,np.array(DATA1),label=names1_out[n],color='red')
                    plt.plot(lat,np.array(DATA2),label=names2_out[n],color='blue')
                    plt.plot(lat,np.array(DATA2)*1.15,linestyle='--',color='black')
                    plt.plot(lat,np.array(DATA2)*0.85,linestyle='--',color='black')
                    pl.title(names1_out[n]+' vs '+names2_out[n],fontsize=fh)
                    ax2.tick_params(labelsize=fh)
                    plt.legend(fontsize=fh)
                    plt.xlabel('Latitude',fontsize=fh)
                    plt.ylabel('Feedback [W/m2/K]', fontsize=fh)
    
                plt.tight_layout(rect=[0, 0.03, 1, 0.95])
                plt.savefig(figdir+'zonalmean-clear-sky-linear-test-'+exp_cntl+'-'+used_models[imod]+'-'+case_stamp+'.png',bbox_inches='tight')
                ### test --------------------------------------------------------------------------
    
                #################################################################################
                #### Adjusted cloud feedback calculation
                #################################################################################
                # calculate cloud masking 
                dic_kern['dLW_t_psum_mask'] = dic_kern['dLW_t_clr_psum'] - dic_kern['dLW_t_psum']
                dic_kern['dLW_q_psum_mask'] = dic_kern['dLW_q_clr_psum'] - dic_kern['dLW_q_psum']
                dic_kern['dSW_q_psum_mask'] = dic_kern['dSW_q_clr_psum'] - dic_kern['dSW_q_psum']
                dic_kern['dSW_alb_mask'] = dic_kern['dSW_alb_clr'] - dic_kern['dSW_alb']
                dic_kern['dLW_t_psum_mask'].setAxisList(AXL3d)
                dic_kern['dLW_q_psum_mask'].setAxisList(AXL3d)
                dic_kern['dSW_q_psum_mask'].setAxisList(AXL3d)
                dic_kern['dSW_alb_mask'].setAxisList(AXL3d)
    
                # adjusted CRE
                dic_kern['dLW_adj'] = dic_kern['dLW_t_psum_mask'] + dic_kern['dLW_q_psum_mask']
                dic_kern['dSW_adj'] = dic_kern['dSW_q_psum_mask'] + dic_kern['dSW_alb_mask']
                dic_kern['net_adj'] = dic_kern['dLW_adj'] + dic_kern['dSW_adj']
                dic_kern['dLW_adj'].setAxisList(AXL3d)
                dic_kern['dSW_adj'].setAxisList(AXL3d)
                dic_kern['net_adj'].setAxisList(AXL3d)
    
                dic_kern['SWCRE_ano_grd_adj'] = dic_grd['SWCRE_ano_grd'] + dic_kern['dSW_adj']
                dic_kern['LWCRE_ano_grd_adj'] = dic_grd['LWCRE_ano_grd'] + dic_kern['dLW_adj']
                dic_kern['netCRE_ano_grd_adj'] = dic_grd['netCRE_ano_grd'] + dic_kern['dSW_adj'] + dic_kern['dLW_adj']
                dic_kern['SWCRE_ano_grd_adj'].setAxisList(AXL3d)
                dic_kern['LWCRE_ano_grd_adj'].setAxisList(AXL3d)
                dic_kern['netCRE_ano_grd_adj'].setAxisList(AXL3d)
    
                ############## get cloudy residual term
                dic_kern['dLW_cld_sum'] = dic_kern['dLW_t_psum'] + dic_kern['dLW_q_psum'] + dic_kern['LWCRE_ano_grd_adj']
                dic_kern['dSW_cld_sum'] = dic_kern['dSW_alb'] + dic_kern['dSW_q_psum'] + dic_kern['SWCRE_ano_grd_adj']
                dic_kern['dLW_cld_sum'].setAxisList(AXL3d)
                dic_kern['dSW_cld_sum'].setAxisList(AXL3d)
    
                dic_kern['dLW_cld_dir'] = -1. * dic_grd['rlut_ano_grd']
                dic_kern['dSW_cld_dir'] = dic_grd['rsnt_ano_grd']
                dic_kern['dLW_cld_dir'].setAxisList(AXL3d)
                dic_kern['dSW_cld_dir'].setAxisList(AXL3d)
    
                # difference between direct and kernel-calculation
                dic_kern['dLW_cld_dir_sum'] = dic_kern['dLW_cld_dir'] - dic_kern['dLW_cld_sum']
                dic_kern['dSW_cld_dir_sum'] = dic_kern['dSW_cld_dir'] - dic_kern['dSW_cld_sum']
                dic_kern['dLW_cld_dir_sum'].setAxisList(AXL3d)
                dic_kern['dSW_cld_dir_sum'].setAxisList(AXL3d)
    
                dic_kern['dnet_cld_dir_sum'] = dic_kern['dLW_cld_dir_sum'] + dic_kern['dSW_cld_dir_sum']
                dic_kern['dnet_cld_dir_sum'].setAxisList(AXL3d)
    

                names = ['dLW_cld_dir','dSW_cld_dir','dLW_cld_sum','dSW_cld_sum','dLW_cld_dir_sum','dSW_cld_dir_sum']
                for name in names:
                    DATA = dic_kern[name]
                    result = get_fdbk(DATA,AXL2d,ts_ano_grd_avg,ts_ano_grd_ann,name,exp_cntl,dic_gfdbk,dic_feedback)
 

                names = ['dLW_t_psum_mask','dSW_alb_mask','dLW_q_psum_mask','dSW_q_psum_mask','dLW_adj','dSW_adj','SWCRE_ano_grd','LWCRE_ano_grd','netCRE_ano_grd','SWCRE_ano_grd_adj','LWCRE_ano_grd_adj','net_adj','netCRE_ano_grd_adj','dLW_cld_sum','dSW_cld_sum','dLW_cld_dir','dSW_cld_dir','dLW_cld_dir_sum','dSW_cld_dir_sum','dnet_cld_dir_sum']
                for name in names:
                    if name in ['SWCRE_ano_grd','LWCRE_ano_grd','netCRE_ano_grd']:
                        DATA = dic_grd[name]
                    else:
                        DATA = dic_kern[name]
                    result = get_fdbk(DATA,AXL2d,ts_ano_grd_avg,ts_ano_grd_ann,name,exp_cntl,dic_gfdbk,dic_feedback)
   
                ### test --------------------------------------------------------------------------------
                # plot spatial distribution of other non-cloud feedbacks and adjusted CRE feedback
                # added by YiQin on Aug 6, 2020
                fig = plt.figure(figsize=(15,12))
                plt.suptitle(case_stamp,fontsize=fh,y=0.97)
                
                plot_names = ['planck','lapserate','q','alb',\
                'netCRE_ano_grd_adj','SWCRE_ano_grd_adj','LWCRE_ano_grd_adj']
                plot_names_out = ['Planck','LR','WV','ALB','netCRE_adj','SWCRE_adj','LWCRE_adj']

                # lat-lon
                for n,name in enumerate(plot_names):
                    DATA = dic_gfdbk[name+'_gfdbk']
    
                    bounds = np.arange(-10,12,2)/4.
                    cmap = pl.cm.RdBu_r
                    bounds2 = np.append(np.append(-500,bounds),500)
                    norm = mpl.colors.BoundaryNorm(bounds2,cmap.N)
    
                    ax1 = fig.add_subplot(4,2,n+1,projection=ccrs.Robinson(central_longitude=180.))
                    im1 = ax1.contourf(lon,lat,DATA,bounds,transform=ccrs.PlateCarree(),cmap=cmap,norm=norm,extend='both')
                    ax1.coastlines()
                    ax1.set_global()
                    DATA.setAxisList(AXL2d)
                    avgDATA = cdutil.averager(DATA,axis='xy',weights='weighted')
                    pl.title(plot_names_out[n] + ' ['+str(np.round(avgDATA,3))+']',fontsize=fh)
                    cb = plt.colorbar(im1,orientation='vertical',drawedges=True,ticks=bounds)
    
                    cb.set_label('W/m$^2$/K')
                plt.tight_layout(rect=[0, 0.03, 1, 0.95])
                plt.savefig(figdir+'lat-lon-NonCloudFeedback-adjCRE-'+exp_cntl+'-'+used_models[imod]+'-'+case_stamp+'.png',bbox_inches='tight')
                ### test --------------------------------------------------------------------------------

                # ------------------------------------------------------------
                # output global mean values 
                # ------------------------------------------------------------
                # June 6, 2020: add clear-sky feedback output
                # June 20, 2020: add water vapor SW/LW feedback
                outputs = ['t_feedback','planck_feedback','lapserate_feedback','q_feedback','alb_feedback','dLW_adj_feedback','dSW_adj_feedback','net_adj_feedback','SWCRE_ano_grd_feedback','LWCRE_ano_grd_feedback','netCRE_ano_grd_feedback','SWCRE_ano_grd_adj_feedback', 'LWCRE_ano_grd_adj_feedback', 'netCRE_ano_grd_adj_feedback','dLW_cld_dir_sum_feedback','dSW_cld_dir_sum_feedback','dnet_cld_dir_sum_feedback','t_clr_feedback','planck_clr_feedback','lapserate_clr_feedback','q_clr_feedback','alb_clr_feedback','q_sw_feedback','q_lw_feedback','q_clr_sw_feedback','q_clr_lw_feedback']
    
                # output variables to csv file
                df_all = pd.DataFrame([dic_feedback.get(key) for key in outputs],index=['T','Planck','LR','WV','ALB','dLW_adj','dSW_adj','dnet_adj','SWCRE','LWCRE','netCRE','SWCRE_adj','LWCRE_adj','netCRE_adj','LW_resd','SW_resd','net_resd','T_clr','Planck_clr','LR_clr','WV_clr','ALB_clr','WV_SW','WV_LW','WV_clr_SW','WV_clr_LW'],columns=[used_models[imod]])

                if get_tropo_Mark_method:
                    df_all.to_csv(outdir+'FDBK_'+exp_cntl+'_'+used_models[imod]+'_'+case_stamp+'_VaryingTrop.csv')
                else:
                    df_all.to_csv(outdir+'FDBK_'+exp_cntl+'_'+used_models[imod]+'_'+case_stamp+'.csv')
    
                # ------------------------------------------------------------
                # output spatial pattern
                # ------------------------------------------------------------
                outputs = ['t_gfdbk','planck_gfdbk','lapserate_gfdbk','q_gfdbk','alb_gfdbk','dLW_adj_gfdbk','dSW_adj_gfdbk','net_adj_gfdbk','SWCRE_ano_grd_gfdbk','LWCRE_ano_grd_gfdbk','netCRE_ano_grd_gfdbk','SWCRE_ano_grd_adj_gfdbk','LWCRE_ano_grd_adj_gfdbk','netCRE_ano_grd_adj_gfdbk','dLW_cld_dir_sum_gfdbk','dSW_cld_dir_sum_gfdbk','dnet_cld_dir_sum_gfdbk','t_clr_gfdbk','planck_clr_gfdbk','lapserate_clr_gfdbk','q_clr_gfdbk','alb_clr_gfdbk','q_sw_gfdbk','q_lw_gfdbk','q_clr_sw_gfdbk','q_clr_lw_gfdbk']
    
                if get_tropo_Mark_method:
                    out1 = cdms.open(outdir+'lat-lon-gfdbk-'+exp_cntl+'-'+used_models[imod]+'-'+case_stamp+'-VaryingTrop.nc','w')
                else:
                    out1 = cdms.open(outdir+'lat-lon-gfdbk-'+exp_cntl+'-'+used_models[imod]+'-'+case_stamp+'.nc','w')
    
                for output in outputs:
                    print(output)
                    DATA = dic_gfdbk[output]
                    DATA.id = output
                    out1.write(DATA)
                out1.close()

