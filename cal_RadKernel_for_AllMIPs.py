
## this script is used to calculate radiative kernel following Soden et al. and Shell et al.
## some trival but important things need to be taken care:
## 1. albedo should be in percent (%) unit;
## 2. pay attention to the masking, which should be level < p_tropopause;

## created: March 12, 2020 by Yi Qin
## modified: July 7 -- including relative humidity feedback following Held and Shell (2012)
## Oct 10, 2020 -- change ts to tas


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
import DefinedFunction as RD
import genutil
import numpy.ma as ma
from genutil import statistics
from scipy.interpolate import interp1d

import CMIP6_utils as CU
sys.path.append('/work/share/')
import PMC_utils
sys.path.append('/work/zelinka1/scripts/')
sys.path.append('/work/caldwell19/cdat/kernels/')
import MDZ_utils as MU

import psutil 
import gc
import copy

# Oct 7, 2020: move read data as a function
def read_data_e3sm(var,var2d,var3d,direc_data1,direc_data2,exp1,exp2,yrS_4d,monS_2d,yrE_4d,monE_2d,dic_invar):

    for svar in var:
        if svar in var:
            print(" =================== we are processing E3SM amip data", svar, " locally ====================")

            f1 = cdms.open(direc_data1+svar+'_'+exp1+'_'+yrS_4d+monS_2d+'-'+yrE_4d+monE_2d+'.nc','r')
            dic_invar[svar+'_pi'] = f1(svar)[:nyears*12,:,:]
            f1.close()
 
            f2 = cdms.open(direc_data2+svar+'_'+exp2+'_'+yrS_4d+monS_2d+'-'+yrE_4d+monE_2d+'.nc','r')
            dic_invar[svar+'_ab'] = f2(svar)[:nyears*12,:,:]
            f2.close()
 
            # reverse lev direction
            if svar in var3d:
                dic_invar[svar+'_pi'] = dic_invar[svar+'_pi'][:,::-1,:,:]
                dic_invar[svar+'_ab'] = dic_invar[svar+'_ab'][:,::-1,:,:]
 
            stop_here = False
        else:
            stop_here = True

        if stop_here: ### we don't get enough data to do further processing. just skip out this loop.
            print('stop_here is', stop_here)
            continue
   
    ### get SWCF, LWCF and their anomalies
    dic_invar['SWCRE_ano']  = dic_invar['rsutcs_ano'] - dic_invar['rsut_ano']
    dic_invar['LWCRE_ano']  = dic_invar['rlutcs_ano'] - dic_invar['rlut_ano']
    dic_invar['netCRE_ano'] = dic_invar['SWCRE_ano']  + dic_invar['LWCRE_ano']

    delterm = ['rsutcs_pi','rsutcs_ab','rsut_pi','rsut_ab','rlutcs_pi','rlutcs_ab','rlut_pi','rlut_ab','ps_ab','ps_ano']
    dic_invar = delete_vars(delterm,dic_invar)

    AXL2d = dic_invar['rsus_pi'].getAxisList()

    ## get albedo
    dic_invar['alb_pi'] = dic_invar['rsus_pi']/dic_invar['rsds_pi'] * 100.
    dic_invar['alb_ab'] = dic_invar['rsus_ab']/dic_invar['rsds_ab'] *100.
    dic_invar['alb_pi'] = MV.masked_outside(dic_invar['alb_pi'], 0.0, 100.)
    dic_invar['alb_ab'] = MV.masked_outside(dic_invar['alb_ab'], 0.0, 100.)
    dic_invar['alb_ano'] = dic_invar['alb_ab'] - dic_invar['alb_pi']
    dic_invar['alb_pi'].setAxisList(AXL2d)
    dic_invar['alb_ab'].setAxisList(AXL2d)
    dic_invar['alb_ano'].setAxisList(AXL2d)

    delterm = ['rsus_pi','rsds_pi','rsus_ab','rsds_ab','alb_pi','alb_ab','rsdt_pi']
    dic_invar = delete_vars(delterm,dic_invar)
   
    # Oct 08, 2020: we don't need alb_clr_ano. I misunderstood the clear-sky albedo feedback.
#    dic_invar['alb_clr_pi'] = dic_invar['rsuscs_pi']/dic_invar['rsdscs_pi'] * 100.
#    dic_invar['alb_clr_ab'] = dic_invar['rsuscs_ab']/dic_invar['rsdscs_ab'] *100.
#    dic_invar['alb_clr_pi'] = MV.masked_outside(dic_invar['alb_clr_pi'], 0.0, 100.)
#    dic_invar['alb_clr_ab'] = MV.masked_outside(dic_invar['alb_clr_ab'], 0.0, 100.)
#    dic_invar['alb_clr_ano'] = dic_invar['alb_clr_ab'] - dic_invar['alb_clr_pi'] 
#    dic_invar['alb_clr_pi'].setAxisList(AXL2d)
#    dic_invar['alb_clr_ab'].setAxisList(AXL2d)
#    dic_invar['alb_clr_ano'].setAxisList(AXL2d)

#    delterm = ['rsuscs_pi','rsdscs_pi','rsuscs_ab','rsdscs_ab','alb_clr_pi','alb_clr_ab','ts_pi','ts_ab','rsus_ano','rsds_ano','rsuscs_ano','rsdscs_ano']
    delterm = ['rsuscs_pi','rsdscs_pi','rsuscs_ab','rsdscs_ab','ts_pi','ts_ab','rsus_ano','rsds_ano','rsuscs_ano','rsdscs_ano','tas_pi','tas_ano']

    dic_invar = delete_vars(delterm,dic_invar)

    return dic_invar
 

def read_data(var,var2d,var3d,nyears,direc_data,project_new,exp_new,comp,freq,used_inst,used_model,used_ripf,project_cntl,exp_cntl,phase,startyear,dic_invar):

    # read 2D variables
    for svar in var:
        print(svar)

#        if svar in var2d:
        if svar in var:
            sub_name = 'Reading '+svar
            print_memory_status(sub_name)

            output = RD.get_data(svar,nyears,direc_data,project_new,exp_new,comp,freq,used_inst,used_model,used_ripf,project_cntl,exp_cntl,phase,startyear)
            sub_name = 'Reading '+svar
            print_memory_status(sub_name)

        else: 
            output = RD.get_data_cloud3d(usedin,svar,nyears,direc_data,project_new,exp_new,comp,freq,used_inst,used_model,used_ripf,project_cntl,exp_cntl,phase,startyear)
            sub_name = 'Reading '+svar
            print_memory_status(sub_name)

        if output == -1:
            print("============Warning==============This model",used_model,"does not have all suitable data. So let us skip it now")
            stop_here = True
            break

        dic_invar[svar+'_pi'] = output[0]
        dic_invar[svar+'_ab'] = output[1]
        stop_here = False

        dic_invar[svar+'_ano'] = dic_invar[svar+'_ab'] - dic_invar[svar+'_pi']
        dic_invar[svar+'_ano'].setAxisList(dic_invar[svar+'_pi'].getAxisList())

        if stop_here: ### we don't get enough data to do further processing. just skip out this loop.
            print('stop_here is', stop_here)
            continue

    print('output getrefcount is',sys.getrefcount(output))
   
    ### get SWCF, LWCF and their anomalies
    dic_invar['SWCRE_ano']  = dic_invar['rsutcs_ano'] - dic_invar['rsut_ano']
    dic_invar['LWCRE_ano']  = dic_invar['rlutcs_ano'] - dic_invar['rlut_ano']
    dic_invar['netCRE_ano'] = dic_invar['SWCRE_ano']  + dic_invar['LWCRE_ano']

    delterm = ['rsutcs_pi','rsutcs_ab','rsut_pi','rsut_ab','rlutcs_pi','rlutcs_ab','rlut_pi','rlut_ab','ps_ab','ps_ano']
    dic_invar = delete_vars(delterm,dic_invar)

    AXL2d = dic_invar['rsus_pi'].getAxisList()

    ## get albedo
    dic_invar['alb_pi'] = dic_invar['rsus_pi']/dic_invar['rsds_pi'] * 100.
    dic_invar['alb_ab'] = dic_invar['rsus_ab']/dic_invar['rsds_ab'] *100.
    dic_invar['alb_pi'] = MV.masked_outside(dic_invar['alb_pi'], 0.0, 100.)
    dic_invar['alb_ab'] = MV.masked_outside(dic_invar['alb_ab'], 0.0, 100.)
    dic_invar['alb_ano'] = dic_invar['alb_ab'] - dic_invar['alb_pi']
    dic_invar['alb_pi'].setAxisList(AXL2d)
    dic_invar['alb_ab'].setAxisList(AXL2d)
    dic_invar['alb_ano'].setAxisList(AXL2d)

    delterm = ['rsus_pi','rsds_pi','rsus_ab','rsds_ab','alb_pi','alb_ab','rsdt_pi']
    dic_invar = delete_vars(delterm,dic_invar)
   
    # Oct 08, 2020: we don't need clear-sky albedo for clear-sky albedo feedback!
#    dic_invar['alb_clr_pi'] = dic_invar['rsuscs_pi']/dic_invar['rsdscs_pi'] * 100.
#    dic_invar['alb_clr_ab'] = dic_invar['rsuscs_ab']/dic_invar['rsdscs_ab'] *100.
#    dic_invar['alb_clr_pi'] = MV.masked_outside(dic_invar['alb_clr_pi'], 0.0, 100.)
#    dic_invar['alb_clr_ab'] = MV.masked_outside(dic_invar['alb_clr_ab'], 0.0, 100.)
#    dic_invar['alb_clr_ano'] = dic_invar['alb_clr_ab'] - dic_invar['alb_clr_pi'] 
#    dic_invar['alb_clr_pi'].setAxisList(AXL2d)
#    dic_invar['alb_clr_ab'].setAxisList(AXL2d)
#    dic_invar['alb_clr_ano'].setAxisList(AXL2d)

#    delterm = ['rsuscs_pi','rsdscs_pi','rsuscs_ab','rsdscs_ab','alb_clr_pi','alb_clr_ab','ts_pi','ts_ab','rsus_ano','rsds_ano','rsuscs_ano','rsdscs_ano']
    delterm = ['rsuscs_pi','rsdscs_pi','rsuscs_ab','rsdscs_ab','ts_pi','ts_ab','rsus_ano','rsds_ano','rsuscs_ano','rsdscs_ano','tas_pi','tas_ab']

    dic_invar = delete_vars(delterm,dic_invar)

    return dic_invar
 

# Oct 4, 2020: chagne Mark's qsat method as a function
def qsat_blend_Mark(avgta,lev_4d):
    wsl=MU.qsat_water(avgta,lev_4d*100.)
    wsl_plus1=MU.qsat_water(avgta+1,lev_4d*100.)
    wsi=PMC_utils.qsat_ice(avgta,lev_4d*100.)
    wsi_plus1=PMC_utils.qsat_ice(avgta+1,lev_4d*100.)
    qsl=wsl/(1+wsl) # convert from mixing ratio (kg/kg) to specific humidity
    qsl_plus1=wsl_plus1/(1+wsl_plus1)
    qsi=wsi/(1+wsi) # convert from mixing ratio (kg/kg) to specific humidity
    qsi_plus1=wsi_plus1/(1+wsi_plus1)   
    del wsl,wsl_plus1,wsi,wsi_plus1

    # Compute blend of qsi and qsl between -40 and 0
    blend = (avgta-233)*qsl/40 + (273-avgta)*qsi/40
    blend_plus1 = (avgta-233)*qsl_plus1/40 + (273-avgta)*qsi_plus1/40

    qs0 = np.where((avgta>233) & (avgta<273), blend, qsi)#[0]
    qs1 = np.where((avgta>233) & (avgta<273), blend_plus1, qsi_plus1)#[0]
    
    qs0 = np.where(avgta >= 273, qsl, qs0)#[0]
    qs1 = np.where(avgta >= 273, qsl_plus1, qs1)#[0]
    del blend, blend_plus1,qsi,qsi_plus1,qsl,qsl_plus1

    qs0 = np.float32(qs0)
    qs1 = np.float32(qs1)

    return qs0, qs1

# Oct 4, 2020: get feedback
#@profile
def get_fdbk(invar1,invar2,outvar,dic_invar,AXL4d, AXL3d, dp_use_method='SPC'):
    for ivar,svar in enumerate(invar1):
        ovar = outvar[ivar]
        print('get_fdbk: invar1, invar2, outvar: ',invar1[ivar],invar2[ivar],outvar[ivar])
        dic_invar[ovar] = dic_invar[svar] * dic_invar[invar2[ivar]]
        dic_invar[ovar].setAxisList(AXL4d)

        # vertical integral (sum) -- first get weights: dp
        if dp_use_method == 'DynPs' or dp_use_method =='SPC':
            dic_invar[ovar+'_psum'] = VertSum(dic_invar[ovar],dic_invar['dp_4d'])
            dic_invar[ovar+'_psum'].setAxisList(AXL3d)
        else:
            dic_invar[ovar+'_psum'] = MV.sum(dic_invar[ovar]*dic_invar['dp_4d']/100., axis=1)
            dic_invar[ovar+'_psum'].setAxisList(AXL3d)

    return dic_invar

# Oct 3, 2020: delete vars
#@profile
def delete_vars(delterm,dic_invar):
    sub_name = 'before delete item'
    print_memory_status(sub_name)

    for svar in delterm:
        if svar in dic_invar.keys():
            print(svar,'dic_invar[svar] getrefcount is', sys.getrefcount(dic_invar[svar]))
            print('delete ',svar)
            del dic_invar[svar]
            gc.collect()
    print('after delete, dic_invar.keys() are ',dic_invar.keys())

    sub_name = 'after delete item'
    print_memory_status(sub_name)

    return dic_invar
    
 

# Oct 2, 2020: vertical interpolation with loop to test the memory usage
#@profile
def Vert_RegridStd_Loop(DATA,stdlevs,std_plev):

    rawlevs = DATA.getLevel()[:]
    if DATA.getLevel().units == 'Pa':
        rawlevs = rawlevs/100.

    DATA2 = np.where(DATA.mask,np.nan,DATA)
    nlat = DATA2.shape[2]
    nlon = DATA2.shape[3]
    ntime = DATA2.shape[0]

    newdata = np.empty((DATA2.shape[0],len(stdlevs),nlat,nlon))

    for itime in range(ntime):
        for ilat,ilon in zip(range(nlat),range(nlon)):
            newdata[itime,:,ilat,ilon] = np.interp(stdlevs.data,rawlevs,DATA2[itime,:,ilat,ilon])

    newdata = np.float32(newdata)
    newdata=MV.masked_where(np.isnan(newdata),newdata)
    newdata = cdms.asVariable(newdata)

    newdata.setAxis(0,DATA.getAxis(0))
    newdata.setAxis(1,std_plev)
    newdata.setAxis(2,DATA.getAxis(2))
    newdata.setAxis(3,DATA.getAxis(3))
    del DATA2,rawlevs,nlat,nlon,ntime

    return newdata1

 
# Oct 2, 2020: move vertical interpolation as the function
#@profile
def Vert_RegridStd(DATA,stdlevs,std_plev):
    rawlevs = DATA.getLevel()[:]
    if DATA.getLevel().units == 'Pa':
        rawlevs = rawlevs/100.
    
    # "steal" from Mark's script: load_TOA_kernels.py -- March 12
    # Oct 3, 2020: MV.where has more memory leak than np.where. so, use np.where here.
#    DATA2=MV.where(DATA.mask,np.nan,DATA)
    DATA2=np.where(DATA.mask,np.nan,DATA)

    f = interp1d(rawlevs,DATA2,axis=1,copy=False)
    newdata = f(stdlevs.data)

    newdata = np.float32(newdata)
    newdata = cdms.asVariable(newdata)
    newdata=np.ma.masked_where(np.isnan(newdata),newdata,copy=False)

    newdata.setAxis(0,DATA.getAxis(0))
    newdata.setAxis(1,std_plev)
    newdata.setAxis(2,DATA.getAxis(2))
    newdata.setAxis(3,DATA.getAxis(3))

    del f,DATA2,rawlevs

    return newdata


# Oct 2, 2020: add memory diagnostics functions in MB and % separately
#@profile
def print_memory_status(sub_name):
    mem1 = memory_usage_psutil_Percent()
    mem2 = memory_usage_psutil_MB()
    print('---------------------------------------------------------------')
    print(sub_name+' is done.')
    print('Currently using '+str(np.round(mem1,6))+'%, '+str(mem2)+' GB')
    print('---------------------------------------------------------------')
    del mem1, mem2
  
def memory_usage_psutil_MB():
    # return the memory usage in MB
    import psutil
    process = psutil.Process(os.getpid())
    mem = process.memory_info()[0]/1.024e9
    return mem

def memory_usage_psutil_Percent():
    # return the memory usage in percentage like top
    process = psutil.Process(os.getpid())
    mem = process.memory_percent()
    return mem


# Sep 22, 2020: get vertical intergal at midpoints 
#@profile
def VertSum(varin,dp_4d):
    # get midpoint var first
    var1 = (varin[:,:-1,:,:] + varin[:,1:,:,:])/2.
    # vertical integral
    outvar = MV.sum(var1 * dp_4d/100., axis=1)
    # set axis
    outvar.setAxisList(varin[:,0,:,:].getAxisList())

    return outvar 

# defined function is here
#@profile
def get_feedback(exp_cntl,DATA,tas_ano_grd_ann,AXL2d):
    if exp_cntl == 'amip':
        DATA_ann = DATA
        newdata1 = MV.average(DATA_ann,axis=0)/MV.average(tas_ano_grd_ann)
        newdata1.setAxisList(AXL2d)
        newdata2 = cdutil.averager(newdata1,axis='xy',weights='weighted')
        newdata3 = MV.average(DATA_ann,axis=0)
        newdata3.setAxisList(AXL2d)
        newdata4 = cdutil.averager(newdata3,axis='xy',weights='weighted')
        del DATA_ann
    elif exp_cntl == 'piControl':
        cdutil.setTimeBoundsMonthly(DATA)
        DATA_ann = cdutil.YEAR(DATA)
        slope,intercept = genutil.statistics.linearregression(DATA_ann,x = tas_ano_grd_ann)
        newdata1 = slope
        newdata2 = cdutil.averager(slope,axis='xy',weights='weighted')
        newdata3 = intercept
        newdata4 = cdutil.averager(intercept,axis='xy',weights='weighted')
        del slope, intercept, DATA_ann
    return newdata1,newdata2,newdata3,newdata4


#@profile
def main():
    
    cdms.axis.latitude_aliases.append("Y")
    cdms.axis.longitude_aliases.append("X")
    
    ########## MAIN SUBROUTINE STARTS HERE ....
    
#    start = time.time()
    
    # June 6: add output directory for consistent use
    outdir = '../post-data/RadKernel_OtherAMIPs_TestTropMethod/'
    figdir = '../figure/RadKernel_OtherAMIPs_TestTropMethod/'
    
    #outdir = '../post-data/RadKernel_OtherAMIPs_Update_v1/'
    #figdir = '../figure/RadKernel_OtherAMIPs_Update_v1/'
    
    RD.make_dir(outdir)
    RD.make_dir(figdir)
    
    # Oct 6, 2020: add switch for plot
    plotting = False
    
    ### Mar 11: add CMIP5 models
    phases = ['CMIP5','CMIP6']
    #phases = ['CMIP5']
    
    consider_amip = True
    e3sm_cmip = False
    
    # Aug 31, 2020: merge processing E3SM amip simulation here
    e3sm_amip = False
    
    # Sep 20, 2020: consider SUNDOWN
    sundown = True
    
    # Sep 21, 2020: all regrid to standard 17 levels
    RegridVert_method = 'StdVert' #'OldVert', 'StdVert' 
    
    # Sep 25, 2020
    get_netRH_method = 'Mark'
    
    # Sep 24, 2020: add mask kernel
    do_mask_kernel = 'ToZero' # 'OnlyMask'
    
    # Sep 24, 2020: use Mark's kernel data
    kernel_source = 'Mark' # 'Raw','Mark'
    
    
    # Sep 21, 2020: method to get pressure thickness
    dp_use_methods = ['SPC','DynPs'] # # 'DynPs', 'pTrop', 'Fixed','SPC','Mark'  # notion: pTrop is exactly same as Fixed.
    
    get_tropo_methods = ['Mark','Soden'] # 'NoLimit'; 'Soden'; 'Mark'
    
    dlogq_method = 'Mixed' # 'Mark','Mixed'
    
    ta_methods = ['taavg','tapi']  # 'taavg','tapi'
    qsat_methods = ['Mark','Yi'] # 'Mark', 'Yi'
    log_method ='Mark'  # 'Mark','Yi','angie'
    
    
    # Oct 5, 2020: sensitivity tests with varying one case at one time
    test_list = [dp_use_methods[0],get_tropo_methods[0],ta_methods[0],qsat_methods[0]]
    
    count = 0
    
    for i,l in  enumerate([dp_use_methods, get_tropo_methods, ta_methods, qsat_methods]):
        for j in range(0,len(l)):
            test_list1 = copy.deepcopy(test_list)
            test_list1[i] = l[j]
    
            if test_list1 == test_list:
                count += 1
                if count > 1:
                    continue
    
            print('test_list1 = ', test_list1)
    
            dp_use_method = test_list1[0]
            get_tropo_method = test_list1[1]
            ta_method = test_list1[2]
            qsat_method = test_list1[3]
    
            dlogq_out = ta_method+'.'+qsat_method+'.'+log_method
     
            # -------------------------------------------------------------------
            #AddString = get_tropo_method+'Trop_'+dp_use_method+'_'+RegridVert_method+'_MaskKernel'+do_mask_kernel+'_KernelFrom'+kernel_source+'_netRHFrom'+get_netRH_method+'_dlogqMethodFrom'+dlogq_out
            
            AddString = get_tropo_method+'Trop_'+dp_use_method+'_'+RegridVert_method+'_MaskKernel'+do_mask_kernel+'_KernelFrom'+kernel_source+'_netRHFrom'+get_netRH_method+'_dlogqMethodFrom'+dlogq_out+'_NoExecTest_update.clrsky.albedo.fdbk_update.to.tas'
    
            print('AddString = ',AddString)
    
            # ---------------------
            #----------------------
            # ---------------------
            if e3sm_cmip:
                consider_amip = False
                phases = ['CMIP6']
            
            if e3sm_amip:
                phases = ['CMIP6']
                direc_data='/work/qin4/Data/E3SM_data/'
            
                fname1 = '20200428.DECKv1b_amip1-CFMIP.ne30_oEC.cori-knl-L'
                fname2 = '20200428.DECKv1b_amip1.plus4K-CFMIP.ne30_oEC.cori-knl-L'
            
                direc_data1 = direc_data+'/'+fname1+'/'
                direc_data2 = direc_data+'/'+fname2+'/'
                
                yrS = 1979
                yrE = 2014
                monS = 1
                monE = 12
                
                yrS_4d='{:04d}'.format(yrS)
                yrE_4d='{:04d}'.format(yrE)
                monS_2d='{:02d}'.format(monS)
                monE_2d='{:02d}'.format(monE) 
                
                nyears_e3sm = yrE - yrS + 1
                
                exp1='amip'
                exp2='amip_p4K'
            
            for phase in phases:
                print('----------Hi, here is processing ',phase,' Data----------------')
            
                direc_kernel='/work/qin4/script/cloud-radiative-kernels/'
                direc_data='/p/user_pub/xclim/'+phase+'/'
            
                if e3sm_amip:
                    amips = ['amip-p4K']
                else:
                    if phase == 'CMIP5':
            #            amips = ['amip4K','amip4xCO2','amipFuture']
                        amips = ['amip4K']
                    else:
            #           amips = ['amip-p4K','amip-4xCO2','amip-future4K','amip-m4K']
                        amips = ['amip-p4K']
            
                for iamip,samip in enumerate(amips):
                    print('--------------- amips is ', samip,'------------------------')
                    ### Jan 23: merge amip and cmip experiments in the same script. here.
                    All_project_cntl = ['CMIP','CMIP']
                    All_exp_cntl = ['piControl','amip']
                    All_project_new = ['CMIP','CFMIP']
                    if phase == 'CMIP5':
                        All_exp_new = ['abrupt4xCO2',samip]
                    else:
                        All_exp_new = ['abrupt-4xCO2',samip]
                
                    if phase == 'CMIP5':
                        All_nyears=[150,5]
                    else:
                        All_nyears=[150,36]
                
                    All_startyear = [1,1]
                    
                    comp='atmos'
                    freq='mon'
                    # Sep 22, 2020: add ps as an input variable to get pressure thickness later 
                    var2d = ["ts","rlut","rsut","rlutcs","rsutcs","rsdt","rsus","rsds","rsuscs","rsdscs","ps","tas"]
            #        var3d = ['ta','hus']
                    var3d = ['hus','ta']
            
                    var = var2d + var3d
                
                    usedin = 'RadKernel'
                
                    if e3sm_amip:
                        used_insts = ['E3SM-project']
                        used_models = ['E3SM-1-0']
                        used_ripf = ['r1i1p1f1']
                        All_nyears = [150,nyears_e3sm]
                        print('used_models',used_models)
                        print('model number is',len(used_models))
                    else:
                        # get model list
                        used_insts,used_models,used_ripf = RD.find_model_list(direc_data,All_project_cntl,All_exp_cntl,All_project_new,All_exp_new,comp,freq,var,consider_amip)
                        print('used_insts',used_insts)
                        print('used_models',used_models)
                        print('used_ripf',used_ripf)
                        print('model number is',len(used_models))
                
#                    for icase in range(1,2): # amip
                    for icase in range(0,1): # coupled
            
                        project_cntl = All_project_cntl[icase]
                        exp_cntl = All_exp_cntl[icase]
                        project_new = All_project_new[icase]
                        exp_new = All_exp_new[icase]
                        nyears = All_nyears[icase]
                        startyear = All_startyear[icase]
                    
                        print(project_cntl, exp_cntl, project_new, exp_new, nyears)

                        AddTag = str(startyear)+'yr-'+str(startyear+nyears-1)+'yr'

                        dir1 = dir()

                        print('--------------------------------------------------')
                        print('preliminary scope variables are ',dir())
                        print('--------------------------------------------------')
            
                        # read model data
                        for imod in range(len(used_models)):
            #            for imod in range(4,len(used_models)): # Start from MIROC5
            #            for imod in range(0,2):
            
                            print(imod,used_models[imod])
            
                            # special process for CanESM2 and HadGEM2-ES
                            if used_models[imod] == 'CanESM2' and 'amip' in All_exp_new[icase]:
                                used_models[imod] = 'CanAM4'
                            elif used_models[imod] == 'HadGEM2-ES' and 'amip' in All_exp_new[icase]:
                                used_models[imod] = 'HadGEM2-A'
            
#                            if used_models[imod] not in ['MPI-ESM-LR','MPI-ESM-MR','MRI-CGCM3']:
                            if used_models[imod] not in ['MPI-ESM-LR']:
#                            if used_models[imod] not in ['MRI-CGCM3']:
                                continue
           
#                            if used_models[imod] not in ['CNRM-CM5','CanESM2', 'CCSM4', 'IPSL-CM5A-LR', 'MIROC5', 'MPI-ESM-LR', 'MRI-CGCM3',\
#                            'BCC-CSM2-MR','CanESM5' 'CNRM-CM6-1' 'IPSL-CM6A-LR','GFDL-CM4']:
#                                continue
            
                            # June 2: only process E3SM-1-0
                            if e3sm_cmip:
                                if used_models[imod] != 'E3SM-1-0' or exp_cntl == 'amip':
                                    continue
                
                            if phase == 'CMIP6' and used_models[imod] == 'CESM2' and nyears == 36:
                                nyears = 35 # CESM2 only has the data to 2013 not 2014. So cut the nyear to 35 here to avoid triggering error.
                   
                            # Oct 3, 2020: define a big dictionary to save data
                            dic_invar = {}

                            if e3sm_amip:
                                dic_invar = read_data_e3sm(var,var2d,var3d,direc_data1,direc_data2,exp1,exp2,yrS_4d,monS_2d,yrE_4d,monE_2d,dic_invar)
                            else:
                                dic_invar = read_data(var,var2d,var3d,nyears,direc_data,project_new,exp_new,comp,freq,used_insts[imod],used_models[imod],used_ripf[imod],project_cntl,exp_cntl,phase,startyear,dic_invar)

                            gc.collect(2)

                            sub_name = 'Reading all data'
                            print_memory_status(sub_name)
            
                            ##================================================================================###
                            ###------------------------- Read kernel data -----------------------###
            
                            if kernel_source == 'Raw':
                                kernel_dir = "/work/qin4/Data/Huang_kernel_data/"
            
                                ### attention: the kernel data's latitude is from north to south [90....-90]
                                ### we need to reverse it.
                                f1 = cdms.open(kernel_dir+"RRTMG_ts_toa_cld_highR.nc")
                                dic_invar['ts_KernCld'] = f1('lwkernel')[:,::-1,:]
                                f1.close()
                                f1 = cdms.open(kernel_dir+"RRTMG_ts_toa_clr_highR.nc")
                                dic_invar['ts_KernClr'] = f1('lwkernel')[:,::-1,:]
                                f1.close()
                    
                                f1 = cdms.open(kernel_dir+"RRTMG_t_toa_cld_highR.nc")
                                dic_invar['t_KernCld'] = f1('lwkernel')[:,:,::-1,:]
                                f1.close()
                                f1 = cdms.open(kernel_dir+"RRTMG_t_toa_clr_highR.nc")
                                dic_invar['t_KernClr'] = f1('lwkernel')[:,:,::-1,:]
                                f1.close()
                    
                                f1 = cdms.open(kernel_dir+"RRTMG_wv_sw_toa_cld_highR.nc")
                                dic_invar['wv_sw_KernCld'] = f1('swkernel')[:,:,::-1,:]
                                f1.close()
                                f1 = cdms.open(kernel_dir+"RRTMG_wv_sw_toa_clr_highR.nc")
                                dic_invar['wv_sw_KernClr'] = f1('swkernel')[:,:,::-1,:]
                                f1.close()
                    
                                f1 = cdms.open(kernel_dir+"RRTMG_wv_lw_toa_cld_highR.nc")
                                dic_invar['wv_lw_KernCld'] = f1('lwkernel')[:,:,::-1,:]
                                f1.close()
                                f1 = cdms.open(kernel_dir+"RRTMG_wv_lw_toa_clr_highR.nc")
                                dic_invar['wv_lw_KernClr'] = f1('lwkernel')[:,:,::-1,:]
                                f1.close()
                    
                                f1 = cdms.open(kernel_dir+"RRTMG_alb_toa_cld_highR.nc")
                                dic_invar['alb_KernCld'] = f1('swkernel')[:,::-1,:]
                                f1.close()
                                f1 = cdms.open(kernel_dir+"RRTMG_alb_toa_clr_highR.nc")
                                dic_invar['alb_KernClr'] = f1('swkernel')[:,::-1,:]
                                f1.close()
            
                            elif kernel_source == 'Mark':
            
                                # Sep 24, 2020: read Mark's kernel data and further processing
                                kernel_dir = '/work/zelinka1/kernels/'
            
                                f1 = cdms.open(kernel_dir+'TOA_kernels.nc')
            
                                dic_invar['t_KernCld'] = f1('big_LW_Tkern')[3,:,:,:,:]
                                dic_invar['t_KernClr'] = f1('big_LWclr_Tkern')[3,:,:,:,:]
                                dic_invar['wv_lw_KernCld'] = f1('big_LW_qkern')[3,:,:,:,:]
                                dic_invar['wv_lw_KernClr'] = f1('big_LWclr_qkern')[3,:,:,:,:]
                                dic_invar['wv_sw_KernCld'] = f1('big_SW_qkern')[3,:,:,:,:]
                                dic_invar['wv_sw_KernClr'] = f1('big_SWclr_qkern')[3,:,:,:,:]
            
                                dic_invar['alb_KernCld'] = f1('big_SW_akern_surf')[3,:,:,:]
                                dic_invar['alb_KernClr'] = f1('big_SWclr_akern_surf')[3,:,:,:]
                                dic_invar['ts_KernCld'] = f1('big_LW_Tkern_surf')[3,:,:,:]
                                dic_invar['ts_KernClr'] = f1('big_LWclr_Tkern_surf')[3,:,:,:]
            
                                print('t_KernCld from Mark shape is',dic_invar['t_KernCld'].shape)
                                print('t_KernCld from Mark minmax is',genutil.minmax(dic_invar['t_KernCld']))
            
                                f1.close()

                                print('f1 getrefcount is',sys.getrefcount(f1))
                                print('f1 getsizeof is',sys.getsizeof(f1))
                                del f1
            
                            sub_name = 'Reading kernel data'
                            print_memory_status(sub_name)
            
                            # =========================================================================================================#
                            # Vertical Regridding 
                            # =========================================================================================================#
                            if RegridVert_method == 'StdVert':
                                stdlevs = np.array([100000, 92500, 85000, 70000, 60000, 50000, 40000, 30000, 25000, 20000, 15000, 10000, 7000, 5000, 3000, 2000, 1000])/100.
                                std_plev = cdms.createAxis(stdlevs)
                                std_plev.designateLevel()
                                std_plev.units = 'hPa'
                                std_plev.id = 'plev'
                                
                                these_vars = ['ta_ano','hus_ano','ta_pi','ta_ab','hus_pi','hus_ab','t_KernCld','t_KernClr','wv_sw_KernCld','wv_sw_KernClr','wv_lw_KernCld','wv_lw_KernClr']
                                sub_name = 'Start Vertical Regrid'
                                print_memory_status(sub_name)
            
                                for ivar in these_vars:
                                    if kernel_source == 'Mark' and 'Kern' in ivar:
                                        dic_invar[ivar+'_vert'] = dic_invar[ivar]
                                    else:
                                        dic_invar[ivar+'_vert'] = Vert_RegridStd(dic_invar[ivar],stdlevs,std_plev)
            
            #                        dic_invar[ivar+'_vert'] = Vert_RegridStd_Loop(dic_invar[ivar],stdlevs,std_plev)
            #                        newdata = DATA.pressureRegrid(std_plev)
            
                                    print('Vertical regrid ',ivar,'from shape',dic_invar[ivar].shape,'to shape',dic_invar[ivar+'_vert'].shape)
                                    del dic_invar[ivar]
            
                            else:
                                if len(dic_invar['t_KernCld'].getLevel()[:]) > len(dic_invar['ta_ano'].getLevel()[:]):
                                    HAX = dic_invar['t_KernCld'].getAxisList()
                                    stdlevs=dic_invar['ta_ano'].getLevel()[:]
                                    if dic_invar['ta_ano'].getLevel().units == 'Pa':
                                        stdlevs=stdlevs/100.
                                    std_plev = cdms.createAxis(stdlevs)
                                    std_plev.designateLevel()
                                    std_plev.units = "hPa"
                                    std_plev.id = "plev"
                                    ## vertically regrid kernel pressure to model pressure, because model pressure grid is less than kernel one
                                    ## model: 19 levels; kernel: 24 levels
                                    these_vars = ['t_KernCld','t_KernClr','wv_sw_KernCld','wv_sw_KernClr','wv_lw_KernCld','wv_lw_KernClr']
                                else:
                                    HAX = dic_invar['ta_pi'].getAxisList()
                                    stdlevs=dic_invar['t_KernCld'].getLevel()[:]
                                    if dic_invar['t_KernCld'].getLevel().units == 'Pa':
                                        stdlevs = stdlevs/100.
                                    std_plev = cdms.createAxis(stdlevs)
                                    std_plev.designateLevel()
                                    std_plev.units = "hPa"
                                    std_plev.id = "plev"
                                    ## vertically regrid model pressure to kernel pressure
                                    these_vars = ['ta_ano','hus_ano','ta_pi','ta_ab','hus_pi','hus_ab']
            
                                for ivar in these_vars:
                                    dic_invar[ivar+'_vert'] = Vert_RegridStd(dic_invar[ivar],stdlevs,std_plev)
                                    print('Vertical regrid ',ivar,'from shape',dic_invar[ivar].shape,'to shape',dic_invar[ivar+'_vert'].shape)
                                    del dic_invar[ivar]
            
                            del stdlevs, std_plev

                            sub_name = 'Vertical regrid '
                            print_memory_status(sub_name)
            
                            # =========================================================================================================#
                            # expand kernel data into several years * 12 month
                            # =========================================================================================================#
                            stdtime=dic_invar['ta_ano_vert'].getTime()[:]
                            kern_time = cdms.createAxis(stdtime)
                            kern_time.designateTime()
                            kern_time.units = dic_invar['ta_ano_vert'].getTime().units
                            kern_time.id = "time"
                            del stdtime
                
                            these_vars = ['ts_KernCld','ts_KernClr','alb_KernCld','alb_KernClr',\
                            't_KernCld_vert','t_KernClr_vert','wv_sw_KernCld_vert','wv_sw_KernClr_vert','wv_lw_KernCld_vert','wv_lw_KernClr_vert']
                            HAX = dic_invar['t_KernCld_vert'].getAxisList()
                
                            for ivar in these_vars:
                                DATA = dic_invar[ivar]
                                if len(DATA.shape) == 3: # (time,lat,lon)
                                    newdata2 = np.tile(DATA,(nyears,1,1))
                                    newdata2 = np.float32(newdata2)
                                    newdata2.setAxis(0,kern_time)
                                    newdata2.setAxis(1,HAX[2])
                                    newdata2.setAxis(2,HAX[3])
                                elif len(DATA.shape) == 4: # (time,lev,lat,lon)
                                    newdata2 = np.tile(DATA,(nyears,1,1,1))
                                    newdata2 = np.float32(newdata2)
                                    newdata2.setAxis(0,kern_time)
                                    newdata2.setAxis(1,HAX[1])
                                    newdata2.setAxis(2,HAX[2])
                                    newdata2.setAxis(3,HAX[3])
                                dic_invar[ivar+'_mon'] = newdata2
                                print(ivar,'expanding from shape',DATA.shape,'to shape',newdata2.shape)
                                del DATA,newdata2,dic_invar[ivar]

                            del kern_time,HAX
            
                            sub_name = 'Expanding kernel data to 12*nyears data '
                            print_memory_status(sub_name)
            
                            # =========================================================================================================#
                            # horizontally regrid the input data to kernel grid
                            # =========================================================================================================#
                            these_vars = ['tas_ano','alb_ano',\
                            'ta_ano_vert','hus_ano_vert','ta_pi_vert','ta_ab_vert','hus_pi_vert','hus_ab_vert',\
                            'rsdt_ab','rsdt_ano','rsutcs_ano','rlutcs_ano','SWCRE_ano','LWCRE_ano','netCRE_ano','rlut_ano','rsut_ano','ps_pi']
                
                            kernel_grid = dic_invar['ts_KernCld_mon'].getGrid() # get new grid info
                
                            for ivar in these_vars:
                                DATA = dic_invar[ivar]
                                newdata = DATA.regrid(kernel_grid,regridTool='esmf',regridMethod='linear')
                                print(ivar,'horizonally regrid from shape',DATA.shape,'to shape', newdata.shape) #(12*nyears, 64, 128)
                                dic_invar[ivar+'_grd'] = newdata
                                del newdata, DATA, dic_invar[ivar]
            
                            del kernel_grid
                            sub_name = 'Horizontally regrid the input data to kernel grid '
                            print_memory_status(sub_name)

                            #################################################################################
                            #### Temperature feedback calculation
                            #################################################################################
                
                            lat  = dic_invar['t_KernCld_vert_mon'].getLatitude()[:]
                            lon  = dic_invar['t_KernCld_vert_mon'].getLongitude()[:]
                            time = dic_invar['t_KernCld_vert_mon'].getTime()[:]
                            lev  = dic_invar['t_KernCld_vert_mon'].getLevel()[:]
            
                            nlat = len(lat)
                            nlon = len(lon)
                            nlev = len(lev)
                            ntime = len(time)
                
                            # March 18: change the raw data to get axis from t_KernCld_vert_mon to ta_pi_vert_grd
                            # however, I don't know why the axis getting from t_KernCld_vert_mon cannot be used for two specific models: CCSM4 and CM4 in processing coupled simulation. So strange.
                            AXL4d = dic_invar['ta_pi_vert_grd'].getAxisList()
                            AXL3d = dic_invar['ta_pi_vert_grd'][:,0,:,:].getAxisList()
                            AXL2d = dic_invar['ta_pi_vert_grd'][0,0,:,:].getAxisList()
              
                            # =========================================================================================================#
                            # get tropopause height 
                            # =========================================================================================================#
                            # April 27, 2020: update the calculation of tropopause pressure by using the time-varying one. 
                            # refer to Mark et al. (2020)
                
                            if get_tropo_method == 'NoLimit':
                                lev_4d = np.transpose(np.tile(lev,(ntime,nlat,nlon,1)),(0,3,1,2))
                                lev_4d = np.float32(lev_4d)
                                # temporially set the tropopause height is the minimum level
                                dic_invar['p_tropopause'] = cdms.asVariable(np.transpose(np.tile(lev_4d[:,-1,:,:],(nlev,1,1,1)),(1,0,2,3)))
                                dic_invar['p_tropopause'].setAxisList(AXL4d)
            
                            elif get_tropo_method == 'Soden':
                                lev_4d = np.transpose(np.tile(lev,(ntime,nlat,nlon,1)),(0,3,1,2))
                                lev_4d = np.float32(lev_4d)
                                # crude tropopause estimate: 100 hPa in the tropics, lowering with 
                                # cosine to 300 hPa at the poles.
                                pi = 3.1415926
                                x = np.cos(lat*pi/100.)
                                p_tropopause_zonalmean = 300-200*x
                                p_tropopause = np.float32(np.transpose(np.tile(p_tropopause_zonalmean,(ntime,nlev,nlon,1)),(0,1,3,2)))
                                dic_invar['p_tropopause'] = cdms.asVariable(p_tropopause)
                                dic_invar['p_tropopause'].setAxisList(AXL4d)
                                del x, p_tropopause_zonalmean, p_tropopause
                
                            elif get_tropo_method == 'Mark':
                                lev_4d = np.transpose(np.tile(lev,(ntime,nlat,nlon,1)),(0,3,1,2))
                                lev_4d = np.float32(lev_4d)
                                # I don't know why Mark uses the p_tropopause derived from ta_ab_vert_grd. I will try it here.
                                dic_invar['p_tropopause_3d'] = RD.get_tropopause_pressure(dic_invar['ta_ab_vert_grd'])
                                p_tropopause = np.float32(np.transpose(np.tile(dic_invar['p_tropopause_3d'],(nlev,1,1,1)),(1,0,2,3)))
                                dic_invar['p_tropopause'] = cdms.asVariable(p_tropopause)
                                print('p_tropopause shape is',dic_invar['p_tropopause'].shape,'minmax is',genutil.minmax(dic_invar['p_tropopause']))
                                del p_tropopause
                
                            # --------------------------------------------------------------------------------------------
                            # set the temperature change to zero in the stratosphere (mask out stratosphere)
                            if get_tropo_method == 'Soden' or get_tropo_method == 'Mark':
                                dic_invar['ta_ano_vert_grd'] = MV.masked_where(lev_4d<=dic_invar['p_tropopause'],dic_invar['ta_ano_vert_grd']) 
            
                            sub_name = 'get tropopause height '
                            print_memory_status(sub_name)
            
                            # =========================================================================================================#
                            # mask kernel method
                            # =========================================================================================================#
                            if do_mask_kernel == 'ToZero':
                                these_vars = ['ts_KernCld_mon','ts_KernClr_mon','alb_KernCld_mon','alb_KernClr_mon',\
                                't_KernCld_vert_mon','t_KernClr_vert_mon','wv_lw_KernCld_vert_mon','wv_lw_KernClr_vert_mon',\
                                'wv_sw_KernCld_vert_mon','wv_sw_KernClr_vert_mon']
            
                                for svar in these_vars:
                                    dic_invar[svar][dic_invar[svar].mask] = 0.0
            
                            elif do_mask_kernel == 'OnlyMask':
                                if get_tropo_method == 'Soden' or get_Tropo_method == 'Mark':
                                    these_vars = ['t_KernCld_vert_mon','t_KernClr_vert_mon','wv_lw_KernCld_vert_mon','wv_lw_KernClr_vert_mon',\
                                    'wv_sw_KernCld_vert_mon','wv_sw_KernClr_vert_mon']
                                    for ivar in these_vars:
                                        dic_invar[ivar] = MV.masked_where(lev_4d<=dic_invar['p_tropopause'],dic_invar[ivar])
                    
                                # ------------------------------------------------------------
                                # July 7th
                                # important! wv_lw_KernCld_vert_mon and wv_sw_KernCld_vert_mon don't exactly have the same effective value in each level.
                                # Use masked_where to keep them the same. wv_sw_KernCld_vert_mon has missing value at the lowest level, but wv_lw_KernCld_vert_mon not. 
                                # I don't know how to avoid it using a more effective way. Just keep the "stupid" method here. -- July 7
                                # ------------------------------------------------------------
                                dic_invar['wv_lw_KernCld_vert_mon'] = MV.masked_where(dic_invar['wv_sw_KernCld_vert_mon'].mask==True,dic_invar['wv_lw_KernCld_vert_mon'])
                                dic_invar['wv_lw_KernClr_vert_mon'] = MV.masked_where(dic_invar['wv_sw_KernClr_vert_mon'].mask==True,dic_invar['wv_lw_KernClr_vert_mon'])
                                dic_invar['t_KernCld_vert_mon'] = MV.masked_where(dic_invar['wv_sw_KernCld_vert_mon'].mask==True,dic_invar['t_KernCld_vert_mon'])
                                dic_invar['t_KernClr_vert_mon'] = MV.masked_where(dic_invar['wv_sw_KernClr_vert_mon'].mask==True,dic_invar['t_KernClr_vert_mon'])
                                # ------------------------------------------------------------
            
                            sub_name = 'do mask kernel '
                            print_memory_status(sub_name)
            
                            # =========================================================================================================#
                            # get pressure thickness
                            # =========================================================================================================#
                            psurf = 1000.
                            ptop = min(lev)
            
                            if dp_use_method == 'DynPs':
                                # Sep 21, 2020: use dynamic surface pressure (ps_pi_grd) rather than the fixed psurf
                                dic_invar['dp_4d'] = RD.get_dp_DynPs(dic_invar['ps_pi_grd'],ptop,lev,dic_invar['p_tropopause'])
                            elif dp_use_method == 'pTrop':
                                # Sep 20, 2020: update to include p_tropopause
                                dic_invar['dp_4d'] = RD.get_dp(psurf,ptop,lev,dic_invar['p_tropopause'])
                            elif dp_use_method == 'Fixed':
                                # the oldest method
                                dp = RD.get_dp_old(psurf, ptop, lev)
                                dic_invar['dp_4d'] = np.float32(np.transpose(np.tile(np.tile(dp,[1,1,1,1]),[ntime,nlat,nlon,1]),(0,3,1,2)))
                            elif dp_use_method == 'SPC': # Sep 24, 2020
                                trop_wts, atm_wts = CU.get_weights_SPC(dic_invar['ps_pi_grd'], dic_invar['p_tropopause_3d'], dic_invar['ta_pi_vert_grd'])
                                dic_invar['dp_4d'] = trop_wts * 100. # convert to thickness to be consistent with my codes. 
                                del trop_wts, atm_wts,dic_invar['p_tropopause_3d']
                            elif dp_use_method == 'Mark': # Sep 24, 2020
                                trop_wts, atm_wts = CU.get_pressure_weights2new(lev,dic_invar['ps_pi_grd'], dic_invar['p_tropopause_3d'])
                                dic_invar['dp_4d'] = trop_wts * 100. # convert to thickness to be consistent with my codes. 
                                del trop_wts, atm_wts,dic_invar['p_tropopause_3d']
            
                            dic_invar['dp_4d'] = np.float32(dic_invar['dp_4d'])
                            del dic_invar['ps_pi_grd'],psurf,ptop
            
                            sub_name = 'Getting pressure thickness '
                            print_memory_status(sub_name)
            
                            # =========================================================================================================#
                            # set SUNDOWN case
                            # =========================================================================================================#
                            if sundown:
                                SUNDOWN = dic_invar['rsdt_ab_grd']
            
                            # =========================================================================================================#
                            # TS feedback
                            # =========================================================================================================#
                            dic_invar['tas_ano_grd'][dic_invar['tas_ano_grd'].mask] = 0
                            dic_invar['dLW_ts'] = dic_invar['ts_KernCld_mon'] * dic_invar['tas_ano_grd']
                            dic_invar['dLW_ts'].setAxisList(AXL3d)
                
                            # for clear-sky
                            dic_invar['dLW_ts_clr'] = dic_invar['ts_KernClr_mon'] * dic_invar['tas_ano_grd']
                            dic_invar['dLW_ts_clr'].setAxisList(AXL3d)
                            del dic_invar['ts_KernClr_mon']
            
                            sub_name = 'Getting TS feedback '
                            print_memory_status(sub_name)
            
                            # =========================================================================================================#
                            # get total water vapor kernel
                            # =========================================================================================================#
                            dic_invar['wv_KernCld_vert_mon'] = dic_invar['wv_sw_KernCld_vert_mon'] + dic_invar['wv_lw_KernCld_vert_mon']
                            dic_invar['wv_KernClr_vert_mon'] = dic_invar['wv_sw_KernClr_vert_mon'] + dic_invar['wv_lw_KernClr_vert_mon']
                            dic_invar['wv_KernCld_vert_mon'].setAxisList(AXL4d)
                            dic_invar['wv_KernClr_vert_mon'].setAxisList(AXL4d)
            
                            dic_invar['Tq_KernCld_vert_mon'] = dic_invar['t_KernCld_vert_mon'] + dic_invar['wv_KernCld_vert_mon']
                            dic_invar['Tq_KernClr_vert_mon'] = dic_invar['t_KernClr_vert_mon'] + dic_invar['wv_KernClr_vert_mon']
                            dic_invar['Tq_KernCld_vert_mon'].setAxisList(AXL4d)
                            dic_invar['Tq_KernClr_vert_mon'].setAxisList(AXL4d)
            
                            # =========================================================================================================#
                            # TA feedback
                            # =========================================================================================================#
                            # get air temperature feedback
                            dic_invar['ta_ano_vert_grd'][dic_invar['ta_ano_vert_grd'].mask]=0
    
                            invar1 = ['t_KernCld_vert_mon','t_KernClr_vert_mon','wv_KernCld_vert_mon','wv_KernClr_vert_mon']
                            invar2 = ['ta_ano_vert_grd','ta_ano_vert_grd','ta_ano_vert_grd','ta_ano_vert_grd']
                            outvar = ['dLW_ta','dLW_ta_clr','dLW_ta_fxRH','dLW_ta_clr_fxRH']
            
                            dic_invar = get_fdbk(invar1,invar2,outvar,dic_invar,AXL4d, AXL3d, dp_use_method='SPC')
                            print(dic_invar.keys())
              
                            dic_invar['dLW_t_clr_psum'] = dic_invar['dLW_ta_clr_psum'] + dic_invar['dLW_ts_clr']
                            dic_invar['dLW_t_psum'] = dic_invar['dLW_ta_psum'] + dic_invar['dLW_ts']
            
                            dic_invar['dLW_t_clr_psum'].setAxisList(AXL3d)
                            dic_invar['dLW_t_psum'].setAxisList(AXL3d)
            
                            # =========================================================================================================#
                            ### Plotting temperature test 
                            # =========================================================================================================#
                            if plotting:
                                fig = plt.figure(figsize=(18,12))
                                fh = 20
                                plt.suptitle('clear-sky linear test',fontsize=fh)
                                
                                bounds = np.arange(-12,14,2)
                                cmap = pl.cm.RdBu_r
                                bounds2 = np.append(np.append(-500,bounds),500)
                                norm = mpl.colors.BoundaryNorm(bounds2,cmap.N)
                    
                                names = ['ta_ano_vert_grd','ta_ano_vert_grd','t_KernCld_vert_mon','t_KernClr_vert_mon','dLW_ta','dLW_ta_clr']
                    
                                for n,name in enumerate(names):
                                    DATA = MV.average(MV.average(dic_invar[name],axis=0),axis=2)
                                    # zonal-mean structure
                                    ax1 = fig.add_subplot(3,2,n+1)
                                    im1 = ax1.contourf(lat,lev,DATA)
                                    pl.title(name)
                                    cb = plt.colorbar(im1,orientation='vertical',drawedges=True)
                                    plt.gca().invert_yaxis()
                    
                                plt.tight_layout()
                                plt.savefig(figdir+'lat-lon-clear-sky-linear-test-temp-'+phase+'-'+exp_new+'-'+used_models[imod]+'-'+AddString+'_'+AddTag+'.png',bbox_inches='tight')
                                fig.clf()
                                plt.close()
                                gc.collect()
            
                            # =========================================================================================================#
                            # delete unnecessary vars 
                            # =========================================================================================================#
                            delterm = ['dLW_ta','dLW_ta_clr','dLW_ta_fxRH','dLW_ta_clr_fxRH']
                            dic_invar = delete_vars(delterm,dic_invar)
            
                            sub_name = 'Getting TA feedback '
                            print_memory_status(sub_name)
            
                            #################################################################################
                            #### Planck feedback calculation
                            #################################################################################
                            # extend to 4-d
                            dic_invar['tas_ano_grd_4d']   = cdms.asVariable(np.float32(np.transpose(np.tile(dic_invar['tas_ano_grd'],(nlev,1,1,1)),(1,0,2,3))))
                
                            # mask stratosphere
                            if get_tropo_method == 'Soden' or get_tropo_method == 'Mark':
                                dic_invar['tas_ano_grd_4d']   = MV.masked_where(lev_4d <=dic_invar['p_tropopause'],dic_invar['tas_ano_grd_4d'])
                
                            # mask ts_ano_grd_4d where ta_ano_vert_grd is True
                            dic_invar['tas_ano_grd_4d']   = MV.masked_where(dic_invar['ta_ano_vert_grd'].mask==True,dic_invar['tas_ano_grd_4d'])
                            dic_invar['tas_ano_grd_4d'][dic_invar['tas_ano_grd_4d'].mask]=0
            
                            # note: fxRH Planck feedback: T_kernel * uniform warming anomaly + Q_kernel * uniform warming anomaly
                            invar1 = ['t_KernCld_vert_mon','t_KernClr_vert_mon','Tq_KernCld_vert_mon','Tq_KernClr_vert_mon']
                            invar2 = ['tas_ano_grd_4d','tas_ano_grd_4d','tas_ano_grd_4d','tas_ano_grd_4d']
                            outvar = ['dLW_planck','dLW_planck_clr','dLW_planck_fxRH','dLW_planck_clr_fxRH']
            
                            dic_invar = get_fdbk(invar1,invar2,outvar,dic_invar,AXL4d, AXL3d, dp_use_method='SPC')
            
                            # add ts feedback
                            dic_invar['dLW_planck_psum'] = dic_invar['dLW_planck_psum'] + dic_invar['dLW_ts']
                            dic_invar['dLW_planck_clr_psum'] = dic_invar['dLW_planck_clr_psum'] + dic_invar['dLW_ts_clr']
                            dic_invar['dLW_planck_fxRH_psum'] = dic_invar['dLW_planck_fxRH_psum'] + dic_invar['dLW_ts']
                            dic_invar['dLW_planck_clr_fxRH_psum'] = dic_invar['dLW_planck_clr_fxRH_psum'] + dic_invar['dLW_ts_clr']
            
                            dic_invar['dLW_planck_psum'].setAxisList(AXL3d)
                            dic_invar['dLW_planck_clr_psum'].setAxisList(AXL3d)
                            dic_invar['dLW_planck_fxRH_psum'].setAxisList(AXL3d)
                            dic_invar['dLW_planck_clr_fxRH_psum'].setAxisList(AXL3d)
            
                            delterm = outvar
                            dic_invar = delete_vars(delterm,dic_invar)
            
                            sub_name = 'Getting Planck feedback '
                            print_memory_status(sub_name)
            
                            #################################################################################
                            #### Lapse rate feedback calculation
                            #################################################################################
                
                            # difference between ta and ts
                            dic_invar['dt_ano'] = dic_invar['ta_ano_vert_grd'] - dic_invar['tas_ano_grd_4d']
                            dic_invar['dt_ano'] = MV.masked_where(dic_invar['ta_ano_vert_grd'].mask==True,dic_invar['dt_ano'])
                            dic_invar['dt_ano'].setAxisList(AXL4d)
                
                            # mask stratosphere
                            if get_tropo_method == 'Soden' or get_tropo_method == 'Mark':
                                dic_invar['dt_ano']   = MV.masked_where(lev_4d <=dic_invar['p_tropopause'],dic_invar['dt_ano'])
                        
                            dic_invar['dt_ano'][dic_invar['dt_ano'].mask]=0
            
                            # note: fxRH Lapse rate feedback: T_kernel * LR warming anomaly + Q_kernel * LR warming anomaly
                            invar1 = ['t_KernCld_vert_mon','t_KernClr_vert_mon','Tq_KernCld_vert_mon','Tq_KernClr_vert_mon']
                            invar2 = ['dt_ano','dt_ano','dt_ano','dt_ano']
                            outvar = ['dLW_lapserate','dLW_lapserate_clr','dLW_lapserate_fxRH','dLW_lapserate_clr_fxRH']
            
                            dic_invar = get_fdbk(invar1,invar2,outvar,dic_invar,AXL4d, AXL3d, dp_use_method='SPC')
            
                            # ---------------------
                            delterm = outvar+['dt_ano','tas_ano_grd_4d']
                            dic_invar = delete_vars(delterm,dic_invar)
                
                            sub_name = 'Getting LapseRate feedback '
                            print_memory_status(sub_name)
            
                            #################################################################################
                            #### Albedo feedback calculation
                            #################################################################################
                
                            dic_invar['alb_ano_grd'][dic_invar['alb_ano_grd'].mask]=0
            
                            invar1 = ['alb_KernCld_mon','alb_KernClr_mon']
                            invar2 = ['alb_ano_grd','alb_ano_grd']
                            outvar = ['dSW_alb','dSW_alb_clr']
            
                            for ivar,svar in enumerate(invar1):
                                ovar = outvar[ivar]
                                dic_invar[ovar] = dic_invar[invar1[ivar]] * dic_invar[invar2[ivar]]
            
                            sub_name = 'Getting Albedo feedback '
                            print_memory_status(sub_name)
            
                            #################################################################################
                            #### Water vapor feedback calculation
                            #################################################################################
            
                            print('-------------------We use the dlogq_method from ', dlogq_method, '------------------------------')
    
                            # ===============================================================#
                            # dlogq method
                            # ===============================================================#
                            if dlogq_method == 'Mixed':
                                # ---------------------------------------------------------------#
                                #  ta_method
                                # ---------------------------------------------------------------#
                                if ta_method == 'taavg':
                                    avgta = (dic_invar['ta_pi_vert_grd'] + dic_invar['ta_ab_vert_grd'])/2.0
                                elif ta_method == 'tapi':
                                    avgta = dic_invar['ta_pi_vert_grd']
            
                                # ---------------------------------------------------------------#
                                #  qsat_method 
                                # ---------------------------------------------------------------#
                                if qsat_method == 'Mark':
                                    qs0,qs1 = qsat_blend_Mark(avgta,lev_4d)
                                    rh0=dic_invar['hus_pi_vert_grd']/qs0
                                    q1k=rh0*qs1
                                elif qsat_method == 'Yi':
                                    qs0 = RD.r_star_GG(lev_4d*100.,avgta)
                                    qs0 = np.float32(qs0)
                                    ta1k = avgta+1.0
                                    qs1 = RD.r_star_GG(lev_4d*100.,ta1k)
                                    qs1 = np.float32(qs1)
                                    rh0 = dic_invar['hus_pi_vert_grd']/qs0
                                    q1k = rh0*qs1
                                del avgta
            
                                # ---------------------------------------------------------------#
                                #  log_method
                                # ---------------------------------------------------------------#
                                if log_method == 'Mark':
                                    dic_invar['dlogq1k'] = MV.log(qs1) - MV.log(qs0)
                                    dlogq = MV.log(dic_invar['hus_ab_vert_grd']) - MV.log(dic_invar['hus_pi_vert_grd'])
                                elif log_method == 'Yi':
                                    dic_invar['dlogq1k'] = MV.log(q1k) - MV.log(dic_invar['hus_pi_vert_grd'])
                                    dlogq = MV.log(dic_invar['hus_ab_vert_grd']) - MV.log(dic_invar['hus_pi_vert_grd'])
                                elif log_method == 'angie':
                                    dic_invar['dlogq1k'] = (q1k-dic_invar['hus_pi_vert_grd'])/dic_invar['hus_pi_vert_grd']
                                    dlogq = dic_invar['hus_ano_vert_grd']/dic_invar['hus_pi_vert_grd']
            
                                # ---------------------------------------------------------------#
                                # get dlogq2 here
                                # ---------------------------------------------------------------#
                                dic_invar['dlogq2'] = np.ma.true_divide(dlogq,dic_invar['dlogq1k'])
                                dic_invar['dlogq2'].setAxisList(dic_invar['hus_pi_vert_grd'].getAxisList())    
            
                                print('dlogq2 shape is',dic_invar['dlogq2'].shape,', minmax is',genutil.minmax(dic_invar['dlogq2']))
              
                                del qs0,qs1,rh0,q1k,dlogq
                                delterm = ['dlogq1k','hus_pi_vert_grd','hus_ab_vert_grd','ta_pi_vert_grd']
                                dic_invar = delete_vars(delterm,dic_invar)
            
                            elif dlogq_method == 'Mark':
                                q_norm_flag = ''
                                blend_flag = ''
                                lognorm_flag = ''
                                # note: Mark's function needs lev's unit Pa. I should change my lev's unit hPa to Pa. so multiply 100.
                                norm_hus_anom = CU.ta_normalized_qv(lev*100,dic_invar['ta_pi_vert_grd'],dic_invar['ta_ab_vert_grd'],dic_invar['hus_pi_vert_grd'],dic_invar['hus_ab_vert_grd'],q_norm_flag,blend_flag,lognorm_flag)
                                dic_invar['dlogq2'] = norm_hus_anom
                                dic_invar['dlogq2'].setAxisList(AXL4d)
                                print('dlogq2 shape is',dic_invar['dlogq2'].shape,', minmax is',genutil.minmax(dic_invar['dlogq2']))
            
                                del norm_hus_anom
                                delterm = ['hus_pi_vert_grd','hus_ab_vert_grd','ta_pi_vert_grd']
                                dic_invar = delete_vars(delterm,dic_invar)
            
                            # ----------------------------------------------------------------#
                            # mask the statrosphere
                            # ----------------------------------------------------------------#
                            if get_tropo_method  == 'Soden' or get_tropo_method == 'Mark':
                                dic_invar['dlogq2'] = MV.masked_where(lev_4d<=dic_invar['p_tropopause'],dic_invar['dlogq2']) 
                            del lev_4d
            
                            dic_invar['dlogq2'][dic_invar['dlogq2'].mask] = 0
    
    
                            # ----------------------------------------------------------------#
                            # get water vapor feedback here
                            # ----------------------------------------------------------------#
                            invar1 = ['wv_lw_KernCld_vert_mon','wv_sw_KernCld_vert_mon','wv_lw_KernClr_vert_mon','wv_sw_KernClr_vert_mon',\
                            'wv_lw_KernCld_vert_mon','wv_sw_KernCld_vert_mon','wv_lw_KernClr_vert_mon','wv_sw_KernClr_vert_mon']
                            invar2 = ['dlogq2','dlogq2','dlogq2','dlogq2',\
                            'ta_ano_vert_grd','ta_ano_vert_grd','ta_ano_vert_grd','ta_ano_vert_grd']
            
                            outvar = ['dLW_q','dSW_q','dLW_q_clr','dSW_q_clr',\
                            'dLW_q_fxRH','dSW_q_fxRH','dLW_q_clr_fxRH','dSW_q_clr_fxRH']
            
                            # note: fxRH LW/SW WV feedback: Q_kernel * ta warming anomaly
                            dic_invar = get_fdbk(invar1,invar2,outvar,dic_invar,AXL4d, AXL3d, dp_use_method='SPC')
    
                            # ---------------------------------------------------------------#
                            # get final RH feedback
                            # ---------------------------------------------------------------#
                            # July 7: final RH feedback related to RH change
                            # RH feedback = default water vapor feedback - (water vapor kernel * uniform warming anomaly (Ts) - water vapor kernel * lapse rate temperature anomaly)
                            # RH feedback = default water vapor feedback - water vapor kernel * atmospheric temperature change
            
                            if get_netRH_method == 'Raw':
                                invar1 = ['dLW_q','dSW_q','dLW_q_clr','dSW_q_clr']
                                outvar = ['dLW_netRH','dSW_netRH','dLW_clr_netRH','dSW_clr_netRH']
            
                                for ivar,svar in enumerate(invar1):
                                    ovar = outvar[ivar]
                                    dic_invar[ovar] = dic_invar[svar] - dic_invar[svar+'_fxRH']
                                    dic_invar[ovar].setAxisList(AXL4d)
            
                                    if dp_use_method == 'DynPs' or dp_use_method =='SPC':
                                        dic_invar[ovar+'_psum'] = VertSum(dic_invar[ovar],dic_invar['dp_4d'])
                                    else:
                                        dic_invar[ovar+'_psum'] = MV.sum(dic_invar[ovar]*dic_invar['dp_4d']/100., axis=1)
            
                                    dic_invar[ovar+'_psum'].setAxisList(AXL3d)
            
                            elif get_netRH_method == 'Mark':
            
                                invar1 = ['dLW_q','dSW_q','dLW_q_clr','dSW_q_clr']
                                outvar = ['dLW_netRH','dSW_netRH','dLW_clr_netRH','dSW_clr_netRH']
            
                                for ivar,svar in enumerate(invar1):
                                    ovar = outvar[ivar]
                                    dic_invar[ovar+'_psum'] = dic_invar[svar+'_psum'] - dic_invar[svar+'_fxRH_psum']
                                    dic_invar[ovar+'_psum'].setAxisList(AXL3d)
            
                            # ---------------------------------------------------------------#
                            # set effect of SUNDOWN
                            # ---------------------------------------------------------------#
                            if sundown:
                                for svar in ['dSW_alb','dSW_alb_clr','dSW_q_psum','dSW_q_clr_psum',\
                                'dSW_q_fxRH_psum','dSW_q_clr_fxRH_psum','dSW_netRH_psum','dSW_clr_netRH_psum']:
                                    dic_invar[svar][SUNDOWN==0] = 0.
                            del SUNDOWN
             
                            # ---------------------------------------------------------------#
                            ### Plotting Water vapor test
                            # ---------------------------------------------------------------#
                            if plotting:
                                fig = plt.figure(figsize=(18,12))
                                fh = 20
                                plt.suptitle('clear-sky linear test',fontsize=fh)
        
                                names = ['dlogq2','wv_lw_KernCld_vert_mon','wv_sw_KernCld_vert_mon',\
                                'wv_lw_KernClr_vert_mon','wv_sw_KernClr_vert_mon','dLW_q','dSW_q','dLW_q_clr','dSW_q_clr']
                                for n,name in enumerate(names):
                                    DATA = MV.average(MV.average(dic_invar[name],axis=0),axis=2)
                                    # zonal-mean structure
                                    ax1 = fig.add_subplot(4,3,n+1)
                                    im1 = ax1.contourf(lat,lev,DATA)
                                    pl.title(name,fontsize=fh)
                                    cb = plt.colorbar(im1,orientation='vertical',drawedges=True)
                                    plt.gca().invert_yaxis()
                    
                                plt.tight_layout()
                                plt.savefig(figdir+'lat-lon-clear-sky-linear-test-wv-'+phase+'-'+exp_new+'-'+used_models[imod]+'-'+AddString+'_'+AddTag+'.png',bbox_inches='tight')
                                fig.clf()
                                plt.close()
                                gc.collect()
            
                            # ---------------------------------------------------------------#
                            # delete unnecessary vars 
                            # ---------------------------------------------------------------#
                            delterm = ['dlogq2','dLW_q','dSW_q','dLW_q_clr','dSW_q_clr',\
                            'dLW_q_fxRH','dLW_q_clr_fxRH','dSW_q_fxRH','dSW_q_clr_fxRH',\
                            't_KernCld_vert_mon','t_KernClr_vert_mon','ts_KernCld_mon','ts_KernClr_mon',\
                            'wv_lw_KernCld_vert_mon','wv_sw_KernCld_vert_mon','wv_lw_KernClr_vert_mon','wv_sw_KernClr_vert_mon']
                            dic_invar = delete_vars(delterm,dic_invar)
            
                            if get_netRH_method == 'Raw':
                                delterm = ['dLW_netRH','dSW_netRH','dSW_clr_netRH','dSW_clr_netRH']
                                dic_invar = delete_vars(delterm,dic_invar)
            
                            sub_name = 'Getting water vapor feedback '
                            print_memory_status(sub_name)
             
                            #################################################################################
                            #### Adjusted cloud feedback calculation
                            #################################################################################
                
                            # ---------------------------------------------------------------#
                            # calculate cloud masking 
                            # ---------------------------------------------------------------#
    
                            invar1 = ['dLW_t_clr_psum','dLW_q_clr_psum','dSW_q_clr_psum','dSW_alb_clr']
                            invar2 = ['dLW_t_psum','dLW_q_psum','dSW_q_psum','dSW_alb']
            
                            for ivar,svar in enumerate(invar1):
                                dic_invar[invar2[ivar]+'_mask'] = dic_invar[invar1[ivar]] - dic_invar[invar2[ivar]]
                                dic_invar[invar2[ivar]+'_mask'].setAxisList(AXL3d)
               
                            # ---------------------------------------------------------------#
                            # adjusted CRE
                            # ---------------------------------------------------------------#
    
                            invar1 = ['dLW_t_psum_mask','dSW_q_psum_mask','SWCRE_ano_grd','LWCRE_ano_grd']
                            invar2 = ['dLW_q_psum_mask','dSW_alb_mask','dSW_adj','dLW_adj']
                            outvar = ['dLW_adj','dSW_adj','SWCRE_ano_grd_adj','LWCRE_ano_grd_adj']
                            for ivar,svar in enumerate(invar1):
                                dic_invar[outvar[ivar]] = dic_invar[invar1[ivar]] + dic_invar[invar2[ivar]]
                                dic_invar[outvar[ivar]].setAxisList(AXL3d)
            
                            dic_invar['net_adj'] = dic_invar['dLW_adj'] + dic_invar['dSW_adj']
                            dic_invar['net_adj'].setAxisList(AXL3d)
            
                            dic_invar['netCRE_ano_grd_adj'] = dic_invar['netCRE_ano_grd'] + dic_invar['net_adj']
                            dic_invar['netCRE_ano_grd_adj'].setAxisList(AXL3d)
            
                            sub_name = 'Getting adjusted cloud feedback '
                            print_memory_status(sub_name)
            
                            # ---------------------------------------------------------------#
                            # get cloudy residual term
                            # ---------------------------------------------------------------#
    
                            # get sum of kernel effect 
                            dic_invar['dLW_cld_sum'] = dic_invar['dLW_t_psum'] + dic_invar['dLW_q_psum'] + dic_invar['LWCRE_ano_grd_adj']
                            dic_invar['dSW_cld_sum'] = dic_invar['dSW_alb'] + dic_invar['dSW_q_psum'] + dic_invar['SWCRE_ano_grd_adj']
                            dic_invar['dLW_cld_sum'].setAxisList(AXL3d)
                            dic_invar['dSW_cld_sum'].setAxisList(AXL3d)
    
                            dic_invar['dnet_cld_sum'] = dic_invar['dLW_cld_sum'] + dic_invar['dSW_cld_sum']
                            dic_invar['dnet_cld_sum'].setAxisList(AXL3d)
            
                            dic_invar['dLW_clr_sum'] = dic_invar['dLW_t_clr_psum'] + dic_invar['dLW_q_clr_psum']
                            dic_invar['dSW_clr_sum'] = dic_invar['dSW_alb_clr'] + dic_invar['dSW_q_clr_psum']
                            dic_invar['dLW_clr_sum'].setAxisList(AXL3d)
                            dic_invar['dSW_clr_sum'].setAxisList(AXL3d)
    
                            dic_invar['dnet_clr_sum'] = dic_invar['dLW_clr_sum'] + dic_invar['dSW_clr_sum']
                            dic_invar['dnet_clr_sum'].setAxisList(AXL3d)
                
                            # get the TOA anomalies from direct model data
                            dic_invar['dLW_cld_dir'] = -1. * dic_invar['rlut_ano_grd']
            #                dic_invar['dSW_cld_dir'] = dic_invar['rsdt_ano_grd'] - dic_invar['rsut_ano_grd']
                            dic_invar['dSW_cld_dir'] = -1. * dic_invar['rsut_ano_grd']
                            dic_invar['dLW_cld_dir'].setAxisList(AXL3d)
                            dic_invar['dSW_cld_dir'].setAxisList(AXL3d)
    
                            dic_invar['dnet_cld_dir'] = dic_invar['dLW_cld_dir'] + dic_invar['dSW_cld_dir']
                            dic_invar['dnet_cld_dir'].setAxisList(AXL3d)
            
                            dic_invar['dLW_clr_dir'] = -1. * dic_invar['rlutcs_ano_grd']
                            dic_invar['dSW_clr_dir'] = -1. * dic_invar['rsutcs_ano_grd']
                            dic_invar['dLW_clr_dir'].setAxisList(AXL3d)
                            dic_invar['dSW_clr_dir'].setAxisList(AXL3d)
    
                            dic_invar['dnet_clr_dir'] = dic_invar['dLW_clr_dir'] + dic_invar['dSW_clr_dir']
                            dic_invar['dnet_clr_dir'].setAxisList(AXL3d)
                
                            # difference between direct and kernel-calculation
                            dic_invar['dLW_cld_dir_sum'] = dic_invar['dLW_cld_dir'] - dic_invar['dLW_cld_sum']
                            dic_invar['dSW_cld_dir_sum'] = dic_invar['dSW_cld_dir'] - dic_invar['dSW_cld_sum']
                            dic_invar['dLW_cld_dir_sum'].setAxisList(AXL3d)
                            dic_invar['dSW_cld_dir_sum'].setAxisList(AXL3d)
            
                            dic_invar['dLW_clr_dir_sum'] = dic_invar['dLW_clr_dir'] - dic_invar['dLW_clr_sum']
                            dic_invar['dSW_clr_dir_sum'] = dic_invar['dSW_clr_dir'] - dic_invar['dSW_clr_sum']
                            dic_invar['dLW_clr_dir_sum'].setAxisList(AXL3d)
                            dic_invar['dSW_clr_dir_sum'].setAxisList(AXL3d)
                
                            dic_invar['dnet_cld_dir_sum'] = dic_invar['dLW_cld_dir_sum'] + dic_invar['dSW_cld_dir_sum']
                            dic_invar['dnet_cld_dir_sum'].setAxisList(AXL3d)
            
                            dic_invar['dnet_clr_dir_sum'] = dic_invar['dLW_clr_dir_sum'] + dic_invar['dSW_clr_dir_sum']
                            dic_invar['dnet_clr_dir_sum'].setAxisList(AXL3d)
            
            
                            # ---------------------------------------------------------------#
                            # delete unnecessary vars
                            # ---------------------------------------------------------------#
                            delterm = ['rlut_ano_grd','rsut_ano_grd','rlutcs_ano_grd','rsutcs_ano_grd']
                            dic_invar = delete_vars(delterm,dic_invar)
            
                            sub_name = 'Getting Cloudy residual term '
                            print_memory_status(sub_name)
            
                            #############################################################################################################################
                            ################################ move lat-lon feedback and global mean feedback here ########################################
                            # get annual global-mean surface temperature anomaly time series
                            tas_ano_grd_avg = cdutil.averager(dic_invar['tas_ano_grd'],axis='xy',weights='weighted')
                            cdutil.setTimeBoundsMonthly(tas_ano_grd_avg)
                            dic_invar['tas_ano_grd_ann'] = cdutil.YEAR(tas_ano_grd_avg)
            
                            del tas_ano_grd_avg
               
                            # July 6: add global surface temperature anomaly standarized by global-mean sfc temperature anomaly
                            Allname1 = ['tas_ano_grd','dLW_ts','dLW_ta_psum','dLW_ts_clr','dLW_ta_clr_psum','dLW_ta_fxRH_psum','dLW_ta_clr_fxRH_psum',\
                            'dLW_planck_psum','dLW_planck_clr_psum','dLW_planck_fxRH_psum','dLW_planck_clr_fxRH_psum',\
                            'dLW_lapserate_psum','dLW_lapserate_clr_psum','dLW_lapserate_fxRH_psum','dLW_lapserate_clr_fxRH_psum',\
                            'dSW_alb','dSW_alb_clr',\
                            'dSW_cld_dir','dLW_cld_dir',\
                            'dSW_clr_dir','dLW_clr_dir',\
                            'dSW_cld_sum','dLW_cld_sum',\
                            'dSW_clr_sum','dLW_clr_sum',\
                            'dnet_cld_dir','dnet_clr_dir',\
                            'dnet_cld_sum','dnet_clr_sum']
                            Allname2 = ['tas_ano_grd','ts','ta','ts_clr','ta_clr','t_fxRH','t_clr_fxRH',\
                            'planck','planck_clr','planck_fxRH','planck_clr_fxRH',\
                            'lapserate','lapserate_clr','lapserate_fxRH','lapserate_clr_fxRH',\
                            'alb','alb_clr',\
                            'dSW_cld_dir','dLW_cld_dir',\
                            'dSW_clr_dir','dLW_clr_dir',\
                            'dSW_cld_sum','dLW_cld_sum',\
                            'dSW_clr_sum','dLW_clr_sum',\
                            'dnet_cld_dir','dnet_clr_dir',\
                            'dnet_cld_sum','dnet_clr_sum']
            
                            for iname,name in enumerate(Allname1):
                                print('Regresion calculation for',name)
                                newdata1,newdata2,newdata3,newdata4 = get_feedback(exp_cntl,dic_invar[name],dic_invar['tas_ano_grd_ann'],AXL2d)
                                dic_invar[Allname2[iname]+'_gfdbk'] = newdata1
                                dic_invar[Allname2[iname]+'_feedback'] = newdata2
                                dic_invar[Allname2[iname]+'_gforcing'] = newdata3
                                dic_invar[Allname2[iname]+'_forcing'] = newdata4
            
                                del newdata1, newdata2, newdata3, newdata4
            
            
                            # ---------------------------------------------------------------#
                            #  get final t_feedbacks 
                            # ---------------------------------------------------------------#
    
                            outname = ['gfdbk','feedback','gforcing','forcing']
                            for iname,name in enumerate(outname):
                                dic_invar['t_'+name] = dic_invar['ts_'+name] + dic_invar['ta_'+name]
                                dic_invar['t_clr_'+name] = dic_invar['ts_clr_'+name] + dic_invar['ta_clr_'+name]
            
                                if name not in ['feedback','forcing']:
                                    dic_invar['t_'+name].setAxisList(AXL2d)
                                    dic_invar['t_clr_'+name].setAxisList(AXL2d)
            
                            # ---------------------------------------------------------------#
                            # get final LW and SW water vapor feedback
                            # ---------------------------------------------------------------#
                            Allname1 = ['dSW_q_psum','dSW_q_clr_psum','dLW_q_psum','dLW_q_clr_psum',\
                            'dSW_q_fxRH_psum','dSW_q_clr_fxRH_psum','dLW_q_fxRH_psum','dLW_q_clr_fxRH_psum',\
                            'dLW_netRH_psum','dSW_netRH_psum','dLW_clr_netRH_psum','dSW_clr_netRH_psum']
            
                            Allname2 = ['q_sw','q_sw_clr','q_lw','q_lw_clr',\
                            'q_sw_fxRH','q_sw_clr_fxRH','q_lw_fxRH','q_lw_clr_fxRH',\
                            'netRH_lw','netRH_sw','netRH_lw_clr','netRH_sw_clr']
            
                            for iname,name in enumerate(Allname1):
                                print('Regresion calculation for',name)
                                newdata1,newdata2,newdata3,newdata4 = get_feedback(exp_cntl,dic_invar[name],dic_invar['tas_ano_grd_ann'],AXL2d)
                                dic_invar[Allname2[iname]+'_gfdbk'] = newdata1
                                dic_invar[Allname2[iname]+'_feedback'] = newdata2
                                dic_invar[Allname2[iname]+'_gforcing'] = newdata3
                                dic_invar[Allname2[iname]+'_forcing'] = newdata4
            
                                del newdata1, newdata2, newdata3, newdata4
            
                            # ---------------------------------------------------------------#
                            #  get net q_feedbacks and netRH_feedbacks 
                            # ---------------------------------------------------------------#
                            outname = ['gfdbk','feedback','gforcing','forcing']
                            for iname,name in enumerate(outname):
                                dic_invar['q_'+name] = dic_invar['q_sw_'+name] + dic_invar['q_lw_'+name]
                                dic_invar['q_clr_'+name] = dic_invar['q_sw_clr_'+name] + dic_invar['q_lw_clr_'+name]
            
                                dic_invar['q_fxRH_'+name] = dic_invar['q_sw_fxRH_'+name] + dic_invar['q_lw_fxRH_'+name]
                                dic_invar['q_clr_fxRH_'+name] = dic_invar['q_sw_clr_fxRH_'+name] + dic_invar['q_lw_clr_fxRH_'+name]
            
                                dic_invar['netRH_'+name] = dic_invar['netRH_sw_'+name] + dic_invar['netRH_lw_'+name]
                                dic_invar['netRH_clr_'+name] = dic_invar['netRH_sw_clr_'+name] + dic_invar['netRH_lw_clr_'+name]
            
                                if name not in ['feedback','forcing']:
                                    dic_invar['q_'+name].setAxisList(AXL2d)
                                    dic_invar['q_clr_'+name].setAxisList(AXL2d)
                                    dic_invar['q_fxRH_'+name].setAxisList(AXL2d)
                                    dic_invar['q_clr_fxRH_'+name].setAxisList(AXL2d)
                                    dic_invar['netRH_'+name].setAxisList(AXL2d)
                                    dic_invar['netRH_clr_'+name].setAxisList(AXL2d)
            
                            # =============================================================================================#
                            # print summary table here
                            # =============================================================================================#
    
                            print('------ Hi, summary all feedbacks except for cloud feedback------------')
                            print('Planck feedback: ',              dic_invar['planck_feedback'],'W/m2/K')
                            print('Lapse rate feedback: ',          dic_invar['lapserate_feedback'],'W/m2/K')
                            print('Lapse rate + Planck feedback: ', dic_invar['lapserate_feedback']+dic_invar['planck_feedback'],'W/m2/K')
                            print('Temperature feedback: ',         dic_invar['t_feedback'],'W/m2/K')
                            print('Water vapor feedback: ',         dic_invar['q_feedback'],'W/m2/K')
                            print("Surface albedo feedback: ",      dic_invar['alb_feedback'], "W/m2/K")
            
                            print('fixedRH Planck feedback: ',              dic_invar['planck_fxRH_feedback'],'W/m2/K')
                            print('fixedRH Lapse rate feedback: ',          dic_invar['lapserate_fxRH_feedback'],'W/m2/K')
                            print('fixedRH Lapse rate + Planck feedback: ', dic_invar['lapserate_fxRH_feedback']+dic_invar['planck_fxRH_feedback'],'W/m2/K')
                            print('fixedRH Water vapor feedback: ',         dic_invar['q_fxRH_feedback'],'W/m2/K')
                            print('fixedRH RH feedback: ',                  dic_invar['netRH_feedback'],'W/m2/K')
            
                            print('--------clear-sky component----------------------------------------------')
                            print('clr Planck feedback: ',                  dic_invar['planck_clr_feedback'],'W/m2/K')
                            print('clr Lapse rate feedback: ',              dic_invar['lapserate_clr_feedback'],'W/m2/K')
                            print('clr Lapse rate + clr Planck feedback: ', dic_invar['lapserate_clr_feedback']+dic_invar['planck_clr_feedback'],'W/m2/K')
                            print('clr Temperature feedback: ',             dic_invar['t_clr_feedback'],'W/m2/K')
                            print('clr Water vapor feedback: ',             dic_invar['q_clr_feedback'],'W/m2/K')
                            print("clr Surface albedo feedback: ",          dic_invar['alb_clr_feedback'],"W/m2/K")
            
                            print('fixedRH clr Planck feedback: ',              dic_invar['planck_clr_fxRH_feedback'],'W/m2/K')
                            print('fixedRH clr Lapse rate feedback: ',          dic_invar['lapserate_clr_fxRH_feedback'],'W/m2/K')
                            print('fixedRH clr Lapse rate + Planck feedback: ', dic_invar['lapserate_clr_fxRH_feedback'] + dic_invar['planck_clr_fxRH_feedback'],'W/m2/K')
                            print('fixedRH clr Water vapor feedback: ',         dic_invar['q_clr_fxRH_feedback'],'W/m2/K')
                            print('fixedRH clr RH feedback: ',                  dic_invar['netRH_clr_feedback'],'W/m2/K')
            
                            print('--------------- quick check here --------------------------')
                            print('fixedRH (Lapse rate + Planck) - (Lapse rate + Planck): ',dic_invar['lapserate_fxRH_feedback'] + dic_invar['planck_fxRH_feedback'] - dic_invar['lapserate_feedback'] - dic_invar['planck_feedback'],'W/m2/K')
                            print('fixedRH WV: ',dic_invar['q_fxRH_feedback'],'W/m2/K')
                            print('fixedRH Ta: ',dic_invar['t_fxRH_feedback'],'W/m2/K')
            
                            print('clr fixedRH (Lapse rate + Planck) - (Lapse rate + Planck): ',dic_invar['lapserate_clr_fxRH_feedback'] + dic_invar['planck_clr_fxRH_feedback'] - dic_invar['lapserate_clr_feedback'] - dic_invar['planck_clr_feedback'],'W/m2/K')
                            print('clr fixedRH WV: ',dic_invar['q_clr_fxRH_feedback'],'W/m2/K')
                            print('clr fixedRH Ta: ',dic_invar['t_clr_fxRH_feedback'],'W/m2/K')
            
                            print('Lapse rate + Planck + WV feedback: ',dic_invar['lapserate_feedback'] + dic_invar['planck_feedback'] + dic_invar['q_feedback'],'W/m2/K')
                            print('fixedRH Lapse rate + Planck + WV feedback: ',dic_invar['lapserate_fxRH_feedback'] + dic_invar['planck_fxRH_feedback'] + dic_invar['netRH_feedback'],'W/m2/K')
            
                            print('clr Lapse rate + Planck + WV feedback: ',dic_invar['lapserate_clr_feedback'] + dic_invar['planck_clr_feedback'] + dic_invar['q_clr_feedback'],'W/m2/K')
                            print('fixedRH clr Lapse rate + Planck + WV feedback: ',dic_invar['lapserate_clr_fxRH_feedback'] + dic_invar['planck_clr_fxRH_feedback'] + dic_invar['netRH_clr_feedback'],'W/m2/K')
            
                            print('sum of all clear-sky sw feedbacks',dic_invar['dSW_clr_sum_feedback'],'W/m2/K')
                            print('sum of all clear-sky lw feedbacks',dic_invar['dLW_clr_sum_feedback'],'W/m2/K')
            
                            print('TOA direct clear-sky sw radiation feedback',dic_invar['dSW_clr_dir_feedback'],'W/m2/K')
                            print('TOA direct clear-sky lw radiation feedback',dic_invar['dLW_clr_dir_feedback'],'W/m2/K')
            
                            sub_name = 'Getting final regression on tas_ano '
                            print_memory_status(sub_name)
             
            
                            # =============================================================================================#
                            # clear-sky test
                            # =============================================================================================#
    
                            # ---------------------------------------------------------------#
                            # plot lat-lon figure as Shell et al. (2006) Figure2
                            # get lat-lon data
                            # ---------------------------------------------------------------#
           
                            Allname1 = ['dLW_t_psum_mask','dSW_alb_mask','dLW_q_psum_mask','dSW_q_psum_mask','dLW_adj','dSW_adj','SWCRE_ano_grd','LWCRE_ano_grd','netCRE_ano_grd','SWCRE_ano_grd_adj','LWCRE_ano_grd_adj','net_adj','netCRE_ano_grd_adj','dLW_cld_dir_sum','dSW_cld_dir_sum','dnet_cld_dir_sum','dLW_clr_dir_sum','dSW_clr_dir_sum','dnet_clr_dir_sum']
                            Allname2 = Allname1
                            for iname,name in enumerate(Allname1):
                                print('Regresion calculation for',name)
                                newdata1,newdata2,newdata3,newdata4 = get_feedback(exp_cntl,dic_invar[name],dic_invar['tas_ano_grd_ann'],AXL2d)
                                dic_invar[Allname2[iname]+'_gfdbk'] = newdata1
                                dic_invar[Allname2[iname]+'_feedback'] = newdata2
                                dic_invar[Allname2[iname]+'_gforcing'] = newdata3
                                dic_invar[Allname2[iname]+'_forcing'] = newdata4
            
                                del newdata1, newdata2, newdata3, newdata4
            
            
                            sub_name = 'Getting clear-sky linearity test '
                            print_memory_status(sub_name)
             
                            # ---------------------------------------------------------------#
                            # Plotting Lat-Lon feedback test 
                            # ---------------------------------------------------------------#
                            if plotting:
    
                                fig = plt.figure(figsize=(18,9))
                                fh = 10
                                plt.suptitle('clear-sky linear test',fontsize=fh)
                                
                                names = ['ta_gfdbk','ts_gfdbk','t_gfdbk','q_lw_gfdbk','alb_gfdbk','q_sw_gfdbk','t_clr_gfdbk','q_lw_clr_gfdbk','alb_clr_gfdbk','q_sw_clr_gfdbk','dLW_clr_sum_gfdbk','dSW_clr_sum_gfdbk','dLW_clr_dir_gfdbk','dSW_clr_dir_gfdbk','dLW_clr_dir_sum_gfdbk','dSW_clr_dir_sum_gfdbk']
                
                                # lat-lon
                                for n,name in enumerate(names):
                                    if 'cld_psum' in name:
                                        DATA = MV.average(dic_invar[name],axis=0)
                                    else:
                                        DATA = dic_invar[name]
                
                                    if name == 'dLW_t_gfdbk' or name == 'dLW_t_clr_gfdbk' or name == 'dLW_ta_psum_gfdbk' or name == 'dLW_ts_gfdbk':
                                        bounds = np.arange(-5,0.5,0.5)
                                    elif name == 'dLW_q_psum_gfdbk' or name == 'dLW_q_clr_psum_gfdbk':
                                        bounds = np.arange(0,5.5,0.5)
                                    elif 'dir_sum' in name:
                                        bounds = np.arange(-0.5,0.6,0.1)
                                    elif name == 'dSW_clr_sum':
                                        bounds = np.arange(-9,9.5,1)
                                    else:
                                        bounds = np.arange(-2.5,3.0,0.5)
                
                                    cmap = pl.cm.RdBu_r
                                    bounds2 = np.append(np.append(-500,bounds),500)
                                    norm = mpl.colors.BoundaryNorm(bounds2,cmap.N)
                
                                    ax1 = fig.add_subplot(5,4,n+1,projection=ccrs.Robinson(central_longitude=180.))
                                    im1 = ax1.contourf(lon,lat,DATA,bounds,transform=ccrs.PlateCarree(),cmap=cmap,norm=norm,extend='both')
                                    ax1.coastlines()
                                    ax1.set_global()
                                    DATA.setAxisList(AXL2d)
                                    avgDATA = cdutil.averager(DATA,axis='xy',weights='weighted')
                                    pl.title(name + ' ['+str(np.round(avgDATA,3))+' '+str(np.round(np.min(DATA),3))+' '+str(np.round(np.max(DATA),3))+']',fontsize=fh)
                                    cb = plt.colorbar(im1,orientation='vertical',drawedges=True,ticks=bounds)
                                    del avgDATA
                
                                    cb.set_label('W/m$^2$/K')
                
                                plt.tight_layout()
                                plt.savefig(figdir+'lat-lon-clear-sky-linear-test-'+phase+'-'+exp_new+'-'+used_models[imod]+'-'+AddString+'_'+AddTag+'.png',bbox_inches='tight')
                                fig.clf()
                                plt.close()
                                gc.collect()
                
                                # ---------------------------------------------------------------#
                                ## Plotting zonal-mean feedback test 
                                # ---------------------------------------------------------------#
        
                                fig = plt.figure(figsize=(18,9))
                                fh = 10
                                plt.suptitle('clear-sky linear test',fontsize=fh)
                 
                                names1 = ['dLW_clr_sum_gfdbk','dSW_clr_sum_gfdbk']
                                names2 = ['dLW_clr_dir_gfdbk','dSW_clr_dir_gfdbk']
                                for n,name in enumerate(names1):
                                    DATA1 = MV.average(dic_invar[names1[n]],axis=1)
                                    DATA2 = MV.average(dic_invar[names2[n]],axis=1)
                
                                    ax2 = fig.add_subplot(1,2,n+1)
                                    plt.plot(lat,np.array(DATA1),label=names1[n],color='red')
                                    plt.plot(lat,np.array(DATA2),label=names2[n],color='blue')
                                    plt.plot(lat,np.array(DATA2)*1.15,linestyle='--',color='black')
                                    plt.plot(lat,np.array(DATA2)*0.85,linestyle='--',color='black')
                                    pl.title(names1[n]+' and '+names2[n],fontsize=fh)
                                    plt.legend()
                
                                    del DATA1, DATA2
                
                                plt.tight_layout()
                                plt.savefig(figdir+'zonalmean-clear-sky-linear-test-'+phase+'-'+exp_new+'-'+used_models[imod]+'-'+AddString+'_'+AddTag+'.png',bbox_inches='tight')
        
                                fig.clf()
                                plt.close()
                                gc.collect()
                
                                # ---------------------------------------------------------------#
                                # Plotting lat-lon feedback and clear-sky difference test 
                                # ---------------------------------------------------------------#
                                fig = plt.figure(figsize=(18,9))
                                fh = 10
                                plt.suptitle('clear-sky linear test',fontsize=fh)
                
                                names = ['dLW_t_psum_mask','dSW_alb_mask','dLW_q_psum_mask','dSW_q_psum_mask','dLW_adj','dSW_adj','SWCRE_ano_grd','LWCRE_ano_grd','netCRE_ano_grd','SWCRE_ano_grd_adj','LWCRE_ano_grd_adj','net_adj','netCRE_ano_grd_adj','dLW_cld_sum','dSW_cld_sum','dLW_cld_dir','dSW_cld_dir','dLW_cld_dir_sum','dSW_cld_dir_sum','dnet_cld_dir_sum']
                                # lat-lon
                                for n,name in enumerate(names):
                                    DATA = dic_invar[name+'_gfdbk']
                
                                    bounds = np.arange(-10,12,2)/4.
                                    cmap = pl.cm.RdBu_r
                                    bounds2 = np.append(np.append(-500,bounds),500)
                                    norm = mpl.colors.BoundaryNorm(bounds2,cmap.N)
                
                                    ax1 = fig.add_subplot(5,4,n+1,projection=ccrs.Robinson(central_longitude=180.))
                                    im1 = ax1.contourf(lon,lat,DATA,bounds,transform=ccrs.PlateCarree(),cmap=cmap,norm=norm,extend='both')
                                    ax1.coastlines()
                                    ax1.set_global()
                                    DATA.setAxisList(AXL2d)
                                    avgDATA = cdutil.averager(DATA,axis='xy',weights='weighted')
                                    pl.title(name + ' ['+str(np.round(avgDATA,3))+']',fontsize=fh)
                                    cb = plt.colorbar(im1,orientation='vertical',drawedges=True,ticks=bounds)
                                    del avgDATA
                
                                    cb.set_label('W/m$^2$/K')
                                plt.tight_layout()
                
                                plt.savefig(figdir+'lat-lon-adjust-CRE-'+phase+'-'+exp_new+'-'+used_models[imod]+'-'+AddString+'_'+AddTag+'.png',bbox_inches='tight')
                                
                                fig.clf()
                                plt.close()
                                gc.collect()
            
                            # ===============================================================#
                            # All plotting tests are done here. start deleting vars ....
                            # ===============================================================#
    
                            delterm = ['dLW_t_clr_psum','dLW_t_psum','dLW_q_clr_psum','dLW_q_psum','dSW_q_clr_psum','dSW_q_psum','dSW_alb_clr','dSW_alb',\
                            'dLW_t_psum_mask','dLW_q_psum_mask','dSW_q_psum_mask','dSW_alb_mask',\
                            'SWCRE_ano_grd','LWCRE_ano_grd','netCRE_ano_grd','dSW_adj','dLW_adj','net_adj',\
                            'LWCRE_ano_grd_adj','SWCRE_ano_grd_adj','netCRE_ano_grd_adj',\
                            'dLW_cld_dir','dLW_cld_sum','dSW_cld_dir','dSW_cld_sum','dLW_cld_dir_sum','dSW_cld_dir_sum','dnet_cld_dir_sum',\
                            'rsutcs_ano_grd','rsdt_ano_grd','rlutcs_ano_grd']
                            dic_invar = delete_vars(delterm,dic_invar)
            
                            sub_name = 'Getting final plotting '
                            print_memory_status(sub_name)
            
                            # ===============================================================#
                            # Final output to CSV files and NC files
                            # ===============================================================#
    
                            # ---------------------------------------------------------------#
                            # output global mean values 
                            # ---------------------------------------------------------------#
           
                            invars1 = ['t','planck','lapserate','q','alb',\
                            'dLW_adj','dSW_adj','net_adj',\
                            'SWCRE_ano_grd','LWCRE_ano_grd','netCRE_ano_grd',\
                            'SWCRE_ano_grd_adj','LWCRE_ano_grd_adj','netCRE_ano_grd_adj',\
                            'dLW_cld_dir_sum','dSW_cld_dir_sum','dnet_cld_dir_sum',\
                            't_clr','planck_clr','lapserate_clr','q_clr','alb_clr',\
                            'q_sw','q_lw',\
                            'q_sw_clr','q_lw_clr',\
                            'planck_fxRH','lapserate_fxRH','netRH',\
                            'planck_clr_fxRH','lapserate_clr_fxRH','netRH_clr',\
                            'dLW_clr_sum','dSW_clr_sum','dnet_clr_sum',\
                            'dLW_clr_dir','dSW_clr_dir','dnet_clr_dir',\
                            'dLW_cld_sum','dSW_cld_sum','dnet_cld_sum',\
                            'dLW_cld_dir','dSW_cld_dir','dnet_cld_dir',\
                            'dLW_clr_dir_sum','dSW_clr_dir_sum','dnet_clr_dir_sum']
            
                            outputs = []
                            outputs_forcing = []
                            for ivar,svar in enumerate(invars1):
                                outputs.append(dic_invar[svar+'_feedback'])
                                outputs_forcing.append(dic_invar[svar+'_forcing'])
    
                            final_index = ['T','Planck','LR','WV','ALB','dLW_adj','dSW_adj','dnet_adj','SWCRE','LWCRE','netCRE','SWCRE_adj','LWCRE_adj','netCRE_adj','LW_resd','SW_resd','net_resd','T_clr','Planck_clr','LR_clr','WV_clr','ALB_clr','WV_SW','WV_LW','WV_clr_SW','WV_clr_LW','Planck_fxRH','LR_fxRH','RH','Planck_clr_fxRH','LR_clr_fxRH','RH_clr','LW_clr_sum','SW_clr_sum','net_clr_sum','LW_clr_dir','SW_clr_dir','net_clr_dir','LW_cld_sum','SW_cld_sum','net_cld_sum','LW_cld_dir','SW_cld_dir','net_cld_dir','LW_clr_resd','SW_clr_resd','net_clr_resd']

                            # output variables to csv file
                            df_all = pd.DataFrame(outputs,index=final_index,columns=[used_models[imod]])
                            df_all.to_csv(outdir+'FDBK_'+phase+'_'+exp_new+'_'+used_models[imod]+'_'+AddString+'_'+AddTag+'.csv')

                            df_all_forcing = pd.DataFrame(outputs_forcing,index=final_index,columns=[used_models[imod]])
                            df_all_forcing.to_csv(outdir+'FDBK_forcing_'+phase+'_'+exp_new+'_'+used_models[imod]+'_'+AddString+'_'+AddTag+'.csv')

                            del final_index, df_all, df_all_forcing, outputs, outputs_forcing
            
                            # ---------------------------------------------------------------#
                            # output spatial pattern
                            # ---------------------------------------------------------------#
                            invars2 = invars1 + ['tas_ano_grd']
             
                            out1 = cdms.open(outdir+'lat-lon-gfdbk-'+phase+'-'+exp_new+'-'+used_models[imod]+'_'+AddString+'_'+AddTag+'.nc','w')
                            out2 = cdms.open(outdir+'lat-lon-gfdbk-forcing-'+phase+'-'+exp_new+'-'+used_models[imod]+'_'+AddString+'_'+AddTag+'.nc','w')
              
                            for ivar,svar in enumerate(invars2):
                                print('Writing',svar,'_gfdbk to file')
                                DATA = dic_invar[svar+'_gfdbk']
                                DATA.id = svar
                                out1.write(DATA)
                                del DATA,dic_invar[svar+'_gfdbk']
                
                            for ivar,svar in enumerate(invars2): 
                                print('Writing',svar,'_gforcing to file')
                                DATA = dic_invar[svar+'_gforcing']
                                DATA.id = svar
                                out2.write(DATA)
                                del DATA, dic_invar[svar+'_gforcing']
             
                            out1.close()
                            out2.close()
                            del out1,out2
            
                            sub_name = 'All things for this model:'+used_models[imod]
                            print_memory_status(sub_name)

                            # ====================================================================#
                            # mannually delete unnecessary variables here
                            # ====================================================================#
                            print('dic_invar.keys() are',dic_invar.keys())
                            del dic_invar,AXL2d, AXL3d, AXL4d,delterm
                            del iname,name,outname,Allname1,Allname2,invar1,invar2,invars1,invars2
                            del lat,lon,time,lev,nlat,nlon,ntime,nlev,outvar,ovar,svar,these_vars,ivar
                            del kernel_dir,sub_name
                            gc.collect(2)
    
                            # ====================================================================#
                            # print final memory check and left variables 
                            # ====================================================================#
            
                            dir2 = dir()
                            print('--------------------------------------------------')
                            print('final scope variables are ',dir())
                            print('--------------------------------------------------')
            
                            for xx in dir2:
                                if xx not in ['dir1','dir2','imod','used_models'] and xx not in dir1:
                                    print('additional variable is',xx,'getrefcount is',sys.getrefcount(locals()[xx]))
#                                    del locals()[xx]
            
                            sub_name = 'Also delete those unnecesasry variables to skim the code!'
                            print_memory_status(sub_name)
    
                            print('final variables afte deleting ', dir())
    
                            print('==============================================================')
                            print(used_models[imod], 'feedback analysis is done! Congrats!')
                            sub_name = used_models[imod]+' feedback analysis'
                            print_memory_status(sub_name)
                            print('==============================================================')
                            exit()
    
                        print('==============================================================')
                        print(exp_new,' processing is done!')
                        sub_name = exp_new+ ', All models  processing'
                        print_memory_status(sub_name)
                        print('==============================================================')
                        exit()
    
                    print('==============================================================')
                    print(samip,' processing is done!')
                    sub_name = samip+', all '+exp_new+', all models processing'
                    print_memory_status(sub_name)
                    print('==============================================================')
                    exit()
    
                print('==============================================================')
                print(phase,' processing is done!')
                sub_name = phase+', All amips and all models processing'
                print_memory_status(sub_name)
                print('==============================================================')
    
            exit() 

# call main function 
main()



