#usr/bin/env python
# coding: utf-8

"""
# This script demonstrates how to compute the cloud feedback using cloud radiative kernels for a 
# short (2-year) period of MPI-ESM-LR using the difference between amipFuture and amip runs.
# One should difference longer periods for more robust results -- these are just for demonstrative purposes

# Additionally, this script demonstrates how to compute the Zelinka et al (2013) decomposition of cloud
feedback into components due to changes in cloud amount, altitude, optical depth, and a residual, and 
does this separately for low and non-low clouds following Zelinka et al (2016).

# Data that are used in this script:
# 1. model clisccp field
# 2. model rsuscs field
# 3. model rsdscs field
# 4. model tas field
# 5. cloud radiative kernels

# This script written by Mark Zelinka (zelinka1@llnl.gov) on 30 October 2018

References:
Zelinka, M. D., S. A. Klein, and D. L. Hartmann, 2012: Computing and Partitioning Cloud Feedbacks Using 
    Cloud Property Histograms. Part I: Cloud Radiative Kernels. J. Climate, 25, 3715-3735. 
    doi:10.1175/JCLI-D-11-00248.1.

Zelinka, M. D., S. A. Klein, and D. L. Hartmann, 2012: Computing and Partitioning Cloud Feedbacks Using 
    Cloud Property Histograms. Part II: Attribution to Changes in Cloud Amount, Altitude, and Optical Depth. 
    J. Climate, 25, 3736-3754. doi:10.1175/JCLI-D-11-00249.1.

Zelinka, M.D., S.A. Klein, K.E. Taylor, T. Andrews, M.J. Webb, J.M. Gregory, and P.M. Forster, 2013: 
    Contributions of Different Cloud Types to Feedbacks and Rapid Adjustments in CMIP5. 
    J. Climate, 26, 5007-5027. doi: 10.1175/JCLI-D-12-00555.1.
    
Zelinka, M. D., C. Zhou, and S. A. Klein, 2016: Insights from a Refined Decomposition of Cloud Feedbacks, 
    Geophys. Res. Lett., 43, 9259-9269, doi:10.1002/2016GL069917.  
"""

# July 1, 2020 modified by Yi Qin
# change into a function to be used for E3SM; 
# change 'exec' into dictionary type
# Aug 12, 2021 -- make it work on FISCCP1_COSP not only isccp variable name.
# Aug 21, 2021 -- change caselist if 'coupled' in case_stamp to run linear regression.
# Mar 11, 2023 -- replace CDAT
 
#IMPORT STUFF:
#=====================
import xarray as xr
import numpy as np
import pylab as pl
import matplotlib as mpl
mpl.use('Agg')
import pandas as pd
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import os
import sys
sys.path.append('../')
import cases_lookup as CL
from PlotDefinedFunction import linearregression_nd,area_averager,weighted_annual_mean, save_big_dataset
from get_mip_data import read_mip_data,read_amip_data,read_pickle,write_pickle,read_e3sm_data,read_e3smdiag_data

###########################################################################
# MAIN ROUTINE FOLLOWS
##########################################################################
def CloudRadKernel_MIP(direc_kernel,case_stamp,outdir,figdir,filenames,tslice):

    outfile_gm_sw = outdir+'decomp_global_mean_sw_'+case_stamp+'.csv'
    outfile_gm_lw = outdir+'decomp_global_mean_lw_'+case_stamp+'.csv'
    outfile_map = outdir+'global_cloud_feedback_'+case_stamp+'.nc'

    if os.path.isfile(outfile_gm_sw) and os.path.isfile(outfile_gm_lw) and os.path.isfile(outfile_map):
        print('CloudRadKernel result is already there.')
        return

    # Read E3SM data 
    Vars = ['clisccp','rsuscs','rsdscs','tas']
    dics_invar = read_mip_data(Vars,filenames,tslice)

    nyears = int(dics_invar[list(dics_invar.keys())[0]].shape[0]/12)

    calculation(direc_kernel,dics_invar,nyears,outfile_gm_sw,outfile_gm_lw,outfile_map,case_stamp)

    return
 
def CloudRadKernel_E3SM(direc_kernel,direc_data,case_stamp,yearS,yearE,fname1,fname2,outdir,figdir,UseE3smDiagOutput=False,grd_info=None,num_years_per_file=None):

    outfile_gm_sw = outdir+'decomp_global_mean_sw_'+case_stamp+'.csv'
    outfile_gm_lw = outdir+'decomp_global_mean_lw_'+case_stamp+'.csv'
    outfile_map = outdir+'global_cloud_feedback_'+case_stamp+'.nc'

    if os.path.isfile(outfile_gm_sw) and os.path.isfile(outfile_gm_lw) and os.path.isfile(outfile_map):
        print('CloudRadKernel result is already there.')
        return

    # Read E3SM data 
    Vars = ['clisccp','rsuscs','rsdscs','tas']
    nyears = int(yearE-yearS+1)

    if UseE3smDiagOutput:
        dics_invar = read_e3smdiag_data(Vars,direc_data,case_stamp,yearS,yearE,fname1,fname2,grd_info,num_years_per_file)
    else:
        dics_invar = read_e3sm_data(Vars,direc_data,case_stamp,yearS,yearE,fname1,fname2)

    calculation(direc_kernel,dics_invar,nyears,outfile_gm_sw,outfile_gm_lw,outfile_map,case_stamp)

    return
    
# ========================================================================================================
def read_kernel(direc_kernel,nyears):

    # Load in the Zelinka et al 2012 kernels:
    f=xr.open_dataset(direc_kernel+'cloud_kernels2.nc')
    LWkernel=f['LWkernel']
    SWkernel=f['SWkernel']
    f.close()
    
    albcs=np.arange(0.0,1.5,0.5) # the clear-sky albedos over which the kernel is computed
    
    # LW kernel does not depend on albcs, just repeat the final dimension over longitudes:
    LWkernel_map=xr.DataArray(np.transpose(np.tile(np.tile(LWkernel[:,:,:,:,0],(1,nyears,1,1,1)),(144,1,1,1,1)),[1,2,3,4,0]))

    # Define the cloud kernel axis attributes
    lats=LWkernel.lat.values
    lons=np.arange(1.25,360,2.5)

    return lats,lons,LWkernel_map,SWkernel


# ========================================================================================================
def cal_sfc_albedo(dics_invar):
    # Calculate surface albedo 
    albcs1=xr.DataArray(dics_invar['rsuscs_pi']/dics_invar['rsdscs_pi'].where((~np.isnan(dics_invar['rsdscs_pi'])) & (dics_invar['rsdscs_pi']!=0)), coords=dics_invar['rsuscs_pi'].coords)
    albcs1=xr.where(albcs1>1.,1,albcs1) # where(condition, x, y) is x where condition is true, y otherwise
    albcs1=xr.where(albcs1<0.,0,albcs1)
    del dics_invar['rsuscs_pi'], dics_invar['rsdscs_ab']

    dics_invar['albcs1'] = albcs1

    # Compute tas anomaly and global annual mean 
    anomtas = dics_invar['tas_ano']

    # Define new time coordinate
    newtime = pd.date_range(start='1850-01-01', periods=dics_invar['tas_pi'].shape[0], freq='MS')
    anomtas = anomtas.assign_coords({'time':("time",newtime)})

    # Get annual mean 
    anomtas_ann = weighted_annual_mean(anomtas.time,anomtas)

    # Get global annual mean 
    anomtas_ann_gm = area_averager(anomtas_ann)
    
    # Compute global annual mean tas anomalies
    avgdtas = area_averager(anomtas_ann.mean(axis=0)) # (scalar)
    print('avgdtas = ',avgdtas.values)

    dics_invar['anomtas_ann_gm'] = avgdtas 

    return dics_invar 

# ========================================================================================================
def calculation(direc_kernel,dics_invar,nyears,outfile_gm_sw,outfile_gm_lw,outfile_map,case_stamp):

    # Calculate surface albedo and global climo surface temperature anomaly
    dics_invar = cal_sfc_albedo(dics_invar)
    print('Surface albedo calculation is done.')

    # Read kernel and its coordinates 
    lats,lons,LWkernel_map,SWkernel = read_kernel(direc_kernel,nyears)
    print('Reading kernel data is done.')

    # Regrid model data to kernel's lats and lons
    dics_in = {}
    for svar in dics_invar.keys():
        if len(dics_invar[svar].shape) > 2: 
            dics_in[svar]  = dics_invar[svar].interp(lat=lats,lon=lons)
        else:
            dics_in[svar] = dics_invar[svar] 
    del dics_invar 
    print('Regrid model data to kernel grid is done.')

    # Use control albcs to map SW kernel to appropriate longitudes
    newtime = pd.date_range(start='1850-01-01', periods=dics_in['tas_pi'].shape[0], freq='MS')
    dics_in['albcs1'] = dics_in['albcs1'].assign_coords({'time':("time",newtime)})
    avgalbcs1 = dics_in['albcs1'].groupby('time.month').mean('time')
    SWkernel_map_tmp = map_SWkern_to_lon(SWkernel,avgalbcs1)
    SWkernel_map = xr.DataArray(np.tile(SWkernel_map_tmp,(nyears,1,1,1,1)), coords=dics_in['clisccp_pi'].coords)
    del SWkernel_map_tmp
    print('Map SW kernel is done.')

    # The sun is down if every bin of the SW kernel is zero:
    sundown=np.ma.sum(np.ma.sum(SWkernel_map,axis=2),axis=1)  #12*nyears,90,144
    night=np.where(sundown==0)
    print("Data processing is done.")
    print()

    # Prepare input data for later feedback calculation 
    dics_in['dcld_dT'] = dics_in['clisccp_ano'] 
    dics_in['c1'] = dics_in['clisccp_pi']
    dics_in['Klw'] = LWkernel_map
    dics_in['Ksw'] = SWkernel_map
    dics_in['night'] = night
    dics_in['anomtas_ann_gm'] = dics_in['anomtas_ann_gm']

    # Compute cloud feedbacks and their breakdown into components
    print('Start computing cloud feedbacks ...')

    # Define a python dictionary containing the sections of the histogram to consider
    # These are the same as in Zelinka et al, GRL, 2016
    sections = ['ALL','HI680','LO680']
    Psections=[slice(0,7),slice(2,7),slice(0,2)]
    sec_dic=dict(zip(sections,Psections))
    
    df_sw_all = pd.DataFrame()
    df_lw_all = pd.DataFrame()

    dic_out1 = {}
    for sec in sections:
        print ('Using '+sec+' CTP bins')
        choose=sec_dic[sec]
        LC = len(np.ones(100)[choose])
    
        ## Preallocation of arrays:
        #LWcld_tot=nanarray((12*nyears,90,144))
        #LWcld_amt=nanarray((12*nyears,90,144))
        #LWcld_alt=nanarray((12*nyears,90,144))
        #LWcld_tau=nanarray((12*nyears,90,144))
        #LWcld_err=nanarray((12*nyears,90,144))
        #SWcld_tot=nanarray((12*nyears,90,144))
        #SWcld_amt=nanarray((12*nyears,90,144))
        #SWcld_alt=nanarray((12*nyears,90,144))
        #SWcld_tau=nanarray((12*nyears,90,144))
        #SWcld_err=nanarray((12*nyears,90,144))
        #dc_star=nanarray((12*nyears,7,LC,90,144))
        #dc_prop=nanarray((12*nyears,7,LC,90,144))
    
        #for mm in range(12*nyears):
        #    dcld_dT = dics_in['dcld_dT'][mm,:,choose,:]
        #    c1 = dics_in['c1'][mm,:,choose,:]
        #    c2 = c1 + dcld_dT
        #    Klw = dics_in['Klw'][mm,:,choose,:]
        #    Ksw = dics_in['Ksw'][mm,:,choose,:]

        #    print('mm=',mm, dcld_dT.shape, c1.shape, c2.shape, Klw.shape, Ksw.shape) 

        #    # The following performs the amount/altitude/optical depth decomposition of
        #    # Zelinka et al., J Climate (2012b), as modified in Zelinka et al., J. Climate (2013)
        #    (LWcld_tot[mm,:],LWcld_amt[mm,:],LWcld_alt[mm,:],LWcld_tau[mm,:],LWcld_err[mm,:],SWcld_tot[mm,:],SWcld_amt[mm,:],SWcld_alt[mm,:],SWcld_tau[mm,:],SWcld_err[mm,:],dc_star[mm,:],dc_prop[mm,:]) = KT_decomposition_4D(c1,c2,Klw,Ksw)

        dcld_dT = dics_in['dcld_dT'][:,:,choose,:]
        c1      = dics_in['c1'][:,:,choose,:]
        c2      = c1 + dcld_dT
        Klw     = dics_in['Klw'][:,:,choose,:]
        Ksw     = dics_in['Ksw'][:,:,choose,:]

        output = KT_decomposition_general(c1, c2, Klw, Ksw)

        LWcld_tot, LWcld_amt, LWcld_alt, LWcld_tau, LWcld_err, SWcld_tot, SWcld_amt, SWcld_alt, SWcld_tau, SWcld_err  = output['LWcld_tot'], output['LWcld_amt'], output['LWcld_alt'], output['LWcld_tau'], output['LWcld_err'], output['SWcld_tot'], output['SWcld_amt'], output['SWcld_alt'], output['SWcld_tau'], output['SWcld_err']


        # Set the SW cloud feedbacks to zero in the polar night
        # Do this since they may come out of previous calcs as undefined, but should be zero:
        night = dics_in['night']
        SWcld_tot[night]=0
        SWcld_amt[night]=0
        SWcld_alt[night]=0
        SWcld_tau[night]=0
        SWcld_err[night]=0

        AX3 = dics_in['c1'][:,0,0,:].coords #[coord_time,coord_lats,coord_lons]

        # July 7, 2020 save variables into dictionary
        names = [\
        'SWcld_tot','SWcld_amt','SWcld_alt','SWcld_tau','SWcld_err',\
        'LWcld_tot','LWcld_amt','LWcld_alt','LWcld_tau','LWcld_err',\
        ]
        variables = [\
        SWcld_tot,SWcld_amt,SWcld_alt,SWcld_tau,SWcld_err,\
        LWcld_tot,LWcld_amt,LWcld_alt,LWcld_tau,LWcld_err,\
        ]

        dic_all = {}
        for n,name in enumerate(names):
            DATA_anom = variables[n].filled(fill_value=np.nan) # convert masked array to array with nan as missing value
            DATA_anom = xr.DataArray(variables[n], coords=AX3)
            DATA_am = weighted_annual_mean(DATA_anom.time,DATA_anom)

            if 'coupled' in case_stamp:
                slope,intercept = linearregression_nd(DATA_am, x = np.reshape(dics_in['anomtas_ann_gm'].values, (dics_in['anomtas_ann_gm'].shape[0],1,1)))
                #print('slope=',slope.shape, np.min(slope), np.max(slope), 'intercept=',intercept.shape, np.min(intercept), np.max(intercept))
            else:
                slope = DATA_anom.groupby('time.month').mean('time').mean(axis=0)/dics_in['anomtas_ann_gm'].mean()

            dic_out1[str(sec)+'_'+str(name)] = slope 
            dic_out1[str(sec)+'_'+str(name)].attrs['long_name'] = str(sec)+'_'+str(name)

            avgDATA = area_averager(slope).values
            if 'SW' in name:
                df_sw_tmp = pd.DataFrame([[sec,'CTP bins',name,str(np.round(avgDATA,3))]],columns=['type','bin','decomp',case_stamp])
                df_sw_all = pd.concat([df_sw_all,df_sw_tmp])
            elif 'LW' in name:
                df_lw_tmp = pd.DataFrame([[sec,'CTP bins',name,str(np.round(avgDATA,3))]],columns=['type','bin','decomp',case_stamp])
                df_lw_all = pd.concat([df_lw_all,df_lw_tmp])

   
    print(df_lw_all.head())
    print(df_sw_all.head())
    
    # Save global mean values
    df_lw_all.to_csv(outfile_gm_lw)
    df_sw_all.to_csv(outfile_gm_sw)

    # Save nc file
    comment = 'created by cal_CloudRadKernel_E3SM.py; Author: Yi Qin (yi.qin@pnnl.gov)'
    save_big_dataset(dic_out1,outfile_map,comment)

###########################################################################
# HELPFUL FUNCTIONS FOLLOW
###########################################################################

def YEAR(data):
    """
    Compute annual means without forcing it to be Jan through Dec
    """
    A=data.shape[0]
    anndata0=nanarray(data.shape)
    cnt=-1
    for i in np.arange(0,A,12):
        if len(data[i:i+12])==12: # only take full 12-month periods
            cnt+=1
            anndata0[cnt] = np.nanmean(data[i:i+12],0)
    B=cnt+1
    anndata = anndata0[:B]
    if type(anndata)!=float:
        anndata  = xr.DataArray(anndata, coords=data[:B*12:12].coords)

    return anndata

###########################################################################
def add_cyclic(data):
    # Add Cyclic point around 360 degrees longitude:
    lons=data.lon[:]
    dx=np.gradient(lons)[-1]
    lons_new = np.arange(0,np.max(lons),dx)
    data2 = data.assign_coords({"lon": lons_new})
    return data2

###########################################################################
def nanarray(vector):
    # this generates a masked array with the size given by vector
    # example: vector = (90,144,28)
    # similar to this=NaN*ones(x,y,z) in matlab
    #this = np.full(vector,np.nan)

    this = np.zeros(vector)
    this = np.ma.masked_where(this==0,this)

    return this

###########################################################################
def map_SWkern_to_lon(Ksw,albcsmap):
    from scipy.interpolate import interp1d
    ## Map each location's clear-sky surface albedo to the correct albedo bin
    # Ksw is size 12,7,7,lats,3
    # albcsmap is size A,lats,lons
    albcs=np.arange(0.0,1.5,0.5) 
    A=albcsmap.shape[0]
    TT=Ksw.shape[1]
    PP=Ksw.shape[2]
    lenlat=Ksw.shape[3]
    lenlon=albcsmap.shape[2]
    SWkernel_map=nanarray((A,TT,PP,lenlat,lenlon))
    for M in range(A):
        MM=M
        while MM>11:
            MM=MM-12
        for LA in range(lenlat):
            #print('M=',M, 'LA=', LA)
            alon=albcsmap[M,LA,:] 
            # interp1d can't handle mask but it can deal with NaN (?)
            alon2=alon
            if np.count_nonzero(~np.isnan(alon2))>1: # at least 1 unmasked value
                if len(np.where(Ksw[MM,:,:,LA,:]>0))==0:
                    SWkernel_map[M,:,:,LA,:] = 0
                else:
                    f = interp1d(albcs,Ksw[MM,:,:,LA,:],axis=2)
                    ynew = f(alon2.data)
                    SWkernel_map[M,:,:,LA,:] = ynew
            else:
                continue

    return SWkernel_map

###########################################################################
def KT_decomposition_4D(c1,c2,Klw,Ksw):

    # this function takes in a (tau,CTP,lat,lon) matrix and performs the 
    # decomposition of Zelinka et al 2013 doi:10.1175/JCLI-D-12-00555.1

    # reshape to be (CTP,tau,lat,lon)
    # This is inefficient but done purely because Mark can't think in tau,CTP space
    c1 = np.ma.transpose(c1.to_masked_array(),(1,0,2,3)) # control cloud fraction histogram
    c2 = np.ma.transpose(c2.to_masked_array(),(1,0,2,3)) # perturbed cloud fraction histogram
    Klw = np.ma.transpose(Klw.to_masked_array(),(1,0,2,3)) # LW Kernel histogram
    Ksw = np.ma.transpose(Ksw.to_masked_array(),(1,0,2,3)) # SW Kernel histogram

    P=c1.shape[0]
    T=c1.shape[1]

    c=c1
    sum_c=np.tile(np.ma.sum(np.ma.sum(c,1),0),(P,T,1,1))                                  # Eq. B2
    dc = c2-c1 
    sum_dc=np.tile(np.ma.sum(np.ma.sum(dc,1),0),(P,T,1,1))
    dc_prop = c*(sum_dc/sum_c)
    dc_star = dc - dc_prop                                                          # Eq. B1

    # LW components
    Klw0 = np.tile(np.ma.sum(np.ma.sum(Klw*c/sum_c,1),0),(P,T,1,1))                       # Eq. B4
    Klw_prime = Klw - Klw0                                                          # Eq. B3
    this=np.ma.sum(Klw_prime*np.tile(np.ma.sum(c/sum_c,0),(P,1,1,1)),1)                   # Eq. B7a
    Klw_p_prime=np.transpose(np.tile(np.tile(this,(1,1,1,1)),(T,1,1,1)),[1,0,2,3])         # Eq. B7b 
    that=np.transpose(np.tile(np.tile(np.ma.sum(c/sum_c,1),(1,1,1,1)),(T,1,1,1)),[1,0,2,3])   # Eq. B8a
    Klw_t_prime = np.tile(np.ma.sum(Klw_prime*that,0),(P,1,1,1))                       # Eq. B8b
    Klw_resid_prime = Klw_prime - Klw_p_prime - Klw_t_prime                         # Eq. B9
    dRlw_true = np.ma.sum(np.ma.sum(Klw*dc,1),0)                                          # LW total
    dRlw_prop = Klw0[0,0,:,:]*sum_dc[0,0,:,:]                                       # LW amount component
    dRlw_dctp = np.ma.sum(np.ma.sum(Klw_p_prime*dc_star,1),0)                             # LW altitude component
    dRlw_dtau = np.ma.sum(np.ma.sum(Klw_t_prime*dc_star,1),0)                             # LW optical depth component
    dRlw_resid = np.ma.sum(np.ma.sum(Klw_resid_prime*dc_star,1),0)                        # LW residual
    dRlw_sum = dRlw_prop + dRlw_dctp + dRlw_dtau + dRlw_resid                       # sum of LW components -- should equal LW total

    # SW components
    Ksw0 = np.tile(np.ma.sum(np.ma.sum(Ksw*c/sum_c,1),0),(P,T,1,1))                       # Eq. B4
    Ksw_prime = Ksw - Ksw0                                                          # Eq. B3
    this=np.ma.sum(Ksw_prime*np.tile(np.ma.sum(c/sum_c,0),(P,1,1,1)),1)                   # Eq. B7a 
    Ksw_p_prime=np.transpose(np.tile(np.tile(this,(1,1,1,1)),(T,1,1,1)),[1,0,2,3])         # Eq. B7b  
    that=np.transpose(np.tile(np.tile(np.ma.sum(c/sum_c,1),(1,1,1,1)),(T,1,1,1)),[1,0,2,3])   # Eq. B8a
    Ksw_t_prime = np.tile(np.ma.sum(Ksw_prime*that,0),(P,1,1,1))                       # Eq. B8b
    Ksw_resid_prime = Ksw_prime - Ksw_p_prime - Ksw_t_prime                         # Eq. B9
    dRsw_true = np.ma.sum(np.ma.sum(Ksw*dc,1),0)                                          # SW total
    dRsw_prop = Ksw0[0,0,:,:]*sum_dc[0,0,:,:]                                       # SW amount component
    dRsw_dctp = np.ma.sum(np.ma.sum(Ksw_p_prime*dc_star,1),0)                             # SW altitude component
    dRsw_dtau = np.ma.sum(np.ma.sum(Ksw_t_prime*dc_star,1),0)                             # SW optical depth component
    dRsw_resid = np.ma.sum(np.ma.sum(Ksw_resid_prime*dc_star,1),0)                        # SW residual
    dRsw_sum = dRsw_prop + dRsw_dctp + dRsw_dtau + dRsw_resid                       # sum of SW components -- should equal SW total

    dc_star = np.ma.transpose(dc_star,(1,0,2,3)) 
    dc_prop = np.ma.transpose(dc_prop,(1,0,2,3)) 

    return (dRlw_true,dRlw_prop,dRlw_dctp,dRlw_dtau,dRlw_resid,dRsw_true,dRsw_prop,dRsw_dctp,dRsw_dtau,dRsw_resid,dc_star,dc_prop)

###########################################################################
def KT_decomposition_general(c1,c2,Klw,Ksw):
    """
    this function takes in a (month,TAU,CTP,lat,lon) matrix and performs the 
    decomposition of Zelinka et al 2013 doi:10.1175/JCLI-D-12-00555.1
    """
    
    # To help with broadcasting, move month axis to the end so that TAU,CTP are first
    c1 = np.array(np.moveaxis(c1.to_masked_array(),0,-1))
    c2 = np.array(np.moveaxis(c2.to_masked_array(),0,-1))
    Klw = np.moveaxis(Klw.to_masked_array(),0,-1)
    Ksw = np.moveaxis(Ksw.to_masked_array(),0,-1)
    
    sum_c=np.ma.sum(np.ma.sum(c1,0),0)                              # Eq. B2
    dc = c2-c1 
    sum_dc=np.ma.sum(np.ma.sum(dc,0),0)
    dc_prop = c1*(sum_dc/sum_c)
    dc_star = dc - dc_prop                                          # Eq. B1

    # LW components
    Klw0 = np.ma.sum(np.ma.sum(Klw*c1/sum_c,0),0)                   # Eq. B4
    Klw_prime = Klw - Klw0                                          # Eq. B3
    B7a = np.ma.sum(c1/sum_c,1,keepdims=True)                       # need to keep this as [TAU,1,...]
    Klw_p_prime = np.ma.sum(Klw_prime*B7a,0)                        # Eq. B7
    Klw_t_prime = np.ma.sum(Klw_prime*np.ma.sum(c1/sum_c,0),1)      # Eq. B8   
    Klw_resid_prime = Klw_prime - np.expand_dims(Klw_p_prime,0) - np.expand_dims(Klw_t_prime,1)        # Eq. B9
    dRlw_true = np.ma.sum(np.ma.sum(Klw*dc,1),0)                    # LW total
    dRlw_prop = Klw0*sum_dc                                         # LW amount component
    dRlw_dctp = np.ma.sum(Klw_p_prime*np.ma.sum(dc_star,0),0)       # LW altitude component
    dRlw_dtau = np.ma.sum(Klw_t_prime*np.ma.sum(dc_star,1),0)       # LW optical depth component
    dRlw_resid = np.ma.sum(np.ma.sum(Klw_resid_prime*dc_star,1),0)  # LW residual
    dRlw_sum = dRlw_prop + dRlw_dctp + dRlw_dtau + dRlw_resid       # sum of LW components -- should equal LW total

    # SW components
    Ksw0 = np.ma.sum(np.ma.sum(Ksw*c1/sum_c,0),0)                   # Eq. B4
    Ksw_prime = Ksw - Ksw0                                          # Eq. B3
    B7a = np.ma.sum(c1/sum_c,1,keepdims=True)                       # need to keep this as [TAU,1,...]
    Ksw_p_prime = np.ma.sum(Ksw_prime*B7a,0)                        # Eq. B7
    Ksw_t_prime = np.ma.sum(Ksw_prime*np.ma.sum(c1/sum_c,0),1)      # Eq. B8  
    Ksw_resid_prime = Ksw_prime - np.expand_dims(Ksw_p_prime,0) - np.expand_dims(Ksw_t_prime,1)        # Eq. B9 
    dRsw_true = np.ma.sum(np.ma.sum(Ksw*dc,1),0)                    # SW total
    dRsw_prop = Ksw0*sum_dc                                         # SW amount component
    dRsw_dctp = np.ma.sum(Ksw_p_prime*np.ma.sum(dc_star,0),0)       # SW altitude component
    dRsw_dtau = np.ma.sum(Ksw_t_prime*np.ma.sum(dc_star,1),0)       # SW optical depth component
    dRsw_resid = np.ma.sum(np.ma.sum(Ksw_resid_prime*dc_star,1),0)  # SW residual
    dRsw_sum = dRsw_prop + dRsw_dctp + dRsw_dtau + dRsw_resid       # sum of SW components -- should equal SW total

    # Set SW fields to zero where the sun is down
    RR = Ksw0.mask
    dRsw_true = np.ma.where(RR,0,dRsw_true)
    dRsw_prop = np.ma.where(RR,0,dRsw_prop)
    dRsw_dctp = np.ma.where(RR,0,dRsw_dctp)
    dRsw_dtau = np.ma.where(RR,0,dRsw_dtau)
    dRsw_resid = np.ma.where(RR,0,dRsw_resid)

    # Move month axis back to the beginning
    output={}
    output['LWcld_tot'] = np.moveaxis(dRlw_true,-1,0)
    output['LWcld_amt'] = np.moveaxis(dRlw_prop,-1,0)
    output['LWcld_alt'] = np.moveaxis(dRlw_dctp,-1,0)
    output['LWcld_tau'] = np.moveaxis(dRlw_dtau,-1,0)
    output['LWcld_err'] = np.moveaxis(dRlw_resid,-1,0)
    output['SWcld_tot'] = np.moveaxis(dRsw_true,-1,0)
    output['SWcld_amt'] = np.moveaxis(dRsw_prop,-1,0)
    output['SWcld_alt'] = np.moveaxis(dRsw_dctp,-1,0)
    output['SWcld_tau'] = np.moveaxis(dRsw_dtau,-1,0)
    output['SWcld_err'] = np.moveaxis(dRsw_resid,-1,0)
    #output['dc_star'] = MV.array(np.moveaxis(dc_star,-1,0))
    #output['dc_prop'] = MV.array(np.moveaxis(dc_prop,-1,0))    
    
    return output    

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if __name__ == "__main__":

#    direc_kernel = '../CloudRadKernel_input/'
#    #direc_data = '/compyfs/qiny108/diag_feedback_E3SM_postdata/'
#    #direc_data = '/compyfs/qiny108/colla/diag_feedback_E3SM_postdata/'
#    direc_data = '/p/user_pub/climate_work/qin4/From_Compy/compyfs_dir/colla/diag_feedback_E3SM_postdata/'
#
#    case_stamps = [\
#    'v2test'
#    ]
#
#    for case_stamp in case_stamps:
#
#        fname1,_,_ = CL.get_lutable(case_stamp,'amip')
#        fname2,_,_ = CL.get_lutable(case_stamp,'amip4K')
#
#        outdir = './'
#        figdir = './'
#
#        yearS = 2
#        yearE = 3
#
#        CloudRadKernel_E3SM(direc_kernel,direc_data,case_stamp,yearS,yearE,fname1,fname2,outdir,figdir)

    # ------------------------------
    # Input requirements for MIP data
    # ------------------------------
    model = 'GFDL-CM4'  
    institution = 'NOAA-GFDL'
    variant = 'r1i1p1f1'

    ff = 'filenames_'+model+'_'+variant+'.pickle'
    filenames = read_pickle(ff)
    #print('filenames=',filenames)

    case_stamp = model+'_'+variant
    outdir = './'
    figdir = './'
    direc_kernel = '../CloudRadKernel_input/'

    tslice = slice('1980-01-01','1981-12-31')

    CloudRadKernel_MIP(direc_kernel,case_stamp,outdir,figdir,filenames,tslice)

