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
from PlotDefinedFunction import linearregression_nd

###########################################################################
# HELPFUL FUNCTIONS FOLLOW
###########################################################################
def weighted_temporal_mean(time, obs):
  """
  weight by days in each month
  """
  # Determine the month length
  month_length = time.dt.days_in_month

  # Calculate the weights
  wgts = month_length.groupby("time.year") / month_length.groupby("time.year").sum()

  # Make sure the weights in each year add up to 1
  np.testing.assert_allclose(wgts.groupby("time.year").sum(xr.ALL_DIMS), 1.0)

  # Setup our masking for nan values
  cond = obs.isnull()
  ones = xr.where(cond, 0.0, 1.0)

  # Calculate the numerator
  obs_sum = (obs * wgts).resample(time="AS").sum(dim="time")

  # Calculate the denominator
  ones_out = (ones * wgts).resample(time="AS").sum(dim="time")

  # Return the weighted average
  return obs_sum / ones_out


# ----------------------------------------------------
def area_averager(data_plot_xr):
    '''
    calculate weighted area mean
    input data is xarray DataArray
    '''
    weights = np.cos(np.deg2rad(data_plot_xr.lat))
    weights.name = "weights"
    # available in xarray version 0.15 and later
    data_weighted = data_plot_xr.weighted(weights)

    weighted_mean = data_weighted.mean(("lat", "lon"))

    return weighted_mean

# ------------------------------------------------------------------------
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
# MAIN ROUTINE FOLLOWS
##########################################################################

###########################################################################
# Part 1: Read in data, regrid, compute anomalies, map kernels to lat/lon
# This is identical to the first part of apply_cloud_kernels_v2.py
###########################################################################
def CloudRadKernel(direc_kernel,direc_data,case_stamp,yearS,yearE,fname1,fname2,outdir,figdir):

    if os.path.isfile(outdir+'global_cloud_feedback_'+case_stamp+'.nc'):
        print('CloudRadKernel is already there.')
        return

    yearS_4d = "{:04d}".format(yearS)
    yearE_4d = "{:04d}".format(yearE)
    nyears = yearE - yearS + 1

    direc_data1 = direc_data+'/'+fname1+'/'
    direc_data2 = direc_data+'/'+fname2+'/'

    exp1='FC5'
    exp2='FC5_4K'
    
    used_models = 'E3SM-1-0'
    
    yrS=yearS
    yrE=yearE
    monS=1
    monE=12
    
    yrS_4d='{:04d}'.format(yrS)
    yrE_4d='{:04d}'.format(yrE)
    monS_2d='{:02d}'.format(monS)
    monE_2d='{:02d}'.format(monE)
    
    # Load in the Zelinka et al 2012 kernels:
    # ------------------------------------------------------------------------------------
    f=xr.open_dataset(direc_kernel+'cloud_kernels2.nc')
    LWkernel=f['LWkernel']
    SWkernel=f['SWkernel']
    f.close()
    
    #LWkernel=MV.masked_where(np.isnan(LWkernel),LWkernel)
    #SWkernel=MV.masked_where(np.isnan(SWkernel),SWkernel)
    
    albcs=np.arange(0.0,1.5,0.5) # the clear-sky albedos over which the kernel is computed
    
    # LW kernel does not depend on albcs, just repeat the final dimension over longitudes:
    LWkernel_map=xr.DataArray(np.transpose(np.tile(np.tile(LWkernel[:,:,:,:,0],(1,nyears,1,1,1)),(144,1,1,1,1)),[1,2,3,4,0]))

    # Define the cloud kernel axis attributes
    lats=LWkernel.lat.values
    lons=np.arange(1.25,360,2.5)

    print('Start reading clisccp ...')
    # ------------------------------------------------------------------------------------
    # Load in clisccp from control and perturbed simulation

    if not os.path.isfile(direc_data1+'clisccp_'+exp1+'_'+yrS_4d+monS_2d+'-'+yrE_4d+monE_2d+'.nc'):
        svar_in = 'FISCCP1_COSP'
    else:
        svar_in = 'clisccp'
        
    f=xr.open_dataset(direc_data1+svar_in+'_'+exp1+'_'+yrS_4d+monS_2d+'-'+yrE_4d+monE_2d+'.nc')
    clisccp1=f[svar_in]#.transpose('time','cosp_tau','cosp_prs','lat','lon') # the old order is: (time, CTP, TAU, LAT, LON)
    f.close()
    f=xr.open_dataset(direc_data2+svar_in+'_'+exp2+'_'+yrS_4d+monS_2d+'-'+yrE_4d+monE_2d+'.nc')
    clisccp2=f[svar_in]#.transpose('time','cosp_tau','cosp_prs','lat','lon')
    f.close()

    # Make sure clisccp is in percent  
    sumclisccp1 = clisccp1.sum(axis=2).sum(axis=1)
    sumclisccp2 = clisccp2.sum(axis=2).sum(axis=1)
    if np.max(sumclisccp1) <= 1.:
        print('Changing clisccp in percent...')
        clisccp1 = clisccp1*100.        
    if np.max(sumclisccp2) <= 1.:
        clisccp2 = clisccp2*100.
    
    # Compute clisccp anomalies
    anomclisccp = xr.DataArray(clisccp2 - clisccp1, coords=clisccp1.coords)
    
    clisccp1_grd = clisccp1.interp(lat=lats,lon=lons).transpose('time','cosp_tau','cosp_prs','lat','lon')
    del clisccp1
    clisccp2_grd = clisccp2.interp(lat=lats,lon=lons).transpose('time','cosp_tau','cosp_prs','lat','lon')
    del clisccp2
    anomclisccp_grd = anomclisccp.interp(lat=lats,lon=lons).transpose('time','cosp_tau','cosp_prs','lat','lon')
    del anomclisccp

    print('Start reading albedo ...')
    # ------------------------------------------------------------------------------------
    # Compute clear-sky surface albedo
    f=xr.open_dataset(direc_data1+'rsuscs_'+exp1+'_'+yrS_4d+monS_2d+'-'+yrE_4d+monE_2d+'.nc')
    rsuscs1 = f['rsuscs']
    f.close()
    f=xr.open_dataset(direc_data1+'rsdscs_'+exp1+'_'+yrS_4d+monS_2d+'-'+yrE_4d+monE_2d+'.nc')
    rsdscs1 = f['rsdscs']
    f.close()

    albcs1=xr.DataArray(rsuscs1/rsdscs1, coords=rsuscs1.coords)
    albcs1=xr.where(albcs1>1.,1,albcs1) # where(condition, x, y) is x where condition is true, y otherwise
    albcs1=xr.where(albcs1<0.,0,albcs1)

    albcs1_grd = albcs1.interp(lat=lats,lon=lons)
    del albcs1, rsuscs1, rsdscs1 
 
    print('Start reading surface air temperature ...')
    # ------------------------------------------------------------------------------------
    # Load surface air temperature
    f=xr.open_dataset(direc_data1+'tas_'+exp1+'_'+yrS_4d+monS_2d+'-'+yrE_4d+monE_2d+'.nc',decode_times=False)
    tas1 = f['tas']
    f.close()
    f=xr.open_dataset(direc_data2+'tas_'+exp2+'_'+yrS_4d+monS_2d+'-'+yrE_4d+monE_2d+'.nc',decode_times=False)
    tas2 = f['tas']
    f.close()

    # Compute tas anomaly and global annual mean 
    anomtas = xr.DataArray(tas2 - tas1, coords=tas2.coords)

    # Define new time coordinate
    newtime = pd.date_range(start='1850-01-01', periods=tas1.shape[0], freq='MS')
    anomtas = anomtas.assign_coords({'time':("time",newtime)})

#    anomtas_ann = anomtas.groupby('time.year').mean('time')
    anomtas_ann = weighted_temporal_mean(anomtas.time,anomtas)

    anomtas_ann_gm = area_averager(anomtas_ann)

    del(tas1,tas2)
    
    # Compute global annual mean tas anomalies
    avgdtas = area_averager(anomtas_ann.mean(axis=0)) # (scalar)
    print('avgdtas = ',avgdtas.values)
   
    # =======================================================================================
    # Use control albcs to map SW kernel to appropriate longitudes
    # Jan 09, 2021: follow Mark's method -- use climatological control albedo to map SWkernel_map
    # rather than the 150-yr annual cycle control albedo
    albcs1_grd = albcs1_grd.assign_coords({'time':("time",newtime)})
    avgalbcs1 = albcs1_grd.groupby('time.month').mean('time')
    SWkernel_map_tmp = map_SWkern_to_lon(SWkernel,avgalbcs1)
    SWkernel_map = xr.DataArray(np.tile(SWkernel_map_tmp,(nyears,1,1,1,1)), coords=clisccp1_grd.coords)
    del SWkernel_map_tmp

    # The sun is down if every bin of the SW kernel is zero:
    sundown=np.ma.sum(np.ma.sum(SWkernel_map,axis=2),axis=1)  #12*nyears,90,144
    night=np.where(sundown==0)
    print("data processing is done. Please continue.")
    
    ###########################################################################
    # Part 2: Compute cloud feedbacks and their breakdown into components
    ###########################################################################         
    print('Start computing cloud feedbacks ...')

    # Define a python dictionary containing the sections of the histogram to consider
    # These are the same as in Zelinka et al, GRL, 2016
    sections = ['ALL','HI680','LO680']
    Psections=[slice(0,7),slice(2,7),slice(0,2)]

    sec_dic=dict(zip(sections,Psections))
    
    df_sw_all = pd.DataFrame()
    df_lw_all = pd.DataFrame()
    

    dic_out1 = {}
    #<qinyi 2021-02-25 #------------------
    # add output of monthly radiation anomalies caused by different cloud properties
    dic_out2 = {}

    for sec in sections:
        print ('Using '+sec+' CTP bins')
        choose=sec_dic[sec]
        LC = len(np.ones(100)[choose])
    
        # Preallocation of arrays:
        LWcld_tot=nanarray((12*nyears,90,144))
        LWcld_amt=nanarray((12*nyears,90,144))
        LWcld_alt=nanarray((12*nyears,90,144))
        LWcld_tau=nanarray((12*nyears,90,144))
        LWcld_err=nanarray((12*nyears,90,144))
        SWcld_tot=nanarray((12*nyears,90,144))
        SWcld_amt=nanarray((12*nyears,90,144))
        SWcld_alt=nanarray((12*nyears,90,144))
        SWcld_tau=nanarray((12*nyears,90,144))
        SWcld_err=nanarray((12*nyears,90,144))
        dc_star=nanarray((12*nyears,7,LC,90,144))
        dc_prop=nanarray((12*nyears,7,LC,90,144))
    
        for mm in range(12*nyears):
            dcld_dT = anomclisccp_grd[mm,:,choose,:]
            c1 = clisccp1_grd[mm,:,choose,:]
            c2 = c1 + dcld_dT
            Klw = LWkernel_map[mm,:,choose,:]
            Ksw = SWkernel_map[mm,:,choose,:]

            # The following performs the amount/altitude/optical depth decomposition of
            # Zelinka et al., J Climate (2012b), as modified in Zelinka et al., J. Climate (2013)
            (LWcld_tot[mm,:],LWcld_amt[mm,:],LWcld_alt[mm,:],LWcld_tau[mm,:],LWcld_err[mm,:],SWcld_tot[mm,:],SWcld_amt[mm,:],SWcld_alt[mm,:],SWcld_tau[mm,:],SWcld_err[mm,:],dc_star[mm,:],dc_prop[mm,:]) = KT_decomposition_4D(c1,c2,Klw,Ksw)

        # Set the SW cloud feedbacks to zero in the polar night
        # Do this since they may come out of previous calcs as undefined, but should be zero:
        SWcld_tot[night]=0
        SWcld_amt[night]=0
        SWcld_alt[night]=0
        SWcld_tau[night]=0
        SWcld_err[night]=0

        AX3 = albcs1_grd.coords #[coord_time,coord_lats,coord_lons]
        AX= albcs1_grd[0,:,:].coords #[coord_lats,coord_lons]

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
            DATA_am = weighted_temporal_mean(DATA_anom.time,DATA_anom)

            if 'coupled' in case_stamp:
                slope,intercept = linearregression_nd(DATA_am, x = np.reshape(anomtas_ann_gm.values, (anomtas_ann_gm.shape[0],1,1)))
                #print('slope=',slope.shape, np.min(slope), np.max(slope), 'intercept=',intercept.shape, np.min(intercept), np.max(intercept))
            else:
                slope = DATA_anom.groupby('time.month').mean('time').mean(axis=0)/anomtas_ann_gm.mean()

            dic_out1[str(sec)+'_'+str(name)] = slope 
            dic_out1[str(sec)+'_'+str(name)].attrs['long_name'] = str(sec)+'_'+str(name)

            avgDATA = area_averager(slope).values
            if 'SW' in name:
                df_sw_tmp = pd.DataFrame([[sec,'CTP bins',name,str(np.round(avgDATA,3))]],columns=['type','bin','decomp',used_models])
                df_sw_all = pd.concat([df_sw_all,df_sw_tmp])
            elif 'LW' in name:
                df_lw_tmp = pd.DataFrame([[sec,'CTP bins',name,str(np.round(avgDATA,3))]],columns=['type','bin','decomp',used_models])
                df_lw_all = pd.concat([df_lw_all,df_lw_tmp])

   
    print(df_lw_all.head())
    print(df_sw_all.head())
    
    df_lw_all.to_csv(outdir+'decomp_global_mean_lw_'+case_stamp+'_xr.csv')
    df_sw_all.to_csv(outdir+'decomp_global_mean_sw_'+case_stamp+'_xr.csv')

    
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if __name__ == "__main__":

    direc_kernel = '../CloudRadKernel_input/'
    direc_data = '/compyfs/qiny108/diag_feedback_E3SM_postdata/'
    #direc_data = '/compyfs/qiny108/colla/diag_feedback_E3SM_postdata/'

    case_stamps = [\
    'v2test_coupled'
    ]

    for case_stamp in case_stamps:

        fname1,_,_ = CL.get_lutable(case_stamp,'amip')
        fname2,_,_ = CL.get_lutable(case_stamp,'amip4K')

        outdir = './'
        figdir = './'

        exp1 = 'FC5'
        exp2 = 'FC5_4K'

        yearS = 2
        yearE = 6

        CloudRadKernel(direc_kernel,direc_data,case_stamp,yearS,yearE,fname1,fname2,outdir,figdir)
