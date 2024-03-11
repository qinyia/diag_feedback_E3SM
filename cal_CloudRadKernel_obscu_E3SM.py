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
 
#IMPORT STUFF:
#=====================
import cdms2 as cdms
import cdutil
import MV2 as MV
import numpy as np
import pylab as pl
import matplotlib as mpl
mpl.use('Agg')
import sys
import pandas as pd
import genutil

## qinyi 
import cartopy.crs as ccrs
import cartopy
import matplotlib.pyplot as plt
import os
 
###########################################################################
# HELPFUL FUNCTIONS FOLLOW
###########################################################################
def compute_fbk(ctl, fut, DT):
    DR = fut - ctl
    fbk = DR / DT
    baseline = ctl
    return fbk, baseline

###########################################################################
def obscuration_feedback_terms_general(L_R_bar0,dobsc_fbk,dunobsc_fbk,dobsc_cov_fbk,Klw,Ksw):
    """
    Estimate unobscured low cloud feedback, 
    the low cloud feedback arising solely from changes in obscuration by upper-level clouds,
    and the covariance term
    
    This function takes in a (month,tau,CTP,lat,lon) matrix
   
    Klw and Ksw contain just the low bins
    
    the following terms are generated in obscuration_terms():
    dobsc = L_R_bar0 * F_prime
    dunobsc = L_R_prime * F_bar
    dobsc_cov = (L_R_prime * F_prime) - climo(L_R_prime * F_prime)
    """
    
    Klw_low = Klw
    Ksw_low = Ksw
    L_R_bar0 = 100*L_R_bar0
    dobsc_fbk = 100*dobsc_fbk
    dunobsc_fbk = 100*dunobsc_fbk
    dobsc_cov_fbk = 100*dobsc_cov_fbk
    
    LWdobsc_fbk = MV.sum(MV.sum(Klw_low * dobsc_fbk,1),1)
    LWdunobsc_fbk = MV.sum(MV.sum(Klw_low * dunobsc_fbk,1),1)
    LWdobsc_cov_fbk = MV.sum(MV.sum(Klw_low * dobsc_cov_fbk,1),1)    
    
    SWdobsc_fbk = MV.sum(MV.sum(Ksw_low * dobsc_fbk,1),1)
    SWdunobsc_fbk = MV.sum(MV.sum(Ksw_low * dunobsc_fbk,1),1)
    SWdobsc_cov_fbk = MV.sum(MV.sum(Ksw_low * dobsc_cov_fbk,1),1)    
    
    ###########################################################################
    # Further break down the true and apparent low cloud-induced radiation anomalies into components
    ###########################################################################
    # No need to break down dobsc_fbk, as that is purely an amount component.
        
    # Break down dunobsc_fbk:
    C_ctl = L_R_bar0
    dC = dunobsc_fbk
    C_fut = C_ctl + dC
    
    obsc_fbk_output = KT_decomposition_4D(C_ctl,C_fut,Klw_low,Ksw_low)        
   
    obsc_fbk_output['LWdobsc_fbk'] = LWdobsc_fbk
    obsc_fbk_output['LWdunobsc_fbk'] = LWdunobsc_fbk
    obsc_fbk_output['LWdobsc_cov_fbk'] = LWdobsc_cov_fbk
    obsc_fbk_output['SWdobsc_fbk'] = SWdobsc_fbk
    obsc_fbk_output['SWdunobsc_fbk'] = SWdunobsc_fbk
    obsc_fbk_output['SWdobsc_cov_fbk'] = SWdobsc_cov_fbk
    
    return obsc_fbk_output
        
###########################################################################
def obscuration_terms3(c1,c2):
    """
    USE THIS VERSION FOR DIFFERENCES OF 2 CLIMATOLOGIES (E.G. AMIP4K, 2xCO2 SLAB RUNS)
    
    Compute the components required for the obscuration-affected low cloud feedback
    These are the terms shown in Eq 4 of Scott et al (2020) DOI: 10.1175/JCLI-D-19-1028.1
    L_prime = dunobsc + dobsc + dobsc_cov, where
    dunobsc = L_R_prime * F_bar     (delta unobscured low clouds, i.e., true low cloud feedback)
    dobsc = L_R_bar * F_prime       (delta obscuration by upper level clouds)
    dobsc_cov = (L_R_prime * F_prime) - climo(L_R_prime * F_prime)  (covariance term)
    """
    # c is [mo,tau,ctp,lat,lon]
    # c is in percent
    
    AX = c2.getAxisList()
    
    c1=MV.masked_where(c2.mask,c1)
    c2=MV.masked_where(c1.mask,c2)
    
    # SPLICE c1 and c2:
    # MAKE SURE c1 and c2 are the same size!!!
    if c1.shape != c2.shape:
        raise RuntimeError('c1 and c2 are NOT the same size!!!')
        
    c12=np.ma.append(c1,c2,axis=0)    
    
    midpt=len(c1)
           
    U12 = MV.sum(MV.sum(c12[:,:,2:,:],1),1)/100.
    
    L12 = c12[:,:,:2,:]/100.  
    
    F12 = 1. - U12
    F12=MV.masked_less(F12,0)

    F12b = MV.array(np.expand_dims(np.expand_dims(F12,axis=1),axis=1))
    F12b=MV.masked_where(L12[:,:1,:1,:].mask,F12b)
    
    L_R12 = L12/F12b
    sum_L_R12 = MV.sum(MV.sum(L_R12,1),1)
    sum_L_R12b = MV.array(np.expand_dims(np.expand_dims(sum_L_R12,axis=1),axis=1))
    sum_L_R12c = np.broadcast_to(sum_L_R12b,L_R12.shape)
    this = MV.masked_outside(sum_L_R12c,0,1)
    L_R12 = MV.masked_where(this.mask,L_R12)
    
    L_R12 = MV.masked_where(sum_L_R12c>1,L_R12)
        
    L_R_prime,L_R_bar = monthly_anomalies(L_R12)
    F_prime,F_bar = monthly_anomalies(F12b)    
    L_prime,L_bar = monthly_anomalies(L12)
    
    # Cannot have negative cloud fractions:
    L_R_bar[L_R_bar<0]=0
    F_bar[F_bar<0]=0    
    
    rep_L_bar = tile_uneven(L_bar,L12)    
    rep_L_R_bar = tile_uneven(L_R_bar,L_R12)            
    rep_F_bar = tile_uneven(F_bar,F12b)
    
    # Cannot have negative cloud fractions:
    L_R_bar[L_R_bar<0]=0
    F_bar[F_bar<0]=0

    dobsc = rep_L_R_bar * F_prime
    dunobsc = L_R_prime * rep_F_bar
    prime_prime = (L_R_prime * F_prime)

    dobsc_cov,climo_prime_prime = monthly_anomalies(prime_prime)   
    
    # Re-scale these anomalies by 2, since we have computed all anomalies w.r.t. 
    # the ctl+pert average rather than w.r.t. the ctl average
    dobsc*=2
    dunobsc*=2
    dobsc_cov*=2
    
    return(rep_L_R_bar[midpt:],dobsc[midpt:],dunobsc[midpt:],dobsc_cov[midpt:])    

###########################################################################

def do_obscuration_calcs(CTL, FUT, Klw, Ksw, DT):
    (L_R_bar, dobsc, dunobsc, dobsc_cov) = obscuration_terms3(CTL, FUT)

    # Get unobscured low-cloud feedbacks and those caused by change in obscuration
    ZEROS = np.zeros(L_R_bar.shape)
    dummy, L_R_bar_base = compute_fbk(L_R_bar, L_R_bar, DT)
    dobsc_fbk, dummy = compute_fbk(ZEROS, dobsc, DT)
    dunobsc_fbk, dummy = compute_fbk(ZEROS, dunobsc, DT)
    dobsc_cov_fbk, dummy = compute_fbk(ZEROS, dobsc_cov, DT)
    obsc_output = obscuration_feedback_terms_general(
        L_R_bar_base, dobsc_fbk, dunobsc_fbk, dobsc_cov_fbk, Klw, Ksw
    )

    return obsc_output
 
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
            anndata0[cnt] = MV.average(data[i:i+12],0)
    B=cnt+1
    anndata = anndata0[:B]
    if type(anndata)!=float:
        anndata.setAxisList(data[:B*12:12].getAxisList())

    return anndata

###########################################################################
def add_cyclic(data):
    # Add Cyclic point around 360 degrees longitude:
    lons=data.getLongitude()[:]
    dx=np.gradient(lons)[-1]
    data2 = data(longitude=(0, dx+np.max(lons)), squeeze=True)    
    return data2

###########################################################################
def nanarray(vector):
    # this generates a masked array with the size given by vector
    # example: vector = (90,144,28)
    # similar to this=NaN*ones(x,y,z) in matlab
    this=MV.zeros(vector)
    this=MV.masked_where(this==0,this)

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
            try:
                alon2=MV.where(alon.mask,np.nan,alon)   
            except:
                alon2=alon
            if np.ma.count(alon2)>1: # at least 1 unmasked value
                if len(np.where(Ksw[MM,:,:,LA,:]>0))==0:
#                if len(np.where(Ksw[MM,:,:,LA,:]>0))==0:
                    SWkernel_map[M,:,:,LA,:] = 0
                else:
                    f = interp1d(albcs,Ksw[MM,:,:,LA,:],axis=2)
                    ynew = f(alon2.data)
                    ynew=MV.masked_where(alon2.mask,ynew)
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
    c1 = MV.transpose(c1,(1,0,2,3)) # control cloud fraction histogram
    c2 = MV.transpose(c2,(1,0,2,3)) # perturbed cloud fraction histogram
    Klw = MV.transpose(Klw,(1,0,2,3)) # LW Kernel histogram
    Ksw = MV.transpose(Ksw,(1,0,2,3)) # SW Kernel histogram

    P=c1.shape[0]
    T=c1.shape[1]

    c=c1
    sum_c=np.tile(MV.sum(MV.sum(c,1),0),(P,T,1,1))                                  # Eq. B2
    dc = c2-c1 
    sum_dc=np.tile(MV.sum(MV.sum(dc,1),0),(P,T,1,1))
    dc_prop = c*(sum_dc/sum_c)
    dc_star = dc - dc_prop                                                          # Eq. B1

    # LW components
    Klw0 = np.tile(MV.sum(MV.sum(Klw*c/sum_c,1),0),(P,T,1,1))                       # Eq. B4
    Klw_prime = Klw - Klw0                                                          # Eq. B3
    this=MV.sum(Klw_prime*np.tile(MV.sum(c/sum_c,0),(P,1,1,1)),1)                   # Eq. B7a
    Klw_p_prime=np.tile(np.tile(this,(1,1,1,1)),(T,1,1,1))(order=[1,0,2,3])         # Eq. B7b 
    that=np.tile(np.tile(MV.sum(c/sum_c,1),(1,1,1,1)),(T,1,1,1))(order=[1,0,2,3])   # Eq. B8a
    Klw_t_prime = np.tile(MV.sum(Klw_prime*that,0),(P,1,1,1))                       # Eq. B8b
    Klw_resid_prime = Klw_prime - Klw_p_prime - Klw_t_prime                         # Eq. B9
    dRlw_true = MV.sum(MV.sum(Klw*dc,1),0)                                          # LW total
    dRlw_prop = Klw0[0,0,:,:]*sum_dc[0,0,:,:]                                       # LW amount component
    dRlw_dctp = MV.sum(MV.sum(Klw_p_prime*dc_star,1),0)                             # LW altitude component
    dRlw_dtau = MV.sum(MV.sum(Klw_t_prime*dc_star,1),0)                             # LW optical depth component
    dRlw_resid = MV.sum(MV.sum(Klw_resid_prime*dc_star,1),0)                        # LW residual
    dRlw_sum = dRlw_prop + dRlw_dctp + dRlw_dtau + dRlw_resid                       # sum of LW components -- should equal LW total

    # SW components
    Ksw0 = np.tile(MV.sum(MV.sum(Ksw*c/sum_c,1),0),(P,T,1,1))                       # Eq. B4
    Ksw_prime = Ksw - Ksw0                                                          # Eq. B3
    this=MV.sum(Ksw_prime*np.tile(MV.sum(c/sum_c,0),(P,1,1,1)),1)                   # Eq. B7a 
    Ksw_p_prime=np.tile(np.tile(this,(1,1,1,1)),(T,1,1,1))(order=[1,0,2,3])         # Eq. B7b  
    that=np.tile(np.tile(MV.sum(c/sum_c,1),(1,1,1,1)),(T,1,1,1))(order=[1,0,2,3])   # Eq. B8a
    Ksw_t_prime = np.tile(MV.sum(Ksw_prime*that,0),(P,1,1,1))                       # Eq. B8b
    Ksw_resid_prime = Ksw_prime - Ksw_p_prime - Ksw_t_prime                         # Eq. B9
    dRsw_true = MV.sum(MV.sum(Ksw*dc,1),0)                                          # SW total
    dRsw_prop = Ksw0[0,0,:,:]*sum_dc[0,0,:,:]                                       # SW amount component
    dRsw_dctp = MV.sum(MV.sum(Ksw_p_prime*dc_star,1),0)                             # SW altitude component
    dRsw_dtau = MV.sum(MV.sum(Ksw_t_prime*dc_star,1),0)                             # SW optical depth component
    dRsw_resid = MV.sum(MV.sum(Ksw_resid_prime*dc_star,1),0)                        # SW residual
    dRsw_sum = dRsw_prop + dRsw_dctp + dRsw_dtau + dRsw_resid                       # sum of SW components -- should equal SW total

    dc_star = MV.transpose(dc_star,(1,0,2,3)) 
    dc_prop = MV.transpose(dc_prop,(1,0,2,3)) 

    return (dRlw_true,dRlw_prop,dRlw_dctp,dRlw_dtau,dRlw_resid,dRsw_true,dRsw_prop,dRsw_dctp,dRsw_dtau,dRsw_resid,dc_star,dc_prop)
 
###########################################################################
# MAIN ROUTINE FOLLOWS
##########################################################################

###########################################################################
# Part 1: Read in data, regrid, compute anomalies, map kernels to lat/lon
# This is identical to the first part of apply_cloud_kernels_v2.py
###########################################################################
def CloudRadKernel(direc_kernel,direc_data,case_stamp,yearS,yearE,fname1,fname2,outdir,figdir):

    if os.path.isfile(outdir+'global_cloud_feedback_'+case_stamp+'.nc') and os.path.isfile(outdir+'decomp_global_mean_lw_'+case_stamp+'.csv'):
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
    f=cdms.open(direc_kernel+'cloud_kernels2.nc')
    LWkernel=f('LWkernel')
    SWkernel=f('SWkernel')
    f.close()
    
    LWkernel=MV.masked_where(np.isnan(LWkernel),LWkernel)
    SWkernel=MV.masked_where(np.isnan(SWkernel),SWkernel)
    
    albcs=np.arange(0.0,1.5,0.5) # the clear-sky albedos over which the kernel is computed
    
    # LW kernel does not depend on albcs, just repeat the final dimension over longitudes:
    LWkernel_map=np.tile(np.tile(LWkernel[:,:,:,:,0],(1,nyears,1,1,1)),(144,1,1,1,1))(order=[1,2,3,4,0])

    # Define the cloud kernel axis attributes
    lats=LWkernel.getLatitude()[:]
    lons=np.arange(1.25,360,2.5)
    grid = cdms.createGenericGrid(lats,lons)

    # =======================================================================================
    print('Start reading clisccp ...')
    # Load in clisccp from control and perturbed simulation

    #<qinyi 2021-08-12 #------------------
    if not os.path.isfile(direc_data1+'clisccp_'+exp1+'_'+yrS_4d+monS_2d+'-'+yrE_4d+monE_2d+'.nc'):
        svar_in = 'FISCCP1_COSP'
    else:
        svar_in = 'clisccp'
    #<qinyi 2021-08-12 #------------------
        
    f=cdms.open(direc_data1+svar_in+'_'+exp1+'_'+yrS_4d+monS_2d+'-'+yrE_4d+monE_2d+'.nc','r')
    clisccp1=f(svar_in,order="02134") # the old order is: (time, CTP, TAU, LAT, LON)
    f.close()
    f=cdms.open(direc_data2+svar_in+'_'+exp2+'_'+yrS_4d+monS_2d+'-'+yrE_4d+monE_2d+'.nc','r')
    clisccp2=f(svar_in,order="02134")
    f.close()

    # Make sure clisccp is in percent  
    sumclisccp1=MV.sum(MV.sum(clisccp1,2),1)
    sumclisccp2=MV.sum(MV.sum(clisccp2,2),1)   
    if np.max(sumclisccp1) <= 1.:
        clisccp1 = clisccp1*100.        
    if np.max(sumclisccp2) <= 1.:
        clisccp2 = clisccp2*100.
    
    # Compute clisccp anomalies
    anomclisccp = clisccp2 - clisccp1
    
    clisccp1 = add_cyclic(clisccp1)
    clisccp1_grd = clisccp1.regrid(grid,regridTool="esmf",regridMethod = "linear")
    del clisccp1
    clisccp2 = add_cyclic(clisccp2)
    clisccp2_grd = clisccp2.regrid(grid,regridTool="esmf",regridMethod = "linear")
    del clisccp2
    anomclisccp = add_cyclic(anomclisccp)
    anomclisccp_grd = anomclisccp.regrid(grid,regridTool="esmf",regridMethod = "linear")
    del anomclisccp

    # =======================================================================================
    print('Start reading albedo ...')

    # Compute clear-sky surface albedo
    f=cdms.open(direc_data1+'rsuscs_'+exp1+'_'+yrS_4d+monS_2d+'-'+yrE_4d+monE_2d+'.nc','r')
    rsuscs1 = f('rsuscs')
    f.close()
    f=cdms.open(direc_data1+'rsdscs_'+exp1+'_'+yrS_4d+monS_2d+'-'+yrE_4d+monE_2d+'.nc','r')
    rsdscs1 = f('rsdscs')
    f.close()

    albcs1=rsuscs1/rsdscs1
    albcs1.setAxisList(rsuscs1.getAxisList())
    albcs1=MV.where(albcs1>1.,1,albcs1) # where(condition, x, y) is x where condition is true, y otherwise
    albcs1=MV.where(albcs1<0.,0,albcs1)

    albcs1 = add_cyclic(albcs1)
    albcs1_grd = albcs1.regrid(grid,regridTool="esmf",regridMethod = "linear")
    del albcs1, rsuscs1, rsdscs1 
 
    # =======================================================================================
    print('Start reading surface air temperature ...')

    # Load surface air temperature
    f=cdms.open(direc_data1+'tas_'+exp1+'_'+yrS_4d+monS_2d+'-'+yrE_4d+monE_2d+'.nc','r')
    tas1 = f('tas')
    f.close()
    f=cdms.open(direc_data2+'tas_'+exp2+'_'+yrS_4d+monS_2d+'-'+yrE_4d+monE_2d+'.nc','r')
    tas2 = f('tas')

    time = f['time']
    calendar = time.getCalendar()
    f.close()

    # Compute tas anomaly and global annual mean 
    anomtas = tas2 - tas1
    anomtas = cdms.asVariable(anomtas)
    anomtas.setAxisList(tas2.getAxisList())
    print('anomtas shape is',anomtas.shape)

    cdutil.setTimeBoundsMonthly(anomtas,1)
    anomtas_ann = cdutil.YEAR(anomtas)
    print('anomtas_ann shape is',anomtas_ann.shape)
    anomtas_ann_gm = cdutil.averager(anomtas_ann,axis='xy',weights='weighted')
    print('anomtas_ann_gm shape is ',anomtas_ann_gm.shape)

    del(tas1,tas2)
    
    # Compute global annual mean tas anomalies
    avgdtas = cdutil.averager(MV.average(anomtas_ann,axis=0), axis='xy', weights='weighted') # (scalar)
    print('avgdtas = ',avgdtas)
    print('MV.average(anomtas_ann_gm)=',MV.average(anomtas_ann_gm))
   
    # =======================================================================================
    # Use control albcs to map SW kernel to appropriate longitudes
    # Jan 09, 2021: follow Mark's method -- use climatological control albedo to map SWkernel_map
    # rather than the 150-yr annual cycle control albedo
    cdutil.setTimeBoundsMonthly(albcs1_grd,1)
    avgalbcs1 = cdutil.ANNUALCYCLE.climatology(albcs1_grd)
    SWkernel_map_tmp = map_SWkern_to_lon(SWkernel,avgalbcs1)
    SWkernel_map = np.tile(SWkernel_map_tmp,(nyears,1,1,1,1))
    SWkernel_map.setAxisList(clisccp1_grd.getAxisList())
    print('SWkernel_map.shape=',SWkernel_map.shape)
    del SWkernel_map_tmp

    # The sun is down if every bin of the SW kernel is zero:
    sundown=MV.sum(MV.sum(SWkernel_map,axis=2),axis=1)  #12*nyears,90,144
    night=np.where(sundown==0)
    print('sundown.shape=',sundown.shape)

    print("data processing is done. Please continue.")
    
    ###########################################################################
    # Part 2: Compute cloud feedbacks and their breakdown into components
    ###########################################################################         
    print('Start computing cloud feedbacks ...')

    # Define a python dictionary containing the sections of the histogram to consider
    # These are the same as in Zelinka et al, GRL, 2016
    sections = ['ALL','HI680','LO680']
    Psections=[slice(0,7),slice(2,7),slice(0,2)]

#    sections = ['ALL','HI680','LO680','HI560','LO560']
#    Psections=[slice(0,7),slice(2,7),slice(0,2),slice(3,7),slice(0,3)]

    sec_dic=dict(zip(sections,Psections))
    
    df_sw_all = pd.DataFrame()
    df_lw_all = pd.DataFrame()
    df_net_all = pd.DataFrame()
    

    value = 0
    cdms.setNetcdfShuffleFlag(value) ## where value is either 0 or 1
    cdms.setNetcdfDeflateFlag(value) ## where value is either 0 or 1
    cdms.setNetcdfDeflateLevelFlag(value) ## where value is a integer between 0 and 9 included

    out1 = cdms.open(outdir+'global_cloud_feedback_'+case_stamp+'.nc','w')

    #<qinyi 2021-02-25 #------------------
    # add output of monthly radiation anomalies caused by different cloud properties
    out2 = cdms.open(outdir+'global_cloud_anomaly_'+case_stamp+'.nc','w')

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

        # July 7, 2020 save variables into dictionary
        names = ['LWcld_tot','LWcld_amt','LWcld_alt','LWcld_tau','LWcld_err',\
        'SWcld_tot','SWcld_amt','SWcld_alt','SWcld_tau','SWcld_err',\
        'NETcld_tot','NETcld_amt','NETcld_alt','NETcld_tau','NETcld_err']
        variables = [LWcld_tot,LWcld_amt,LWcld_alt,LWcld_tau,LWcld_err,\
        SWcld_tot,SWcld_amt,SWcld_alt,SWcld_tau,SWcld_err,\
        SWcld_tot+LWcld_tot,SWcld_amt+LWcld_amt,SWcld_alt+LWcld_alt,SWcld_tau+LWcld_tau,SWcld_err+LWcld_err,\
        ]

        dic_all = {}
        for n,name in enumerate(names):
            dic_all[name] = variables[n]

        AX3 = albcs1_grd.getAxisList() #[coord_time,coord_lats,coord_lons]
        AX= albcs1_grd[0,:,:].getAxisList() #[coord_lats,coord_lons]

        # ==============================================
        # Plot Maps
        lons=albcs1_grd.getLongitude()[:]
        lats=albcs1_grd.getLatitude()[:]
        LON, LAT = np.meshgrid(lons,lats)
    
        # SW
        fig=plt.figure(figsize=(18,12)) # this creates and increases the figure size
        plt.suptitle(case_stamp+': '+sec+' CTP bins',fontsize=16,y=0.95)
        bounds = np.arange(-4,4.5,0.5)
        cmap = pl.cm.RdBu_r
        bounds2 = np.append(np.append(-500,bounds),500) # This is only needed for norm if colorbar is extended
        norm = mpl.colors.BoundaryNorm(bounds2, cmap.N) # make sure the colors vary linearly even if the desired color boundaries are at varied intervals
        names = ['SWcld_tot','SWcld_amt','SWcld_alt','SWcld_tau','SWcld_err']
        for n,name in enumerate(names):
            ax1 = fig.add_subplot(3,2,n+1,projection=ccrs.Robinson(central_longitude=180.))
            DATA_anom = dic_all[name]
            DATA_anom.setAxisList(AX3)

            cdutil.setTimeBoundsMonthly(DATA_anom,1)
            DATA_am = cdutil.YEAR(DATA_anom)

            if 'coupled' in case_stamp:
                print('DATA_am.shape=',DATA_am.shape)
                print('anomtas_ann_gm.shape=',anomtas_ann_gm.shape)
                slope, intercept = genutil.statistics.linearregression(DATA_am,x=anomtas_ann_gm)
            else:
                slope = MV.average(cdutil.ANNUALCYCLE.climatology(DATA_anom),axis=0)/MV.average(anomtas_ann_gm)

            DATA = slope
            DATA.id = str(sec)+"_"+str(name)
            DATA.long_name = str(sec)+"_"+str(name)
            out1.write(DATA)

            #im1 = ax1.contourf(LON,LAT,DATA,bounds,transform=ccrs.PlateCarree(),cmap=cmap,norm=norm,extend='both',corner_mask = False)
            im1 = ax1.pcolormesh(LON,LAT,np.round(DATA,3),vmin=min(bounds),vmax=max(bounds),transform=ccrs.PlateCarree(),cmap=cmap,norm=norm)

            ax1.coastlines()
            ax1.set_global()
            DATA.setAxisList(AX)
            avgDATA = cdutil.averager(DATA, axis='xy', weights='weighted')
            plt.title(name+' ['+str(np.round(avgDATA,3))+']',fontsize=14)
            cb = plt.colorbar(im1,orientation='vertical',drawedges=True,ticks=bounds)
            cb.set_label('W/m$^2$/K')
            
            df_sw_tmp = pd.DataFrame([[sec,'CTP bins',name,str(np.round(avgDATA,3))]],columns=['type','bin','decomp',used_models])
            df_sw_all = pd.concat([df_sw_all,df_sw_tmp])
        
        plt.savefig(figdir+'SW_'+sec+'_cld_fbk_maps-'+used_models+'-'+case_stamp+'.png', bbox_inches='tight')
 
        # LW
        fig=plt.figure(figsize=(18,12)) # this creates and increases the figure size
        plt.suptitle(case_stamp+': '+sec+' CTP bins',fontsize=16,y=0.95)
        bounds = np.arange(-4,4.5,0.5)
        cmap = plt.cm.RdBu_r
        bounds2 = np.append(np.append(-500,bounds),500) # This is only needed for norm if colorbar is extended
        norm = mpl.colors.BoundaryNorm(bounds2, cmap.N) # make sure the colors vary linearly even if the desired color boundaries are at varied intervals
        names = ['LWcld_tot','LWcld_amt','LWcld_alt','LWcld_tau','LWcld_err']
        for n,name in enumerate(names):
            print('varname=',name)
            ax1 = fig.add_subplot(3,2,n+1,projection=ccrs.Robinson(central_longitude=180.))
            DATA_anom = dic_all[name]
            DATA_anom.setAxisList(AX3)

            cdutil.setTimeBoundsMonthly(DATA_anom,1)
            DATA_am = cdutil.YEAR(DATA_anom)

            if 'coupled' in case_stamp:
                print('DATA_am.shape=',DATA_am.shape)
                print('anomtas_ann_gm.shape=',anomtas_ann_gm.shape)
                slope, intercept = genutil.statistics.linearregression(DATA_am,x=anomtas_ann_gm)
            else:
                slope = MV.average(cdutil.ANNUALCYCLE.climatology(DATA_anom),axis=0)/MV.average(anomtas_ann_gm)

            DATA = slope
            DATA.id = str(sec)+"_"+str(name)
            DATA.long_name = str(sec)+"_"+str(name)
            out1.write(DATA)

            #im1 = ax1.contourf(LON,LAT,np.round(DATA,3),bounds,transform=ccrs.PlateCarree(),cmap=cmap,norm=norm,extend='both',corner_mask = False)
            im1 = ax1.pcolormesh(LON,LAT,np.round(DATA,3),vmin=min(bounds),vmax=max(bounds),transform=ccrs.PlateCarree(),cmap=cmap,norm=norm)


            ax1.coastlines()
            ax1.set_global()
            DATA.setAxisList(AX)
            avgDATA = cdutil.averager(DATA, axis='xy', weights='weighted')
            pl.title(name+' ['+str(np.round(avgDATA,3))+']',fontsize=14)
            cb = plt.colorbar(im1,orientation='vertical',drawedges=True,ticks=bounds)
            cb.set_label('W/m$^2$/K')
            
            df_lw_tmp = pd.DataFrame([[sec,'CTP bins',name,str(np.round(avgDATA,3))]],columns=['type','bin','decomp',used_models])
            df_lw_all = pd.concat([df_lw_all,df_lw_tmp])
              
        plt.savefig(figdir+'LW_'+sec+'_cld_fbk_maps-'+used_models+'-'+case_stamp+'.png', bbox_inches='tight')
    
   
        # NET
        fig=plt.figure(figsize=(18,12)) # this creates and increases the figure size
        plt.suptitle(case_stamp+': '+sec+' CTP bins',fontsize=16,y=0.95)
        bounds = np.arange(-4,4.5,0.5)
        cmap = plt.cm.RdBu_r
        bounds2 = np.append(np.append(-500,bounds),500) # This is only needed for norm if colorbar is extended
        norm = mpl.colors.BoundaryNorm(bounds2, cmap.N) # make sure the colors vary linearly even if the desired color boundaries are at varied intervals
        names = ['NETcld_tot','NETcld_amt','NETcld_alt','NETcld_tau','NETcld_err']
        for n,name in enumerate(names):
            print('varname=',name)
            ax1 = fig.add_subplot(3,2,n+1,projection=ccrs.Robinson(central_longitude=180.))
            DATA_anom = dic_all[name]
            DATA_anom.setAxisList(AX3)

            cdutil.setTimeBoundsMonthly(DATA_anom,1)
            DATA_am = cdutil.YEAR(DATA_anom)

            if 'coupled' in case_stamp:
                print('DATA_am.shape=',DATA_am.shape)
                print('anomtas_ann_gm.shape=',anomtas_ann_gm.shape)
                slope, intercept = genutil.statistics.linearregression(DATA_am,x=anomtas_ann_gm)
            else:
                slope = MV.average(cdutil.ANNUALCYCLE.climatology(DATA_anom),axis=0)/MV.average(anomtas_ann_gm)

            DATA = slope
            DATA.id = str(sec)+"_"+str(name)
            DATA.long_name = str(sec)+"_"+str(name)
            out1.write(DATA)

            #im1 = ax1.contourf(LON,LAT,np.round(DATA,3),bounds,transform=ccrs.PlateCarree(),cmap=cmap,norm=norm,extend='both',corner_mask = False)
            im1 = ax1.pcolormesh(LON,LAT,np.round(DATA,3),vmin=min(bounds),vmax=max(bounds),transform=ccrs.PlateCarree(),cmap=cmap,norm=norm)


            ax1.coastlines()
            ax1.set_global()
            DATA.setAxisList(AX)
            avgDATA = cdutil.averager(DATA, axis='xy', weights='weighted')
            pl.title(name+' ['+str(np.round(avgDATA,3))+']',fontsize=14)
            cb = plt.colorbar(im1,orientation='vertical',drawedges=True,ticks=bounds)
            cb.set_label('W/m$^2$/K')
            
            df_net_tmp = pd.DataFrame([[sec,'CTP bins',name,str(np.round(avgDATA,3))]],columns=['type','bin','decomp',used_models])
            df_net_all = pd.concat([df_net_all,df_net_tmp])
              
        plt.savefig(figdir+'NET_'+sec+'_cld_fbk_maps-'+used_models+'-'+case_stamp+'.png', bbox_inches='tight')
    

    out1.close()
    out2.close()
    
    print(df_lw_all.head())
    print(df_sw_all.head())
    print(df_net_all.head())
    
    df_lw_all.to_csv(outdir+'decomp_global_mean_lw_'+case_stamp+'.csv')
    df_sw_all.to_csv(outdir+'decomp_global_mean_sw_'+case_stamp+'.csv')
    df_net_all.to_csv(outdir+'decomp_global_mean_net_'+case_stamp+'.csv')


    ###########################################################################
    # Compute obscuration feedback components
    ###########################################################################  
    sec='LO680' # this should already be true, but just in case...
    PP=sec_dic[sec]  
    PP=sec_dic[sec]   
    print('Get Obscuration Terms')
    CTL,FUT = clisccp1_grd,clisccp1_grd+anomclisccp_grd
    LWK = LWkernel_map[:,:,PP,:]
    SWK = SWkernel_map[:,:,PP,:]
    dTs = MV.average(anomtas_ann_gm)
    obsc_output={}
    obsc_output[sec] = do_obscuration_calcs(CTL,FUT,LWK,SWK,dTs)
    
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if __name__ == "__main__":

    import cases_lookup as CL 

    RadKernel_dir = '/qfs/people/qiny108/diag_feedback_E3SM/CloudRadKernel_input/'
    direc_data = '/compyfs/qiny108/diag_feedback_E3SM_postdata/'
    case_stamp = 'v2.OutTend'
    yearS = 2
    yearE = 6
    fname1 = CL.get_lutable(case_stamp,'amip')
    fname2 = CL.get_lutable(case_stamp,'amip4K')
    outdir = './test1003/'
    figdir = './test1003/'
    
    result = CloudRadKernel(RadKernel_dir,direc_data,case_stamp,yearS,yearE,fname1,fname2,outdir,figdir)

