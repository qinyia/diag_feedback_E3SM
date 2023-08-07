#****************************************************************
#
#    Filename: cal_RadKern_regime_new.py 
#
#    Author: Yi Qin - qin4@llnl.gov
#    Description: Used to do radiative kernel analysis. 
#    Input: tas - Surface air temperature [K]
#           rlut - TOA outgoing longwave radiation [W/m2]
#           rlutcs - TOA clear-sky outgoing longwave radiation [W/m2]
#           rsutcs - TOA clear-sky upward shortwave radiation [W/m2]
#           rsus - surface upward shortwave radiation [W/m2]
#           rsds - surface downward shortwave radiation [W/m2]
#           ps - surface pressure [Pa]
#           rsdt - TOA downward shortwave radiation [W/m2]
#           psl - sea level pressure [Pa]
#           ta - atmospheric temperature [K]
#           hus - atmospheric humidity [kg/kg]
#           OMEGA - vertical velocity [Pa/s]
#           Z3 - geopotential height [m]
#    Output: 
#           Spatial map and global mean of Planck, lapse rate, water vapor, cloud feedbacks 
#
#    Create: 2021-11-09 10:04:40
#    Last Modified: 2023-08-02 14:13:40
#****************************************************************

import numpy as np
import pandas as pd
import netCDF4
import xarray as xr
from global_land_mask import globe
import cartopy.crs as ccrs
import os
import sys
sys.path.append('../')
import psutil 
import glob
import pickle 
import csv
import matplotlib.pyplot as plt 
import matplotlib as mpl
mpl.use('Agg')

import cases_lookup as CL
import time
from loguru import logger
from PlotDefinedFunction import area_averager

### Horizontally regrid data into one defined grid first 
do_hor_regrid = True

### define uniform horizontal grids 
if do_hor_regrid:
#    Dlats = np.arange(-90,92.5,2.5)
#    Dlons = np.arange(0,360,2.5)

    Dlats = np.arange(-87.5,90,2.5)
    Dlons = np.arange(1.25,360,2.5)

### Define standard pressure levels  [hPa]
Dlevs = np.array([100000, 92500, 85000, 70000, 60000, 50000, 40000, 30000, 25000, 20000, 15000, 10000, 7000, 5000, 3000, 2000, 1000])/100.

### read necessary variables from models [month,(level),lat,lon]
var2d = ["tas","rlut","rsut","rlutcs","rsutcs","rsus","rsds","ps","rsdt","psl","ts"]
var3d = ["ta","hus","OMEGA","Z3"]

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def RadKernel(kernel_dir,direc_data,case_stamp,yearS,yearE,fname1,fname2,outdir,figdir,exp1,exp2):

    outfile_map = "RadKern_map_"+case_stamp+".nc"
    outfile_gm  = "RadKern_gm_"+case_stamp+".csv" 
    if os.path.isfile(outdir+outfile_map) and os.path.isfile(outdir+outfile_gm):
        print('RadKenel anlaysis is done.')
        return 

    logger.remove()
    fmt = "<green>{time:YYYY-MM-DD HH:mm:ss}</green> | <cyan>{level}</cyan> | {message} |{elapsed}"
    logger.add(sys.stdout, format=fmt)
    #logger.add("NewRegress",format=fmt)
    ###################

    yearS_4d = "{:04d}".format(yearS)
    yearE_4d = "{:04d}".format(yearE)
    nyears = yearE - yearS + 1

    print(direc_data,fname1)
    direc_data1 = direc_data+'/'+fname1+'/'
    direc_data2 = direc_data+'/'+fname2+'/'

    monS = 1
    monE = 12
    monS_2d='{:02d}'.format(monS)
    monE_2d='{:02d}'.format(monE) 

    Vars = var2d + var3d

    my_timer = Timer()
    #=============================================================
    # read model's data
    #=============================================================
    dic_mod = read_data_e3sm(Vars,var2d,var3d,direc_data1,direc_data2,exp1,exp2,yearS_4d,monS_2d,yearE_4d,monE_2d,nyears)
    
    sub_name = 'Read model data'
    print_memory_status(sub_name)

    print(dic_mod.keys())

    #=============================================================
    # Regridding model's data
    #=============================================================
    dic_mod_fn = {}
    for svar in dic_mod.keys():
        data = dic_mod[svar]

        if do_hor_regrid:
            # horizontal regrid model data to (Dlats,Dlons)
            data_grd = data.interp(lat=Dlats, lon=Dlons,kwargs={"fill_value": "extrapolate"})
        else:
            data_grd = data
    
        if len(data.shape) == 3:
            dic_mod_fn[svar] = data_grd 
        elif len(data.shape) == 4:
            # vertical regrid model data to Dlevs
            if data_grd.coords['lev'].max().values > 1300: # pressure in Pa
                logger.info(f'convert level from Pa to hPa')
                data_grd = data_grd.assign_coords({"lev":data_grd.coords['lev'].values/100})
            data_grd_vrt = data_grd.interp(lev=Dlevs,kwargs={"fill_value": "extrapolate"})
            dic_mod_fn[svar] = data_grd_vrt 
       
    sub_name = 'Regrid model data'
    print_memory_status(sub_name)

    del dic_mod

    #=============================================================
    # get coordinates
    #=============================================================
    ntime,nlev,nlat,nlon = dic_mod_fn['ta_pi'].shape
    logger.info(f'ntime={ntime},nlev={nlev},nlat={nlat},nlon={nlon}')
    
    #=============================================================
    # calculate annual global-mean surface temp anomaly
    #=============================================================
    dtas = dic_mod_fn['tas_ab'] - dic_mod_fn['tas_pi']
    # define an identifiable time
    time_coord = pd.date_range("1850-01", periods=dtas.shape[0], name="time",freq='MS')
    dtas = dtas.assign_coords({"time":time_coord})
    dtas_ann = dtas.groupby('time.year').mean('time') #[year,lat,lon]
    dtas_avg = area_averager(dtas_ann) #[years]

    if do_hor_regrid:
        dtas_ann_grd = dtas_ann.interp(lat=Dlats,lon=Dlons,kwargs={"fill_value": "extrapolate"})
    else:
        dtas_ann_grd = dtas_ann

    logger.info(f'dtas_avg.shape={dtas_avg.shape},dtags_avg_minmax={dtas_avg.min().values},{dtas_avg.max().values},dtas_avg={np.mean(dtas_avg.data)}')

    #=============================================================
    # define coordinates
    #=============================================================
    timec = dic_mod_fn['tas_pi'].coords['time'].data
    latc = dic_mod_fn['tas_pi'].coords['lat'].data
    lonc = dic_mod_fn['tas_pi'].coords['lon'].data

    coords_3d = {"time": timec, "lat": ("lat",latc), "lon": ("lon",lonc)}
    dims_3d   = ["time","lat","lon"]
    coords_4d = {"time": timec, "lev": ("lev",Dlevs), "lat": ("lat",latc), "lon": ("lon",lonc)}
    dims_4d   = ["time","lev","lat","lon"]
    
    # temporal coordinate for dp4d 
    #Dlevs_mid = (Dlevs[:-1] + Dlevs[1:])/2.
    Dlevs_mid = Dlevs[:-1]
    coords_4d_dp = {"time": timec, "lev": ("lev",Dlevs_mid), "lat": ("lat",latc), "lon": ("lon",lonc)}
    
    #=============================================================
    # expand tas to verticals
    #=============================================================
    tas_pi = dic_mod_fn['tas_pi']
    tas_pi_4d = xr.DataArray(np.transpose(np.tile(tas_pi.data,(nlev,1,1,1)),(1,0,2,3)),coords=coords_4d,dims=dims_4d)

    tas_ab = dic_mod_fn['tas_ab']
    tas_ab_4d = xr.DataArray(np.transpose(np.tile(tas_ab.data,(nlev,1,1,1)),(1,0,2,3)),coords=coords_4d,dims=dims_4d)

    # mask tas_ano_grd_4d where ta_ano_vert_grd is True
    dic_mod_fn['tas4d_pi'] = tas_pi_4d.where((dic_mod_fn['ta_pi'].notnull()) & (tas_pi_4d.notnull()))
    dic_mod_fn['tas4d_ab'] = tas_ab_4d.where((dic_mod_fn['ta_ab'].notnull()) & (tas_ab_4d.notnull()))

    dic_mod_fn['ta_pi'] = dic_mod_fn['ta_pi'].where((dic_mod_fn['ta_pi'].notnull()) & (tas_pi_4d.notnull()))
    dic_mod_fn['ta_ab'] = dic_mod_fn['ta_ab'].where((dic_mod_fn['ta_ab'].notnull()) & (tas_ab_4d.notnull()))

    dic_mod_fn['dt_ab'] = dic_mod_fn['ta_ab'] - dic_mod_fn['tas4d_ab']
    dic_mod_fn['dt_pi'] = dic_mod_fn['ta_pi'] - dic_mod_fn['tas4d_pi']

    sub_name = 'Expand tas to verticals'
    print_memory_status(sub_name)

    #=============================================================
    # calculate q1k
    #=============================================================
#    dic_mod_fn['q1k_pi'] = cal_q1k_Mark(dic_mod_fn)
    dic_mod_fn['q1k_pi'] = cal_q1k_Yi(dic_mod_fn)
    
    sub_name = 'Calculate q1k'
    print_memory_status(sub_name)

    #=============================================================
    # calculate dlogq2 
    #=============================================================
    dlogq2_pi,dlogq2_ab = cal_dlogq2_separate(dic_mod_fn)
    dic_mod_fn['dlogq2_pi'] = dlogq2_pi
    dic_mod_fn['dlogq2_ab'] = dlogq2_ab

    sub_name = 'calculate dlogq2'
    print_memory_status(sub_name)

    #=============================================================
    # calculate tropopause height
    #=============================================================
    p_tropo = get_tropopause_pressure(dic_mod_fn['ta_ab'])
    logger.info(f'p_tropo.shape={p_tropo.shape}')
    tropo4d = np.transpose(np.tile(p_tropo,(nlev,1,1,1)),(1,0,2,3)) # [time,lev,lat,lon]
    # set as xarray DataArray
    tropo4d = xr.DataArray(tropo4d, coords=coords_4d, dims=dims_4d)
    
    sub_name = 'Getting Tropopause Pressure'
    print_memory_status(sub_name)

    #=============================================================
    # read radiative kernel data
    #=============================================================
    model_lat = dic_mod_fn['tas_pi'].coords['lat'].values
    model_lon = dic_mod_fn['tas_pi'].coords['lon'].values
    dic_kern = read_kernel_data(kernel_dir, nyears, coords_3d,dims_3d,coords_4d,dims_4d,model_lat,model_lon)
    
    sub_name = 'Read radiative kernel data'
    print_memory_status(sub_name)

    #=============================================================
    # mask out stratosphere 
    #=============================================================
    dic_mod_tropo, dic_kern_tropo = apply_mask_tropo(Dlevs,dic_mod_fn,dic_kern,tropo4d,coords_4d,dims_4d)
    
    del dic_mod_fn, dic_kern

    sub_name = 'Mask out stratosphere'
    print_memory_status(sub_name)

    #=============================================================
    # get pressure thickness
    #=============================================================
    trop_wts, atm_wts = get_weights_SPC(dic_mod_tropo['ps_pi'], p_tropo, dic_mod_tropo['ta_pi'])
    dp4d = trop_wts * 100. # convert to per hPa
    dp4d = xr.DataArray(dp4d,coords=coords_4d_dp, dims=dims_4d) 
    dic_mod_tropo['dp4d_pi'] = dp4d # add into dictionary to be sorted by omega

    sub_name = 'Get pressure thickness'
    print_memory_status(sub_name)

    #=============================================================
    # calculate Kernel * xx_ab and Kernel * xx_pi first 
    #=============================================================
    dic_rad = cal_Kern_rad_separate(dic_mod_tropo, dic_kern_tropo,dic_mod_tropo['dp4d_pi'],dic_mod_tropo['rsdt_ab'])

    sub_name = 'calculate kernel*xx_ab and kernel*xx_pi first'
    print_memory_status(sub_name)

    dic_rad_wap = dic_rad
    del dic_mod_tropo, dic_kern_tropo

    time_coord = pd.date_range("1850-01", periods=dtas_avg.shape[0]*12, name="time",freq='MS')

    #=============================================================
    # calculate anomalies (xx_ab - xx_pi)
    #=============================================================
    dic_rad_wap_ano = get_anomalies(dic_rad_wap)

    sub_name = 'get anomalies: kernel*xx_ab minus kernel*xx_pi'
    print_memory_status(sub_name)

    #=============================================================
    # regression against global-mean surface air temperature anomaly [yearly]
    #=============================================================
    dic_rad_perK = regress_tas(case_stamp,dic_rad_wap_ano,dtas_avg,dtas_ann_grd,time_coord)

    sub_name = 'Regression onto tas'
    print_memory_status(sub_name)

    #=============================================================
    # get final outputs
    #=============================================================
    dic_final = get_outputs(dic_rad_perK)

    sub_name = 'Get outputs'
    print_memory_status(sub_name)
    
    #=============================================================
    # print regional/global values for test                    
    #=============================================================
    prints_gm(dic_final)
    plt_figure(dic_final,case_stamp,figdir)

    sub_name = 'Plot figures'
    print_memory_status(sub_name)

    #=============================================================
    # output data into csv and NC files
    #=============================================================
    output_file(outdir,dic_final,case_stamp, outfile_map, outfile_gm)
        
    sub_name = 'Output data'
    print_memory_status(sub_name)

    return None

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def output_file(outdir,dic_final,model_out, outfile_map, outfile_gm):

    # output spatial 
    outname = ['fbk']
    for iname,name in enumerate(outname):
        da = xr.Dataset(
        data_vars = {
        "T":               (('lat','lon'),dic_final['T_'+name].data),
        "Planck":          (('lat','lon'),dic_final['Planck_'+name].data),
        "LR":              (('lat','lon'),dic_final['LR_'+name].data),
        "WV":              (('lat','lon'),dic_final['WV_'+name].data),
        "ALB":             (('lat','lon'),dic_final['ALB_'+name].data),
        "LW_adj":          (('lat','lon'),dic_final['LW_adj_'+name].data),
        "SW_adj":          (('lat','lon'),dic_final['SW_adj_'+name].data),
        "net_adj":         (('lat','lon'),dic_final['net_adj_'+name].data),
        "SWCRE":           (('lat','lon'),dic_final['SWCRE_'+name].data),
        "LWCRE":           (('lat','lon'),dic_final['LWCRE_'+name].data),
        "netCRE":          (('lat','lon'),dic_final['netCRE_'+name].data),
        "SWCRE_adj":       (('lat','lon'),dic_final['SWCRE_adj_'+name].data),
        "LWCRE_adj":       (('lat','lon'),dic_final['LWCRE_adj_'+name].data),
        "netCRE_adj":      (('lat','lon'),dic_final['netCRE_adj_'+name].data),
        "SW_resd":         (('lat','lon'),dic_final['SW_cld_dir_sum_'+name].data),
        "LW_resd":         (('lat','lon'),dic_final['LW_cld_dir_sum_'+name].data),
        "net_resd":        (('lat','lon'),dic_final['net_cld_dir_sum_'+name].data),
        "T_clr":           (('lat','lon'),dic_final['T_clr_'+name].data),
        "Planck_clr":      (('lat','lon'),dic_final['Planck_clr_'+name].data), 
        "LR_clr":          (('lat','lon'),dic_final['LR_clr_'+name].data),
        "WV_clr":          (('lat','lon'),dic_final['WV_clr_'+name].data),
        "ALB_clr":         (('lat','lon'),dic_final['ALB_clr_'+name].data),
        "WV_SW":           (('lat','lon'),dic_final['WV_sw_'+name].data),
        "WV_LW":           (('lat','lon'),dic_final['WV_lw_'+name].data),
        "WV_clr_SW":       (('lat','lon'),dic_final['WV_sw_clr_'+name].data),
        "WV_clr_LW":       (('lat','lon'),dic_final['WV_lw_clr_'+name].data),
        "Planck_fxRH":     (('lat','lon'),dic_final['Planck_fxRH_'+name].data),
        "LR_fxRH":         (('lat','lon'),dic_final['LR_fxRH_'+name].data),
        "RH":              (('lat','lon'),dic_final['netRH_'+name].data),
        "Planck_clr_fxRH": (('lat','lon'),dic_final['Planck_fxRH_clr_'+name].data),
        "LR_clr_fxRH":     (('lat','lon'),dic_final['LR_fxRH_clr_'+name].data),
        "RH_clr":          (('lat','lon'),dic_final['netRH_clr_'+name].data),
        "LW_clr_sum":      (('lat','lon'),dic_final['LW_clr_sum_'+name].data),
        "SW_clr_sum":      (('lat','lon'),dic_final['SW_clr_sum_'+name].data),
        "net_clr_sum":     (('lat','lon'),dic_final['net_clr_sum_'+name].data),
        "LW_clr_dir":      (('lat','lon'),dic_final['LW_clr_dir_'+name].data),
        "SW_clr_dir":      (('lat','lon'),dic_final['SW_clr_dir_'+name].data),
        "net_clr_dir":     (('lat','lon'),dic_final['net_clr_dir_'+name].data),
        "LW_cld_sum":      (('lat','lon'),dic_final['LW_cld_sum_'+name].data),
        "SW_cld_sum":      (('lat','lon'),dic_final['SW_cld_sum_'+name].data),
        "net_cld_sum":     (('lat','lon'),dic_final['net_cld_sum_'+name].data),
        "LW_cld_dir":      (('lat','lon'),dic_final['LW_cld_dir_'+name].data),
        "SW_cld_dir":      (('lat','lon'),dic_final['SW_cld_dir_'+name].data),
        "net_cld_dir":     (('lat','lon'),dic_final['net_cld_dir_'+name].data),
        "LW_clr_resd":     (('lat','lon'),dic_final['LW_clr_dir_sum_'+name].data),
        "SW_clr_resd":     (('lat','lon'),dic_final['SW_clr_dir_sum_'+name].data),
        "net_clr_resd":    (('lat','lon'),dic_final['net_clr_dir_sum_'+name].data),
        },
        coords = {
        "lat": dic_final['ALB_'+name].coords['lat'].values,
        "lon": dic_final['ALB_'+name].coords['lon'].values,
        }
        )
        
        da.to_netcdf(outdir+outfile_map)

    varlist = ['T','Planck','LR','WV','ALB','LW_adj','SW_adj','net_adj','SWCRE','LWCRE','netCRE',
    'SWCRE_adj','LWCRE_adj','netCRE_adj','LW_cld_dir_sum','SW_cld_dir_sum','net_cld_dir_sum',
    'T_clr','Planck_clr','LR_clr','WV_clr','ALB_clr',
    'WV_sw','WV_lw','WV_sw_clr','WV_lw_clr','Planck_fxRH','LR_fxRH','netRH','Planck_fxRH_clr','LR_fxRH_clr','netRH_clr',
    'LW_clr_sum','SW_clr_sum','net_clr_sum','LW_clr_dir','SW_clr_dir','net_clr_dir','LW_cld_sum','SW_cld_sum','net_cld_sum',
    'LW_cld_dir','SW_cld_dir','net_cld_dir','LW_clr_dir_sum','SW_clr_dir_sum','net_clr_dir_sum']

    varlist_out = ['T','Planck','LR','WV','ALB','dLW_adj','dSW_adj','dnet_adj','SWCRE','LWCRE','netCRE',
    'SWCRE_adj','LWCRE_adj','netCRE_adj','LW_resd','SW_resd','net_resd',
    'T_clr','Planck_clr','LR_clr','WV_clr','ALB_clr',
    'WV_SW','WV_LW','WV_clr_SW','WV_clr_LW','Planck_fxRH','LR_fxRH','RH','Planck_clr_fxRH','LR_clr_fxRH','RH_clr',
    'LW_clr_sum','SW_clr_sum','net_clr_sum','LW_clr_dir','SW_clr_dir','net_clr_dir','LW_cld_sum','SW_cld_sum','net_cld_sum',
    'LW_cld_dir','SW_cld_dir','net_cld_dir','LW_clr_resd','SW_clr_resd','net_clr_resd']

    # output global mean 
    with open(outdir+outfile_gm, 'w') as f:
        f.write("%s,%s\n"%('var','E3SM-1-0'))

        for var,var_out in zip(varlist, varlist_out):
            f.write("%s,%s\n"%(var_out,dic_final[var+'_gfbk'].values))

    return None

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def plt_figure(dic_final,case_stamp,figdir):

    # each variable
    fig = plt.figure(figsize=(12,8))
    axl = fig.add_subplot(1,1,1)

    varx = ['Planck','LR','WV','ALB','SWCRE_adj','LWCRE_adj','netCRE_adj','SWCRE','LWCRE','netCRE']
    colors = ['tab:blue','tab:cyan','tab:orange','tab:grey','blue','red','black','blue','red','black']
    ls = ['-','-','-','-','-','-','-','--','--','--']
    for ivar,svar in enumerate(varx):
        svar1 = svar+'_fbk'
        #ax = fig.add_subplot(4,3,ivar+1,projection=ccrs.PlateCarree(180))
        #datap = dic_final[svar1]
        #lons = datap.lon
        #lats = datap.lat
        #im = ax.contourf(lons,lats,datap,transform=ccrs.PlateCarree())
        #ax.coastlines()
        #ax.set_global()
        #fig.colorbar(im,fraction=0.02)
        #ax.set_title(svar)

        datap = dic_final[svar1].mean(axis=1)
        lats = datap.lat.values
        axl.plot(lats,datap,c=colors[ivar],label=svar,ls=ls[ivar])
        axl.axhline(y=0,ls='--',c='grey')
        axl.set_xlabel('Latitude')
        axl.set_ylabel('Feedbacks [W/m2/K]')

    axl.legend(ncol=2)
    fig.savefig(figdir+'Zonalmean-allfbks-'+case_stamp+'.png', dpi = 300)

    #==============================================================
    # clear-sky linearity test 
    fig = plt.figure(figsize=(12,9))
    nrow = 3
    ncol = 1

    var1 = ['LW_clr_sum','SW_clr_sum','net_clr_sum']
    var2 = ['LW_clr_dir','SW_clr_dir','net_clr_dir']
    var3 = ['LW_clr_dir_sum', 'SW_clr_dir_sum', 'net_clr_dir_sum']
    
    for ivar,svar in enumerate(var1):
        ax = fig.add_subplot(nrow,ncol,ivar+1)

        svar1 = var1[ivar]+'_fbk'
        svar2 = var2[ivar]+'_fbk'
        svar3 = var3[ivar]+'_fbk'
        print(svar1, svar2, svar3) 

        data = np.array([dic_final[svar1],dic_final[svar2],dic_final[svar3]])
        labels = ['Kernel','Model','Model - Kernel']

        for ii,da in enumerate(data):
            if ii == 2:
                marker = '.'
                ls = ':'
            else:
                marker = '.'
                ls = '-'
            
            ax.plot(lats,da.mean(axis=1),marker=marker,ls=ls,
            label=svar1.split('_')[0]+': '+labels[ii])

        ax.legend()

        ax.axhline(y=0,ls=':',color='grey')

        ax.set_xlabel('Latitude')
        ax.set_ylabel('Feedback [W/m2/K]')

    fig.savefig(figdir+'ClrSkyLinearityTest-'+case_stamp+'.png', dpi = 300)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def prints_gm(dic_final):
    for svar in dic_final.keys():
        if 'gfbk' in svar:
            logger.debug(f'{svar}.mean = {dic_final[svar]}')

    str_len = 60
    print('------ Hi, summary all feedbacks except for cloud feedback------------')
    print('Planck feedback: '.rjust(str_len),              dic_final['Planck_gfbk'].values,'W/m2/K')
    print('Lapse rate feedback: '.rjust(str_len),          dic_final['LR_gfbk'].values,'W/m2/K')
    print('Lapse rate + Planck feedback: '.rjust(str_len), dic_final['LR_gfbk'].values+dic_final['Planck_gfbk'].values,'W/m2/K')
    print('Temperature feedback: '.rjust(str_len),         dic_final['T_gfbk'].values,'W/m2/K')
    print('Water vapor feedback: '.rjust(str_len),         dic_final['WV_gfbk'].values,'W/m2/K')
    print("Surface albedo feedback: ".rjust(str_len),      dic_final['ALB_gfbk'].values, "W/m2/K")
    
    print('fixedRH Planck feedback: '.rjust(str_len),              dic_final['Planck_fxRH_gfbk'].values,'W/m2/K')
    print('fixedRH Lapse rate feedback: '.rjust(str_len),          dic_final['LR_fxRH_gfbk'].values,'W/m2/K')
    print('fixedRH Lapse rate + Planck feedback: '.rjust(str_len), dic_final['LR_fxRH_gfbk'].values+dic_final['Planck_fxRH_gfbk'].values,'W/m2/K')
    print('fixedRH Water vapor feedback: '.rjust(str_len),         dic_final['WV_fxRH_gfbk'].values,'W/m2/K')
    print('fixedRH RH feedback: '.rjust(str_len),                  dic_final['netRH_gfbk'].values,'W/m2/K')
    
    print('--------clear-sky component----------------------------------------------')
    print('clr Planck feedback: '.rjust(str_len),                  dic_final['Planck_clr_gfbk'].values,'W/m2/K')
    print('clr Lapse rate feedback: '.rjust(str_len),              dic_final['LR_clr_gfbk'].values,'W/m2/K')
    print('clr Lapse rate + clr Planck feedback: '.rjust(str_len), dic_final['LR_clr_gfbk'].values+dic_final['Planck_clr_gfbk'].values,'W/m2/K')
    print('clr Temperature feedback: '.rjust(str_len),             dic_final['T_clr_gfbk'].values,'W/m2/K')
    print('clr Water vapor feedback: '.rjust(str_len),             dic_final['WV_clr_gfbk'].values,'W/m2/K')
    print("clr Surface albedo feedback: ".rjust(str_len),          dic_final['ALB_clr_gfbk'].values,"W/m2/K")
    
    print('fixedRH clr Planck feedback: '.rjust(str_len),              dic_final['Planck_fxRH_clr_gfbk'].values,'W/m2/K')
    print('fixedRH clr Lapse rate feedback: '.rjust(str_len),          dic_final['LR_fxRH_clr_gfbk'].values,'W/m2/K')
    print('fixedRH clr Lapse rate + Planck feedback: '.rjust(str_len), dic_final['LR_fxRH_clr_gfbk'].values+dic_final['Planck_fxRH_clr_gfbk'].values,'W/m2/K')
    print('fixedRH clr Water vapor feedback: '.rjust(str_len),         dic_final['WV_lw_clr_fxRH_gfbk'].values+dic_final['WV_sw_clr_fxRH_gfbk'].values,'W/m2/K')
    print('fixedRH clr RH feedback: '.rjust(str_len),                  dic_final['netRH_clr_gfbk'].values,'W/m2/K')

    print('--------------- quick check here --------------------------')
    print('fixedRH (Lapse rate + Planck) - (Lapse rate + Planck): '.rjust(str_len), dic_final['LR_fxRH_gfbk'].values + dic_final['Planck_fxRH_gfbk'].values - dic_final['LR_gfbk'].values - dic_final['Planck_gfbk'].values,'W/m2/K')
    print('fixedRH WV: '.rjust(str_len),dic_final['WV_fxRH_gfbk'].values,'W/m2/K')
    print('fixedRH Ta: '.rjust(str_len),dic_final['T_fxRH_gfbk'].values,'W/m2/K')
    
    print('clr fixedRH (Lapse rate + Planck) - (Lapse rate + Planck): '.rjust(str_len), dic_final['LR_fxRH_clr_gfbk'].values + dic_final['Planck_fxRH_clr_gfbk'].values - dic_final['LR_clr_gfbk'].values - dic_final['Planck_clr_gfbk'].values,'W/m2/K')
    print('clr fixedRH WV: '.rjust(str_len), dic_final['WV_lw_clr_fxRH_gfbk'].values + dic_final['WV_sw_clr_fxRH_gfbk'].values,'W/m2/K')
    print('clr fixedRH Ta: '.rjust(str_len), dic_final['T_fxRH_clr_gfbk'].values,'W/m2/K')
    
    print('Lapse rate + Planck + WV feedback: '.rjust(str_len), dic_final['LR_gfbk'].values + dic_final['Planck_gfbk'].values + dic_final['WV_gfbk'].values,'W/m2/K')
    print('fixedRH Lapse rate + Planck + WV feedback: '.rjust(str_len), dic_final['LR_fxRH_gfbk'].values + dic_final['Planck_fxRH_gfbk'].values + dic_final['netRH_gfbk'].values,'W/m2/K')
    
    print('clr Lapse rate + Planck + WV feedback: '.rjust(str_len), dic_final['LR_clr_gfbk'].values + dic_final['Planck_clr_gfbk'].values + dic_final['WV_clr_gfbk'].values,'W/m2/K')
    print('fixedRH clr Lapse rate + Planck + WV feedback: '.rjust(str_len), dic_final['LR_fxRH_clr_gfbk'].values + dic_final['Planck_fxRH_clr_gfbk'].values + dic_final['netRH_clr_gfbk'].values,'W/m2/K')
    
    print('sum of all clear-sky sw feedbacks: '.rjust(str_len), dic_final['SW_clr_sum_gfbk'].values,'W/m2/K')
    print('sum of all clear-sky lw feedbacks: '.rjust(str_len), dic_final['LW_clr_sum_gfbk'].values,'W/m2/K')
    
    print('TOA direct clear-sky sw radiation feedback: '.rjust(str_len), dic_final['SW_clr_dir_gfbk'].values,'W/m2/K')
    print('TOA direct clear-sky lw radiation feedback: '.rjust(str_len), dic_final['LW_clr_dir_gfbk'].values,'W/m2/K')

    print('TOA direct clear-sky sw/sum of all clear-sky sw: '.rjust(str_len), (dic_final['SW_clr_sum_gfbk'].values - dic_final['SW_clr_dir_gfbk'].values)/dic_final['SW_clr_dir_gfbk'].values*100., '%')
    print('TOA direct clear-sky lw/sum of all clear-sky lw: '.rjust(str_len), (dic_final['LW_clr_sum_gfbk'].values - dic_final['LW_clr_dir_gfbk'].values)/dic_final['LW_clr_dir_gfbk'].values*100., '%')

    return None

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def get_outputs(dic_rad_perK):
    dic_final = {}

    outname = ['gfbk','fbk']
    for iname,name in enumerate(outname):

        for key in dic_rad_perK.keys():
            dic_final[key] = dic_rad_perK[key]

        # get final water vapor feedbacks
        dic_final['WV_'+name]    = dic_rad_perK['WV_sw_'+name] + dic_rad_perK['WV_lw_'+name]
        dic_final['WV_clr_'+name]    = dic_rad_perK['WV_sw_clr_'+name] + dic_rad_perK['WV_lw_clr_'+name]
        dic_final['WV_fxRH_'+name]    = dic_rad_perK['WV_sw_fxRH_'+name] + dic_rad_perK['WV_lw_fxRH_'+name]

        dic_final['netRH_'+name]    = dic_rad_perK['netRH_sw_'+name] + dic_rad_perK['netRH_lw_'+name]
        dic_final['netRH_clr_'+name]    = dic_rad_perK['netRH_sw_clr_'+name] + dic_rad_perK['netRH_lw_clr_'+name]


    return dic_final

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 2022-03-25: replace the regression from genutil.statistics.linearregression to new function based on xarray data.
def regress_tas(exp_ctl,dic_rad,dtas_avg,dtas_ann,time_coord):

    # rename dimension name from year to time 
    dtas_avg = dtas_avg.rename({'year':'time'})
    dtas_ann = dtas_ann.rename({'year':'time'})

    dic_rad_perK = {}
    for svar in dic_rad.keys():
        logger.debug(f'We are regressing {svar} on tas')

        data = dic_rad[svar]
        data = data.assign_coords({"time":time_coord})
        # get annual mean for each variable 
        data_ann = data.groupby('time.year').mean('time') #[years,nbins]
        # rename dimension name from year to time 
        data_ann = data_ann.rename({'year':'time'})
        logger.debug(f'data_ann.shape={data_ann.shape},dtas_avg.shape={dtas_avg.shape},dtas_ann.shape={dtas_ann.shape}')

        # regress on tas anomaly 
        if 'coupled' not in exp_ctl:
            slope = data_ann.mean(axis=0)/dtas_avg.mean()
        else:
            _,_,slope,intercept,_,_ = lag_linregress_3D(dtas_avg, data_ann, lagx=0, lagy=0)

        # save
        dic_rad_perK[svar+'_fbk'] = xr.DataArray(slope,coords=data_ann[0,:].coords, dims=data_ann[0,:].dims)

        # get bin average and global average
        latS = -90.
        latE = 90.
        tmp = area_averager(slope)

        # save
        dic_rad_perK[svar+'_gfbk'] = tmp
   
    return dic_rad_perK

#=============================================================
def cal_q1k_Yi(dic_mod):
    hus_pi = dic_mod['hus_pi']
    hus_ab = dic_mod['hus_ab']

    avgta = (dic_mod['ta_pi'] + dic_mod['ta_ab'])/2.0

    levs = dic_mod['ta_pi'].coords['lev'].values
    ntime,nlev,nlat,nlon = dic_mod['ta_pi'].shape
    
    tmp = np.tile(levs,(ntime,nlat,nlon,1))
    levnd = np.transpose(tmp,(0,3,1,2))

    qs0 = r_star_GG(levnd*100.,avgta)
    qs0 = np.float32(qs0)
    ta1k = avgta+1.0
    qs1 = r_star_GG(levnd*100.,ta1k)
    qs1 = np.float32(qs1)
    rh0 = hus_pi/qs0
    q1k = rh0*qs1

    #q1k_xr = xr.DataArray(q1k,coords=ta_pi.coords,dims=ta_pi.dims)
    q1k_xr = q1k

    return q1k_xr

def r_star_GG(p, T):
    Rv = 461.5    # Gas constant of water vapor [J kg^-1 K^-1]
    Rd = 287.04   # Gas constant of dry air [J kg^-1 K^-1]
    epsilon = Rd/Rv

    rr = epsilon*e_star_GG(T)/(p - e_star_GG(T))
    hus = rr/(1.+rr)
    return hus

def e_star_GG(T):
   # use the Goff-Gratch equation: unit of T is K, unit of e_star is hPa.
   # to liquid (-49.9 ~ 100 degC)
   T00 = 273.16
   t1 = 10.79574*(1.-T00/T)
   t2 = 5.02800*np.ma.log10(T/T00)
   t3 = 1.50475*1e-4*(1.-10**(-8.2969*(T/T00-1)))
   t4 = 0.42873*1e-3*(10**(4.76955*(1.-T/T00))-1.)
   t5 = 0.78614
   esw = 10**(t1-t2+t3+t4+t5)

   # to ice (-100 ~ 0 degC)
   t1 = -9.09685*(T00/T-1.)
   t2 = -3.56654*np.ma.log10(T00/T)
   t3 = 0.87682*(1-T/T00)
   t4 = 0.78614
   esi = 10**(t1+t2+t3+t4)

   es = np.where(T<T00-49.9,esi,esw)
   es = np.ma.masked_where(T.to_masked_array().mask == True, es*100.) # convert to Pa
   print(np.nanmin(es), np.nanmax(es))
   return es 

#=============================================================
def cal_q1k_Mark(dic_mod):
    hus_pi = dic_mod['hus_pi']
    hus_ab = dic_mod['hus_ab']

    avgta = (dic_mod['ta_pi'] + dic_mod['ta_ab'])/2.0

    levs = dic_mod['ta_pi'].coords['lev'].values
    ntime,nlev,nlat,nlon = dic_mod['ta_pi'].shape
    
    tmp = np.tile(levs,(ntime,nlat,nlon,1))
    levnd = np.transpose(tmp,(0,3,1,2))

    qs0,qs1 = qsat_blend_Mark(avgta,levnd)
    rh0 = hus_pi/qs0
    q1k = rh0*qs1

    #q1k_xr = xr.DataArray(q1k,coords=ta_pi.coords,dims=ta_pi.dims)
    q1k_xr = q1k

    return q1k_xr

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def get_anomalies(dic_mod_wap):
    dic_mod_ano = {}
    for svar in dic_mod_wap.keys():
        if '_pi' in svar:
            var1 = '_'.join(svar.split('_')[:-1])
            if var1 not in ['ps','q1k','dp4d']: # ignore ps because only read ps_pi
                da1 = dic_mod_wap[var1+'_ab']# * dic_mod_wap[var1+'_ab_N']
                da2 = dic_mod_wap[var1+'_pi']# * dic_mod_wap[var1+'_pi_N']
                diff = da1 - da2 
                dic_mod_ano[var1] = diff 

    return dic_mod_ano

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def apply_mask_tropo(Dlevs,dic_mod_fn,dic_kern,tropo4d,coords_4d,dims_4d):
    '''
    mask data above tropopause. Set mask and nan as zero.
    '''
    ntime,nlev,nlat,nlon = dic_mod_fn['ta_pi'].shape
    lev4d = np.transpose(np.tile(Dlevs,(ntime,nlat,nlon,1)),(0,3,1,2)) 
    logger.debug(f'lev4d = {lev4d.shape}, {np.nanmin(lev4d)}, {np.nanmax(lev4d)}')
    logger.debug(f'tropo4d = {tropo4d.shape}, {np.nanmin(tropo4d)}, {np.nanmax(tropo4d)}')

    dic_mod_tropo = {}
    for svar in dic_mod_fn.keys():
        if len(dic_mod_fn[svar].shape) == 4: # [time,lev,lat,lon]
            dic_mod_tropo[svar] = np.ma.masked_where(lev4d<=tropo4d, dic_mod_fn[svar])
            dic_mod_tropo[svar][dic_mod_tropo[svar].mask] = 0.0
            dic_mod_tropo[svar] = xr.DataArray(dic_mod_tropo[svar], coords=coords_4d, dims=dims_4d)
        else:
            dic_mod_tropo[svar] = dic_mod_fn[svar]


    dic_kern_tropo = {}
    for svar in dic_kern.keys():
        if len(dic_kern[svar].shape) == 4: 
            dic_kern_tropo[svar] = np.ma.masked_where(lev4d<=tropo4d, dic_kern[svar])
            dic_kern_tropo[svar][dic_kern_tropo[svar].mask] = 0.0
            dic_kern_tropo[svar] = xr.DataArray(dic_kern_tropo[svar], coords=coords_4d, dims=dims_4d)
            dic_kern_tropo[svar] = dic_kern_tropo[svar].fillna(0) # set nan as zero
        else:
            dic_kern_tropo[svar] = dic_kern[svar]

    return dic_mod_tropo, dic_kern_tropo

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def get_tropopause_pressure(Ta0):
    '''
    April 27, 2021: source code from Mark
    Perform tropopause pressure calculation
    generated via: f2py -c -m tropo tropo.f90
    tropo.f90 came from http://www.inscc.utah.edu/~reichler/research/projects/TROPO/code.txt

    Note: levs need to be at hPa

    - Yi Qin:
    Oct 13, 2021: revised from tropo.f90 to tropo4d.f90 to process 4-D data [time,lev,lat,lon]
    This one will be faster than looping over lons.
    '''
    import tropo4d

    PP = Ta0.coords['lev'].values
    lons = Ta0.coords['lon'].values

    plimu=45000
    pliml=7500
    plimlex=7500
    dofill=0

    # fortran wants this in lon,lat,lev,time
    tropp=np.ma.zeros((Ta0[:,0,:,:].shape)) #[time,lat,lon]
    tropp=np.ma.masked_where(tropp==0,tropp)
    temp=Ta0.transpose("lon","lat","lev","time")
    tp,tperr=tropo4d.tropo4d(temp,PP,plimu,pliml,plimlex,dofill) # returned tp in [lon,lat,time]
    logger.debug(f'temp.shape={temp.shape},{temp.max().values},{temp.min().values},tp.shape={tp.shape},{np.min(tp)/100},{np.max(tp)/100} hPa')

    tropp=np.transpose(tp,(2,1,0)) # move time to the first dim
    tropp = tropp/100. # convert to hPa

    tropp=np.ma.masked_where(tropp>plimu/100.,tropp)
    tropp=np.ma.masked_where(tropp<0.,tropp)

    return tropp

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def read_kernel_data(kernel_dir, nyears, coords_3d,dims_3d, coords_4d,dims_4d, model_lat, model_lon):
    '''
    read radiative kernel [12,(level),lat,lon]
    input: directory of kernel data, number of model years
    output: dictionary saving all kernel data in (time,(lev),lat,lon)
    '''
    dic_kern = {}
    cases = ['cld','clr']
    kern_vars = ['ts','t','wv_sw','wv_lw','alb']
    Types = ['lw','lw','sw','lw','sw']
    for case in cases:
        for ivar,kern_var in enumerate(kern_vars):
    
            varout = kern_var+'_Kern'+case.capitalize()
            Type = Types[ivar]
    
            fname = "RRTMG_"+kern_var+"_toa_"+case+"_highR.nc"
            f1 = xr.open_dataset(kernel_dir+fname)
            data = f1[Type+'kernel']
    
            ## The kernel data's latitude is from north to south [90....-90]
            ## we need to reverse it.
            if len(data.shape) == 3: #(time,lat,lon)
                # horizontally regrid
                data_grd = data[:,::-1,:].interp(lat=model_lat,lon=model_lon,kwargs={"fill_value": "extrapolate"})
                # rename dimension name
                data_grd = data_grd.rename({'month':'time'})
                dic_kern[varout] = data_grd

            else: #(time,lev,lat,lon)
                # horizontally regrid
                data_grd = data[:,:,::-1,:].interp(lat=model_lat,lon=model_lon,kwargs={"fill_value": "extrapolate"})
                # rename dimension name
                data_grd = data_grd.rename({'player':'lev','month':'time'})
                ## vertical regrid to standard levels 
                if data_grd.coords['lev'].max().values > 1300: # pressure in Pa
                    data_grd = data_grd.assign_coords({"lev":data_grd.coords['lev'].values/100})
                data_grd_vrt = data_grd.interp(lev=Dlevs,kwargs={"fill_value": "extrapolate"})
                dic_kern[varout] = data_grd_vrt
    
            logger.debug(f'{varout}, {dic_kern[varout].shape},{dic_kern[varout].min().values},{dic_kern[varout].max().values}')
            
    # mask out arrays when one has missing value at this level but the other does not 
    dic_kern['wv_lw_KernCld'] = dic_kern['wv_lw_KernCld'].where(dic_kern['wv_sw_KernCld'].notnull())
    dic_kern['wv_lw_KernClr'] = dic_kern['wv_lw_KernClr'].where(dic_kern['wv_sw_KernClr'].notnull())
    dic_kern['t_KernCld'] = dic_kern['t_KernCld'].where(dic_kern['wv_sw_KernCld'].notnull())
    dic_kern['t_KernClr'] = dic_kern['t_KernClr'].where(dic_kern['wv_sw_KernClr'].notnull())

    del data, data_grd, data_grd_vrt

    ## expand kernel from [12,(level),lat,lon] to [month,(level),lat,lon]
    dic_kern_mon = {}
    for svar in dic_kern.keys():
        data = dic_kern[svar]
        if len(data.shape) == 3: # (time,lat,lon)
            data1 = xr.DataArray(np.tile(data,(nyears,1,1)),coords=coords_3d,dims=dims_3d)
        elif len(data.shape) == 4: #(time,lev,lat,lon)
            data1 = xr.DataArray(np.tile(data,(nyears,1,1,1)),coords=coords_4d,dims=dims_4d)

        dic_kern_mon[svar] = data1

    ## calculate total water vapor kernel
    dic_kern_mon['wv_KernCld'] = dic_kern_mon['wv_lw_KernCld'] + dic_kern_mon['wv_sw_KernCld']
    dic_kern_mon['wv_KernClr'] = dic_kern_mon['wv_lw_KernClr'] + dic_kern_mon['wv_sw_KernClr']

    ## calculate T+Q kernel
    dic_kern_mon['Tq_KernCld'] = dic_kern_mon['t_KernCld'] + dic_kern_mon['wv_KernCld']
    dic_kern_mon['Tq_KernClr'] = dic_kern_mon['t_KernClr'] + dic_kern_mon['wv_KernClr']

    return dic_kern_mon

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def get_weights_SPC(ps, trop, field):
    """
    Jul 06, 2021: copied from /home/zelinka1/scripts/CMIP6_utils.py to reduce the dependency on other py files.
    from Stephen Po-Chedley
    wts = get_weights(ps, trop, field)
    Function returns the layer weights (layer thickness) given the upper
    integration boundary (e.g., the tropopause, trop) and the lower bound
    (e.g., surface pressure, ps). Uses the input field to get pressure level
    information.
         ps - response field [time, lat, lon]
         trop - radiative kernel [time, lat, lon]
         field - field that will be integrated [time, plev, lat, lon]
         Dlevs - vertical pressure levels [plev]
    """     
            
    PAL = np.ma.zeros(ps.shape)
    PAT = np.ma.ones(ps.shape)          
    trop_wts = np.ma.zeros(field.shape - np.ma.array([0,1,0,0]))
    atm_wts = np.ma.zeros(field.shape - np.ma.array([0,1,0,0]))
    plev = field.coords['lev'].values
    
    # make sure pressures are in Pa.  
    if np.ma.max(plev)<=2000:
        plev=100*plev         
    if np.ma.max(ps)<=2000:
        ps=100*ps
    if np.ma.max(trop)<=2000:
        trop=100*trop
            
    if plev[0] < plev[1]:
        raise ValueError('This script assumes that plev[0] > plev[1].')
    for i in range(len(plev)-1):
        # allocate slice
        sfield = np.ma.zeros(ps.shape)
        tfield = np.ma.zeros(ps.shape)
        # get first level above surface
        p = plev[i]
        pp = plev[i+1]
        ISURF = np.ma.greater(ps, p) # 1 where ps>p, 0 where ps<p
        II = np.ma.equal(PAL, 0)
        IA = np.ma.logical_and(ISURF, II)                  
        PAL = np.ma.where(IA, pp, PAL)
        # get first layer below tropopause
        ITROP = np.ma.less(trop, p)
        II = np.ma.greater_equal(trop, pp)
        IB = np.ma.logical_and(ITROP, II)
        PAT = np.ma.where(IB, p, PAT)
        # layers above tropopause or below surface (i.e., zero weight)
        above_trop = np.ma.logical_not(ITROP)
        below_surf = np.ma.logical_not(ISURF)
        IC = np.ma.logical_or(below_surf,above_trop)
        # layers adjacent to both tropopause and surface (i.e., layer is between tropopause and surface)
        ID = np.ma.logical_and(IA, IB)
        # layers not adjacent to tropopause or surface (i.e., full layers)
        IE = np.ma.logical_or(IA, IB)
        IE = np.ma.logical_and(np.ma.logical_not(IC), np.ma.logical_not(IE)) 
        # layers not adjacent to surface (i.e., full layers)
        IF = np.ma.logical_and(np.ma.logical_not(below_surf), np.ma.logical_not(IA)) 
        # TROPOSPHERIC VALUES
        sfield = np.ma.where(IA, ps-PAL, sfield)
        sfield = np.ma.where(IB, PAT-trop, sfield)
        sfield = np.ma.where(IC, 0., sfield)
        sfield = np.ma.where(ID, ps-trop, sfield)
        sfield = np.ma.where(IE, p - pp, sfield)
        # store field and weight by per 100 hPa (1/100 for Pa to hPa and 1/100 for per *100* hPa)
        trop_wts[:, i, :, :] = sfield / 10000.
        
        # ATMOSPHERIC VALUES
        tfield = np.ma.where(IA, ps-PAL, tfield)
        tfield = np.ma.where(below_surf, 0., tfield)
        tfield = np.ma.where(IF, p - pp, tfield)
        # store field and weight by per 100 hPa (1/100 for Pa to hPa and 1/100 for per *100* hPa)
        atm_wts[:, i, :, :] = tfield / 10000.
        
    return trop_wts,atm_wts

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def VertSum(varin,dp4d):
    # get midpoint var first
    var1 = (varin[:,:-1,:].data + varin[:,1:,:].data)/2.
    # vertical integral
    outvar = np.ma.sum(var1 * dp4d/100., axis=1)

    return outvar 

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def qsat_water(t,p):
    """
    Jul 7, 2021: copied from /home/zelinka1/scripts/MDZ_utils.py to reduce the dependency on other py files 

    saturation mixing ratio w/respect to liquid water.
    I checked that this code gives answers very similar 
    to its counterpart in CAM5. This code uses the Goff 
    + Gratch (1946) method, which is the default
    version used in CAM. If 4d, uses a pre-compiled F90 
    version of the exact same method as used for 3d (for 
    efficiency). 

    INPUTS:
    t - temperature (K)
    p - pressure (Pa)

    OUTPUTS:
    qs - saturation mixing ratio (kg water/kg air)
    """

    # The fortran version has stopped working for reasons I don't understand
    """
    #USE FORTRAN VERSION IF 4D:
    if len(t.shape)==4:

        #MAKE SURE DATATYPE IS WHAT F90 EXPECTED:
        if t.dtype!=np.float32:
            t=t.astype(np.float32)
            print('recasting t to float32!')
        if p.dtype!=np.float32:
            t=t.astype(np.float32)
            print('recasting p to float32!')

        #CONVERT MASKED STUFF TO -999.
        t=t.filled(-999.)
        p=p.filled(-999.)

        #ACTUALLY RUN THE FORTRAN CODE
        import qsat_water4d
        qsat_water = qsat_water4d.qsat_water4d(t,p)

        #CONVERT STUFF BACK TO MASKED ARRAYS
        qsat_water=np.ma.masked_less(qsat_water,0.)

    else:
    """
    #DEFINE THE LOG FUNCTION DEPENDING ON THE DATATYPE OF THE INPUT.
    #=======================================
    #<qinyi 2021-10-27 #------------------
    #if type(t)==type(np.arange(2)):
    #    from MV2 import log10
    #elif type(t)==type(np.arange(2)):
    #    from numpy import log10

    from numpy import log10
    #>qinyi 2021-10-27 #------------------

    epsqs=0.6218847083189506

    ps = 1013.246
    ts = 373.16
    e1 = 11.344*(1.0 - t/ts)
    e2 = -3.49149*(ts/t - 1.0)
    f1 = -7.90298*(ts/t - 1.0)
    f2 = 5.02808*log10(ts/t)
    f3 = -1.3816*(10.0**e1 - 1.0)/10000000.0
    f4 = 8.1328*(10.0**e2 - 1.0)/1000.0
    f5 = log10(ps)
    f  = f1 + f2 + f3 + f4 + f5
    es = (10.0**f)*100.0

    qsat_water = epsqs*es/(p-(1.-epsqs)*es) # saturation w/respect to liquid only

    return qsat_water

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def qsat_ice(t,p):
    """
    Jul 6, 2021: Yi copied from /home/share/PMC_utils.py

    saturation mixing ratio w/respect to ice. Uses the 
    Goff + Gratch (1946) method, which is the default
    version used in CAM. This ver looks different than 
    qsat_water because it is pulled from CAM5.3.01, which
    has been updated.

    INPUTS:
    t - temperature (K)
    p - pressure (Pa)

    OUTPUTS:
    qs - saturation mixing ratio (kg ice/kg air)
    """
    #<qinyi 2021-10-27 #------------------
    #if type(t)==type(np.arange(2)):
    #    from MV2 import log10
    #elif type(t)==type(np.arange(2)) or type(t)==type(273.):
    #    from numpy import log10

    from numpy import log10
    #>qinyi 2021-10-27 #------------------

    epsqs=0.6218847083189506
    h2otrip=273.16 # Triple point temperature of water (K)

    #Compute es (saturation vap pres) in Pa:
    es = 10.0**(-9.09718*(h2otrip/t-1.)-3.56654*
                   log10(h2otrip/t)+0.876793*(1.-t/h2otrip)+
                   log10(6.1071))*100.

    qsat_ice = epsqs*es/(p-(1.-epsqs)*es)

    return qsat_ice

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def get_qsl(avgta, lev_4d):
    print_memory_status('start get_qsl')

    wsl=qsat_water(avgta,lev_4d*100.)
    print_memory_status('get wsl')

    wsl_plus1=qsat_water(avgta+1,lev_4d*100.)
    print_memory_status('get wsl_plus1')

    qsl=wsl/(1+wsl) # convert from mixing ratio (kg/kg) to specific humidity
    print_memory_status('get qsl')

    qsl_plus1=wsl_plus1/(1+wsl_plus1)
    print_memory_status('get qsl_plus1')

    return qsl, qsl_plus1

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def get_qsi(avgta, lev_4d):
    print_memory_status('start get_qsi')

    wsi=qsat_ice(avgta,lev_4d*100.)
    print_memory_status('get wsi')

    wsi_plus1=qsat_ice(avgta+1,lev_4d*100.)
    print_memory_status('get wsi_plus1')

    qsi=wsi/(1+wsi) # convert from mixing ratio (kg/kg) to specific humidity
    print_memory_status('get qsi')

    qsi_plus1=wsi_plus1/(1+wsi_plus1)   
    print_memory_status('get qsi_plus1')

    return qsi, qsi_plus1

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def get_qs0(avgta, lev_4d):

    print_memory_status('start get_qs0')

    qsl,_ = get_qsl(avgta, lev_4d)
    qsi,_ = get_qsi(avgta, lev_4d)
    print_memory_status('get qsl and qsi')

    # Compute blend of qsi and qsl between -40 and 0
    blend = (avgta-233)*qsl/40 + (273-avgta)*qsi/40
    print_memory_status('get blend')

    my_timer = Timer()

    inds = avgta > 233
    inds &= avgta < 273
    qs0 = blend.where(inds,qsi)
    print_memory_status('get qs0 A')
    inds2 = avgta >= 273
    qs0 = qs0.where(~inds2,qsl)
    print_memory_status('get qs0 B')

    #qs0 = np.where((avgta>233) & (avgta<273), blend, qsi)#[0]
    #qs0 = np.where(avgta >= 273, qsl, qs0)#[0]

    return qs0

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def get_qs1(avgta, lev_4d):

    _,qsl_plus1 = get_qsl(avgta, lev_4d)
    _,qsi_plus1 = get_qsi(avgta, lev_4d)

    # Compute blend of qsi and qsl between -40 and 0
    blend_plus1 = (avgta-233)*qsl_plus1/40 + (273-avgta)*qsi_plus1/40

    inds = avgta > 233
    inds &= avgta < 273
    qs1 = blend_plus1.where(inds,qsi_plus1)

    inds2 = avgta >= 273
    qs1 = qs1.where(~inds2,qsl_plus1)

    return qs1

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Oct 4, 2020: change Mark's qsat method as a function
def qsat_blend_Mark(avgta,lev_4d):

    qs0 = get_qs0(avgta, lev_4d)
    qs1 = get_qs1(avgta, lev_4d)

    return qs0, qs1

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def print_memory_status(sub_name):
    mem1 = memory_usage_psutil_Percent()
    mem2 = memory_usage_psutil_MB()
    logger.info(f'{sub_name} is done. Currently using {str(np.round(mem1,6))}%, {str(mem2)} GB')
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

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def delete_vars(delterm,dic_final):
    sub_name = 'before delete item'
    print_memory_status(sub_name)

    for svar in delterm:
        if svar in dic_final.keys():
            logger.debug(f'deleting {svar}')
            del dic_final[svar]

    sub_name = 'after delete item'
    print_memory_status(sub_name)

    return dic_final

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def save_big_dataset(dic_mod,outfile):
    '''
    create a big dataset based on all variables in a dictionary and save to netcdf file.
    '''
    datalist = []
    for svar in dic_mod.keys():
        data = xr.DataArray(dic_mod[svar],name=svar)
        datalist.append(data)

    data_big = xr.merge(datalist)

    #data_big.to_netcdf(outfile,encoding={'time':{'dtype': 'i4'},'bin':{'dtype':'i4'}})
    data_big.to_netcdf(outfile)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def read_pickle(filename):
    pickle_in = open(filename,"rb")
    dic = pickle.load(pickle_in)
    return dic

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def write_pickle(filename,dic):
    pickle_out = open(filename,"wb")
    pickle.dump(dic,pickle_out)
    pickle_out.close()

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
class Timer:
    '''
    Example:
    my_timer = Timer() # Start Timer
    .... do something....
    time_hhmmss = my_timer.get_time_hhmmss()
    print("Time elapsed: %s" % time_hhmmss)
    '''
    def __init__(self):
        self.start = time.time()

    def restart(self):
        self.start = time.time()

    def get_time_hhmmss(self):
        end = time.time()
        m, s = divmod(end - self.start, 60)
        h, m = divmod(m, 60)
        time_str = "%02d:%02d:%02d" % (h, m, s)
        return time_str

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def make_dir(outdir_out):
    '''
    Aug 30, 2020: make a directory if it does not exist
    Nov 15, 2020: make it work for recursive directory
    '''
    try:
#        os.mkdir(outdir_out)
        os.makedirs(outdir_out,exist_ok=True)
    except OSError:
        print("Creation of the directory %s failed" % outdir_out, ", since it already exists.")
    else:
        print("Successfully created the directory %s " % outdir_out)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def lag_linregress_3D(x, y, lagx=0, lagy=0):
    """
    Input: Two xr.Datarrays of any dimensions with the first dim being time. 
    Thus the input data could be a 1D time series, or for example, have three 
    dimensions (time,lat,lon). 
    Datasets can be provided in any order, but note that the regression slope 
    and intercept will be calculated for y with respect to x.
    Output: Covariance, correlation, regression slope and intercept, p-value, 
    and standard error on regression between the two datasets along their 
    aligned time dimension.  
    Lag values can be assigned to either of the data, with lagx shifting x, and
    lagy shifting y, with the specified lag amount. 
    """ 
    #1. Ensure that the data are properly alinged to each other. 
    x,y = xr.align(x,y)
    
    #2. Add lag information if any, and shift the data accordingly
    if lagx!=0:
    
        # If x lags y by 1, x must be shifted 1 step backwards. 
        # But as the 'zero-th' value is nonexistant, xr assigns it as invalid 
        # (nan). Hence it needs to be dropped
        x   = x.shift(time = -lagx).dropna(dim='time')
    
        # Next important step is to re-align the two datasets so that y adjusts
        # to the changed coordinates of x
        x,y = xr.align(x,y)
    
    if lagy!=0:
        y   = y.shift(time = -lagy).dropna(dim='time')
        x,y = xr.align(x,y)
    
    #3. Compute data length, mean and standard deviation along time axis: 
    n = y.notnull().sum(dim='time')
    xmean = x.mean(axis=0)
    ymean = y.mean(axis=0)
    xstd  = x.std(axis=0)
    ystd  = y.std(axis=0)
    
    #4. Compute covariance along time axis
    cov   =  np.ma.sum((x - xmean)*(y - ymean), axis=0)/(n)
    cov   =  xr.DataArray(cov, coords=ymean.coords)
    
    #5. Compute correlation along time axis
    cor   = cov/(xstd*ystd)
    
    #6. Compute regression slope and intercept:
    slope     = cov/(xstd**2)
    intercept = ymean - xmean*slope  
    
    #7. Compute P-value and standard error
    #Compute t-statistics
    tstats = cor*np.ma.sqrt(n-2)/np.ma.sqrt(1-cor**2)
    stderr = slope/tstats
    
    from scipy.stats import t
    pval   = t.sf(tstats, n-2)*2
    pval   = xr.DataArray(pval, dims=cor.dims, coords=cor.coords)
    
    return cov,cor,slope,intercept,pval,stderr

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def cal_Kern_rad_separate(dic_mod_wap,dic_kern_wap,dp4d,rsdt_ab):

    dic_rad = {}
    for case in ['pi','ab']:
        #----------------------------------------------------------
        # Temperature Feedback
        #----------------------------------------------------------
        kern1 = ['ts_KernCld', 'ts_KernClr','t_KernCld', 't_KernClr', 'wv_KernCld',  'wv_KernClr',  't_KernCld',  't_KernClr',  'Tq_KernCld',   'Tq_KernClr'     ]
        var1 =  ['tas',         'tas',      'ta',        'ta',        'ta',           'ta',         'tas4d',      'tas4d',      'tas4d',        'tas4d'      ]
        out1  = ['TS',         'TS_clr',    'T',         'T_clr',     'T_fxRH',      'T_fxRH_clr',  'Planck',     'Planck_clr', 'Planck_fxRH',  'Planck_fxRH_clr']
    
        kern1_in = [kern for kern in kern1]
        var1_in = [var+'_'+case for var in var1]
        out1_in = [out+'_'+case for out in out1]
        dic_rad = get_mult_latlon(dic_rad,kern1_in,var1_in,out1_in,dic_kern_wap,dic_mod_wap,dp4d)

              
        # add TS feedback to get total temperature and Planck feedbacks 
        dic_rad['T'+'_'+case] = dic_rad['T'+'_'+case] + dic_rad['TS'+'_'+case]
        dic_rad['T_clr'+'_'+case] = dic_rad['T_clr'+'_'+case] + dic_rad['TS_clr'+'_'+case]
        dic_rad['T_fxRH'+'_'+case] = dic_rad['T_fxRH'+'_'+case] + dic_rad['TS'+'_'+case]
        dic_rad['T_fxRH_clr'+'_'+case] = dic_rad['T_fxRH_clr'+'_'+case] + dic_rad['TS_clr'+'_'+case]
    
        dic_rad['Planck'+'_'+case] = dic_rad['Planck'+'_'+case] + dic_rad['TS'+'_'+case]
        dic_rad['Planck_clr'+'_'+case] = dic_rad['Planck_clr'+'_'+case] + dic_rad['TS_clr'+'_'+case]
        dic_rad['Planck_fxRH'+'_'+case] = dic_rad['Planck_fxRH'+'_'+case] + dic_rad['TS'+'_'+case]
        dic_rad['Planck_fxRH_clr'+'_'+case] = dic_rad['Planck_fxRH_clr'+'_'+case] + dic_rad['TS_clr'+'_'+case]

        #----------------------------------------------------------
        # Lapse rate Feedback
        #----------------------------------------------------------
        kern1 = ['t_KernCld', 't_KernClr', 'Tq_KernCld',  'Tq_KernClr' ]
        var1 =  ['dt',        'dt',         'dt',         'dt'         ]
        out1 =  ['LR',        'LR_clr',    'LR_fxRH',     'LR_fxRH_clr']

        kern1_in = [kern for kern in kern1]
        var1_in = [var+'_'+case for var in var1]
        out1_in = [out+'_'+case for out in out1]
        dic_rad = get_mult_latlon(dic_rad,kern1_in,var1_in,out1_in,dic_kern_wap,dic_mod_wap,dp4d)

        #----------------------------------------------------------
        # Albedo Feedback
        #----------------------------------------------------------
        kern1 = ['alb_KernCld', 'alb_KernClr']
        var1  = ['alb',         'alb'    ]
        out1 =  ['ALB',         'ALB_clr' ]
        
        kern1_in = [kern for kern in kern1]
        var1_in = [var+'_'+case for var in var1]
        out1_in = [out+'_'+case for out in out1]
        dic_rad = get_mult_latlon(dic_rad,kern1_in,var1_in,out1_in,dic_kern_wap,dic_mod_wap,dp4d)

        #----------------------------------------------------------
        # Water Vapor Feedback
        #----------------------------------------------------------
        kern1 = ['wv_lw_KernCld', 'wv_sw_KernCld', 'wv_lw_KernClr', 'wv_sw_KernClr', 'wv_lw_KernCld', 'wv_sw_KernCld', 'wv_lw_KernClr', 'wv_sw_KernClr' ]
        var1  = ['dlogq2',        'dlogq2',        'dlogq2',        'dlogq2',        'ta',            'ta',            'ta',            'ta'        ]
        out1 =  ['WV_lw',         'WV_sw',         'WV_lw_clr',     'WV_sw_clr',     'WV_lw_fxRH',    'WV_sw_fxRH',    'WV_lw_clr_fxRH','WV_sw_clr_fxRH']
        
        # note: fxRH LW/SW WV feedback: Q_kernel * ta warming anomaly
        kern1_in = [kern for kern in kern1]
        var1_in = [var+'_'+case for var in var1]
        out1_in = [out+'_'+case for out in out1]
        dic_rad = get_mult_latlon(dic_rad,kern1_in,var1_in,out1_in,dic_kern_wap,dic_mod_wap,dp4d)

        # ---------------------------------------------------------------#
        # get final RH feedback
        # ---------------------------------------------------------------#
        # final RH feedback related to RH change
        # RH feedback = default water vapor feedback - (water vapor kernel * uniform warming anomaly (Ts) - water vapor kernel * lapse rate temperature anomaly)
        # RH feedback = default water vapor feedback - water vapor kernel * atmospheric temperature change
        invar1 = ['WV_lw','WV_sw','WV_lw_clr','WV_sw_clr']
        outvar = ['netRH_lw','netRH_sw','netRH_lw_clr','netRH_sw_clr']
        
        for ivar,svar in enumerate(invar1):
            ovar = outvar[ivar]
            dic_rad[ovar+'_'+case] = dic_rad[svar+'_'+case] - dic_rad[svar+'_fxRH'+'_'+case]

        SUNDOWN = rsdt_ab
        for key in dic_rad.keys():
            if 'sw' in key or 'ALB' in key:
                logger.debug(f'Doing SUNDOWN for {key}')
                dic_rad[key] = dic_rad[key].where(SUNDOWN != 0.0, 0.0)

        print(dic_rad.keys())

        #----------------------------------------------------------
        # Adjusted cloud feedback     
        #----------------------------------------------------------
        # calculate cloud masking 
        invar1 = ['T_clr','WV_lw_clr','WV_sw_clr','ALB_clr']
        invar2 = ['T','WV_lw','WV_sw','ALB']
        
        for ivar,svar in enumerate(invar1):
            var1 = svar
            var2 = invar2[ivar]
            dic_rad[var2+'_mask'+'_'+case] = dic_rad[var1+'_'+case] - dic_rad[var2+'_'+case]
        
        dic_rad['SWCRE'+'_'+case] = dic_mod_wap['SWCRE'+'_'+case]
        dic_rad['LWCRE'+'_'+case] = dic_mod_wap['LWCRE'+'_'+case]
        dic_rad['netCRE'+'_'+case] = dic_mod_wap['netCRE'+'_'+case]
        dic_rad['rlut'+'_'+case] =   dic_mod_wap['rlut'+'_'+case]
        dic_rad['rsut'+'_'+case] =   dic_mod_wap['rsut'+'_'+case]
        dic_rad['rlutcs'+'_'+case] = dic_mod_wap['rlutcs'+'_'+case]
        dic_rad['rsutcs'+'_'+case] = dic_mod_wap['rsutcs'+'_'+case]

        # calculate adjusted CRE
        invar1 = ['T_mask','WV_sw_mask','SWCRE','LWCRE']
        invar2 = ['WV_lw_mask','ALB_mask','SW_adj','LW_adj']
        outvar = ['LW_adj','SW_adj','SWCRE_adj','LWCRE_adj']
        for ivar,svar in enumerate(invar1):
            var1 = svar
            var2 = invar2[ivar]
            ovar = outvar[ivar]
            dic_rad[ovar+'_'+case] = dic_rad[var1+'_'+case] + dic_rad[var2+'_'+case]

        dic_rad['net_adj'+'_'+case] = dic_rad['LW_adj'+'_'+case] + dic_rad['SW_adj'+'_'+case]
        dic_rad['netCRE_adj'+'_'+case] = dic_rad['netCRE'+'_'+case] + dic_rad['net_adj'+'_'+case]
 
        # calculate cloudy residual term 
        # get sum of kernel effect 
        dic_rad['LW_cld_sum'+'_'+case]  = dic_rad['T'+'_'+case] + dic_rad['WV_lw'+'_'+case] + dic_rad['LWCRE_adj'+'_'+case]
        dic_rad['SW_cld_sum'+'_'+case]  = dic_rad['ALB'+'_'+case] + dic_rad['WV_sw'+'_'+case] + dic_rad['SWCRE_adj'+'_'+case]
        dic_rad['net_cld_sum'+'_'+case] = dic_rad['LW_cld_sum'+'_'+case] + dic_rad['SW_cld_sum'+'_'+case]
                
        dic_rad['LW_clr_sum'+'_'+case]  = dic_rad['T_clr'+'_'+case] + dic_rad['WV_lw_clr'+'_'+case]
        dic_rad['SW_clr_sum'+'_'+case]  = dic_rad['ALB_clr'+'_'+case] + dic_rad['WV_sw_clr'+'_'+case]
        dic_rad['net_clr_sum'+'_'+case] = dic_rad['LW_clr_sum'+'_'+case] + dic_rad['SW_clr_sum'+'_'+case]
                    
        # get the TOA anomalies from direct model data
        dic_rad['LW_cld_dir'+'_'+case]  = -1. * dic_rad['rlut'+'_'+case]
        dic_rad['SW_cld_dir'+'_'+case]  = -1. * dic_rad['rsut'+'_'+case]
        dic_rad['net_cld_dir'+'_'+case] = dic_rad['LW_cld_dir'+'_'+case] + dic_rad['SW_cld_dir'+'_'+case]

        dic_rad['LW_clr_dir'+'_'+case]  = -1. * dic_rad['rlutcs'+'_'+case]
        dic_rad['SW_clr_dir'+'_'+case]  = -1. * dic_rad['rsutcs'+'_'+case]
        dic_rad['net_clr_dir'+'_'+case] = dic_rad['LW_clr_dir'+'_'+case] + dic_rad['SW_clr_dir'+'_'+case]
                    
        # difference between direct and kernel-calculation
        dic_rad['LW_cld_dir_sum'+'_'+case]  = dic_rad['LW_cld_dir'+'_'+case] - dic_rad['LW_cld_sum'+'_'+case]
        dic_rad['SW_cld_dir_sum'+'_'+case]  = dic_rad['SW_cld_dir'+'_'+case] - dic_rad['SW_cld_sum'+'_'+case]
        dic_rad['net_cld_dir_sum'+'_'+case] = dic_rad['LW_cld_dir_sum'+'_'+case] + dic_rad['SW_cld_dir_sum'+'_'+case]
                
        dic_rad['LW_clr_dir_sum'+'_'+case]  = dic_rad['LW_clr_dir'+'_'+case] - dic_rad['LW_clr_sum'+'_'+case]
        dic_rad['SW_clr_dir_sum'+'_'+case]  = dic_rad['SW_clr_dir'+'_'+case] - dic_rad['SW_clr_sum'+'_'+case]
        dic_rad['net_clr_dir_sum'+'_'+case] = dic_rad['LW_clr_dir_sum'+'_'+case] + dic_rad['SW_clr_dir_sum'+'_'+case]

    return dic_rad

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def get_mult_latlon(dic_rad, kern1,var1,out1,dic_kern_wap,dic_mod_wap,dp4d):
    '''
    kernel * xx_ab / xx_pi
    '''
    for ivar,svar in enumerate(kern1):
        a1 = dic_kern_wap[svar].to_masked_array() # [time,lat,lon]
        a2 = dic_mod_wap[var1[ivar]].to_masked_array() 
        #a1[a1.mask] = 0
        #a2[a2.mask] = 0
        outvar = out1[ivar]
        a3 = a1 * a2
        logger.debug(f'kernel={svar}, data={var1[ivar]}, {a1.shape}') 

        if len(a1.shape) == 3:
            dic_rad[outvar] = xr.DataArray(a3,coords=dic_kern_wap[svar].coords, dims=dic_kern_wap[svar].dims)
        elif len(a1.shape) == 4: #[time,lev,lat,lon]
            # vertical integral (sum)
            tmp = VertSum(a3, dp4d)
            dic_rad[outvar] = xr.DataArray(tmp,coords=dic_kern_wap[svar][:,0,:].coords, dims=dic_kern_wap[svar][:,0,:].dims)

    return dic_rad

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def cal_dlogq2_separate(dic_mod_wap):

    q1k =    dic_mod_wap['q1k_pi']
    hus_pi = dic_mod_wap['hus_pi']
    hus_ab = dic_mod_wap['hus_ab']

    dlogq1k = np.ma.log(q1k) - np.ma.log(hus_pi)
    dlogq   = np.ma.log(hus_ab) - np.ma.log(hus_pi)

    dlogq2_pi = np.ma.log(hus_pi)/dlogq1k
    dlogq2_ab = np.ma.log(hus_ab)/dlogq1k

    #dlogq2  = dlogq/dlogq1k
    #dlogq2_xr = xr.DataArray(dlogq2,coords=hus_pi.coords,dims=hus_pi.dims)

    dlogq2_pi_xr = xr.DataArray(dlogq2_pi,coords=hus_pi.coords,dims=hus_pi.dims)
    dlogq2_ab_xr = xr.DataArray(dlogq2_ab,coords=hus_pi.coords,dims=hus_pi.dims)

    return dlogq2_pi_xr, dlogq2_ab_xr


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def prints(varname,data):
    return print(varname,' = ',data.shape, np.nanmin(data),np.nanmax(data),type(data),'number of NaN =',np.sum(np.isnan(data)))
    #return print(varname,' = ',data.shape, np.min(data),np.max(data),type(data),'number of NaN =',np.sum(np.isnan(data)),'number of mask =',np.sum(data.mask))

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def read_data_e3sm(var,var2d,var3d,direc_data1,direc_data2,exp1,exp2,yrS_4d,monS_2d,yrE_4d,monE_2d,nyears):

    dic_invar = {}
    for svar in var:
        if svar in var:
            if not os.path.isfile(direc_data1+svar+'_'+exp1+'_'+yrS_4d+monS_2d+'-'+yrE_4d+monE_2d+'.nc'):
                if svar == 'ta':
                    svar_in = 'T'
                elif svar == 'hus':
                    svar_in = 'Q'
                else:
                    svar_in = svar
            else:
                svar_in = svar

            print(" =================== we are processing E3SM amip data", svar, " locally ====================")

            f1 = xr.open_dataset(direc_data1+svar_in+'_'+exp1+'_'+yrS_4d+monS_2d+'-'+yrE_4d+monE_2d+'.nc')
            dic_invar[svar+'_pi'] = f1[svar_in][:nyears*12,:,:]
            f1.close()
 
            f2 = xr.open_dataset(direc_data2+svar_in+'_'+exp2+'_'+yrS_4d+monS_2d+'-'+yrE_4d+monE_2d+'.nc')
            dic_invar[svar+'_ab'] = f2[svar_in][:nyears*12,:,:]
            f2.close()
 
            # reverse lev direction
            if svar in var3d:
                dic_invar[svar+'_pi'] = dic_invar[svar+'_pi'][:,::-1,:,:]
                dic_invar[svar+'_ab'] = dic_invar[svar+'_ab'][:,::-1,:,:]

            dic_invar[svar+'_ano'] = xr.DataArray(dic_invar[svar+'_ab'] - dic_invar[svar+'_pi'], coords=dic_invar[svar+'_pi'].coords)
            stop_here = False
        else:
            stop_here = True

        if stop_here: ### we don't get enough data to do further processing. just skip out this loop.
            print('stop_here is', stop_here)
            continue
   
    ### get SWCF, LWCF and their anomalies
    #dic_invar['SWCRE_ano']  = dic_invar['rsutcs_ano'] - dic_invar['rsut_ano']
    #dic_invar['LWCRE_ano']  = dic_invar['rlutcs_ano'] - dic_invar['rlut_ano']
    #dic_invar['netCRE_ano'] = dic_invar['SWCRE_ano']  + dic_invar['LWCRE_ano']

    dic_invar['SWCRE_pi']  = dic_invar['rsutcs_pi'] - dic_invar['rsut_pi']
    dic_invar['LWCRE_pi']  = dic_invar['rlutcs_pi'] - dic_invar['rlut_pi']
    dic_invar['netCRE_pi'] = dic_invar['SWCRE_pi']  + dic_invar['LWCRE_pi']
    dic_invar['SWCRE_ab']  = dic_invar['rsutcs_ab'] - dic_invar['rsut_ab']
    dic_invar['LWCRE_ab']  = dic_invar['rlutcs_ab'] - dic_invar['rlut_ab']
    dic_invar['netCRE_ab'] = dic_invar['SWCRE_ab']  + dic_invar['LWCRE_ab']


    ## get albedo
    dic_invar['alb_pi']  = np.true_divide(dic_invar['rsus_pi'],dic_invar['rsds_pi']) * 100.
    dic_invar['alb_ab']  = np.true_divide(dic_invar['rsus_ab'],dic_invar['rsds_ab']) * 100.
    # use np.ma.masked_invalid to mask nan although the warning is still there...
    #dic_mod['alb_pi']  = xr.DataArray(np.ma.masked_outside(np.ma.masked_invalid(dic_mod['alb_pi']), 0.0, 100.),coords=dic_mod['rsus_pi'].coords)
    #dic_mod['alb_ab']  = xr.DataArray(np.ma.masked_outside(np.ma.masked_invalid(dic_mod['alb_ab']), 0.0, 100.),coords=dic_mod['rsus_ab'].coords)
   
    dic_invar['alb_pi']  = xr.DataArray(np.ma.masked_outside(dic_invar['alb_pi'], 0.0, 100.),coords=dic_invar['rsus_pi'].coords)
    dic_invar['alb_ab']  = xr.DataArray(np.ma.masked_outside(dic_invar['alb_ab'], 0.0, 100.),coords=dic_invar['rsus_ab'].coords)

    delterm = ["rsus_pi","rsds_pi",\
               "rsus_ab","rsds_ab","ps_ab"]
    dic_invar = delete_vars(delterm,dic_invar)

    return dic_invar

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if __name__ == "__main__":

    RadKernel_dir = '/qfs/people/qiny108/diag_feedback_E3SM/Huang_kernel_data/'
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
    
    RadKernel(RadKernel_dir,direc_data,case_stamp,yearS,yearE,fname1,fname2,outdir,figdir,exp1,exp2)

