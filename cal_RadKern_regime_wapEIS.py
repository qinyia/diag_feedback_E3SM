#****************************************************************
#
#    Filename: cal_RadKern_omega_space_tropo4d_weighted_amip_2.5x2.5.py
#
#    Author: Yi Qin - qin4@llnl.gov
#    Description: the last version filename is cal_RadKern_omega_space_tropo4d_weighted_cmip_rawgrid_NewRegress.py
#                 copied from feedback2.llnl.gov  on 2022-05-04
#    Input: 
#    Output: 
#    Create: 2021-11-09 10:04:40
#    Last Modified: 2021-11-09 10:04:40
#    2021-11-09: change binedges from (binS,binE,width) to (binS,binE+width,width)
#                change horizontal resolution from 2x2.5 to 2.5x2.5 -- consistent with previous cdat version
#    2021-11-11: prescribe the start date of amip experiment. If not, some models' start date for amip might
#                be 1950 or something. This is not reasonable to compare with its amip4K from 1979.
#    2021-11-12: do not regrid model data to uniform grid at first. Regrid kernel data to model grid.
#    2022-03-09: output dCRE/dTlocal to get the contributors to amip vs cmip cloud feedback differences.
#    2022-03-24: calculate epoch difference -- output '_loc_epoch' variables
#    2022-03-25: replace genutil.statistics.linearregression by lag_linregress_3D() function. 
#    2022-05-05: modify it to fit to E3SM raw output data
#****************************************************************

import numpy as np
import pandas as pd
import netCDF4
import xarray as xr
from global_land_mask import globe
import os
import sys
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

sorting = True
sorting_method = 'Yi'
sorting_var = 'wapEIS' #'wap'

epoch_method = False

if sorting:
    tag = 'binraw'
else:
    tag = 'latlon'

do_hor_regrid = True

# define bins: 10 hPa/dy wide bins
width = 10 
binS = -100
binE = 100

binedges = np.arange(binS,binE+width,width) # range105~115 is implicitly included while sorting.
if sorting_method == 'Yi':
    bincenters = binedges
else:
    bincenters = np.arange(binS-width/2.,binE+width/2.+width,width)

### define uniform horizontal grids 
if do_hor_regrid:
    Dlats = np.arange(-90,92.5,2.5)
    Dlons = np.arange(1.25,360,2.5)

Dlevs = np.array([100000, 92500, 85000, 70000, 60000, 50000, 40000, 30000, 25000, 20000, 15000, 10000, 7000, 5000, 3000, 2000, 1000])/100.

kernel_dir = '/home/qin4/Data/Huang_kernel_data/'

### read necessary variables from models [month,(level),lat,lon]
var2d = ["tas","rlut","rsut","rlutcs","rsutcs","rsus","rsds","ps","rsdt","psl","ts"]
var3d = ["ta","hus","OMEGA","Z3"]

var = var2d + var3d

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def RadKernel_regime(kernel_dir,direc_data,case_stamp,yearS,yearE,fname1,fname2,outdir,figdir,exp1,exp2):

    outfile = "RadKern_regime_wapEIS_"+tag+"_"+case_stamp+".nc"

    if os.path.isfile(outdir+outfile):
        print('RadKenel is already there.')
        return 

    logger.remove()
    fmt = "<green>{time:YYYY-MM-DD HH:mm:ss}</green> | <cyan>{level}</cyan> | {message} |{elapsed}"
    logger.add(sys.stdout, format=fmt)
    logger.add("NewRegress",format=fmt)
    ###################

    yearS_4d = "{:04d}".format(yearS)
    yearE_4d = "{:04d}".format(yearE)
    nyears = yearE - yearS + 1

    direc_data1 = direc_data+'/'+fname1+'/'
    direc_data2 = direc_data+'/'+fname2+'/'

    monS = 1
    monE = 12
    monS_2d='{:02d}'.format(monS)
    monE_2d='{:02d}'.format(monE) 

    #=============================================================
    # read model's data
    #=============================================================
 
    dic_mod = read_data_e3sm(var,var2d,var3d,direc_data1,direc_data2,exp1,exp2,yearS_4d,monS_2d,yearE_4d,monE_2d,nyears)

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
            data_grd = data.interp(lat=Dlats, lon=Dlons)
        else:
            #<qinyi 2021-11-27 #------------------
            # HadGEM2-ES has different horizontal resolution for wap and other vars.
            if model in ['HadGEM2-ES','HadGEM2-A']:
                lats_here = dic_mod['tas_pi'].coords['lat'].values
                lons_here = dic_mod['tas_pi'].coords['lon'].values
                logger.info(f'Hello, we will do horizontally regrid for HadGEM2-ES specially.')
                data_grd = data.interp(lat=lats_here,lon=lons_here)
            else:
                data_grd = data
    
        if len(data.shape) == 3:
            dic_mod_fn[svar] = data_grd 
        elif len(data.shape) == 4:
            # vertical regrid model data to Dlevs
            #data_grd = data_grd.rename({'plev':'lev'}) # rename from plev to lev
            if data_grd.coords['lev'].max().values > 1300: # pressure in Pa
                logger.info(f'convert level from Pa to hPa')
                data_grd = data_grd.assign_coords({"lev":data_grd.coords['lev'].values/100})
            data_grd_vrt = data_grd.interp(lev=Dlevs)
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
        dtas_ann_grd = dtas_ann.interp(lat=Dlats,lon=Dlons)
    else:
        dtas_ann_grd = dtas_ann

    logger.info(f'dtas_avg.shape={dtas_avg.shape},{dtas_avg.min().values},{dtas_avg.max().values},dtas_avg={np.mean(dtas_avg.data)}')

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

    print('tas_ab_4d',tas_ab_4d[0,:,50,100].values)
    print('ta_ab',dic_mod_fn['ta_ab'][0,:,50,100].values)

    # mask tas_ano_grd_4d where ta_ano_vert_grd is True
    dic_mod_fn['tas4d_pi'] = tas_pi_4d.where((dic_mod_fn['ta_pi'].notnull()) & (tas_pi_4d.notnull()))
    dic_mod_fn['tas4d_ab'] = tas_ab_4d.where((dic_mod_fn['ta_ab'].notnull()) & (tas_ab_4d.notnull()))

    dic_mod_fn['ta_pi'] = dic_mod_fn['ta_pi'].where((dic_mod_fn['ta_pi'].notnull()) & (tas_pi_4d.notnull()))
    dic_mod_fn['ta_ab'] = dic_mod_fn['ta_ab'].where((dic_mod_fn['ta_ab'].notnull()) & (tas_ab_4d.notnull()))

    print('tas4d_ab',dic_mod_fn['tas4d_ab'][0,:,50,100].values)

    dic_mod_fn['dt_ab'] = dic_mod_fn['ta_ab'] - dic_mod_fn['tas4d_ab']
    dic_mod_fn['dt_pi'] = dic_mod_fn['ta_pi'] - dic_mod_fn['tas4d_pi']

    print('dt_ab',dic_mod_fn['dt_ab'][0,:,50,100].values)
    
    sub_name = 'Expand tas to verticals'
    print_memory_status(sub_name)

    #=============================================================
    # calculate q1k
    #=============================================================
    dic_mod_fn['q1k_pi'] = cal_q1k(dic_mod_fn)
    
    sub_name = 'Calculate q1k'
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
    print('dic_mod_tropo.keys()=',dic_mod_tropo.keys())

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
    # calculate dlogq2 
    #=============================================================
    dlogq2_pi,dlogq2_ab = cal_dlogq2_separate(dic_mod_tropo)
    dic_mod_tropo['dlogq2_pi'] = dlogq2_pi
    dic_mod_tropo['dlogq2_ab'] = dlogq2_ab

    sub_name = 'calculate dlogq2'
    print_memory_status(sub_name)
 
    #=============================================================
    # calculate Kernel * xx_ab and Kernel * xx_pi first 
    #=============================================================
    dic_rad = cal_Kern_rad_separate(dic_mod_tropo, dic_kern_tropo,dic_mod_tropo['dp4d_pi'],dic_mod_tropo['rsdt_ab'])
    print(dic_rad.keys())

    sub_name = 'calculate kernel*xx_ab and kernel*xx_pi first'
    print_memory_status(sub_name)
 
    #=============================================================
    # sort model data and kernel data into omega bins 
    #=============================================================
    if sorting:
        if sorting_var == 'wapEIS':
            dic_rad_wap = sort_var_wapEIS(dic_mod_tropo,dic_rad,case_stamp,outdir)
    else:
        dic_rad_wap = dic_rad

    sub_name = 'sort model and kernel data into omega bins'
    print_memory_status(sub_name)

    del dic_mod_tropo, dic_kern_tropo

    if sorting_var == 'wapEIS':
        #=============================================================
        # calculate dynamic, thermodynamic and co-variance
        #=============================================================
        for svar in dic_rad_wap.keys():
            dics_var_comp = {}
        
            if '_pi' in svar:
                svarh = '_'.join(svar.split('_')[:-1])
        
                svarh_pi = svar
                svarh_ab = svarh+'_ab'
        
                data1,data1_N = dic_rad_wap[svarh_pi]
                data2,data2_N = dic_rad_wap[svarh_ab]
        
                tot = data2*data2_N - data1*data1_N
                thermo = (data2-data1)*data1_N
                dyn = (data2_N-data1_N)*data1
                cov = (data2_N-data1_N)*(data2-data1)
                sum = thermo + dyn + cov 
        
                # here need to keep NaN consistent because calculation
                # of dyn term
                # only needs data1, which can be non-nan but data2 
                # can be nan. This will lead to the inconsistency b/t
                # therm and dyn terms.
                dyn = xr.where(np.isnan(thermo),np.nan,dyn)
        
                dics_var_comp[svarh+'_tot'] = tot/dtas_avg.mean()
                dics_var_comp[svarh+'_thermo'] = thermo/dtas_avg.mean()
                dics_var_comp[svarh+'_dyn'] = dyn/dtas_avg.mean()
                dics_var_comp[svarh+'_cov'] = cov/dtas_avg.mean()
                dics_var_comp[svarh+'_sum'] = sum/dtas_avg.mean()
        
                save_big_dataset(dics_var_comp,outdir+"/middata/Regime_Sorted_Components_"+svarh+"_"+case_stamp+".nc")
        
        logger.info(f'dics_var_comp.keys() = {dics_var_comp.keys()}')

    #=============================================================
    # calculate anomalies (xx_ab - xx_pi)
    #=============================================================
    dic_rad_wap_ano = get_anomalies(dic_rad_wap)

    sub_name = 'get anomalies in omega bins'
    print_memory_status(sub_name)

    print(dic_rad_wap_ano.keys())
   
    #=============================================================
    # regression against global-mean surface air temperature anomaly [yearly]
    #=============================================================
    dic_rad_perK = regress_tas('amip',dic_rad_wap_ano,dtas_avg,dtas_ann_grd,time_coord)

    sub_name = 'Regression onto tas'
    print_memory_status(sub_name)

    save_big_dataset(dic_rad_perK,outdir+outfile)

    return -1

    #=============================================================
    # get final outputs
    #=============================================================
    dic_final = get_outputs(dic_rad_perK)

    sub_name = 'Get outputs'
    print_memory_status(sub_name)
    
    #=============================================================
    # print regional/global values for test                    
    #=============================================================
    #prints_gm(dic_final)
    plt_figure(dic_final,model_out,figdir)

    sub_name = 'Plot figures'
    print_memory_status(sub_name)

    #=============================================================
    # output data into csv and NC files
    #=============================================================
    if sorting:
        output_file_bin(outdir,dic_final,model_out)
    else:
        output_file(outdir,dic_final,model_out)
        
    sub_name = 'Output data'
    print_memory_status(sub_name)

    return None

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def output_file_bin(outdir,dic_final,model_out):

    # output spatial 
    outname = ['fbk','avg']
    for iname,name in enumerate(outname):
        da = xr.Dataset(
        data_vars = {
        "T":               (('bin'),dic_final['T_'+name]),
        "Planck":          (('bin'),dic_final['Planck_'+name]),
        "LR":              (('bin'),dic_final['LR_'+name]),
        "WV":              (('bin'),dic_final['WV_'+name]),
        "ALB":             (('bin'),dic_final['ALB_'+name]),
        "LW_adj":          (('bin'),dic_final['LW_adj_'+name]),
        "SW_adj":          (('bin'),dic_final['SW_adj_'+name]),
        "net_adj":         (('bin'),dic_final['net_adj_'+name]),
        "SWCRE":           (('bin'),dic_final['SWCRE_'+name]),
        "LWCRE":           (('bin'),dic_final['LWCRE_'+name]),
        "netCRE":          (('bin'),dic_final['netCRE_'+name]),
        "SWCRE_adj":       (('bin'),dic_final['SWCRE_adj_'+name]),
        "LWCRE_adj":       (('bin'),dic_final['LWCRE_adj_'+name]),
        "netCRE_adj":      (('bin'),dic_final['netCRE_adj_'+name]),
        "SW_resd":         (('bin'),dic_final['SW_cld_dir_sum_'+name]),
        "LW_resd":         (('bin'),dic_final['LW_cld_dir_sum_'+name]),
        "net_resd":        (('bin'),dic_final['net_cld_dir_sum_'+name]),
        "T_clr":           (('bin'),dic_final['T_clr_'+name]),
        "Planck_clr":      (('bin'),dic_final['Planck_clr_'+name]), 
        "LR_clr":          (('bin'),dic_final['LR_clr_'+name]),
        "WV_clr":          (('bin'),dic_final['WV_clr_'+name]),
        "ALB_clr":         (('bin'),dic_final['ALB_clr_'+name]),
        "WV_SW":           (('bin'),dic_final['WV_sw_'+name]),
        "WV_LW":           (('bin'),dic_final['WV_lw_'+name]),
        "WV_clr_SW":       (('bin'),dic_final['WV_sw_clr_'+name]),
        "WV_clr_LW":       (('bin'),dic_final['WV_lw_clr_'+name]),
        "Planck_fxRH":     (('bin'),dic_final['Planck_fxRH_'+name]),
        "LR_fxRH":         (('bin'),dic_final['LR_fxRH_'+name]),
        "RH":              (('bin'),dic_final['netRH_'+name]),
        "Planck_clr_fxRH": (('bin'),dic_final['Planck_fxRH_clr_'+name]),
        "LR_clr_fxRH":     (('bin'),dic_final['LR_fxRH_clr_'+name]),
        "RH_clr":          (('bin'),dic_final['netRH_clr_'+name]),
        "LW_clr_sum":      (('bin'),dic_final['LW_clr_sum_'+name]),
        "SW_clr_sum":      (('bin'),dic_final['SW_clr_sum_'+name]),
        "net_clr_sum":     (('bin'),dic_final['net_clr_sum_'+name]),
        "LW_clr_dir":      (('bin'),dic_final['LW_clr_dir_'+name]),
        "SW_clr_dir":      (('bin'),dic_final['SW_clr_dir_'+name]),
        "net_clr_dir":     (('bin'),dic_final['net_clr_dir_'+name]),
        "LW_cld_sum":      (('bin'),dic_final['LW_cld_sum_'+name]),
        "SW_cld_sum":      (('bin'),dic_final['SW_cld_sum_'+name]),
        "net_cld_sum":     (('bin'),dic_final['net_cld_sum_'+name]),
        "LW_cld_dir":      (('bin'),dic_final['LW_cld_dir_'+name]),
        "SW_cld_dir":      (('bin'),dic_final['SW_cld_dir_'+name]),
        "net_cld_dir":     (('bin'),dic_final['net_cld_dir_'+name]),
        "LW_clr_resd":     (('bin'),dic_final['LW_clr_dir_sum_'+name]),
        "SW_clr_resd":     (('bin'),dic_final['SW_clr_dir_sum_'+name]),
        "net_clr_resd":    (('bin'),dic_final['net_clr_dir_sum_'+name]),
        },
        coords = {
        "bin": bincenters,
        }
        )
        
        da.to_netcdf(outdir+'latlon_fbk_'+tag+'_'+model_out+'_'+name+'.nc')

    varlist = ['T','Planck','LR','WV','ALB','LW_adj','SW_adj','net_adj','SWCRE','LWCRE','netCRE',
    'SWCRE_adj','LWCRE_adj','netCRE_adj','SW_cld_dir_sum','LW_cld_dir_sum','net_cld_dir_sum',
    'T_clr','Planck_clr','LR_clr','WV_clr','ALB_clr',
    'WV_sw','WV_lw','WV_sw_clr','WV_lw_clr','Planck_fxRH','LR_fxRH','netRH','Planck_fxRH_clr','LR_fxRH_clr','netRH_clr',
    'LW_clr_sum','SW_clr_sum','net_clr_sum','LW_clr_dir','SW_clr_dir','net_clr_dir','LW_cld_sum','SW_cld_sum','net_cld_sum',
    'LW_cld_dir','SW_cld_dir','net_cld_dir','LW_clr_dir_sum','SW_clr_dir_sum','net_clr_dir_sum']

    # output mean 
    with open(outdir+'gfbk_'+tag+'_'+model_out+'.csv', 'w') as f:
        f.write("%s,%s\n"%('var',model_out))

        for var in varlist:
            f.write("%s,%s\n"%(var,dic_final[var+'_gfbk']))

    return None

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def output_file(outdir,dic_final,model_out):

    # output spatial 
    outname = ['fbk']
    for iname,name in enumerate(outname):
        da = xr.Dataset(
        data_vars = {
        "T":               (('lat','lon'),dic_final['T_'+name]),
        "Planck":          (('lat','lon'),dic_final['Planck_'+name]),
        "LR":              (('lat','lon'),dic_final['LR_'+name]),
        "WV":              (('lat','lon'),dic_final['WV_'+name]),
        "ALB":             (('lat','lon'),dic_final['ALB_'+name]),
        "LW_adj":          (('lat','lon'),dic_final['LW_adj_'+name]),
        "SW_adj":          (('lat','lon'),dic_final['SW_adj_'+name]),
        "net_adj":         (('lat','lon'),dic_final['net_adj_'+name]),
        "SWCRE":           (('lat','lon'),dic_final['SWCRE_'+name]),
        "LWCRE":           (('lat','lon'),dic_final['LWCRE_'+name]),
        "netCRE":          (('lat','lon'),dic_final['netCRE_'+name]),
        "SWCRE_adj":       (('lat','lon'),dic_final['SWCRE_adj_'+name]),
        "LWCRE_adj":       (('lat','lon'),dic_final['LWCRE_adj_'+name]),
        "netCRE_adj":      (('lat','lon'),dic_final['netCRE_adj_'+name]),
        "SW_resd":         (('lat','lon'),dic_final['SW_cld_dir_sum_'+name]),
        "LW_resd":         (('lat','lon'),dic_final['LW_cld_dir_sum_'+name]),
        "net_resd":        (('lat','lon'),dic_final['net_cld_dir_sum_'+name]),
        "T_clr":           (('lat','lon'),dic_final['T_clr_'+name]),
        "Planck_clr":      (('lat','lon'),dic_final['Planck_clr_'+name]), 
        "LR_clr":          (('lat','lon'),dic_final['LR_clr_'+name]),
        "WV_clr":          (('lat','lon'),dic_final['WV_clr_'+name]),
        "ALB_clr":         (('lat','lon'),dic_final['ALB_clr_'+name]),
        "WV_SW":           (('lat','lon'),dic_final['WV_sw_'+name]),
        "WV_LW":           (('lat','lon'),dic_final['WV_lw_'+name]),
        "WV_clr_SW":       (('lat','lon'),dic_final['WV_sw_clr_'+name]),
        "WV_clr_LW":       (('lat','lon'),dic_final['WV_lw_clr_'+name]),
        "Planck_fxRH":     (('lat','lon'),dic_final['Planck_fxRH_'+name]),
        "LR_fxRH":         (('lat','lon'),dic_final['LR_fxRH_'+name]),
        "RH":              (('lat','lon'),dic_final['netRH_'+name]),
        "Planck_clr_fxRH": (('lat','lon'),dic_final['Planck_fxRH_clr_'+name]),
        "LR_clr_fxRH":     (('lat','lon'),dic_final['LR_fxRH_clr_'+name]),
        "RH_clr":          (('lat','lon'),dic_final['netRH_clr_'+name]),
        "LW_clr_sum":      (('lat','lon'),dic_final['LW_clr_sum_'+name]),
        "SW_clr_sum":      (('lat','lon'),dic_final['SW_clr_sum_'+name]),
        "net_clr_sum":     (('lat','lon'),dic_final['net_clr_sum_'+name]),
        "LW_clr_dir":      (('lat','lon'),dic_final['LW_clr_dir_'+name]),
        "SW_clr_dir":      (('lat','lon'),dic_final['SW_clr_dir_'+name]),
        "net_clr_dir":     (('lat','lon'),dic_final['net_clr_dir_'+name]),
        "LW_cld_sum":      (('lat','lon'),dic_final['LW_cld_sum_'+name]),
        "SW_cld_sum":      (('lat','lon'),dic_final['SW_cld_sum_'+name]),
        "net_cld_sum":     (('lat','lon'),dic_final['net_cld_sum_'+name]),
        "LW_cld_dir":      (('lat','lon'),dic_final['LW_cld_dir_'+name]),
        "SW_cld_dir":      (('lat','lon'),dic_final['SW_cld_dir_'+name]),
        "net_cld_dir":     (('lat','lon'),dic_final['net_cld_dir_'+name]),
        "LW_clr_resd":     (('lat','lon'),dic_final['LW_clr_dir_sum_'+name]),
        "SW_clr_resd":     (('lat','lon'),dic_final['SW_clr_dir_sum_'+name]),
        "net_clr_resd":    (('lat','lon'),dic_final['net_clr_dir_sum_'+name]),
        },
        coords = {
        "lat": dic_final['ALB_'+name].coords['lat'].values,
        "lon": dic_final['ALB_'+name].coords['lon'].values,
        }
        )
        
        da.to_netcdf(outdir+'latlon_fbk_'+tag+'_'+model_out+'.nc')

    varlist = ['T','Planck','LR','WV','ALB','LW_adj','SW_adj','net_adj','SWCRE','LWCRE','netCRE',
    'SWCRE_adj','LWCRE_adj','netCRE_adj','SW_cld_dir_sum','LW_cld_dir_sum','net_cld_dir_sum',
    'T_clr','Planck_clr','LR_clr','WV_clr','ALB_clr',
    'WV_sw','WV_lw','WV_sw_clr','WV_lw_clr','Planck_fxRH','LR_fxRH','netRH','Planck_fxRH_clr','LR_fxRH_clr','netRH_clr',
    'LW_clr_sum','SW_clr_sum','net_clr_sum','LW_clr_dir','SW_clr_dir','net_clr_dir','LW_cld_sum','SW_cld_sum','net_cld_sum',
    'LW_cld_dir','SW_cld_dir','net_cld_dir','LW_clr_dir_sum','SW_clr_dir_sum','net_clr_dir_sum']

    # output mean 
    with open(outdir+'gfbk_'+tag+'_'+model_out+'.csv', 'w') as f:
        f.write("%s,%s\n"%('var',model_out))

        for var in varlist:
            f.write("%s,%s\n"%(var,dic_final[var+'_gfbk']))

    return None

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def plt_figure(dic_final,model,figdir):

    # each variable
    fig = plt.figure(figsize=(12,9))
    ax = fig.add_subplot(1,1,1)

    varx = ['Planck','LR','WV','ALB','SWCRE_adj','LWCRE_adj','netCRE_adj','SWCRE','LWCRE','netCRE']
    colors = ['tab:blue','tab:cyan','tab:orange','tab:grey','blue','red','black','blue','red','black']
    for ivar,svar in enumerate(varx):
        svar1 = svar+'_fbk'
        if svar in ['SWCRE','LWCRE','netCRE']:
            ls = ':'
        else:
            ls = '-'

        ax.plot(bincenters, dic_final[svar1],marker='.',label=svar,ls=ls,color=colors[ivar])


    ax.legend()
    ax.axhline(y=0,ls=':',color='grey')

    ax.set_xlabel('omega [hPa/day]')
    ax.set_ylabel('Feedback [W/m2/K]')
    fig.savefig(figdir+'test-allfbks-'+model+'.png', dpi = 300)

    #==============================================================
    # clear-sky linearity test 
    var1 = ['LW_clr_sum','SW_clr_sum','net_clr_sum']
    var2 = ['LW_clr_dir','SW_clr_dir','net_clr_dir']
    var3 = ['LW_clr_dir_sum', 'SW_clr_dir_sum', 'net_clr_dir_sum']
    
    for ivar,svar in enumerate(var1):
        svar1 = var1[ivar]+'_fbk'
        svar2 = var2[ivar]+'_fbk'
        svar3 = var3[ivar]+'_fbk'

        fig = plt.figure(figsize=(12,9))
        ax = fig.add_subplot(1,1,1)

        data = np.array([dic_final[svar1],dic_final[svar2],dic_final[svar3]])
        labels = ['Kernel','Model','Model - Kernel']

        for ii,da in enumerate(data):
            if ii == 2:
                marker = '.'
                ls = ':'
            else:
                marker = '.'
                ls = '-'
            
            ax.plot(bincenters,da,marker=marker,ls=ls,
            label=svar1.split('_')[0]+': '+labels[ii])

            if ii == 2:
                for jj,binc in enumerate(bincenters):
                    ax.text(binc,da[jj],np.round(da[jj]/data[1][jj]*100,1))

        ax.legend()

        ax.axhline(y=0,ls=':',color='grey')

        ax.set_xlabel('omega [hPa/day]')
        ax.set_ylabel('Feedback [W/m2/K]')

        fig.savefig(figdir+'test-'+svar1.split('_')[0]+'-'+model+'.png', dpi = 300)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def prints_gm(dic_final):
    for svar in dic_final.keys():
        if 'gfbk' in svar:
            logger.debug(f'{svar}.mean = {dic_final[svar]}')

    print('------ Hi, summary all feedbacks except for cloud feedback------------')
    print('Planck feedback: ',              dic_final['Planck_gfbk'],'W/m2/K')
    print('Lapse rate feedback: ',          dic_final['LR_gfbk'],'W/m2/K')
    print('Lapse rate + Planck feedback: ', dic_final['LR_gfbk']+dic_final['Planck_gfbk'],'W/m2/K')
    print('Temperature feedback: ',         dic_final['T_gfbk'],'W/m2/K')
    print('Water vapor feedback: ',         dic_final['WV_gfbk'],'W/m2/K')
    print("Surface albedo feedback: ",      dic_final['ALB_gfbk'], "W/m2/K")
    
    print('fixedRH Planck feedback: ',              dic_final['Planck_fxRH_gfbk'],'W/m2/K')
    print('fixedRH Lapse rate feedback: ',          dic_final['LR_fxRH_gfbk'],'W/m2/K')
    print('fixedRH Lapse rate + Planck feedback: ', dic_final['LR_fxRH_gfbk']+dic_final['Planck_fxRH_gfbk'],'W/m2/K')
    print('fixedRH Water vapor feedback: ',         dic_final['WV_fxRH_gfbk'],'W/m2/K')
    print('fixedRH RH feedback: ',                  dic_final['netRH_gfbk'],'W/m2/K')
    
    print('--------clear-sky component----------------------------------------------')
    print('clr Planck feedback: ',                  dic_final['Planck_clr_gfbk'],'W/m2/K')
    print('clr Lapse rate feedback: ',              dic_final['LR_clr_gfbk'],'W/m2/K')
    print('clr Lapse rate + clr Planck feedback: ', dic_final['LR_clr_gfbk']+dic_final['Planck_clr_gfbk'],'W/m2/K')
    print('clr Temperature feedback: ',             dic_final['T_clr_gfbk'],'W/m2/K')
    print('clr Water vapor feedback: ',             dic_final['WV_clr_gfbk'],'W/m2/K')
    print("clr Surface albedo feedback: ",          dic_final['ALB_clr_gfbk'],"W/m2/K")
    
    print('fixedRH clr Planck feedback: ',              dic_final['Planck_fxRH_clr_gfbk'],'W/m2/K')
    print('fixedRH clr Lapse rate feedback: ',          dic_final['LR_fxRH_clr_gfbk'],'W/m2/K')
    print('fixedRH clr Lapse rate + Planck feedback: ', dic_final['LR_fxRH_clr_gfbk']+dic_final['Planck_fxRH_clr_gfbk'],'W/m2/K')
    print('fixedRH clr Water vapor feedback: ',         dic_final['WV_lw_clr_fxRH_gfbk']+dic_final['WV_sw_clr_fxRH_gfbk'],'W/m2/K')
    print('fixedRH clr RH feedback: ',                  dic_final['netRH_clr_gfbk'],'W/m2/K')

    print('--------------- quick check here --------------------------')
    print('fixedRH (Lapse rate + Planck) - (Lapse rate + Planck): ',dic_final['LR_fxRH_gfbk'] + dic_final['Planck_fxRH_gfbk'] - dic_final['LR_gfbk'] - dic_final['Planck_gfbk'],'W/m2/K')
    print('fixedRH WV: ',dic_final['WV_fxRH_gfbk'],'W/m2/K')
    print('fixedRH Ta: ',dic_final['T_fxRH_gfbk'],'W/m2/K')
    
    print('clr fixedRH (Lapse rate + Planck) - (Lapse rate + Planck): ',dic_final['LR_fxRH_clr_gfbk'] + dic_final['Planck_fxRH_clr_gfbk'] - dic_final['LR_clr_gfbk'] - dic_final['Planck_clr_gfbk'],'W/m2/K')
    print('clr fixedRH WV: ',dic_final['WV_lw_clr_fxRH_gfbk']+dic_final['WV_sw_clr_fxRH_gfbk'],'W/m2/K')
    print('clr fixedRH Ta: ',dic_final['T_fxRH_clr_gfbk'],'W/m2/K')
    
    print('Lapse rate + Planck + WV feedback: ',dic_final['LR_gfbk'] + dic_final['Planck_gfbk'] + dic_final['WV_gfbk'],'W/m2/K')
    print('fixedRH Lapse rate + Planck + WV feedback: ',dic_final['LR_fxRH_gfbk'] + dic_final['Planck_fxRH_gfbk'] + dic_final['netRH_gfbk'],'W/m2/K')
    
    print('clr Lapse rate + Planck + WV feedback: ',dic_final['LR_clr_gfbk'] + dic_final['Planck_clr_gfbk'] + dic_final['WV_clr_gfbk'],'W/m2/K')
    print('fixedRH clr Lapse rate + Planck + WV feedback: ',dic_final['LR_fxRH_clr_gfbk'] + dic_final['Planck_fxRH_clr_gfbk'] + dic_final['netRH_clr_gfbk'],'W/m2/K')
    
    print('sum of all clear-sky sw feedbacks',dic_final['SW_clr_sum_gfbk'],'W/m2/K')
    print('sum of all clear-sky lw feedbacks',dic_final['LW_clr_sum_gfbk'],'W/m2/K')
    
    print('TOA direct clear-sky sw radiation feedback',dic_final['SW_clr_dir_gfbk'],'W/m2/K')
    print('TOA direct clear-sky lw radiation feedback',dic_final['LW_clr_dir_gfbk'],'W/m2/K')

    print('TOA direct clear-ksy sw/sum of all clear-sky sw',(dic_final['SW_clr_sum_gfbk']-dic_final['SW_clr_dir_gfbk'])/dic_final['SW_clr_dir_gfbk']*100., '%')
    print('TOA direct clear-ksy lw/sum of all clear-sky lw',(dic_final['LW_clr_sum_gfbk']-dic_final['LW_clr_dir_gfbk'])/dic_final['LW_clr_dir_gfbk']*100., '%')

    return None
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def get_outputs(dic_rad_perK):
    dic_final = {}

    outname = ['gfbk','fbk','avg']
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
        if exp_ctl == 'amip':
            slope = data_ann.mean(axis=0)/dtas_avg.mean()
        else:
            _,_,slope,intercept,_,_ = lag_linregress_3D(dtas_avg, data_ann, lagx=0, lagy=0)

        # save
        dic_rad_perK[svar+'_fbk'] = xr.DataArray(slope,coords=data_ann[0,:].coords, dims=data_ann[0,:].dims)

        # get bin average and global average
        if sorting:
            if sorting_var == 'wapEIS':
                tmp = np.ma.sum(slope[:-1]) # ignore the last one for 'glob' if sorting_var = 'wapEIS'
            else:
                tmp = np.ma.sum(slope)
        else:
            latS = -90.
            latE = 90.
            tmp = area_averager(slope)

        # save
        dic_rad_perK[svar+'_gfbk'] = tmp

        print(svar, tmp)

        # also save time-averaged anomaly in each bin 
        tmp1 = data_ann.mean(axis=0)
        dic_rad_perK[svar+'_avg'] = xr.DataArray(tmp1,coords=data_ann[0,:].coords, dims=data_ann[0,:].dims)

        if not sorting and epoch_method: # on lat-lon space
            # ============================================================================
            # 2022-03-09: regress on local tas anomaly rather than global mean tas anomaly
            # get annual mean for tas
            if exp_ctl == 'amip':
                slope = data_ann.mean(axis=0)/dtas_ann.mean(axis=0)
            else:
                _,_,slope,intercept,_,_ = lag_linregress_3D(dtas_ann, data_ann, lagx=0, lagy=0)
            # save
            dic_rad_perK[svar+'_fbk_loc'] = xr.DataArray(slope,coords=data_ann[0,:].coords, dims=data_ann[0,:].dims)
    
            ntt = 20
            # 2022-03-24: get epoch difference (last 20yr minus first 20yr) for coupled local feedback.
            if exp_ctl != 'amip':
                da1 = data_ann[-ntt:,:].mean(axis=0)
                da2 = data_ann[:ntt,:].mean(axis=0)
                t1 = dtas_avg[-ntt:].mean()
                t2 = dtas_avg[:ntt].mean()
                slope = (da1-da2)/(t1-t2)
                # save
                dic_rad_perK[svar+'_fbk_epoch'] = xr.DataArray(slope,coords=data_ann[0,:].coords, dims=data_ann[0,:].dims)
     
                da1 = data_ann[-ntt:,:].mean(axis=0)
                da2 = data_ann[:ntt,:].mean(axis=0)
                t1 = dtas_ann[-ntt:,:].mean(axis=0)
                t2 = dtas_ann[:ntt,:].mean(axis=0)
                slope = (da1-da2)/(t1-t2)
                # save
                dic_rad_perK[svar+'_fbk_loc_epoch'] = xr.DataArray(slope,coords=data_ann[0,:].coords, dims=data_ann[0,:].dims)
            #=============================================================
    
            # ============================================================================
            # 2022-03-10: regress lat-lon tas anomaly on global mean tas anomaly
            if exp_ctl =='amip':
                slope = dtas_ann.mean(axis=0)/dtas_avg.mean()
            else:
                _,_,slope,intercept,_,_ = lag_linregress_3D(dtas_avg, dtas_ann, lagx=0, lagy=0)
    
            # save
            dic_rad_perK['tas_fbk'] = xr.DataArray(slope,coords=dtas_ann[0,:].coords, dims=dtas_ann[0,:].dims)
    
            # 2022-03-24: get tas anomaly from epoch 
            if exp_ctl !='amip':
                da1 = dtas_ann[-ntt:,:].mean(axis=0)
                da2 = dtas_ann[:ntt,:].mean(axis=0)
                t1 = dtas_avg[-ntt:].mean()
                t2 = dtas_avg[:ntt].mean()
                slope = (da1-da2)/(t1-t2)
                # save 
                dic_rad_perK['tas_fbk_epoch'] = xr.DataArray(slope,coords=dtas_ann[0,:].coords, dims=dtas_ann[0,:].dims)
            # ============================================================================
   
    return dic_rad_perK

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def cal_Kern_rad(dic_mod_wap_ano, dic_kern_wap,dp4d,rsdt_ab):

    dic_rad = {}
    #----------------------------------------------------------
    # Temperature Feedback
    #----------------------------------------------------------
    kern1 = ['ts_KernCld', 'ts_KernClr','t_KernCld', 't_KernClr', 'wv_KernCld',  'wv_KernClr',  't_KernCld',  't_KernClr',  'Tq_KernCld',   'Tq_KernClr'     ]
    var1 =  ['tas_ano',    'tas_ano',   'ta_ano',    'ta_ano',    'ta_ano',      'ta_ano',      'tas4d_ano',  'tas4d_ano',  'tas4d_ano',    'tas4d_ano'      ]
    out1  = ['TS',         'TS_clr',    'T',         'T_clr',     'T_fxRH',      'T_fxRH_clr',  'Planck',     'Planck_clr', 'Planck_fxRH',  'Planck_fxRH_clr']

    kern1_in = [kern for kern in kern1]
    dic_rad = get_mult(dic_rad,kern1_in,var1,out1,dic_kern_wap,dic_mod_wap_ano,dp4d)
           
    # add TS feedback to get total temperature and Planck feedbacks 
    dic_rad['T'] = dic_rad['T'] + dic_rad['TS']
    dic_rad['T_clr'] = dic_rad['T_clr'] + dic_rad['TS_clr']
    dic_rad['T_fxRH'] = dic_rad['T_fxRH'] + dic_rad['TS']
    dic_rad['T_fxRH_clr'] = dic_rad['T_fxRH_clr'] + dic_rad['TS_clr']

    dic_rad['Planck'] = dic_rad['Planck'] + dic_rad['TS']
    dic_rad['Planck_clr'] = dic_rad['Planck_clr'] + dic_rad['TS_clr']
    dic_rad['Planck_fxRH'] = dic_rad['Planck_fxRH'] + dic_rad['TS']
    dic_rad['Planck_fxRH_clr'] = dic_rad['Planck_fxRH_clr'] + dic_rad['TS_clr']

    #----------------------------------------------------------
    # Lapse rate Feedback
    #----------------------------------------------------------
    kern1 = ['t_KernCld', 't_KernClr', 'Tq_KernCld',  'Tq_KernClr'     ]
    var1 =  ['dt_ano',    'dt_ano',    'dt_ano',      'dt_ano'         ]
    out1 =  ['LR',        'LR_clr',    'LR_fxRH',     'LR_fxRH_clr']

    kern1_in = [kern for kern in kern1]
    dic_rad = get_mult(dic_rad,kern1_in,var1,out1,dic_kern_wap,dic_mod_wap_ano,dp4d)
    
    #----------------------------------------------------------
    # Albedo Feedback
    #----------------------------------------------------------
    kern1 = ['alb_KernCld', 'alb_KernClr']
    var1  = ['alb_ano',     'alb_ano'    ]
    out1 =  ['ALB',         'ALB_clr' ]
    
    kern1_in = [kern for kern in kern1]
    dic_rad = get_mult(dic_rad,kern1_in,var1,out1,dic_kern_wap,dic_mod_wap_ano,dp4d)

    #----------------------------------------------------------
    # Water Vapor Feedback
    #----------------------------------------------------------
    kern1 = ['wv_lw_KernCld', 'wv_sw_KernCld', 'wv_lw_KernClr', 'wv_sw_KernClr', 'wv_lw_KernCld', 'wv_sw_KernCld', 'wv_lw_KernClr', 'wv_sw_KernClr' ]
    var1  = ['dlogq2',        'dlogq2',        'dlogq2',        'dlogq2',        'ta_ano',        'ta_ano',        'ta_ano',        'ta_ano'        ]
    out1 =  ['WV_lw',         'WV_sw',         'WV_lw_clr',     'WV_sw_clr',     'WV_lw_fxRH',    'WV_sw_fxRH',    'WV_lw_clr_fxRH','WV_sw_clr_fxRH']
    
    # note: fxRH LW/SW WV feedback: Q_kernel * ta warming anomaly
    kern1_in = [kern for kern in kern1]
    dic_rad = get_mult(dic_rad,kern1_in,var1,out1,dic_kern_wap,dic_mod_wap_ano,dp4d)

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
        dic_rad[ovar] = dic_rad[svar] - dic_rad[svar+'_fxRH']

    # SUNDOWN is not necessary because we don't consider regions outside of tropics. 
    if not sorting:
        SUNDOWN = rsdt_ab
        for key in dic_rad.keys():
            if 'sw' in key or 'ALB' in key:
                logger.debug(f'Doing SUNDOWN for {key}')
                dic_rad[key] = dic_rad[key].where(SUNDOWN != 0.0, 0.0)

    #----------------------------------------------------------
    # Adjusted cloud feedback     
    #----------------------------------------------------------
    # calculate cloud masking 
    invar1 = ['T_clr','WV_lw_clr','WV_sw_clr','ALB_clr']
    invar2 = ['T','WV_lw','WV_sw','ALB']
    
    for ivar,svar in enumerate(invar1):
        var1 = svar
        var2 = invar2[ivar]
        dic_rad[var2+'_mask'] = dic_rad[var1] - dic_rad[var2]
    
    dic_rad['SWCRE'] = dic_mod_wap_ano['SWCRE_ano']
    dic_rad['LWCRE'] = dic_mod_wap_ano['LWCRE_ano']
    dic_rad['netCRE'] = dic_mod_wap_ano['netCRE_ano']
    dic_rad['rlut'] =   dic_mod_wap_ano['rlut_ano']
    dic_rad['rsut'] =   dic_mod_wap_ano['rsut_ano']
    dic_rad['rlutcs'] = dic_mod_wap_ano['rlutcs_ano']
    dic_rad['rsutcs'] = dic_mod_wap_ano['rsutcs_ano']

    # calculate adjusted CRE
    invar1 = ['T_mask','WV_sw_mask','SWCRE','LWCRE']
    invar2 = ['WV_lw_mask','ALB_mask','SW_adj','LW_adj']
    outvar = ['LW_adj','SW_adj','SWCRE_adj','LWCRE_adj']
    for ivar,svar in enumerate(invar1):
        var1 = svar
        var2 = invar2[ivar]
        ovar = outvar[ivar]
        dic_rad[ovar] = dic_rad[var1] + dic_rad[var2]

    dic_rad['net_adj'] = dic_rad['LW_adj'] + dic_rad['SW_adj']
    dic_rad['netCRE_adj'] = dic_rad['netCRE'] + dic_rad['net_adj']
 
    # calculate cloudy residual term 
    # get sum of kernel effect 
    dic_rad['LW_cld_sum']  = dic_rad['T'] + dic_rad['WV_lw'] + dic_rad['LWCRE_adj']
    dic_rad['SW_cld_sum']  = dic_rad['ALB'] + dic_rad['WV_sw'] + dic_rad['SWCRE_adj']
    dic_rad['net_cld_sum'] = dic_rad['LW_cld_sum'] + dic_rad['SW_cld_sum']
            
    dic_rad['LW_clr_sum']  = dic_rad['T_clr'] + dic_rad['WV_lw_clr']
    dic_rad['SW_clr_sum']  = dic_rad['ALB_clr'] + dic_rad['WV_sw_clr']
    dic_rad['net_clr_sum'] = dic_rad['LW_clr_sum'] + dic_rad['SW_clr_sum']
                
    # get the TOA anomalies from direct model data
    dic_rad['LW_cld_dir']  = -1. * dic_rad['rlut']
    dic_rad['SW_cld_dir']  = -1. * dic_rad['rsut']
    dic_rad['net_cld_dir'] = dic_rad['LW_cld_dir'] + dic_rad['SW_cld_dir']

    dic_rad['LW_clr_dir']  = -1. * dic_rad['rlutcs']
    dic_rad['SW_clr_dir']  = -1. * dic_rad['rsutcs']
    dic_rad['net_clr_dir'] = dic_rad['LW_clr_dir'] + dic_rad['SW_clr_dir']
                
    # difference between direct and kernel-calculation
    dic_rad['LW_cld_dir_sum']  = dic_rad['LW_cld_dir'] - dic_rad['LW_cld_sum']
    dic_rad['SW_cld_dir_sum']  = dic_rad['SW_cld_dir'] - dic_rad['SW_cld_sum']
    dic_rad['net_cld_dir_sum'] = dic_rad['LW_cld_dir_sum'] + dic_rad['SW_cld_dir_sum']
            
    dic_rad['LW_clr_dir_sum']  = dic_rad['LW_clr_dir'] - dic_rad['LW_clr_sum']
    dic_rad['SW_clr_dir_sum']  = dic_rad['SW_clr_dir'] - dic_rad['SW_clr_sum']
    dic_rad['net_clr_dir_sum'] = dic_rad['LW_clr_dir_sum'] + dic_rad['SW_clr_dir_sum']

    return dic_rad

#=============================================================
def cal_q1k(dic_mod):
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
def cal_dlogq2_bin(dic_mod_wap):

    ta_pi = dic_mod_wap['ta_pi']# * dic_mod_wap['ta_pi_N']
    ta_ab = dic_mod_wap['ta_ab']# * dic_mod_wap['ta_ab_N']
    avgta = (ta_pi + ta_ab)/2.0
    #avgta = ta_pi

    hus_pi = dic_mod_wap['hus_pi']# *dic_mod_wap['hus_pi_N']
    hus_ab = dic_mod_wap['hus_ab']# *dic_mod_wap['hus_ab_N']

    levs = dic_mod_wap['ta_pi'].coords['lev'].values
    ntime,nlev,nbin = dic_mod_wap['ta_pi'].shape
    
    levnd = np.transpose(np.tile(levs,(ntime,nbin,1)),(0,2,1))

    qs0,qs1 = qsat_blend_Mark(avgta,levnd)
    rh0 = hus_pi/qs0
    q1k = rh0*qs1 

    dlogq1k = np.ma.log(q1k) - np.ma.log(hus_pi)
    dlogq   = np.ma.log(hus_ab) - np.ma.log(hus_pi)
    dlogq2  = dlogq/dlogq1k

    dif1 = q1k - hus_pi
    dif2 = hus_ab - hus_pi

    dic_tmp = {}
    dic_tmp['qs0'] = qs0
    dic_tmp['qs1'] = qs1
    dic_tmp['rh0'] = rh0
    dic_tmp['q1k'] = q1k
    dic_tmp['dif1'] = dif1
    dic_tmp['dif2'] = dif2
    dic_tmp['hus_pi'] = hus_pi
    dic_tmp['hus_ab'] = hus_ab
    dic_tmp['q1k'] = q1k
    dic_tmp['hus_pi'] = hus_pi
    dic_tmp['hus_ab'] = hus_ab
    dic_tmp['dlogq1k'] = dlogq1k
    dic_tmp['dlogq'] = dlogq
    dic_tmp['dlogq2'] = dlogq2

    save_big_dataset(dic_tmp,"test-q1k-dlogq2-bin.nc")

    dlogq2_xr = xr.DataArray(dlogq2,coords=hus_pi.coords,dims=hus_pi.dims)

    return dlogq2_xr

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def cal_dlogq2(dic_mod_wap):

    q1k =    dic_mod_wap['q1k_pi']# *dic_mod_wap['q1k_pi_N']
    hus_pi = dic_mod_wap['hus_pi']# *dic_mod_wap['hus_pi_N']
    hus_ab = dic_mod_wap['hus_ab']# *dic_mod_wap['hus_ab_N']

    dlogq1k = np.ma.log(q1k) - np.ma.log(hus_pi)
    dlogq   = np.ma.log(hus_ab) - np.ma.log(hus_pi)
    dlogq2  = dlogq/dlogq1k

    dlogq2_xr = xr.DataArray(dlogq2,coords=hus_pi.coords,dims=hus_pi.dims)

    return dlogq2_xr

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def cal_dlogq2_wts(dic_mod_wap):

    q1k =    dic_mod_wap['q1k_pi']# *dic_mod_wap['hus_pi_N']
    hus_pi = dic_mod_wap['hus_pi']# *dic_mod_wap['hus_pi_N']
    hus_ab = dic_mod_wap['hus_ab']# *dic_mod_wap['hus_ab_N']

    dlogq1k = np.ma.log(q1k)*dic_mod_wap['hus_pi_N'] - np.ma.log(hus_pi)*dic_mod_wap['hus_pi_N']
    dlogq   = np.ma.log(hus_ab)*dic_mod_wap['hus_ab_N'] - np.ma.log(hus_pi)*dic_mod_wap['hus_pi_N']

    dlogq2  = dlogq/dlogq1k

    dic_tmp = {}
    dic_tmp['q1k'] = q1k
    dic_tmp['hus_pi'] = hus_pi
    dic_tmp['hus_ab'] = hus_ab
    dic_tmp['dlogq1k'] = xr.DataArray(dlogq1k,coords=hus_pi.coords,dims=hus_pi.dims)
    dic_tmp['dlogq'] = xr.DataArray(dlogq,coords=hus_pi.coords,dims=hus_pi.dims)
    dic_tmp['dlogq2'] = xr.DataArray(dlogq2,coords=hus_pi.coords,dims=hus_pi.dims)

    save_big_dataset(dic_tmp,"test-dlogq2-"+tag+".nc")

    dlogq2_xr = xr.DataArray(dlogq2,coords=hus_pi.coords,dims=hus_pi.dims)

    return dlogq2_xr

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def get_anomalies(dic_mod_wap):
    dic_mod_ano = {}
    for svar in dic_mod_wap.keys():
        if '_pi' in svar:
            var1 = '_'.join(svar.split('_')[:-1])
            if var1 not in ['ps','q1k','dp4d']: # ignore ps because only read ps_pi
                da1,da1_N = dic_mod_wap[var1+'_pi']
                da2,da2_N = dic_mod_wap[var1+'_ab']
                diff = da2*da2_N - da1*da1_N
                dic_mod_ano[var1+'_ano'] = diff 

    ## temperature deviation from vertical uniform warming 
    #dt_ano = dic_mod_ano['ta_ano'] - dic_mod_ano['tas4d_ano']
    #dt_ano  = np.ma.masked_where(dic_mod_ano['ta_ano'].to_masked_array().mask, dt_ano)
    #dic_mod_ano['dt_ano'] = xr.DataArray(dt_ano, coords=dic_mod_ano['ta_ano'].coords, dims=dic_mod_ano['ta_ano'].dims)

    return dic_mod_ano

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def apply_mask_tropo(Dlevs,dic_mod_fn,dic_kern,tropo4d,coords_4d,dims_4d):
    '''
    mask data above tropopause 
    ??? Should I set mask as zero here????
    '''
    ntime,nlev,nlat,nlon = dic_mod_fn['ta_pi'].shape
    lev4d = np.transpose(np.tile(Dlevs,(ntime,nlat,nlon,1)),(0,3,1,2)) 
    logger.debug(f'lev4d = {lev4d.shape}, {np.nanmin(lev4d)}, {np.nanmax(lev4d)}')
    logger.debug(f'tropo4d = {tropo4d.shape}, {np.nanmin(tropo4d)}, {np.nanmax(tropo4d)}')

    dic_mod_tropo = {}
    #for svar in dic_mod_fn.keys():
    #    if len(dic_mod_fn[svar].shape) == 4: # [time,lev,lat,lon]
    #        dic_mod_tropo[svar] = np.ma.masked_where(lev4d<=tropo4d, dic_mod_fn[svar])
    #        #dic_mod_tropo[svar][dic_mod_tropo[svar].mask] = 0.0
    #        dic_mod_tropo[svar] = xr.DataArray(dic_mod_tropo[svar], coords=coords_4d, dims=dims_4d)
    #    else:
    #        dic_mod_tropo[svar] = dic_mod_fn[svar]
    dic_mod_tropo = dic_mod_fn


    dic_kern_tropo = {}
    for svar in dic_kern.keys():
        if len(dic_kern[svar].shape) == 4: 
            dic_kern_tropo[svar] = np.ma.masked_where(lev4d<=tropo4d, dic_kern[svar])
            #dic_kern_tropo[svar][dic_kern_tropo[svar].mask] = 0.0
            dic_kern_tropo[svar] = xr.DataArray(dic_kern_tropo[svar], coords=coords_4d, dims=dims_4d)
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
    logger.debug(f'PP={PP}')

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
def sort_Kern_data_yi(dic_mod,dic_kern):
    '''
    sort kernel data and variables 
    output: model and kernel data in wap bins: [time,(lev),nbins]
    '''
    # get lats and lons
    ctl_wap = dic_mod['OMEGA_pi'].sel(lev=500,method='nearest')
    fut_wap = dic_mod['OMEGA_ab'].sel(lev=500,method='nearest')

    # Make sure wap is in hPa/day
    ctl_wap = 36*24*ctl_wap # Pa/s -> hPa/day
    fut_wap = 36*24*fut_wap 

    lats = ctl_wap.coords['lat']
    lons = ctl_wap.coords['lon']

    # define tropical regions
    latS = -30; latE = 30
    lats_tr = lats.sel(lat=slice(latS,latE))

    ## mask land and only select tropical data
    ctl_wap_ocean,_ = mask_land(lons,lats_tr,ctl_wap.sel(lat=slice(latS,latE)))
    fut_wap_ocean,_ = mask_land(lons,lats_tr,fut_wap.sel(lat=slice(latS,latE)))

    # sort model variables into bins
    dic_mod_wap = {}
    for svar in dic_mod.keys():
        logger.debug(f'=======>>>>>>> svar ={svar} is doing bony_sorting')
        data,_ = mask_land(lons,lats_tr, dic_mod[svar].sel(lat=slice(latS,latE)))
        if '_pi' in svar:
            omega = ctl_wap_ocean
        else:
            omega = fut_wap_ocean

        DATA = sort_var_by_omega(width,binedges,omega,data)

        # define coordinates in bins
        if len(DATA.shape) == 2:
            coords_bin = {"time":ctl_wap.coords['time'].values, "bin": ("bin",bincenters)}
            dims_bin   = ['time','bin']
            DATA_out = DATA[:,:]
        else:
            coords_bin = {"time":ctl_wap.coords['time'].values, "lev":data.coords['lev'].values,"bin": ("bin",bincenters)}
            dims_bin   = ['time','lev','bin']
            DATA_out = DATA[:,:,:]

        DATA_xr = xr.DataArray(DATA_out,coords=coords_bin, dims=dims_bin) # ignore land
        dic_mod_wap[svar] = DATA_xr

    # sort kernel data into bins
    dic_kern_wap = {}
    for svar in dic_kern.keys():
        logger.debug(f'=======>>>>>>> svar ={svar} is doing bony_sorting_part2')
        data,_ = mask_land(lons,lats_tr, dic_kern[svar].sel(lat=slice(latS,latE)))
        omega = ctl_wap_ocean

        DATA = sort_var_by_omega(width,binedges,omega,data)

        # define coordinates in bins
        if len(DATA.shape) == 2:
            coords_bin = {"time":ctl_wap.coords['time'], "bin": ("bin",bincenters)}
            dims_bin   = ['time','bin']
            DATA_out = DATA[:,:]
        else:
            coords_bin = {"time":ctl_wap.coords['time'], "lev":data.coords['lev'],"bin": ("bin",bincenters)}
            dims_bin   = ['time','lev','bin']
            DATA_out = DATA[:,:,:]

        DATA_xr = xr.DataArray(DATA_out,coords=coords_bin, dims=dims_bin) # ignore land
        dic_kern_wap[svar] = DATA_xr
       
    return dic_mod_wap, dic_kern_wap

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
                data_grd = data[:,::-1,:].interp(lat=model_lat,lon=model_lon)
                # rename dimension name
                data_grd = data_grd.rename({'month':'time'})
                dic_kern[varout] = data_grd

            else: #(time,lev,lat,lon)
                # horizontally regrid
                data_grd = data[:,:,::-1,:].interp(lat=model_lat,lon=model_lon)
                # rename dimension name
                data_grd = data_grd.rename({'player':'lev','month':'time'})
                ## vertical regrid to standard levels 
                if data_grd.coords['lev'].max().values > 1300: # pressure in Pa
                    data_grd = data_grd.assign_coords({"lev":data_grd.coords['lev'].values/100})
                data_grd_vrt = data_grd.interp(lev=Dlevs)
                dic_kern[varout] = data_grd_vrt
    
            logger.debug(f'{varout}, {dic_kern[varout].shape},{dic_kern[varout].min().values},{dic_kern[varout].max().values}')
            
    #print('before wv_lw_KernCld',dic_kern['wv_lw_KernCld'][0,:,50,100].values)
    #print('before wv_sw_KernCld',dic_kern['wv_sw_KernCld'][0,:,50,100].values)
    #print('before wv_lw_KernClr',dic_kern['wv_lw_KernClr'][0,:,50,100].values)
    #print('before wv_sw_KernClr',dic_kern['wv_sw_KernClr'][0,:,50,100].values)
    #print('before t_KernCld',dic_kern['t_KernCld'][0,:,50,100].values)
    #print('before t_KernClr',dic_kern['t_KernClr'][0,:,50,100].values)

    # mask out arrays when one has missing value at this level but the other does not 
    dic_kern['wv_lw_KernCld'] = dic_kern['wv_lw_KernCld'].where(dic_kern['wv_sw_KernCld'].notnull())
    dic_kern['wv_lw_KernClr'] = dic_kern['wv_lw_KernClr'].where(dic_kern['wv_sw_KernClr'].notnull())
    dic_kern['t_KernCld'] = dic_kern['t_KernCld'].where(dic_kern['wv_sw_KernCld'].notnull())
    dic_kern['t_KernClr'] = dic_kern['t_KernClr'].where(dic_kern['wv_sw_KernClr'].notnull())

    #print('after wv_lw_KernCld',dic_kern['wv_lw_KernCld'][0,:,50,100].values)
    #print('after wv_sw_KernCld',dic_kern['wv_sw_KernCld'][0,:,50,100].values)
    #print('after wv_lw_KernClr',dic_kern['wv_lw_KernClr'][0,:,50,100].values)
    #print('after wv_sw_KernClr',dic_kern['wv_sw_KernClr'][0,:,50,100].values)
    #print('after t_KernCld',dic_kern['t_KernCld'][0,:,50,100].values)
    #print('after t_KernClr',dic_kern['t_KernClr'][0,:,50,100].values)

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
def Convert_HybridHeight_to_p19(wap,plevel,plev_nm):
    '''
    March 18, 2021: convert hybrid height coordinate to pressure coordinate and finally to standard 19 pressure levels.
    '''
    import metpy
    from metpy.units import units
    import metpy.calc
    import copy

    # get height level
    height = plevel

    # use metpy to convert height to pressure
    height = height * units.meter
    pressure = metpy.calc.height_to_pressure_std(height)
    nan_idx = np.argwhere(np.isnan(pressure))
    # first index with nan
    nan_idx1 = nan_idx[0][0]

    wap_pres = copy.deepcopy(wap)

    wap_pres = wap_pres.rename({plev_nm:'plev'}) # rename from plev to lev
    wap_pres = wap_pres.assign_coords({"plev":pressure})
    wap_pres.coords['plev']['units'] = 'hPa'

    # interp to vertical pressure levels - p19
    wap_interp = wap_pres[:,:nan_idx1,:].interp(plev=Dlevs)

    return wap_interp

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def height_to_pressure_std(height):
    import metpy.constants as mpconst
    from metpy.units import units

    # The following variables are constants for a standard atmosphere
    t0 = units.Quantity(288., 'kelvin')
    p0 = units.Quantity(1013.25, 'hPa')
    gamma = units.Quantity(6.5, 'K/km')

    return p0 * (1 - (gamma / t0) * height) ** (mpconst.g / (mpconst.Rd * gamma))

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def read_model_data(phase,model,filenames,varlist,explist,t1,t2):
    dic_mod = {}
    for svar in varlist:
        for iexp,exp in enumerate(explist):
            # get branch_time_in_parent
            if exp in ['abrupt-4xCO2','abrupt4xCO2']:
                _,branch_time = get_data_fut(model,exp,filenames[exp],svar,t1,t2)

            # read data
            if exp in ['abrupt-4xCO2','abrupt4xCO2']:
                data,_ = get_data_fut(model,exp,filenames[exp],svar,t1,t2)
            elif exp in ['piControl']:
                logger.debug(f'{exp} branch_time={branch_time}')
                data   = get_data_ctl(phase,model,exp,filenames[exp],svar,t1,t2,branch_time_list=branch_time)
            else:
                data   = get_data_ctl_amip(model,exp,filenames[exp],svar,t1,t2)

            #<qinyi 2021-03-20 #------------------
            # it looks like the amip-p4K wap value is too small for GISS-E2-1-G.
            #if svar == 'wap' and model == 'GISS-E2-1-G' and exp == 'amip-p4K':
            if svar == 'wap' and model == 'GISS-E2-1-G' and exp == 'amip':
                logger.debug(f'We are modifying this data for {model} {exp} {svar}')
                data = data * 1000.
            #>qinyi 2021-03-20 #------------------

            # set the same time axis
            if exp in ['piControl','amip']:
                exp_id = 'pi'
                data = data.assign_coords(coords)
            if exp in ['abrupt-4xCO2','abrupt4xCO2','amip-p4K','amip4K']:
                exp_id = 'ab'
                coords = data.coords

            logger.debug(f'filenames={filenames[exp][svar]}')
            logger.debug(f'{svar} {exp_id} {exp} {data.shape}')

            dic_mod[svar+'_'+exp_id] = data

            sub_name = 'Reading '+svar
            print_memory_status(sub_name)

    if len(varlist) != 1: # only read tas 
        ### get SWCRE, LWCRE and netCRE
        dic_mod['SWCRE_pi']  = dic_mod['rsutcs_pi'] - dic_mod['rsut_pi']
        dic_mod['LWCRE_pi']  = dic_mod['rlutcs_pi'] - dic_mod['rlut_pi']
        dic_mod['netCRE_pi'] = dic_mod['SWCRE_pi']  + dic_mod['LWCRE_pi']

        dic_mod['SWCRE_ab']  = dic_mod['rsutcs_ab'] - dic_mod['rsut_ab']
        dic_mod['LWCRE_ab']  = dic_mod['rlutcs_ab'] - dic_mod['rlut_ab']
        dic_mod['netCRE_ab'] = dic_mod['SWCRE_ab']  + dic_mod['LWCRE_ab']

        ## get albedo
        dic_mod['alb_pi']  = np.true_divide(dic_mod['rsus_pi'],dic_mod['rsds_pi']) * 100.
        dic_mod['alb_ab']  = np.true_divide(dic_mod['rsus_ab'],dic_mod['rsds_ab']) * 100.
        # use np.ma.masked_invalid to mask nan although the warning is still there...
        #dic_mod['alb_pi']  = xr.DataArray(np.ma.masked_outside(np.ma.masked_invalid(dic_mod['alb_pi']), 0.0, 100.),coords=dic_mod['rsus_pi'].coords)
        #dic_mod['alb_ab']  = xr.DataArray(np.ma.masked_outside(np.ma.masked_invalid(dic_mod['alb_ab']), 0.0, 100.),coords=dic_mod['rsus_ab'].coords)

        dic_mod['alb_pi']  = xr.DataArray(np.ma.masked_outside(dic_mod['alb_pi'], 0.0, 100.),coords=dic_mod['rsus_pi'].coords)
        dic_mod['alb_ab']  = xr.DataArray(np.ma.masked_outside(dic_mod['alb_ab'], 0.0, 100.),coords=dic_mod['rsus_ab'].coords)

        delterm = ["rsus_pi","rsds_pi",\
                   "rsus_ab","rsds_ab","ps_ab"]
        dic_mod = delete_vars(delterm,dic_mod)

    return dic_mod
 
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def mask_land(lons,lats,data):
    lons_here = np.where(lons>180,lons-360,lons)
    lon_grid,lat_grid = np.meshgrid(lons_here,lats)
    globe_land_mask = globe.is_land(lat_grid,lon_grid)
    if len(data.shape) == len(globe_land_mask.shape):
        globe_land_mask_nd = globe_land_mask
    elif len(data.shape) - len(globe_land_mask.shape) == 1: # have time dimension
        globe_land_mask_nd = np.tile(globe_land_mask,(data.shape[0],1,1))
    elif len(data.shape) - len(globe_land_mask.shape) == 2: # have ctp and tau dimension / time and lev
        globe_land_mask_nd = np.tile(globe_land_mask,(data.shape[0],data.shape[1],1,1))
    elif len(data.shape) - len(globe_land_mask.shape) == 3: # have time, ctp and tau dimension
        globe_land_mask_nd = np.tile(globe_land_mask,(data.shape[0],data.shape[1],data.shape[2],1,1))

    data_ocean = np.ma.masked_where(globe_land_mask_nd==True,data)
    data_land = np.ma.masked_where(globe_land_mask_nd==False,data)

    data_ocean = xr.DataArray(data_ocean, coords=data.coords, dims=data.dims)
    data_land = xr.DataArray(data_land, coords=data.coords, dims=data.dims)

    return data_ocean,data_land

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def get_area_wts(data):
    # Create map of area weights -- from Mark Zelinka
    lats = data.coords['lat'].values
    A,B,C = data.shape
    coslat = np.cos(lats[:]*np.pi/180)
    coslat2 = coslat/np.sum(coslat)
    area_wts = np.ma.array(np.moveaxis(np.tile(coslat2/C,[A,C,1]),1,2)) # summing this over lat and lon = 1
    return area_wts

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
    var1 = (varin[:,:-1,:] + varin[:,1:,:])/2.
    # vertical integral
    outvar = np.ma.sum(var1 * dp4d/100., axis=1)

    return outvar 

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def get_mult(dic_rad, kern1,var1,out1,dic_kern_wap,dic_mod_wap,dp4d):
    '''
    kernel * (xx_ab - xx_pi)
    '''
    for ivar,svar in enumerate(kern1):
        print('svar=',svar,var1[ivar])
        a1 = dic_kern_wap[svar].to_masked_array() # [time,nbins]
        a2 = dic_mod_wap[var1[ivar]].to_masked_array() 
        a1[a1.mask] = 0
        a2[a2.mask] = 0
        outvar = out1[ivar]
        a3 = a1 * a2
        logger.debug(f'kernel={svar}, data={var1[ivar]}, {a1.shape}') 

        if sorting:
            if len(a1.shape) == 2:
                dic_rad[outvar] = xr.DataArray(a3,coords=dic_kern_wap[svar].coords, dims=dic_kern_wap[svar].dims)
            elif len(a1.shape) == 3: #[time,lev,nbins]
                # vertical integral (sum)
                tmp = VertSum(a3, dp4d)
                dic_rad[outvar] = xr.DataArray(tmp,coords=dic_kern_wap[svar][:,0,:].coords, dims=dic_kern_wap[svar][:,0,:].dims)
        else:
            if len(a1.shape) == 3:
                dic_rad[outvar] = xr.DataArray(a3,coords=dic_kern_wap[svar].coords, dims=dic_kern_wap[svar].dims)
            elif len(a1.shape) == 4: #[time,lev,lat,lon]
                # vertical integral (sum)
                tmp = VertSum(a3, dp4d)
                dic_rad[outvar] = xr.DataArray(tmp,coords=dic_kern_wap[svar][:,0,:].coords, dims=dic_kern_wap[svar][:,0,:].dims)

            #tmp = a3
            #dic_rad[outvar] = xr.DataArray(tmp,coords=dic_kern_wap[svar].coords, dims=dic_kern_wap[svar].dims)
 
    return dic_rad

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

    print_memory_status('start get_qs1')

    _,qsl_plus1 = get_qsl(avgta, lev_4d)
    _,qsi_plus1 = get_qsi(avgta, lev_4d)
    print_memory_status('get qsl_plus1 and qsi_plus1')

    # Compute blend of qsi and qsl between -40 and 0
    blend_plus1 = (avgta-233)*qsl_plus1/40 + (273-avgta)*qsi_plus1/40
    print_memory_status('get blend_plus1')

    my_timer = Timer()

    inds = avgta > 233
    inds &= avgta < 273
    qs1 = blend_plus1.where(inds,qsi_plus1)
    print_memory_status('get qs1 A')
    inds2 = avgta >= 273
    qs1 = qs1.where(~inds2,qsl_plus1)
    print_memory_status('get qs1 B')

    #qs1 = np.where((avgta>233) & (avgta<273), blend_plus1, qsi_plus1)#[0]
    #qs1 = np.where(avgta >= 273, qsl_plus1, qs1)#[0]

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
def get_files(datadir,phase,MIP,inst,model,exp,ripf,mon,gr,version,var):
    '''
    Example: 
    datadir = '/p/css03/esgf_publish'
    phase = 'CMIP6'
    MIP = 'CMIP'
    inst = 'NOAA-GFDL'
    model = 'GFDL-CM4'
    exp = 'piControl'
    ripf = 'r1i1p1f1'
    mon = 'Amon'
    var = 'ts'
    gr = 'gr1'
    version = 'v20180701'

    files = get_files(datadir,phase,MIP,inst,model,exp,ripf,mon,gr,version,var)
    f3 = xr.open_mfdataset(files,decode_times=False)
    ts = f3[var]
    '''

    fname = var+'_'+mon+'_'+model+'_'+exp+'_'+ripf+'_'+gr+'_'
    base_url = datadir+'/'+phase+'/'+MIP+'/'+inst+'/'+model+'/'+exp+'/'+ripf+'/'+mon+'/'+var+'/'+gr+'/'+version+'/'+fname+'*'
    logger.debug(f'base_url = {base_url}')
    filenames = glob.glob(base_url)

    return filenames 

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def save_big_dataset(dic_mod,outfile):
    '''
    create a big dataset based on all variables in a dictionary and save to netcdf file.
    '''
    datalist = []
    for svar in dic_mod.keys():
        data = xr.DataArray(dic_mod[svar],name=svar)
        datalist.append(data)

    data_big = xr.merge(datalist,compat='override')

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
def get_tslice(start,nyears):
    endyear = start+nyears-1
    t1 = (start-1)*12
    t2 = endyear*12
    return t1,t2,endyear

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
def sort_var_by_omega(dbin,bins,omega,fdbk):

    ntime = fdbk.shape[0]
    nlev = fdbk.shape[1]

    if len(fdbk.shape) == 4: # [time,lev,lat,lon]
        omega_scaled = np.transpose(np.tile(omega,(nlev,1,1,1)),(1,0,2,3))
        data_bin_avg = np.ma.empty((len(bins),ntime,nlev))
    else:
        omega_scaled = omega
        data_bin_avg = np.ma.empty((len(bins),ntime))

    logger.debug(f'omega_scaled.shape={omega_scaled.shape}, fdbk.shape={fdbk.shape}')

    for ibin,binx in enumerate(bins):
#         print('ibin=',ibin,'binx=',binx)

        if ibin == 0:
            tmp_more = np.ma.masked_where(omega_scaled >= binx+dbin/2.0, fdbk)
        elif ibin == len(bins)-1:
            tmp_more = np.ma.masked_where(omega_scaled < binx-dbin/2.0, fdbk)
        else:
            tmp_less = np.ma.masked_where(omega_scaled < binx-dbin/2.0, fdbk)
            tmp_more = np.ma.masked_where(omega_scaled >= binx+dbin/2.0, tmp_less)

        tmp_more = xr.DataArray(tmp_more, coords=fdbk.coords, dims = fdbk.dims)

        data_bin_avg[ibin,:] = area_averager(tmp_more) #cdutil.averager(tmp_more,axis='xy',weights='weighted')

    data_out = np.moveaxis(data_bin_avg, 0, -1) # move bin to the last dimension

    return data_out

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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
def dt2cal(dt):
    """
    Reference: https://stackoverflow.com/questions/13648774/get-year-month-or-day-from-numpy-datetime64
    Convert array of datetime64 to a calendar array of year, month, day, hour,
    minute, seconds, microsecond with these quantites indexed on the last axis.

    Parameters
    ----------
    dt : datetime64 array (...)
        numpy.ndarray of datetimes of arbitrary shape

    Returns
    -------
    cal : uint32 array (..., 7)
        calendar array with last axis representing year, month, day, hour,
        minute, second, microsecond
    """

    # allocate output 
    out = np.empty(dt.shape + (7,), dtype="u4")
    # decompose calendar floors
    Y, M, D, h, m, s = [dt.astype(f"M8[{x}]") for x in "YMDhms"]
    out[..., 0] = Y + 1970 # Gregorian Year
    out[..., 1] = (M - Y) + 1 # month
    out[..., 2] = (D - M) + 1 # dat
    out[..., 3] = (dt - D).astype("m8[h]") # hour
    out[..., 4] = (dt - h).astype("m8[m]") # minute
    out[..., 5] = (dt - m).astype("m8[s]") # second
    out[..., 6] = (dt - s).astype("m8[us]") # microsecond
    return out

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

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def sort_var_by_regime(omega,EIS,fdbk):

    do_test = True

    fvalue = np.nan

    ntime = fdbk.shape[0]

    lons = fdbk.coords['lon'].data
    lats = fdbk.coords['lat'].data

    # define omega and EIS bins 
    domega  = 10
    #omega_bins = np.arange(-100,110,domega)
    omega_bins = np.arange(-95,95+domega,domega)
    omega_filter = 5
    omega_left = 11 # number of bins smaller than 15 hPa/day (exclude 15)
    omega_right = 9 # number of bins larger than 15 hPa/day (include 15)
    dEIS = 4
    #EIS_bins = np.arange(1,5+dEIS,dEIS)
    EIS_bins = np.arange(5,9+dEIS,dEIS)
    print('omega_bins = ',omega_bins)
    print('EIS_bins = ',EIS_bins)

    # define regime names 
    regimes = [None]*(omega_left+omega_right*len(EIS_bins))
    ii = 0 
    for omega_bin in omega_bins:
        if omega_bin <= omega_filter:
            omega_bin_str = str(abs(omega_bin))+'lhs'
            regimes[ii] = omega_bin_str
            ii += 1
        else:
            omega_bin_str = str(abs(omega_bin))+'rhs'

            for ix,EIS_bin in enumerate(EIS_bins):
                if EIS_bin == EIS_bins[0]:
                    jj = ii
                    ii = ii+1
                else:
                    jj = int(omega_right*ix+ii-1)
                
                EIS_bin_str = str(abs(EIS_bin))+'pos'

                regimes[jj] = omega_bin_str+'_'+EIS_bin_str

    print(regimes)

    regimes_others = ['troplnd','30S60S','30N60N','60S90S','60N90N','glb']
    regimes_all = regimes + regimes_others
    print(regimes_all,len(regimes_all))

    nregime = len(regimes_all)

    # process 4-dimension data 
    if len(fdbk.shape) == 4: # [time,lev,lat,lon]
        nlev = fdbk.shape[1]
        nlat = fdbk.shape[2]
        nlon = fdbk.shape[3]

        omega_scaled = np.transpose(np.tile(omega,(nlev,1,1,1)),(1,0,2,3))
        EIS_scaled = np.transpose(np.tile(EIS,(nlev,1,1,1)),(1,0,2,3))
        omega_scaled = xr.DataArray(omega_scaled, coords=fdbk.coords, dims=fdbk.dims)
        EIS_scaled = xr.DataArray(EIS_scaled, coords=fdbk.coords, dims=fdbk.dims)

        mask_regime_avg = np.ma.empty((nregime,ntime,nlev))
        data_regime_avg = np.ma.empty((nregime,ntime,nlev))
        if do_test:
            data_regime = np.ma.empty((nregime,ntime,nlev,nlat,nlon))
    else:
        nlat = fdbk.shape[1]
        nlon = fdbk.shape[2]

        omega_scaled = omega
        EIS_scaled = EIS

        mask_regime_avg = np.ma.empty((nregime,ntime))
        data_regime_avg = np.ma.empty((nregime,ntime))
        if do_test:
            data_regime = np.ma.empty((nregime,ntime,nlat,nlon))
            
    mask_regime_avg[:] = np.nan
    data_regime_avg[:] = np.nan
    if do_test:
        data_regime[:] = np.nan

    logger.debug(f'omega_scaled.shape={omega_scaled.shape}, fdbk.shape={fdbk.shape}')

    # ================== tropical regimes ======================================
    latS = -30
    latE = 30

    ## mask land and only select tropical data
    omega_scaled_ocn,_ = mask_land(lons,lats,omega_scaled)
    EIS_scaled_ocn,_ = mask_land(lons,lats,EIS_scaled)

    omega_mask = omega_scaled_ocn.where((omega_scaled_ocn['lat']>=latS) & (omega_scaled_ocn['lat']<=latE))
    EIS_mask = EIS_scaled_ocn.where((EIS_scaled_ocn['lat']>=latS) & (EIS_scaled_ocn['lat']<=latE))

    ii = 0
    for omega_bin in omega_bins:
        if omega_bin == omega_bins[0]:
            tmp = xr.where((omega_mask < omega_bin+domega/2) , fdbk,fvalue)
            print('case00', ii, omega_bin,area_averager(tmp)[0])

            tmp_mask = xr.where((omega_mask < omega_bin+domega/2) ,1.0,0.0)

            tmp = xr.DataArray(tmp, coords=fdbk.coords, dims = fdbk.dims)
            if do_test:
                data_regime[ii,:] = tmp
            data_regime_avg[ii,:] = area_averager(tmp)
            mask_regime_avg[ii,:] = area_averager(tmp_mask)

            ii += 1
        elif omega_bin == omega_bins[-1]:
            for ix,EIS_bin in enumerate(EIS_bins):
                if EIS_bin == EIS_bins[0]:
                    tmp = xr.where((omega_mask >= omega_bin-domega/2) & (EIS_mask<EIS_bin+dEIS/2), fdbk,fvalue)
                    print('case10',ii, omega_bin, omega_bin-domega/2, EIS_bin,area_averager(tmp).values[0])

                    tmp_mask = xr.where((omega_mask >= omega_bin-domega/2) & (EIS_mask<EIS_bin+dEIS/2),1,0)

                elif EIS_bin == EIS_bins[-1]:
                    tmp = xr.where((omega_mask >= omega_bin-domega/2) & (EIS_mask>=EIS_bin-dEIS/2), fdbk,fvalue)
                    print('case11',ii, omega_bin, omega_bin-domega/2, EIS_bin,area_averager(tmp).values[0])

                    tmp_mask = xr.where((omega_mask >= omega_bin-domega/2) & (EIS_mask>=EIS_bin-dEIS/2),1,0)

                else:
                    tmp = xr.where((omega_mask >= omega_bin-domega/2) & (EIS_mask>=EIS_bin-dEIS/2) & (EIS_mask<EIS_bin+dEIS/2), fdbk,fvalue)
                    print('case12',ii, omega_bin, omega_bin-domega/2, EIS_bin,area_averager(tmp).values[0])

                    tmp_mask = xr.where((omega_mask >= omega_bin-domega/2) & (EIS_mask>=EIS_bin-dEIS/2) & (EIS_mask<EIS_bin+dEIS/2),1,0)

                tmp = xr.DataArray(tmp, coords=fdbk.coords, dims = fdbk.dims)

                if EIS_bin == EIS_bins[0]:
                    jj = ii 
                    ii += 1
                else:
                    jj = int(omega_right*ix+ii-1)

                if do_test:
                    data_regime[jj,:] = tmp
                data_regime_avg[jj,:] = area_averager(tmp)
                mask_regime_avg[jj,:] = area_averager(tmp_mask)

                print('jj=',jj)

        else:
            if omega_bin <= omega_filter:
                tmp = xr.where((omega_mask >= omega_bin-domega/2) & (omega_mask <omega_bin+domega/2), fdbk,fvalue)
                print('case20', ii, omega_bin, omega_bin-domega/2, omega_bin+domega/2,area_averager(tmp).values[0])
                tmp_mask = xr.where((omega_mask >= omega_bin-domega/2) & (omega_mask <omega_bin+domega/2),1,0)

                tmp = xr.DataArray(tmp, coords=fdbk.coords, dims = fdbk.dims)
                if do_test:
                    data_regime[ii,:] = tmp
                data_regime_avg[ii,:] = area_averager(tmp) 
                mask_regime_avg[ii,:] = area_averager(tmp_mask) 


                ii += 1

            else:
                for ix,EIS_bin in enumerate(EIS_bins):
                    if EIS_bin == EIS_bins[0]:
                        tmp = xr.where((omega_mask >= omega_bin-domega/2) & (omega_mask <omega_bin+domega/2) & (EIS_mask<EIS_bin+dEIS/2) , fdbk,fvalue)

                        print('case30',ii, omega_bin, omega_bin-domega/2, omega_bin+domega/2, EIS_bin,area_averager(tmp).values[0])
                        tmp_mask = xr.where((omega_mask >= omega_bin-domega/2) & (omega_mask <omega_bin+domega/2) & (EIS_mask<EIS_bin+dEIS/2),1,0)


                    elif EIS_bin == EIS_bins[-1]:
                        tmp = xr.where((omega_mask >= omega_bin-domega/2) & (omega_mask <omega_bin+domega/2) & (EIS_mask>=EIS_bin-dEIS/2) , fdbk,fvalue)

                        print('case31',ii, omega_bin, omega_bin-domega/2, omega_bin+domega/2, EIS_bin,area_averager(tmp).values[0])
                        tmp_mask = xr.where((omega_mask >= omega_bin-domega/2) & (omega_mask <omega_bin+domega/2) & (EIS_mask>=EIS_bin-dEIS/2),1,0)

                    else:
                        tmp = xr.where((omega_mask >= omega_bin-domega/2) & (omega_mask <omega_bin+domega/2) & (EIS_mask>=EIS_bin-dEIS/2) & (EIS_mask<EIS_bin+dEIS/2), fdbk,fvalue)

                        print('case32',ii, omega_bin, EIS_bin,area_averager(tmp).values[0])

                        tmp_mask = xr.where((omega_mask >= omega_bin-domega/2) & (omega_mask <omega_bin+domega/2) & (EIS_mask>=EIS_bin-dEIS/2) & (EIS_mask<EIS_bin+dEIS/2),1,0)


                    if EIS_bin == EIS_bins[0]:
                        jj = ii
                        ii += 1
                    else:
                        jj = int(omega_right*ix+ii-1)

                    if do_test:
                        data_regime[jj,:] = tmp
                    data_regime_avg[jj,:] = area_averager(tmp)
                    mask_regime_avg[jj,:] = area_averager(tmp_mask)
                    print('jj=',jj)

    print(data_regime_avg[:,0])
    print(mask_regime_avg[:,0],mask_regime_avg[:-6,0].sum())

    # ================== Other regimes ======================================

    for ireg,regime in enumerate(regimes_others):
        print('ireg=',ireg,'regime=',regime)

        if regime == 'troplnd' :
            _,fdbk_lnd = mask_land(lons,lats,fdbk)
            tmp = fdbk_lnd.where((fdbk['lat']>=latS) & (fdbk['lat']<=latE),fvalue).fillna(fvalue)

        if regime == '60N90N':
            tmp = fdbk.where((fdbk['lat']>60) & (fdbk['lat']<=90),fvalue)

        if regime == '30N60N':
            tmp = fdbk.where((fdbk['lat']>30) & (fdbk['lat']<=60), fvalue)

        if regime == '60S90S':
            tmp = fdbk.where((fdbk['lat']>=-90) & (fdbk['lat']<-60), fvalue)

        if regime == '30S60S':
            tmp = fdbk.where((fdbk['lat']>=-60) & (fdbk['lat']<-30), fvalue)

        if regime == 'glb':
            tmp = fdbk

        tmp_mask = xr.where(np.isnan(tmp),0,1)

        tmp = xr.DataArray(tmp, coords=fdbk.coords, dims = fdbk.dims)

        if do_test:
            data_regime[ireg+len(regimes),:] = tmp

        data_regime_avg[ireg+len(regimes),:] = area_averager(tmp) 
        mask_regime_avg[ireg+len(regimes),:] = area_averager(tmp_mask)


    # test: print regime names and regime averaged values
    print("========== check ==============")
    for ireg,regime in enumerate(regimes_all):
        print(regime,data_regime_avg[ireg,0],mask_regime_avg[ireg,0])
        
    # print the sum of all separated regimes 
    print('Real=',data_regime_avg[-1,0], mask_regime_avg[-1,0])
    print('Check sum of regimes', np.nansum(data_regime_avg[:-1,0]*mask_regime_avg[:-1,0]),np.nansum(mask_regime_avg[:-1,0]))

    # final adjustment 
    data_avg = np.moveaxis(data_regime_avg, 0, -1) # move bin to the last dimension
    mask_avg = np.moveaxis(mask_regime_avg, 0, -1)

    if do_test:
        data = np.moveaxis(data_regime, 0, -1) # move bin to the last dimension

        # redefine coordinates
    
    if len(fdbk.shape) == 4:
        coords = {"time": fdbk.coords['time'].data, "lev": fdbk.coords['lev'].data, "lat": fdbk.coords['lat'].data, "lon": fdbk.coords['lon'].data,"regime":range(len(regimes_all))}
        dims = ['time','lev','lat','lon','regime']
        coords_avg = {"time": fdbk.coords['time'].data, "lev": fdbk.coords['lev'].data,"regime":range(len(regimes_all))}
        dims_avg = ['time','lev','regime']
    else:
        coords = {"time": fdbk.coords['time'].data, "lat": fdbk.coords['lat'].data, "lon": fdbk.coords['lon'].data,"regime":range(len(regimes_all))}
        dims = ['time','lat','lon','regime']
        coords_avg = {"time": fdbk.coords['time'].data,"regime":range(len(regimes_all))}
        dims_avg = ['time','regime']

    if do_test:
        data = xr.DataArray(data, coords=coords,dims=dims)

    data_avg = xr.DataArray(data_avg, coords=coords_avg,dims=dims_avg)
    mask_avg = xr.DataArray(mask_avg, coords=coords_avg,dims=dims_avg)

    print(regimes_all)

    if do_test:
        return data_avg,mask_avg,data
    else:
        return data_avg,mask_avg

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
        kern1 = ['t_KernCld', 't_KernClr', 'Tq_KernCld',  'Tq_KernClr'     ]
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

        # SUNDOWN is not necessary because we don't consider regions outside of tropics. 
        if not sorting:
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
    kernel * (xx_ab - xx_pi)
    '''
    for ivar,svar in enumerate(kern1):
        print('svar=',svar,var1[ivar])
        a1 = dic_kern_wap[svar].to_masked_array() # [time,nbins]
        a2 = dic_mod_wap[var1[ivar]].to_masked_array() 
        a1[a1.mask] = 0
        a2[a2.mask] = 0
        outvar = out1[ivar]
        a3 = a1 * a2
        logger.debug(f'kernel={svar}, data={var1[ivar]}, {a1.shape}') 

        if len(a1.shape) == 3:
            dic_rad[outvar] = xr.DataArray(a3,coords=dic_kern_wap[svar].coords, dims=dic_kern_wap[svar].dims)
        elif len(a1.shape) == 4: #[time,lev,lat,lon]
            # vertical integral (sum)
            tmp = VertSum(a3, dp4d)
            dic_rad[outvar] = xr.DataArray(tmp,coords=dic_kern_wap[svar][:,0,:].coords, dims=dic_kern_wap[svar][:,0,:].dims)

        #tmp = a3
        #dic_rad[outvar] = xr.DataArray(tmp,coords=dic_kern_wap[svar].coords, dims=dic_kern_wap[svar].dims)
 
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
def sort_var_wapEIS(dic_mod,dic_rad,case_stamp,outdir):
    '''
    sort kernel data and variables 
    output: model and kernel data in defined regimes: [time,(lev),nregimes]
    '''
    # read wap 
    ctl_wap = dic_mod['OMEGA_pi'].sel(lev=500,method='nearest')
    fut_wap = dic_mod['OMEGA_ab'].sel(lev=500,method='nearest')


    # calculate EIS 
    ctl_EIS,_ = calc_EIS(dic_mod['ts_pi'], dic_mod['psl_pi'], dic_mod['ta_pi'].sel(lev=700,method='nearest'), dic_mod['Z3_pi'].sel(lev=700,method='nearest'))
    fut_EIS,_ = calc_EIS(dic_mod['ts_ab'], dic_mod['psl_ab'], dic_mod['ta_ab'].sel(lev=700,method='nearest'), dic_mod['Z3_ab'].sel(lev=700,method='nearest'))
    print('ctl_EIS.shape=',ctl_EIS.shape,ctl_EIS.min().values,ctl_EIS.max().values)

    # Make sure wap is in hPa/day
    ctl_wap = 36*24*ctl_wap # Pa/s -> hPa/day
    fut_wap = 36*24*fut_wap 

    # get lats and lons
    lats = ctl_wap.coords['lat']
    lons = ctl_wap.coords['lon']

    # sort model variables into regimes 
    dic_rad_wap = {}
    for svar in dic_rad.keys():
    #for svar in ['TS_pi','TS_ab']:

        logger.debug(f'=======>>>>>>> svar ={svar} is doing sort based on regime definition')
        data_glb = dic_rad[svar]

        if '_pi' in svar:
            omega = ctl_wap
            EIS = ctl_EIS
        else:
            omega = fut_wap
            EIS = fut_EIS

        if svar in ['T_mask_pi','T_mask_ab','WV_lw_mask_pi','WV_lw_mask_ab','WV_sw_mask_pi','WV_sw_mask_ab','LW_adj_pi','LW_adj_ab','SW_adj_pi','SW_adj_ab','LWCRE_pi','LWCRE_ab','SWCRE_pi','SWCRE_ab','netCRE_pi','netCRE_ab','LWCRE_adj_pi','LWCRE_adj_ab','SWCRE_adj_pi','SWCRE_adj_ab','netCRE_adj_pi','netCRE_adj_ab']:
            DATA_avg,mask_avg,DATA = sort_var_by_regime(omega,EIS,data_glb)
            dich = {}
            dich[svar] = DATA
            dich[svar+'_avg'] = DATA_avg
            dich['mask_avg'] = mask_avg
            print('DATA.shape=',DATA.shape, 'DATA_avg.shape=',DATA_avg.shape)
            save_big_dataset(dich,outdir+"/middata/Regime_Sorted_"+svar+"_"+case_stamp+".nc")
        else:
            DATA_avg,mask_avg,DATA = sort_var_by_regime(omega,EIS,data_glb) 

        dic_rad_wap[svar] = (DATA_avg,mask_avg)

    return dic_rad_wap

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

            #dic_invar[svar+'_ano'] = dic_invar[svar+'_ab'] - dic_invar[svar+'_pi']
 
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

if __name__ == "__main__":

    RadKernel_dir = '/qfs/people/qiny108/diag_feedback_E3SM/Huang_kernel_data/'
    direc_data = '/compyfs/qiny108/diag_feedback_E3SM_postdata/'
    case_stamp = 'v1'
    yearS = 2
    yearE = 6
    fname1 = CL.get_lutable(case_stamp,'amip')
    fname2 = CL.get_lutable(case_stamp,'amip4K')
    outdir = './'
    figdir = './'
    exp1 = 'FC5'
    exp2 = 'FC5_4K'
    
    RadKernel_regime(RadKernel_dir,direc_data,case_stamp,yearS,yearE,fname1,fname2,outdir,figdir,exp1,exp2)

