
#IMPORT STUFF:
#=====================
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import os
import pandas as pd
from global_land_mask import globe 
import sys
sys.path.append('../')
import cases_lookup as CL
 
###########################################################################
# HELPFUL FUNCTIONS FOLLOW
###########################################################################

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


########## MAIN SUBROUTINE STARTS HERE ....
def sort_var_by_temperature(dbin,bins,ref_var,bin_var):
    # ==========================================

    ntime = ref_var.shape[0]
    nlevs = ref_var.shape[1]
    nlats = ref_var.shape[2]
    nlons = ref_var.shape[3]

    #data_bin = np.ma.empty((len(bins),ntime,nlevs,nlats,nlons))
    #data_bin_avg = np.ma.empty((len(bins),ntime))
    data_bin_avg = np.ma.empty((len(bins)))

    for ibin,binx in enumerate(bins):
        if ibin == 0:
            tmp_more = bin_var.where(ref_var < binx+dbin/2.0)
        elif ibin == len(bins)-1:
            tmp_more = bin_var.where(ref_var >= binx-dbin/2.0)
        else:
            tmp_less = bin_var.where(ref_var >= binx-dbin/2.0)
            tmp_more = tmp_less.where(ref_var < binx+dbin/2.0)
            print('binx=',binx, [tmp_less.min().values,tmp_less.max().values], [tmp_more.min().values,tmp_more.max().values])

        #data_bin[ibin,:,:,:,:] = tmp_more

        # time and level average
        tmp_more_avg = tmp_more.mean(axis=0).mean(axis=0)
        print('tmp_more_avg.shape=',tmp_more_avg.shape)

        # mask land regions
        lons = tmp_more_avg.lon.data
        lats = tmp_more_avg.lat.data
        lons_here = np.where(lons>180,lons-360,lons)
        lon_grid,lat_grid = np.meshgrid(lons_here,lats)
        globe_land_mask = globe.is_land(lat_grid,lon_grid)
     
        tmp_more_avg_mask = tmp_more_avg.where(globe_land_mask==False)

        data_bin_avg[ibin] = area_averager(tmp_more_avg_mask.sel(lat=slice(-80,-30)))

    #return data_bin,data_bin_avg
    return data_bin_avg


def cal_LCF(direc_data,case_stamp,yearS,yearE,fname1,fname2,outdir,figdir):

    # ============= sorting setting ==================
    dbin = 4 # set the temperature bin size is 4 K
    minbin = 210
    maxbin = 290
    bins = range(minbin, maxbin, dbin)
    nbin = int((maxbin - minbin)/dbin)

    if os.path.isfile(outdir+'LCF_binned_by_temperature_'+case_stamp+'_'+str(nbin)+'bins-OnlyOcean.csv'):
        print('cal_LCF', case_stamp, 'output is ready. Please continue. ')

    else:

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
    
        lats_spc = np.arange(-90,92.5,2.5)
        lons_spc = np.arange(1.25,360,2.5)
    
        # read CLDLIQ
        print('Reading CLDLIQ....')
        f1=xr.open_dataset(direc_data1+'CLDLIQ_'+exp1+'_'+yrS_4d+monS_2d+'-'+yrE_4d+monE_2d+'.nc')
        clw = f1['CLDLIQ']
        clw1 = clw.interp(lat=lats_spc,lon=lons_spc)
        print('clw1.shape = ',clw1.shape)
        f1.close()
        
        # read CLDICE
        print('Reading CLDICE....')
        f1=xr.open_dataset(direc_data1+'CLDICE_'+exp1+'_'+yrS_4d+monS_2d+'-'+yrE_4d+monE_2d+'.nc')
        cli = f1['CLDICE']
        cli1 = cli.interp(lat=lats_spc,lon=lons_spc)
        print('cli1.shape = ',cli1.shape)
        f1.close()
    
        # read ta
        print('Reading ta....')
        #<qinyi 2021-08-12 #------------------
        if not os.path.isfile(direc_data1+'ta_'+exp1+'_'+yrS_4d+monS_2d+'-'+yrE_4d+monE_2d+'.nc'):
            svar_in = 'T'
        else:
            svar_in = 'ta'
        #>qinyi 2021-08-12 #------------------

        f1=xr.open_dataset(direc_data1+svar_in+'_'+exp1+'_'+yrS_4d+monS_2d+'-'+yrE_4d+monE_2d+'.nc')
        ta=f1[svar_in]
        ta1 = ta.interp(lat=lats_spc,lon=lons_spc)
        print('ta1.shape = ',ta1.shape)
    
        # calculate Liquid Condensate Fraction (LCF) at each vertical level, each month, each model
        limiter = 1e-7
        clt1 = clw1 + cli1
        clt1 = clt1.where(clt1>=limiter)
        LCF1 = clw1 /clt1
    
        del clw1, cli1, clt1 
    
        # ===============================================
        # compositing procedure
        # ===============================================
   
        print('Start sorting LCF by temperature ...')
        binned_LCF = sort_var_by_temperature(dbin,bins,ta1,LCF1)
        print('binned_LCF.shape = ',binned_LCF.shape)
    
        # save data_new into panda csv file for each model
        df_each = pd.DataFrame({'temp':bins,'bin':binned_LCF})
        df_each.to_csv(outdir+'LCF_binned_by_temperature_'+case_stamp+'_'+str(nbin)+'bins-OnlyOcean.csv')
    
        print(df_each)


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if __name__ == "__main__":

    direc_data = '/compyfs/qiny108/diag_feedback_E3SM_postdata/'

    case_stamps = [\
    'v2test'
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

        cal_LCF(direc_data,case_stamp,yearS,yearE,fname1,fname2,outdir,figdir)

