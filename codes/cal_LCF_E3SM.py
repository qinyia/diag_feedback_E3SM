
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
import matplotlib.pyplot as plt
import os
import pandas as pd
import cdtime
import time
import ReadData as RD
import genutil
import numpy.ma as ma
from global_land_mask import globe 
from joblib import Parallel,delayed
import multiprocessing
 
###########################################################################
# HELPFUL FUNCTIONS FOLLOW
###########################################################################

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
            tmp_more = MV.masked_where(ref_var >= binx+dbin/2.0, bin_var)
        elif ibin == len(bins)-1:
            tmp_more = MV.masked_where(ref_var < binx-dbin/2.0, bin_var)
        else:
            tmp_less = MV.masked_where(ref_var < binx-dbin/2.0, bin_var)
            tmp_more = MV.masked_where(ref_var >= binx+dbin/2.0, tmp_less)
            print('binx=',binx, genutil.minmax(tmp_less), genutil.minmax(tmp_more))

        #data_bin[ibin,:,:,:,:] = tmp_more

        # time and level average
        tmp_more_avg = MV.average(MV.average(tmp_more,axis=0),axis=0)
        print('tmp_more_avg.shape=',tmp_more_avg.shape)

        # mask land regions
        lons = tmp_more_avg.getLongitude()[:]
        lats = tmp_more_avg.getLatitude()[:]
        lons_here = np.where(lons>180,lons-360,lons)
        lon_grid,lat_grid = np.meshgrid(lons_here,lats)
        globe_land_mask = globe.is_land(lat_grid,lon_grid)
     
        tmp_more_avg_mask = MV.masked_where(globe_land_mask==True,tmp_more_avg)

        data_bin_avg[ibin] = cdutil.averager(tmp_more_avg_mask.subRegion(latitude=(-80,-30)), axis='xy',weights='weighted')

    #return data_bin,data_bin_avg
    return data_bin_avg


def cal_LCF(direc_data,case_stamp,yearS,yearE,fname1,fname2,outdir,figdir):

    # ============= sorting setting ==================
    dbin = 4 # set the temperature bin size is 4 K
    minbin = 210
    maxbin = 290
    bins = range(minbin, maxbin, dbin)
    nbin = np.int((maxbin - minbin)/dbin)

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
    
        lats = np.arange(-90,92.5,2.5)
        lons = np.arange(1.25,360,2.5)
        grid = cdms.createGenericGrid(lats,lons)
    
    
        # read CLDLIQ
        print('Reading CLDLIQ....')
        f1=cdms.open(direc_data1+'CLDLIQ_'+exp1+'_'+yrS_4d+monS_2d+'-'+yrE_4d+monE_2d+'.nc','r')
        clw = f1('CLDLIQ') 
        clw1 = clw.regrid(grid,regridTool='esmf',regridMethod='linear')
        print('clw1.shape = ',clw1.shape, genutil.minmax(clw1))
        f1.close()
        
        # read CLDICE
        print('Reading CLDICE....')
        f1=cdms.open(direc_data1+'CLDICE_'+exp1+'_'+yrS_4d+monS_2d+'-'+yrE_4d+monE_2d+'.nc','r')
        cli = f1('CLDICE') 
        cli1 = cli.regrid(grid,regridTool='esmf',regridMethod='linear')
        print('cli1.shape = ',cli1.shape, genutil.minmax(cli1))
        f1.close()
    
        # read ta
        print('Reading ta....')
        #<qinyi 2021-08-12 #------------------
        if not os.path.isfile(direc_data1+'ta_'+exp1+'_'+yrS_4d+monS_2d+'-'+yrE_4d+monE_2d+'.nc'):
            svar_in = 'T'
        else:
            svar_in = 'ta'
        #>qinyi 2021-08-12 #------------------

        f1=cdms.open(direc_data1+svar_in+'_'+exp1+'_'+yrS_4d+monS_2d+'-'+yrE_4d+monE_2d+'.nc','r')
        ta=f1(svar_in) 
        ta1 = ta.regrid(grid,regridTool='esmf',regridMethod='linear')
        print('ta1.shape = ',ta1.shape, genutil.minmax(ta1))
    
        # calculate Liquid Condensate Fraction (LCF) at each vertical level, each month, each model
        limiter = 1e-7
        clt1 = clw1 + cli1
        clt1 = MV.masked_less(clt1,limiter)
        LCF1 = clw1 /clt1
    
        del clw1, cli1, clt1 
    
        # ===============================================
        # compositing procedure
        # ===============================================
    
        #binned_LCF_4d,binned_LCF = sort_var_by_temperature(dbin,bins,ta1,LCF1)
        #print('binned_LCF_4d.shape = ',binned_LCF_4d.shape)
    
        print('Start sorting LCF by temperature ...')
        binned_LCF = sort_var_by_temperature(dbin,bins,ta1,LCF1)
        print('binned_LCF.shape = ',binned_LCF.shape)
    
    
        ## output binned_LCF_4d
        #out = cdms.open(figdir+'LCF-binned-by-temperature-with-lat-amip-'+case_stamp+'-'+str(nbin)+'bins-OnlyOcean.nc','w')
        #out.write(binned_LCF_4d)
        #out.close()
    
        # save data_new into panda csv file for each model
    
        df_each = pd.DataFrame({'temp':bins,'bin':binned_LCF})
        df_each.to_csv(outdir+'LCF_binned_by_temperature_'+case_stamp+'_'+str(nbin)+'bins-OnlyOcean.csv')
    
        print(df_each)

