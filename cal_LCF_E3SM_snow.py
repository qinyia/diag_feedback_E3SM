
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
import cases_lookup as CL
 
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
    nbin = int((maxbin - minbin)/dbin)

    outfile = outdir+'LCF_binned_by_temperature_'+case_stamp+'_'+str(nbin)+'bins-OnlyOcean-CTLandP4K.csv'

    if os.path.isfile(outfile):
        print('====> cal_LCF', case_stamp, 'output is ready. Please continue. ')

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
    
        # ======================== read data =========================================
        Vars = ['CLDLIQ','CLDICE','SNOWQM','T']
        data_CTL = {}
        for svar in Vars:
            print('Reading amip '+svar+' ....')
            f1=cdms.open(direc_data1+svar+'_'+exp1+'_'+yrS_4d+monS_2d+'-'+yrE_4d+monE_2d+'.nc','r')
            clw = f1(svar)
            clw1 = clw.regrid(grid,regridTool='esmf',regridMethod='linear')
            print('clw1.shape = ',clw1.shape, genutil.minmax(clw1))
            f1.close()

            data_CTL[svar] = clw 

        data_P4K = {}
        for svar in Vars:
            print('Reading amip4K '+svar+' ....')
            f1=cdms.open(direc_data2+svar+'_'+exp2+'_'+yrS_4d+monS_2d+'-'+yrE_4d+monE_2d+'.nc','r')
            clw = f1(svar)
            clw1 = clw.regrid(grid,regridTool='esmf',regridMethod='linear')
            print('clw1.shape = ',clw1.shape, genutil.minmax(clw1))
            f1.close()

            data_P4K[svar] = clw

        print(data_CTL.keys())
        print(data_P4K.keys())

        # calculate Liquid Condensate Fraction (LCF) at each vertical level, each month, each model
        limiter = 1e-7
        clt_CTL = data_CTL['CLDLIQ'] + data_CTL['CLDICE'] + data_CTL['SNOWQM']
        clt_CTL = MV.masked_less(clt_CTL,limiter)
        LCF_CTL = data_CTL['CLDLIQ'] /clt_CTL
    
        clt_P4K = data_P4K['CLDLIQ'] + data_P4K['CLDICE'] + data_P4K['SNOWQM']
        clt_P4K = MV.masked_less(clt_P4K,limiter)
        LCF_P4K = data_P4K['CLDLIQ'] /clt_P4K

        # ===============================================
        # compositing procedure
        # ===============================================
    
        print('Start sorting LCF by temperature ...')
        binned_LCF_CTL = sort_var_by_temperature(dbin,bins,data_CTL['T'],LCF_CTL)
        binned_LCF_P4K = sort_var_by_temperature(dbin,bins,data_P4K['T'],LCF_P4K)

        print('binned_LCF_CTL.shape = ',binned_LCF_CTL.shape)
        print('binned_LCF_P4K.shape = ',binned_LCF_P4K.shape)
   
        # save data_new into panda csv file for each model
    
        df_each = pd.DataFrame({'temp':bins,'bin_CTL':binned_LCF_CTL, 'bin_P4K':binned_LCF_P4K})
        df_each.to_csv(outfile)

        print(df_each)

#####################################################################
### LCF profile over 30S-80S
#####################################################################
def plot_LCF(datadir,cases_here):
    print('====> plot_LCF starts ..........')

    nbin = 20
    # E3SM
    data_all = np.zeros((nbin,len(cases_here)))

    df_all = pd.DataFrame()
    for icase,case in enumerate(cases_here):

        infile = outdir+'LCF_binned_by_temperature_'+case_stamp+'_'+str(nbin)+'bins-OnlyOcean-CTLandP4K.csv'

        df1 = pd.read_csv(infile,index_col=0)
        df1.index = df1.loc[:,'temp']

        df_all[case+'_CTL'] = df1.loc[:,'bin_CTL']
        df_all[case+'_P4K'] = df1.loc[:,'bin_P4K']

    print(df_all)

    #----------------------------------------------------------
    # start plotting ...
    #----------------------------------------------------------
    fig = plt.figure(figsize=(9,6))
    ax = fig.add_subplot(1,1,1)

    for col in df_all.columns:
        print('col=',col)
        ax.plot(df_all.index, df_all.loc[:,col], label=col)

    ax.set_xlabel('Temperature [K]')
    ax.set_ylabel('Liquid Condensate Fraction') 

    ax.axvline(x=273.15,ls='--',c='grey',lw=1.0)
    ax.axhline(y=0.50,ls='--',c='grey',lw=1.0)

    ax.grid(ls=':',alpha=0.7)

    ax.legend()

    fig.savefig(figdir+'LCF-vs-temperature_'+cases_here[-1]+'.png',bbox_inches='tight',dpi=300)
    plt.close(fig)
           
    print('------------------------------------------------')
    print('plot_LCF is done!')
    print('------------------------------------------------')


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if __name__ == "__main__":

    direc_data = '/compyfs/qiny108/diag_feedback_E3SM_postdata/'

    case_stamps = [\
    'v2.bk.trig_ull'
    ]

    for case_stamp in case_stamps:

        fname1 = CL.get_lutable(case_stamp,'amip')
        fname2 = CL.get_lutable(case_stamp,'amip4K')

        outdir = './'
        figdir = './'

        exp1 = 'FC5'
        exp2 = 'FC5_4K'

        yearS = 2
        yearE = 6

        cal_LCF(direc_data,case_stamp,yearS,yearE,fname1,fname2,outdir,figdir)

    plot_LCF(outdir, case_stamps)

