
#IMPORT STUFF:
#=====================
import numpy as np
import xarray as xr
import os
import pandas as pd
from global_land_mask import globe 
import sys
sys.path.append('../')
import cases_lookup as CL
from PlotDefinedFunction import area_averager 
import xrw
from loguru import logger
from get_mip_data import read_mip_data,read_amip_data,read_pickle,write_pickle,read_e3sm_data,read_e3smdiag_data

# ============= sorting setting ==================
dbin = 4 # set the temperature bin size is 4 K
minbin = 210
maxbin = 290
bins = range(minbin, maxbin, dbin)
nbin = int((maxbin - minbin)/dbin)

### Define standard pressure levels  [Pa]
Dlevs = np.array([100000, 92500, 85000, 70000, 60000, 50000, 40000, 30000, 25000, 20000, 15000, 10000, 7000, 5000, 3000, 2000, 1000, 500, 100])

# Main 
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

        data_bin_avg[ibin] = area_averager(tmp_more_avg_mask.sel(lat=slice(-80,-30))).values

    #return data_bin,data_bin_avg
    return data_bin_avg

def cal_LCF_MIP(case_stamp,outdir,figdir,filenames,tslice):

    Vars = ['clw','cli','ta']

    outfile = outdir+'LCF_binned_by_temperature_'+case_stamp+'_'+str(nbin)+'bins-OnlyOcean.csv'

    if os.path.isfile(outfile):
        print('cal_LCF', case_stamp, 'output is ready. Please continue. ')
        return 

    # read model's data
    dic_mod = read_mip_data(Vars,filenames,tslice)

    calculation(dic_mod, outfile)

    return 


def cal_LCF_E3SM(direc_data,case_stamp,yearS,yearE,fname1,fname2,outdir,figdir,UseE3smDiagOutput=False,grd_info=None,num_years_per_file=None):

    Vars = ['clw','cli','ta']

    outfile = outdir+'LCF_binned_by_temperature_'+case_stamp+'_'+str(nbin)+'bins-OnlyOcean.csv'

    if os.path.isfile(outfile):
        print('cal_LCF', case_stamp, 'output is ready. Please continue. ')
        return 

    nyears = yearE - yearS + 1 

    # read model's data
    if UseE3smDiagOutput:
        dic_mod = read_e3smdiag_data(Vars,direc_data,case_stamp,yearS,yearE,fname1,fname2,grd_info,num_years_per_file)
    else:
        dic_mod = read_e3sm_data(Vars,direc_data,case_stamp,yearS,yearE,fname1,fname2)

    calculation(dic_mod, outfile)

def calculation(dic_mod, outfile):

    # Regrid data
    dic_in = {}
    lats_spc = np.arange(-90,92.5,2.5)
    lons_spc = np.arange(1.25,360,2.5)
    for svar in dic_mod.keys():
        dic_in[svar] = dic_mod[svar].interp(lat=lats_spc,lon=lons_spc)
    
    ta1 = dic_in['ta_pi']
    if 'plev' in list(ta1.coords):
        ta1 = ta1.rename({'plev':'lev'})

    # calculate Liquid Condensate Fraction (LCF) at each vertical level, each month, each model
    limiter = 1e-7
    clt1 = dic_in['clw_pi'] + dic_in['cli_pi']
    clt1 = clt1.where(clt1>=limiter)
    LCF1 = dic_in['clw_pi'] /clt1

    del dic_in, clt1

    # ===============================================
    # compositing procedure
    # ===============================================
   
    print('Start sorting LCF by temperature ...')
    binned_LCF = sort_var_by_temperature(dbin,bins,ta1,LCF1)
    print('binned_LCF.shape = ',binned_LCF.shape)
    
    # save data_new into panda csv file for each model
    df_each = pd.DataFrame({'temp':bins,'bin':binned_LCF})
    df_each.to_csv(outfile)
    
    print(df_each)


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if __name__ == "__main__":

    if False: 
        direc_data = '/p/user_pub/climate_work/qin4/From_Compy/compyfs_dir/colla/diag_feedback_E3SM_postdata/'
    
        case_stamp = 'v2test'
    
        fname1,_,_ = CL.get_lutable(case_stamp,'amip')
        fname2,_,_ = CL.get_lutable(case_stamp,'amip4K')
    
        outdir = './'
        figdir = './'
    
        yearS = 2
        yearE = 3
    
        cal_LCF_E3SM(direc_data,case_stamp,yearS,yearE,fname1,fname2,outdir,figdir)


    # ------------------------------
    # Input requirements for MIP data
    # ------------------------------
    #model = 'GFDL-CM4'  
    #institution = 'NOAA-GFDL'

    model = 'CESM2'
    institution = 'NCAR'
    variant = 'r1i1p1f1'

    # ------------ get filenames ----------------------------------
    ff = 'filenames_'+model+'_'+variant+'.pickle'

    filenames = read_pickle(ff)
    #print('filenames=',filenames)

    case_stamp = model+'_'+variant
    outdir = './'
    figdir = './'

    tslice = slice("1980-01-01","1981-12-31")

    cal_LCF_MIP(case_stamp,outdir,figdir,filenames,tslice)



