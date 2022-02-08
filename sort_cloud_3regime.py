#****************************************************************
#
#    Filename: sort_cloud_3regime.py
#
#    Author: Yi Qin - qin4@llnl.gov
#    Description: Use monthly EIS and omega500 to sort cloud regimes: Sc, Cu and Ascent
#    Input: 
#    Output: 
#    Create: 2022-02-06 15:39:52
#    Last Modified: 2022-02-07 15:39:52
#****************************************************************

#=====================
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
from genutil import statistics
import numpy.ma as ma
from global_land_mask import globe 
import cal_EIS as calEIS
import PlotDefinedFunction as PDF
#################################################################################
def sort_cloud_3regime(direc_data,case_stamp,yearS,yearE,fname1,fname2,outdir,figdir,exp1,exp2):

    outfile = outdir+'global_cloud_3regime_'+case_stamp+'.nc'

    if os.path.isfile(outfile):
        print('sort_cloud_3regime', case_stamp, 'output is ready. Please continue. ')
    else:
        yearS_4d = "{:04d}".format(yearS)
        yearE_4d = "{:04d}".format(yearE)
        nyears = yearE - yearS + 1
    
        direc_data1 = direc_data+'/'+fname1+'/'
        direc_data2 = direc_data+'/'+fname2+'/'
    
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

        var3d = ['CLOUD','CLDLIQ','CLDICE','T','Q','OMEGA']
        var2d = ['psl','ts','tas','CLDLOW']
        vara = ['T700','OMEGA500','Z700','OMEGA700']

        var = var2d + var3d + vara

        # =============================================
        # read 2D variables
        # =============================================
        dic_all = {}

        for svar in var:
            print(svar)

            if svar in vara:
                if svar == 'T700':
                    svarin = 'T'
                elif svar == 'OMEGA500':
                    svarin = 'OMEGA'
                elif svar == 'OMEGA700':
                    svarin = 'OMEGA'
                elif svar == 'Z700':
                    svarin = 'Z3'
            else:
                svarin = svar 

            f1 = cdms.open(direc_data+fname1+'/'+svarin+'_'+exp1+'_'+yearS_4d+'01-'+yearE_4d+'12.nc')
            pi_raw = f1(svarin)
            f1.close()

            f2 = cdms.open(direc_data+fname2+'/'+svarin+'_'+exp2+'_'+yearS_4d+'01-'+yearE_4d+'12.nc')
            ab_raw = f2(svarin)
            f2.close()

            if svar in vara:
                # by default, the level for raw model ouput is in unit of hPa. Be cautious.
                if svar in ['Z700','T700']:
                    spec_lev = 700
                elif svar in ['OMEGA500']:
                    spec_lev = 500
                elif svar in ['OMEGA700']:
                    spec_lev = 700

                pi_raw = pi_raw.pressureRegrid(cdms.createAxis([spec_lev]))[:,0,:,:]
                ab_raw = ab_raw.pressureRegrid(cdms.createAxis([spec_lev]))[:,0,:,:]
            elif svar in var3d:
                # reverse lev direction
                pi_raw = pi_raw[:,::-1,:,:]
                ab_raw = ab_raw[:,::-1,:,:]

            #----------------------------------------------------------
            # regrid data                 
            #----------------------------------------------------------
            pi_raw_grd = pi_raw.regrid(grid,regridTool='esmf',regridMethod='linear')
            ab_raw_grd = ab_raw.regrid(grid,regridTool='esmf',regridMethod='linear')

            print('pi_raw_grd.shape = ',pi_raw_grd.shape,genutil.minmax(pi_raw_grd))
            print('ab_raw_grd.shape = ',ab_raw_grd.shape,genutil.minmax(ab_raw_grd))

            dic_all[svar+'_ano'] = ab_raw_grd - pi_raw_grd
            dic_all[svar+'_pi'] = pi_raw_grd
            dic_all[svar+'_ab'] = ab_raw_grd

            dic_all[svar+'_ano'].setAxisList(pi_raw_grd.getAxisList())
            dic_all[svar+'_pi'].setAxisList(pi_raw_grd.getAxisList())
            dic_all[svar+'_ab'].setAxisList(pi_raw_grd.getAxisList())

            del(pi_raw, ab_raw, pi_raw_grd, ab_raw_grd)

        print(dic_all.keys())

        #----------------------------------------------------------
        # calculate EIS
        #----------------------------------------------------------
        EIS_pi,LTS_pi = calEIS.calc_EIS(dic_all['ts_pi'], dic_all['psl_pi'], dic_all['T700_pi'], dic_all['Z700_pi'])
        EIS_ab,LTS_ab = calEIS.calc_EIS(dic_all['ts_ab'], dic_all['psl_ab'], dic_all['T700_ab'], dic_all['Z700_ab'])

        dic_all['EIS_pi'] = EIS_pi
        dic_all['EIS_ab'] = EIS_ab
        dic_all['EIS_ano'] = EIS_ab - EIS_pi

        AXL = dic_all['ts_pi'].getAxisList()
        dic_all['EIS_ano'].setAxisList(AXL)
        dic_all['EIS_pi'].setAxisList(AXL)
        dic_all['EIS_ab'].setAxisList(AXL)

        #----------------------------------------------------------
        # sort variables in three cloud regimes
        #----------------------------------------------------------
        value = 0
        cdms.setNetcdfShuffleFlag(value) ## where value is either 0 or 1
        cdms.setNetcdfDeflateFlag(value) ## where value is either 0 or 1
        cdms.setNetcdfDeflateLevelFlag(value) ## where value is a integer between 0 and 9 included

        fout = cdms.open(outdir+'global_cloud_3regime_'+case_stamp+'.nc','w')

        var_sort = ['CLOUD','CLDLIQ','CLDICE','T','Q','OMEGA']

        for svar in var_sort:
            print('var_sort=',svar)
            for case in ['pi','ab']:

                latS = -30
                latE = 30 
                data = dic_all[svar+'_'+case].subRegion(lat=(latS,latE))
                EIS = dic_all['EIS_'+case].subRegion(lat=(latS,latE))
                omega500 = dic_all['OMEGA500_'+case].subRegion(lat=(latS,latE))
                omega700 = dic_all['OMEGA700_'+case].subRegion(lat=(latS,latE))

                # change omega from Pa/s to hPa/day
                omega500 = omega500 * 864.
                omega700 = omega700 * 864.

                lons = data.getLongitude()[:]
                lats = data.getLatitude()[:]

                #----------------------------------------------------------
                # first, mask land 
                #----------------------------------------------------------
                data_ocn = PDF.mask_land(lons,lats,data,land=True)
                data_ocn.setAxisList(data.getAxisList())

                #----------------------------------------------------------
                # sort variables in each cloud regime 
                #----------------------------------------------------------
                nbins = 3
                if len(data.shape) == 3: #(time,lat,lon)
                    ntime,nlat,nlon = data.shape[0],data.shape[1],data.shape[2]
                    freq_bins = MV.zeros((nbins,ntime,nlat,nlon)) 
                    freq_bins_avg = MV.zeros((nbins,ntime)) 

                    data_bins = MV.zeros((nbins,ntime,nlat,nlon)) 
                    data_bins_avg = MV.zeros((nbins,ntime)) 

                elif len(data.shape) == 4: #(time,lev,lat,lon)
                    ntime,nlev,nlat,nlon = data.shape[0],data.shape[1],data.shape[2],data.shape[3]
                    freq_bins = MV.zeros((nbins,ntime,nlev,nlat,nlon)) 
                    freq_bins_avg = MV.zeros((nbins,ntime,nlev)) 

                    data_bins = MV.zeros((nbins,ntime,nlev,nlat,nlon)) 
                    data_bins_avg = MV.zeros((nbins,ntime,nlev)) 

                # Sc 
                MASK_Sc = (omega500 >= 10 ) & (omega700 >= 10) & (EIS > 3)
#                MASK_Sc = (omega500 >= 10 ) & (EIS > 3)
                # Cu
                MASK_Cu = (omega500 >= 10 ) & (omega700 >= 10) & (EIS <= 3) 
#                MASK_Cu = (omega500 >= 10 ) & (EIS <= 3) 
                # Ascent
                MASK_As = (omega500 < 0 )

                for ibin,MASK in zip(range(nbins),[MASK_Sc,MASK_Cu,MASK_As]):

                    if len(data.shape) == 4:
                        MASK_4d = np.transpose(np.tile(MASK,(nlev,1,1,1)),(1,0,2,3))
                        print('MASK_4d.shape=',MASK_4d.shape)
    
                        freq_bins[ibin,:,:,:,:] = MASK_4d
                        freq_bins_avg[ibin,:,:] = MV.sum(MASK_4d)
    
                        data_bins[ibin,:,:,:,:] = MV.masked_where(MASK_4d==False,data_ocn)
                        data_bins_avg[ibin,:,:] = cdutil.averager(MV.masked_where(MASK_4d==False,data_ocn),axis='xy',weights='weighted')
    
                    else:
    
                        MASK_4d = MASK 
                        freq_bins[ibin,:,:,:] = MASK_4d
                        freq_bins_avg[ibin,:] = MV.sum(MASK_4d)
    
                        data_bins[ibin,:,:,:] = MV.masked_where(MASK_4d==False,data_ocn)
                        data_bins_avg[ibin,:] = cdutil.averager(MV.masked_where(MASK_4d==False,data_ocn),axis='xy',weights='weighted')

                #----------------------------------------------------------
                # define coordinates 
                #----------------------------------------------------------
                lons = cdms.createAxis(data.getLongitude()[:])
                lons.id="lon"
                lons.units="degrees_E"
                lons.designateLongitude()

                lats = cdms.createAxis(data.getLatitude()[:])
                lats.id="lat"
                lats.units="degrees_N"
                lats.designateLatitude()

                regimes = cdms.createAxis([0,1,2])
                regimes.id="regime"
                regimes.units="none"

                time = data.getTime()

                if len(data.shape) == 4:
                    levs = data.getLevel()
                    data.setAxis(0,time)
                    data.setAxis(1,levs)
                    data.setAxis(2,lats)
                    data.setAxis(3,lons)
                else:
                    data.setAxis(0,time)
                    data.setAxis(1,lats)
                    data.setAxis(2,lons)

                freq_bins.setAxis(0,regimes)
                freq_bins.setAxis(1,time)

                freq_bins_avg.setAxis(0,regimes)
                freq_bins_avg.setAxis(1,time)

                data_bins.setAxis(0,regimes)
                data_bins.setAxis(1,time)

                data_bins_avg.setAxis(0,regimes)
                data_bins_avg.setAxis(1,time)

                if len(data.shape) == 4:
                    levs = data.getLevel()
                    freq_bins.setAxis(2,levs)
                    freq_bins.setAxis(3,lats)
                    freq_bins.setAxis(4,lons)

                    data_bins.setAxis(2,levs)
                    data_bins.setAxis(3,lats)
                    data_bins.setAxis(4,lons)

                    freq_bins_avg.setAxis(2,levs)
                    data_bins_avg.setAxis(2,levs)
                else:
                    freq_bins.setAxis(2,lats)
                    freq_bins.setAxis(3,lons)

                    data_bins.setAxis(2,lats)
                    data_bins.setAxis(3,lons)

                #fout = cdms.open('test_global_cloud_regime_'+svar+'_'+case+'.nc','w')
                #=============================================================
                # write into NC files
                #=============================================================
                #DATA = freq_bins
                #DATA.id = svar+'_'+case+'_freq'
                #fout.write(DATA)

                #DATA = data_bins
                #DATA.id = svar+'_'+case
                #fout.write(DATA)

                DATA = freq_bins_avg
                DATA.id = svar+'_'+case+'_freq_avg'
                fout.write(DATA)

                DATA = data_bins_avg
                DATA.id = svar+'_'+case+'_avg'
                fout.write(DATA)

                #fout.close()

        fout.close()

