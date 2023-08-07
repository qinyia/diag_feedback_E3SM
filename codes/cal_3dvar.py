#****************************************************************
#
#    Filename: cal_3d.py
#
#    Author: Yi Qin - qin4@llnl.gov
#    Description: 
#    Input: 
#    Output: 
#    Create: 2022-07-10 17:01:37
#    Last Modified: 2022-07-10 17:01:37
#****************************************************************

#IMPORT STUFF:
#=====================
import numpy as np
import os
import pandas as pd
import xarray as xr
from PlotDefinedFunction import linearregression_nd, area_averager, weighted_annual_mean
import sys
sys.path.append("../")
import cases_lookup as CL

# ============================================================ 
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

###########################################################################
# HELPFUL FUNCTIONS FOLLOW
###########################################################################
def cal_3dvar(direc_data,case_stamp,yearS,yearE,fname1,fname2,outdir,figdir,exp1,exp2):

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
    
    latspc = np.arange(-90,92.5,2.5)
    lonspc = np.arange(1.25,360,2.5)

    var2d = [
    'ts', 'TGCLDLWP','EIS',
    'CLDLOW', 'CLDMED', 'CLDHGH',
    ]

    var3d = [
    'CLOUD', 'CLDLIQ', 
    ]

    var = var2d + var3d

    # =============================================
    # read 2D variables
    # =============================================
    for svarh in var:
        print('svarh=',svarh)

        outfile = outdir+'global_'+svarh+'_'+case_stamp+'.nc'
    
        if os.path.isfile(outfile):
            print('cal_3dvar', case_stamp, svarh, 'output is ready. Please continue. ')
            continue 

        dic_all = {}
        for svar in ['tas',svarh]:
            print(svar)

            if svar in ['SWCRE','LWCRE','netCRE','EIS','PRECT']:
                if svar == 'netCRE':
                    svar_here = ['rsutcs','rsut','rlutcs','rlut']
                elif svar == 'SWCRE':
                    svar_here = ['rsutcs','rsut']
                elif svar == 'LWCRE':
                    svar_here = ['rlutcs','rlut']
                elif svar == 'EIS':
                    svar_here = ['Z3','T','OMEGA','psl','ts']
                elif svar == 'PRECT':
                    svar_here = ['PRECC','PRECL']


                dics = {}
                for svarh in svar_here:
                    f1 = xr.open_dataset(direc_data+fname1+'/'+svarh+'_'+exp1+'_'+yearS_4d+'01-'+yearE_4d+'12.nc')
                    pi_raw1 = f1[svarh]
                    f1.close()
                    f2 = xr.open_dataset(direc_data+fname2+'/'+svarh+'_'+exp2+'_'+yearS_4d+'01-'+yearE_4d+'12.nc')
                    ab_raw1 = f2[svarh]
                    f2.close()

                    if svarh in ['Z3','T']:
                        spec_lev = 700 
                        pi_raw1 = pi_raw1.interp(lev=[spec_lev])[:,0,:,:]
                        ab_raw1 = ab_raw1.interp(lev=[spec_lev])[:,0,:,:]
                    elif svarh in ['OMEGA']:
                        spec_lev = 500 
                        pi_raw1 = pi_raw1.interp(lev=[spec_lev])[:,0,:,:]
                        ab_raw1 = ab_raw1.interp(lev=[spec_lev])[:,0,:,:]

                    pi_raw1_grd = pi_raw1.interp(lat=latspc,lon=lonspc)
                    ab_raw1_grd = ab_raw1.interp(lat=latspc,lon=lonspc)

                    print(pi_raw1_grd.shape, ab_raw1_grd.shape)

                    dics[svarh] = [pi_raw1_grd,ab_raw1_grd]

                # calculate CRE
                if svar == 'netCRE':
                    cre_pi = (dics['rsutcs'][0]-dics['rsut'][0]) + (dics['rlutcs'][0]-dics['rlut'][0])
                    cre_ab = (dics['rsutcs'][1]-dics['rsut'][1]) + (dics['rlutcs'][1]-dics['rlut'][1])
                elif svar == 'SWCRE':
                    cre_pi = dics['rsutcs'][0]-dics['rsut'][0] 
                    cre_ab = dics['rsutcs'][1]-dics['rsut'][1] 
                elif svar == 'LWCRE':
                    cre_pi = dics['rlutcs'][0]-dics['rlut'][0]
                    cre_ab = dics['rlutcs'][1]-dics['rlut'][1]
                elif svar == 'EIS':
                    cre_pi,_ = calc_EIS(dics['ts'][0], dics['psl'][0], dics['T'][0], dics['Z3'][0])
                    cre_ab,_ = calc_EIS(dics['ts'][1], dics['psl'][1], dics['T'][1], dics['Z3'][1])
                elif svar == 'PRECT':
                    cre_pi = dics['PRECC'][0] + dics['PRECL'][0]
                    cre_ab = dics['PRECC'][1] + dics['PRECL'][1]

                pi_raw_grd = xr.DataArray(cre_pi, coords=dics[list(dics.keys())[0]][0].coords)
                ab_raw_grd = xr.DataArray(cre_ab, coords=dics[list(dics.keys())[0]][0].coords)

                # reverse lev direction
                if len(pi_raw_grd.shape) == 4:
                    pi_raw_grd = pi_raw_grd[:,::-1,:,:]
                    ab_raw_grd = ab_raw_grd[:,::-1,:,:]
            else:
                f1 = xr.open_dataset(direc_data+fname1+'/'+svar+'_'+exp1+'_'+yearS_4d+'01-'+yearE_4d+'12.nc')
                pi_raw = f1[svar]
                f1.close()

                f2 = xr.open_dataset(direc_data+fname2+'/'+svar+'_'+exp2+'_'+yearS_4d+'01-'+yearE_4d+'12.nc')
                ab_raw = f2[svar]
                f2.close()

                # reverse lev direction
                if len(pi_raw.shape) == 4:
                    pi_raw = pi_raw[:,::-1,:,:]
                    ab_raw = ab_raw[:,::-1,:,:]

                #----------------------------------------------------------
                # regrid data                 
                #----------------------------------------------------------
                pi_raw_grd = pi_raw.interp(lat=latspc,lon=lonspc)
                ab_raw_grd = ab_raw.interp(lat=latspc,lon=lonspc)

            print('pi_raw_grd.shape = ',pi_raw_grd.shape)
            print('ab_raw_grd.shape = ',ab_raw_grd.shape)

            dic_all[svar+'_ano'] = ab_raw_grd - pi_raw_grd
            dic_all[svar+'_pi'] = pi_raw_grd
            dic_all[svar+'_ab'] = ab_raw_grd

            del(pi_raw_grd, ab_raw_grd)

        # Define new time coordinate
        newtime = pd.date_range(start='1850-01-01', periods=dic_all[list(dic_all.keys())[0]].shape[0], freq='MS')

        #----------------------------------------------------------
        # calculate global mean surface air temperature anomaly  
        #----------------------------------------------------------
        dic_all['tas_ano'] = dic_all['tas_ano'].assign_coords({'time':("time",newtime)})
        anomtas = weighted_annual_mean(dic_all['tas_ano'].time, dic_all['tas_ano']) #(nyears, 90, 144)
        avgdtas = area_averager(anomtas.mean(axis=0)) # (scalar)
        print('avgdtas = ',avgdtas)

        # get time-series of annual mean
        dic_all2 = {}
        for svar in dic_all.keys():
            if '_ano' in svar:
                dic_all[svar] = dic_all[svar].assign_coords({'time':("time",newtime)})
                dic_all2[svar+'_ann'] = weighted_annual_mean(dic_all[svar].time, dic_all[svar])

        tas_ano_ann_gm = area_averager(dic_all2['tas_ano_ann'])

        #----------------------------------------------------------
        # calculate climatological control state 
        #----------------------------------------------------------
        dic_out = {}
        for svar in dic_all.keys():
            if '_pi' in svar:
                dic_out[svar+'_clim'] = dic_all[svar].mean(axis=0)
            if '_ab' in svar:
                dic_out[svar+'_clim'] = dic_all[svar].mean(axis=0)
            elif '_ano' in svar:
                if 'coupled' not in case_stamp:
                    dic_out[svar+'_clim'] = dic_all[svar].mean(axis=0)/avgdtas
                else:
                    print('doing regression...')
                    dic_out[svar+'_clim'],intercept = linearregression_nd(dic_all2[svar+'_ann'], x=np.reshape(tas_ano_ann_gm,(dic_all2[svar+'_ann'].shape[0],1,1)))

        print('dic_out.keys() = ',dic_out.keys())

        #----------------------------------------------------------
        # save data into file     
        #----------------------------------------------------------
        description = '_pi_clim is the climatological control state, _ab_clim is the climatological warming state and xxx_ano_clim is the anomaly normalized by global mean tas anomaly'

        data_vars = {}
        for svar in dic_out.keys():
            print('svar = ', svar)
            tmp = dic_out[svar] 
            if len(tmp.shape) == 2:
                data_vars[svar] = (('lat','lon'),tmp.data)
            elif len(tmp.shape) == 3:
                data_vars[svar] = (('lev','lat','lon'),tmp.data)
                levspc = tmp.lev.values

        if svarh in var3d: 
            da = xr.Dataset(
            data_vars = data_vars, 
            coords = {
            "lev": levspc,
            "lat": latspc,
            "lon": lonspc, 
            },
            attrs=dict(description=description),
            )
        elif svarh in var2d: 
            da = xr.Dataset(
            data_vars = data_vars, 
            coords = {
            "lat": latspc,
            "lon": lonspc, 
            },
            attrs=dict(description=description),
            )

        da.to_netcdf(outdir+outfile)
        
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if __name__ == "__main__":

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
    
    cal_3dvar(direc_data,case_stamp,yearS,yearE,fname1,fname2,outdir,figdir,exp1,exp2)



