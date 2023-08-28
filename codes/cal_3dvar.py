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
from get_mip_data import read_mip_data,read_amip_data,read_pickle,write_pickle,read_e3sm_data

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
def cal_3dvar(direc_data,case_stamp,yearS,yearE,fname1,fname2,outdir,figdir):
   
    latspc = np.arange(-90,92.5,2.5)
    lonspc = np.arange(1.25,360,2.5)

    var2d = [
    'ts', 'TGCLDLWP','TGCLDIWP',
    'EIS',
    'CLDLOW', 'CLDMED', 'CLDHGH','CLDTOT',
    ]

    var3d = [
    'CLOUD', 'CLDLIQ', 'CLDICE',
    ]

    var = var2d + var3d

    # =============================================
    # read 2D variables
    # =============================================
    for svarh in var:
        print('svarh = ',svarh)

        outfile = outdir+'global_'+svarh+'_'+case_stamp+'.nc'
    
        if os.path.isfile(outfile):
            print('cal_3dvar', case_stamp, svarh, 'output is ready. Please continue. ')
            continue 

        if svarh == 'netCRE':
            svar_here = ['rsutcs','rsut','rlutcs','rlut']
        elif svarh == 'SWCRE':
            svar_here = ['rsutcs','rsut']
        elif svarh == 'LWCRE':
            svar_here = ['rlutcs','rlut']
        elif svarh == 'EIS':
            svar_here = ['Z3','T','OMEGA','psl','ts']
        elif svarh == 'PRECT':
            svar_here = ['PRECC','PRECL']
        else:
            svar_here = [svarh]

        svar_here += ['tas']

        dics = read_e3sm_data(svar_here,direc_data,case_stamp,yearS,yearE,fname1,fname2)

        # Regrid 
        for key in dics.keys():
            dics[key] = dics[key].interp(lat=latspc,lon=lonspc)

        # calculate CRE
        if svarh == 'netCRE':
            cre_pi = (dics['rsutcs_pi']-dics['rsut_pi']) + (dics['rlutcs_pi']-dics['rlut_pi'])
            cre_ab = (dics['rsutcs_ab']-dics['rsut_ab']) + (dics['rlutcs_ab']-dics['rlut_ab'])
        elif svarh == 'SWCRE':
            cre_pi = dics['rsutcs_pi']-dics['rsut_pi'] 
            cre_ab = dics['rsutcs_ab']-dics['rsut_ab'] 
        elif svarh == 'LWCRE':
            cre_pi = dics['rlutcs_pi']-dics['rlut_pi']
            cre_ab = dics['rlutcs_ab']-dics['rlut_ab']
        elif svarh == 'EIS':
            cre_pi,_ = calc_EIS(dics['ts_pi'], dics['psl_pi'], dics['T_pi'].interp(lev=[700]), dics['Z3_pi'].interp(lev=[700]))
            cre_ab,_ = calc_EIS(dics['ts_ab'], dics['psl_ab'], dics['T_ab'].interp(lev=[700]), dics['Z3_ab'].interp(lev=[700]))
        elif svarh == 'PRECT':
            cre_pi = dics['PRECC_pi'] + dics['PRECL_pi']
            cre_ab = dics['PRECC_ab'] + dics['PRECL_ab']
        else:
            cre_pi = dics[svarh+'_pi']
            cre_ab = dics[svarh+'_ab']

        pi_raw_grd = xr.DataArray(cre_pi, coords=dics['tas_pi'].coords)
        ab_raw_grd = xr.DataArray(cre_ab, coords=dics['tas_pi'].coords)

        # reverse lev direction
        if len(pi_raw_grd.shape) == 4:
            pi_raw_grd = pi_raw_grd[:,::-1,:,:]
            ab_raw_grd = ab_raw_grd[:,::-1,:,:]

        dic_all = {}
        dic_all[svarh+'_pi'] = pi_raw_grd
        dic_all[svarh+'_ab'] = ab_raw_grd 
        dic_all[svarh+'_ano'] = ab_raw_grd - pi_raw_grd

        dic_all['tas_pi'] = dics['tas_pi']
        dic_all['tas_ab'] = dics['tas_ab']
        dic_all['tas_ano'] = dics['tas_ano']

        del(pi_raw_grd, ab_raw_grd)

        # Define new time coordinate
        newtime = pd.date_range(start='1850-01-01', periods=dic_all[list(dic_all.keys())[0]].shape[0], freq='MS')

        #----------------------------------------------------------
        # calculate global mean surface air temperature anomaly  
        #----------------------------------------------------------
        dic_all['tas_ano'] = dic_all['tas_ano'].assign_coords({'time':("time",newtime)})
        anomtas = weighted_annual_mean(dic_all['tas_ano'].time, dic_all['tas_ano']) #(nyears, 90, 144)
        avgdtas = area_averager(anomtas.mean(axis=0)) # (scalar)
        print('avgdtas = ',avgdtas.values)

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

        print('dic_out.keys() = ',list(dic_out.keys()))

        #----------------------------------------------------------
        # save data into file     
        #----------------------------------------------------------
        description = '_pi_clim is the climatological control state, _ab_clim is the climatological warming state and xxx_ano_clim is the anomaly normalized by global mean tas anomaly'

        data_vars = {}
        for svar in dic_out.keys():
            print('Saving svar = ', svar)
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

        da.to_netcdf(outfile)
        
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if __name__ == "__main__":

    #direc_data = '/compyfs/qiny108/diag_feedback_E3SM_postdata/'
    #direc_data = '/compyfs/qiny108/colla/diag_feedback_E3SM_postdata/'
    direc_data = '/p/user_pub/climate_work/qin4/From_Compy/compyfs_dir/colla/diag_feedback_E3SM_postdata/'

    case_stamp = 'v2test'
    yearS = 2
    yearE = 3
    fname1,_,_ = CL.get_lutable(case_stamp,'amip')
    fname2,_,_ = CL.get_lutable(case_stamp,'amip4K')
    outdir = './'
    figdir = './'
    exp1 = 'FC5'
    exp2 = 'FC5_4K'
    
    cal_3dvar(direc_data,case_stamp,yearS,yearE,fname1,fname2,outdir,figdir)

