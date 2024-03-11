"""
APRP Functions

input: total cloud cover and SW radiative fluxes at TOA and SFC for clear- and all-sky conditions
       -standard CMIP nomenclature: clt,rsdt,rsut,rsutcs,rsds,rsus,rsdscs,rsuscs
       -flag to do forward, backward, or avg of forward / backward calcuations (the default)

output: TOA SW anomalies due to changes in:
        -surface albedo (for all-, clear-, and overcast-sky conditions)
        -clouds (total change and contributions from changing cloud cover, scattering, and absorption)
        -non-cloud atmosphere (e.g., from changes in water vapor, aerosols, ozone)

Equation numbers throughout refer to Taylor et al. (2007)

Reference:
Taylor, K. E. et al. (2007), Estimating shortwave radiative forcing and response in 
    climate models, J. Clim., 20(11), 2530-2543, doi:10.1175/JCLI4143.1.  

Contact: Mark Zelinka (zelinka1@llnl.gov)
"""
 
#IMPORT STUFF:
#=====================
import xarray as xr
import numpy as np
import cases_lookup_ppe as CLP
import os

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
 
###########################################################################
def albedo(c,a_clr,a_oc,mu_clr,mu_cld,ga_clr,ga_cld):

    mu_oc=mu_clr*mu_cld # Eq. 14
    ga_oc=1-(1-ga_clr)*(1-ga_cld) # Eq. 13
    A_clr=(mu_clr*ga_clr) + ((mu_clr*a_clr*(1-ga_clr)**2)/(1-(a_clr*ga_clr))) # Eq. 7
    A_oc= (mu_oc*ga_oc)   + ((mu_oc*a_oc*(1-ga_oc)**2)/(1-(a_oc*ga_oc))) # Eq. 7
    A=(1-c)*A_clr + c*A_oc # Eq. 15

    return A 

    
###########################################################################
def parameters(SWupsfccs,SWdnsfccs,SWdn,SWupcs,SWupsfcoc,SWdnsfcoc,SWupoc):
    
    # clear sky parameters
    a_clr=SWupsfccs/SWdnsfccs # albedo
    Q=SWdnsfccs/SWdn # ratio of incident sfc flux to TOA insolation
    mu_clr=SWupcs/SWdn + Q*(1-a_clr) # Eq. 9
    ga_clr=(mu_clr-Q)/(mu_clr-a_clr*Q) # Eq. 10

    # overcast parameters
    a_oc=SWupsfcoc/SWdnsfcoc # albedo
    Q=SWdnsfcoc/SWdn # ratio of incident sfc flux to TOA insolation
    mu_oc=SWupoc/SWdn + Q*(1-a_oc) # Eq. 9
    ga_oc=(mu_oc-Q)/(mu_oc-a_oc*Q) # Eq. 10

    # cloud parameters
    mu_cld=mu_oc/mu_clr  # Eq. 14 sometimes this is greater than 1??
    ga_cld=(ga_oc-1)/(1-ga_clr)+1  # Eq. 13

    return (a_clr,mu_clr,ga_clr,a_oc,mu_cld,ga_cld) 

###########################################################################
def APRP(CTL,PERT,flag=''):

    # Get stuff out of the dictionary, give control values suffix of 1, perturbed suffix of 2:
    clt1,clt2 = CTL['clt'],PERT['clt']
    rsdt1,rsdt2 = CTL['rsdt'],PERT['rsdt']
    rsut1,rsut2 = CTL['rsut'],PERT['rsut']
    rsutcs1,rsutcs2 = CTL['rsutcs'],PERT['rsutcs']
    rsds1,rsds2 = CTL['rsds'],PERT['rsds']
    rsus1,rsus2 = CTL['rsus'],PERT['rsus']
    rsdscs1,rsdscs2 = CTL['rsdscs'],PERT['rsdscs']
    rsuscs1,rsuscs2 = CTL['rsuscs'],PERT['rsuscs']
    
    rlnt1,rlnt2 = CTL['FLNT'],PERT['FLNT']
    rlntc1,rlntc2 = CTL['FLNTC'],PERT['FLNTC']
    
    # Make sure the cld fractions are expressed as fraction and not percent
    if np.max(clt1)>1.:
        clt1=clt1/100.
    if np.max(clt2)>1.:
        clt2=clt2/100.

    # Derive overcast conditions
    rsutoc1=(1/clt1)*(rsut1-(1-clt1)*rsutcs1)
    rsdsoc1=(1/clt1)*(rsds1-(1-clt1)*rsdscs1)
    rsusoc1=(1/clt1)*(rsus1-(1-clt1)*rsuscs1)
    rsutoc2=(1/clt2)*(rsut2-(1-clt2)*rsutcs2)
    rsdsoc2=(1/clt2)*(rsds2-(1-clt2)*rsdscs2)
    rsusoc2=(1/clt2)*(rsus2-(1-clt2)*rsuscs2)

    # Mask these where values are unphysical 
    rsdsoc1=xr.where(rsdsoc1 > rsds1,np.nan,rsdsoc1)   
    rsusoc1=xr.where(rsusoc1 > rsus1,np.nan,rsusoc1)
    rsutoc1=xr.where(rsutoc1 < 0,np.nan,rsutoc1)
    rsdsoc1=xr.where(rsdsoc1 < 0,np.nan,rsdsoc1)
    rsusoc1=xr.where(rsusoc1 < 0,np.nan,rsusoc1)

    rsdsoc2=xr.where(rsdsoc2 > rsds2,np.nan,rsdsoc2)
    rsusoc2=xr.where(rsusoc2 > rsus2,np.nan,rsusoc2)
    rsutoc2=xr.where(rsutoc2 < 0,np.nan,rsutoc2)
    rsdsoc2=xr.where(rsdsoc2 < 0,np.nan,rsdsoc2)
    rsusoc2=xr.where(rsusoc2 < 0,np.nan,rsusoc2)    

    ## NOW THE FORMAL APRP CALCULATIONS:
    a_clr1,mu_clr1,ga_clr1,a_oc1,mu_cld1,ga_cld1 = \
        parameters(rsuscs1,rsdscs1,rsdt1,rsutcs1,rsusoc1,rsdsoc1,rsutoc1) # control
    a_clr2,mu_clr2,ga_clr2,a_oc2,mu_cld2,ga_cld2 = \
        parameters(rsuscs2,rsdscs2,rsdt2,rsutcs2,rsusoc2,rsdsoc2,rsutoc2) # perturbed    

    ## Taylor et al. (2007) Eqn. 12b:
    A_1=albedo(clt1,a_clr1,a_oc1,mu_clr1,mu_cld1,ga_clr1,ga_cld1)
    A_2=albedo(clt2,a_clr2,a_oc2,mu_clr2,mu_cld2,ga_clr2,ga_cld2)

   
    # Forward PRP calculation:
    dA_amt_cld_fwd=     albedo(clt2,a_clr1,a_oc1,mu_clr1,mu_cld1,ga_clr1,ga_cld1)-A_1 # 16b.3
    dA_a_clr_fwd=       albedo(clt1,a_clr2,a_oc1,mu_clr1,mu_cld1,ga_clr1,ga_cld1)-A_1 # 16a.1
    dA_a_oc_fwd=        albedo(clt1,a_clr1,a_oc2,mu_clr1,mu_cld1,ga_clr1,ga_cld1)-A_1 # 16a.2
    dA_abs_noncld_fwd=  albedo(clt1,a_clr1,a_oc1,mu_clr2,mu_cld1,ga_clr1,ga_cld1)-A_1 # 16c.1
    dA_abs_cld_fwd=     albedo(clt1,a_clr1,a_oc1,mu_clr1,mu_cld2,ga_clr1,ga_cld1)-A_1 # 16b.1
    dA_scat_noncld_fwd= albedo(clt1,a_clr1,a_oc1,mu_clr1,mu_cld1,ga_clr2,ga_cld1)-A_1 # 16c.2
    dA_scat_cld_fwd=    albedo(clt1,a_clr1,a_oc1,mu_clr1,mu_cld1,ga_clr1,ga_cld2)-A_1 # 16b.2

    # Backward PRP calculation:
    dA_amt_cld_bwd=    A_2-albedo(clt1,a_clr2,a_oc2,mu_clr2,mu_cld2,ga_clr2,ga_cld2)
    dA_a_clr_bwd=      A_2-albedo(clt2,a_clr1,a_oc2,mu_clr2,mu_cld2,ga_clr2,ga_cld2)
    dA_a_oc_bwd=       A_2-albedo(clt2,a_clr2,a_oc1,mu_clr2,mu_cld2,ga_clr2,ga_cld2)
    dA_abs_noncld_bwd= A_2-albedo(clt2,a_clr2,a_oc2,mu_clr1,mu_cld2,ga_clr2,ga_cld2)
    dA_abs_cld_bwd=    A_2-albedo(clt2,a_clr2,a_oc2,mu_clr2,mu_cld1,ga_clr2,ga_cld2)
    dA_scat_noncld_bwd=A_2-albedo(clt2,a_clr2,a_oc2,mu_clr2,mu_cld2,ga_clr1,ga_cld2)
    dA_scat_cld_bwd=   A_2-albedo(clt2,a_clr2,a_oc2,mu_clr2,mu_cld2,ga_clr2,ga_cld1)

    if flag=='': # do forward and backward PRP (default)
        dA_amt_cld =    0.5*(dA_amt_cld_fwd + dA_amt_cld_bwd)
        dA_a_clr =      0.5*(dA_a_clr_fwd + dA_a_clr_bwd) 
        dA_a_oc =       0.5*(dA_a_oc_fwd + dA_a_oc_bwd)
        dA_abs_noncld = 0.5*(dA_abs_noncld_fwd + dA_abs_noncld_bwd)
        dA_abs_cld =    0.5*(dA_abs_cld_fwd + dA_abs_cld_bwd)
        dA_scat_noncld =0.5*(dA_scat_noncld_fwd + dA_scat_noncld_bwd)
        dA_scat_cld =   0.5*(dA_scat_cld_fwd + dA_scat_cld_bwd)        
        RSDT = 0.5*(rsdt1+rsdt2)
    elif flag=='forward': # do forward-only PRP  
        dA_amt_cld = dA_amt_cld_fwd
        dA_a_clr = dA_a_clr_fwd
        dA_a_oc = dA_a_oc_fwd
        dA_abs_noncld = dA_abs_noncld_fwd
        dA_abs_cld = dA_abs_cld_fwd
        dA_scat_noncld = dA_scat_noncld_fwd
        dA_scat_cld = dA_scat_cld_fwd 
        RSDT = rsdt1
    elif flag=='backward': # do backward-only PRP
        dA_amt_cld = dA_amt_cld_bwd
        dA_a_clr = dA_a_clr_bwd
        dA_a_oc = dA_a_oc_bwd
        dA_abs_noncld = dA_abs_noncld_bwd
        dA_abs_cld = dA_abs_cld_bwd
        dA_scat_noncld = dA_scat_noncld_bwd
        dA_scat_cld = dA_scat_cld_bwd
        RSDT = rsdt2

        
    # if the cld fraction is less than 2%, set fields to be zero
    dA_amt_cld=xr.where(clt1<0.02,0.,dA_amt_cld)
    dA_amt_cld=xr.where(clt2<0.02,0.,dA_amt_cld)
    dA_a_oc=xr.where(clt1<0.02,0.,dA_a_oc)
    dA_a_oc=xr.where(clt2<0.02,0.,dA_a_oc)
    dA_abs_cld=xr.where(clt1<0.02,0.,dA_abs_cld)
    dA_abs_cld=xr.where(clt2<0.02,0.,dA_abs_cld)
    dA_scat_cld=xr.where(clt1<0.02,0.,dA_scat_cld)
    dA_scat_cld=xr.where(clt2<0.02,0.,dA_scat_cld)

    dA_a=dA_a_clr+dA_a_oc
    dA_cld=dA_abs_cld+dA_scat_cld+dA_amt_cld
    dA_noncld=dA_abs_noncld+dA_scat_noncld

    ## TOA SW Anomalies due to Surface Albedo Anomalies
    sfc_alb=-dA_a*RSDT
    sfc_alb_clr=-dA_a_clr*RSDT   
    sfc_alb_oc=-dA_a_oc*RSDT

    ## TOA SW Anomalies due to Cloud Anomalies
    cld=-dA_cld*RSDT
    cld_amt=-dA_amt_cld*RSDT    
    cld_scat=-dA_scat_cld*RSDT
    cld_abs=-dA_abs_cld*RSDT

    ## TOA SW Anomalies due to Non-cloud Anomalies
    noncld=-dA_noncld*RSDT
    noncld_scat=-dA_scat_noncld*RSDT  
    noncld_abs=-dA_abs_noncld*RSDT
    
    # set fields to zero when incoming solar radiation is zero
    sfc_alb=xr.where(RSDT<0.1,0.,sfc_alb)
    cld=xr.where(RSDT<0.1,0.,cld)
    noncld=xr.where(RSDT<0.1,0.,noncld)
    sfc_alb_clr=xr.where(RSDT<0.1,0.,sfc_alb_clr)
    sfc_alb_oc=xr.where(RSDT<0.1,0.,sfc_alb_oc)
    cld_amt=xr.where(RSDT<0.1,0.,cld_amt)
    cld_scat=xr.where(RSDT<0.1,0.,cld_scat)
    cld_abs=xr.where(RSDT<0.1,0.,cld_abs)
    noncld_scat=xr.where(RSDT<0.1,0.,noncld_scat)
    noncld_abs=xr.where(RSDT<0.1,0.,noncld_abs)

    coords = rlnt1.coords
    cldlw = xr.DataArray((rlnt2-rlntc2) - (rlnt1-rlntc1), coords=coords)
    cldlw *= -1 
    cldnet = xr.DataArray(cld + cldlw, coords=coords)
    noncldlw = xr.DataArray(rlntc2 - rlntc1, coords=coords)
    noncldlw *= -1 
    noncldnet = xr.DataArray(noncld + noncldlw, coords=coords)
    print('noncldnet.shape=',noncldnet.shape)

    # store in a dataset:
    TIME = sfc_alb.time
    LAT = sfc_alb.lat
    LON = sfc_alb.lon
    DS = xr.Dataset(
    {
        'tas':(('time','lat','lon'), PERT['tas'].data-CTL['tas'].data), 
        'sfc_alb':(('time','lat','lon'),sfc_alb.data),
        'sfc_alb_clr':(('time','lat','lon'),sfc_alb_clr.data),
        'sfc_alb_oc':(('time','lat','lon'),sfc_alb_oc.data),
        'cld':(('time','lat','lon'),cld.data),
        'cld_amt':(('time','lat','lon'),cld_amt.data),
        'cld_scat':(('time','lat','lon'),cld_scat.data),
        'cld_abs':(('time','lat','lon'),cld_abs.data),
        'noncld':(('time','lat','lon'),noncld.data),
        'noncld_scat':(('time','lat','lon'),noncld_scat.data),
        'noncld_abs':(('time','lat','lon'),noncld_abs.data),
        'cldlw':(('time','lat','lon'),cldlw.data),
        'cldnet':(('time','lat','lon'),cldnet.data),
        'noncldlw':(('time','lat','lon'),noncldlw.data),
        'noncldnet':(('time','lat','lon'),noncldnet.data),
    },
    coords={'time': TIME,'lat': LAT,'lon': LON},
    ) 
    DS.lat.attrs["axis"] = "Y"
    DS.lon.attrs["axis"] = "X"
    # output = DS.bounds.add_missing_bounds()
    output = DS

    return output


def cal_APRP(direc_data,case_stamp,yearS,yearE,fname1,fname2,outdir,figdir,exp1,exp2):

    print('Hello')

    # Check whether file exists
    # -------------------------------------------------------
    outfile = "APRP_"+case_stamp+".nc"
    if os.path.isfile(outdir+outfile):
        print('APRP is already there.')
        return

    # Format start and end years/months
    # -------------------------------------------------------
    yearS_4d = "{:04d}".format(yearS)
    yearE_4d = "{:04d}".format(yearE)
    nyears = yearE - yearS + 1

    monS = 1
    monE = 12
    monS_2d='{:02d}'.format(monS)
    monE_2d='{:02d}'.format(monE) 

    # Data directory
    # -------------------------------------------------------
    direc_data2 = direc_data+'/'+fname2+'/'

    # Mapping variable names
    # -------------------------------------------------------
    variables = ['clt','rsdt','rsut','rsutcs','rsds','rsus','rsdscs','rsuscs','tas','FLNT', 'FLNTC']
    variables_E3SM = ['CLDTOT', 'rsdt','rsut','rsutcs','rsds','rsus','rsdscs','rsuscs','tas','rlut', 'rlutcs'] 


    # Read data
    # -------------------------------------------------------
    casenames = [
        fname1, 
        fname2, 
    ]
    casetags = [
        'PD_CTL',
        'P4K_CTL',
    ]

    exps = [
        exp1,
        exp2,
    ]
    
    DATA={}
    for iexp,exp in enumerate(casetags):
        casenameh = casenames[iexp]
    
        DATA[exp]={}
    
        for ivar,svar in enumerate(variables):

            fnamef = direc_data+'/'+casenameh+'/'+variables_E3SM[ivar]+'_'+exps[iexp]+'_'+yearS_4d+'01-'+yearE_4d+'12.nc'
            f1 = xr.open_dataset(fnamef)
            data = f1[variables_E3SM[ivar]]
            f1.close()
    
            DATA[exp][svar] = data 
    
    result = APRP(DATA['PD_CTL'],DATA['P4K_CTL'])

    # Normalized by tas
    # -------------------------------------------------------
    avgdtas = area_averager(result['tas']).mean(dim='time')
    print('avgdtas=',avgdtas)

    outdata = result.mean(dim='time')/avgdtas

    # output time averaged dataset
    # -------------------------------------------------------
    outdata.to_netcdf(outdir+outfile)


    # Calculate and output global average 
    # -------------------------------------------------------
    output1_avg_json = {}
    f1 = xr.open_dataset(outdir+outfile)
        
    # get global average
    output1_avg_json[case_stamp] = {} 
    for key in list(f1.data_vars):
        data = f1[key]
        output1_avg_json[case_stamp][key] = area_averager(data).data.tolist()
        
    print(output1_avg_json.keys())
    
    # save global average file 
    #AR.save_json(output1_avg_json,outdir+'APRP_'+case_stamp+'_gm.json')
    
    f1.close()


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if __name__ == "__main__":

    direc_data = '/compyfs/qiny108/diag_feedback_E3SM_postdata_ppe/'
    case_stamp = 'BASE'
    yearS = 2010
    yearE = 2014
    fname1 = CLP.get_lutable(case_stamp,'PD_FR')
    fname2 = CLP.get_lutable(case_stamp,'P4K_FR')
    outdir = './'
    figdir = './'
    exp1 = 'FC5'
    exp2 = 'FC5_4K'
    
    cal_APRP(direc_data,case_stamp,yearS,yearE,fname1,fname2,outdir,figdir,exp1,exp2)

