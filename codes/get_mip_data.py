import xarray as xr 
import pickle
import pandas as pd
import numpy as np
import glob 

# Standard vertical levels  [Pa]
Dlevs = [100000, 92500, 85000, 70000, 60000, 50000, 40000, 30000, 25000, 20000, 15000, 10000, 7000, 5000, 3000, 2000, 1000, 500, 100]

# ==================================================================================
def read_amip_data(filename,var,tslice, lev=None):
    '''
    Refer to the same function from Mark Zelinka's assessed-cloud-fbks package
    '''
    # load in cmip data using the appropriate function for the experiment/mip
    print(var, filename[var][0]+"*.nc")
    try:
        f=xr.open_dataset(filename[var][0]+"*.nc")
    except:
        f=xr.open_mfdataset(filename[var][0]+"*.nc",data_vars="minimal")

    if lev:
        data = f[var].sel(time=tslice,level=lev).squeeze()
    else:
        data = f[var].sel(time=tslice).squeeze()

    # Check level coordinate
    if len(data.shape) == 4: # and 'plev' not in list(data.coords):  # (time,lev,lat,lon)
        if 'plev' in list(data.coords): # standard levels - not required further calculation
            return data 
        else:
            varlist = list(f.data_vars)

        if ('a' in varlist or 'ap' in varlist) and 'b' in varlist : # hybrid sigma coordinate
            try:
                ap = f['ap']
            except: 
                ap = f['a']
            bp = f['b']

            if 'p0' in varlist:
                p0 = f['p0']

            # Read ps 
            try: 
                fps = xr.open_dataset(filename['ps'][0]+'*.nc')
            except:
                fps = xr.open_mfdataset(filename['ps'][0]+'*.nc')
            ps = fps['ps'].sel(time=tslice)
            fps.close()


            coords={'time':ps.coords['time'],'lev':ap.coords['lev'],'lat':ps.coords['lat'],'lon':ps.coords['lon']}

            ps_4d = xr.DataArray(np.transpose(np.tile(ps,(ap.shape[0],1,1,1)),(1,0,2,3)), coords=coords)
            ap_4d = xr.DataArray(np.transpose(np.tile(ap,(ps.shape[0],ps.shape[1],ps.shape[2],1)),(0,3,1,2)),coords=coords)
            bp_4d = xr.DataArray(np.transpose(np.tile(bp,(ps.shape[0],ps.shape[1],ps.shape[2],1)),(0,3,1,2)),coords=coords)

            if 'p0' in varlist: 
                level_true = ap_4d * p0.values + bp_4d * ps_4d
            else:
                level_true = ap_4d + bp_4d * ps_4d

            print('level_true=',level_true[0,:,50,50].values)
            print('data=',data[0,:,50,50].values)

            data_new = logLinearInterpolation(data,level_true,levels=Dlevs,axis='lev')

            print('level_new=',data_new.lev.values)
            print('data_new=',data_new[0,:,50,50].values)

            #np.savez('test_LCF.npz',level_true=level_true[0,:,50,50].values, level_new=data_new.lev.values,data=data[0,:,50,50].values,data_new=data_new[0,:,50,50].values)
        else:
            data_new = data

    else:
        data_new = data 

    f.close()

    return data_new

# ==================================================================================
def read_mip_data(Vars,filenames,tslice):
    '''
    Refer to get_CRK_data function from Mark Zelinka's assessed-cld-fbks package
    '''
    
    # Load in monthly mean from control and perturbed simulation
    dic_all = {}
    for svar in Vars:
        pi_raw = read_amip_data(filenames['amip'],svar,tslice)
        ab_raw = read_amip_data(filenames['amip-p4K'],svar,tslice)

        dic_all[svar+'_ano'] = xr.DataArray(ab_raw - pi_raw, coords=pi_raw.coords)
        dic_all[svar+'_pi'] = pi_raw
        dic_all[svar+'_ab'] = ab_raw

        del(pi_raw, ab_raw)
                   
    return dic_all

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


# ========================================================================================================
def read_e3sm_data(Vars,direc_data,case_stamp,yearS,yearE,fname1,fname2):
    '''
    Read E3SM data
    ''' 

    exp1 = 'FC5'
    exp2 = 'FC5_4K'

    direc_data1 = direc_data+'/'+fname1+'/'
    direc_data2 = direc_data+'/'+fname2+'/'
    
    nyears = yearE - yearS + 1

    monS=1
    monE=12

    yrS_4d='{:04d}'.format(yearS)
    yrE_4d='{:04d}'.format(yearE)
    monS_2d='{:02d}'.format(monS)
    monE_2d='{:02d}'.format(monE)

    dics_invar = {}
    for svar in Vars:
        print('Start reading '+svar+'...')
        if svar == 'clisccp':
            svar_in = 'FISCCP1_COSP'
        elif svar == 'ta':
            svar_in = 'T'
        elif svar == 'hus':
            svar_in = 'Q'
        elif svar == 'clw':
            svar_in = 'CLDLIQ'
        elif svar == 'cli':
            svar_in = 'CLDICE'
        else:
            svar_in = svar
        
        f=xr.open_dataset(direc_data1+svar_in+'_'+exp1+'_'+yrS_4d+monS_2d+'-'+yrE_4d+monE_2d+'.nc')
        data_pi=f[svar_in]
        f.close()
        f=xr.open_dataset(direc_data2+svar_in+'_'+exp2+'_'+yrS_4d+monS_2d+'-'+yrE_4d+monE_2d+'.nc')
        data_ab=f[svar_in]
        f.close()

        # reverse lev direction
        if svar in ['ta','hus']:
            data_pi = data_pi[:,::-1,:,:]
            data_ab = data_ab[:,::-1,:,:]

        if svar == 'clisccp':
            # Make sure clisccp is in percent  
            sumclisccp1 = data_pi.sum(axis=2).sum(axis=1)
            sumclisccp2 = data_ab.sum(axis=2).sum(axis=1)
            if np.max(sumclisccp1) <= 1.:
                print('Changing clisccp in percent...')
                data_pi = data_pi*100.        
            if np.max(sumclisccp2) <= 1.:
                data_ab = data_ab*100.

            data_pi = data_pi.transpose('time','cosp_tau','cosp_prs','lat','lon')
            data_ab = data_ab.transpose('time','cosp_tau','cosp_prs','lat','lon')
    
        # Compute clisccp anomalies
        data_ano = xr.DataArray(data_ab - data_pi, coords=data_pi.coords)

        dics_invar[svar+'_pi'] = data_pi
        dics_invar[svar+'_ab'] = data_ab
        dics_invar[svar+'_ano'] = data_ano

    return dics_invar 

# ========================================================================================================
def read_e3smdiag_data(Vars,direc_data,case_stamp,yearS,yearE,fname1,fname2,grd_info,num_years_per_file):
    '''
    Read E3SM data from e3sm_diag output
    ''' 

    direc_data1 = direc_data+'/'+fname1+'/post/atm/'+grd_info+'/ts/monthly/'+str(num_years_per_file)+'yr/'
    direc_data2 = direc_data+'/'+fname2+'/post/atm/'+grd_info+'/ts/monthly/'+str(num_years_per_file)+'yr/'
    print('direc_data1 = ', direc_data1)
    print('direc_data2 = ', direc_data2)
    print()

    
    dics_invar = {}
    for svar in Vars:
        print('Start reading '+svar+'...')
        if svar == 'clisccp':
            svar_in = 'FISCCP1_COSP'
        elif svar == 'ta':
            svar_in = 'T'
        elif svar == 'hus':
            svar_in = 'Q'
        elif svar == 'clw':
            svar_in = 'CLDLIQ'
        elif svar == 'cli':
            svar_in = 'CLDICE'
        elif svar == 'tas':
            svar_in = 'TREFHT'
        elif svar == 'rlut':
            svar_in = 'FLUT'
        elif svar == 'rlutcs':
            svar_in = 'FLUTC'
        elif svar == 'rsut':
            svar_in = 'rsut'
        elif svar == 'rsutcs':
            svar_in = 'rsutcs'
        elif svar == 'rsds':
            svar_in = 'FSDS'
        elif svar == 'rsdscs':
            svar_in = 'FSDSC'
        elif svar == 'rsus':
            svar_in = 'rsus'
        elif svar == 'rsuscs':
            svar_in = 'rsuscs'
        elif svar == 'ps':
            svar_in = 'PS'
        elif svar == 'psl':
            svar_in = 'PSL'
        elif svar == 'rsdt':
            svar_in = 'SOLIN'
        elif svar == 'ts':
            svar_in = 'TS'
        else:
            svar_in = svar
        
        if len(glob.glob1(direc_data1,svar_in+"_*.nc")) == 0: 
            print('No direct output for '+svar_in+':')
            if svar_in == 'rsutcs':
                print('    Doing calculation for '+svar_in)
                svar_tmp = 'SOLIN'
                f=xr.open_mfdataset(direc_data1+svar_tmp+'_*.nc')
                data1_pi=f[svar_tmp]
                f.close()
                f=xr.open_mfdataset(direc_data2+svar_tmp+'_*.nc')
                data1_ab=f[svar_tmp]
                f.close()

                svar_tmp = 'FSNTC'
                f=xr.open_mfdataset(direc_data1+svar_tmp+'_*.nc')
                data2_pi=f[svar_tmp]
                f.close()
                f=xr.open_mfdataset(direc_data2+svar_tmp+'_*.nc')
                data2_ab=f[svar_tmp]
                f.close()

                data_pi = xr.DataArray(data1_pi - data2_pi, coords=data1_pi.coords)
                data_ab = xr.DataArray(data1_ab - data2_ab, coords=data1_ab.coords)
            elif svar_in == 'rsut':
                print('    Doing calculation for '+svar_in)
                svar_tmp = 'SOLIN'
                f=xr.open_mfdataset(direc_data1+svar_tmp+'_*.nc')
                data1_pi=f[svar_tmp]
                f.close()
                f=xr.open_mfdataset(direc_data2+svar_tmp+'_*.nc')
                data1_ab=f[svar_tmp]
                f.close()

                svar_tmp = 'FSNT'
                f=xr.open_mfdataset(direc_data1+svar_tmp+'_*.nc')
                data2_pi=f[svar_tmp]
                f.close()
                f=xr.open_mfdataset(direc_data2+svar_tmp+'_*.nc')
                data2_ab=f[svar_tmp]
                f.close()

                data_pi = xr.DataArray(data1_pi - data2_pi, coords=data1_pi.coords)
                data_ab = xr.DataArray(data1_ab - data2_ab, coords=data1_ab.coords)

            elif svar_in == 'rsuscs':
                print('    Doing calculation for '+svar_in)

                svar_tmp = 'FSDSC'
                f=xr.open_mfdataset(direc_data1+svar_tmp+'_*.nc')
                data1_pi=f[svar_tmp]
                f.close()
                f=xr.open_mfdataset(direc_data2+svar_tmp+'_*.nc')
                data1_ab=f[svar_tmp]
                f.close()

                svar_tmp = 'FSNSC'
                f=xr.open_mfdataset(direc_data1+svar_tmp+'_*.nc')
                data2_pi=f[svar_tmp]
                f.close()
                f=xr.open_mfdataset(direc_data2+svar_tmp+'_*.nc')
                data2_ab=f[svar_tmp]
                f.close()

                data_pi = xr.DataArray(data1_pi - data2_pi, coords=data1_pi.coords)
                data_ab = xr.DataArray(data1_ab - data2_ab, coords=data1_ab.coords)

            elif svar_in == 'rsus':
                print('    Doing calculation for '+svar_in)

                svar_tmp = 'FSDS'
                f=xr.open_mfdataset(direc_data1+svar_tmp+'_*.nc')
                data1_pi=f[svar_tmp]
                f.close()
                f=xr.open_mfdataset(direc_data2+svar_tmp+'_*.nc')
                data1_ab=f[svar_tmp]
                f.close()

                svar_tmp = 'FSNS'
                f=xr.open_mfdataset(direc_data1+svar_tmp+'_*.nc')
                data2_pi=f[svar_tmp]
                f.close()
                f=xr.open_mfdataset(direc_data2+svar_tmp+'_*.nc')
                data2_ab=f[svar_tmp]
                f.close()

                data_pi = xr.DataArray(data1_pi - data2_pi, coords=data1_pi.coords)
                data_ab = xr.DataArray(data1_ab - data2_ab, coords=data1_ab.coords)

            else:
                print('ERROR: We should do the calculation for', svar_in, '. Please check.')
                exit()

        else:
            f=xr.open_mfdataset(direc_data1+svar_in+'_*.nc')
            data_pi=f[svar_in]
            f.close()
            f=xr.open_mfdataset(direc_data2+svar_in+'_*.nc')
            data_ab=f[svar_in]
            f.close()

        # reverse lev direction
        if svar in ['ta','hus']:
            data_pi = data_pi[:,::-1,:,:]
            data_ab = data_ab[:,::-1,:,:]

        if svar == 'clisccp':
            # Make sure clisccp is in percent  
            sumclisccp1 = data_pi.sum(axis=2).sum(axis=1)
            sumclisccp2 = data_ab.sum(axis=2).sum(axis=1)
            if np.max(sumclisccp1) <= 1.:
                print('Changing clisccp in percent...')
                data_pi = data_pi*100.        
            if np.max(sumclisccp2) <= 1.:
                data_ab = data_ab*100.

            data_pi = data_pi.transpose('time','cosp_tau','cosp_prs','lat','lon')
            data_ab = data_ab.transpose('time','cosp_tau','cosp_prs','lat','lon')
    
        # Compute clisccp anomalies
        data_ano = xr.DataArray(data_ab - data_pi, coords=data_pi.coords)

        dics_invar[svar+'_pi'] = data_pi
        dics_invar[svar+'_ab'] = data_ab
        dics_invar[svar+'_ano'] = data_ano

    return dics_invar 


# ============================================================================================
def logLinearInterpolation(
        A, P, levels=[100000, 92500, 85000, 70000, 60000, 50000, 40000,
                      30000, 25000, 20000, 15000, 10000, 7000, 5000,
                      3000, 2000, 1000], axis='z'):
    """
    Converted from cdutil.logLinearRegression by Yi Qin (yi.qin@pnnl.gov).
    
    Description: Log-linear interpolation to convert a field from sigma levels to pressure levels. 
                 Values below surface are masked.

    :param A: array on sigma levels
    :param P: pressure field from TOP (level 0) to BOTTOM (last level)
    :param levels: pressure levels to interplate to (same units as P), default levels are:
            [100000, 92500, 85000, 70000, 60000, 50000, 40000, 30000, 25000, 20000, 15000, 10000, 7000, 5000,
            3000, 2000, 1000]
    :type levels: list
    :param axis: axis over which to do the linear interpolation
    :type axis: str

    .. note::
        P and levels must have same units
    :returns: array on pressure levels (levels)
    :Example:
        .. doctest:: vertical_logLinearInterpolation
            >>> A=logLinearInterpolation(A,P) # interpolate A using pressure field P over the default levels
    """

    try:
        nlev = len(levels)  # Number of pressure levels
    except BaseException:
        nlev = 1  # if only one level len(levels) would breaks
        levels = [levels, ]
    order = list(A.coords)
    lev_idx = [idim for idim,dim in enumerate(order) if 'lev' in dim]
    lev_nm = order[lev_idx[0]]

    # The input pressure field needs to be TOP to BOTTOM
    if P[0,-1,0,0] < P[0,0,0,0]:
        print('Reverse pressure field into TOP to BOTTOM as logLinearRegression required')
        P = P[:,::-1,:]
        A = A[:,::-1,:]
    else:
        print('Top level=',P[0,-1,0,0].values, 'Bottom level=',P[0,0,0,0].values)

    A = A.transpose(lev_nm,...)
    P = P.transpose(lev_nm,...)
    sh = list(P.shape)
    nsigma = sh[0]  # number of sigma levels
    sh[0] = nlev
    t = np.ma.zeros(sh, dtype=np.float32)
    sh2 = P[0].shape
    prev = -1

    t = find_level_loop(A.to_masked_array(),P.to_masked_array(),levels,nlev,nsigma,sh2,t)

    ax = A.coords
    lvl = levels 
    try:
        lvl.units = P.units
    except BaseException:
        pass

    try:
        t.units = P.units
    except BaseException:
        pass

    coords = {'lev':lvl,'time':A.coords['time'], 'lat':A.coords['lat'], 'lon':A.coords['lon'],}
    t = xr.DataArray(t, coords=coords).transpose('time',lev_nm,'lat','lon')

    return t

# ============================================================================================
def find_level_loop(A,P,levels,nlev,nsigma,sh2,t):

    for ilev in range(nlev):  # loop through pressure levels
        print('Interpolating to targeted level =',levels[ilev])
        lev = levels[ilev]  # get value for the level
        Pabv = np.ones(sh2, dtype=np.float32)
        Aabv = -1 * Pabv  # Array on sigma level Above
        Abel = -1 * Pabv  # Array on sigma level Below
        Pbel = -1 * Pabv  # Pressure on sigma level Below
        Pabv = -1 * Pabv  # Pressure on sigma level Above
        Peq = np.ma.masked_equal(Pabv,-1)  # Area where Pressure == levels
        for i in range(1, nsigma):  # loop from second sigma level to last one
            a = np.ma.greater_equal(
                P[i],
                lev)  # Where is the pressure greater than lev
            b = np.ma.less_equal(
                P[i - 1],
                lev)  # Where is the pressure less than lev
            # Now looks if the pressure level is in between the 2 sigma levels
            # If yes, sets Pabv, Pbel and Aabv, Abel
            a = np.ma.logical_and(a, b)
            Pabv = xr.where(a, P[i], Pabv)  # Pressure on sigma level Above
            Aabv = xr.where(a, A[i], Aabv)  # Array on sigma level Above
            Pbel = xr.where(a, P[i - 1], Pbel)  # Pressure on sigma level Below
            Abel = xr.where(a, A[i - 1], Abel)  # Array on sigma level Below
            Peq = np.ma.where(np.ma.equal(P[i], lev), A[i], Peq)
    
        val = np.ma.masked_where(
            np.ma.equal(Pbel, -1), np.ones(Pbel.shape) * lev) # set to missing value if no data below lev if there is
    
        tl = np.log(val / Pbel) / np.log(Pabv / Pbel) * (Aabv - Abel) + Abel  # Interpolation

        if ((Peq.mask is None) or (Peq.mask is np.ma.nomask)):
            tl = Peq
        else:
            tl = xr.where(1 - Peq.mask, Peq, tl)

        t[ilev] = tl.astype(np.float32)

    return t

