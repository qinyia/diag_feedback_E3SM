#****************************************************************
#
#    Filename: cal_webb_decomposition.py
#
#    Author: Yi Qin - qin4@llnl.gov

#    Description: calculate the decomposition of adjusted CRE feedback by Webb et al. (2006) method.
#                 source codes from Mark Zelinka (zelinka1@llnl.gov)
#    Input: SW, LW and net adjusted CRE feedbacks [lat, lon]
#    Output: decomposed adjusted CRE feedbacks

#    Create: 2021-08-05 16:07:36
#    Last Modified: 2021-08-05 16:07:36
#****************************************************************

#IMPORT STUFF:
#=====================
import numpy as np
import xarray as xr
import sys
import os
from global_land_mask import globe 
sys.path.append('../')
import cases_lookup as CL
from PlotDefinedFunction import save_big_dataset
 
###########################################################################
# HELPFUL FUNCTIONS FOLLOW
###########################################################################

# ========================================================================
def cal_webb_decomp(direc_data,case_stamp,outdir):

    outfile = outdir+'RadKern_map_'+case_stamp+'-webb-decomp.nc'

    if os.path.isfile(outfile):
        print('cal_webb_decomp ', case_stamp, 'output is ready. Please continue. ')
        return
    else:
        lats_spc = np.arange(-90,92.5,2.5)
        lons_spc = np.arange(1.25,360,2.5)

        var = ['SWCRE_adj', 'LWCRE_adj', 'netCRE_adj']
        var_out = ['SWCRE','LWCRE','netCRE']

        # =============================================
        # read variables
        # =============================================
        dic_all = {}
        for ivar,svar in enumerate(var):
            if svar in var:
                print(svar)
                svar_out = var_out[ivar]

                f1 = xr.open_dataset(direc_data+'/RadKern_map_'+case_stamp+'.nc')
                pi_raw = f1[svar]
                f1.close()

                dic_all[svar_out] = pi_raw

            else:
                print('we dont find this variable:',svar,' in your file. Please check it!')
 
        print(dic_all.keys())

        # =============================================
        # do Webb Decomposition
        # =============================================
        lo_mask, nonlo_mask = webb_decomposition(dic_all['SWCRE'], dic_all['LWCRE'])
        print('lo_mask.shape=',lo_mask.shape)
        print('nonlo_mask.shape=',nonlo_mask.shape)

        dic_all_mask = {}

        dic_all_mask['SWCRE_lo']     = dic_all['SWCRE'].where(lo_mask == True)
        dic_all_mask['SWCRE_lo'].attrs['long_name'] = 'SW low cloud feedback'

        dic_all_mask['SWCRE_nonlo']  = dic_all['SWCRE'].where(nonlo_mask == True)
        dic_all_mask['SWCRE_nonlo'].attrs['long_name'] = 'SW non-low cloud feedback'

        dic_all_mask['LWCRE_lo']     = dic_all['LWCRE'].where(lo_mask == True)
        dic_all_mask['LWCRE_lo'].attrs['long_name'] = 'LW low cloud feedback'

        dic_all_mask['LWCRE_nonlo']  = dic_all['LWCRE'].where(nonlo_mask == True)
        dic_all_mask['LWCRE_nonlo'].attrs['long_name'] = 'LW non-low cloud feedback'

        dic_all_mask['netCRE_lo']    = dic_all['netCRE'].where(lo_mask == True)
        dic_all_mask['netCRE_lo'].attrs['long_name'] = 'NET low cloud feedback'

        dic_all_mask['netCRE_nonlo'] = dic_all['netCRE'].where(nonlo_mask == True)
        dic_all_mask['netCRE_nonlo'].attrs['long_name'] = 'NET non-low cloud feedback'


        # set units
        for svar in dic_all_mask.keys():
            dic_all_mask[svar].attrs['units'] = 'W/m2/K'

        # =============================================
        # save data into file     
        # =============================================
        comment = 'created by cal_webb_decomposition.py; Author: Yi Qin (yi.qin@pnnl.gov)'

        save_big_dataset(dic_all_mask,outfile,comment)

 # =========================================================================================
def webb_decomposition(sw,lw):
    # function to perform decomposition of Webb et al 2006. This is from Mark:)
    import numpy as np
    
    # input: maps of sw and lw cloud-induced radiation anomalies [lat,lon]
    
    tan225 = np.tan(22.5*np.pi/180.)
    SWpos = np.ma.greater(sw , 0)
    SWneg = np.ma.less(sw , 0)
    LWpos = np.ma.greater(lw , 0)
    LWneg = np.ma.less(lw , 0)
    
    # Classes A(S+LN) and E(S-LN)
    AE = np.ma.less_equal(np.abs(lw), tan225*np.abs(sw)) 
    # 1 where np.abs(lw)>np.tan(22.5)*np.abs(sw), 0 elsewhere
    A = np.ma.logical_and(AE, SWpos)
    E = np.ma.logical_and(AE, SWneg)
    
    # Classes D(S-L-) and H(S+L+)
    #Here the values of KSC and KLC are of the same sign,
    #abs(lw) >= tan(22.5)*abs(sw) and abs(sw) >= tan(22.5)*abs(lw)
    DH1 = np.ma.greater_equal(np.abs(lw), tan225*np.abs(sw)) 
    DH2 = np.ma.greater_equal(np.abs(sw), tan225*np.abs(lw)) 
    DH = np.ma.logical_and(DH1, DH2)
    D0 = np.ma.logical_and(DH, SWneg)
    D = np.ma.logical_and(D0, LWneg)
    H0 = np.ma.logical_and(DH, SWpos)
    H = np.ma.logical_and(H0, LWpos)
    
    #Classes B(S+L-) and F(S-L+) comprise the two sectors on the line KSC=-KLC and contain
    # values where KLC and KSC are of comparable magnitude but opposite sign.
    F0 = np.ma.logical_and(DH, SWneg)
    F = np.ma.logical_and(F0, LWpos)
    B0 = np.ma.logical_and(DH, SWpos)
    B = np.ma.logical_and(B0, LWneg)
    
    #Classes C(SNL-) and G(SNL+) contain values of sw 
    #which are relatively neutral compared with lw.
    CG = np.ma.less_equal(np.abs(sw), tan225*np.abs(lw)) 
    C = np.ma.logical_and(CG, LWneg)
    G = np.ma.logical_and(CG, LWpos)
    
    summ=A+B+C+D+E+F+G+H
    intsum = A.astype(int)+B.astype(int)+C.astype(int)+D.astype(int)+\
    		 E.astype(int)+F.astype(int)+G.astype(int)+H.astype(int)
    
    duplicates = np.ma.count(np.ma.where(intsum>1))
    if duplicates!=0:
        print('Some locations assigned to more than one category')
    
    noassign = np.ma.count(np.ma.where(intsum==0))
    if noassign!=0:
        print('Some locations not assigned to any category')
    
    # Soden and Vecchi (2011) aggregation:
    lo_mask = A+E
    hi_mask = B+C+F+G
    mix_mask = D+H
    nonlo_mask = B+C+D+F+G+H
    
    #return(A,B,C,D,E,F,G,H) 
    return lo_mask, nonlo_mask

# =========================================================================================

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if __name__ == "__main__":

    #direc_data = '/compyfs/qiny108/diag_feedback_E3SM_postdata/'
    #direc_data = '/qfs/people/qiny108/colla/diag_feedback_E3SM/data/'
    direc_data = '/p/user_pub/climate_work/qin4/From_Compy/home_dir/colla/diag_feedback_E3SM/codes/'

    case_stamps = [\
    'v2test'
    ]

    for case_stamp in case_stamps:

        fname1,_,_ = CL.get_lutable(case_stamp,'amip')
        fname2,_,_ = CL.get_lutable(case_stamp,'amip4K')

        outdir = './'

        exp1 = 'FC5'
        exp2 = 'FC5_4K'

        cal_webb_decomp(direc_data,case_stamp,outdir)

    model = 'GFDL-CM4'  
    variant = 'r1i1p1f1'
    case_stamp  = model+'_'+variant 
    cal_webb_decomp(direc_data,case_stamp,outdir)
