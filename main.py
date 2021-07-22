
import os
import cal_global_radiation_feedback_E3SM as GRF
import cal_RadKernel_E3SM as RK
import cal_CloudRadKernel_E3SM as CRK

# ----------------------------
# all required environments for this script:
# python (anaconda)
# latest CDAT
# cartopy
# NCO

# ----------------------------Hi, all modifications please start here --------------------------------------------

## convert raw model output to desired data format for later analysis
## if you run this script at the first time, set PreProcess = True; otherwise, PreProcess = False
PreProcess = True
# calculate pure radiative feedbacks (Cess feedbacks)
Global_RadFeedback = True
# calculate radiative kernel feedbacks (Soden et al., 2005), including non-cloud feedbacks like Planck, lapse rate, water vapor, albedo, adjusted cloud feedbacks.
RadKernel = True
# calculate cloud radiative kernel feedbacks (Zelinka et al., 2012, 2016), including the contribution of cloud properties (amount, altitude, altitude) on cloud feedbacks.
CloudRadKernel = True

# give one stamp for each pair experiment, like v1, v2, v3..... Any names you like.
case_stamp = [\
'P4K-syrah',\
#'amip-p4K',\
]
# set the control case names
run_id1s = [\
'20210101.F2010C5-CMIP6-LR.ne30_oECv3.syrah.1024',\
#'20200428.DECKv1b_amip1-CFMIP.ne30_oEC.cori-knl-L',\
]
# set the p4K case names
run_id2s = [\
'20210101.F2010C5-CMIP6-LR-p4K.ne30_oECv3.syrah.1024',\
#'20200428.DECKv1b_amip1.plus4K-CFMIP.ne30_oEC.cori-knl-L',\
]
# set the regrid map from ne30np4 to lat-lon
#rgr_map = '~zender/data/maps/map_ne30np4_to_cmip6_72x144_aave.20181001.nc'
rgr_map = '~zender/data/maps/map_ne30np4_to_fv129x256_aave.20150901.nc'
# set start year
yearS = 1
# set end year
yearE = 2

# set input directory 1 --- the directory before casename in the whole directory
#datadir_in1 = '/global/cscratch1/sd/qinyi/E3SM_predata/'
#datadir_in1= '/global/cscratch1/sd/qinyi/E3SM_simulations/'
datadir_in1 = '/global/cscratch1/sd/qinyi/From_syrah/'
# set input directory 2 --- the directory after casename in the whole directory
#datadir_in2 = 'archive/atm/hist/h0/'
#datadir_in2 = 'archive/atm/hist/'
datadir_in2 = 'atm/hist/'

# set output directory for necessary variables after post-processing E3SM raw data -- caution: this data might be large.
outdir_out = '/global/cscratch1/sd/qinyi/diag_feedback_E3SM_postdata/'


### NOTION: if you work on Cori, you should not change the below directories. If not, you should download data.
# set RadKernel input kernel directory
RadKernel_dir = '/global/project/projectdirs/mp193/www/qinyi/DATA/Huang_kernel_data/'

# ---------------------------- all modifications please stop here --------------------------------------------

# current directory
curdir = os. getcwd()

# set CloudRadKernel input kernel directory
CloudRadKernel_dir = curdir+'/CloudRadKernel_input/'

# set final output directory 
outdir_final = curdir+'/data/'
# set output figure directory
figdir = curdir+'/figure/'

# set the case component
exp1 = 'FC5'
exp2 = 'FC5_4K'

# --------------
try:
    os.mkdir(outdir_out)
except OSError:
    print("Creation of the directory %s failed" % outdir_out)
else:
    print("Successfully created the directory %s " % outdir_out)

try:
    os.mkdir(outdir_final)
except OSError:
    print("Creation of the directory %s failed" % outdir_final)
else:
    print("Successfully created the directory %s " % outdir_final)

try:
    os.mkdir(figdir)
except OSError:
    print("Creation of the directory %s failed" % figdir)
else:
    print("Successfully created the directory %s " % figdir)
 
# ----------------------------------------------------------------------------
for icase,case in enumerate(case_stamp):


    run_id1 = run_id1s[icase]
    run_id2 = run_id2s[icase]
    
    # Part 1: pre-process model output
    if PreProcess:
        os.system('sh get_data.sh '+ run_id1 + ' '+run_id2+' '+rgr_map+' '+str(yearS)+' '+str(yearE)+' '+datadir_in1+' '+datadir_in2+' '+outdir_out)

    # Part 2: run global radiative feedback analysis
    direc_data = outdir_out
    if Global_RadFeedback:
        result = GRF.Global_RadFeedback(direc_data,case_stamp[icase],yearS,yearE,run_id1,run_id2,outdir_final)
    
    # special treatment to amip-4xCO2: it does not need RadKernel and CloudRadKernel analysis
    if case != 'amip-4xCO2':
        # Part 3: run radiative kernel analysis
        if RadKernel:
            result = RK.RadKernel(RadKernel_dir,direc_data,case_stamp[icase],yearS,yearE,run_id1,run_id2,outdir_final,figdir,exp1,exp2)
        
        # Part 4: run cloud radiative kernel analysis: decomposition 
        if CloudRadKernel:
            result = CRK.CloudRadKernel(CloudRadKernel_dir,direc_data,case_stamp[icase],yearS,yearE,run_id1,run_id2,outdir_final,figdir)

