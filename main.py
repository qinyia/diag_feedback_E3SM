
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
# if you run this script at the first time, set PreProcess = True; otherwise, PreProcess = False
PreProcess = False
Global_RadFeedback = True
RadKernel = True
CloudRadKernel = True

# give one stamp for each pair experiment, like v1, v2, v3....
case_stamp = [\
'TEST',\
#'amip-p4K',\
#'amip-4xCO2',\
#'amip-future4K'\
]
# set the control case names
run_id1s = [\
'cori-haswell.20190513.F2010C5-CMIP6-LR.ne30_oECv3',\
#'20200428.DECKv1b_amip1-CFMIP.ne30_oEC.cori-knl-L',\
#'20200428.DECKv1b_amip1-CFMIP.ne30_oEC.cori-knl-L',\
#'20200428.DECKv1b_amip1-CFMIP.ne30_oEC.cori-knl-L',\
]
# set the p4K case names
run_id2s = [\
'cori-haswell.20190513.F2010C5-CMIP6-LR.plus4K.ne30_oECv3',\
#'20200428.DECKv1b_amip1.plus4K-CFMIP.ne30_oEC.cori-knl-L',\
#'20201010.DECKv1b_amip1.4xCO2-CFMIP.ne30_oEC.cori-knl-L',\
#'20201011.DECKv1b_amip1.future4K-CFMIP.ne30_oEC.cori-knl-L',\
]
# set the regrid map from ne30np4 to lat-lon
rgr_map = '~zender/data/maps/map_ne30np4_to_cmip6_72x144_aave.20181001.nc'
# set start year
yearS = 1
# set end year
yearE = 2

# set input directory 1 --- the directory before casename in the whole directory
datadir_in1 = '/global/cscratch1/sd/qinyi/E3SM_predata/'
#datadir_in1= '/global/cscratch1/sd/qinyi/E3SM_simulations/'
# set input directory 2 --- the directory after casename in the whole directory
datadir_in2 = 'archive/atm/hist/'

# set output directory for necessary variables
outdir_out = '/global/cscratch1/sd/qinyi/diag_feedback_E3SM_postdata/'


### NOTION: if you work on Cori, you should not change the below directories.If not, you should download data.
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
        result = GRF.Global_RadFeedback(CloudRadKernel_dir,direc_data,case_stamp[icase],yearS,yearE,run_id1,run_id2,outdir_final)
    
    # Part 3: run radiative kernel analysis
    if RadKernel:
        result = RK.RadKernel(RadKernel_dir,direc_data,case_stamp[icase],yearS,yearE,run_id1,run_id2,outdir_final,figdir,exp1,exp2)
    
    # Part 4: run cloud radiative kernel analysis: decomposition 
    if CloudRadKernel:
        result = CRK.CloudRadKernel(CloudRadKernel_dir,direc_data,case_stamp[icase],yearS,yearE,run_id1,run_id2,outdir_final,figdir)



