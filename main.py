
import os
import cal_global_radiation_feedback_E3SM as GRF
import cal_RadKernel_E3SM as RK
import cal_CloudRadKernel_E3SM as CRK
import cal_LCF_E3SM as LCF
import cal_cloud_E3SM as CLOUD
import cal_webb_decomposition as WD
import sys
import allplots as AP

# ----------------------------
# all required environments for this script:
# python (anaconda)
# latest CDAT
# cartopy
# NCO

# ----------------------------Hi, all modifications please start here --------------------------------------------
# if you run this script at the first time, set PreProcess = True; otherwise, PreProcess = False
PreProcess = True
regrid_SE2FV = "True" # True: need regrid from SE2FV; otherwise, no need regrid.

if PreProcess:
    RunDiag = True
else:
    RunDiag = False

cal_types = [
'RadFeedback',
'RadKernel',
'Webb_Decomp',
'CloudRadKernel',
'cal_LCF',
'cal_cloud',
]

# give one stamp for each pair experiment, like v1, v2, v3....
case_stamp = [\
'F2010-p4K-all',\
]

# set the control case names
run_id1s = [\
'20210101.F2010C5-CMIP6-LR.ne30_oECv3.syrah.1024',\
]

# set the p4K case names
run_id2s = [\
'20210105.F2010C5-CMIP6-LR-p4K-all.ne30_oECv3.syrah.1024',\
]

# set start year
yearS1 = 1 #101
yearS2 = 1
# set end year
yearE1 = 5 #250
yearE2 = 5 #150

# set input directory 1 --- the directory before casename in the whole directory
datadir_in1= '/p/lustre2/qin4/E3SM_simulations/'
#datadir_in1 = '/p/lustre2/qin4/Data_cori/'
# set input directory 2 --- the directory after casename in the whole directory
datadir_in2 = 'run/'
#datadir_in2 = 'archive/atm/hist/'

# set output directory for necessary variables after post-processing E3SM raw data
outdir_out = '/p/lustre2/qin4/diag_feedback_E3SM_postdata/'

### NOTION: if you work on Cori, you should not change the below directories. If not, you should download data.
# set RadKernel input kernel directory
RadKernel_dir = '/p/lustre2/qin4/Data_cori/Huang_kernel_data/'

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

# ---------Create Directories--------------------------------------------
for outdir in [outdir_out, outdir_final, figdir]:
    AP.make_dir(outdir)

# ----------------------------------------------------------------------------
for icase,case in enumerate(case_stamp):
    run_id1 = run_id1s[icase]
    run_id2 = run_id2s[icase]

    #################################################################
    # add option to select component -- for processing E3SMv2 data,
    # which uses 'eam.hx' rather than 'cam.hx'
    #################################################################
    if 'v2rc1c' in case or 'v2_coupled' in case:
        comp = 'eam.h0'
        #rgr_map = '/p/lustre2/qin4/Data_cori/map_ne30pg2_to_cmip6_180x360_aave.20200201.nc'
        rgr_map = '/p/lustre2/qin4/Data_cori/map_ne30pg2_to_cmip6_72x144_aave.20210812.nc'
    else:
        comp = 'cam.h0'
        rgr_map = '/p/lustre2/qin4/Data_cori/map_ne30np4_to_fv129x256_aave.20150901.nc'
        
    #################################################################
    # pre-process model output to get necessary input files
    #################################################################
    if PreProcess:
        # Part -1: archive data from run directory to archive/atm/hist directory
        datadir_in0 = 'case_scripts/'
        os.system('sh archive_data.sh '+run_id1 + ' '+run_id2+' '+datadir_in1+' '+datadir_in0)

        # Part 0: regrid model output if your data in on SE grid. Otherwise, this step can be skipped. 
        os.system('sh regrid_data.sh '+ run_id1 + ' '+run_id2+' '+rgr_map+' '+str(yearS1)+' '+str(yearE1)+' '+str(yearS2)+' '+str(yearE2)+' '+datadir_in1+' '+datadir_in2+' '+outdir_out+' '+comp+' '+regrid_SE2FV)
        # Part 1: pre-process model output to generate necessary files for further diagnostics
        os.system('sh get_data.sh '+ run_id1 + ' '+run_id2+' '+rgr_map+' '+str(yearS1)+' '+str(yearE1)+' '+str(yearS2)+' '+str(yearE2)+' '+datadir_in1+' '+datadir_in2+' '+outdir_out+' '+comp)
        exit()

    #################################################################
    # run diagnostics
    #################################################################
    if RunDiag:
        direc_data = outdir_out

        dics_cal = AP.get_cal_dics(direc_data, case_stamp[icase], yearS2, yearE2, run_id1, run_id2, outdir_final,
                          RadKernel_dir, figdir, exp1, exp2,
                          CloudRadKernel_dir)

        for key in dics_cal:
            if key in cal_types:
                dics_cal[key]()
            
            
            
