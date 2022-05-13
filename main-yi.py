
import os
import cal_global_radiation_feedback_E3SM as GRF
import cal_RadKernel_E3SM as RK
import cal_CloudRadKernel_E3SM as CRK
import cal_LCF_E3SM as LCF
import cal_cloud_E3SM as CLOUD
import cal_webb_decomposition as WD
import sys
import allplots as AP
import cases_lookup as CL

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
    RunDiag = False
else:
    RunDiag = True

#########RunDiag = True

cal_types = [
#'RadFeedback',
#'RadKernel',
#'Webb_Decomp',
#'CloudRadKernel',
#'cal_LCF',
#'cal_cloud',
#'cal_EIS',
##'sort_cloud_regime',
#'sort_cloud_3regime',
'RadKernel_regime',
]

# give one stamp for each pair experiment, like v1, v2, v3....
case_stamp = [\
#'v2_coupled',
#'v1',\
#'v2',\
#'v1.All',\
#'v2.bk.clubb',\
#'v2.bk.clubb.MG',\
#'v2.bk.clubb.MG.ZM',\
#'v2.bk.clubb.MG.ZM.gust',\
#'v2.bk.clubb.MG.ZM.gust.trig',\
#'v2.bk.clubb.MG.ZM.gust.trig.gw',\
'v2.OutTend',\
#'v2.bk.trig',\
'v2.bk.MG',\
#'v2.bk.ZM',\
#'v2.bk.so4',\
'v2.bk.MG_Berg',\
'v2.bk.MG_mincdnc',\
#'v1.SSTv1',\
#'v2.SSTv2',\
#'v1.SSTv1.150yr',\
#'v1.pert1e-14',\
#'v2.bk.clubb.trig',\
#'v2.bk.clubb.trig.MG',\
#'F2010-p4Ka.v2.Allv1nml',\
#'F2010-p4Ka.v2.Allv1nml.ipdf',\
]

# set start year
yearS1 = 2 #101
yearS2 = 2 #1
# set end year
yearE1 = 6 #250
yearE2 = 6 #150

# set output directory for necessary variables after post-processing E3SM raw data
outdir_out = '/compyfs/qiny108/diag_feedback_E3SM_postdata/'

### NOTION: if you work on Cori, you should not change the below directories. If not, you should download data.
# set RadKernel input kernel directory
RadKernel_dir = '/qfs/people/qiny108/diag_feedback_E3SM/Huang_kernel_data/'

# ---------------------------- all main modifications please stop here --------------------------------------------

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
    #run_id1 = run_id1s[icase]
    #run_id2 = run_id2s[icase]

    run_id1 = CL.get_lutable(case,'amip')
    run_id2 = CL.get_lutable(case,'amip4K')
    print(case)
    print(run_id1)
    print(run_id2)

    #################################################################
    # add option to select component -- for processing E3SMv2 data,
    # which uses 'eam.hx' rather than 'cam.hx'
    #################################################################
    if 'v2rc1c' in case or 'v2_coupled' in case or 'v2' in case or 'All' in case:
        # set input directory 1 --- the directory before casename in the whole directory
        datadir_in1= '/compyfs/qiny108/E3SMv2_simulations/'
        # set input directory 2 --- the directory after casename in the whole directory
        datadir_in2 = 'archive/atm/hist/'

        comp = 'eam.h0'
        rgr_map = '/qfs/people/zender/data/maps/map_ne30pg2_to_cmip6_180x360_aave.20200201.nc'

    else:
        # set input directory 1 --- the directory before casename in the whole directory
        datadir_in1= '/compyfs/qiny108/E3SM_simulations/'
        # set input directory 2 --- the directory after casename in the whole directory
        datadir_in2 = 'archive/atm/hist/'

        comp = 'cam.h0'
        rgr_map = "/qfs/people/zender/data/maps/map_ne30np4_to_cmip6_180x360_aave.20181001.nc"
        
    #################################################################
    # pre-process model output to get necessary input files
    #################################################################
    if PreProcess:
#        # Part -1: archive data from run directory to archive/atm/hist directory
#        datadir_in0 = 'case_scripts/'
#        os.system('sh archive_data.sh '+run_id1 + ' '+run_id2+' '+datadir_in1+' '+datadir_in0)
#
#        # Part 0: regrid model output if your data in on SE grid. Otherwise, this step can be skipped. 
#        os.system('sh regrid_data.sh '+ run_id1 + ' '+run_id2+' '+rgr_map+' '+str(yearS1)+' '+str(yearE1)+' '+str(yearS2)+' '+str(yearE2)+' '+datadir_in1+' '+datadir_in2+' '+outdir_out+' '+comp+' '+regrid_SE2FV)

        # Part 1: pre-process model output to generate necessary files for further diagnostics
        #os.system('sh get_data.sh '+ run_id1 + ' '+run_id2+' '+rgr_map+' '+str(yearS1)+' '+str(yearE1)+' '+str(yearS2)+' '+str(yearE2)+' '+datadir_in1+' '+datadir_in2+' '+outdir_out+' '+comp)
        #os.system('sh get_data_tendency.sh '+ run_id1 + ' '+run_id2+' '+rgr_map+' '+str(yearS1)+' '+str(yearE1)+' '+str(yearS2)+' '+str(yearE2)+' '+datadir_in1+' '+datadir_in2+' '+outdir_out+' '+comp)
        os.system('sh get_data_tendency_MG.sh '+ run_id1 + ' '+run_id2+' '+rgr_map+' '+str(yearS1)+' '+str(yearE1)+' '+str(yearS2)+' '+str(yearE2)+' '+datadir_in1+' '+datadir_in2+' '+outdir_out+' '+comp)



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
            
            
            
