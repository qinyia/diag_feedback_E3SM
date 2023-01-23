
import os
import sys
sys.path.append('./codes/')
import cases_lookup as CL
import PlotDefinedFunction as PDF
import allplots as AP
import cal_global_radiation_feedback_E3SM as GRF
import cal_RadKernel_E3SM as RK
import cal_CloudRadKernel_E3SM as CRK
import cal_LCF_E3SM as LCF
import cal_webb_decomposition as WD
import generate_html as gh

# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# ----------------------------Hi, all modifications please start here --------------------------------------------
machine = 'cori'
#machine = 'compy'

# model output directory www_dir and webpage directory run_dir
if machine == 'compy':
    www_dir = "/compyfs/www/qiny108/"
    run_dir = "/compyfs/qiny108/" 
elif machine == 'cori':
    www_dir = "/global/project/projectdirs/mp193/www/qinyi/"
    #run_dir = "/global/cscratch1/sd/qinyi/"
    run_dir = "/global/project/projectdirs/mp193/www/qinyi/diagfbk_testdata/"

e3sm_version = 1 # E3SM version 

PreProcess = True  # True: prepare input data for feedback calculation, including regrid and reorganize data
COSP_output = True # True: you have COSP output; False: no COSP output

RunDiag = False # True: run feedback calculation

GetFigure = False # True: run figure plotting and generate webpage

if GetFigure:
    # -------------------------
    # If you are not on cori machine, please set it as False. The comparison with other CMIP models are not 
    # supported on other machines currently. 
    Add_otherCMIPs = False ## include option about whether adding results from other CMIP models


# give one shorname for each pair experiment, like v1, v2, v3....
case_short = [\
'v1',\
#'v2',\
]

# give the reference case: the reference case is used to compare with the case and generate difference maps. 
# Note: 
# 1. The length of "ref_case_short" must be the same as "case_short".
# 2. Any cases in ref_case_short should be in "case_short" first.
# 3. For each case, the corresponding reference cases can be any lengths.
# 4. If you don't have any reference cases for one case, just set it as [].
ref_case_short = [\
[],\
#['v1'],\
]

# set start and end years: sometime, your start and end years are different for control (CTL) and warming (P4K) exps. 
yearS_CTL,yearE_CTL = 2,3
yearS_P4K,yearE_P4K = 2,3 

# set model output data directory 
if e3sm_version == 2: # E3SM version 2
    # set input directory 1 --- the directory before casename in the whole directory
    datadir_in1= run_dir+'/E3SMv2_simulations/'
    # set input directory 2 --- the directory after casename in the whole directory
    datadir_in2 = 'archive/atm/hist/'

    comp = 'eam.h0'
    if machine == 'compy':
        rgr_map = '/qfs/people/zender/data/maps/map_ne30pg2_to_cmip6_180x360_aave.20200201.nc'
    elif machine == 'cori':
        rgr_map = '/global/cfs/cdirs/e3sm/zender/maps/map_ne30pg2_to_cmip6_180x360_aave.20200201.nc'

elif e3sm_version == 1: # E3SM version 1
    # set input directory 1 --- the directory before casename in the whole directory
    datadir_in1= run_dir+'/E3SM_simulations/'
    # set input directory 2 --- the directory after casename in the whole directory
    datadir_in2 = 'archive/atm/hist/'

    comp = 'cam.h0'
    if machine == 'compy':
        rgr_map = "/qfs/people/zender/data/maps/map_ne30np4_to_cmip6_180x360_aave.20181001.nc"
    elif machine == 'cori':
        rgr_map = "/global/cfs/cdirs/e3sm/zender/maps/map_ne30np4_to_cmip6_180x360_aave.20181001.nc"

 
# set output directory for necessary variables after post-processing E3SM raw data
outdir_out = run_dir+'/diag_feedback_E3SM_postdata/'

### NOTION: if you work on Cori, you should not change the below directories. If not, you should download data.
# set RadKernel input kernel directory
if machine == 'compy':
    RadKernel_dir = '/qfs/people/qiny108/diag_feedback_E3SM/Huang_kernel_data/'
elif machine == 'cori':
    RadKernel_dir = '/global/project/projectdirs/mp193/www/qinyi/DATA/Huang_kernel_data/'

# ---------------------------- all main modifications please stop here --------------------------------------------
# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

# RunDiag types 
if COSP_output: 
    cal_types = [
    'RadFeedback',
    'RadKernel',
    'Webb_Decomp',
    'CloudRadKernel',
    'cal_LCF',
    'cal_3dvar',
    #'RadKernel_regime',
    ]
else:
    cal_types = [
    'RadFeedback',
    'RadKernel',
    'Webb_Decomp',
    'cal_LCF',
    'cal_cloud',
    'cal_3dvar',
    #'RadKernel_regime',
    ]

# current directory
curdir = os. getcwd()

# set CloudRadKernel input kernel directory
CloudRadKernel_dir = curdir+'/CloudRadKernel_input/'

# set final output directory
outdir_final = curdir+'/data/'
# set output figure directory
figdir = curdir+'/figure/'

# set the case tag for control and warming experiments. Dont modify it.
exp1 = 'FC5'
exp2 = 'FC5_4K'

# ---------Create Directories--------------------------------------------
for outdir in [outdir_out, outdir_final, figdir]:
    AP.make_dir(outdir)

# ----------------------------------------------------------------------------
for icase,case in enumerate(case_short):

    run_id1 = CL.get_lutable(case,'amip')
    run_id2 = CL.get_lutable(case,'amip4K')
    print(case)
    print(run_id1)
    print(run_id2)

    #################################################################
    # pre-process model output to get necessary input files
    #################################################################
    if PreProcess:
        # Part -1: archive data from run directory to archive/atm/hist directory
        datadir_in0 = 'case_scripts/'
        os.system('sh archive_data.sh '+run_id1 + ' '+run_id2+' '+datadir_in1+' '+datadir_in0)
        # Part 1: regrid data  
        os.system('sh get_data_simple.sh '+ run_id1 + ' '+run_id2+' '+rgr_map+' '+str(yearS_CTL)+' '+str(yearE_CTL)+' '+str(yearS_P4K)+' '+str(yearE_P4K)+' '+datadir_in1+' '+datadir_in2+' '+outdir_out+' '+comp)

    #continue

    #################################################################
    # run diagnostics
    #################################################################
    if RunDiag:
        direc_data = outdir_out
    
        dics_cal = AP.get_cal_dics(direc_data, case_short[icase], yearS_P4K, yearE_P4K, run_id1, run_id2, outdir_final,
                          RadKernel_dir, figdir, exp1, exp2,
                          CloudRadKernel_dir)

        for key in dics_cal:
            if key in cal_types:
                dics_cal[key]()
 
#################################################################
# Generate plots and get webpage 
#################################################################
if GetFigure:

    # ---------------- please set all plot types you want -----------------------------------------------------------------
    ## choose which type figures you want to plot. If not, just comment them out.
    if COSP_output:
        plot_types = [
        'CRE_globalmean',                   # scatter plot of global mean CRE feedbacks
        'RadKernel_globalmean',             # scatter plot of global mean RadKernel feedback: non-cloud and adjusted CRE feedbacks
        'RadKernel_zonalmean',              # zonal mean plot of adjusted CRE feedback
        'CldRadKernel_globalmean',          # scatter plot of global mean CldRadKernel feedback: decomposition into low and non-low clouds and amount, altitude, optical depth.
        'CldRadKernel_zonalmean',           # zonal mean plot of CldRadKernel feedback
        'RadKernel_latlon',                 # lat-lon plot of RadKernel feedback for each case
        'CldRadKernel_latlon',              # lat-lon plot of CldRadKernel feedback for each case
        'CldRadKernel_latlon_dif',          # lat-lon plot of CldRadKernel feedback difference between case and reference case
        'RadKernel_latlon_dif',             # lat-lon plot of RadKernel feedback difference between case and reference case
        'tas_latlon',                       # lat-lon plot of surface air temperature and the difference between case and reference case
        'LCF',                              # Temperature - Liquid Condensate Fraction
        'zm_CLOUD',                         # zonal mean plot of cloud varaibles difference 
        'latlon_CLOUD',                     # lat-lon plot of cloud varaibles difference
        'webb_decomp',                      # decomposition of adjusted CRE feedback into low and non-low clouds
        ]
    else:
        plot_types = [
        'CRE_globalmean',                   # scatter plot of global mean CRE feedbacks
        'RadKernel_globalmean',             # scatter plot of global mean RadKernel feedback: non-cloud and adjusted CRE feedbacks
        'RadKernel_zonalmean',              # zonal mean plot of adjusted CRE feedback
        #'CldRadKernel_globalmean',          # scatter plot of global mean CldRadKernel feedback: decomposition into low and non-low clouds and amount, altitude, optical depth.
        #'CldRadKernel_zonalmean',           # zonal mean plot of CldRadKernel feedback
        'RadKernel_latlon',                 # lat-lon plot of RadKernel feedback for each case
        #'CldRadKernel_latlon',              # lat-lon plot of CldRadKernel feedback for each case
        #'CldRadKernel_latlon_dif',          # lat-lon plot of CldRadKernel feedback difference between case and reference case
        'RadKernel_latlon_dif',             # lat-lon plot of RadKernel feedback difference between case and reference case
        'tas_latlon',                       # lat-lon plot of surface air temperature and the difference between case and reference case
        'LCF',                              # Temperature - Liquid Condensate Fraction
        'zm_CLOUD',                         # zonal mean plot of cloud varaibles difference 
        'latlon_CLOUD',                     # lat-lon plot of cloud varaibles difference
        'webb_decomp',                      # decomposition of adjusted CRE feedback into low and non-low clouds
        ]

    
    # ---------------- please set other optional setting for figure: start -------------------------------------------------
    colors = PDF.get_color('tab10',len(case_short)) #['tab:red','tab:blue','tab:cyan','tab:orange','tab:purple','tab:green']
        
    linewidths = [2]
    linestyles = ['--']
    linewidths.extend([3]*(len(case_short)-1))
    linestyles.extend(['-']*(len(case_short)-1))
    
    fh = 15     # font size
    fh1 = 13    # font size for legend
    s1 = 120    # marker size for E3SMv2
    s2 = 100    # marker size for other CMIP models
    a1 = 1      # apparency for markers
    
    # advanced setting: if plot_CldRadKernel_zonalmean, control the number of cases you want to show in one figure.
    # for example, if you would like to show the first three cases, then first 6 cases, and all cases, pls set ncase = [3,6,7]
    # generally, if you want all lines in the same plot, just set ncase = [len(cases)]
    ncase = [len(case_short)]
    #if len(case_short) == 1:
    #    ncase = [1]
    #else:
    #    ncase = range(2,len(case_short)+1)
    
    #print('ncase=',list(ncase))

    # add region ranges: [latS, latE, lonS, lonE] to help generate regional figures 
    regions = [-90,90,0,360] 
    #regions = [10,70,220,310]
    #regions = [24,55.5,-140,-70]

    # ----------- set up directories for necessary data --------------------------
    if machine == 'compy':
        datadir_CMIPs = '/compyfs/qiny108/diag_feedback_otherCMIPs/'
    elif machine == 'cori':
        datadir_CMIPs = '/global/project/projectdirs/mp193/www/qinyi/DATA/'
    
    # -- data for E3SMv1 [dont modify data in this directory.]
    datadir_v1 = datadir_CMIPs+'E3SMv1_data/'
    # -- data for other CMIP models from CRE feedback [dont' modify data in it.]
    datadir_Ringer = datadir_CMIPs+'RadFeedback/'
    # -- data for other CMIP models for RadKernel feedback [don't modify data in it.]
    datadir_RadKernel = datadir_CMIPs+'RadKernel/'
    # -- data for other CMIP models for CldRadKernel feedabck [ don't modify data in it.]
    datadir_CldRadKernel = datadir_CMIPs+'CldRadKernel/'
    
    # ----------- please set all these following directories and your prefered styles: start ----------------------
    ## main directory. pls modify it based on your current script directory.
    datadir = os.getcwd()
    
    ## data directory for E3SMv2
    ## [it includes all data that you want to be plotted. If main.py runs successfully, this directory would be enough for further plot.]
    datadir_v2 = datadir+'/data/'
    
    casedir = datadir+'/'+case_short[-1]+'/'
    AP.make_dir(casedir)
    
    ## figure directory
    figdir = datadir+'/'+case_short[-1]+'/figure/'
    ## csv directory
    csvdir = datadir+'/'+case_short[-1]+'/csvfile/'
    ## viewer directory
    viewdir = datadir+'/'+case_short[-1]+'/viewer/'
    
    ## web file directory, like on compy or nersc
    if machine == 'compy':
        webdir = www_dir+"/diag_feedback/"+case_short[-1]+"/"
    elif machine == 'cori':
        webdir = www_dir+"/diag_feedback/"+case_short[-1]+"/"
    
    ## create figure directory if it does not exist
    AP.make_dir(figdir)
    AP.make_dir(csvdir)
    AP.make_dir(viewdir)
    
    if machine in ['compy','cori']:
        AP.make_dir(webdir)
        os.system("cp -p "+datadir+"/viewer/11.css "+viewdir)
    
    # the following lines might be removed later....
    Add_amipFuture = False
    highlight_CESM2 = False
    lw_CESM2 = 2
    ls_CESM2 = ':'
    lc_CESM2 = 'blue'
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    # get dictionary of all plot type lists
    dics_plots = AP.get_plot_dics(case_short,ref_case_short,Add_otherCMIPs,datadir_v2, datadir_v1, s1, s2, fh, fh1, a1, colors, figdir, ncase, linestyles, linewidths,Add_amipFuture,highlight_CESM2,lw_CESM2, ls_CESM2, lc_CESM2, datadir_Ringer, datadir_RadKernel, datadir_CldRadKernel,regions)
    
    for key in dics_plots:
        if key in plot_types:
            pd2html = dics_plots[key]()
            # save pandas dataframe to csv file
            pd2html.to_csv(csvdir+"pd2html_"+key+".csv")
    
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # generate html file
    if machine in ['compy','cori']:
        gh.generate_html(casedir,webdir)
    
    print('Well Done.')

