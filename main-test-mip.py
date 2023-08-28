
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
from get_mip_data import read_pickle,write_pickle

# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# ----------------------------Hi, all modifications please start here --------------------------------------------
#machine = 'cori'
#machine = 'compy'
machine = 'feedback2'

# model output directory run_dir and webpage directory www_dir
if machine == 'compy':
    www_dir = "/compyfs/www/qiny108/"
elif machine == 'cori':
    www_dir = "/global/project/projectdirs/mp193/www/qinyi/"
elif machine == 'feedback2':
    www_dir = None

hasCOSP = True   # True: include COSP output 
RunDiag = True # True: run feedback calculation
GetFigure = True # True: run figure plotting and generate webpage

# ------------------------------
# Input requirements for MIP data
# ------------------------------
model = 'GFDL-CM4'  
institution = 'NOAA-GFDL'
variant = 'r1i1p1f1'
path = '/p/css03/esgf_publish/CMIP6'
phase = 'CMIP6'
exps = ['amip','amip-p4K']

# All variables needed in this package 
Vars = [
"rlut","rsdt","rsut","rlutcs","rsutcs","ts","tas",
"rsus","rsds","rsdscs","rsuscs","ps","psl",
"ta","hus",
"clisccp",
]

tslice = slice("1980-01-01","1981-12-31")

# ------------ get all needed variables' directory and save them for later use ----------------------------------
ff = 'filenames_'+model+'_'+variant+'.pickle'

if os.path.isfile(ff):
    print(f'filenames of {model} is ready. Great!')
else:
    get_var_dirs(ff,exps,Vars,model,variant,phase)

filenames = read_pickle(ff)
#print('filenames=',filenames)

# give one shorname for each pair experiment, like v1, v2, v3....
case_short = [\
model+'_'+variant,\
]

# give the reference case: the reference case is used to compare with the case and generate difference maps. 
# Note: 
# 1. The length of "ref_case_short" must be the same as "case_short".
# 2. Any cases in ref_case_short should be in "case_short" first.
# 3. For each case, the corresponding reference cases can be any lengths.
# 4. If you don't have any reference cases for one case, just set it as [].
ref_case_short = [\
[],\
]

if machine in ['cori','compy']:
    Add_otherCMIPs = True ## include option about whether adding results from other CMIP models
else:
    Add_otherCMIPs = False

### NOTION: if you work on Cori, you should not change the below directories. If not, you should download data.
# set RadKernel input kernel directory
if machine == 'compy':
    RadKernel_dir = '/qfs/people/qiny108/diag_feedback_E3SM/Huang_kernel_data/'
elif machine == 'cori':
    RadKernel_dir = '/global/project/projectdirs/mp193/www/qinyi/DATA/Huang_kernel_data/'
elif machine == 'feedback2':
    RadKernel_dir = '/p/user_pub/climate_work/qin4/From_Compy/home_dir/diag_feedback_E3SM/Huang_kernel_data/'


# RunDiag types 
if RunDiag:
    cal_types = [
    'RadFeedback',
    'RadKernel',
    'Webb_Decomp',
    'CloudRadKernel',
#    'cal_LCF',
    ]

if GetFigure:
    # ---------------- please set all plot types you want -----------------------------------------------------------------
    ## choose which type figures you want to plot. If not, just comment them out.
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
#    'LCF',                              # Temperature - Liquid Condensate Fraction
    'webb_decomp',                      # decomposition of adjusted CRE feedback into low and non-low clouds
    ]

# ---------------------------- all main modifications please stop here --------------------------------------------
# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

# current directory
curdir = os. getcwd()

# set CloudRadKernel input kernel directory
CloudRadKernel_dir = curdir+'/CloudRadKernel_input/'

# set final output directory
outdir_final = curdir+'/data/'
# set output figure directory
diagfigdir = curdir+'/diagfigure/'

# ---------Create Directories--------------------------------------------
for outdir in [outdir_final, diagfigdir]:
    AP.make_dir(outdir)

# ----------------------------------------------------------------------------
for icase,case in enumerate(case_short):

    #################################################################
    # run diagnostics
    #################################################################
    if RunDiag:
    
        dics_cal = AP.get_cal_dics_MIP(case_short[icase], tslice, filenames, outdir_final,
                          RadKernel_dir, diagfigdir, 
                          CloudRadKernel_dir)

        for key in dics_cal:
            if key in cal_types:
                if not hasCOSP and 'CloudRadKernel' in key: 
                    continue

                dics_cal[key]()
 
#################################################################
# Generate plots and get webpage 
#################################################################
    
if GetFigure:
    # ----------- set up directories for necessary data --------------------------
    if machine == 'compy':
        datadir_CMIPs = '/compyfs/qiny108/diag_feedback_otherCMIPs/'
    elif machine == 'cori':
        datadir_CMIPs = '/global/project/projectdirs/mp193/www/qinyi/DATA/'
    else:
        datadir_CMIPs = None
    
    if datadir_CMIPs != None:
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
    
    ## data directory for all diagnostic outputs
    ## [it includes all data that you want to be plotted. If main.py runs successfully, this directory would be enough for further plot.]
    datadir_in = datadir+'/data/'
    
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
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    # get dictionary of all plot type lists
    dics_plots = AP.get_plot_dics(case_short,ref_case_short,figdir,datadir_in,Add_otherCMIPs)
    
    for key in dics_plots:
        if key in plot_types:
            # skip CldRadKernel if no COSP output
            if not hasCOSP and 'CloudRadKernel' in key: 
                continue

            pd2html = dics_plots[key]()
            # save pandas dataframe to csv file
            pd2html.to_csv(csvdir+"pd2html_"+key+".csv")
    
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # generate html file
    if machine in ['compy','cori']:
        gh.generate_html(casedir,webdir)
    else:
        gh.generate_html(casedir,casedir)
   
    print('Well Done.')



