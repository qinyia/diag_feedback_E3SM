
#****************************************************************
#
#    Filename: main-plot.py
#
#    Author: Yi Qin - qin4@llnl.gov
#    Description: This script is used to plot all kinds of figures for E3SMv2 developing versions via comparing with
#                 E3SMv1-amip, E3SMv1-piControl, and other available CMIP5 and CMIP6 models. 
#    Input: 
#    Output: 
#    Create: 2020-08-04 
#    Last Modified: 2021-02-18 16:14:27

#    Dec 29, 2020: add plot for radiative forcing, especially for amip-4xCO2 experiment
#    Jul 18, 2021: add calculation of LCF
#    Jul 27, 2021: add plot_latlon_CLOUD; change plot_CLOUD to plot_zm_CLOUD
#    Aug 05, 2021: add plot_webb_decomp
#    Feb 01, 2022: add plot_CLOUD_profile

#****************************************************************
import sys
import os
import PlotDefinedFunction as PDF
import allplots as AP
import generate_html as gh

########################################## Modification starts here #################################
machine = 'compy'

## notion: if you also want to compare with default E3SM-1-0, please add 'v1_coupled' and 'v1_amip4K' below.

cases = ['F2010-p4K-all', 'F2010-p4K-all.IC']
ref_casesA = [[], ['F2010-p4K-all']]

# -------------------------
# If you are not on cori or LC machines, please set it as False. The comparison with other CMIP models are not 
# supported on other machines currently. 
Add_otherCMIPs = False ## include option about whether adding results from other CMIP models

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
'LCF',                              # Temperature - Liquid Condensate Fraction
'zm_CLOUD',                         # zonal mean plot of cloud varaibles difference 
'latlon_CLOUD',                     # lat-lon plot of cloud varaibles difference
'webb_decomp',                      # decomposition of adjusted CRE feedback into low and non-low clouds
'CLOUD_profile',                    # vertical cloud profile in different regions
#'NRMSE_RadKern',                   # NRMSE and spatial correlation (COR) evolution with incremental denial experiments [note: RadKernel_latlon_dif should run first.]
]

# ---------------- please set other optional setting for figure: start -------------------------------------------------
colors = PDF.get_color('tab10',len(cases)) #['tab:red','tab:blue','tab:cyan','tab:orange','tab:purple','tab:green']
    
linewidths = [2]
linestyles = ['--']
linewidths.extend([3]*(len(cases)-1))
linestyles.extend(['-']*(len(cases)-1))

fh = 15     # font size
fh1 = 13    # font size for legend
s1 = 120    # marker size for E3SMv2
s2 = 100    # marker size for other CMIP models
a1 = 1      # apparency for markers

# advanced setting: if plot_CldRadKernel_zonalmean, control the number of cases you want to show in one figure.
# for example, if you would like to show the first three cases, then first 6 cases, and all cases, pls set ncase = [3,6,7]
# generally, if you want all lines in the same plot, just set ncase = [len(cases)]
ncase = [len(cases)]
#if len(cases) == 1:
#    ncase = [1]
#else:
#    ncase = range(2,len(cases)+1)

print('ncase=',list(ncase))

#%%%%%%%%%%%%%%%%%% Stop modification here !!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# ----------- set up directories for necessary data --------------------------
if machine == 'LC':
    datadir_CMIPs = '/p/lustre2/qin4/Data_cori/'
elif machine == 'compy':
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

casedir = datadir+'/'+cases[-1]+'/'
AP.make_dir(casedir)

## figure directory
figdir = datadir+'/'+cases[-1]+'/figure/'
## csv directory
csvdir = datadir+'/'+cases[-1]+'/csvfile/'
## viewer directory
viewdir = datadir+'/'+cases[-1]+'/viewer/'

## web file directory, like on compy or nersc
if machine == 'compy':
    webdir = "/compyfs/www/qiny108/diag_feedback/"+cases[-1]+"/"
elif machine == 'cori':
    webdir = "/global/project/projectdirs/mp193/www/qinyi/diag_feedback/"+cases[-1]+"/"

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

if True:
    # get dictionary of all plot type lists
    dics_plots = AP.get_plot_dics(cases,ref_casesA,Add_otherCMIPs,datadir_v2, datadir_v1, s1, s2, fh, fh1, a1, colors, figdir, ncase, linestyles, linewidths,Add_amipFuture,highlight_CESM2,lw_CESM2, ls_CESM2, lc_CESM2, datadir_Ringer, datadir_RadKernel, datadir_CldRadKernel)

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

