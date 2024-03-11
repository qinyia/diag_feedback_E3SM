
import os
import sys
import cases_lookup_ppe as CLP
import PlotDefinedFunction as PDF
import allplots as AP
import cal_global_radiation_feedback_E3SM as GRF
import cal_RadKernel_E3SM as RK
import cal_CloudRadKernel_E3SM as CRK
import cal_LCF_E3SM as LCF
import cal_cloud_E3SM as CLOUD
import cal_webb_decomposition as WD
import generate_html as gh

#for yrE in [2010,2011,2012,2013,2014]: # loop over end years for calculating cloud feedback
for yrE in [2012]:
    # XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    # ----------------------------Hi, all modifications please start here --------------------------------------------
    machine = 'compy'
    
    # model output directory www_dir and webpage directory run_dir
    if machine == 'compy':
        www_dir = "/compyfs/www/qiny108/"
        run_dir = "/compyfs/qiny108/" 
    elif machine == 'cori':
        www_dir = "/global/project/projectdirs/mp193/www/qinyi/"
        run_dir = "/global/cscratch1/sd/qinyi/"
    
    e3sm_version = 2 # E3SM version 
    
    PreProcess = False  # True: prepare input data for feedback calculation, including regrid and reorganize data
    COSP_output = True # True: you have COSP output; False: no COSP output
    
    RunDiag = True # True: run feedback calculation
    
    GetFigure = False # True: run figure plotting and generate webpage
    
    # give one shorname for each pair experiment, like v1, v2, v3....
    case_short = [\
    #'BASE_FR',
    'BASE',
#    'nomincdnc',
#    'prc_exp1',
#    'berg',
#    'prc_exp',
#    'prc_coef1',
    'c1',
    'gamma_coef',
    'c8',
#    'accre_enhan',
#    'ice_deep',
#    'clubb_tk1',
#    'dp1',
#    'ice_sed_ai',
#    'so4_sz',
#    'prc_exp1_2',
#    'prc_exp_2',
#    'prc_coef1_2',
#    'ice_deep_2',
#    'dp1_2',
#    'ice_sed_ai_2',
#    'so4_sz_2',
#    'nomincdnc_2',
#    'nomincdnc_3',
#    'prc_exp1_3',
#    'prc_exp_3',
#    'prc_coef1_3',
    'c1_2',
    'c1_3',
    'gamma_coef_2',
    'gamma_coef_3', 
    'gamma_coefb',
    'c8_2',
    'c8_3',
#    'accre_enhan_2',
#    'accre_enhan_3',
#    'wsub',
#    'wsub_2',
#    'berg_2',
#    'berg_3',
#    'wsub_3',
#    'ice_deep_3',
#    'clubb_tk1_2',
#    'clubb_tk1_3',
#    'dp1_3',
#    'ice_sed_ai_3',
#    'so4_sz_3',
#    'BASE_01',
#    'BASE_02',
#    'BASE_03',
#    'nomincdnc.prc_exp1_2',
#    'nomincdnc.prc_exp1_3',
#    'nomincdnc.prc_exp_2',
#    'nomincdnc.prc_exp_3',
#    'nomincdnc.prc_coef1',
#    'nomincdnc.prc_coef1_3',
#    'nomincdnc.prc_exp1_v1',
#    'nomincdnc.berg',
#    'nomincdnc.ice_deep_2',
#    'nomincdnc.clubb_tk1_3',
#    'nomincdnc.dp1_3',
#    'nomincdnc.ice_sed_ai_3',
#    'nomincdnc.so4_sz_3',
#    'nomincdnc.accre_enhan_2',
#    'nomincdnc.c1_2',
#    'nomincdnc.c8_3',
#    'nomincdnc.gamma_coef',
#    'nomincdnc.wsub_3',
    'c1_4',
    'c8_4',
    'gamma_coef_4',
    'gamma_coefb_2',
    'gamma_coefb_3',
#    'wsub_4',
#    'BASE.vars',
    ]
    
    # give the reference case: the reference case is used to compare with the case and generate difference maps. 
    # Note: 
    # 1. The length of "ref_case_short" must be the same as "case_short".
    # 2. Any cases in ref_case_short should be in "case_short" first.
    # 3. For each case, the corresponding reference cases can be any lengths.
    # 4. If you don't need any reference cases for one case, just set it as [].
    ref_case_short = [
    [],
    ['BASE'],
    ['BASE'],
    ['BASE'],
    ['BASE'],
    ['BASE'],
    ['BASE'],
    ['BASE'],
    ]
    
    # set start and end years: sometime, your start and end years are different for control (CTL) and warming (P4K) exps. 
    #yearS_CTL,yearE_CTL = 2,6 #101,250
    #yearS_P4K,yearE_P4K = 2,6 #1,150
    
    # set model output data directory 
    if e3sm_version == 2: # E3SM version 2
        # set input directory 1 --- the directory before casename in the whole directory
        #datadir_in1= run_dir+'/E3SMv2_simulations/'
        datadir_in1= run_dir+'/E3SMv2_simulations/PPE/'
    
        # set input directory 2 --- the directory after casename in the whole directory
        #datadir_in2 = 'archive/atm/hist/'
        datadir_in2 = 'run/'
    
        comp = 'eam.h0'
        rgr_map = '/qfs/people/zender/data/maps/map_ne30pg2_to_cmip6_180x360_aave.20200201.nc'
    
    elif e3sm_version == 1: # E3SM version 1
        # set input directory 1 --- the directory before casename in the whole directory
        datadir_in1= run_dir+'/E3SM_simulations/'
        # set input directory 2 --- the directory after casename in the whole directory
        datadir_in2 = 'archive/atm/hist/'
    
        comp = 'cam.h0'
        rgr_map = "/qfs/people/zender/data/maps/map_ne30np4_to_cmip6_180x360_aave.20181001.nc"
     
    # set output directory for necessary variables after post-processing E3SM raw data
    outdir_out = run_dir+'/diag_feedback_E3SM_postdata_ppe/'
    
    ### NOTION: if you work on Cori, you should not change the below directories. If not, you should download data.
    # set RadKernel input kernel directory
    RadKernel_dir = '/qfs/people/qiny108/diag_feedback_E3SM/Huang_kernel_data/'
    
    # ---------------------------- all main modifications please stop here --------------------------------------------
    # XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    
    # RunDiag types 
    if COSP_output: 
        cal_types = [
#        'RadFeedback',
#        'RadKernel',
#        'Webb_Decomp',
#        'CloudRadKernel',
#        'cal_LCF',
#        'cal_APRP',
        #'cal_cloud',
        #'cal_EIS',
        #'cal_pr',
        'cal_3dvar',
        #'sort_cloud_regime',
        #'sort_cloud_3regime',
        #'RadKernel_regime',
        ]
    else:
        cal_types = [
        #'RadFeedback',
        #'RadKernel',
        #'Webb_Decomp',
        #'cal_LCF',
        #'cal_cloud',
        #'cal_EIS',
        #'cal_pr',
        'cal_3dvar',
        #'sort_cloud_3regime',
        #'RadKernel_regime',
        ]
    
    # current directory
    curdir = os. getcwd()
    
    # set CloudRadKernel input kernel directory
    CloudRadKernel_dir = curdir+'/CloudRadKernel_input/'
    
    # set final output directory
    #outdir_final = curdir+'/data/'
    outdir_final = curdir+'/data_ppe_2010to'+str(yrE)+'/'
    
    # set output figure directory
    #figdir = curdir+'/figure/'
    figdir = curdir+'/figure_ppe/'
    
    # set the case tag for control and warming experiments. Dont modify it.
    exp1 = 'FC5'
    exp2 = 'FC5_4K'
    
    # ---------Create Directories--------------------------------------------
    for outdir in [outdir_out, outdir_final, figdir]:
        AP.make_dir(outdir)
    
    # ----------------------------------------------------------------------------
    for icase,case in enumerate(case_short):
    
    
        if '_FR' in case:
            yearS_CTL = 2010
            yearE_CTL = 2018
            yearS_P4K = 2010
            yearE_P4K = 2018
            run_id1 = CLP.get_lutable(case,'PD')
            run_id2 = CLP.get_lutable(case,'P4K')
        else:
            yearS_CTL = 2010
            yearE_CTL = yrE #2014
            yearS_P4K = 2010
            yearE_P4K = yrE #2014
            run_id1 = CLP.get_lutable(case,'PD_FR')
            run_id2 = CLP.get_lutable(case,'P4K_FR')
    
        print(case)
        print(run_id1,yearS_CTL,yearE_CTL)
        print(run_id2,yearS_P4K,yearE_P4K)
    
        #################################################################
        # pre-process model output to get necessary input files
        #################################################################
        if PreProcess:
            # Part -1: archive data from run directory to archive/atm/hist directory
            datadir_in0 = 'case_scripts/'
            #os.system('sh archive_data.sh '+run_id1 + ' '+run_id2+' '+datadir_in1+' '+datadir_in0)
        
#            os.system('sh get_data_select.sh '+ run_id1 + ' '+run_id2+' '+rgr_map+' '+str(yearS_CTL)+' '+str(yearE_CTL)+' '+str(yearS_P4K)+' '+str(yearE_P4K)+' '+datadir_in1+' '+datadir_in2+' '+outdir_out+' '+comp)

            os.system('sh get_data_spec.sh '+ run_id1 + ' '+run_id2+' '+rgr_map+' '+str(yearS_CTL)+' '+str(yearE_CTL)+' '+str(yearS_P4K)+' '+str(yearE_P4K)+' '+datadir_in1+' '+datadir_in2+' '+outdir_out+' '+comp)

    
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
        # -------------------------
        # If you are not on cori or LC machines, please set it as False. The comparison with other CMIP models are not 
        # supported on other machines currently. 
        Add_otherCMIPs = False ## include option about whether adding results from other CMIP models
    
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
            'APRP_latlon',                       # lat-lon plot of APRP results
    #        'CldRadKernel_latlon_dif',          # lat-lon plot of CldRadKernel feedback difference between case and reference case
    #        'RadKernel_latlon_dif',             # lat-lon plot of RadKernel feedback difference between case and reference case
            'tas_latlon',                       # lat-lon plot of surface air temperature and the difference between case and reference case
            'LCF',                              # Temperature - Liquid Condensate Fraction
    #        'zm_CLOUD',                         # zonal mean plot of cloud varaibles difference 
    #        'latlon_CLOUD',                     # lat-lon plot of cloud varaibles difference
            #'webb_decomp',                      # decomposition of adjusted CRE feedback into low and non-low clouds
            #'CLOUD_profile',                    # vertical cloud profile in different regions
            #'NRMSE_RadKern',                   # NRMSE and spatial correlation (COR) evolution with incremental denial experiments [note: RadKernel_latlon_dif should run first.]
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
            #'CLOUD_profile',                    # vertical cloud profile in different regions
            #'NRMSE_RadKern',                   # NRMSE and spatial correlation (COR) evolution with incremental denial experiments [note: RadKernel_latlon_dif should run first.]
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
        #datadir_v2 = datadir+'/data/'
        #datadir_v2 = datadir+'/data_ppe/'
        datadir_v2 = outdir_final
       
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
