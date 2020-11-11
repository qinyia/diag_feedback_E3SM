#!/usr/bin/env python
# coding: utf-8

'''
This script is used to plot all kinds of figures for E3SMv2 developing versions via comparing with
E3SMv1-amip, E3SMv1-piControl, and other available CMIP5 and CMIP6 models. 

created by: Yi Qin on Aug 4, 2020
modified: Nov 10, 2020

'''
import cdms2 as cdms
import cdutil
import MV2 as MV
import numpy as np
import pylab as pl
import matplotlib as mpl
mpl.use('Agg')
import sys

import cartopy.crs as ccrs
import cartopy
import matplotlib.pyplot as plt
import os
import pandas as pd
import cdtime
import genutil
from matplotlib import ticker
from genutil import statistics
import scipy.stats
import PlotDefinedFunction as PDF


############################################################################################
####Modification ends here #################################################################

# ---------------- please set all these following directories and your prefered styles: start ----------------------

# -------------- set up directories for all kinds of necessary data
## main directory. pls modify it based on your current script directory. 
datadir = '/global/homes/q/qinyi/scripts/diag_feedback_E3SM/'

## data directory for E3SMv2 
## [it includes all data that you want to be plotted. If main.py runs successfully, this directory would be enough for further plot.]
datadir_v2 = datadir+'data/'
## figure directory
figdir = datadir+'figure/'

## notion: if you also want to compare with default E3SM-1-0, please add 'v1_coupled' and 'v1_amip4K' below.
cases = ['v1_coupled','v1_amip4K','TEST']
colors = ['tab:red','tab:blue','tab:orange']
linewidths = [2,2,2]
linestyles = ['-','-','-']

## include option about whether adding results from other CMIP models 
Add_otherCMIPs = True

# ----------------- please set all these following directories and your prefered styles: End!!!!! ----------------------

# ---------------- please set all plot types you want -----------------------------------------------------------------
## choose which type figures you want to plot.
### scatter plot: global mean CRE feedback
plot_CRE_globalmean = True
### scatter plot: global mean RadKernel feedback --> non-cloud feedback and adjusted CRE feedback
plot_RadKernel_globalmean = True
### scatter plot: global mean CldRadKernel feedback --> decomposition of cloud feedback
plot_RadKernel_zonalmean = True
### zonal mean plot: RadKernel feedback --> adjusted CRE feedback
plot_CldRadKernel_globalmean = True
### zonal mean plot: CldRadKernel feedback --> decomposition of cloud feedback
plot_CldRadKernel_zonalmean = True

### scatter plot: amip-p4K vs Future-4K -- just CRE feedback [cannot be used if you don't have amip-future4K!!!!!]
plot_CRE_globalmean_P4KvsFuture = False

# ---------------- please set other optional setting for figure: start -------------------------------------------------

# optional settings starts----------------
# marker size for E3SMv2
s1 = 200
# marker size for other CMIP models
s2 = 100
# apparency for markers
a1 = 0.4

# advanced setting: iff plot_CldRadKernel_zonalmean = True, control the number of cases you want to show in one figure.
# for example, if you would like to show the first three cases, then first 6 cases, and all cases, pls set ncase = [3,6,7]
# generally, if you want all lines in the same plot, just set ncase = [len(cases)]
#ncase = [3,6,len(cases)]
ncase = [len(cases)]

# ---------------- please set other optional setting for figure: end -------------------------------------------------

####Modification ends here #################################################################
############################################################################################

# ----------- set up directories for necessary data --------------------------
datadir_CMIPs = '/global/project/projectdirs/mp193/www/qinyi/DATA/'
# -- data for E3SMv1 [dont modify data in this directory.]
datadir_v1 = datadir_CMIPs+'E3SMv1_data/'
# -- data for other CMIP models from CRE feedback [dont' modify data in it.]
datadir_Ringer = datadir_CMIPs+'RadFeedback/'
# -- data for other CMIP models for RadKernel feedback [don't modify data in it.]
datadir_RadKernel = datadir_CMIPs+'RadKernel/'
# -- data for other CMIP models for CldRadKernel feedabck [ don't modify data in it.]
datadir_CldRadKernel = datadir_CMIPs+'CldRadKernel/'


Add_amipFuture = False

####################################################################################
### 1. bar plot for global mean CRE feedback: including E3SMv1 piControl and amip
####################################################################################

if plot_CRE_globalmean:
    print('ScatterPlot-CRE-feedback starts ........')

    df_all = pd.DataFrame()

    if Add_otherCMIPs:
        # read other CMIP5 and CMIP6 models
        # Oct 19, 2020: reduce model lists to fit both amip-p4K and amipFuture

        exp_cntl = [['piControl','amip'],['piControl','amip']]
        exp_new = [['abrupt4xCO2','amip4K'],['abrupt-4xCO2','amip-p4K']]
        
        prefix = 'global_mean_features'
        suffix1 = '*.csv'
        suffix2 = '*.csv'
        
        models_all,cmip5_models,cmip6_models = PDF.get_intersect_withripf(exp_cntl,exp_new,prefix,suffix1,suffix2,datadir_Ringer)

        exp_cntl = [['piControl','amip'],['piControl','amip']]
        exp_new = [['abrupt4xCO2','amipFuture'],['abrupt-4xCO2','amip-future4K']]
       
        models_all_future,cmip5_models_future,cmip6_models_future = PDF.get_intersect_withripf(exp_cntl,exp_new,prefix,suffix1,suffix2,datadir_Ringer)
        
        #print('models_all',models_all,len(models_all))
        #print('cmip5_models', cmip5_models,len(cmip5_models))
        #print('cmip6_models', cmip6_models,len(cmip6_models))

        #print('models_all_future',models_all_future,len(models_all_future))
        #print('cmip5_models_future', cmip5_models_future,len(cmip5_models_future))
        #print('cmip6_models_future', cmip6_models_future,len(cmip6_models_future))

        # ---- amip4K ---------------------
        df_p4K = pd.DataFrame()
        for model in models_all:
            #print(model)
            if model in cmip5_models:

                if model == 'CanESM2_r1i1p1':
                    model_amip = 'CanAM4_r1i1p1'
                elif model == 'HadGEM2-ES_r1i1p1':
                    model_amip = 'HadGEM2-A_r1i1p1'
                else:
            	    model_amip = model

                filename = datadir_Ringer+'global_mean_features_CMIP5_amip4K_'+model_amip+'.csv'
            else:
                filename = datadir_Ringer+'global_mean_features_CMIP6_amip-p4K_'+model+'.csv'
        
            df = pd.read_csv(filename,index_col=0)
            df.index = df.loc[:,'varname']
            df2 = df.loc[:,'anomaly_perK']
            df_p4K[model] = df2

        # ---- amipFuture ---------------------
        df_future = pd.DataFrame()
        for model in models_all_future:
            if model in cmip5_models_future:

                if model == 'CanESM2_r1i1p1':
                    model_amip = 'CanAM4_r1i1p1'
                elif model == 'HadGEM2-ES_r1i1p1':
                    model_amip = 'HadGEM2-A_r1i1p1'
                else:
            	    model_amip = model

                filename = datadir_Ringer+'global_mean_features_CMIP5_amipFuture_'+model_amip+'.csv'
            else:
                filename = datadir_Ringer+'global_mean_features_CMIP6_amip-future4K_'+model+'.csv'
        
            df = pd.read_csv(filename,index_col=0)
            df.index = df.loc[:,'varname']
            df2 = df.loc[:,'anomaly_perK']
            df_future[model] = df2

    # read amip
    for icase,case in enumerate(cases):
        if case == 'v1_coupled':
            # read v1-coupled 
            df_coupled = pd.read_csv(datadir_v1+'global_mean_features_CMIP6_abrupt-4xCO2_E3SM-1-0_r1i1p1f1.csv',index_col=0)
            df_coupled.index = df_coupled.loc[:,'varname']
            df_all['v1_coupled'] = df_coupled.loc[:,'anomaly_perK']
        elif case == 'v1_amip4K':
            # read v1-amip
            df_amip = pd.read_csv(datadir_v1+'global_mean_features_CMIP6_amip-p4K_E3SM-1-0_r1i1p1f1.csv',index_col=0)
            df_amip.index = df_coupled.loc[:,'varname']
            df_all['v1_amip4K'] = df_amip.loc[:,'anomaly_perK']
        else:   
            df1 = pd.read_csv(datadir_v2+'global_mean_features_'+case+'_E3SM-1-0'+'.csv',index_col=0)
            df1.index = df1.loc[:,'varname']
    
            df2 = df1.loc[:,'anomaly_perK']
            df_all[case] = df2
    

    # start plotting 
    fig = plt.figure(figsize=(12,12))
    ax = fig.add_subplot(1,1,1)
    fh = 15
    
    drop_index = ['ts','SWCLR','LWCLR']
    df_plot = df_all.drop(index=drop_index)

    if Add_otherCMIPs:
        df_p4K_plot = df_p4K.drop(index=drop_index)
        df_future_plot = df_future.drop(index=drop_index)
   
    x = np.arange(1,len(df_plot.index)+1,1)
    
    for idx,index in enumerate(df_plot.index):
        for icol,column in enumerate(df_plot.columns):
            if column == 'v1_coupled':
                L1 = ax.scatter(x[idx],df_plot.loc[index,column].tolist(),s=s1,alpha=a1,label=column,color=colors[icol],marker='x')
            elif column == 'v1_amip4K':
                ax.scatter(x[idx],df_plot.loc[index,column].tolist(),s=s1,alpha=a1,label=column,color=colors[icol],marker='x')
            else:
                ax.scatter(x[idx],df_plot.loc[index,column].tolist(),s=s1,alpha=a1,label=column,color=colors[icol])
    
        if Add_otherCMIPs:
            # add other CMIP models
            for icol,column in enumerate(df_p4K_plot.columns):
#                ax.scatter(x[idx]-0.2,df_p4K_plot.loc[index,column].tolist(),edgecolor='none',facecolor='tab:blue',alpha=a1,s=s2)
                ax.scatter(x[idx]-0.2,df_p4K_plot.loc[index,column].tolist(),edgecolor='none',facecolor='grey',alpha=a1,s=s2)
            # ensemble mean
#            L2 = ax.scatter(x[idx]-0.2,df_p4K_plot.loc[index,:].mean().tolist(),edgecolor='black',facecolor='blue',s=s2)
            L2 = ax.scatter(x[idx]-0.2,df_p4K_plot.loc[index,:].mean().tolist(),color='black',s=s2)

#            Add_amipFuture = True
            if Add_amipFuture:
                for icol,column in enumerate(df_future_plot.columns):
                    ax.scatter(x[idx]+0.2,df_future_plot.loc[index,column].tolist(),edgecolor='none',facecolor='tab:orange',alpha=a1,s=s2)
                # ensemble mean
#                L3 = ax.scatter(x[idx]+0.2,df_future_plot.loc[index,:].mean().tolist(),edgecolor='black',facecolor='orange',s=s2)
                L3 = ax.scatter(x[idx]+0.2,df_future_plot.loc[index,:].mean().tolist(),color='red',s=s2)

            
        ax.tick_params(labelsize=fh)
        ax.set_ylabel('W/m$^2$/K',fontsize=fh)
        if idx==0:
            if Add_amipFuture:
                legend1 = ax.legend([L2,L3],['amip4K','amipFuture'],fontsize=fh,loc='upper left')
                ax.legend(fontsize=fh)
                ax.add_artist(legend1) 
            else:
                legend1 = ax.legend([L2],['amip4K'],fontsize=fh,loc='upper left')
                ax.legend(fontsize=fh)
                ax.add_artist(legend1) 
    
    plt.xticks(x,df_plot.index)
    
    ax.grid(which='major', linestyle=':', linewidth='1.0', color='grey')
    fig.savefig(figdir+'ScatterPlot-CRE-feedback.pdf')
    plt.show()
    del(df_all,df_plot)
    print('------------------------------------------------')
    print('ScatterPlot-CRE-feedback is done!')
    print('------------------------------------------------')

####################################################################################
### 1.1 scatter plot for global mean CRE feedback: p4K vs future4K
####################################################################################

if plot_CRE_globalmean_P4KvsFuture:
    print('ScatterPlot-CRE-feedback P4K.vs.Future starts ........')

    df_all = pd.DataFrame()

    if Add_otherCMIPs:
        # read other CMIP5 and CMIP6 models

        models_cmip6 = ['BCC-CSM2-MR','CNRM-CM6-1','IPSL-CM6A-LR',\
        'MRI-ESM2-0','CESM2','GFDL-CM4','CanESM5']
        models_cmip5 = ['bcc-csm1-1','CNRM-CM5','IPSL-CM5A-LR','IPSL-CM5B-LR',\
        'MIROC5','MPI-ESM-LR','MPI-ESM-MR','MRI-CGCM3','CanAM4']
        models_all = models_cmip6 + models_cmip5


        df_p4K = pd.DataFrame()

        for model in models_all:
            if model in models_cmip5:
                filename = datadir_Ringer+'global_mean_features_CMIP5_amip4K_'+model+'_r1i1p1.csv'
            else:
                if model == 'CNRM-CM6-1':
                    filename = datadir_Ringer+'global_mean_features_CMIP6_amip-p4K_'+model+'_r1i1p1f2.csv'
                elif model == 'CanESM5':
                    filename = datadir_Ringer+'global_mean_features_CMIP6_amip-p4K_'+model+'_r1i1p2f1.csv'
                else:
                    filename = datadir_Ringer+'global_mean_features_CMIP6_amip-p4K_'+model+'_r1i1p1f1.csv'
        
            df = pd.read_csv(filename,index_col=0)
            df.index = df.loc[:,'varname']
            df2 = df.loc[:,'anomaly_perK']
            df_p4K[model] = df2

        df_future = pd.DataFrame()
        for model in models_all:
            if model in models_cmip5:
                filename = datadir_Ringer+'global_mean_features_CMIP5_amipFuture_'+model+'_r1i1p1.csv'
            else:
                if model == 'CNRM-CM6-1':
                    filename = datadir_Ringer+'global_mean_features_CMIP6_amip-future4K_'+model+'_r1i1p1f2.csv'
                elif model == 'CanESM5':
                    filename = datadir_Ringer+'global_mean_features_CMIP6_amip-future4K_'+model+'_r1i1p2f1.csv'
                else:
                    filename = datadir_Ringer+'global_mean_features_CMIP6_amip-future4K_'+model+'_r1i1p1f1.csv'
        
            df = pd.read_csv(filename,index_col=0)
            df.index = df.loc[:,'varname']
            df2 = df.loc[:,'anomaly_perK']
            df_future[model] = df2


    # read amip
    for icase,case in enumerate(cases):
        if case == 'v1_coupled':
            # read v1-coupled 
            df_coupled = pd.read_csv(datadir_v1+'global_mean_features_piControl_E3SM-1-0.csv',index_col=0)
            df_coupled.index = df_coupled.loc[:,'varname']
            df_all['v1_coupled'] = df_coupled.loc[:,'anomaly_perK']
        elif case == 'v1_amip4K':
            # read v1-amip
            df_amip = pd.read_csv(datadir_v1+'global_mean_features_amip_E3SM-1-0.csv',index_col=0)
            df_amip.index = df_coupled.loc[:,'varname']
            df_all['v1_amip4K'] = df_amip.loc[:,'anomaly_perK']
        else:   
            df1 = pd.read_csv(datadir_v2+'global_mean_features_'+case+'_E3SM-1-0'+'.csv',index_col=0)
            df1.index = df1.loc[:,'varname']
    
            df2 = df1.loc[:,'anomaly_perK']
            df_all[case] = df2
    
    # start plotting 
    
    fig = plt.figure(figsize=(12,12))
    fh = 15
    a1 = 0.6
    
    
    drop_index = ['ts','SWCLR','LWCLR']
    df_plot = df_all.drop(index=drop_index)

    if Add_otherCMIPs:
        df_p4K_plot = df_p4K.drop(index=drop_index)
        df_future_plot = df_future.drop(index=drop_index)

    for idx,index in enumerate(df_plot.index):

        if index in ['FTOA']:
    	    valmin = -2
    	    valmax = 0
        elif index in ['SWCRE','LWCRE','netCRE']:
    	    valmin = -1
    	    valmax = 1
        elif index in ['FTOA','FTOACLR']:
            valmin = -2.5
            valmax = -0.5
        elif index in ['FSNT','FSNTCLR']:
            valmin = -0.5
            valmax = 1.5
        elif index in ['FLNT','FLNTCLR']:
            valmin = -3
            valmax = -1

        ax = fig.add_subplot(3,3,idx+1)

        # plot E3SM itself
        ax.scatter(df_plot.loc[index,'amip-p4K'],df_plot.loc[index,'amip-future4K'],s=s1-100,alpha=a1,label='E3SM',color='red',marker='*')

        if Add_otherCMIPs:
            for icol,column in enumerate(df_p4K_plot.columns):
                if model in models_cmip5:
                    L1 = ax.scatter(df_p4K_plot.loc[index,column],df_future_plot.loc[index,column],s=s1-100,alpha=a1,label=column,color='tab:blue')
                else:
                    L1 = ax.scatter(df_p4K_plot.loc[index,column],df_future_plot.loc[index,column],s=s1-100,alpha=a1,label=column,color='tab:red')

        ax.plot([valmin,valmax],[valmin,valmax],ls='--',lw=3,color='grey')

        ax.tick_params(labelsize=fh)
        ax.set_ylabel(index+' Future4K [W/m$^2$/K]',fontsize=fh)
        ax.set_xlabel(index+' P4K [W/m$^2$/K]',fontsize=fh)

        ax.grid(which='major', linestyle=':', linewidth='1.0', color='grey')

#        if idx==0:
#            ax.legend(fontsize=fh)
    
#    plt.xticks(x,df_plot.index)
    
    fig.tight_layout()
    fig.savefig(figdir+'ScatterPlot-CRE-P4KvsFuture-feedback.pdf')
    plt.show()
    del(df_all,df_plot)
    print('------------------------------------------------')
    print('ScatterPlot-CRE-P4KvsFuture-feedback is done!')
    print('------------------------------------------------')


####################################################################
### 2. bar plot for radiative feedback based on Radiative kernel 
####################################################################

if plot_RadKernel_globalmean:
    print('ScatterPlot-RadKernel-Feedback starts........')

    df_all = pd.DataFrame()

    # read other CMIP5&6 models
    if Add_otherCMIPs:

        phases = ['CMIP5','CMIP6']

        exp_cntl = [['piControl','amip'],['piControl','amip']]
        exp_new = [['abrupt4xCO2','amip4K'],['abrupt-4xCO2','amip-p4K']]
        
        prefix = 'FDBK'
        suffix1 = '*1yr-150yr.csv'
        suffix2 = '*.csv'
        
        models_all,cmip5_models,cmip6_models = PDF.get_intersect(exp_cntl,exp_new,prefix,suffix1,suffix2,datadir_RadKernel)
        #print('models_all',models_all,len(models_all))
        #print('cmip5_models',cmip5_models,len(cmip5_models))
        #print('cmip6_models',cmip6_models,len(cmip6_models))

        models = [cmip5_models, cmip6_models]

        df_others = pd.DataFrame()

        for iphase,phase in enumerate(phases):
            if phase == 'CMIP5':
                suffix = '_Latest-Oct18_1yr-27yr'
            else:
                suffix = '_Latest-Oct18_1yr-36yr'

            for imodel,model in enumerate(models[iphase]):

                if model == 'CanESM2':
                    model_amip = 'CanAM4'
                elif model == 'HadGEM2-ES':
                    model_amip = 'HadGEM2-A'
                else:
            	    model_amip = model

                df = pd.read_csv(datadir_RadKernel+'FDBK_'+phase+'_'+exp_new[iphase][1]+'_'+model_amip+suffix+'.csv',index_col=0)
                df2 = df.iloc[:,0]
                df_others[model] = df2


    # E3SM
    for icase,case in enumerate(cases):
        if case == 'v1_coupled':
            # read v1-coupled 
            df_coupled = pd.read_csv(datadir_v1+'FDBK_CMIP6_abrupt-4xCO2_E3SM-1-0_Latest-Oct18_1yr-150yr.csv',index_col=0)
            df_all['v1_coupled'] = df_coupled.iloc[:,0]
        elif case == 'v1_amip4K':
            # read v1-amip
            df_amip = pd.read_csv(datadir_v1+'FDBK_CMIP6_amip-p4K_E3SM-1-0_Latest-Oct18_1yr-36yr.csv',index_col=0)
            df_all['v1_amip4K'] = df_amip.iloc[:,:]
        else:    
            df1 = pd.read_csv(datadir_v2+'FDBK_CMIP6_'+case+'_E3SM-1-0'+'.csv',index_col=0)
            df2 = df1.loc[:,'E3SM-1-0']
            df_all[case] = df2
        
    # start plotting 
    
    fig = plt.figure(figsize=(12,12))
    ax = fig.add_subplot(1,1,1)
    fh = 15
    a1 = 0.6
    
    
    drop_index = ['T','dLW_adj','dSW_adj','dnet_adj','LW_resd','SW_resd','net_resd',\
    'T_clr','Planck_clr','LR_clr','WV_clr','ALB_clr','WV_clr_SW','WV_clr_LW','WV_SW','WV_LW',\
    'SWCRE','LWCRE','netCRE','Planck_clr_fxRH','LR_clr_fxRH','RH_clr','LW_clr_sum','SW_clr_sum',\
    'net_clr_sum','LW_clr_dir','SW_clr_dir','net_clr_dir','LW_cld_sum','SW_cld_sum','net_cld_sum',\
    'LW_cld_dir','SW_cld_dir','net_cld_dir','LW_clr_resd','SW_clr_resd','net_clr_resd']
    

    df_plot = df_all.drop(index=drop_index)
    x = np.arange(1,len(df_plot.index)+1,1)

    if Add_otherCMIPs:
        df_others_plot = df_others.drop(index=drop_index)
    
    for idx,index in enumerate(df_plot.index):
        for icol,column in enumerate(df_plot.columns):
            if column == 'v1_coupled':
                L1 = ax.scatter(x[idx],df_plot.loc[index,column].tolist(),s=s1,alpha=a1,label=column,color=colors[icol],marker='x')
            elif column == 'v1_amip4K':
                ax.scatter(x[idx],df_plot.loc[index,column].tolist(),s=s1,alpha=a1,label=column,color=colors[icol],marker='x')
            else:
                ax.scatter(x[idx],df_plot.loc[index,column].tolist(),s=s1,alpha=a1,label=column,color=colors[icol])
            
        # other CMIP models
        if Add_otherCMIPs:
            for icol,column in enumerate(df_others_plot.columns):
                ax.scatter(x[idx]-0.2, df_others_plot.loc[index,column].tolist(),s=s2,edgecolor='none',facecolor='grey',alpha=a1)
            # ensemble mean
            L2 = ax.scatter(x[idx]-0.2, df_others_plot.loc[index,:].mean(),s=s2,edgecolor='black',facecolor='black')

        ax.tick_params(labelsize=fh)
        ax.set_ylabel('W/m$^2$/K',fontsize=fh)
        if idx == 0:
            legend1 = ax.legend([L2],['amip4K'],fontsize=fh,loc='upper left')
            ax.legend(fontsize=fh)
            ax.add_artist(legend1) 

    ax.grid(which='major', linestyle=':', linewidth='1.0', color='grey')
    degrees = 70
    plt.xticks(x,df_plot.index,rotation=degrees)
    ax.set_title('Radiative Kernel feedback',fontsize=fh)
    
    fig.savefig(figdir+'ScatterPlot-RadKernel-Feedback.pdf')
    print('------------------------------------------------')
    print('ScatterPlot-RadKernel-Feedback is done!')
    print('------------------------------------------------')


#########################################################################
### 3. zonal mean plot of radiative feedback based on radiative kernel
#########################################################################

if plot_RadKernel_zonalmean:
    print('plot_RadKernel_zonalmean starts ..........')

    # variables = ['SWCRE_ano_grd_gfdbk','LWCRE_ano_grd_gfdbk','netCRE_ano_grd_gfdbk',\
    #             'SWCRE_ano_grd_adj_gfdbk','LWCRE_ano_grd_adj_gfdbk','netCRE_ano_grd_adj_gfdbk']
    # variables_out = ['SWCRE','LWCRE','netCRE','SWCRE_adj','LWCRE_adj','netCRE_adj']
    
    variables = ['SWCRE_ano_grd_adj','LWCRE_ano_grd_adj','netCRE_ano_grd_adj']
    variables_out = ['SWCRE_adj','LWCRE_adj','netCRE_adj']
    
    nlat = 73
    nlon = 144
    
    # generate figure based on case categories
    for ii in ncase:
        fig = plt.figure(figsize=(18,12))
        fh = 15
        num1 = 0
        a1 = 0.7

        # add other CMIP models
        if Add_otherCMIPs:

            phases = ['CMIP5','CMIP6']

            exp_cntl = [['piControl','amip'],['piControl','amip']]
            exp_new = [['abrupt4xCO2','amip4K'],['abrupt-4xCO2','amip-p4K']]
            
            prefix = 'FDBK'
            suffix1 = '*1yr-150yr.csv'
            suffix2 = '*.csv'
            
            models_all,cmip5_models,cmip6_models = PDF.get_intersect(exp_cntl,exp_new,prefix,suffix1,suffix2,datadir_RadKernel)
            #print('models_all',models_all,len(models_all))
            #print('cmip5_models',cmip5_models,len(cmip5_models))
            #print('cmip6_models',cmip6_models,len(cmip6_models))
            models = [cmip5_models, cmip6_models]
            model_list = cmip5_models + cmip6_models

        for ivar,svar in enumerate(variables):

            # get other CMIP models
            if Add_otherCMIPs:
                data_others = np.zeros((nlat,len(model_list)))
                for iphase,phase in enumerate(phases):
                    if phase == 'CMIP5':
                        suffix = '_Latest-Oct18_1yr-27yr'
                    else:
                        suffix = '_Latest-Oct18_1yr-36yr'
                    
                    # define the array to save the overall data for cross-model correlation
                    data1 = np.zeros((nlat,len(models[iphase])))
                    for imodel,model in enumerate(models[iphase]):
                        if model == 'CanESM2':
                            model_amip = 'CanAM4'
                        elif model == 'HadGEM2-ES':
                            model_amip = 'HadGEM2-A'
                        else:
                    	    model_amip = model

                        f1 = cdms.open(datadir_RadKernel+'lat-lon-gfdbk-'+phase+'-'+exp_new[iphase][1]+'-'+model_amip+suffix+'.nc')
                        tmp1 = f1(svar)
                    
                        lats = tmp1.getLatitude()[:]
                        data1[:,imodel] = MV.average(tmp1,axis=1)
                    if iphase == 0:
                        data_others[:,:len(cmip5_models)] = data1
                    else:
                        data_others[:,len(cmip5_models):] = data1

 
            # E3SM
            data_all = np.zeros((nlat,len(cases[:ii])))
            for icase,case in enumerate(cases[:ii]):
                if case == 'v1_coupled':
                    f1 = cdms.open(datadir_v1+'lat-lon-gfdbk-CMIP6-abrupt-4xCO2-E3SM-1-0_Latest-Oct18_1yr-150yr.nc')
                elif case == 'v1_amip4K':
                    f1 = cdms.open(datadir_v1+'lat-lon-gfdbk-CMIP6-amip-p4K-E3SM-1-0_Latest-Oct18_1yr-36yr.nc')
                else:
                    f1 = cdms.open(datadir_v2+'lat-lon-gfdbk-CMIP6-'+case+'-E3SM-1-0'+'.nc')
    
                data = f1(svar)
                lats = data.getLatitude()[:]
    
                ######################################################
                # Compute weights and take weighted average over latitude dimension
                clat = np.cos(np.deg2rad(lats))
                clat1 = clat/MV.sum(clat)
                clat1[0] = 0.
    
                clats = np.zeros(len(clat1))
    
                for ilat in range(len(clat1)):
                    clats[ilat] = np.sum(clat1[:ilat])
                clats[0] = 0.
    
                #Needs = [-90,-50,-30,-15,0, 15,30,50,90]
                Needs = [-85, -55, -35, -15, 15, 35, 55, 85]
    
                N = [i for i in range(len(lats)) if lats[i] in Needs]
                spec_lats = Needs
                spec_clats = list(np.array(clats)[N])

#                print('lats=',lats)
#                print('clats=',clats)
#                print('spec_lats=',spec_lats)
#                print('spec_clats=',spec_clats)

                ######################################################
    
                # get zonal mean
                data_all[:,icase] = MV.average(data,axis=1)    
    
            # start plotting ...
            ax = fig.add_subplot(2,2,num1+1)
    
            ax.set_prop_cycle(color=colors,lw=linewidths,linestyle=linestyles)
    
            L1 = ax.plot(clats,data_all,alpha=a1)
            
            # plot other CMIP models
            if Add_otherCMIPs:
                L2 = ax.plot(clats,np.average(data_others,axis=1),lw=3,label='ENS-MEAN',color='grey',linestyle='-')
                ax.fill_between(clats, np.amax(data_others,axis=1),np.amin(data_others,axis=1),alpha=0.2,color='grey')
    
            plt.xticks(spec_clats,spec_lats,fontsize=fh)
            ax.set_xlim((0,1))
            plt.yticks(fontsize=fh)
    
            ax.set_title(variables_out[ivar],fontsize=fh)
            ax.set_ylabel('W/m$^2$/K',fontsize=fh)
    
            ax.tick_params(axis='both', which='major', labelsize=fh)
#            ax.set_xlim((-90,90))
    
        #     ax.set_ylim((-3,2))
            ax.grid(which='major', linestyle=':', linewidth='1.0', color='grey')
    
            ax.axhline(y=0,color='grey',linestyle='--',lw=2)
            if ivar == len(variables)-1:
                ax.legend(L1,cases,fontsize=fh,bbox_to_anchor=(1.04,0), loc='lower left')
    
            num1 += 1
    
        plt.tight_layout()
        fig.savefig(figdir+'Zonal-mean-Cloud-RadKernel-Feedback-'+str(np.round(ii,0))+'.pdf')


    print('------------------------------------------------')
    print('plot_RadKernel_zonalmean is done!')
    print('------------------------------------------------')


###################################################################
### 4. bar plot of cloud feedback based on cloud radiative kernel
###################################################################

if plot_CldRadKernel_globalmean:

    print('ScatterPlot-Cloud-feedback-Decomposition starts .........')
    # other CMIP models
    if Add_otherCMIPs:

        phases = ['CMIP5','CMIP6']
        exp_cntl = [['piControl','amip'],['piControl','amip']]
        exp_new = [['abrupt4xCO2','amip4K'],['abrupt-4xCO2','amip-p4K']]
        
        prefix = 'decomp_global_mean_lw'
        suffix1 = '*1yr-150yr.csv'
        suffix2 = '*.csv'
        
        models_all,cmip5_models,cmip6_models = PDF.get_intersect(exp_cntl,exp_new,prefix,suffix1,suffix2,datadir_CldRadKernel)
        #print('models_all',models_all,len(models_all))
        #print('cmip5_models',cmip5_models,len(cmip5_models))
        #print('cmip6_models',cmip6_models,len(cmip6_models))
        models = [cmip5_models, cmip6_models]
        model_list = cmip5_models + cmip6_models

        df_LW_others = pd.DataFrame()
        df_SW_others = pd.DataFrame()

        for iphase,phase in enumerate(phases):
            for imodel,model in enumerate(models[iphase]):
                df = pd.read_csv(datadir_CldRadKernel+'decomp_global_mean_lw_'+phase+'_'+exp_new[iphase][1]+'_'+model+'.csv',index_col=0)
                df2 = df.loc[:,model]
                df_LW_others[model] = df2

                df = pd.read_csv(datadir_CldRadKernel+'decomp_global_mean_sw_'+phase+'_'+exp_new[iphase][1]+'_'+model+'.csv',index_col=0)
                df2 = df.loc[:,model]
                df_SW_others[model] = df2

    # E3SM
    df_LW_all = pd.DataFrame()
    df_SW_all = pd.DataFrame()
    
    for icase,case in enumerate(cases):
        if case == 'v1_coupled':
            df1 = pd.read_csv(datadir_v1+'decomp_global_mean_lw_CMIP6_abrupt-4xCO2_E3SM-1-0_1yr-150yr.csv',index_col=0)
            df2 = df1.loc[:,'E3SM-1-0']
            df_LW_all[case] = df2
            
            df1 = pd.read_csv(datadir_v1+'decomp_global_mean_sw_CMIP6_abrupt-4xCO2_E3SM-1-0_1yr-150yr.csv',index_col=0)
            df2 = df1.loc[:,'E3SM-1-0']
            df_SW_all[case] = df2
        elif case == 'v1_amip4K':
            df1 = pd.read_csv(datadir_v1+'decomp_global_mean_lw_CMIP6_amip-p4K_E3SM-1-0.csv',index_col=0)
            df2 = df1.loc[:,'E3SM-1-0']
            df_LW_all[case] = df2
            
            df1 = pd.read_csv(datadir_v1+'decomp_global_mean_sw_CMIP6_amip-p4K_E3SM-1-0.csv',index_col=0)
            df2 = df1.loc[:,'E3SM-1-0']
            df_SW_all[case] = df2
        else:
            df1 = pd.read_csv(datadir_v2+'decomp_global_mean_lw_'+case+'_E3SM-1-0'+'.csv',index_col=0)
            df2 = df1.loc[:,'E3SM-1-0']
            df_LW_all[case] = df2
    
            df1 = pd.read_csv(datadir_v2+'decomp_global_mean_sw_'+case+'_E3SM-1-0'+'.csv',index_col=0)
            df2 = df1.loc[:,'E3SM-1-0']
            df_SW_all[case] = df2
        
    # get net cloud feedback 
    df_net_all = df_SW_all + df_LW_all
    if Add_otherCMIPs:
        df_net_others = df_SW_others + df_LW_others
        
    # start plotting 
    
    fig, axes = plt.subplots(nrows=3,ncols=1,figsize=(12,15))
    
    titles = ['All Cloud CTP bins', "Non-Low Cloud CTP bins", "Low Cloud CTP bins"]
    labels = ['Total','Amount','Altitude','Optical Depth','Residual']
    wc = 1 # circle
    wf = 0 # filled 
    
    fh = 15
    a1 = 0.7
    
    w = 0.25
    w2 = 0.08
    w3 = 0.08
    
    
    x = np.asarray([1,2,3,4,5])
    
    for ii in range(3): ## loop for each panel, ii reprents All, Non-Low and Low Cloud 
        for icol,column in enumerate(df_LW_all.columns):
            y1 = df_LW_all.iloc[ii*5:(ii+1)*5,icol]
            if column == 'v1_coupled':
                axes[ii].scatter(x-w+w2,y1.values.tolist(),marker='x',s=s1,color=colors[icol],linewidths=wf,alpha=a1,label=column)
            elif column == 'v1_amip4K':
                axes[ii].scatter(x-w+w2,y1.values.tolist(),marker='x',s=s1,color=colors[icol],linewidths=wf,alpha=a1,label=column)
            else:
#                L1 = axes[ii].scatter(x-w+w2,y1.values.tolist(),marker='o',s=s2,color=colors[icol],linewidths=wf,alpha=a1,label=column)
                axes[ii].scatter(x-w+w2,y1.values.tolist(),marker='o',s=s2,color=colors[icol],linewidths=wf,alpha=a1,label=column)

               
        if ii == 0:
            axes[ii].legend(fontsize=fh)
            
        for icol,column in enumerate(df_net_all.columns):
            y1 = df_net_all.iloc[ii*5:(ii+1)*5,icol]
            if column == 'v1_coupled':
                axes[ii].scatter(x+w2,y1.values.tolist(),marker='x',s=s1,color=colors[icol],linewidths=wf,alpha=a1,label=column)
            elif column == 'v1_amip4K':
                axes[ii].scatter(x+w2,y1.values.tolist(),marker='x',s=s1,color=colors[icol],linewidths=wf,alpha=a1,label=column)
            else:
#                L2 = axes[ii].scatter(x+w2,y1.values.tolist(),marker='o',s=s2,color=colors[icol],linewidths=wf,alpha=a1,label=column)
                axes[ii].scatter(x+w2,y1.values.tolist(),marker='o',s=s2,color=colors[icol],linewidths=wf,alpha=a1,label=column)

    
        for icol,column in enumerate(df_SW_all.columns):
            y1 = df_SW_all.iloc[ii*5:(ii+1)*5,icol]
            if column == 'v1_coupled':
                axes[ii].scatter(x+w+w2,y1.values.tolist(),marker='x',s=s1,color=colors[icol],linewidths=wf,alpha=a1,label=column)
            elif column == 'v1_amip4K':
                axes[ii].scatter(x+w+w2,y1.values.tolist(),marker='x',s=s1,color=colors[icol],linewidths=wf,alpha=a1,label=column)
            else:
 #               L3 = axes[ii].scatter(x+w+w2,y1.values.tolist(),marker='o',s=s2,color=colors[icol],linewidths=wf,alpha=a1,label=column)
                axes[ii].scatter(x+w+w2,y1.values.tolist(),marker='o',s=s2,color=colors[icol],linewidths=wf,alpha=a1,label=column)
   
        # CMIP - other models
        if Add_otherCMIPs:
            a2 = 0.3
            s3 = 100
            for icol,column in enumerate(df_LW_others.columns):
                y1 = df_LW_others.iloc[ii*5:(ii+1)*5,icol]
                axes[ii].scatter(x-w+w2-w3,y1.values.tolist(),marker='o',s=s3,color='red',linewidths=wf,alpha=a2,label='_nolegend_')

            for icol,column in enumerate(df_net_others.columns):
                y1 = df_net_others.iloc[ii*5:(ii+1)*5,icol]
                axes[ii].scatter(x+w2-w3,y1.values.tolist(),marker='o',s=s3,color='grey',linewidths=wf,alpha=a2,label='_nolegend_')
      
            for icol,column in enumerate(df_SW_others.columns):
                y1 = df_SW_others.iloc[ii*5:(ii+1)*5,icol]
                axes[ii].scatter(x+w+w2-w3,y1.values.tolist(),marker='o',s=s3,color='blue',linewidths=wf,alpha=a2,label='_nolegend_')
            
            L1 = axes[ii].scatter(x-w+w2-w3,df_LW_others.iloc[ii*5:(ii+1)*5,:].mean(axis=1),marker='o',s=s3,color='red',linewidths=wf,alpha=1.0,label='_nolegend_')
            L2 = axes[ii].scatter(x+w2-w3,df_net_others.iloc[ii*5:(ii+1)*5,:].mean(axis=1),marker='o',s=s3,color='grey',linewidths=wf,alpha=1.0,label='_nolegend_')
            L3 = axes[ii].scatter(x+w+w2-w3,df_SW_others.iloc[ii*5:(ii+1)*5,:].mean(axis=1),marker='o',s=s3,color='blue',linewidths=wf,alpha=1.0,label='_nolegend_')

            plt.legend((L1,L2,L3),["LW","NET","SW"],scatterpoints=1,bbox_to_anchor=(1,1),loc="best",fontsize=15)


        axes[ii].grid(which='major', linestyle=':', linewidth='1.0', color='grey')
        
        axes[ii].axhline(0, color="grey", linestyle="-",linewidth=2)
        axes[ii].set_ylabel('W/$m^2$/K',fontsize=fh)
        axes[ii].set_title(titles[ii],fontsize=fh)
        axes[ii].set_ylim([-0.5,1.5])
        
        major_ticks = np.arange(-1.0,1.5,0.5)
        minor_ticks = np.arange(-1.0,1.5,0.25)
        axes[ii].set_yticks(major_ticks)
        axes[ii].set_yticks(minor_ticks,minor=True)
        
        axes[ii].tick_params(axis='both', which='major', labelsize=fh)
        
        axes[ii].grid(axis='y',c='grey',linestyle='-.',which = 'major')
        
        if ii == 2:
            plt.xticks(x,['Total','Amount','Altitude','Optical Depth','Residual'])
        else:
            axes[ii].set_xticklabels("")
            
    plt.tight_layout()
    fig.savefig(figdir+'ScatterPlot-Cloud-feedback-Decomposition.pdf')

    print('------------------------------------------------')
    print('ScatterPlot-Cloud-feedback-Decomposition is done!')
    print('------------------------------------------------')


#####################################################################
### 5. zonal mean cloud feedback based on cloud radiative kernel
#####################################################################
if plot_CldRadKernel_zonalmean:

    print('ZonalMean-Cloud-feedback-Decomposition starts ........')

    nlat = 90
    nlon = 144

    if Add_otherCMIPs:
        phases = ['CMIP5','CMIP6']
        exp_cntl = [['piControl','amip'],['piControl','amip']]
        exp_new = [['abrupt4xCO2','amip4K'],['abrupt-4xCO2','amip-p4K']]
        
        prefix = 'decomp_global_mean_lw'
        suffix1 = '*1yr-150yr.csv'
        suffix2 = '*.csv'
        
        models_all,cmip5_models,cmip6_models = PDF.get_intersect(exp_cntl,exp_new,prefix,suffix1,suffix2,datadir_CldRadKernel)
        #print('models_all',models_all,len(models_all))
        #print('cmip5_models',cmip5_models,len(cmip5_models))
        #print('cmip6_models',cmip6_models,len(cmip6_models))
        models = [cmip5_models, cmip6_models]
        model_list = cmip5_models + cmip6_models
 
    
    # --- define variables 
    levs = ['ALL','HI680','LO680',]
    components = ['NET','SW','LW']
    decomps = ['tot','amt','alt','tau']
    
    # generate figure based on case categories

    for ii in ncase:
        for component in components:
            for lev in levs:
                fig = plt.figure(figsize=(18,12))
                num1 = 0
                fh = 15
    
                for decomp in decomps:
                    varname = lev+'_'+component+'cld_'+decomp
                    if component == 'NET':
                        varSW = lev+'_SWcld_'+decomp
                        varLW = lev+'_LWcld_'+decomp
    
                    # add other CMIP models
                    if Add_otherCMIPs:
            
                        data_others = np.zeros((nlat,len(model_list)))
                        avgdata_others = np.zeros(len(model_list))

                        for iphase,phase in enumerate(phases):
                            data1 = np.zeros((nlat,len(models[iphase])))
                            avgdata1 = np.zeros((len(models[iphase])))
                            for imodel,model in enumerate(models[iphase]):
                                f1 = cdms.open(datadir_CldRadKernel+'global_cloud_feedback_'+phase+'_'+exp_new[iphase][1]+'_'+model+'.nc')
                                if component in ['LW','SW']:
                                    tmp1 = f1(varname)
                                else:
                                    dataSW = f1(varSW)
                                    dataLW = f1(varLW)
                                    tmp1 = dataSW + dataLW
                                    tmp1.setAxisList(dataSW.getAxisList())

                                data1[:,imodel] = MV.average(tmp1,axis=1)
                                avgdata1[imodel] = cdutil.averager(tmp1,axis='xy',weights='weighted')
                            if iphase == 0:
                                data_others[:,:len(cmip5_models)] = data1
                                avgdata_others[:len(cmip5_models)] = avgdata1
                            else:
                                data_others[:,len(cmip5_models):] = data1
                                avgdata_others[len(cmip5_models):] = avgdata1
 
                    # E3SM 
                    data_all = np.zeros((nlat,len(cases[:ii])))
                    avgdata = np.zeros(len(cases[:ii]))
    
                    for icase,case in enumerate(cases[:ii]):
                        if case == 'v1_coupled':
                            f1 = cdms.open(datadir_v1+'global_cloud_feedback_CMIP6_abrupt-4xCO2_E3SM-1-0_1yr-150yr.nc')
                        elif case == 'v1_amip4K':
                            f1 = cdms.open(datadir_v1+'global_cloud_feedback_CMIP6_amip-p4K_E3SM-1-0.nc')
                        else:
                            f1 = cdms.open(datadir_v2+'global_cloud_feedback_'+case+'_E3SM-1-0'+'.nc')
    
                        if component in ['LW','SW']:
                            data = f1(varname)
                        else:
    
                            dataSW = f1(varSW)
                            dataLW = f1(varLW)
                            data = dataSW + dataLW
                            data.setAxisList(dataSW.getAxisList())
    
                        lats = data.getLatitude()[:]
    
                        ######################################################
                        # Compute weights and take weighted average over latitude dimension
                        clat = np.cos(np.deg2rad(lats))
                        clat1 = clat/MV.sum(clat)
                        clat1[0] = 0.
    
                        clats = np.zeros(len(clat1))
    
                        for ilat in range(len(clat1)):
                            clats[ilat] = np.sum(clat1[:ilat])
                        clats[0] = 0.
    
#                        Needs = [-90,-50,-30,-15,0, 15,30,50,90]
                        Needs = [-85, -55, -35, -15, 15, 35, 55, 85]

                        N = [i for i in range(len(lats)) if lats[i] in Needs]
                        spec_lats = Needs
                        spec_clats = list(np.array(clats)[N])

#                        print('lats=',lats)
#                        print('clats=',clats)
#                        print('spec_lats=',spec_lats)
#                        print('spec_clats=',spec_clats)

                        ######################################################
    
                        # get zonal mean
                        data_all[:,icase] = MV.average(data,axis=1)
                        # get global mean
                        avgdata[icase] = cdutil.averager(data,axis='xy',weights='weighted')


                    # start plotting ...
                    ax = fig.add_subplot(2,2,num1+1)
    
                    ax.set_prop_cycle(color=colors,lw=linewidths,linestyle=linestyles)
    
                    L1 = ax.plot(clats,data_all)

                    # plot other CMIP models
                    if Add_otherCMIPs:
                        L2 = ax.plot(clats,np.average(data_others,axis=1),lw=3,label='ENS-MEAN',color='grey',linestyle='-')
                        ax.fill_between(clats, np.amax(data_others,axis=1),np.amin(data_others,axis=1),alpha=0.2,color='grey')

                    plt.xticks(spec_clats,spec_lats,fontsize=fh)
                    ax.set_xlim((0,1))

                    plt.yticks(fontsize=fh)
                    ax.set_title(varname,fontsize=fh)
                    ax.set_ylabel('W/m$^2$/K',fontsize=fh)
    
        #             ax.set_ylim((-2,2))
                    ax.grid(which='major', linestyle=':', linewidth='1.0', color='grey')
    
                    ax.axhline(y=0,color='grey',linestyle='--',lw=2)
    
                    ax.tick_params(axis='both', which='major', labelsize=fh)
                    #ax.set_xlim((-90,90))
    
                    if 'tau' in varname:
                        ax.legend(L1,cases,fontsize=fh,bbox_to_anchor=(1.04,0), loc='lower left')
    
                    num1 += 1
    
                plt.tight_layout()
                fig.savefig(figdir+'ZonalMean-Cloud-feedback-Decomposition-'+lev+'-'+component+'-'+str(np.round(ii,0))+'.pdf')

    print('------------------------------------------------')
    print('ZonalMean-Cloud-feedback-Decomposition is done!')
    print('------------------------------------------------')





