#!/usr/bin/env python
# coding: utf-8

'''
This script is used to plot all kinds of figures for E3SMv2 developing versions via comparing with
E3SMv1-amip, E3SMv1-piControl, and other available CMIP5 and CMIP6 models. 

created by: Yi Qin on Aug 4, 2020
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

## Modification starts here ################################################################################################################
# please set all these following directories and your prefered styles: start ----------------------

# -------------- set up directories for all kinds of necessary data
# main directory. pls modify it based on your current script directory. 
datadir = '/global/homes/q/qinyi/E3SM_scripts/diag_feedback_E3SM/'

# -- data for E3SMv2 [it includes all data that you want to be plotted. If main.py runs successfully, this directory would be enough for further plot.]
datadir_v2 = datadir+'data/'
# figure directory
figdir = datadir+'figure/'

cases = ['v1_coupled','amip-p4K','amip-future4K']
colors = ['tab:red','tab:blue','tab:orange']
linewidths = [2, 2, 2]
linestyles = ['-','-','-']

#cases = ['v1_coupled','v1_amip', 'alpha3_0','alpha2', 'alpha2_4', 'alpha2_min10nc', 'alpha2_v1p']
#colors = ['black', 'brown','red','m','b','g','c']
#linewidths = [2, 2, 2, 2, 2, 2, 2]
#linestyles = ['-','-','-','-.','-.','-.','-']

# --- include option about whether adding results from other CMIP models 
#Add_otherCMIPs = False

# please set all these following directories and your prefered styles: End!!!!! ----------------------

#----------------- choose which type figures you want to plot.
# scatter plot: global mean CRE feedback
plot_CRE_globalmean = True
# scatter plot: global mean RadKernel feedback --> non-cloud feedback and adjusted CRE feedback
plot_RadKernel_globalmean = True
# scatter plot: global mean CldRadKernel feedback --> decomposition of cloud feedback
plot_RadKernel_zonalmean = True
# zonal mean plot: RadKernel feedback --> adjusted CRE feedback
plot_CldRadKernel_globalmean = True
# zonal mean plot: CldRadKernel feedback --> decomposition of cloud feedback
plot_CldRadKernel_zonalmean = True
# scatter plot: amip-p4K vs Future-4K -- just CRE feedback
plot_CRE_globalmean_P4KvsFuture = True



# ----------- optional settings ----------------
# marker size for E3SMv2
s1 = 200
# marker size for other CMIP models
s2 = 150
# apparency for markers
a1 = 0.6

# advanced setting: iff plot_CldRadKernel_zonalmean = True, control the number of cases you want to show in one figure.
# for example, if you would like to show the first three cases, then first 6 cases, and all cases, pls set ncase = [3,6,7]
# generally, if you want all lines in the same plot, just set ncase = [len(cases)]
#ncase = [3,6,len(cases)]
ncase = [len(cases)]

# ---------- optional settings ----------------
####Modification ends here ##############################################################################################################

#------------------------------------------------------------------------------------------
# ------------------------- set up directories for necessary data --------------------------
datadir_CMIPs = datadir+'plot_data/'
# -- data for E3SMv1 [dont modify data in this directory.]
datadir_v1 = datadir_CMIPs+'E3SMv1_data/'
# -- data for other CMIP models from CRE feedback [dont' modify data in it.]
datadir_Ringer = datadir_CMIPs+'Ringer2014_forcing_update/'
# -- data for other CMIP models for RadKernel feedback [don't modify data in it.]
datadir_RadKernel = datadir_CMIPs+'RadKernel-Jul6/'
# -- data for other CMIP models for CldRadKernel feedabck [ don't modify data in it.]
datadir_CldRadKernel = datadir_CMIPs+'CloudRadKernel-2/'

####################################################################################
### 1. bar plot for global mean CRE feedback: including E3SMv1 piControl and amip
####################################################################################

Add_otherCMIPs = True

if plot_CRE_globalmean:
    print('ScatterPlot-CRE-feedback starts ........')

    df_all = pd.DataFrame()

    if Add_otherCMIPs:
        # read other CMIP5 and CMIP6 models

	    # here use E3SM-FC5 as the fake of E3SM-1-0 in AMIP simulation
#        models_cmip6 = ['BCC-CSM2-MR','CNRM-CM6-1','IPSL-CM6A-LR',\
#        'MRI-ESM2-0','GISS-E2-1-G','CESM2','GFDL-CM4','CanESM5']
#        models_cmip5 = ['bcc-csm1-1','CNRM-CM5','IPSL-CM5A-LR','IPSL-CM5B-LR',\
#        'FGOALS-g2','MIROC5','MPI-ESM-LR','MPI-ESM-MR','MRI-CGCM3','CanAM4']

        # Oct 19, 2020: reduce model lists to fit both amip-p4K and amipFuture
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
        elif case == 'v1_amip':
            # read v1-amip
            df_amip = pd.read_csv(datadir_v1+'global_mean_features_amip_E3SM-1-0.csv',index_col=0)
            df_amip.index = df_coupled.loc[:,'varname']
            df_all['v1_amip'] = df_amip.loc[:,'anomaly_perK']
        else:   
            df1 = pd.read_csv(datadir_v2+'global_mean_features_'+case+'_E3SM-1-0'+'.csv',index_col=0)
            df1.index = df1.loc[:,'varname']
    
            df2 = df1.loc[:,'anomaly_perK']
            df_all[case] = df2
    
    # start plotting 
    
    fig = plt.figure(figsize=(12,12))
    ax = fig.add_subplot(1,1,1)
    fh = 15
    a1 = 0.6
    
    
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
            elif column == 'v1_amip':
                ax.scatter(x[idx],df_plot.loc[index,column].tolist(),s=s1,alpha=a1,label=column,color=colors[icol],marker='x')
            else:
                ax.scatter(x[idx],df_plot.loc[index,column].tolist(),s=s1,alpha=a1,label=column,color=colors[icol])
    
        if Add_otherCMIPs:
            # add other CMIP models
            for icol,column in enumerate(df_p4K_plot.columns):
                ax.scatter(x[idx]-0.2,df_p4K_plot.loc[index,column].tolist(),edgecolor='none',facecolor='tab:blue',alpha=a1,s=s2)
            # ensemble mean
            ax.scatter(x[idx]-0.2,df_p4K_plot.loc[index,:].mean().tolist(),edgecolor='black',facecolor='blue',s=s2)

            for icol,column in enumerate(df_future_plot.columns):
                ax.scatter(x[idx]+0.2,df_future_plot.loc[index,column].tolist(),edgecolor='none',facecolor='tab:orange',alpha=a1,s=s2)
            # ensemble mean
            ax.scatter(x[idx]+0.2,df_future_plot.loc[index,:].mean().tolist(),edgecolor='black',facecolor='orange',s=s2)

            
        ax.tick_params(labelsize=fh)
        ax.set_ylabel('W/m$^2$/K',fontsize=fh)
        if idx==0:
            ax.legend(fontsize=fh)
    
    plt.xticks(x,df_plot.index)
    
    ax.grid(which='major', linestyle=':', linewidth='1.0', color='grey')
    fig.savefig(figdir+'ScatterPlot-CRE-feedback.pdf')
    plt.show()
    del(df_all,df_plot)
    print('ScatterPlot-CRE-feedback is done!')
    exit()

####################################################################################
### 1. scatter plot for global mean CRE feedback: p4K vs future4K
####################################################################################

if plot_CRE_globalmean_P4KvsFuture:
    print('ScatterPlot-CRE-feedback P4K.vs.Future starts ........')

    df_all = pd.DataFrame()

    if Add_otherCMIPs:
        # read other CMIP5 and CMIP6 models

	    # here use E3SM-FC5 as the fake of E3SM-1-0 in AMIP simulation
#        models_cmip6 = ['BCC-CSM2-MR','CNRM-CM6-1','IPSL-CM6A-LR',\
#        'MRI-ESM2-0','GISS-E2-1-G','CESM2','GFDL-CM4','CanESM5']
#        models_cmip5 = ['bcc-csm1-1','CNRM-CM5','IPSL-CM5A-LR','IPSL-CM5B-LR',\
#        'FGOALS-g2','MIROC5','MPI-ESM-LR','MPI-ESM-MR','MRI-CGCM3','CanAM4']

        # Oct 19, 2020: reduce model lists to fit both amip-p4K and amipFuture
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
        elif case == 'v1_amip':
            # read v1-amip
            df_amip = pd.read_csv(datadir_v1+'global_mean_features_amip_E3SM-1-0.csv',index_col=0)
            df_amip.index = df_coupled.loc[:,'varname']
            df_all['v1_amip'] = df_amip.loc[:,'anomaly_perK']
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

        print(df_p4K_plot)
        print(df_future_plot)

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
                print(index,column)
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
    print('ScatterPlot-CRE-P4KvsFuture-feedback is done!')



####################################################################
### 2. bar plot for radiative feedback based on Radiative kernel 
####################################################################

Add_otherCMIPs = False
if plot_RadKernel_globalmean:
    print('ScatterPlot-RadKernel-Feedback starts........')

    df_all = pd.DataFrame()

    # read other CMIP5&6 models
    if Add_otherCMIPs:
        df_others = pd.DataFrame()
        phases = ['CMIP5','CMIP6']
        cmip5_models = ['bcc-csm1-1', 'CCSM4', 'CNRM-CM5', 'IPSL-CM5A-LR', 'IPSL-CM5B-LR', 'MIROC5', 'MPI-ESM-LR', 'MPI-ESM-MR', 'MRI-CGCM3']
        cmip6_models = ['BCC-CSM2-MR', 'CESM2', 'CNRM-CM6-1', 'GFDL-CM4', 'GISS-E2-1-G','IPSL-CM6A-LR', 'MRI-ESM2-0']
        models = [cmip5_models, cmip6_models]
        exps = ['amip', 'piControl']

        for iphase,phase in enumerate(phases):
            for imodel,model in enumerate(models[iphase]):
                df = pd.read_csv(datadir_RadKernel+'FDBK_'+phase+'_'+exps[0]+'_'+model+'.csv',index_col=0)
                df2 = df.iloc[:26,0]
                df_others[model] = df2

    # E3SM
    for icase,case in enumerate(cases):
        if case == 'v1_coupled':
            # read v1-coupled 
            df_coupled = pd.read_csv(datadir_v1+'FDBK_CMIP6_abrupt-4xCO2_E3SM-1-0_Latest-Oct18_1yr-150yr.csv',index_col=0)
            df_all['v1_coupled'] = df_coupled.iloc[:,0]
        elif case == 'v1_amip':
            # read v1-amip
            df_amip = pd.read_csv(datadir_v1+'FDBK_CMIP6_amip_E3SM-1-0.csv',index_col=0)
            df_all['v1_amip'] = df_amip.iloc[:26,:]
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
            elif column == 'v1_amip':
                ax.scatter(x[idx],df_plot.loc[index,column].tolist(),s=s1,alpha=a1,label=column,color=colors[icol],marker='x')
            else:
                ax.scatter(x[idx],df_plot.loc[index,column].tolist(),s=s1,alpha=a1,label=column,color=colors[icol])
            
        # other CMIP models
        if Add_otherCMIPs:
            for icol,column in enumerate(df_others_plot.columns):
                ax.scatter(x[idx]-0.2, df_others_plot.loc[index,column].tolist(),s=s2,edgecolor='none',facecolor='grey',alpha=a1)
            # ensemble mean
            ax.scatter(x[idx]-0.2, df_others_plot.loc[index,:].mean(),s=s2,edgecolor='black',facecolor='black')

        ax.tick_params(labelsize=fh)
        ax.set_ylabel('W/m$^2$/K',fontsize=fh)
        if idx == 0:
            ax.legend(fontsize=fh)

    ax.grid(which='major', linestyle=':', linewidth='1.0', color='grey')
    degrees = 70
    plt.xticks(x,df_plot.index,rotation=degrees)
    ax.set_title('CRE feedback',fontsize=fh)
    
    fig.savefig(figdir+'ScatterPlot-RadKernel-Feedback.pdf')
    print('ScatterPlot-RadKernel-Feedback is done!')


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
            cmip5_models = ['bcc-csm1-1', 'CCSM4', 'CNRM-CM5', 'IPSL-CM5A-LR', 'IPSL-CM5B-LR', 'MIROC5', 'MPI-ESM-LR', 'MPI-ESM-MR', 'MRI-CGCM3']
            cmip6_models = ['BCC-CSM2-MR', 'CESM2', 'CNRM-CM6-1', 'GFDL-CM4', 'GISS-E2-1-G','IPSL-CM6A-LR', 'MRI-ESM2-0']
            models = [cmip5_models, cmip6_models]
            exps = ['amip', 'piControl']
            model_list = cmip5_models + cmip6_models

        for ivar,svar in enumerate(variables):

            # get other CMIP models
            if Add_otherCMIPs:
                data_others = np.zeros((nlat,len(model_list)))
                for iphase,phase in enumerate(phases):
                    # define the array to save the overall data for cross-model correlation
                    data1 = np.zeros((nlat,len(models[iphase])))
                    for imodel,model in enumerate(models[iphase]):
                        f1 = cdms.open(datadir_RadKernel+'lat-lon-gfdbk-'+phase+'-'+exps[0]+'-'+model+'.nc')
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
                elif case == 'v1_amip':
                    f1 = cdms.open(datadir_v1+'lat-lon-gfdbk-CMIP6-amip-E3SM-1-0.nc')
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
    
                Needs = [-75,-45,-15,0, 15,45,75]
    
                N = [i for i in range(len(lats)) if lats[i] in Needs]
                spec_lats = Needs
                spec_clats = list(np.array(clats)[N])
                ######################################################
    
                # get zonal mean
                data_all[:,icase] = MV.average(data,axis=1)    
    
            # start plotting ...
            ax = fig.add_subplot(2,2,num1+1)
    
    #        ax.set_prop_cycle(color=colors,lw=[2,2,2,2,2,2,2],linestyle=['-','-','-','-.','-.','-.','-'])
            ax.set_prop_cycle(color=colors,lw=linewidths,linestyle=linestyles)
    
            L1 = ax.plot(lats,data_all,alpha=a1)
            
            # plot other CMIP models
            if Add_otherCMIPs:
                L2 = ax.plot(lats,np.average(data_others,axis=1),lw=3,label='ENS-MEAN',color='grey',linestyle='-')
                ax.fill_between(lats, np.amax(data_others,axis=1),np.amin(data_others,axis=1),alpha=0.2)
    
        #     plt.xticks(spec_clats,spec_lats,fontsize=fh)
            plt.yticks(fontsize=fh)
    
            ax.set_title(variables_out[ivar],fontsize=fh)
            ax.set_ylabel('W/m$^2$/K',fontsize=fh)
    
            ax.tick_params(axis='both', which='major', labelsize=fh)
            ax.set_xlim((-90,90))
    
        #     ax.set_ylim((-3,2))
            ax.grid(which='major', linestyle=':', linewidth='1.0', color='grey')
    
            ax.axhline(y=0,color='grey',linestyle='--',lw=2)
            if ivar == len(variables)-1:
                ax.legend(L1,cases,fontsize=fh,bbox_to_anchor=(1.04,0), loc='lower left')
    
            num1 += 1
    
        plt.tight_layout()
        fig.savefig(figdir+'Zonal-mean-Cloud-RadKernel-Feedback-'+str(np.round(ii,0))+'.pdf')

    print('plot_RadKernel_zonalmean is done!')

###################################################################
### 4. bar plot of cloud feedback based on cloud radiative kernel
###################################################################

if plot_CldRadKernel_globalmean:

    print('ScatterPlot-Cloud-feedback-Decomposition starts .........')
    # other CMIP models
    if Add_otherCMIPs:
        phases = ['CMIP5','CMIP6']
        cmip5_models = ['MIROC5', 'MPI-ESM-LR', 'MRI-CGCM3']
        cmip6_models = ['CESM2','E3SM-1-0', 'GFDL-CM4', 'IPSL-CM6A-LR', 'MRI-ESM2-0']
        models = [cmip5_models, cmip6_models]
        exps = ['amip', 'piControl']

        df_LW_others = pd.DataFrame()
        df_SW_others = pd.DataFrame()

        for iphase,phase in enumerate(phases):
            for imodel,model in enumerate(models[iphase]):
                df = pd.read_csv(datadir_CldRadKernel+'decomp_global_mean_lw_'+phase+'_'+exps[0]+'_'+model+'.csv',index_col=0)
                df2 = df.loc[:,model]
                df_LW_others[model] = df2

                df = pd.read_csv(datadir_CldRadKernel+'decomp_global_mean_sw_'+phase+'_'+exps[0]+'_'+model+'.csv',index_col=0)
                df2 = df.loc[:,model]
                df_SW_others[model] = df2

    # E3SM
    df_LW_all = pd.DataFrame()
    df_SW_all = pd.DataFrame()
    
    for icase,case in enumerate(cases):
        if case == 'v1_coupled':
            df1 = pd.read_csv(datadir_v1+'decomp_global_mean_lw_CMIP6_piControl_E3SM-1-0.csv',index_col=0)
            df2 = df1.loc[:,'E3SM-1-0']
            df_LW_all[case] = df2
            
            df1 = pd.read_csv(datadir_v1+'decomp_global_mean_sw_CMIP6_piControl_E3SM-1-0.csv',index_col=0)
            df2 = df1.loc[:,'E3SM-1-0']
            df_SW_all[case] = df2
        elif case == 'v1_amip':
            df1 = pd.read_csv(datadir_v1+'decomp_global_mean_lw_CMIP6_amip_E3SM-1-0.csv',index_col=0)
            df2 = df1.loc[:,'E3SM-1-0']
            df_LW_all[case] = df2
            
            df1 = pd.read_csv(datadir_v1+'decomp_global_mean_sw_CMIP6_amip_E3SM-1-0.csv',index_col=0)
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
            elif column == 'v1_amip':
                axes[ii].scatter(x-w+w2,y1.values.tolist(),marker='x',s=s1,color=colors[icol],linewidths=wf,alpha=a1,label=column)
            else:
                L1 = axes[ii].scatter(x-w+w2,y1.values.tolist(),marker='o',s=s2,color=colors[icol],linewidths=wf,alpha=a1,label=column)

               
        if ii == 0:
            axes[ii].legend(fontsize=fh)
            
        for icol,column in enumerate(df_net_all.columns):
            y1 = df_net_all.iloc[ii*5:(ii+1)*5,icol]
            if column == 'v1_coupled':
                axes[ii].scatter(x+w2,y1.values.tolist(),marker='x',s=s1,color=colors[icol],linewidths=wf,alpha=a1,label=column)
            elif column == 'v1_amip':
                axes[ii].scatter(x+w2,y1.values.tolist(),marker='x',s=s1,color=colors[icol],linewidths=wf,alpha=a1,label=column)
            else:
                L2 = axes[ii].scatter(x+w2,y1.values.tolist(),marker='o',s=s2,color=colors[icol],linewidths=wf,alpha=a1,label=column)
    
        for icol,column in enumerate(df_SW_all.columns):
            y1 = df_SW_all.iloc[ii*5:(ii+1)*5,icol]
            if column == 'v1_coupled':
                axes[ii].scatter(x+w+w2,y1.values.tolist(),marker='x',s=s1,color=colors[icol],linewidths=wf,alpha=a1,label=column)
            elif column == 'v1_amip':
                axes[ii].scatter(x+w+w2,y1.values.tolist(),marker='x',s=s1,color=colors[icol],linewidths=wf,alpha=a1,label=column)
            else:
                L3 = axes[ii].scatter(x+w+w2,y1.values.tolist(),marker='o',s=s2,color=colors[icol],linewidths=wf,alpha=a1,label=column)
    
        # CMIP - other models
        if Add_otherCMIPs:
            a2 = 0.3
            s3 = 100
            for icol,column in enumerate(df_LW_others.columns):
                y1 = df_LW_others.iloc[ii*5:(ii+1)*5,icol]
                axes[ii].scatter(x-w+w2-w3,y1.values.tolist(),marker='o',s=s3,color='grey',linewidths=wf,alpha=a2,label='_nolegend_')

            for icol,column in enumerate(df_net_others.columns):
               y1 = df_net_others.iloc[ii*5:(ii+1)*5,icol]
               axes[ii].scatter(x+w2-w3,y1.values.tolist(),marker='o',s=s3,color='grey',linewidths=wf,alpha=a2,label='_nolegend_')
      
            for icol,column in enumerate(df_SW_others.columns):
               y1 = df_SW_others.iloc[ii*5:(ii+1)*5,icol]
               axes[ii].scatter(x+w+w2-w3,y1.values.tolist(),marker='o',s=s3,color='grey',linewidths=wf,alpha=a2,label='_nolegend_')
            
            axes[ii].scatter(x-w+w2-w3,df_LW_others.iloc[ii*5:(ii+1)*5,:].mean(axis=1),marker='o',s=s3,color='black',linewidths=wf,alpha=1.0,label='_nolegend_')
            axes[ii].scatter(x+w2-w3,df_net_others.iloc[ii*5:(ii+1)*5,:].mean(axis=1),marker='o',s=s3,color='black',linewidths=wf,alpha=1.0,label='_nolegend_')
            axes[ii].scatter(x+w+w2-w3,df_SW_others.iloc[ii*5:(ii+1)*5,:].mean(axis=1),marker='o',s=s3,color='black',linewidths=wf,alpha=1.0,label='_nolegend_')

        axes[ii].grid(which='major', linestyle=':', linewidth='1.0', color='grey')
    
        plt.legend((L1,L2,L3),["LW","NET","SW"],scatterpoints=1,bbox_to_anchor=(1,1),loc="best",fontsize=15)
        
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

    print('ScatterPlot-Cloud-feedback-Decomposition is done!')

#####################################################################
### 5. zonal mean cloud feedback based on cloud radiative kernel
#####################################################################
if plot_CldRadKernel_zonalmean:

    print('ZonalMean-Cloud-feedback-Decomposition starts ........')

    nlat = 90
    nlon = 144
    
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
                        phases = ['CMIP5','CMIP6']
                        cmip5_models = ['MIROC5', 'MPI-ESM-LR', 'MRI-CGCM3']
                        cmip6_models = ['CESM2','E3SM-1-0', 'GFDL-CM4', 'IPSL-CM6A-LR', 'MRI-ESM2-0']
                        models = [cmip5_models, cmip6_models]
                        exps = ['amip', 'piControl']
                        model_list = cmip5_models + cmip6_models
            
                        data_others = np.zeros((nlat,len(model_list)))
                        avgdata_others = np.zeros(len(model_list))

                        for iphase,phase in enumerate(phases):
                            data1 = np.zeros((nlat,len(models[iphase])))
                            avgdata1 = np.zeros((len(models[iphase])))
                            for imodel,model in enumerate(models[iphase]):
                                f1 = cdms.open(datadir_CldRadKernel+'global_cloud_feedback_'+phase+'_'+exps[0]+'_'+model+'.nc')
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
                            f1 = cdms.open(datadir_v1+'global_cloud_feedback_CMIP6_piControl_E3SM-1-0.nc')
                        elif case == 'v1_amip':
                            f1 = cdms.open(datadir_v1+'global_cloud_feedback_CMIP6_amip_E3SM-1-0.nc')
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
    
                        Needs = [-75,-45,-15,15,45,75]
    
                        N = [i for i in range(len(lats)) if lats[i] in Needs]
                        spec_lats = Needs
                        spec_clats = list(np.array(clats)[N])
                        ######################################################
    
                        # get zonal mean
                        data_all[:,icase] = MV.average(data,axis=1)
                        # get global mean
                        avgdata[icase] = cdutil.averager(data,axis='xy',weights='weighted')

                    # start plotting ...
                    ax = fig.add_subplot(2,2,num1+1)
    
                    ax.set_prop_cycle(color=colors,lw=linewidths,linestyle=linestyles)
    
                    L1 = ax.plot(lats,data_all)

                    # plot other CMIP models
                    if Add_otherCMIPs:
                        L2 = ax.plot(lats,np.average(data_others,axis=1),lw=3,label='ENS-MEAN',color='grey',linestyle='-')
                        ax.fill_between(lats, np.amax(data_others,axis=1),np.amin(data_others,axis=1),alpha=0.2)

        #             plt.xticks(spec_clats,spec_lats,fontsize=fh)
                    plt.yticks(fontsize=fh)
                    ax.set_title(varname,fontsize=fh)
                    ax.set_ylabel('W/m$^2$/K',fontsize=fh)
    
        #             ax.set_ylim((-2,2))
                    ax.grid(which='major', linestyle=':', linewidth='1.0', color='grey')
    
                    ax.axhline(y=0,color='grey',linestyle='--',lw=2)
    
                    ax.tick_params(axis='both', which='major', labelsize=fh)
                    ax.set_xlim((-90,90))
    
                    if 'tau' in varname:
                        ax.legend(L1,cases,fontsize=fh,bbox_to_anchor=(1.04,0), loc='lower left')
    
                    num1 += 1
    
                plt.tight_layout()
                fig.savefig(figdir+'ZonalMean-Cloud-feedback-Decomposition-'+lev+'-'+component+'-'+str(np.round(ii,0))+'.pdf')
    print('ZonalMean-Cloud-feedback-Decomposition is done!')

    exit()



