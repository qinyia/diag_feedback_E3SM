
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
import copy

import cal_global_radiation_feedback_E3SM as GRF
import cal_RadKernel_E3SM as RK
import cal_CloudRadKernel_E3SM as CRK
import cal_LCF_E3SM as LCF
import cal_cloud_E3SM as CLOUD
import cal_webb_decomposition as WD


#############################################################################################
def make_dir(outdir):
    if os.path.isdir(outdir):
        print("----------"+outdir+' exists.-------------')
    else:
        os.mkdir(outdir)
        print("----------Successfully created the directory %s " % outdir)

#############################################################################################

def get_cal_dics(direc_data, case_stamp, yearS2, yearE2, run_id1, run_id2, outdir_final,
                  RadKernel_dir, figdir, exp1, exp2, 
                  CloudRadKernel_dir):

    dics_cal = {}

    print(direc_data, case_stamp, yearS2, yearE2, run_id1, run_id2, outdir_final)

    my_cal = calculation(direc_data, case_stamp, yearS2, yearE2, run_id1, run_id2, outdir_final,
                  RadKernel_dir, figdir, exp1, exp2,
                  CloudRadKernel_dir)

    dics_cal['RadFeedback']         = my_cal.cal_Global_RadFeedback
    dics_cal['RadKernel']           = my_cal.cal_RadKernel
    dics_cal['Webb_Decomp']         = my_cal.cal_webb_decomp
    dics_cal['CloudRadKernel']      = my_cal.cal_CloudRadKernel
    dics_cal['cal_LCF']             = my_cal.cal_LCF
    dics_cal['cal_cloud']           = my_cal.cal_cloud

    return dics_cal

class calculation:
    def __init__(self,direc_data, case_stamp, yearS2, yearE2, run_id1, run_id2, outdir_final,
                      RadKernel_dir, figdir, exp1, exp2, 
                      CloudRadKernel_dir ):

        self.direc_data = direc_data
        self.case_stamp = case_stamp
        self.yearS2 = yearS2
        self.yearE2 = yearE2
        self.run_id1 = run_id1
        self.run_id2 = run_id2
        self.outdir_final = outdir_final
        self.RadKernel_dir = RadKernel_dir
        self.figdir = figdir
        self.exp1 = exp1
        self.exp2 = exp2
        self.CloudRadKernel_dir = CloudRadKernel_dir 

    def cal_Global_RadFeedback(self):
        result = GRF.Global_RadFeedback(self.direc_data, self.case_stamp, self.yearS2, self.yearE2, self.run_id1, self.run_id2, self.outdir_final,self.exp1,self.exp2)
        
    def cal_RadKernel(self):
        result = RK.RadKernel(self.RadKernel_dir,self.direc_data,self.case_stamp,self.yearS2,self.yearE2,self.run_id1,self.run_id2,self.outdir_final,self.figdir,self.exp1,self.exp2)

    def cal_webb_decomp(self):
        result = WD.cal_webb_decomp(self.outdir_final,self.case_stamp,self.yearS2,self.yearE2,self.outdir_final,self.figdir)

    def cal_CloudRadKernel(self):
        result = CRK.CloudRadKernel(self.CloudRadKernel_dir,self.direc_data,self.case_stamp,self.yearS2,self.yearE2,self.run_id1,self.run_id2,self.outdir_final,self.figdir)

    def cal_LCF(self):
        result = LCF.cal_LCF(self.direc_data,self.case_stamp,self.yearS2,self.yearE2,self.run_id1,self.run_id2,self.outdir_final,self.figdir)

    def cal_cloud(self):
        result = CLOUD.cal_cloud(self.direc_data,self.case_stamp,self.yearS2,self.yearE2,self.run_id1,self.run_id2,self.outdir_final,self.figdir,self.exp1,self.exp2)


#############################################################################################

def get_plot_dics(cases,ref_casesA,Add_otherCMIPs,datadir_v2, datadir_v1, s1, s2, fh, fh1, a1, colors, figdir,ncase, linestyles, linewidths,Add_amipFutue,highlight_CESM2,lw_CESM2,ls_CESM2,lc_CESM2, datadir_Ringer, datadir_RadKernel, datadir_CldRadKernel):
    '''
    Aug 20, 201: get the plot dictionary for all plot types.
    '''
    dics_plots = {}

    my_plot = plots(cases,ref_casesA,Add_otherCMIPs,datadir_v2, datadir_v1, s1, s2, fh, fh1, a1, colors, figdir,ncase, linestyles, linewidths,Add_amipFutue,highlight_CESM2,lw_CESM2,ls_CESM2,lc_CESM2, datadir_Ringer, datadir_RadKernel, datadir_CldRadKernel)

    dics_plots['CRE_globalmean']             = my_plot.plot_CRE_globalmean
    dics_plots['RadKernel_globalmean']       = my_plot.plot_RadKernel_globalmean
    dics_plots['CldRadKernel_globalmean']    = my_plot.plot_CldRadKernel_globalmean

    dics_plots['RadKernel_zonalmean']        = my_plot.plot_RadKernel_zonalmean
    dics_plots['CldRadKernel_zonalmean']     = my_plot.plot_CldRadKernel_zonalmean

    dics_plots['RadKernel_latlon']           = my_plot.plot_RadKernel_latlon
    dics_plots['CldRadKernel_latlon']        = my_plot.plot_CldRadKernel_latlon
    dics_plots['tas_latlon']                 = my_plot.plot_tas_latlon

    dics_plots['LCF']                        = my_plot.plot_LCF 
    dics_plots['zm_CLOUD']                   = my_plot.plot_zm_CLOUD
    dics_plots['latlon_CLOUD']               = my_plot.plot_latlon_CLOUD
    dics_plots['webb_decomp']                = my_plot.plot_webb_decomp
    
    dics_plots['CRE_globalmean_P4KvsFuture'] = my_plot.plot_CRE_globalmean_P4KvsFuture
    dics_plots['RadForcing_globalmean']      = my_plot.plot_RadForcing_globalmean

    return dics_plots


class plots:
    def __init__(self, cases,ref_casesA,Add_otherCMIPs,datadir_v2, datadir_v1, s1, s2, fh, fh1, a1, colors,figdir,ncase, linestyles, linewidths,Add_amipFuture,highlight_CESM2,lw_CESM2,ls_CESM2,lc_CESM2, datadir_Ringer, datadir_RadKernel, datadir_CldRadKernel):
        self.cases = cases
        self.ref_casesA = ref_casesA
        self.Add_otherCMIPs = Add_otherCMIPs
        self.datadir_v2 = datadir_v2
        self.datadir_v1 = datadir_v1
        self.s1 = s1
        self.s2 = s2
        self.fh = fh
        self.fh1 = fh1
        self.a1 = a1
        self.colors = colors
        self.figdir = figdir
        self.ncase = ncase
        self.linestyles = linestyles
        self.linewidths = linewidths 
        self.Add_amipFuture = Add_amipFuture
        self.highlight_CESM2 = highlight_CESM2
        self.lw_CESM2 = lw_CESM2
        self.ls_CESM2 = ls_CESM2
        self.lc_CESM2 = lc_CESM2
        self.datadir_Ringer = datadir_Ringer
        self.datadir_RadKernel = datadir_RadKernel
        self.datadir_CldRadKernel = datadir_CldRadKernel

    ####################################################################################
    ### bar plot for global mean CRE feedback: including E3SMv1 piControl and amip
    ####################################################################################
    def plot_CRE_globalmean(self):
        print('ScatterPlot-CRE-feedback starts ........')
    
        cases_here = copy.deepcopy(self.cases)
        if 'amip-4xCO2' in self.cases:
            cases_here.remove('amip-4xCO2')
    
        df_all = pd.DataFrame()
    
        if self.Add_otherCMIPs:
            # read other CMIP5 and CMIP6 models
            # Oct 19, 2020: reduce model lists to fit both amip-p4K and amipFuture
    
            exp_cntl = [['piControl','amip'],['piControl','amip']]
            exp_new = [['abrupt4xCO2','amip4K'],['abrupt-4xCO2','amip-p4K']]
            
            prefix = 'global_mean_features'
            suffix1 = '*.csv'
            suffix2 = '*.csv'
            
            models_all,cmip5_models,cmip6_models = PDF.get_intersect_withripf(exp_cntl,exp_new,prefix,suffix1,suffix2,self.datadir_Ringer)
    
            exp_cntl = [['piControl','amip'],['piControl','amip']]
            exp_new = [['abrupt4xCO2','amipFuture'],['abrupt-4xCO2','amip-future4K']]
           
            models_all_future,cmip5_models_future,cmip6_models_future = PDF.get_intersect_withripf(exp_cntl,exp_new,prefix,suffix1,suffix2,self.datadir_Ringer)
            
            #print('models_all',models_all,len(models_all))
            #print('cmip5_models', cmip5_models,len(cmip5_models))
            #print('cmip6_models', cmip6_models,len(cmip6_models))
    
            #print('models_all_future',models_all_future,len(models_all_future))
            #print('cmip5_models_future', cmip5_models_future,len(cmip5_models_future))
            #print('cmip6_models_future', cmip6_models_future,len(cmip6_models_future))
    
            # ---- amip4K ---------------------
            df_p4K = pd.DataFrame()
            for model in models_all:
                if model in cmip5_models:
    
                    if model == 'CanESM2_r1i1p1':
                        model_amip = 'CanAM4_r1i1p1'
                    elif model == 'HadGEM2-ES_r1i1p1':
                        model_amip = 'HadGEM2-A_r1i1p1'
                    else:
                	    model_amip = model
    
                    filename = self.datadir_Ringer+'global_mean_features_CMIP5_amip4K_'+model_amip+'.csv'
                else:
                    filename = self.datadir_Ringer+'global_mean_features_CMIP6_amip-p4K_'+model+'.csv'
            
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
    
                    filename = self.datadir_Ringer+'global_mean_features_CMIP5_amipFuture_'+model_amip+'.csv'
                else:
                    filename = self.datadir_Ringer+'global_mean_features_CMIP6_amip-future4K_'+model+'.csv'
            
                df = pd.read_csv(filename,index_col=0)
                df.index = df.loc[:,'varname']
                df2 = df.loc[:,'anomaly_perK']
                df_future[model] = df2
    
        # read amip
        for icase,case in enumerate(cases_here):
            if case == 'v1_coupled':
                # read v1-coupled 
                df_coupled = pd.read_csv(self.datadir_v1+'global_mean_features_CMIP6_abrupt-4xCO2_E3SM-1-0_r1i1p1f1.csv',index_col=0)
                df_coupled.index = df_coupled.loc[:,'varname']
                df_all['v1_coupled'] = df_coupled.loc[:,'anomaly_perK']
            elif case == 'v1_amip4K':
                # read v1-amip
                df_amip = pd.read_csv(self.datadir_v1+'global_mean_features_CMIP6_amip-p4K_E3SM-1-0_r2i1p1f1.csv',index_col=0)
                df_amip.index = df_amip.loc[:,'varname']
                df_all['v1_amip4K'] = df_amip.loc[:,'anomaly_perK']
            elif case == 'v1_future4K':
                # read v1-amip
                df_amip = pd.read_csv(self.datadir_v1+'global_mean_features_CMIP6_amip-future4K_E3SM-1-0_r2i1p1f1.csv',index_col=0)
                df_amip.index = df_amip.loc[:,'varname']
                df_all['v1_future4K'] = df_amip.loc[:,'anomaly_perK']
            elif case == 'amip-4xCO2':   
                continue
            else:
                df1 = pd.read_csv(self.datadir_v2+'global_mean_features_'+case+'.csv',index_col=0)
                df1.index = df1.loc[:,'varname']
        
                df2 = df1.loc[:,'anomaly_perK']
                df_all[case] = df2
        
        # start plotting 
        fig = plt.figure(figsize=(18,12))
        ax = fig.add_subplot(1,1,1)
        
        if 'ts' in df_all.index:
            drop_index = ['ts','SWCLR','LWCLR']
        else:
            drop_index = ['tas','SWCLR','LWCLR']
    
        df_plot = df_all.drop(index=drop_index)
    
        if self.Add_otherCMIPs:
            df_p4K_plot = df_p4K.drop(index=drop_index)
            df_future_plot = df_future.drop(index=drop_index)
       
        x = np.arange(1,len(df_plot.index)+1,1)
        
        for idx,index in enumerate(df_plot.index):
            for icol,column in enumerate(df_plot.columns):
                if column == 'v1_coupled' or column == 'v2_coupled':
                    L1 = ax.scatter(x[idx],df_plot.loc[index,column].tolist(),s=self.s1,alpha=self.a1,label=column,color=self.colors[icol],marker='x')
                elif column == 'v1_amip4K':
                    ax.scatter(x[idx],df_plot.loc[index,column].tolist(),s=self.s1,alpha=self.a1,label=column,color=self.colors[icol],marker='x')
                elif column == 'v1_future4K':
                    ax.scatter(x[idx],df_plot.loc[index,column].tolist(),s=self.s1,alpha=self.a1,label=column,color=self.colors[icol],marker='x')
                else:
                    ax.scatter(x[idx],df_plot.loc[index,column].tolist(),s=self.s1,alpha=self.a1,label=column,color=self.colors[icol])
        
            if self.Add_otherCMIPs:
                # add other CMIP models
                for icol,column in enumerate(df_p4K_plot.columns):
                    if self.highlight_CESM2 and 'CESM2' in column:
                        ax.scatter(x[idx]-0.2,df_p4K_plot.loc[index,column].tolist(),edgecolor='none',facecolor=self.lc_CESM2,alpha=1.0,s=self.s2,marker='x',\
                        label=column.split('_')[0]+'_amip-p4K')
                    else:
                        ax.scatter(x[idx]-0.2,df_p4K_plot.loc[index,column].tolist(),edgecolor='none',facecolor='grey',alpha=self.a1,s=self.s2)
    
                # ensemble mean
                L2 = ax.scatter(x[idx]-0.2,df_p4K_plot.loc[index,:].mean().tolist(),color='black',s=self.s2)
    
                if self.Add_amipFuture:
                    for icol,column in enumerate(df_future_plot.columns):
                        if self.highlight_CESM2 and 'CESM2' in column:
                            ax.scatter(x[idx]+0.2,df_future_plot.loc[index,column].tolist(),edgecolor='none',facecolor=self.lc_CESM2,alpha=1.0,s=self.s2,marker='x',\
                            label=column.split('_')[0]+'_amip-future4K')
                        else:
                            ax.scatter(x[idx]+0.2,df_future_plot.loc[index,column].tolist(),edgecolor='none',facecolor='grey',alpha=self.a1,s=self.s2)
    
                    # ensemble mean
                    L3 = ax.scatter(x[idx]+0.2,df_future_plot.loc[index,:].mean().tolist(),color='red',s=self.s2)
    
                
            ax.tick_params(labelsize=self.fh)
            ax.set_ylabel('W/m$^2$/K',fontsize=self.fh)
            if idx==0:
                if self.Add_amipFuture:
                    if Add_otherCMIPs:
                        legend1 = ax.legend([L2,L3],['amip4K','amipFuture'],fontsize=self.fh1,loc='upper left')
                        ax.legend(fontsize=self.fh1)
                        ax.add_artist(legend1) 
                    else:
                        ax.legend(fontsize=self.fh1)
                else:
                    if self.Add_otherCMIPs:
                        legend1 = ax.legend([L2],['amip4K'],fontsize=self.fh1,loc='upper left')
                        ax.legend(fontsize=self.fh1)
                        ax.add_artist(legend1) 
                    else:
                        ax.legend(fontsize=self.fh1)
        
        plt.xticks(x,df_plot.index)
        ax.set_title('Radiative Feedback',fontsize=self.fh)
        
        ax.set_ylim(-2.5,2.5)
        ax.grid(which='major', linestyle=':', linewidth='1.0', color='grey')
        fig.savefig(self.figdir+'ScatterPlot-CRE-feedback-'+self.cases[-1]+'.png',bbox_inches='tight',dpi=300)
        plt.close(fig)
    
        del(df_all,df_plot)
        print('------------------------------------------------')
        print('ScatterPlot-CRE-feedback is done!')
        print('------------------------------------------------')
    
    ####################################################################
    ### bar plot for radiative feedback based on Radiative kernel 
    ####################################################################
    def plot_RadKernel_globalmean(self):
        print('ScatterPlot-RadKernel-Feedback starts........')
    
        cases_here = copy.deepcopy(self.cases)
        if 'amip-4xCO2' in self.cases:
            cases_here.remove('amip-4xCO2')
    
        df_all = pd.DataFrame()
    
        # read other CMIP5&6 models
        if self.Add_otherCMIPs:
    
            phases = ['CMIP5','CMIP6']
    
            exp_cntl = [['piControl','amip'],['piControl','amip']]
            exp_new = [['abrupt4xCO2','amip4K'],['abrupt-4xCO2','amip-p4K']]
            
            prefix = 'FDBK'
            suffix1 = '*1yr-150yr.csv'
            suffix2 = '*.csv'
            
            models_all,cmip5_models,cmip6_models = PDF.get_intersect(exp_cntl,exp_new,prefix,suffix1,suffix2,self.datadir_RadKernel)
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
    
                    df = pd.read_csv(self.datadir_RadKernel+'FDBK_'+phase+'_'+exp_new[iphase][1]+'_'+model_amip+suffix+'.csv',index_col=0)
                    df2 = df.iloc[:,0]
                    df_others[model] = df2
    
    
        # E3SM
        for icase,case in enumerate(cases_here):
            if case == 'v1_coupled':
                # read v1-coupled 
                df_coupled = pd.read_csv(self.datadir_v1+'FDBK_CMIP6_abrupt-4xCO2_E3SM-1-0_r1i1p1f1_1yr-150yr.csv',index_col=0)
                df_all['v1_coupled'] = df_coupled.iloc[:,0]
            elif case == 'v1_amip4K':
                # read v1-amip
                df_amip = pd.read_csv(self.datadir_v1+'FDBK_CMIP6_amip-p4K_E3SM-1-0_r2i1p1f1_1yr-36yr.csv',index_col=0)
                df_all['v1_amip4K'] = df_amip.iloc[:,0]
            elif case == 'v1_future4K':
                # read v1-amip
                df_amip = pd.read_csv(self.datadir_v1+'FDBK_CMIP6_amip-future4K_E3SM-1-0_r2i1p1f1_1yr-36yr.csv',index_col=0)
                df_all['v1_future4K'] = df_amip.iloc[:,0]
            elif case == 'amip-4xCO2':
                continue
            else:    
                df1 = pd.read_csv(self.datadir_v2+'FDBK_CMIP6_'+case+'.csv',index_col=0)
                if 'a4SST' in case or 'amip-p4K-CESM2' in case:
                    MODEL = 'CESM2'
                else:
                    MODEL = 'E3SM-1-0'
                
                df2 = df1.loc[:,MODEL]
    
                df_all[case] = df2
            
        # start plotting 
        
        fig = plt.figure(figsize=(18,12))
        ax = fig.add_subplot(1,1,1)
        
        drop_index = ['T','dLW_adj','dSW_adj','dnet_adj','LW_resd','SW_resd','net_resd',\
        'T_clr','Planck_clr','LR_clr','WV_clr','ALB_clr','WV_clr_SW','WV_clr_LW','WV_SW','WV_LW',\
        'SWCRE','LWCRE','netCRE','Planck_clr_fxRH','LR_clr_fxRH','RH_clr','LW_clr_sum','SW_clr_sum',\
        'net_clr_sum','LW_clr_dir','SW_clr_dir','net_clr_dir','LW_cld_sum','SW_cld_sum','net_cld_sum',\
        'LW_cld_dir','SW_cld_dir','net_cld_dir','LW_clr_resd','SW_clr_resd','net_clr_resd']
        
    
        df_plot = df_all.drop(index=drop_index)
        x = np.arange(1,len(df_plot.index)+1,1)
    
        if self.Add_otherCMIPs:
            df_others_plot = df_others.drop(index=drop_index)
        
        for idx,index in enumerate(df_plot.index):
            for icol,column in enumerate(df_plot.columns):
                if column == 'v1_coupled' or column == 'v2_coupled':
                    L1 = ax.scatter(x[idx],df_plot.loc[index,column].tolist(),s=self.s1,alpha=self.a1,label=column,color=self.colors[icol],marker='x')
                elif column == 'v1_amip4K':
                    ax.scatter(x[idx],df_plot.loc[index,column].tolist(),s=self.s1,alpha=self.a1,label=column,color=self.colors[icol],marker='x')
                elif column == 'v1_future4K':
                    ax.scatter(x[idx],df_plot.loc[index,column].tolist(),s=self.s1,alpha=self.a1,label=column,color=self.colors[icol],marker='x')
                else:
                    ax.scatter(x[idx],df_plot.loc[index,column].tolist(),s=self.s1,alpha=self.a1,label=column,color=self.colors[icol])
                
            # other CMIP models
            if self.Add_otherCMIPs:
                for icol,column in enumerate(df_others_plot.columns):
                    if self.highlight_CESM2 and 'CESM2' in column:
                        ax.scatter(x[idx], df_others_plot.loc[index,column].tolist(),s=self.s2,edgecolor='none',facecolor=self.lc_CESM2,alpha=1, marker='X',\
                        label = column.split('_')[0]+'_amip-p4K')
                    else:
                        ax.scatter(x[idx]-0.2, df_others_plot.loc[index,column].tolist(),s=self.s2,edgecolor='none',facecolor='grey',alpha=self.a1)
    
                # ensemble mean
                L2 = ax.scatter(x[idx]-0.2, df_others_plot.loc[index,:].mean(),s=self.s2,edgecolor='black',facecolor='black')
    
            ax.tick_params(labelsize=self.fh)
            ax.set_ylabel('W/m$^2$/K',fontsize=self.fh)
            if idx == 0:
                if self.Add_otherCMIPs:
                    legend1 = ax.legend([L2],['amip4K'],fontsize=self.fh1,loc='upper left')
                    ax.legend(fontsize=self.fh1)
                    ax.add_artist(legend1) 
                else:
                    ax.legend(fontsize=self.fh1)
    
        ax.grid(which='major', linestyle=':', linewidth='1.0', color='grey')
        degrees = 0
        plt.xticks(x,df_plot.index,rotation=degrees)
        ax.set_title('Radiative Kernel feedback',fontsize=self.fh)
        
        fig.savefig(self.figdir+'ScatterPlot-RadKernel-Feedback-'+self.cases[-1]+'.png',bbox_inches='tight',dpi=300)
        plt.close(fig)
    
        print('------------------------------------------------')
        print('ScatterPlot-RadKernel-Feedback is done!')
        print('------------------------------------------------')
    
    #########################################################################
    ### zonal mean plot of radiative feedback based on radiative kernel
    #########################################################################
    def plot_RadKernel_zonalmean(self):
        print('plot_RadKernel_zonalmean starts ..........')
    
        cases_here = copy.deepcopy(self.cases)
        if 'amip-4xCO2' in self.cases:
            cases_here.remove('amip-4xCO2')
    
        variables = ['SWCRE_ano_grd_adj','LWCRE_ano_grd_adj','netCRE_ano_grd_adj']
        variables_out = ['SW Cloud Feedback','LW Cloud Feedback','NET Cloud Feedback']
        
        nlat = 73
        nlon = 144
        
        # generate figure based on case categories
        for ii in self.ncase:
            fig = plt.figure(figsize=(18,9))
            num1 = 0
    
            # add other CMIP models
            if self.Add_otherCMIPs:
    
                phases = ['CMIP5','CMIP6']
    
                exp_cntl = [['piControl','amip'],['piControl','amip']]
                exp_new = [['abrupt4xCO2','amip4K'],['abrupt-4xCO2','amip-p4K']]
                
                prefix = 'FDBK'
                suffix1 = '*1yr-150yr.csv'
                suffix2 = '*.csv'
                
                models_all,cmip5_models,cmip6_models = PDF.get_intersect(exp_cntl,exp_new,prefix,suffix1,suffix2,self.datadir_RadKernel)
                #print('models_all',models_all,len(models_all))
                #print('cmip5_models',cmip5_models,len(cmip5_models))
                #print('cmip6_models',cmip6_models,len(cmip6_models))
                models = [cmip5_models, cmip6_models]
                model_list = cmip5_models + cmip6_models
    
            for ivar,svar in enumerate(variables):
    
                # get other CMIP models
                if self.Add_otherCMIPs:
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
    
                            f1 = cdms.open(self.datadir_RadKernel+'lat-lon-gfdbk-'+phase+'-'+exp_new[iphase][1]+'-'+model_amip+suffix+'.nc')
                            tmp1 = f1(svar)
                        
                            lats = tmp1.getLatitude()[:]
                            data1[:,imodel] = MV.average(tmp1,axis=1)
                        if iphase == 0:
                            data_others[:,:len(cmip5_models)] = data1
                        else:
                            data_others[:,len(cmip5_models):] = data1
    
     
                # E3SM
                data_all = np.zeros((nlat,len(cases_here[:ii])))
                for icase,case in enumerate(cases_here[:ii]):
                    if case == 'v1_coupled':
                        f1 = cdms.open(self.datadir_v1+'lat-lon-gfdbk-CMIP6-abrupt-4xCO2-E3SM-1-0_r1i1p1f1_1yr-150yr.nc')
                    elif case == 'v1_amip4K':
                        f1 = cdms.open(self.datadir_v1+'lat-lon-gfdbk-CMIP6-amip-p4K-E3SM-1-0_r2i1p1f1_1yr-36yr.nc')
                    elif case == 'v1_future4K':
                        f1 = cdms.open(self.datadir_v1+'lat-lon-gfdbk-CMIP6-amip-future4K-E3SM-1-0_r2i1p1f1_1yr-36yr.nc')
    
                    elif case == 'amip-4xCO2':
                        continue
                    else:
                        f1 = cdms.open(self.datadir_v2+'lat-lon-gfdbk-CMIP6-'+case+'.nc')
        
                    data = f1(svar)
                    lats = data.getLatitude()[:]
        
                    ######################################################
                    # Compute weights and take weighted average over latitude dimension
                    clat = np.cos(np.deg2rad(lats))
                    clat1 = clat/MV.sum(clat)
                    clat1[0] = 0.
    
                    clats = np.zeros(len(clat1))
        
                    for ilat in range(len(clat1)):
                        clats[ilat] = np.sum(clat1[:ilat+1])
    #                clats[0] = 0.
        
                    #Needs = [-90,-50,-30,-15,0, 15,30,50,90]
                    Needs = [-85, -55, -35, -15, 15, 35, 55, 85]
        
                    N = [i for i in range(len(lats)) if lats[i] in Needs]
                    spec_lats = Needs
                    spec_clats = list(np.array(clats)[N])
                    
    #                print('clat1=',clat1)
    #                print('lats=',lats)
    #                print('clats=',clats)
    #                print('spec_lats=',spec_lats)
    #                print('spec_clats=',spec_clats)
    
                    ######################################################
        
                    # get zonal mean
                    data_all[:,icase] = MV.average(data,axis=1)    
        
                # start plotting ...
                ax = fig.add_subplot(2,2,num1+1)
        
                ax.set_prop_cycle(color=self.colors,lw=self.linewidths,linestyle=self.linestyles)
        
                L1 = ax.plot(clats,data_all,alpha=self.a1)
                
                # highlight CESM2
                if self.highlight_CESM2:
                    data_others = pd.DataFrame(data_others,columns = model_list)
                    for column in data_others.columns:
                        if column == 'CESM2':
                            L3 = ax.plot(clats,data_others.loc[:,column],lw=self.lw_CESM2,ls=self.ls_CESM2,color=self.lc_CESM2)
    
                # plot other CMIP models
                if self.Add_otherCMIPs:
                    L2 = ax.plot(clats,np.average(data_others,axis=1),lw=3,label='ENS-MEAN',color='grey',linestyle='-')
                    ax.fill_between(clats, np.amax(data_others,axis=1),np.amin(data_others,axis=1),alpha=0.2,color='grey')
        
                plt.xticks(spec_clats,spec_lats,fontsize=self.fh)
                ax.set_xlim((0,1))
                plt.yticks(fontsize=self.fh)
        
                ax.set_title(variables_out[ivar],fontsize=self.fh)
                ax.set_ylabel('W/m$^2$/K',fontsize=self.fh)
        
                ax.tick_params(axis='both', which='major', labelsize=self.fh)
    #            ax.set_xlim((-90,90))
        
                ax.set_ylim((-2,2))
                ax.grid(which='major', linestyle=':', linewidth='1.0', color='grey')
        
                ax.axhline(y=0,color='grey',linestyle='--',lw=2)
                if ivar == len(variables)-1:
                    if self.highlight_CESM2:
                        ax.legend(L1+L3,cases_here+['CESM2-amip4K'],fontsize=self.fh1,bbox_to_anchor=(1.04,0), loc='lower left')
                    else:
                        #ax.legend(L1,cases_here,fontsize=self.fh1,bbox_to_anchor=(1.04,0), loc='lower left')
                        ax.legend(L1,cases_here,fontsize=self.fh1,loc='best')
    
                num1 += 1
        
            plt.tight_layout()
            fig.savefig(self.figdir+'Zonal-mean-Cloud-RadKernel-Feedback-'+str(np.round(ii,0))+'-'+self.cases[-1]+'.png',bbox_inches='tight',dpi=300)
            plt.close(fig)
    
    
        print('------------------------------------------------')
        print('plot_RadKernel_zonalmean is done!')
        print('------------------------------------------------')
    
    
    ###################################################################
    ### 4. bar plot of cloud feedback based on cloud radiative kernel
    ###################################################################
    
    def plot_CldRadKernel_globalmean(self):
    
        cases_here = copy.deepcopy(self.cases)
        if 'amip-4xCO2' in self.cases:
            cases_here.remove('amip-4xCO2')
    
        print('ScatterPlot-Cloud-feedback-Decomposition starts .........')
        # other CMIP models
        if self.Add_otherCMIPs:
    
            phases = ['CMIP5','CMIP6']
            exp_cntl = [['piControl','amip'],['piControl','amip']]
            exp_new = [['abrupt4xCO2','amip4K'],['abrupt-4xCO2','amip-p4K']]
    #        exp_cntl = [['amip','amip'],['amip','amip']]
    #        exp_new = [['amip4K','amip4K'],['amip-p4K','amip-p4K']]
           
            prefix = 'decomp_global_mean_lw'
            suffix1 = '*1yr-150yr.csv'
            suffix2 = '*.csv'
            
            models_all,cmip5_models,cmip6_models = PDF.get_intersect(exp_cntl,exp_new,prefix,suffix1,suffix2,self.datadir_CldRadKernel)
            models = [cmip5_models, cmip6_models]
            model_list = cmip5_models + cmip6_models
    
    #        print('models_all=',models_all)
    #        print('cmip5_models=',cmip5_models)
    #        print('cmip6_models=',cmip6_models)
    
            df_LW_others = pd.DataFrame()
            df_SW_others = pd.DataFrame()
    
            for iphase,phase in enumerate(phases):
                for imodel,model in enumerate(models[iphase]):
                    df = pd.read_csv(self.datadir_CldRadKernel+'decomp_global_mean_lw_'+phase+'_'+exp_new[iphase][1]+'_'+model+'.csv',index_col=0)
                    df2 = df.loc[:,model]
                    df_LW_others[model] = df2
    
                    df = pd.read_csv(self.datadir_CldRadKernel+'decomp_global_mean_sw_'+phase+'_'+exp_new[iphase][1]+'_'+model+'.csv',index_col=0)
                    df2 = df.loc[:,model]
                    df_SW_others[model] = df2
    
        # E3SM
        for ii in self.ncase:
    
            df_LW_all = pd.DataFrame()
            df_SW_all = pd.DataFrame()
            
            for icase,case in enumerate(cases_here[:ii]):
                if case == 'v1_coupled':
                    df1 = pd.read_csv(self.datadir_v1+'decomp_global_mean_lw_CMIP6_abrupt-4xCO2_E3SM-1-0_r1i1p1f1_1yr-150yr.csv',index_col=0)
                    df2 = df1.loc[:,'E3SM-1-0']
                    df_LW_all[case] = df2
                    
                    df1 = pd.read_csv(self.datadir_v1+'decomp_global_mean_sw_CMIP6_abrupt-4xCO2_E3SM-1-0_r1i1p1f1_1yr-150yr.csv',index_col=0)
                    df2 = df1.loc[:,'E3SM-1-0']
                    df_SW_all[case] = df2
                elif case == 'v1_amip4K':
                    df1 = pd.read_csv(self.datadir_v1+'decomp_global_mean_lw_CMIP6_amip-p4K_E3SM-1-0_r2i1p1f1_1yr-36yr.csv',index_col=0)
                    df2 = df1.loc[:,'E3SM-1-0']
                    df_LW_all[case] = df2
                    
                    df1 = pd.read_csv(self.datadir_v1+'decomp_global_mean_sw_CMIP6_amip-p4K_E3SM-1-0_r2i1p1f1_1yr-36yr.csv',index_col=0)
                    df2 = df1.loc[:,'E3SM-1-0']
                    df_SW_all[case] = df2
                elif case == 'v1_future4K':
                    df1 = pd.read_csv(self.datadir_v1+'decomp_global_mean_lw_CMIP6_amip-future4K_E3SM-1-0_r2i1p1f1_1yr-36yr.csv',index_col=0)
                    df2 = df1.loc[:,'E3SM-1-0']
                    df_LW_all[case] = df2
                    
                    df1 = pd.read_csv(self.datadir_v1+'decomp_global_mean_sw_CMIP6_amip-future4K_E3SM-1-0_r2i1p1f1_1yr-36yr.csv',index_col=0)
                    df2 = df1.loc[:,'E3SM-1-0']
                    df_SW_all[case] = df2
    
                elif case == 'amip-4xCO2':
                    continue
                else:
                    if 'a4SST' in case or 'amip-p4K-CESM2' in case:
                        MODEL = 'CESM2'
                    else:
                        MODEL = 'E3SM-1-0'
     
                    df1 = pd.read_csv(self.datadir_v2+'decomp_global_mean_lw_'+case+'.csv',index_col=0)
                    df2 = df1.loc[:,MODEL]
                    df_LW_all[case] = df2
            
                    df1 = pd.read_csv(self.datadir_v2+'decomp_global_mean_sw_'+case+'.csv',index_col=0)
                    df2 = df1.loc[:,MODEL]
                    df_SW_all[case] = df2
                
            # get net cloud feedback 
            df_net_all = df_SW_all + df_LW_all
            if self.Add_otherCMIPs:
                df_net_others = df_SW_others + df_LW_others
                
            # ----------------------------------------------------------
            # start plotting 
            fig, axes = plt.subplots(nrows=3,ncols=1,figsize=(12,15))
            
            titles = ['All Cloud CTP bins', "Non-Low Cloud CTP bins", "Low Cloud CTP bins"]
            labels = ['Total','Amount','Altitude','Optical Depth','Residual']
            wc = 1 # circle
            wf = 0 # filled 
            
            w = 0.25
            w2 = 0.08
            w3 = 0.08
        
            x = np.asarray([1,2,3,4,5])
            
            for jj in range(3): ## loop for each panel, jj reprents All, Non-Low and Low Cloud 
                for icol,column in enumerate(df_LW_all.columns):
                    y1 = df_LW_all.iloc[jj*5:(jj+1)*5,icol]
                    if column == 'v1_coupled' or column == 'v2_coupled':
                        axes[jj].scatter(x-w+w2,y1.values.tolist(),marker='x',s=self.s1,color=self.colors[icol],alpha=self.a1,label=column)
                    elif column == 'v1_amip4K':
                        axes[jj].scatter(x-w+w2,y1.values.tolist(),marker='x',s=self.s1,color=self.colors[icol],alpha=self.a1,label=column)
                    elif column == 'v1_future4K':
                        axes[jj].scatter(x-w+w2,y1.values.tolist(),marker='x',s=self.s1,color=self.colors[icol],alpha=self.a1,label=column)
    
                    else:
        #                L1 = axes[jj].scatter(x-w+w2,y1.values.tolist(),marker='o',s=self.s2,color=self.colors[icol],alpha=self.a1,label=column)
                        La = axes[jj].scatter(x-w+w2,y1.values.tolist(),marker='v',s=self.s2,color=self.colors[icol],alpha=self.a1,label=column)
        
                for icol,column in enumerate(df_net_all.columns):
                    y1 = df_net_all.iloc[jj*5:(jj+1)*5,icol]
                    if column == 'v1_coupled' or column == 'v2_coupled':
                        axes[jj].scatter(x+w2,y1.values.tolist(),marker='x',s=self.s1,color=self.colors[icol],alpha=self.a1,label='_nolegend_')
                    elif column == 'v1_amip4K':
                        axes[jj].scatter(x+w2,y1.values.tolist(),marker='x',s=self.s1,color=self.colors[icol],alpha=self.a1,label='_nolegend_')
                    elif column == 'v1_future4K':
                        axes[jj].scatter(x+w2,y1.values.tolist(),marker='x',s=self.s1,color=self.colors[icol],alpha=self.a1,label='_nolegend_')
    
                    else:
        #                L2 = axes[jj].scatter(x+w2,y1.values.tolist(),marker='o',s=self.s2,color=self.colors[icol],alpha=self.a1,label=column)
                        Lb = axes[jj].scatter(x+w2,y1.values.tolist(),marker='o',s=self.s2,color=self.colors[icol],alpha=self.a1,label='_nolegend_')
        
                for icol,column in enumerate(df_SW_all.columns):
                    y1 = df_SW_all.iloc[jj*5:(jj+1)*5,icol]
                    if column == 'v1_coupled' or column == 'v2_coupled':
                        axes[jj].scatter(x+w+w2,y1.values.tolist(),marker='x',s=self.s1,color=self.colors[icol],alpha=self.a1,label='_nolegend_')
                    elif column == 'v1_amip4K':
                        axes[jj].scatter(x+w+w2,y1.values.tolist(),marker='x',s=self.s1,color=self.colors[icol],alpha=self.a1,label='_nolegend_')
                    elif column == 'v1_future4K':
                        axes[jj].scatter(x+w+w2,y1.values.tolist(),marker='x',s=self.s1,color=self.colors[icol],alpha=self.a1,label='_nolegend_')
    
                    else:
         #               L3 = axes[jj].scatter(x+w+w2,y1.values.tolist(),marker='o',s=self.s2,color=self.colors[icol],alpha=self.a1,label=column)
                        Lc = axes[jj].scatter(x+w+w2,y1.values.tolist(),marker='^',s=self.s2,color=self.colors[icol],alpha=self.a1,label='_nolegend_')
           
        
                plt.legend([La,Lb,Lc],["LW","NET","SW"],scatterpoints=1,loc="upper left",fontsize=self.fh1)

                # CMIP - other models
                if self.Add_otherCMIPs:
                    a2 = 0.3
                    s3 = 100
                    for icol,column in enumerate(df_LW_others.columns):
                        y1 = df_LW_others.iloc[jj*5:(jj+1)*5,icol]
                        if self.highlight_CESM2 and 'CESM2' in column:
                            axes[jj].scatter(x-w+w2-w3,y1.values.tolist(),marker='X',s=s3,color='red',alpha=a2,\
                            label=column.split('_')[0]+'_amip-p4K')
                        else:
                            axes[jj].scatter(x-w+w2-w3,y1.values.tolist(),marker='v',s=s3,color='red',alpha=a2,label='_nolegend_')
        
                    for icol,column in enumerate(df_net_others.columns):
                        y1 = df_net_others.iloc[jj*5:(jj+1)*5,icol]
                        if self.highlight_CESM2 and 'CESM2' in column:
                            axes[jj].scatter(x+w2-w3,y1.values.tolist(),marker='X',s=s3,color='grey',alpha=a2,\
                            label=column.split('_')[0]+'_amip-p4K')
                        else:
                            axes[jj].scatter(x+w2-w3,y1.values.tolist(),marker='o',s=s3,color='grey',alpha=a2,label='_nolegend_')
        
                    for icol,column in enumerate(df_SW_others.columns):
                        y1 = df_SW_others.iloc[jj*5:(jj+1)*5,icol]
                        if self.highlight_CESM2 and 'CESM2' in column:
                            axes[jj].scatter(x+w+w2-w3,y1.values.tolist(),marker='X',s=s3,color='blue',alpha=a2,\
                            label=column.split('_')[0]+'_amip-p4K')
                        else:
                            axes[jj].scatter(x+w+w2-w3,y1.values.tolist(),marker='^',s=s3,color='blue',alpha=a2,label='_nolegend_')
        
        
                    L1 = axes[jj].scatter(x-w+w2-w3,df_LW_others.iloc[jj*5:(jj+1)*5,:].mean(axis=1),marker='v',s=s3,color='red',alpha=1.0,label='_nolegend_')
                    L2 = axes[jj].scatter(x+w2-w3,df_net_others.iloc[jj*5:(jj+1)*5,:].mean(axis=1),marker='o',s=s3,color='grey',alpha=1.0,label='_nolegend_')
                    L3 = axes[jj].scatter(x+w+w2-w3,df_SW_others.iloc[jj*5:(jj+1)*5,:].mean(axis=1),marker='^',s=s3,color='blue',alpha=1.0,label='_nolegend_')
        
                    #plt.legend((L1,L2,L3),["LW","NET","SW"],scatterpoints=1,bbox_to_anchor=(1,1),loc="best",fontsize=self.fh1)
        
                if jj == 0:
                    axes[jj].legend(fontsize=self.fh1,ncol=2)
         
                axes[jj].grid(which='major', linestyle=':', linewidth='1.0', color='grey')
                
                axes[jj].axhline(0, color="grey", linestyle="-",linewidth=2)
                axes[jj].set_ylabel('W/$m^2$/K',fontsize=self.fh)
                axes[jj].set_title(titles[jj],fontsize=self.fh)
                axes[jj].set_ylim([-0.5,1.5])
                
                major_ticks = np.arange(-1.0,1.5,0.5)
                minor_ticks = np.arange(-1.0,1.5,0.25)
                axes[jj].set_yticks(major_ticks)
                axes[jj].set_yticks(minor_ticks,minor=True)
                
                axes[jj].tick_params(axis='both', which='major', labelsize=self.fh)
                
                axes[jj].grid(axis='y',c='grey',linestyle='-.',which = 'major')
                
                if jj == 2:
                    plt.xticks(x,['Total','Amount','Altitude','Optical Depth','Residual'])
                else:
                    axes[jj].set_xticklabels("")
                    
            plt.tight_layout()
            fig.savefig(self.figdir+'ScatterPlot-Cloud-feedback-Decomposition-'+str(np.round(ii,0))+'-'+self.cases[-1]+'.png',bbox_inches='tight',dpi=300)
            plt.close(fig)
    
        print('------------------------------------------------')
        print('ScatterPlot-Cloud-feedback-Decomposition is done!')
        print('------------------------------------------------')
    
    #####################################################################
    ### 5. zonal mean cloud feedback based on cloud radiative kernel
    #####################################################################
    def plot_CldRadKernel_zonalmean(self):
    
        print('ZonalMean-Cloud-feedback-Decomposition starts ........')
    
        cases_here = copy.deepcopy(self.cases)
        if 'amip-4xCO2' in self.cases:
            cases_here.remove('amip-4xCO2')
    
    
        nlat = 90
        nlon = 144
    
        if self.Add_otherCMIPs:
            phases = ['CMIP5','CMIP6']
            exp_cntl = [['piControl','amip'],['piControl','amip']]
            exp_new = [['abrupt4xCO2','amip4K'],['abrupt-4xCO2','amip-p4K']]
            
            prefix = 'decomp_global_mean_lw'
            suffix1 = '*1yr-150yr.csv'
            suffix2 = '*.csv'
            
            models_all,cmip5_models,cmip6_models = PDF.get_intersect(exp_cntl,exp_new,prefix,suffix1,suffix2,self.datadir_CldRadKernel)
            #print('models_all',models_all,len(models_all))
            #print('cmip5_models',cmip5_models,len(cmip5_models))
            #print('cmip6_models',cmip6_models,len(cmip6_models))
            models = [cmip5_models, cmip6_models]
            model_list = cmip5_models + cmip6_models
     
        
        # --- define variables 
        levs = ['ALL','HI680','LO680',]
        components = ['NET','SW','LW']
        decomps = ['tot','amt','alt','tau']
    
        levs_out = ['ALL', 'High cloud', 'Low cloud']
        components_out = ['NET','SW','LW']
        decomps_out = ['Total','Amount','Altitude','Optical Depth']
        
        # generate figure based on case categories
    
        for ii in self.ncase:
            #<qinyi 2021-05-19 #------------------
            fig1 = plt.figure(figsize=(24,12))
            num2 = 0 
    
            for icomp,component in enumerate(components):
                for ilev,lev in enumerate(levs):
                    fig = plt.figure(figsize=(18,9))
                    num1 = 0
        
                    for idecomp,decomp in enumerate(decomps):
                        varname = lev+'_'+component+'cld_'+decomp
                        varname_out = components_out[icomp]+' '+levs_out[ilev]+' '+decomps_out[idecomp]+' Feedback'
    
                        if component == 'NET':
                            varSW = lev+'_SWcld_'+decomp
                            varLW = lev+'_LWcld_'+decomp
        
                        # add other CMIP models
                        if self.Add_otherCMIPs:
                
                            data_others = np.zeros((nlat,len(model_list)))
                            avgdata_others = np.zeros(len(model_list))
    
                            for iphase,phase in enumerate(phases):
                                data1 = np.zeros((nlat,len(models[iphase])))
                                avgdata1 = np.zeros((len(models[iphase])))
                                for imodel,model in enumerate(models[iphase]):
                                    f1 = cdms.open(self.datadir_CldRadKernel+'global_cloud_feedback_'+phase+'_'+exp_new[iphase][1]+'_'+model+'.nc')
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
                        data_all = np.zeros((nlat,len(cases_here[:ii])))
                        avgdata = np.zeros(len(cases_here[:ii]))
        
                        for icase,case in enumerate(cases_here[:ii]):
                            if case == 'v1_coupled':
                                f1 = cdms.open(self.datadir_v1+'global_cloud_feedback_CMIP6_abrupt-4xCO2_E3SM-1-0_r1i1p1f1_1yr-150yr.nc')
                            elif case == 'v1_amip4K':
                                f1 = cdms.open(self.datadir_v1+'global_cloud_feedback_CMIP6_amip-p4K_E3SM-1-0_r2i1p1f1_1yr-36yr.nc')
                            elif case == 'v1_future4K':
                                f1 = cdms.open(self.datadir_v1+'global_cloud_feedback_CMIP6_amip-future4K_E3SM-1-0_r2i1p1f1_1yr-36yr.nc')
                            elif case == 'amip-4xCO2':
                                continue
                            else:
                                f1 = cdms.open(self.datadir_v2+'global_cloud_feedback_'+case+'.nc')
        
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
                                clats[ilat] = np.sum(clat1[:ilat+1])
    #                        clats[0] = 0.
        
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
                        ax.set_prop_cycle(color=self.colors,lw=self.linewidths,linestyle=self.linestyles)
    
                        L1 = ax.plot(clats,data_all)
    
                        # plot other CMIP models
                        if self.Add_otherCMIPs:
                            L2 = ax.plot(clats,np.average(data_others,axis=1),lw=3,label='ENS-MEAN',color='grey',linestyle='-')
                            ax.fill_between(clats, np.amax(data_others,axis=1),np.amin(data_others,axis=1),alpha=0.2,color='grey')
    
                        plt.xticks(spec_clats,spec_lats,fontsize=self.fh)
                        ax.set_xlim((0,1))
    
                        plt.yticks(fontsize=self.fh)
                        ax.set_title(varname_out,fontsize=self.fh)
                        ax.set_ylabel('W/m$^2$/K',fontsize=self.fh)
    
                        ax.set_ylim((-2,2))
                        ax.grid(which='major', linestyle=':', linewidth='1.0', color='grey')
                        ax.axhline(y=0,color='grey',linestyle='--',lw=2)
                        ax.tick_params(axis='both', which='major', labelsize=self.fh)
        
                        if 'tau' in varname:
                            #ax.legend(L1,cases_here,fontsize=self.fh1,bbox_to_anchor=(1.04,0), loc='lower left')
                            ax.legend(L1,cases_here,fontsize=self.fh1,loc='best')
    
                        num1 += 1
    
                        #<qinyi 2021-05-19 #------------------
                        if varname in ['LO680_SWcld_amt','LO680_SWcld_tau','HI680_SWcld_amt','HI680_SWcld_tau','HI680_LWcld_amt','HI680_LWcld_alt','HI680_LWcld_tau']:
                            ax1 = fig1.add_subplot(2,4,num2+1)
                            ax1.set_prop_cycle(color=self.colors,lw=self.linewidths,linestyle=self.linestyles)
    
                            L3 = ax1.plot(clats,data_all)
    
                            ax1.set_xticks(spec_clats)
                            ax1.set_xticklabels(spec_lats,fontsize=self.fh)
                            ax1.set_xlim((0,1))
        
                            plt.yticks(fontsize=self.fh)
                            ax1.set_title(varname_out,fontsize=self.fh)
                            ax1.set_ylabel('W/m$^2$/K',fontsize=self.fh)
        
                            ax1.set_ylim((-2,2))
                            #ax1.grid(which='major', linestyle=':', linewidth='1.0', color='grey')
                            #ax1.grid(axis='y')
                            ax1.axhline(y=0,color='grey',linestyle='--',lw=2)
                            ax1.tick_params(axis='both', which='major', labelsize=self.fh)
    
                            if varname == 'HI680_SWcld_amt':
                                ax1.legend(L3,cases_here,fontsize=self.fh1,loc='best')
    
                            xpos = 0.95
                            ypos = 0.05
                            for jj in range(len(avgdata)):
                                ax1.text(xpos,ypos+jj*0.05,cases_here[jj]+': '+'{:.2f}'.format(np.round(avgdata[jj],2)),transform=ax1.transAxes,\
                                color=self.colors[jj],fontsize=self.fh1+5,ha='right',va='center')
    
                            num2 += 1
    
                        #<qinyi 2021-05-19 #------------------
    
                    plt.tight_layout()
                    fig.savefig(self.figdir+'ZonalMean-Cloud-feedback-Decomposition-'+lev+'-'+component+'-'+str(np.round(ii,0))+'-'+self.cases[-1]+'.png',bbox_inches='tight',dpi=300)
                    plt.close(fig)
    
            #<qinyi 2021-05-19 #------------------
            plt.tight_layout()
            fig1.savefig(self.figdir+'ZonalMean-Cloud-feedback-Decomposition-KeyComponents-'+str(np.round(ii,0))+'-'+self.cases[-1]+'.png',bbox_inches='tight',dpi=300)
            plt.close(fig)
            #>qinyi 2021-05-19 #------------------
    
    
        print('------------------------------------------------')
        print('ZonalMean-Cloud-feedback-Decomposition is done!')
        print('------------------------------------------------')
    
    
    #####################################################################
    ### 7. LAT-LON cloud feedback difference based on cloud radiative kernel
    #####################################################################
    def plot_CldRadKernel_latlon(self):
    
        print('LatLon-Cloud-feedback-Decomposition starts ........')
    
        cases_here = copy.deepcopy(self.cases)
        if 'amip-4xCO2' in self.cases:
            cases_here.remove('amip-4xCO2')
    
    
        nlat = 90
        nlon = 144
        
        # --- define variables 
        sections = ['ALL','HI680','LO680']
        components = ['NET','SW','LW']
        decomps = ['tot','amt','alt','tau']
        
        # generate figure based on case categories
        for icomp,component in enumerate(components):
            for sec in sections:
    
                # E3SM 
                data_all = np.zeros((nlat,nlon,len(decomps),len(cases_here)))
                avgdata = np.zeros((len(decomps),len(cases_here)))
                data_all = cdms.asVariable(data_all)
    
                for idecomp,decomp in enumerate(decomps):
                    varname = sec+'_'+component+'cld_'+decomp
                    if component == 'NET':
                        varSW = sec+'_SWcld_'+decomp
                        varLW = sec+'_LWcld_'+decomp
    
                    for icase,case in enumerate(cases_here):
                        if case == 'v1_coupled':
                            f1 = cdms.open(self.datadir_v1+'global_cloud_feedback_CMIP6_abrupt-4xCO2_E3SM-1-0_r1i1p1f1_1yr-150yr.nc')
                        elif case == 'v1_amip4K':
                            f1 = cdms.open(self.datadir_v1+'global_cloud_feedback_CMIP6_amip-p4K_E3SM-1-0_r2i1p1f1_1yr-36yr.nc')
                        elif case == 'v1_future4K':
                            f1 = cdms.open(self.datadir_v1+'global_cloud_feedback_CMIP6_amip-future4K_E3SM-1-0_r2i1p1f1_1yr-36yr.nc')
    
                        elif case == 'amip-4xCO2':
                            continue
                        else:
                            f1 = cdms.open(self.datadir_v2+'global_cloud_feedback_'+case+'.nc')
    
                        if component in ['LW','SW']:
                            data = f1(varname)
                        else:
    
                            dataSW = f1(varSW)
                            dataLW = f1(varLW)
                            data = dataSW + dataLW
                            data.setAxisList(dataSW.getAxisList())
    
                        lats = data.getLatitude()
                        lons = data.getLongitude()
    
                        ######################################################
                        data_all[:,:,idecomp,icase] = data
                        # get global mean
                        avgdata[idecomp,icase] = cdutil.averager(data,axis='xy',weights='weighted')
    
                data_all.setAxis(0,lats)
                data_all.setAxis(1,lons)
    
                print('data_all.shape=',data_all.shape)
                
                # -----------------------------------------------------------------
                # start plotting ...
                # -----------------------------------------------------------------
    
                for icase,case in enumerate(cases_here):
                    ref_cases = self.ref_casesA[icase]
                    if len(ref_cases) == 0:
                        continue
    
                    for ref_case in ref_cases:
                        iref = cases_here.index(ref_case)
                        print('iref = ',iref, ref_case)
    
                        fig=plt.figure(figsize=(18,12)) # this creates and increases the figure size
                        plt.suptitle(sec+' CTP bins ['+cases_here[icase]+' minus '+cases_here[iref]+']',fontsize=self.fh,y=0.95)
                        bounds = np.arange(-3,3.25,0.25)
                        cmap = plt.cm.RdBu_r
                        bounds2 = np.append(np.append(-500,bounds),500) # This is only needed for norm if colorbar is extended
                        norm = mpl.colors.BoundaryNorm(bounds2, cmap.N) # make sure the colors vary linearly even if the desired color boundaries are at varied intervals
                        names = [component+'cld_tot',component+'cld_amt',component+'cld_alt',component+'cld_tau']
        
                        for n,name in enumerate(names):
                            # difference b/t icase and v1-coupled
                            DATA = data_all[:,:,n,icase] - data_all[:,:,n,iref]
       
                            DATA.setAxis(0,lats)
                            DATA.setAxis(1,lons)
        
                            ax1 = fig.add_subplot(3,2,n+1,projection=ccrs.Robinson(central_longitude=180.))
                            im1 = ax1.contourf(lons[:],lats[:],DATA,bounds,transform=ccrs.PlateCarree(),cmap=cmap,norm=norm,extend='both')
                            ax1.coastlines()
                            ax1.set_global()
        
                            avgDATA = avgdata[n,icase] - avgdata[n,iref]

                            # get spatial correlation, NRMSE and RMSE
                            wts = np.cos(np.deg2rad(lats[:]))
                            daa = data_all[:,:,n,icase]
                            dbb = data_all[:,:,n,iref]
                            cor,NRMSE, RMSE = PDF.pattern_cor(daa,dbb,wts,1)
                            print('cor=',cor, 'NRMSE=',NRMSE, 'RMSE=',RMSE)
 
                            plt.title(name+' ['+str(np.round(avgDATA,3))+']\nNRMSE='+str(np.round(NRMSE,2))+', COR='+str(np.round(cor,2)),fontsize=self.fh)

                            cb = plt.colorbar(im1,orientation='vertical',drawedges=True,ticks=bounds[::2])
                            cb.set_label('W/m$^2$/K')
    
        
                        fig.subplots_adjust(top=0.9)
    
                        fig.savefig(self.figdir+'LatLon-Cloud-feedback-Decomposition-'+cases_here[icase]+'.minus.'+cases_here[iref]+'-'+component+'-'+sec+'.png',dpi=300,bbox_inches='tight')
                        plt.close(fig)
                
    
    #####################################################################
    ### 8. LAT-LON tas anomaly based on radiative kernel output
    #####################################################################
    def plot_tas_latlon(self):
        print('plot_tas_latlon starts ..........')
    
        cases_here = copy.deepcopy(self.cases)
        if 'amip-4xCO2' in self.cases:
            cases_here.remove('amip-4xCO2')
    
        variables = ['tas_ano_grd']
        variables_out = ['TAS']
        
        nlat = 73
        nlon = 144
        
        # add other CMIP models
        if self.Add_otherCMIPs:
    
            phases = ['CMIP5','CMIP6']
    
            exp_cntl = [['piControl','amip'],['piControl','amip']]
            exp_new = [['abrupt4xCO2','amip4K'],['abrupt-4xCO2','amip-p4K']]
            
            prefix = 'FDBK'
            suffix1 = '*1yr-150yr.csv'
            suffix2 = '*.csv'
            
            models_all,cmip5_models,cmip6_models = PDF.get_intersect(exp_cntl,exp_new,prefix,suffix1,suffix2,self.datadir_RadKernel)
            models = [cmip5_models, cmip6_models]
            model_list = cmip5_models + cmip6_models
    
        for ivar,svar in enumerate(variables):
    
            # get other CMIP models
            if self.Add_otherCMIPs:
                data_others = np.zeros((nlat,nlon,len(model_list)))
                data_others_zm = np.zeros((nlat,len(model_list)))
    
                for iphase,phase in enumerate(phases):
                    if phase == 'CMIP5':
                        suffix = '_Latest-Oct18_1yr-27yr'
                    else:
                        suffix = '_Latest-Oct18_1yr-36yr'
                    
                    data1 = np.zeros((nlat,nlon,len(models[iphase])))
                    data1_zm = np.zeros((nlat,len(models[iphase])))
    
                    for imodel,model in enumerate(models[iphase]):
                        if model == 'CanESM2':
                            model_amip = 'CanAM4'
                        elif model == 'HadGEM2-ES':
                            model_amip = 'HadGEM2-A'
                        else:
                    	    model_amip = model
    
                        f1 = cdms.open(self.datadir_RadKernel+'lat-lon-gfdbk-'+phase+'-'+exp_new[iphase][1]+'-'+model_amip+suffix+'.nc')
                        tmp1 = f1(svar)
                    
                        lats = tmp1.getLatitude()[:]
                        lons = tmp1.getLongitude()[:]
                        data1[:,:,imodel] = tmp1
                        data1_zm[:,imodel] = MV.average(tmp1,axis=1)
    
                    if iphase == 0:
                        data_others[:,:,:len(cmip5_models)] = data1
                        data_others_zm[:,:len(cmip5_models)] = data1_zm
                    else:
                        data_others[:,:,len(cmip5_models):] = data1
                        data_others_zm[:,len(cmip5_models):] = data1_zm
    
            # E3SM
            data_all = np.zeros((nlat,nlon,len(cases_here)))
            data_all_zm = np.zeros((nlat,len(cases_here)))
    
            for icase,case in enumerate(cases_here):
                if case == 'v1_coupled':
                    f1 = cdms.open(self.datadir_v1+'lat-lon-gfdbk-CMIP6-abrupt-4xCO2-E3SM-1-0_r1i1p1f1_1yr-150yr.nc')
                elif case == 'v1_amip4K':
                    f1 = cdms.open(self.datadir_v1+'lat-lon-gfdbk-CMIP6-amip-p4K-E3SM-1-0_r2i1p1f1_1yr-36yr.nc')
                elif case == 'v1_future4K':
                    f1 = cdms.open(self.datadir_v1+'lat-lon-gfdbk-CMIP6-amip-future4K-E3SM-1-0_r2i1p1f1_1yr-36yr.nc')
                elif case == 'amip-4xCO2':
                    continue
                else:
                    f1 = cdms.open(self.datadir_v2+'lat-lon-gfdbk-CMIP6-'+case+'.nc')
    
                data = f1(svar)
                lats = data.getLatitude()
                lons = data.getLongitude()
    
                # get zonal mean
                data_all[:,:,icase] = data
                data_all_zm[:,icase] = MV.average(data,axis=1)    
    
            #----------------------------------------------------------
            # start plotting ...
            #----------------------------------------------------------
    
            for icase,case in enumerate(cases_here):
                ref_cases = self.ref_casesA[icase]
                if len(ref_cases) == 0:
                    continue
    
                for ref_case in ref_cases:
                    iref = cases_here.index(ref_case)
    
                    fig=plt.figure(figsize=(18,12)) # this creates and increases the figure size
                    nrow = 3; ncol = 1
    
                    cmap0 = plt.cm.Reds
                    bounds0 = np.arange(0,1.7,0.1)
                    bounds01 = np.append(np.append(-500,bounds0),500) # This is only needed for norm if colorbar is extended
                    norm0 = mpl.colors.BoundaryNorm(bounds01, cmap0.N) # make sure the colors vary linearly even if the desired color boundaries are at varied intervals
    
                    cmap1 = plt.cm.RdBu_r
                    bounds1 = np.arange(-0.8,0.85,0.1)
                    bounds11 = np.append(np.append(-500,bounds1),500) # This is only needed for norm if colorbar is extended
                    norm1 = mpl.colors.BoundaryNorm(bounds11, cmap1.N) # make sure the colors vary linearly even if the desired color boundaries are at varied intervals
    
                    # difference b/t icase and reference case
                    DATAd = data_all[:,:,icase] - data_all[:,:,iref]
                    DATAd = cdms.asVariable(DATAd)
    
                    DATAd.setAxis(0,lats)
                    DATAd.setAxis(1,lons)
    
                    DATA_all = [data_all[:,:,icase], data_all[:,:,iref], DATAd]
                    titles = [cases_here[icase], cases_here[iref], cases_here[icase]+' minus '+cases_here[iref]]
                    bounds = [bounds0, bounds0, bounds1]
                    norms = [norm0, norm0, norm1]
                    cmaps = [cmap0, cmap0, cmap1]
    
                    for iDATA,DATA in enumerate(DATA_all):
    
                        ax1 = fig.add_subplot(nrow,ncol,iDATA+1,projection=ccrs.Robinson(central_longitude=180.))
                        im1 = ax1.contourf(lons[:],lats[:],DATA,bounds[iDATA],transform=ccrs.PlateCarree(),cmap=cmaps[iDATA],norm=norms[iDATA],extend='both')
                        ax1.coastlines()
                        ax1.set_global()
        
                        plt.title(titles[iDATA],fontsize=17)
    
                        PDF.make_colorbar(ax1, 'K/K', 13, im1, nbins = 9, orientation='vertical')
    
    #                cb = plt.colorbar(im1,orientation='vertical',drawedges=True,ticks=bounds)
    #                cb.set_label('K/K')
     
                    # add_common_colorbar(fig,im,axes,units,orientation='vertical',nbins=9,fontsize=15)
                    # PDF.add_common_colorbar(fig,im1,ax1,'K/K', orientation='vertical', nbins=9, fontsize=15)
    
                    fig.subplots_adjust(top=0.9)
    
                    fig.savefig(self.figdir+'LatLon-TAS-'+cases_here[icase]+'.minus.'+cases_here[iref]+'.png',bbox_inches='tight',dpi=300)
                    plt.close(fig)
                
    
        print('------------------------------------------------')
        print('plot_tas_latlon is done!')
        print('------------------------------------------------')
    
    
    #####################################################################
    ### 9. LAT-LON RadKernel
    #####################################################################
    def plot_RadKernel_latlon(self):
        print('plot_RadKernel_latlon starts ..........')
        
        cases_here = copy.deepcopy(self.cases)
        if 'amip-4xCO2' in self.cases:
            cases_here.remove('amip-4xCO2')
    
        variables = ['SWCRE_ano_grd_adj','LWCRE_ano_grd_adj','netCRE_ano_grd_adj']
        variables_out = ['SWCRE','LWCRE','netCRE']
        
        nlat = 73
        nlon = 144
        
        # E3SM
        data_all = np.zeros((nlat,nlon,len(variables),len(cases_here)))
        data_all_zm = np.zeros((nlat,len(variables),len(cases_here)))
    
        for ivar,svar in enumerate(variables):
            for icase,case in enumerate(cases_here):
                if case == 'v1_coupled':
                    f1 = cdms.open(self.datadir_v1+'lat-lon-gfdbk-CMIP6-abrupt-4xCO2-E3SM-1-0_r1i1p1f1_1yr-150yr.nc')
                elif case == 'v1_amip4K':
                    f1 = cdms.open(self.datadir_v1+'lat-lon-gfdbk-CMIP6-amip-p4K-E3SM-1-0_r2i1p1f1_1yr-36yr.nc')
                elif case == 'v1_future4K':
                    f1 = cdms.open(self.datadir_v1+'lat-lon-gfdbk-CMIP6-amip-future4K-E3SM-1-0_r2i1p1f1_1yr-36yr.nc')
                elif case == 'amip-4xCO2':
                    continue
                else:
                    f1 = cdms.open(self.datadir_v2+'lat-lon-gfdbk-CMIP6-'+case+'.nc')
    
                data = f1(svar)
                lats = data.getLatitude()
                lons = data.getLongitude()
    
                # get zonal mean
                data_all[:,:,ivar,icase] = data
                data_all_zm[:,ivar,icase] = MV.average(data,axis=1)
        
        #=============================================================
        # start plotting ...
        #=============================================================
        for icase,case in enumerate(cases_here):
            ref_cases = self.ref_casesA[icase]
            if len(ref_cases) == 0:
                continue
    
            for ref_case in ref_cases:
                iref = cases_here.index(ref_case)
    
                #----------------------------------------------------------
                # define figures                         
                #----------------------------------------------------------
                nrow = 3; ncol = 3
                fig = plt.figure(figsize=(ncol*8,nrow*4)) # this creates and increases the figure size
    
                bounds = np.arange(-3,3.25,0.25)

                cmap = plt.cm.RdBu_r
                bounds2 = np.append(np.append(-500,bounds),500) # This is only needed for norm if colorbar is extended
                norm = mpl.colors.BoundaryNorm(bounds2, cmap.N) # make sure the colors vary linearly even if the desired color boundaries are at varied intervals

                case_out = cases_here[icase]
                ref_case_out = cases_here[iref]
    
                for ivar,svar in enumerate(variables):

                    data_plot = [data_all[:,:,ivar,icase], data_all[:,:,ivar,iref], data_all[:,:,ivar,icase] - data_all[:,:,ivar,iref] ]
                    title_plot = [case_out, ref_case_out, case_out+' minus '+ref_case_out]

                    num0 = 0
                    for DATA in data_plot:

                        print('minmax DATA=',genutil.minmax(DATA))
                        DATA = cdms.asVariable(DATA)
                        DATA.setAxis(0,lats)
                        DATA.setAxis(1,lons)
    
                        ax1 = fig.add_subplot(nrow,ncol,ivar*ncol+num0+1,projection=ccrs.Robinson(central_longitude=180.))
                        im1 = ax1.contourf(lons[:],lats[:],DATA,bounds,transform=ccrs.PlateCarree(),cmap=cmap,norm=norm,extend='both')
                        ax1.coastlines()
                        ax1.set_global()
    
                        # global-mean 
                        avgDATA = cdutil.averager(DATA,axis='xy',weights='weighted')

                        if num0 == 2: # only for different plot
                            # get spatial correlation, NRMSE and RMSE
                            wts = np.cos(np.deg2rad(lats[:]))
                            daa = data_all[:,:,ivar,icase]
                            dbb = data_all[:,:,ivar,iref]
                            cor,NRMSE, RMSE = PDF.pattern_cor(daa,dbb,wts,1)
                            print('cor=',cor, 'NRMSE=',NRMSE, 'RMSE=',RMSE)
    
                            ax1.set_title(title_plot[num0]+'\n'+variables_out[ivar]+' ['+str(np.round(avgDATA,2))+']\nNRMSE='+str(np.round(NRMSE,2))+', COR='+str(np.round(cor,2)),fontsize=self.fh)

                            PDF.make_colorbar(ax1, 'W/m$^2$/K', self.fh-5, im1, orientation='vertical')

                        else:
                            ax1.set_title(title_plot[num0]+'\n'+variables_out[ivar]+' ['+str(np.round(avgDATA,2))+']',fontsize=self.fh)

                        num0 += 1
    
                fig.subplots_adjust(top=0.9)
    
                fig.savefig(self.figdir+'LatLon-adjusted-CRE-'+case_out+'.minus.'+ref_case_out+'.png',bbox_inches='tight',dpi=300)
                plt.close(fig)
                
        print('------------------------------------------------')
        print('plot_RadKernel_latlon is done!')
        print('------------------------------------------------')
    
    
    #####################################################################
    ### LCF profile over 30S-80S
    #####################################################################
    def plot_LCF(self):
        print('plot_LCF starts ..........')
    
        cases_here = copy.deepcopy(self.cases)
        if 'amip-4xCO2' in self.cases:
            cases_here.remove('amip-4xCO2')
    
        # ===============================================================        
        # ===============================================================        
        nbin = 20
        # E3SM
        data_all = np.zeros((nbin,len(cases_here)))
    
        df_all = pd.DataFrame()
        for icase,case in enumerate(cases_here):
            if case in ['v1_coupled','v1_amip4K','v1_future4K','amip-4xCO2']:
                continue
            else:
                df1 = pd.read_csv(self.datadir_v2+'LCF_binned_by_temperature_'+case+'_20bins-OnlyOcean.csv',index_col=0)
                df1.index = df1.loc[:,'temp']
    
                df_all[case] = df1.loc[:,'bin']
    
        print(df_all)
    
        #----------------------------------------------------------
        # start plotting ...
        #----------------------------------------------------------
        fig = plt.figure(figsize=(12,9))
        ax = fig.add_subplot(1,1,1)
    
        for col in df_all.columns:
            ax.plot(df_all.index, df_all.loc[:,col], label=col)
    
        ax.set_xlabel('Temperature [K]')
        ax.set_ylabel('Liquid Condensate Fraction') 
    
        ax.legend(fontsize=self.fh-3)
    
        fig.savefig(self.figdir+'LCF-vs-temperature-'+cases_here[-1]+'.png',bbox_inches='tight',dpi=300)
        plt.close(fig)
                
        print('------------------------------------------------')
        print('plot_LCF is done!')
        print('------------------------------------------------')
    
    
    #####################################################################
    ### plot_zm_CLOUD
    #####################################################################
    def plot_zm_CLOUD(self):
        print('plot_zm_CLOUD starts ..........')
    
        cases_here = copy.deepcopy(self.cases)
        if 'amip-4xCO2' in self.cases:
            cases_here.remove('amip-4xCO2')
    
        var1 = ['CLOUD','CLDLIQ','CLDICE']
        var1_tmp = ['cl', 'clw', 'cli'] # for v1_coupled
        var1_range = [[0,1,0.1],[0,50,5],[0,10,1]]
        var1_range_d = [[-0.02,0.02,0.004], [-0.8,0.8,0.1], [-0.8,0.8,0.1]]
    
        var1_units = ['fraction','mg/kg','mg/kg']
    
        # ===============================================================        
        # ===============================================================        
        # E3SM
        nlat = 73
        nlon = 144
        nlev = 72
     
        for icase,case in enumerate(cases_here):
            ref_cases = self.ref_casesA[icase]
            if len(ref_cases) == 0:
                continue
            else:
                for ref_case in ref_cases:
                    for ivar,svar in enumerate(var1):
    
                        # -------------------- generate figures -------------------
                        nrow = 3
                        ncol = 2
                        fig = plt.figure(figsize=(ncol*5,nrow*3))
                        # -------------------- generate figures -------------------
    
                        f1 = cdms.open(self.datadir_v2+'global_cloud_'+case+'.nc')
                        if case == 'v1_coupled':
                            svar_in = var1_tmp[ivar]
                        else:
                            svar_in = var1[ivar]
    
                        tmp1 = f1(svar_in+'_pi_clim')
                        tmp2 = f1(svar_in+'_ano_clim')
    
                        levs = tmp1.getLevel()[:]
                        lats = tmp1.getLatitude()[:]
    
                        data_all_pi = MV.average(tmp1,axis=2)
                        data_all_ano =  MV.average(tmp2,axis=2)
    
                        # convert unit from kg/kg to mg/kg
                        if svar in ['CLDLIQ','CLDICE']:
                            data_all_pi = data_all_pi * 1e6
                            data_all_ano = data_all_ano * 1e6
    
                        if svar in ['CLOUD'] and case == 'v1_coupled':
                            data_all_pi = data_all_pi/100.
                            data_all_ano = data_all_ano/100.
    
                        print('data_all_pi.shape=',data_all_pi.shape)
                        print('data_all_ano.shape=',data_all_ano.shape)
        
                        spec_lats = [-85, -55, -35, -15, 15, 35, 55, 85]
                        spec_lats_str = ['85S','55S','35S','15S','15N','35N','55N','85N']
                        clats,spec_clats = PDF.get_scaled_lat(lats, spec_lats)
    
                        # -------------- read reference P4K data -----------------
                        fref = cdms.open(self.datadir_v2+'global_cloud_'+ref_case+'.nc')
                        if ref_case == 'v1_coupled':
                            svar_in = var1_tmp[ivar]
                        else:
                            svar_in = var1[ivar]
    
                        data_all_pi_ref = MV.average(fref(svar_in+'_pi_clim'),axis=2)
                        data_all_ano_ref = MV.average(fref(svar_in+'_ano_clim'),axis=2)
    
                        if svar in ['CLDLIQ','CLDICE']:
                            data_all_pi_ref = data_all_pi_ref * 1e6
                            data_all_ano_ref = data_all_ano_ref * 1e6
    
                        if svar in ['CLOUD'] and ref_case == 'v1_coupled':
                            data_all_pi_ref = data_all_pi_ref/100.
                            data_all_ano_ref = data_all_ano_ref/100.
    
                        #----------------------------------------------------------
                        # start plotting ...
                        #----------------------------------------------------------
                        data_plot = [data_all_pi_ref, data_all_ano_ref, \
                                    data_all_pi, data_all_ano, \
                                    data_all_pi - data_all_pi_ref, data_all_ano - data_all_ano_ref]
    
                        for idata,data1 in enumerate(data_plot):
                            if idata in [0]:
                                bounds = np.arange(var1_range[ivar][0],var1_range[ivar][1]+var1_range[ivar][2], var1_range[ivar][2])
                                unit = var1_units[ivar]
                                title = ref_case+' CNTL'
                            elif idata in [1]:
                                bounds = np.arange(var1_range_d[ivar][0],var1_range_d[ivar][1]+var1_range_d[ivar][2], var1_range_d[ivar][2])
                                unit = var1_units[ivar]+'/K'
                                title = ref_case
                            elif idata in [2]:
                                bounds = np.arange(var1_range[ivar][0],var1_range[ivar][1]+var1_range[ivar][2], var1_range[ivar][2])
                                unit = var1_units[ivar]
                                title = case+' CNTL'
                            elif idata in [3,4,5]:
                                bounds = np.arange(var1_range_d[ivar][0],var1_range_d[ivar][1]+var1_range_d[ivar][2], var1_range_d[ivar][2])
                                unit = var1_units[ivar]+'/K'
                                title = case
    
                                if idata == 4:
                                    bounds = bounds * 5.
                                    title = case+' CNTL minus '+ref_case+' CNTL'
                                elif idata == 5:
                                    title = case+' minus '+ref_case
    
                            ax = fig.add_subplot(nrow,ncol,idata+1)
                            im1 = ax.contourf(clats, levs, data1,\
                            bounds,\
                            cmap='RdBu_r',\
                            extend='both')
    
                            fig.colorbar(im1,ax=ax,fraction=0.05,label=unit)
    
                            ax.set_xlabel('Latitude')
                            ax.set_ylabel('Pressure [hPa]') 
    
                            ax.set_xticks(spec_clats)
                            ax.set_xticklabels(spec_lats)
    
                            ax.set_ylim(max(levs),min(levs))
                            if idata in [0,1]:
                                ax.set_title(svar,loc='left')
    
                            ax.set_title(title, loc='right')
    
                        fig.tight_layout()
                        fig.savefig(self.figdir+'LatLev-'+svar+'-'+case+'-vs-'+ref_case+'.png',bbox_inches='tight',dpi=300)
                        plt.close(fig)
                
        print('------------------------------------------------')
        print('plot_zm_CLOUD is done!')
        print('------------------------------------------------')
    
    #####################################################################
    ### plot_latlon_CLOUD
    #####################################################################
    def plot_latlon_CLOUD(self):
        print('plot_latlon_CLOUD starts ..........')
    
        cases_here = copy.deepcopy(self.cases)
        if 'amip-4xCO2' in self.cases:
            cases_here.remove('amip-4xCO2')
    
        var1 = ['TGCLDLWP','CLDLOW','TGCLDIWP']
        var1_tmp = ['clwvi', '', 'clivi']
        var1_range = [[0,100,5],[0,100,5],[0,100,5]]
        var1_range_d = [[-5,5,0.5], [-5,5,0.5],[-5,5,0.5]]
        var1_cntl_range_d = [[-20,20,2], [-20,20,2],[-20,20,2]]
    
        var1_out = ['Liquid Water Path','Low Cloud Fraction', 'Ice Water Path']
    
        var1_units = ['g/m$^2$','%','g/m$^2$']
    
        # ===============================================================        
        # ===============================================================        
        # E3SM
     
        for icase,case in enumerate(cases_here):
            ref_cases = self.ref_casesA[icase]
            if len(ref_cases) == 0:
                continue
    
            for ref_case in ref_cases:
    
                for ivar,svar in enumerate(var1):
                    # I don't have CLDLOW for v1_coupled. So, skip it if case or ref_case is it.
                    if svar == 'CLDLOW' and (ref_case == 'v1_coupled' or case == 'v1_coupled'):
                        continue
    
                    # -------------------- generate figures -------------------
                    nrow = 3
                    ncol = 2
                    fig = plt.figure(figsize=(ncol*5,nrow*2.5))
                    # -------------------- generate figures -------------------
    
                    f1 = cdms.open(self.datadir_v2+'global_cloud_'+case+'.nc')
                    if case == 'v1_coupled':
                        svar_in = var1_tmp[ivar]
                    else:
                        svar_in = var1[ivar]
    
                    tmp1 = f1(svar_in+'_pi_clim')
                    tmp2 = f1(svar_in+'_ano_clim')
    
                    lats = tmp1.getLatitude()[:]
                    lons = tmp1.getLongitude()[:]
    
                    data_all_pi = tmp1
                    data_all_ano =  tmp2
    
                    # convert unit from kg/m2 to g/m2
                    if svar in ['TGCLDLWP','TGCLDIWP']:
                        data_all_pi = data_all_pi * 1e3
                        data_all_ano = data_all_ano * 1e3
                    elif svar in ['CLDLOW']:
                        data_all_pi = data_all_pi * 1e2
                        data_all_ano = data_all_ano * 1e2
    
                    if svar in ['CLOUD'] and case == 'v1_coupled':
                        data_all_pi = data_all_pi/100.
                        data_all_ano = data_all_ano/100.
    
    
                    print('data_all_pi.shape=',data_all_pi.shape)
                    print('data_all_ano.shape=',data_all_ano.shape)
        
                    # -------------- read reference p4K data -----------------
                    fref = cdms.open(self.datadir_v2+'global_cloud_'+ref_case+'.nc')
                    if ref_case == 'v1_coupled':
                        svar_in = var1_tmp[ivar]
                    else:
                        svar_in = var1[ivar]
    
                    data_all_pi_ref = fref(svar_in+'_pi_clim')
                    data_all_ano_ref = fref(svar_in+'_ano_clim')
    
                    if svar in ['TGCLDLWP','TGCLDIWP']:
                        data_all_pi_ref = data_all_pi_ref * 1e3
                        data_all_ano_ref = data_all_ano_ref * 1e3
                    elif svar in ['CLDLOW']:
                        data_all_pi_ref = data_all_pi_ref * 1e2
                        data_all_ano_ref = data_all_ano_ref * 1e2
    
                    if svar in ['CLOUD'] and ref_case == 'v1_coupled':
                        data_all_pi_ref = data_all_pi_ref/100.
                        data_all_ano_ref = data_all_ano_ref/100.
    
                    #----------------------------------------------------------
                    # start plotting ...
                    #----------------------------------------------------------
                    data_plot = [data_all_pi_ref, data_all_ano_ref, \
                                data_all_pi, data_all_ano, \
                                data_all_pi - data_all_pi_ref, data_all_ano - data_all_ano_ref]
    
                    for idata,data1 in enumerate(data_plot):
                        if idata in [0,2]:
                            bounds = np.arange(var1_range[ivar][0],var1_range[ivar][1]+var1_range[ivar][2], var1_range[ivar][2])
                            unit = var1_units[ivar]
                        elif idata in [1,3,5]:
                            bounds = np.arange(var1_range_d[ivar][0],var1_range_d[ivar][1]+var1_range_d[ivar][2], var1_range_d[ivar][2])
                            unit = var1_units[ivar]+'/K'
                        elif idata in [4]:
                            bounds = np.arange(var1_cntl_range_d[ivar][0],var1_cntl_range_d[ivar][1]+var1_cntl_range_d[ivar][2], var1_cntl_range_d[ivar][2])
                            unit = var1_units[ivar]
    
                        if idata in [0]:
                            title = ref_case+' CNTL'
                        elif idata in [1]:
                            title = ref_case
                        elif idata in [2]:
                            title = case+' CNTL'
                        elif idata in [3]:
                            title = case
                        elif idata == 4:
                            title = case+' CNTL minus '+ref_case+' CNTL'
                        elif idata == 5:
                            title = case+' minus '+ref_case
    
                        print(svar, genutil.minmax(data1))
    
                        ax = fig.add_subplot(nrow,ncol,idata+1,projection=ccrs.Robinson(central_longitude=180.))
                        im1 = ax.contourf(lons, lats, data1,\
                        bounds,\
                        transform=ccrs.PlateCarree(),\
                        cmap='RdBu_r',\
                        extend='both')
    
                        avgdata1 = cdutil.averager(data1,axis='xy',weights='weighted')
                        print('avgdata1 = ',avgdata1)
    
                        ax.coastlines()
                        ax.set_global()
    
                        fig.colorbar(im1,ax=ax,fraction=0.025,label=unit)
    
    #                    if idata in [0,1]:
    #                        ax.set_title(var1_out[ivar],loc='left')
    
                        ax.set_title(title+'\n'+var1_out[ivar]+' ['+str(np.round(avgdata1,2))+']')
    
                    fig.tight_layout()
                    fig.savefig(self.figdir+'LatLon-'+svar+'-'+case+'-vs-'+ref_case+'.png',bbox_inches='tight',dpi=300)
                    plt.close(fig)
                
        print('------------------------------------------------')
        print('plot_latlon_CLOUD is done!')
        print('------------------------------------------------')
    
    #####################################################################
    ### plot_webb_decomposition
    #####################################################################
    def plot_webb_decomp(self):
        print('plot_webb_decomp starts ..........')
    
        cases_here = copy.deepcopy(self.cases)
        if 'amip-4xCO2' in self.cases:
            cases_here.remove('amip-4xCO2')
    
        variables = ['netCRE_lo','netCRE_nonlo','SWCRE_lo','SWCRE_nonlo','LWCRE_lo','LWCRE_nonlo']
        variables_out = []
        for svar in variables:
            if '_lo' in svar:
                variables_out.append('Low Cloud '+svar.split('_')[0])
            elif '_nonlo' in svar:
                variables_out.append('Non-low Cloud '+svar.split('_')[0])
        print(variables_out)
        
        nlat = 73
        nlon = 144
        
        # E3SM
        for icase,case in enumerate(cases_here):
            if case == 'v1_coupled':
                f1 = cdms.open(self.datadir_v1+'lat-lon-gfdbk-CMIP6-abrupt-4xCO2-E3SM-1-0_r1i1p1f1_1yr-150yr-webb-decomp.nc')
            elif case == 'v1_amip4K':
                f1 = cdms.open(self.datadir_v1+'lat-lon-gfdbk-CMIP6-amip-p4K-E3SM-1-0_r2i1p1f1_1yr-36yr-webb-decomp.nc')
            elif case == 'v1_future4K':
                f1 = cdms.open(self.datadir_v1+'lat-lon-gfdbk-CMIP6-amip-future4K-E3SM-1-0_r2i1p1f1_1yr-36yr-webb-decomp.nc')
            elif case == 'amip-4xCO2':
                continue
            else:
                f1 = cdms.open(self.datadir_v2+'lat-lon-gfdbk-CMIP6-'+case+'-webb-decomp.nc')
    
            nrow = 3; ncol = 2
            fig = plt.figure(figsize=(ncol*8,nrow*4)) # this creates and increases the figure size
    
            for ivar,svar in enumerate(variables):
                data = f1(svar)
                lats = data.getLatitude()
                lons = data.getLongitude()
    
                #----------------------------------------------------------
                # start plotting ...
                #----------------------------------------------------------
                bounds = np.arange(-3,3.25,0.25)
                cmap = plt.cm.RdBu_r
                bounds2 = np.append(np.append(-500,bounds),500) # This is only needed for norm if colorbar is extended
                norm = mpl.colors.BoundaryNorm(bounds2, cmap.N) # make sure the colors vary linearly even if the desired color boundaries are at varied intervals
    
                ax1 = fig.add_subplot(nrow,ncol,ivar+1,projection=ccrs.Robinson(central_longitude=180.))
                im1 = ax1.contourf(lons[:],lats[:],data,bounds,transform=ccrs.PlateCarree(),cmap=cmap,norm=norm,extend='both')
                ax1.coastlines()
                ax1.set_global()
    
                # global-mean 
                avgDATA = cdutil.averager(data,axis='xy',weights='weighted')
        
                case_out = case
    
                ax1.set_title(variables_out[ivar]+' ['+str(np.round(avgDATA,2))+']',fontsize=self.fh)
    
                PDF.make_colorbar(ax1, 'W/m$^2$/K', self.fh-5, im1, orientation='vertical')
    
            fig.subplots_adjust(top=0.9)
            plt.suptitle(case_out,fontsize=self.fh,y=0.95)
    
            fig.savefig(self.figdir+'LatLon-'+case+'-webb-decomp.png',bbox_inches='tight',dpi=300)
            plt.close(fig)
                
        print('------------------------------------------------')
        print('plot_web_decomp is done!')
        print('------------------------------------------------')
    
    
    ####################################################################################
    ### scatter plot for global mean CRE feedback: p4K vs future4K
    ####################################################################################
    def plot_CRE_globalmean_P4KvsFuture(self):
        print('ScatterPlot-CRE-feedback P4K.vs.Future starts ........')
    
        cases_here = copy.deepcopy(self.cases)
        if 'amip-4xCO2' in self.cases:
            cases_here.remove('amip-4xCO2')
    
        df_all = pd.DataFrame()
    
        if self.Add_otherCMIPs:
            # read other CMIP5 and CMIP6 models
    
            models_cmip6 = ['BCC-CSM2-MR','CNRM-CM6-1','IPSL-CM6A-LR',\
            'MRI-ESM2-0','CESM2','GFDL-CM4','CanESM5']
            models_cmip5 = ['bcc-csm1-1','CNRM-CM5','IPSL-CM5A-LR','IPSL-CM5B-LR',\
            'MIROC5','MPI-ESM-LR','MPI-ESM-MR','MRI-CGCM3','CanAM4']
            models_all = models_cmip6 + models_cmip5
    
    
            df_p4K = pd.DataFrame()
    
            for model in models_all:
                if model in models_cmip5:
                    filename = self.datadir_Ringer+'global_mean_features_CMIP5_amip4K_'+model+'_r1i1p1.csv'
                else:
                    if model == 'CNRM-CM6-1':
                        filename = self.datadir_Ringer+'global_mean_features_CMIP6_amip-p4K_'+model+'_r1i1p1f2.csv'
                    elif model == 'CanESM5':
                        filename = self.datadir_Ringer+'global_mean_features_CMIP6_amip-p4K_'+model+'_r1i1p2f1.csv'
                    else:
                        filename = self.datadir_Ringer+'global_mean_features_CMIP6_amip-p4K_'+model+'_r1i1p1f1.csv'
            
                df = pd.read_csv(filename,index_col=0)
                df.index = df.loc[:,'varname']
                df2 = df.loc[:,'anomaly_perK']
                df_p4K[model] = df2
    
            df_future = pd.DataFrame()
            for model in models_all:
                if model in models_cmip5:
                    filename = self.datadir_Ringer+'global_mean_features_CMIP5_amipFuture_'+model+'_r1i1p1.csv'
                else:
                    if model == 'CNRM-CM6-1':
                        filename = self.datadir_Ringer+'global_mean_features_CMIP6_amip-future4K_'+model+'_r1i1p1f2.csv'
                    elif model == 'CanESM5':
                        filename = self.datadir_Ringer+'global_mean_features_CMIP6_amip-future4K_'+model+'_r1i1p2f1.csv'
                    else:
                        filename = self.datadir_Ringer+'global_mean_features_CMIP6_amip-future4K_'+model+'_r1i1p1f1.csv'
            
                df = pd.read_csv(filename,index_col=0)
                df.index = df.loc[:,'varname']
                df2 = df.loc[:,'anomaly_perK']
                df_future[model] = df2
    
    
        # read amip
        for icase,case in enumerate(cases_here):
            if case == 'v1_coupled':
                # read v1-coupled 
                df_coupled = pd.read_csv(self.datadir_v1+'global_mean_features_CMIP6_abrupt-4xCO2_E3SM-1-0_r1i1p1f1.csv',index_col=0)
                df_coupled.index = df_coupled.loc[:,'varname']
                df_all['v1_coupled'] = df_coupled.loc[:,'anomaly_perK']
            elif case == 'v1_amip4K':
                # read v1-amip
                df_amip = pd.read_csv(self.datadir_v1+'global_mean_features_CMIP6_amip-p4K_E3SM-1-0_r2i1p1f1.csv',index_col=0)
                df_amip.index = df_amip.loc[:,'varname']
                df_all['v1_amip4K'] = df_amip.loc[:,'anomaly_perK']
            elif case == 'v1_future4K':
                # read v1-amip
                df_amip = pd.read_csv(self.datadir_v1+'global_mean_features_CMIP6_amip-future4K_E3SM-1-0_r2i1p1f1.csv',index_col=0)
                df_amip.index = df_amip.loc[:,'varname']
                df_all['v1_future4K'] = df_amip.loc[:,'anomaly_perK']
            elif case == 'amip-4xCO2':
                continue
            else:   
                df1 = pd.read_csv(self.datadir_v2+'global_mean_features_'+case+'.csv',index_col=0)
                df1.index = df1.loc[:,'varname']
        
                df2 = df1.loc[:,'anomaly_perK']
                df_all[case] = df2
        
        # start plotting 
        
        fig = plt.figure(figsize=(18,12))
        
        drop_index = ['tas','SWCLR','LWCLR']
        df_plot = df_all.drop(index=drop_index)
    
        if self.Add_otherCMIPs:
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
            ax.scatter(df_plot.loc[index,'amip-p4K'],df_plot.loc[index,'amip-future4K'],s=self.s1,alpha=self.a1,label='E3SM',color='red',marker='*')
    
            if self.Add_otherCMIPs:
                for icol,column in enumerate(df_p4K_plot.columns):
                    if column in models_cmip5:
                        L1 = ax.scatter(df_p4K_plot.loc[index,column],df_future_plot.loc[index,column],s=self.s1-100,alpha=self.a1,label='_nolegend_',color='tab:blue',edgecolor='none')
                    else:
                        if self.highlight_CESM2 and column == 'CESM2':
                            L1 = ax.scatter(df_p4K_plot.loc[index,column],df_future_plot.loc[index,column],s=self.s1,alpha=self.a1,label=column,color='tab:red',marker='X',edgecolor='none')
                        else:
                            L1 = ax.scatter(df_p4K_plot.loc[index,column],df_future_plot.loc[index,column],s=self.s1-100,alpha=self.a1,label='_nolegend_',color='tab:red',edgecolor='none')
    
            ax.plot([valmin,valmax],[valmin,valmax],ls='--',lw=3,color='grey')
    
            ax.tick_params(labelsize=self.fh)
            ax.set_ylabel(index+' Future4K [W/m$^2$/K]',fontsize=self.fh)
            ax.set_xlabel(index+' P4K [W/m$^2$/K]',fontsize=self.fh)
    
            ax.grid(which='major', linestyle=':', linewidth='1.0', color='grey')
    
            if idx == 0:
                ax.legend(fontsize=self.fh1)
    
    #    plt.xticks(x,df_plot.index)
        
        fig.tight_layout()
        fig.savefig(self.figdir+'ScatterPlot-CRE-P4KvsFuture-feedback-'+self.cases[-1]+'.png',bbox_inches='tight',dpi=300)
        plt.close(fig)
    
        del(df_all,df_plot)
        print('------------------------------------------------')
        print('ScatterPlot-CRE-P4KvsFuture-feedback is done!')
        print('------------------------------------------------')
    
    
    ####################################################################################
    ### 6. scatter plot for global mean radiative forcing: including E3SMv1 piControl and amip
    ####################################################################################
    def plot_RadForcing_globalmean(self):
    
        if any(case for case in self.cases if case == 'amip-4xCO2'):
            return 
    
        print('ScatterPlot-RadForcing starts ........')
    
        df_all = pd.DataFrame()
    
        if self.Add_otherCMIPs:
            # read other CMIP5 and CMIP6 models
    
            exp_cntl = [['piControl','amip'],['piControl','amip']]
            exp_new = [['abrupt4xCO2','amip4xCO2'],['abrupt-4xCO2','amip-4xCO2']]
            
            prefix = 'global_mean_features'
            suffix1 = '*.csv'
            suffix2 = '*.csv'
            
            models_all,cmip5_models,cmip6_models = PDF.get_intersect_withripf(exp_cntl,exp_new,prefix,suffix1,suffix2,self.datadir_Ringer)
    
    #        print('models_all',models_all,len(models_all))
    #        print('cmip5_models', cmip5_models,len(cmip5_models))
    #        print('cmip6_models', cmip6_models,len(cmip6_models))
    
            # ---- amip4xCO2 ---------------------
            df_4xCO2 = pd.DataFrame()
            for model in models_all:
                if model in cmip5_models:
    
                    if model == 'CanESM2_r1i1p1':
                        model_amip = 'CanAM4_r1i1p1'
                    elif model == 'HadGEM2-ES_r1i1p1':
                        model_amip = 'HadGEM2-A_r1i1p1'
                    else:
                	    model_amip = model
    
                    filename = self.datadir_Ringer+'global_mean_features_CMIP5_amip4xCO2_'+model_amip+'.csv'
                else:
                    filename = self.datadir_Ringer+'global_mean_features_CMIP6_amip-4xCO2_'+model+'.csv'
            
                df = pd.read_csv(filename,index_col=0)
                df.index = df.loc[:,'varname']
                df2 = df.loc[:,'anomaly']
                df_4xCO2[model] = df2
    
        # ------ read E3SM amip4xCO2 ------------------------
        for icase,case in enumerate(self.cases):
            if case == 'v1_coupled':
                # read v1-coupled 
                df_coupled = pd.read_csv(self.datadir_v1+'global_mean_features_CMIP6_abrupt-4xCO2_E3SM-1-0_r1i1p1f1.csv',index_col=0)
                df_coupled.index = df_coupled.loc[:,'varname']
                df_all['v1_coupled'] = df_coupled.loc[:,'forcing']
            elif case not in ['amip-4xCO2','v1_coupled']:
                continue
            else:   
                df1 = pd.read_csv(self.datadir_v2+'global_mean_features_'+case+'.csv',index_col=0)
                df1.index = df1.loc[:,'varname']
        
                df2 = df1.loc[:,'anomaly']
                df_all[case] = df2
        
    
        # start plotting 
        fig = plt.figure(figsize=(18,12))
        ax = fig.add_subplot(1,1,1)
        
        drop_index = ['tas','SWCLR','LWCLR']
        df_plot = df_all.drop(index=drop_index)
    
        if self.Add_otherCMIPs:
            df_4xCO2_plot = df_4xCO2.drop(index=drop_index)
       
        x = np.arange(1,len(df_plot.index)+1,1)
        
        for idx,index in enumerate(df_plot.index):
            for icol,column in enumerate(df_plot.columns):
                if column == 'v1_coupled':
                    L1 = ax.scatter(x[idx],df_plot.loc[index,column].tolist(),s=self.s1,alpha=self.a1,label=column,color=self.colors[icol],marker='x')
                else:
                    ax.scatter(x[idx],df_plot.loc[index,column].tolist(),s=self.s1,alpha=self.a1,label=column,color=self.colors[icol])
        
            if self.Add_otherCMIPs:
                # add other CMIP models
                for icol,column in enumerate(df_4xCO2_plot.columns):
                    if self.highlight_CESM2 and 'CESM2' in column:
                        ax.scatter(x[idx]-0.2,df_4xCO2_plot.loc[index,column].tolist(),edgecolor='none',facecolor=self.lc_CESM2,alpha=self.a1,s=self.s2,marker='X',\
                        label=column.split('_')[0]+'_amip-p4K')
                    else:
                        ax.scatter(x[idx]-0.2,df_4xCO2_plot.loc[index,column].tolist(),edgecolor='none',facecolor='grey',alpha=self.a1,s=self.s2,\
                        label='_nolegend_')
    
                # ensemble mean
                L2 = ax.scatter(x[idx]-0.2,df_4xCO2_plot.loc[index,:].mean().tolist(),color='black',s=self.s2)
                
            ax.tick_params(labelsize=self.fh)
            ax.set_ylabel('Forcing/Adjustment [W/m$^2$]',fontsize=self.fh)
            if idx==0:
                legend1 = ax.legend([L2],['amip4xCO2'],fontsize=self.fh1,loc='lower right')
                ax.legend(fontsize=self.fh1)
                ax.add_artist(legend1) 
        
        plt.xticks(x,df_plot.index)
        
        ax.grid(which='major', linestyle=':', linewidth='1.0', color='grey')
        fig.savefig(self.figdir+'ScatterPlot-RadForcing-'+self.cases[-1]+'.png',bbox_inches='tight',dpi=300)
        plt.close(fig)
    
        del(df_all,df_plot)
        print('------------------------------------------------')
        print('ScatterPlot-RadForcing is done!')
        print('------------------------------------------------')




