
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
import cal_EIS as calEIS
import sort_cloud_regime as SCR
import sort_cloud_3regime as SCR3

import glob

import LonPivot as LP
#############################################################################################
def prepare_pd2html(outfig,varname,desc,casevscase=''):
    print(casevscase)
    pd_plot = pd.DataFrame([[varname,desc,casevscase,outfig]],columns=['Variables','Description','Case.VS.Case','Plot'])
    pd_plot['Plot'] = pd_plot['Plot'].apply(lambda x: '<a href='+outfig+' target="_blanck">plot</a>')

    return pd_plot

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
    dics_cal['cal_EIS']             = my_cal.cal_EIS
    dics_cal['sort_cloud_regime']   = my_cal.sort_cloud_regime
    dics_cal['sort_cloud_3regime']   = my_cal.sort_cloud_3regime


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

    def cal_EIS(self):
        result = calEIS.cal_EIS(self.direc_data,self.case_stamp,self.yearS2,self.yearE2,self.run_id1,self.run_id2,self.outdir_final,self.figdir,self.exp1,self.exp2)

    def sort_cloud_regime(self):
        result = SCR.sort_cloud_regime(self.direc_data,self.case_stamp,self.yearS2,self.yearE2,self.run_id1,self.run_id2,self.outdir_final,self.figdir,self.exp1,self.exp2)

    def sort_cloud_3regime(self):
        result = SCR3.sort_cloud_3regime(self.direc_data,self.case_stamp,self.yearS2,self.yearE2,self.run_id1,self.run_id2,self.outdir_final,self.figdir,self.exp1,self.exp2)


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

    dics_plots['RadKernel_latlon_dif']           = my_plot.plot_RadKernel_latlon_dif
    dics_plots['CldRadKernel_latlon_dif']        = my_plot.plot_CldRadKernel_latlon_dif

    dics_plots['tas_latlon']                 = my_plot.plot_tas_latlon

    dics_plots['LCF']                        = my_plot.plot_LCF 
    dics_plots['zm_CLOUD']                   = my_plot.plot_zm_CLOUD
    dics_plots['latlon_CLOUD']               = my_plot.plot_latlon_CLOUD
    dics_plots['webb_decomp']                = my_plot.plot_webb_decomp
    dics_plots['CLOUD_profile']              = my_plot.plot_CLOUD_profile
    dics_plots['NRMSE_RadKern']              = my_plot.plot_NRMSE_RadKern
    dics_plots['cal_regionCor']              = my_plot.cal_RadKernel_regional_correspondence
    
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
        fig.savefig(self.figdir+'ScatterPlot-CRE-feedback_'+self.cases[-1]+'.png',bbox_inches='tight',dpi=300)
        plt.close(fig)
        
        #<2021-12-21
        pd_plot = prepare_pd2html('../figure/ScatterPlot-CRE-feedback_'+self.cases[-1]+'.png',
                                  'Radiative feedback [W/m2/K]',
                                  'TOA radiative feedbacks [SWCRE, LWCRE, netCRE, clear-sky SW, clear-sky LW, net TOA, net TOA SW, net TOA LW, clear-sky net TOA, clear-sky net TOA SW, clear-sky net TOA LW]',
                                  '')
        #>2021-12-21
   
        del(df_all,df_plot)
        print('------------------------------------------------')
        print('ScatterPlot-CRE-feedback is done!')
        print('------------------------------------------------')

        return pd_plot 
    
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
    
            # if only add those models with both amip and cmip runs, set do_amip_cmip = True
            do_amip_cmip = False

            phases = ['CMIP5','CMIP6']
    
            if do_amip_cmip:
                exp_cntl = [['piControl','amip'],['piControl','amip']]
                exp_new = [['abrupt4xCO2','amip4K'],['abrupt-4xCO2','amip-p4K']]

                suffix1 = '*_1yr-150yr.csv'
                suffix2 = '*.csv'

            else:
                exp_cntl = [['piControl','piControl'],['piControl','piControl']]
                exp_new = [['abrupt4xCO2','abrupt4xCO2'],['abrupt-4xCO2','abrupt-4xCO2']]
            
                suffix1 = '*_1yr-150yr.csv'
                suffix2 = '*_1yr-150yr.csv'
            
            prefix = 'FDBK'

            models_all,cmip5_models,cmip6_models = PDF.get_intersect(exp_cntl,exp_new,prefix,suffix1,suffix2,self.datadir_RadKernel+'/20220208/')
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
    
                    if do_amip_cmip:
                        if model == 'CanESM2':
                            model_amip = 'CanAM4'
                        elif model == 'HadGEM2-ES':
                            model_amip = 'HadGEM2-A'
                        else:
                	        model_amip = model
    
                        df = pd.read_csv(self.datadir_RadKernel+'FDBK_'+phase+'_'+exp_new[iphase][1]+'_'+model_amip+suffix+'.csv',index_col=0)
                    else:
                        ff = glob.glob(self.datadir_RadKernel+'/20220208/FDBK_'+phase+'_'+exp_new[iphase][1]+'_'+model+'_*_1yr-150yr.csv')
                        if len(ff) == 0:
                            print('We dont find data for ',model,'please check!!!')
                        else:
                            print(ff)
                            df = pd.read_csv(ff[0],index_col=0)

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
        
        fig = plt.figure(figsize=(18,9))
        ax = fig.add_subplot(1,1,1)
        
        drop_index = ['T','dLW_adj','dSW_adj','dnet_adj','LW_resd','SW_resd',\
        #'net_resd',\
        'T_clr','Planck_clr','LR_clr','WV_clr','ALB_clr','WV_clr_SW','WV_clr_LW','WV_SW','WV_LW',\
        'SWCRE','LWCRE','netCRE','Planck_clr_fxRH','LR_clr_fxRH','RH_clr','LW_clr_sum','SW_clr_sum',\
        'net_clr_sum','LW_clr_dir','SW_clr_dir','net_clr_dir','LW_cld_sum','SW_cld_sum',\
        #,'net_cld_sum',\
        'LW_cld_dir','SW_cld_dir','net_cld_dir','LW_clr_resd','SW_clr_resd','net_clr_resd']
        
    
        df_plot = df_all.drop(index=drop_index)
        x = np.arange(1,len(df_plot.index)+1,1)
        print(df_plot.index)

        # redefine column orders 
        indexA = ['Planck','LR','WV','ALB','netCRE_adj','SWCRE_adj','LWCRE_adj','net_cld_sum','net_resd','Planck_fxRH','LR_fxRH','RH']
        xticks = ['Planck','LR','WV','Albedo','Cloud','Cloud$_{sw}$','Cloud$_{lw}$','Total','Residual','Planck\n[fixed RH]','LR\n[fixed RH]','RH']

        # output the dataframe df_plot as a table
        df_plot_flip = df_plot.transpose()
        df_plot_flip = df_plot_flip.reindex(columns=indexA)
        df_plot_flip.columns = xticks
        df_plot_flip.index = [item.split('.')[-1] for item in df_plot_flip.index]
        print(df_plot_flip.round(2))
        df_plot_flip_out = df_plot_flip.round(2)
        df_out = df_plot_flip_out.applymap("{:.2f}".format)
        print(df_out)
        df_out.to_csv(self.datadir_v2+'climate_feedback_table_'+cases_here[-1]+'.csv')
    
        if self.Add_otherCMIPs:
            df_others_plot = df_others.drop(index=drop_index)
        
        for idx,index in enumerate(indexA):
            for icol,column in enumerate(df_plot.columns):
                if column == 'v1_coupled' or column == 'v2_coupled':
                    label = column.split('_')[0]+" [abrupt4xCO2]"
                    if 'v2.NARRM.coupled' in df_plot.columns and column == 'v2_coupled':
                        label = 'v2.LR'
                    L1 = ax.scatter(x[idx],df_plot.loc[index,column].tolist(),s=self.s1,alpha=self.a1,label=label,color=self.colors[icol],marker='x')
                elif column == 'v1_amip4K':
                    ax.scatter(x[idx],df_plot.loc[index,column].tolist(),s=self.s1,alpha=self.a1,label=column,color=self.colors[icol],marker='x')
                elif column == 'v1_future4K':
                    ax.scatter(x[idx],df_plot.loc[index,column].tolist(),s=self.s1,alpha=self.a1,label=column,color=self.colors[icol],marker='x')
                else:
                    if 'gwenergy' in column:
                        label = 'All'
                    elif column == 'v2.NARRM.coupled':
                        label = 'v2.NARRM'
                    else:
                        label = column.split('.')[-1]

                    ax.scatter(x[idx]+icol/20.-len(df_plot.columns)/2.*1/20.,df_plot.loc[index,column].tolist(),s=self.s1,alpha=self.a1,label=label,color=self.colors[icol],edgecolor='grey')
               
            # other CMIP models
            if self.Add_otherCMIPs:
                for icol,column in enumerate(df_others_plot.columns):
                    if self.highlight_CESM2 and 'CESM2' in column:
                        ax.scatter(x[idx], df_others_plot.loc[index,column].tolist(),s=self.s2,edgecolor='none',facecolor=self.lc_CESM2,alpha=1, marker='X',\
                        label = column.split('_')[0]+'_amip-p4K')
                    else:
                        ax.scatter(x[idx]-0.2, df_others_plot.loc[index,column].tolist(),s=self.s2,edgecolor='none',facecolor='grey',alpha=0.3,marker='x')
    
                # ensemble mean
                L2 = ax.scatter(x[idx]-0.2, df_others_plot.loc[index,:].mean(),s=self.s2,edgecolor='black',facecolor='black',marker='X')
    
            ax.tick_params(labelsize=self.fh)
            ax.set_ylabel('Feedback [W/m$^2$/K]',fontsize=self.fh)
            if idx == 0:
                if self.Add_otherCMIPs:
                    if do_amip_cmip:
                        label = 'amip4K'
                    else:
                        label = 'abrupt4xCO2'

                    legend1 = ax.legend([L2],[label],fontsize=self.fh1,loc='upper left')
                    ax.legend(fontsize=self.fh1)
                    ax.add_artist(legend1) 
                else:
                    ax.legend(fontsize=self.fh1)
    
        #ax.grid(which='major', linestyle=':', linewidth='1.0', color='grey')
        ax.axhline(y=0,ls='-',color='black',lw=1)
        degrees = 0

        plt.xticks(x,xticks,rotation=degrees)
        ax.set_title('Global Mean Feedbacks',fontsize=self.fh)

        # Add reference line between Cloud and Cloud_sw
        ax.axvline(x=5.5,ls='--',color='grey')
        # between Cloud_lw and Total 
        ax.axvline(x=7.5,ls='--',color='grey')
        # between Total and Residual
        ax.axvline(x=8.5,ls='--',color='grey')
        # between Residual and Planck_fxRH
        ax.axvline(x=9.5,ls='--',color='grey')
        
        fig.savefig(self.figdir+'ScatterPlot-RadKernel-Feedback_'+self.cases[-1]+'.png',bbox_inches='tight',dpi=300)
        plt.close(fig)

        #<2021-12-22
        pd_plot = prepare_pd2html('../figure/ScatterPlot-RadKernel-Feedback_'+self.cases[-1]+'.png',
                                  'CRE [W/m2/K]',
                                  'Adjusted CRE feedbacks derived from radiative kernel method (SW, LW, NET)',
                                  '')
        #>2021-12-22
    
        print('------------------------------------------------')
        print('ScatterPlot-RadKernel-Feedback is done!')
        print('------------------------------------------------')
        return pd_plot 

    
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
        
        pd_plot_all = pd.DataFrame(columns=['Variables','Description','Case.VS.Case','Plot'])

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
            fig.savefig(self.figdir+'Zonal-mean-Cloud-RadKernel-Feedback_'+str(np.round(ii,0))+'-'+self.cases[-1]+'.png',bbox_inches='tight',dpi=300)
            plt.close(fig)

            pd_plot = prepare_pd2html('../figure/Zonal-mean-Cloud-RadKernel-Feedback_'+str(np.round(ii,0))+'-'+self.cases[-1]+'.png',
                                      'CRE [W/m2/K]',
                                      'Adjusted CRE feedbacks derived from radiative kernel method (SW, LW, NET)',
                                      'ncases = '+str(np.round(ii,0)))
    
            pd_plot_all = pd.merge(pd_plot_all, pd_plot, on =['Variables','Description','Case.VS.Case','Plot'],how='outer')

        print('------------------------------------------------')
        print('plot_RadKernel_zonalmean is done!')
        print('------------------------------------------------')

        return pd_plot_all
    
    
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
        pd_plot_all = pd.DataFrame(columns=['Variables','Description','Case.VS.Case','Plot'])

        for ii in self.ncase:
            print('ii = ',ii)
    
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
                    print('Doing LW...')
                    y1 = df_LW_all.iloc[jj*5:(jj+1)*5,icol]
                    if column == 'v1_coupled' or column == 'v2_coupled':
                        La = axes[jj].scatter(x-w+w2,y1.values.tolist(),marker='v',s=self.s1,color=self.colors[icol],alpha=self.a1,label=column)
                    elif column == 'v1_amip4K':
                        axes[jj].scatter(x-w+w2,y1.values.tolist(),marker='x',s=self.s1,color=self.colors[icol],alpha=self.a1,label=column)
                    elif column == 'v1_future4K':
                        axes[jj].scatter(x-w+w2,y1.values.tolist(),marker='x',s=self.s1,color=self.colors[icol],alpha=self.a1,label=column)
    
                    else:
        #                L1 = axes[jj].scatter(x-w+w2,y1.values.tolist(),marker='o',s=self.s2,color=self.colors[icol],alpha=self.a1,label=column)
                        La = axes[jj].scatter(x-w+w2,y1.values.tolist(),marker='v',s=self.s2,color=self.colors[icol],alpha=self.a1,label=column)
        
                for icol,column in enumerate(df_net_all.columns):
                    print('Doing NET...')
                    y1 = df_net_all.iloc[jj*5:(jj+1)*5,icol]
                    if column == 'v1_coupled' or column == 'v2_coupled':
                        Lb = axes[jj].scatter(x+w2,y1.values.tolist(),marker='o',s=self.s1,color=self.colors[icol],alpha=self.a1,label='_nolegend_')
                    elif column == 'v1_amip4K':
                        axes[jj].scatter(x+w2,y1.values.tolist(),marker='x',s=self.s1,color=self.colors[icol],alpha=self.a1,label='_nolegend_')
                    elif column == 'v1_future4K':
                        axes[jj].scatter(x+w2,y1.values.tolist(),marker='x',s=self.s1,color=self.colors[icol],alpha=self.a1,label='_nolegend_')
    
                    else:
        #                L2 = axes[jj].scatter(x+w2,y1.values.tolist(),marker='o',s=self.s2,color=self.colors[icol],alpha=self.a1,label=column)
                        Lb = axes[jj].scatter(x+w2,y1.values.tolist(),marker='o',s=self.s2,color=self.colors[icol],alpha=self.a1,label='_nolegend_')
        
                for icol,column in enumerate(df_SW_all.columns):
                    print('Doing SW...')
                    y1 = df_SW_all.iloc[jj*5:(jj+1)*5,icol]
                    if column == 'v1_coupled' or column == 'v2_coupled':
                        Lc = axes[jj].scatter(x+w+w2,y1.values.tolist(),marker='^',s=self.s1,color=self.colors[icol],alpha=self.a1,label='_nolegend_')
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
            fig.savefig(self.figdir+'ScatterPlot-Cloud-feedback-Decomposition_'+str(np.round(ii,0))+'-'+self.cases[-1]+'.png',bbox_inches='tight',dpi=300)
            plt.close(fig)

            pd_plot = prepare_pd2html('../figure/ScatterPlot-Cloud-feedback-Decomposition_'+str(np.round(ii,0))+'-'+self.cases[-1]+'.png',
                                      'Cloud feedback [W/m2/K]',
                                      'Cloud feedback components derived from cloud radiative kernel method (ISCCP simulator output): Low clouds and non-low clouds; Cloud amount, altitude and optical depth',
                                      'ncases = '+str(np.round(ii,0)))

            pd_plot_all = pd.merge(pd_plot_all, pd_plot, on =['Variables','Description','Case.VS.Case','Plot'],how='outer')

    
        print('------------------------------------------------')
        print('ScatterPlot-Cloud-feedback-Decomposition is done!')
        print('------------------------------------------------')

        return pd_plot_all
    
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

        pd_plot_all_1 = pd.DataFrame(columns=['Variables','Description','Case.VS.Case','Plot'])
        pd_plot_all = pd.DataFrame(columns=['Variables','Description','Case.VS.Case','Plot'])
    
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
                    fig.savefig(self.figdir+'ZonalMean-Cloud-feedback-Decomposition_'+lev+'-'+component+'-'+str(np.round(ii,0))+'-'+self.cases[-1]+'.png',bbox_inches='tight',dpi=300)
                    plt.close(fig)

                    pd_plot = prepare_pd2html('../figure/ZonalMean-Cloud-feedback-Decomposition_'+lev+'-'+component+'-'+str(np.round(ii,0))+'-'+self.cases[-1]+'.png',
                                      component+'-'+lev+' [W/m2/K]',
                                      components_out[icomp]+' '+levs_out[ilev]+' Feedback',
                                      'ncases = '+str(np.round(ii,0)))

                    pd_plot_all = pd.merge(pd_plot_all, pd_plot, on =['Variables','Description','Case.VS.Case','Plot'],how='outer')

    
            #<qinyi 2021-05-19 #------------------
            plt.tight_layout()
            fig1.savefig(self.figdir+'ZonalMean-Cloud-feedback-Decomposition-KeyComponents_'+str(np.round(ii,0))+'-'+self.cases[-1]+'.png',bbox_inches='tight',dpi=300)
            plt.close(fig)
            #>qinyi 2021-05-19 #------------------

            pd_plot = prepare_pd2html('../figure/ZonalMean-Cloud-feedback-Decomposition-KeyComponents_'+str(np.round(ii,0))+'-'+self.cases[-1]+'.png',
                              'Key components [W/m2/K]',
                              'LO680_SWcld_amt, LO680_SWcld_tau, HI680_SWcld_amt, HI680_SWcld_tau, HI680_LWcld_amt, HI680_LWcld_alt and HI680_LWcld_tau',
                              'ncases = '+str(np.round(ii,0)))

            pd_plot_all_1 = pd.merge(pd_plot_all_1, pd_plot, on =['Variables','Description','Case.VS.Case','Plot'],how='outer')

        
        print('------------------------------------------------')
        print('ZonalMean-Cloud-feedback-Decomposition is done!')
        print('------------------------------------------------')

        return pd_plot_all
    
    
    #####################################################################
    ### LAT-LON cloud feedback difference based on cloud radiative kernel
    #####################################################################
    def plot_CldRadKernel_latlon_dif(self):
    
        print('plot_CldRadKernel_latlon_dif starts ........')
    
        cases_here = copy.deepcopy(self.cases)
        if 'amip-4xCO2' in self.cases:
            cases_here.remove('amip-4xCO2')
    
    
        nlat = 90
        nlon = 144
        
        # define DataFrame to save spatial correlation and NRMSE for further plot
        pdata = pd.DataFrame(columns=['case','ref_case','lev','var','COR','NRMSE'])

        # define DataFrame to save figures to generate html file
        pd_plot_all = pd.DataFrame(columns=['Variables','Description','Case.VS.Case','Plot'])

        # --- define variables 
        sections = ['ALL','HI680','LO680']
        components = ['NET','SW','LW']
        decomps = ['tot','amt','alt','tau']

        sections_out = ['ALL', 'High cloud', 'Low cloud']
        components_out = ['NET','SW','LW']
        
        # generate figure based on case categories
        for icomp,component in enumerate(components):
            for isec,sec in enumerate(sections):
    
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

                            ##pdata = pd.DataFrame(columns=['case','ref_case','var','COR','NRMSE'])
                            pd_data = pd.DataFrame([[cases_here[icase],cases_here[iref],sec,name,cor,NRMSE]],columns=['case','ref_case','lev','var','COR','NRMSE'])
                            pdata = pd.merge(pdata, pd_data, on = ['case','ref_case','lev','var','COR','NRMSE'], how='outer')
                            print(pdata)
 
                            plt.title(name+' ['+str(np.round(avgDATA,3))+']\nNRMSE='+str(np.round(NRMSE,2))+', COR='+str(np.round(cor,2)),fontsize=self.fh)

                            cb = plt.colorbar(im1,orientation='vertical',drawedges=True,ticks=bounds[::2])
                            cb.set_label('W/m$^2$/K')
    
        
                        fig.subplots_adjust(top=0.9)
    
                        fig.savefig(self.figdir+'LatLon-Cloud-feedback-Decomposition_'+cases_here[icase]+'.minus.'+cases_here[iref]+'-'+component+'-'+sec+'.png',dpi=300,bbox_inches='tight')
                        plt.close(fig)

                        pd_plot = prepare_pd2html('../figure/LatLon-Cloud-feedback-Decomposition_'+cases_here[icase]+'.minus.'+cases_here[iref]+'-'+component+'-'+sec+'.png',
                                      component+'-'+sec+' [W/m2/K]',
                                      components_out[icomp]+' '+sections_out[isec]+' feedback',
                                      cases_here[icase]+' .VS. '+cases_here[iref])

                        pd_plot_all = pd.merge(pd_plot_all, pd_plot, on =['Variables','Description','Case.VS.Case','Plot'],how='outer')

        print(pdata)
        # save into csv file for further plot 
        pdata.to_csv(self.datadir_v2+'COR_RMSE_CldRadKernel_latlon_dif_'+cases_here[-1]+'.csv')

        print('------------------------------------------------')
        print('plot_CldRadKernel_latlon_dif is done!')
        print('------------------------------------------------')

        return pd_plot_all

    
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
            pd_plot_all = pd.DataFrame(columns=['Variables','Description','Case.VS.Case','Plot'])

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
    
                    fig.savefig(self.figdir+'LatLon-TAS_'+cases_here[icase]+'.minus.'+cases_here[iref]+'.png',bbox_inches='tight',dpi=300)
                    plt.close(fig)
                    
                    pd_plot = prepare_pd2html('../figure/LatLon-TAS_'+cases_here[icase]+'.minus.'+cases_here[iref]+'.png',
                                      'TAS [K/K]',
                                      'surface air temperature response',
                                      cases_here[icase]+' .VS. '+cases_here[iref])

                    pd_plot_all = pd.merge(pd_plot_all, pd_plot, on =['Variables','Description','Case.VS.Case','Plot'],how='outer')
                
        print('------------------------------------------------')
        print('plot_tas_latlon is done!')
        print('------------------------------------------------')

        return pd_plot_all
    
    #####################################################################
    ### LAT-LON RadKernel difference
    #####################################################################
    def plot_RadKernel_latlon_dif(self):
        print('plot_RadKernel_latlon_dif starts ..........')
        
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
        # define DataFrame to save spatial correlation and NRMSE for further plot
        pdata = pd.DataFrame(columns=['case','ref_case','var','COR','NRMSE'])

        pd_plot_all = pd.DataFrame(columns=['Variables','Description','Case.VS.Case','Plot'])

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
                #<qinyi 2021-12-30 #------------------
                #n = 11
                #new_cmap = PDF.map_white(n,0.1,cmap)
                new_cmap = cmap 

                bounds2 = np.append(np.append(-500,bounds),500) # This is only needed for norm if colorbar is extended
                norm = mpl.colors.BoundaryNorm(bounds2, new_cmap.N) # make sure the colors vary linearly even if the desired color boundaries are at varied intervals

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
                        im1 = ax1.contourf(lons[:],lats[:],DATA,bounds,transform=ccrs.PlateCarree(),cmap=new_cmap,norm=norm,extend='both')
                        #im1 = ax1.contourf(lons[:],lats[:],DATA,bounds,transform=ccrs.PlateCarree(),cmap=new_cmap,extend='both')

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

                            ##pdata = pd.DataFrame(columns=['case','ref_case','var','COR','NRMSE'])
                            pd_data = pd.DataFrame([[cases_here[icase],cases_here[iref],variables_out[ivar],cor,NRMSE]],columns=['case','ref_case','var','COR','NRMSE'])
                            pdata = pd.merge(pdata, pd_data, on = ['case','ref_case','var','COR','NRMSE'], how='outer')
                            print(pdata)
    
                            ax1.set_title(title_plot[num0]+'\n'+variables_out[ivar]+' ['+str(np.round(avgDATA,2))+']\nNRMSE='+str(np.round(NRMSE,2))+', COR='+str(np.round(cor,2)),fontsize=self.fh)

                            PDF.make_colorbar(ax1, 'W/m$^2$/K', self.fh-5, im1, orientation='vertical')

                        else:
                            ax1.set_title(title_plot[num0]+'\n'+variables_out[ivar]+' ['+str(np.round(avgDATA,2))+']',fontsize=self.fh)

                        num0 += 1
    
                fig.subplots_adjust(top=0.9)
    
                fig.savefig(self.figdir+'LatLon-adjusted-CRE_'+case_out+'.minus.'+ref_case_out+'.png',bbox_inches='tight',dpi=300)
                plt.close(fig)

                pd_plot = prepare_pd2html('../figure/LatLon-adjusted-CRE_'+case_out+'.minus.'+ref_case_out+'.png',
                                      'CRE [W/m2/K]',
                                      'Adjusted CRE feedback derived from radiative kernel method (SW, LW, and NET)',
                                      case_out+' .VS. '+ref_case_out)

                pd_plot_all = pd.merge(pd_plot_all, pd_plot, on =['Variables','Description','Case.VS.Case','Plot'],how='outer')

                
        print(pdata)
        # save into csv file for further plot
        pdata.to_csv(self.datadir_v2+'COR_RMSE_RadKernel_latlon_dif_'+cases_here[-1]+'.csv')

        print('------------------------------------------------')
        print('plot_RadKernel_latlon_dif is done!')
        print('------------------------------------------------')

        return pd_plot_all
    
    
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
    
        fig.savefig(self.figdir+'LCF-vs-temperature_'+cases_here[-1]+'.png',bbox_inches='tight',dpi=300)
        plt.close(fig)

        pd_plot = prepare_pd2html('../figure/LCF-vs-temperature_'+cases_here[-1]+'.png',
                                      'LCF [fraction]',
                                      'The ratio of Liquid cloud water to total cloud water',
                                      '')
                
        print('------------------------------------------------')
        print('plot_LCF is done!')
        print('------------------------------------------------')

        return pd_plot
    
    
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
     
        pd_plot_all = pd.DataFrame(columns=['Variables','Description','Case.VS.Case','Plot'])

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
                                if idata == 4:
                                    unit = var1_units[ivar]
                                else:
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
                        fig.savefig(self.figdir+'LatLev_'+svar+'-'+case+'-vs-'+ref_case+'.png',bbox_inches='tight',dpi=300)
                        plt.close(fig)

                        pd_plot = prepare_pd2html('../figure/LatLev_'+svar+'-'+case+'-vs-'+ref_case+'.png',
                                      svar+' ['+var1_units[ivar]+'/K]',
                                      'Zonal mean '+svar,
                                      case+' .VS. '+ref_case)

                        pd_plot_all = pd.merge(pd_plot_all, pd_plot, on =['Variables','Description','Case.VS.Case','Plot'],how='outer')
                
        print('------------------------------------------------')
        print('plot_zm_CLOUD is done!')
        print('------------------------------------------------')

        return pd_plot_all
    
    #####################################################################
    ### plot_latlon_CLOUD
    #####################################################################
    def plot_latlon_CLOUD(self):
        print('plot_latlon_CLOUD starts ..........')
    
        cases_here = copy.deepcopy(self.cases)
        if 'amip-4xCO2' in self.cases:
            cases_here.remove('amip-4xCO2')
    
        var1 = ['TGCLDLWP','CLDLOW','TGCLDIWP','CLDMED','CLDHGH']
        var1_tmp = ['clwvi', '', 'clivi','','']
        var1_range = [[0,100,5],[0,100,5],[0,100,5],[0,100,5],[0,100,5]]
        var1_range_d = [[-5,5,0.5], [-5,5,0.5],[-5,5,0.5],[-5,5,0.5],[-5,5,0.5]]
        var1_cntl_range_d = [[-20,20,2], [-20,20,2],[-20,20,2],[-20,20,2],[-20,20,2]]
    
        var1_out = ['Liquid Water Path','Low Cloud Fraction', 'Ice Water Path','Middle Cloud Fraction', 'High Cloud Fraction']
    
        var1_units = ['g/m$^2$','%','g/m$^2$','%','%']
        var1_units1 = ['g/m2','%','g/m2','%','%']
    
        # ===============================================================        
        # ===============================================================        
        # E3SM
        pd_plot_all = pd.DataFrame(columns=['Variables','Description','Case.VS.Case','Plot'])
 
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
                    elif svar in ['CLDLOW','CLDMED','CLDHGH']:
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
                    elif svar in ['CLDLOW','CLDMED','CLDHGH']:
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
                    fig.savefig(self.figdir+'LatLon_'+svar+'-'+case+'-vs-'+ref_case+'.png',bbox_inches='tight',dpi=300)
                    plt.close(fig)

                    pd_plot = prepare_pd2html('../figure/LatLon_'+svar+'-'+case+'-vs-'+ref_case+'.png',
                                      svar+' ['+var1_units1[ivar]+'/K]',
                                      svar,
                                      case+' .VS. '+ref_case)

                    pd_plot_all = pd.merge(pd_plot_all, pd_plot, on =['Variables','Description','Case.VS.Case','Plot'],how='outer')

                
        print('------------------------------------------------')
        print('plot_latlon_CLOUD is done!')
        print('------------------------------------------------')

        return pd_plot_all
    
    #####################################################################
    ### plot_webb_decomposition
    #####################################################################
    def plot_webb_decomp(self):
        print('plot_webb_decomp starts ..........')
    
        cases_here = copy.deepcopy(self.cases[-1:])
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
        pd_plot_all = pd.DataFrame(columns=['Variables','Description','Case.VS.Case','Plot'])

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
    
            fig.savefig(self.figdir+'LatLon_'+case+'_webb-decomp.png',bbox_inches='tight',dpi=300)
            plt.close(fig)

            pd_plot = prepare_pd2html('../figure/LatLon_'+case+'_webb-decomp.png',
                                      'CRE [W/m2/K]',
                                      'Cloud feedback components derived from Webb method',
                                      case)

            pd_plot_all = pd.merge(pd_plot_all, pd_plot, on =['Variables','Description','Case.VS.Case','Plot'],how='outer')
                
        print('------------------------------------------------')
        print('plot_web_decomp is done!')
        print('------------------------------------------------')

        return pd_plot_all
    
    # Added - Dec 23, 2021
    #####################################################################
    ### LAT-LON RadKernel
    #####################################################################
    def plot_RadKernel_latlon(self):
        print('plot_RadKernel_latlon starts ..........')
        
        cases_here = copy.deepcopy(self.cases[-1:])
        if 'amip-4xCO2' in self.cases:
            cases_here.remove('amip-4xCO2')
    
        variables = ['SWCRE_ano_grd_adj','LWCRE_ano_grd_adj','netCRE_ano_grd_adj']
        variables_out = ['SWCRE','LWCRE','netCRE']
        
        nlat = 73
        nlon = 144
        
        # E3SM
        data_all = np.zeros((nlat,nlon,len(variables),len(cases_here)))
    
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
        
        #=============================================================
        # start plotting ...
        #=============================================================
        pd_plot_all = pd.DataFrame(columns=['Variables','Description','Case.VS.Case','Plot'])

        # only plot the last case -- that is the one we care about now.
        for icase,case in enumerate(cases_here):
            #----------------------------------------------------------
            # define figures                         
            #----------------------------------------------------------
            nrow = 3; ncol = 1
            fig = plt.figure(figsize=(ncol*8,nrow*4)) # this creates and increases the figure size
    
            bounds = np.arange(-3,3.25,0.25)

            cmap = plt.cm.RdBu_r
            bounds2 = np.append(np.append(-500,bounds),500) # This is only needed for norm if colorbar is extended
            norm = mpl.colors.BoundaryNorm(bounds2, cmap.N) # make sure the colors vary linearly even if the desired color boundaries are at varied intervals

            case_out = case
    
            for ivar,svar in enumerate(variables):

                data_plot = [data_all[:,:,ivar,icase]]
                title_plot = [case_out]

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

                    PDF.make_colorbar(ax1, 'W/m$^2$/K', self.fh-5, im1, orientation='vertical')

                    ax1.set_title(title_plot[num0]+'\n'+variables_out[ivar]+' ['+str(np.round(avgDATA,2))+']',fontsize=self.fh)

                    num0 += 1
    
            fig.subplots_adjust(top=0.9)
    
            fig.savefig(self.figdir+'LatLon-adjusted-CRE_'+case_out+'.png',bbox_inches='tight',dpi=300)
            plt.close(fig)

            pd_plot = prepare_pd2html('../figure/LatLon-adjusted-CRE_'+case_out+'.png',
                                  'CRE [W/m2/K]',
                                  'Adjusted CRE feedback derived from radiative kernel method (SW, LW, and NET)',
                                  case_out)

            pd_plot_all = pd.merge(pd_plot_all, pd_plot, on =['Variables','Description','Case.VS.Case','Plot'],how='outer')

                
        print('------------------------------------------------')
        print('plot_RadKernel_latlon is done!')
        print('------------------------------------------------')

        return pd_plot_all
 
    
    #####################################################################
    ### LAT-LON cloud feedback based on cloud radiative kernel
    #####################################################################
    def plot_CldRadKernel_latlon(self):
    
        print('LatLon-Cloud-feedback-Decomposition starts ........')
    
        cases_here = copy.deepcopy(self.cases[-1:])
        if 'amip-4xCO2' in self.cases:
            cases_here.remove('amip-4xCO2')
    
        nlat = 90
        nlon = 144

        pd_plot_all = pd.DataFrame(columns=['Variables','Description','Case.VS.Case','Plot'])

        # --- define variables 
        sections = ['ALL','HI680','LO680']
        components = ['NET','SW','LW']
        decomps = ['tot','amt','alt','tau']

        sections_out = ['ALL', 'High cloud', 'Low cloud']
        components_out = ['NET','SW','LW']
        
        # generate figure based on case categories
        for icomp,component in enumerate(components):
            for isec,sec in enumerate(sections):
    
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
                # only plot the last case -- that is what we care about now.
                for icase,case in enumerate(cases_here):
   
                    fig=plt.figure(figsize=(18,12)) # this creates and increases the figure size
                    nrow = 3
                    ncol = 2
                    plt.suptitle(sec+' CTP bins ['+case+']',fontsize=self.fh,y=0.95)
                    bounds = np.arange(-3,3.25,0.25)
                    cmap = plt.cm.RdBu_r
                    bounds2 = np.append(np.append(-500,bounds),500) # This is only needed for norm if colorbar is extended
                    norm = mpl.colors.BoundaryNorm(bounds2, cmap.N) # make sure the colors vary linearly even if the desired color boundaries are at varied intervals
                    names = [component+'cld_tot',component+'cld_amt',component+'cld_alt',component+'cld_tau']
        
                    for n,name in enumerate(names):
                        DATA = data_all[:,:,n,icase]
       
                        DATA.setAxis(0,lats)
                        DATA.setAxis(1,lons)
        
                        ax1 = fig.add_subplot(nrow,ncol,n+1,projection=ccrs.Robinson(central_longitude=180.))
                        im1 = ax1.contourf(lons[:],lats[:],DATA,bounds,transform=ccrs.PlateCarree(),cmap=cmap,norm=norm,extend='both')
                        ax1.coastlines()
                        ax1.set_global()
        
                        avgDATA = avgdata[n,icase]

                        plt.title(name+' ['+str(np.round(avgDATA,3))+']',fontsize=self.fh)

                        cb = plt.colorbar(im1,orientation='vertical',drawedges=True,ticks=bounds[::2])
                        cb.set_label('W/m$^2$/K')
        
                    fig.subplots_adjust(top=0.9)
    
                    fig.savefig(self.figdir+'LatLon-Cloud-feedback-Decomposition_'+cases_here[icase]+'-'+component+'-'+sec+'.png',dpi=300,bbox_inches='tight')
                    plt.close(fig)

                    pd_plot = prepare_pd2html('../figure/LatLon-Cloud-feedback-Decomposition_'+cases_here[icase]+'-'+component+'-'+sec+'.png',
                                  component+'-'+sec+' [W/m2/K]',
                                  components_out[icomp]+' '+sections_out[isec]+' feedback',
                                  cases_here[icase])

                    pd_plot_all = pd.merge(pd_plot_all, pd_plot, on =['Variables','Description','Case.VS.Case','Plot'],how='outer')

        print('------------------------------------------------')
        print('plot_CldRadKernel_latlon is done!')
        print('------------------------------------------------')

        return pd_plot_all


    #####################################################################
    ### plot_CLOUD_profile
    #####################################################################
    def plot_CLOUD_profile(self):
        print('plot_CLOUD_profile starts ..........')
    
        cases_here = copy.deepcopy(self.cases)
        if 'amip-4xCO2' in self.cases:
            cases_here.remove('amip-4xCO2')
    
        #if len(cases_here) > 2:
        #    print('Hi, you have more than 2 cases to produce cloud profile. It will be very noisy...')
        #    print('Will skip this...')
        #    return

        var1 = ['CLOUD','CLDLIQ','CLDICE','T','Q','OMEGA']
        var1_tmp = ['cl', 'clw', 'cli','T','Q','OMEGA'] # for v1_coupled
        var1_range = [[-0.1,0.35],[-10,40],[-5,10],[190,310],[0,25],[-50,50]]
        var1_units = ['fraction','mg/kg','mg/kg','K','g/kg','hPa/day']
    
        # ===============================================================        
        # ===============================================================        
        # E3SM
        nlat = 73
        nlon = 144
        nlev = 72

        ## regions 
        #regions = ['90S-90N','30S-30N','30S-60S','30N-60N','Sc','Cu','As']
        #regions_bnds = [[-90,90], [-30,30], [-60,-30], [30,60],[-30,30],[-30,30],[-30,30]]
        regions = ['Sc','Cu','As']
        regions_bnds = [[-30,30],[-30,30],[-30,30]]

        lndocns = ['all','lnd','ocn']
        #lndocns = ['ocn']

        all_colors = ['tab:blue','tab:orange','tab:green','tab:red','tab:purple']
        exps = ['pi','ab','ano']
        exps_out = ['CTL','P4K','P4K-CTL[*4]']
     
        pd_plot_all = pd.DataFrame(columns=['Variables','Description','Case.VS.Case','Plot'])

        for lndocn in lndocns:
            for iregion,region in enumerate(regions):

                if lndocn in ['all','lnd'] and region in ['Sc','Cu','As']:
                    continue

                print('=======lndocn = ', lndocn, 'region = ',region,'=========')

                # -------------------- generate figures -------------------
                nrow = 2
                #ncol = len(var1)/nrow
                ncol = 3
                fig = plt.figure(figsize=(ncol*5,nrow*5))
                # -------------------- generate figures -------------------

                plt.suptitle(lndocn+': '+region,fontsize=12)

                for ivar,svar in enumerate(var1):

                    ax = fig.add_subplot(nrow,ncol,ivar+1)

                    for icase,case in enumerate(cases_here):

                        color = all_colors[icase]

                        # additional processing for Sc, Cu and Tropical Ascent regimes
                        if lndocn == 'ocn' and region in ['Sc','Cu','As']:
                            f1 = cdms.open(self.datadir_v2+'global_cloud_3regime_'+case+'.nc')
                            svar_in = var1[ivar]
                        else:
                            f1 = cdms.open(self.datadir_v2+'global_cloud_'+case+'.nc')
                            svar_in = var1[ivar]
                
                        for exp,exp_out in zip(exps,exps_out):

                            if lndocn == 'ocn' and region in ['Sc','Cu','As']:
                                if exp == 'ano':
                                    tmp_pi = f1(svar_in+'_pi_avg')
                                    tmp_ab = f1(svar_in+'_ab_avg')
                                    tmp1 = tmp_ab - tmp_pi
                                    tmp1.setAxisList(tmp_pi.getAxisList())
                                else:
                                    tmp1 = f1(svar_in+'_'+exp+'_avg')
                                
                                levs = tmp1.getLevel()[:]
                                
                                if region == 'Sc':
                                    ir = 0
                                elif region == 'Cu':
                                    ir = 1
                                elif region == 'As':
                                    ir = 2
                                data2 = MV.average(tmp1[ir,:,:],axis=0) #average over time -> (regime,level)

                            else:
                                tmp1 = f1(svar_in+'_'+exp+'_clim')

                                if exp == 'ano':
                                    tmp1 = tmp1 * 4. # scale it larger
                
                                levs = tmp1.getLevel()[:]
                                lats = tmp1.getLatitude()[:]
                                lons = tmp1.getLongitude()[:]
                
                                if lndocn == 'lnd':
                                    data1 = PDF.mask_land(lons,lats,tmp1,land=False)
                                    data1.setAxisList(tmp1.getAxisList())
                                elif lndocn == 'ocn':
                                    data1 = PDF.mask_land(lons,lats,tmp1,land=True)
                                    data1.setAxisList(tmp1.getAxisList())
                                else:
                                    data1 = tmp1 

                                # regime average 
                                data2 = cdutil.averager(data1.subRegion(lat=(regions_bnds[iregion][0],regions_bnds[iregion][1])), axis='xy',weights='weighted')

                            print(svar,case,'data2.shape=',data2.shape,np.min(data2),np.max(data2))
                            
                            # convert unit from kg/kg to mg/kg
                            if svar in ['CLDLIQ','CLDICE']:
                                data2 = data2 * 1e6
                            if svar in ['CLOUD'] and case == 'v1_coupled':
                                data2 = data2/100.
                            if svar in ['Q']:
                                data2 = data2 * 1e3
                            if svar in ['OMEGA']:
                                data2 = data2 * 864.
                            if svar in ['T'] and region in ['Sc','Cu']:
                                # convert into potentail temperature here
                                data2=data2*(1000/levs)**(287/1004)
   
                            #----------------------------------------------------------
                            # start plotting ...
                            #----------------------------------------------------------
                            unit = var1_units[ivar]
                            title = svar 

                            if exp in ['pi']:
                                ls = '-'
                            elif exp in ['ab']:
                                ls = '--'
                            elif exp in ['ano']:
                                ls = '-.'
                            
                            if exp == 'ano' and region in ['Sc','Cu','As']:
                                label = case.split('.')[-1]+': P4K-CTL'
                            else:
                                label=case.split('.')[-1]+': '+exp_out

                            im1 = ax.plot(np.array(data2),levs,label=label,c=color,ls=ls)
    
                            if region in ['Sc','Cu'] and svar == 'T':
                                ax.set_xlabel('$\Theta$ ['+unit+']')
                            else:
                                ax.set_xlabel(svar+' ['+unit+']')
                            ax.set_ylabel('Pressure [hPa]') 
    
                            ax.set_ylim(max(levs),min(levs))
                            if region in ['Sc','Cu']:
                                ax.set_ylim(max(levs),700)

                            ax.axvline(x=0,ls='-',color='grey',lw=0.5)
                            if region in ['Sc','Cu'] and svar == 'T':
                                ax.set_xlim((290,320)) # range for potential temp
                            else:
                                ax.set_xlim((var1_range[ivar][0],var1_range[ivar][1]))

                    ax.legend()

                fig.tight_layout()
                fig.savefig(self.figdir+'RegionalAvg_CloudProfile_'+lndocn+'_'+region+'.png',bbox_inches='tight',dpi=300)
                plt.close(fig)

                pd_plot = prepare_pd2html('../figure/RegionalAvg_CloudProfile_'+lndocn+'_'+region+'.png',
                              'regional average cloud profile',
                              lndocn+' '+region,
                              'none')

                pd_plot_all = pd.merge(pd_plot_all, pd_plot, on =['Variables','Description','Case.VS.Case','Plot'],how='outer')
                
        print('------------------------------------------------')
        print('plot_CLOUD_profile is done!')
        print('------------------------------------------------')

        return pd_plot_all
    
 
    #####################################################################
    ### plot_NRMSE
    #####################################################################
    def plot_NRMSE_RadKern(self):
        print('plot_NRMSE_RadKern starts ..........')
    
        cases_here = copy.deepcopy(self.cases)
        if 'amip-4xCO2' in self.cases:
            cases_here.remove('amip-4xCO2')
 
        df = pd.read_csv(self.datadir_v2+'COR_RMSE_RadKernel_latlon_dif_'+cases_here[-1]+'.csv',index_col=0)
        #print(df)

        df1 = df[(df['ref_case']=='F2010-p4Ka.v1') & (df['case']!='F2010-p4Ka.v1')]
        #print(df1)
        
        df21 = df1[df1['var']=='netCRE']
        df22 = df1[df1['var']=='SWCRE']
        df23 = df1[df1['var']=='LWCRE']
        
        print(df21)
        print(df21.index)
        print(df21.columns)

        pd_plot_all = pd.DataFrame(columns=['Variables','Description','Case.VS.Case','Plot'])

        #### plotting...
        fig = plt.figure(figsize=(18,12))
        nrow = 2
        ncol = 1

        xvalue = np.arange(len(df21['case']))
        print(xvalue)
        
        xticks = [case.split('.')[-1] for case in df21['case']]
        xticks = ['All' if item == 'gwenergy' else item for item in xticks]
        print(xticks)
        
        ii = 0
        for mm in ['COR','NRMSE']:

            ax = fig.add_subplot(nrow,ncol,ii+1)
            
            w = 0.10
            a1 = 0.5
            ww = w+0.02
            
            ax.bar(xvalue,df21[mm],width=w,alpha=a1,label='NET')
            ax.bar(xvalue-ww,df22[mm],width=w,alpha=a1,label='SW')
            ax.bar(xvalue+ww,df23[mm],width=w,alpha=a1,label='LW')
            
            ax.plot(xvalue,df21[mm],marker='o')
            ax.plot(xvalue-ww,df22[mm],marker='o')
            ax.plot(xvalue+ww,df23[mm],marker='o')
            
            ax.set_xticks(xvalue,xticks,fontsize=14)
            ax.tick_params(labelsize=14)
            ax.set_ylabel(mm+' [ref=v1]',fontsize=14)
            ax.legend(fontsize=12)
            ax.grid(axis='y')

            if mm == 'COR':
                ax.set_ylim((0.4,1.0))
            elif mm == 'NRMSE':
                ax.set_ylim((0.4,1.2))
        
            ii = ii + 1
        
        fig.savefig(self.figdir+'plot_NRMSE_COR_RadKern_'+cases_here[-1]+'.png',bbox_inches='tight',dpi=300)

        pd_plot = prepare_pd2html('../figure/plot_NRMSE_COR_RadKern_'+cases_here[-1]+'.png',
                      'spatial correlation (COR) and NRMSE for RadKern',
                      'bar plot',
                      'none')

        pd_plot_all = pd.merge(pd_plot_all, pd_plot, on =['Variables','Description','Case.VS.Case','Plot'],how='outer')
 
        return pd_plot_all

    #####################################################################
    ### cal regional LAT-LON RadKernel correspondence b/t two exps
    ### NOTE: this function is only used for ref_case number == 1 now. 
    #####################################################################
    def cal_RadKernel_regional_correspondence(self):
        print('cal_RadKernel_regional_correspondence starts ..........')
        
        cases_here = copy.deepcopy(self.cases)
        if 'amip-4xCO2' in self.cases:
            cases_here.remove('amip-4xCO2')
    
        variables = ['SWCRE_ano_grd_adj','LWCRE_ano_grd_adj','netCRE_ano_grd_adj']
        variables_out = ['SWCRE','LWCRE','netCRE']
        colors = ['tab:blue','tab:green','tab:orange']
        
        nlat = 73
        nlon = 144
        
        #=============================================================
        # read ALL data in 
        #=============================================================
        data_all = np.zeros((nlat,nlon,len(variables),len(cases_here)))
    
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
    
                data_all[:,:,ivar,icase] = data
        
        #=============================================================
        # calculate REGIONAL NRMSE and spatial correlation
        #=============================================================
        regions = ['90S-90N','30S-60S','30N-60N','WPacific','EPacific','Atlantic']
        regions_bnds = [[-90,90],[-60,-30],[30,60], [-30,30],[-30,30],[-30,30]]
        regions_lonbnds = [[0,360],[0,360],[0,360], [30,180],[185,290],[-70,30]]

        #regions = ['WPacific','EPacific','Atlantic']
        #regions_bnds = [[-30,30],[-30,30],[-30,30]]
        #regions_lonbnds = [[30,180],[185,290],[-70,30]]

        #lndocns = ['all','lnd','ocn']
        lndocns = ['ocn','lnd']
        #lndocns = ['lnd']

        pd_plot_all = pd.DataFrame(columns=['Variables','Description','Case.VS.Case','Plot'])

        for lndocn in lndocns:

            ## save bar plot: fig1 -- correlation; fig2 -- NRMSE
            fig1 = plt.figure(figsize=(18,12))
            fig2 = plt.figure(figsize=(18,12))
            ii = 0

            for iregion,region in enumerate(regions):

                # only calculate global land rather than regional distribution 
                if lndocn == 'lnd' and region != '90S-90N':
                    continue

                ax2 = fig1.add_subplot(2,3,ii+1)
                ax3 = fig2.add_subplot(2,3,ii+1)
                ii = ii + 1

                ## plot spatial distribution
                fig = plt.figure(figsize=(18,21))
                nrow = len(cases_here)
                ncol = 3

                for ivar,svar in enumerate(variables):

                    dics_cor = {}
                    dics_rmse = {}
                    for icase,case in enumerate(cases_here):
                        ref_cases = self.ref_casesA[icase]
                        if len(ref_cases) == 0:
                            continue

                        if len(ref_cases) > 1:
                            print('Warning: we only need one ref_case for this plot. You have ',len(ref_cases),'Exiting...')
                            return 

                        for ref_case in ref_cases:
                            iref = cases_here.index(ref_case)

                            # calculate COR and NRMSE
                            daa = data_all[:,:,ivar,icase]
                            dbb = data_all[:,:,ivar,iref]
    
                            daa = cdms.asVariable(daa)
                            dbb = cdms.asVariable(dbb)
                            daa.setAxis(0,lats)
                            daa.setAxis(1,lons)
                            dbb.setAxis(0,lats)
                            dbb.setAxis(1,lons)
    
                            if lndocn == 'lnd':
                                land_mask = False
                                daa1 = PDF.mask_land(lons[:],lats[:],daa,land=land_mask)
                                dbb1 = PDF.mask_land(lons[:],lats[:],dbb,land=land_mask)
                                daa1.setAxisList(daa.getAxisList())
                                dbb1.setAxisList(dbb.getAxisList())
                            elif lndocn == 'ocn':
                                land_mask = True
                                daa1 = PDF.mask_land(lons[:],lats[:],daa,land=land_mask)
                                dbb1 = PDF.mask_land(lons[:],lats[:],dbb,land=land_mask)
                                daa1.setAxisList(daa.getAxisList())
                                dbb1.setAxisList(dbb.getAxisList())
                            else:
                                daa1 = daa
                                dbb1 = dbb
    
                            if region == 'Atlantic':
                                pivotLon = 290.
                                daa2,tlon = LP.lonPivot(daa1, lons[:], pivotLon)
                                dbb2,_ = LP.lonPivot(dbb1, lons[:], pivotLon)
        
                                tlonc = cdms.createAxis(tlon)
                                tlonc.id="lon"
                                tlonc.units="degrees_E"
                                tlonc.designateLongitude()
        
                                daa2.setAxis(0,lats)
                                dbb2.setAxis(0,lats)
                                daa2.setAxis(1,tlonc)
                                dbb2.setAxis(1,tlonc)
                            else:
                                daa2 = daa1
                                dbb2 = dbb1
    
    
                            daa1_reg = daa2.subRegion(lat=(regions_bnds[iregion][0],regions_bnds[iregion][1]),lon=(regions_lonbnds[iregion][0],regions_lonbnds[iregion][1]))
                            dbb1_reg = dbb2.subRegion(lat=(regions_bnds[iregion][0],regions_bnds[iregion][1]),lon=(regions_lonbnds[iregion][0],regions_lonbnds[iregion][1]))
    
                            wts = np.cos(np.deg2rad(daa1_reg.getLatitude()[:]))
                            cor,NRMSE, RMSE = PDF.pattern_cor(dbb1_reg,daa1_reg,wts,1)
                            print('lndocn=',lndocn,'region=',region,'case=',case,'ref_case=',ref_case,'cor=',cor, 'NRMSE=',NRMSE)

                            dics_cor[case] = cor
                            dics_rmse[case] = NRMSE

                            #### LAT-LON plot
                            levels = np.arange(-3,3.5,0.5)
                            cmap = 'RdBu_r'
                            lon_here = daa1_reg.getLongitude()[:]
                            lat_here = daa1_reg.getLatitude()[:]
                            ax1 = fig.add_subplot(nrow,ncol,icase*ncol+ivar+1,projection=ccrs.PlateCarree())
#                            ax1.contourf(lon_here,lat_here, daa1_reg,transform=ccrs.PlateCarree(),levels=levels,extend='both',cmap=cmap)
                            ax1.pcolormesh(lon_here,lat_here,daa1_reg, transform=ccrs.PlateCarree(),vmin=-3,vmax=3,cmap=cmap)
                            ### !!! set_extent really needs the lon ranging from -180 to 180!!!
                            lon_bounds = lon_here
                            lon_bounds[lon_bounds>180] -= 360
                            print(lon_bounds)
                            ax1.set_extent([min(lon_bounds),max(lon_bounds),min(lat_here),max(lat_here)], crs=ccrs.PlateCarree())
                            ax1.coastlines()
                            ax1.set_title(svar.split('_')[0]+' '+case.split('.')[-1]+' R='+str(cor.round(2))+', NRMSE='+str(NRMSE.round(2)))

                    
                    #### bar plot 
                    w = 0.25
                    a1 = 0.5
                    ww = w+0.02

                    xvalue = np.arange(len(dics_cor.keys()))
                    xticks = [case.split('.')[-1] for case in dics_cor.keys()]
                    xticks = ['All' if item == 'gwenergy' else item for item in xticks]
                    print(xticks)

                    label = svar.split('_')[0]

                    #---- correlation 
                    datap = [dics_cor[case] for case in dics_cor.keys()]
                    ax2.bar(xvalue+ivar*ww,datap, width=w,alpha=a1,label=label,color=colors[ivar])
                    ax2.plot(xvalue+ivar*ww,datap ,marker='o',color=colors[ivar])
                    
                    #---- NRMSE
                    datap = [dics_rmse[case] for case in dics_cor.keys()]
                    ax3.bar(xvalue+ivar*ww,datap, width=w,alpha=a1,label=label,color=colors[ivar])
                    ax3.plot(xvalue+ivar*ww,datap ,marker='o',color=colors[ivar])

                    for axx in [ax2,ax3]:
                        axx.set_xticks(xvalue,xticks,fontsize=14)
                        axx.tick_params(labelsize=12)
                        axx.legend()
    
                        axx.set_title(region+' '+lndocn,fontsize=15)

                    ax2.set_ylabel('spatial correlation',fontsize=14)
                    ax2.set_ylim((-0.1,1))

                    ax3.set_ylabel('NRMSE',fontsize=14)
                    ax3.set_ylim((0,1.6))

                fig.savefig(self.figdir+'plot_regional_cor_'+lndocn+'_'+region+'.png',dpi=300,bbox_inches='tight')
                pd_plot = prepare_pd2html('../figure/'+'plot_regional_cor_'+lndocn+'_'+region+'.png',
                              'regional RadKern cloud feedbacks (all are on the same page)',
                              region+' '+lndocn, 
                              case+' vs '+ref_case)
                pd_plot_all = pd.merge(pd_plot_all, pd_plot, on =['Variables','Description','Case.VS.Case','Plot'],how='outer')

            fig1.savefig(self.figdir+'plot_regional_cor_barplot_'+lndocn+'_'+cases_here[-1]+'.png',dpi=300,bbox_inches='tight')
            pd_plot = prepare_pd2html('../figure/'+'plot_regional_cor_barplot_'+lndocn+'_'+cases_here[-1]+'.png',
                          'spatial correlation (COR) for RadKern cloud feedbacks',
                          'bar plot: '+lndocn,
                          case+' vs '+ref_case)
            pd_plot_all = pd.merge(pd_plot_all, pd_plot, on =['Variables','Description','Case.VS.Case','Plot'],how='outer')

            fig2.savefig(self.figdir+'plot_regional_nrmse_barplot_'+lndocn+'_'+cases_here[-1]+'.png',dpi=300,bbox_inches='tight')
            pd_plot = prepare_pd2html('../figure/'+'plot_regional_nrmse_barplot_'+lndocn+'_'+cases_here[-1]+'.png',
                          'NRMSE for RadKern cloud feedbacks',
                          'bar plot: '+lndocn,
                          case+' vs '+ref_case)
            pd_plot_all = pd.merge(pd_plot_all, pd_plot, on =['Variables','Description','Case.VS.Case','Plot'],how='outer')

        print(pd_plot_all)
 
        return pd_plot_all

        print('------------------------------------------------')
        print('cal_RadKernel_regional_correspondence is done!')
        print('------------------------------------------------')
    
