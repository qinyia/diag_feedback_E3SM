
# Aug 31, 2020: This py file is used to save all functions for plotting. 

import cdms2 as cdms
import cdutil
import MV2 as MV
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import os
import pandas as pd
import cdtime
import time
import genutil
from genutil import statistics
from scipy import stats
import numpy.ma as ma
import scipy as sp
import scipy.interpolate
from global_land_mask import globe

from sklearn import datasets, linear_model
import statsmodels.api as sm
import statsmodels.formula.api as smf
from sklearn.metrics import r2_score

from matplotlib import ticker

#========================================================================================================
def find_model_list(datadir,prefix,suffix1,suffix2,phases,exp_cntl,exp_new):
    '''
    get intersection between coupled and amip model lists
    '''
    # coupled
    for iphase,phase in enumerate(phases):
        if phase == 'CMIP5':
            stream = os.popen('ls '+datadir+prefix+'_'+phase+'_'+exp_new[iphase][0]+suffix1)
        else:
            stream = os.popen('ls '+datadir+prefix+'_'+phase+'_'+exp_new[iphase][0]+suffix1)
        output = stream.readlines()
        # separate by '_'
        newlist_1 = []
        for x in output:
            x1 = x.split('_')
            idx_exp = x1.index(exp_new[iphase][0])
            x2 = x1[idx_exp+1].split('.')
            newlist_1.append(x2[0])

        if phase == 'CMIP5':
            cmip5_cmip_models = newlist_1
        else:
            cmip6_cmip_models = newlist_1

    # amip
    for iphase,phase in enumerate(phases):
        if phase == 'CMIP5':
            stream = os.popen('ls '+datadir+prefix+'_'+phase+'_'+exp_new[iphase][1]+suffix2)
        else:
            stream = os.popen('ls '+datadir+prefix+'_'+phase+'_'+exp_new[iphase][1]+suffix2)
        output = stream.readlines()
        # separate by '_'
        newlist_1 = []
        for x in output:
            x1 = x.split('_')
            idx_exp = x1.index(exp_new[iphase][1])
            x2 = x1[idx_exp+1].split('.')
            newlist_1.append(x2[0])
        if phase == 'CMIP5':
            cmip5_amip_models = newlist_1
        else:
            cmip6_amip_models = newlist_1

#    print(cmip5_cmip_models)
#    print(cmip5_amip_models)
    
    # part 3: get the intersection b/t cmip and amip experiments in both CMIP5 and CMIP6
    cmip5_models = [ x for x in cmip5_amip_models if x in cmip5_cmip_models ]
    cmip6_models = [ x for x in cmip6_amip_models if x in cmip6_cmip_models ]

    for xx in cmip5_cmip_models:
        if xx in ['CanESM2','HadGEM2-ES'] and xx not in cmip5_models:
            cmip5_models.append(xx)

    return cmip5_models, cmip6_models

# ==========================================================================================================================
def find_model_ripf_list(datadir,prefix,suffix1,suffix2,phases,exp_cntl,exp_new):
    '''
    find model list and related ripf list
    '''
    # coupled
    for iphase,phase in enumerate(phases):
        if phase == 'CMIP5':
            stream = os.popen('ls '+datadir+prefix+'_'+phase+'_'+exp_new[iphase][0]+suffix1)
        else:
            stream = os.popen('ls '+datadir+prefix+'_'+phase+'_'+exp_new[iphase][0]+suffix1)
        output = stream.readlines()
        # separate by '_'
        newlist_1 = []
        newlist_2 = []
        for x in output:
            x1 = x.split('_')
            idx_exp = x1.index(exp_new[iphase][0])
            x2 = x1[idx_exp+1].split('.')
            newlist_1.append(x2[0])
            newlist_2.append(x1[idx_exp+2].split('.')[0])

        if phase == 'CMIP5':
            cmip5_cmip_models = newlist_1
            cmip5_cmip_ripfs = newlist_2
        else:
            cmip6_cmip_models = newlist_1
            cmip6_cmip_ripfs = newlist_2

    # amip
    for iphase,phase in enumerate(phases):
        if phase == 'CMIP5':
            stream = os.popen('ls '+datadir+prefix+'_'+phase+'_'+exp_new[iphase][1]+suffix2)
        else:
            stream = os.popen('ls '+datadir+prefix+'_'+phase+'_'+exp_new[iphase][1]+suffix2)
        output = stream.readlines()
        # separate by '_'
        newlist_1 = []
        newlist_2 = []
        for x in output:
            x1 = x.split('_')
            idx_exp = x1.index(exp_new[iphase][1])
            x2 = x1[idx_exp+1].split('.')
            newlist_1.append(x2[0])
            newlist_2.append(x1[idx_exp+2].split('.')[0])
        if phase == 'CMIP5':
            cmip5_amip_models = newlist_1
            cmip5_amip_ripfs = newlist_2
        else:
            cmip6_amip_models = newlist_1
            cmip6_amip_ripfs = newlist_2

    # part 3: get the intersection b/t cmip and amip experiments in both CMIP5 and CMIP6
    cmip5_models = [ x for x in cmip5_amip_models if x in cmip5_cmip_models ]
    cmip6_models = [ x for x in cmip6_amip_models if x in cmip6_cmip_models ]
    
    id1 = [ ix for ix,x in enumerate(cmip5_amip_models) if x in cmip5_cmip_models ]
    id2 = [ ix for ix,x in enumerate(cmip6_amip_models) if x in cmip6_cmip_models ]
    
    
    cmip5_ripfs = [j for ij,j in enumerate(cmip5_amip_ripfs) if ij in id1]
    cmip6_ripfs = [j for ij,j in enumerate(cmip6_amip_ripfs) if ij in id2]
    
#     print(len(id1),len(id2),len(cmip5_models),len(cmip6_models),len(cmip5_ripfs),len(cmip6_ripfs))

    for xx in cmip5_cmip_models:
        if xx in ['CanESM2','HadGEM2-ES','CanAM4','HadGEM2-A'] and xx not in cmip5_models:
            if (xx == 'CanESM2') and ('CanAM4' in cmip5_amip_models or 'CanESM2' in cmip5_amip_models):
                cmip5_models.append(xx)
                cmip5_ripfs.append('r1i1p1')
            elif (xx == 'HadGEM2-ES') and ('HadGEM2-A' in cmip5_amip_models or 'HadGEM2-ES' in cmip5_amip_models):
                cmip5_models.append(xx)
                cmip5_ripfs.append('r1i1p1')

    return cmip5_models, cmip6_models,cmip5_ripfs,cmip6_ripfs


#========================================================================================================
def get_color(colormap,color_nums):
    '''
    get color based on user-defined colormap and number or colors needed.
    '''
    palette = plt.get_cmap(colormap)
    colors = []
    for ic in range(color_nums):
        color = palette(ic)
        colors.append(color)
    return colors

#========================================================================================================
def reg_and_cor(XX,YY):
    '''
    goal: get regression and correlation
    input: XX and YY
    output: 
    xmin, xmax, ymin, ymax --- XXmin, XXmax, predicted YYmin, predicted YYmax
    cor -- correlation 
    pval -- pvalue
    r2 -- r2_score. Slightly different from the general casde. Use XX as YY_pred
    '''
    
    # convert from pandas to array; eleminate the additional column number 
    XX = np.array(XX.iloc[:,0])
    YY = np.array(YY.iloc[:,0])
    
    # regression 
    reg1 = linear_model.LinearRegression().fit(XX.reshape(-1,1),YY.reshape(-1,1))
    YY_pred = reg1.predict(XX.reshape(-1,1))
    reg1_coef = reg1.coef_
    reg1_intercept = reg1.intercept_

    xmin = min(XX)
    xmax = max(XX)
    ymin = xmin * reg1_coef[0,0] + reg1_intercept[0]
    ymax = xmax * reg1_coef[0,0] + reg1_intercept[0]
    
    # correlation 
#    cor,pval = scipy.stats.pearsonr(XX,YY)
    cor,pval = scipy.stats.spearmanr(XX,YY)

    # r2 score --- this one is used to estimate whether they are on 1v1 line, 
    # not the same as the traditional R2 from regression
#    r2 = r2_score(XX,YY)
    # r2_score(yy_true,yy_pred)
    r2 = r2_score(YY,XX)
    
    return xmin,xmax,ymin,ymax,cor,pval,r2

#========================================================================================================
def styling_specific_cell(x,input_ind):
    '''
    highlight some values in the table when some criteria is reached
    '''
    df_styler = pd.DataFrame('', index=x.index, columns=x.columns)
    for i,j in input_ind:
        print('in function',i,j)
#         color = 'background-color: white; color: deepskyblue'
        color = 'color: deepskyblue'
        df_styler.loc[i, j] = color
    return df_styler


#========================================================================================================
def reg_and_cor_coef_intercept(XX,YY):
    '''
    similar to reg_and_cor(XX,YY). Just add two additional outputs: regression coefficient and intercept.
    '''

    # convert from pandas to array; eleminate the additional column number
    XX = np.array(XX.iloc[:,0])
    YY = np.array(YY.iloc[:,0])

    # regression
    reg1 = linear_model.LinearRegression().fit(XX.reshape(-1,1),YY.reshape(-1,1))
    YY_pred = reg1.predict(XX.reshape(-1,1))
    reg1_coef = reg1.coef_
    reg1_intercept = reg1.intercept_

    xmin = min(XX)
    xmax = max(XX)
    ymin = xmin * reg1_coef[0,0] + reg1_intercept[0]
    ymax = xmax * reg1_coef[0,0] + reg1_intercept[0]

    coef = reg1_coef[0,0]
    intercept = reg1_intercept[0]
    
    # correlation
    cor,pval = scipy.stats.pearsonr(XX,YY)
#     cor,pval = scipy.stats.spearmanr(XX,YY)

    # r2 score --- this one is used to estimate whether they are on 1v1 line,
    # not the same as the traditional R2 from regression
#    r2 = r2_score(XX,YY)
    # r2_score(yy_true,yy_pred)
    r2 = r2_score(YY,XX)

    return xmin,xmax,ymin,ymax,cor,pval,r2,coef,intercept

#========================================================================================================
def pd_find_nan_index(df):
    '''
    find index when value is NaN.
    '''
    is_NaN = df.isnull()
    row_has_NaN = is_NaN.any(axis=1)
    rows_with_NaN = df[row_has_NaN]
    rows_idx_with_NaN = list(rows_with_NaN.index)

    return rows_idx_with_NaN

#========================================================================================================
def sm_find_outlier(xnew,ynew,models_all):
    
    # --------- part 1: clear raw data ----------------------------------------
    # first figure out whether there is NaN in raw data.
    # if so, find their index and remove them from both xnew and ynew.
    list_1 = pd_find_nan_index(xnew)
    list_2 = pd_find_nan_index(ynew)
    # merge the NaN index list into a new list -- combined_list
    list_2_items_not_in_list_1 = list(set(list_2) - set(list_1))
    combined_list = list_1 + list_2_items_not_in_list_1
    # drop those NaN index from both xnew and ynew
    xnew1 = xnew.drop(index=combined_list)
    ynew1 = ynew.drop(index=combined_list)
    # drop those NaN index also from index name list: models_all
    models_all_1 = [i for j, i in enumerate(models_all) if i not in combined_list]
#     print(models_all_1)

    # ------- part 2: OLS and influence analysis ------------------------------
    # start get OLS and influence plots
    lm1 = sm.OLS(np.array(xnew1.iloc[:,0]), np.array(ynew1.iloc[:,0])).fit()
    rsquared = lm1.rsquared
    
    infl= lm1.get_influence()
    xx = infl.summary_frame().filter(["hat_diag","student_resid","dffits","cooks_d"])
    xx.index = models_all_1
    n = xnew1.shape[0] # sample size
    p = 1 # number of predictant variables
    D_lim = 4/(n-p-1) # threshold for cook's influence. if greater, this point will be treated as a potential outlier.
    outlier_idx = list(xx[xx['cooks_d'] > D_lim].index)

#     figx,axx = plt.subplots(figsize=(12,8))
#     sm.graphics.influence_plot(lm1, ax=axx, criterion="cooks")

    return xnew1,ynew1,outlier_idx,rsquared

# =================================================================================
def reg_and_cor_nan(models_all,data1,data2,do_kick_outlier): 
    xnew = pd.DataFrame(data1)
    ynew = pd.DataFrame(data2)
    
    xnew1,ynew1,outlier_idx,rsquared = sm_find_outlier(xnew,ynew,models_all)
    
    if do_kick_outlier:
        if len(outlier_idx) > 0:
            print('Hi~ we will drop these models: ',outlier_idx)
            xnew2 = xnew1.drop(index=outlier_idx)
            ynew2 = ynew1.drop(index=outlier_idx)
        else:
            xnew2 = xnew1
            ynew2 = ynew1
    else:
        xnew2 = xnew1
        ynew2 = ynew1
    
#     print('-----------')
#     print(len(xnew2.index),' models are used for regression and correlation. They are: ',xnew2.index)
#     print(len(xnew2.index),' models are used for regression and correlation.')    
    xmin2,xmax2,ymin2,ymax2,cor2,val2,r22,coef,intercept = reg_and_cor_coef_intercept(xnew2,ynew2)
    
    return xmin2,xmax2,ymin2,ymax2,cor2,val2,r22,coef,intercept

# ==========================================================================================================================
def add_regline_cortext(ax0,xnew1,ynew1,model_list,do_kick_outlier,xpos,ypos,c1,ls1,lw1,fh,add_cor=True,add_r2=False,add_regtext=False,
                       suffix='',add_regline=False):
    '''
    add regression lines and correlation as text based on pvalues to add * additionally
    Oct 29, 2020: the default option is to add correlation coefficient. 
    if adding regression equation, adding correlation coefficient is turned off.
    '''
    xmin1,xmax1,ymin1,ymax1,cor1,pval1,r1,coef,intercept = reg_and_cor_nan(model_list,xnew1,ynew1,do_kick_outlier=False)

    if add_regline:
        ax0.plot([xmin1,xmax1],[ymin1,ymax1],lw=lw1,color=c1,ls=ls1)
    
    if add_regtext:
        add_cor = False
        
    if add_cor:
        if add_r2:
            if pval1 < 0.05:
                outstring = suffix+' cor='+str(np.round(cor1,2))+'$^*$, r2='+str(np.round(r1,2))
            else:
                outstring = suffix+' cor='+str(np.round(cor1,2))+', r2='+str(np.round(r1,2))
            
        else:
            if pval1 < 0.05:
                outstring = suffix+' cor='+str(np.round(cor1,2))+'$^*$'
            else:
                outstring = suffix+' cor='+str(np.round(cor1,2))
            
        ax0.text(xpos, ypos,outstring,horizontalalignment='left',verticalalignment='center',
             transform = ax0.transAxes,fontsize=fh,color=c1)
    
    if add_regtext:
        ax0.text(xpos, ypos,'Y='+str(np.round(coef,2))+'*X $+$ '+str(np.round(intercept,2)),horizontalalignment='left',verticalalignment='center',
                 transform = ax0.transAxes,fontsize=fh,color=c1)

# ==========================================================================================================================

def get_intersect(exp_cntl,exp_new,prefix,suffix1,suffix2,datadir):
    '''
    get model list for specific experiment
    '''
    phases = ['CMIP5','CMIP6']
    
    cmip5_models,cmip6_models = \
    find_model_list(datadir,prefix,suffix1,suffix2,phases,exp_cntl,exp_new)
   
    models_all = cmip5_models + cmip6_models

    return models_all,cmip5_models,cmip6_models

# ==========================================================================================================================

def get_intersect_withripf(exp_cntl,exp_new,prefix,suffix1,suffix2,datadir):
    '''
    get model list for specific experiment
    '''
    phases = ['CMIP5','CMIP6']
    
    cmip5_models,cmip6_models,cmip5_ripfs,cmip6_ripfs = \
    find_model_ripf_list(datadir,prefix,suffix1,suffix2,phases,exp_cntl,exp_new)
   
    models_all = cmip5_models + cmip6_models
    ripfs_all = cmip5_ripfs + cmip6_ripfs
    
    models_ripfs_all = [m+'_'+n for m,n in zip(models_all,ripfs_all)]
    cmip5_models_ripfs = [m+'_'+n for m,n in zip(cmip5_models,cmip5_ripfs)]
    cmip6_models_ripfs = [m+'_'+n for m,n in zip(cmip6_models,cmip6_ripfs)]

    return models_ripfs_all,cmip5_models_ripfs,cmip6_models_ripfs


# ==========================================================================================================================
def make_colorbar(ax, units, fh, mappable,nbins=11, **kwargs):
    '''
    Nov 23, 2020: a function to add colorbar fitted well with the figure
    Referred from: https://github.com/pydata/xarray/issues/619
    units: the unit of colorbar
    fh: fontsize of colorbar label font size
    '''
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    import matplotlib as mpl

    divider = make_axes_locatable(ax)
    orientation = kwargs.pop('orientation', 'vertical')
    if orientation == 'vertical':
        loc = 'right'
    elif orientation == 'horizontal':
        loc = 'bottom'
        
    cax = divider.append_axes(loc, '5%', pad='5%', axes_class=mpl.pyplot.Axes)
#    cb = ax.get_figure().colorbar(mappable, cax=cax, orientation=orientation,extend='both')
    cb = ax.get_figure().colorbar(mappable, cax=cax, orientation=orientation)

    cb.set_label(units,fontsize=fh)
    cb.ax.tick_params(labelsize=fh)

    tick_locator = ticker.MaxNLocator(nbins=nbins)
    cb.locator = tick_locator
    cb.update_ticks()
# ==========================================================================================================================


def add_common_colorbar(fig,im,axes,units,orientation='vertical',nbins=9,fontsize=15):
    '''
    Feb 12, 2021: generate a common colorbar for several specific subplots.
    For example, subplots on the same row or the same column.
    -----------------------
    Inputs:
    fig -- 
    im -- return handle. 
    axes --- the list for all subplots.
    units --- unit of colorbar.
    orientation --- direction of colorbar.
    nbins --- number of labeled ticks. default: 9
    fontsize --- label fontsize for colorbar. 
    '''
    pos1 = axes[-1].get_position() # get the original position for the last subplot 
    if orientation == 'vertical':
        pos2 = [pos1.x0 + pos1.width + 0.01, pos1.y0,  pos1.width / 20.0, pos1.height ]   
    else:
        pos2 = [pos1.x0, pos1.y0 - 0.10, pos1.width, pos1.height / 25.0]
        
    cbar_ax = fig.add_axes(pos2)
    cb = fig.colorbar(im,ax=axes, orientation=orientation, cax=cbar_ax)
    cb.set_label(units,fontsize=fontsize)
    cb.ax.tick_params(labelsize=fontsize)
    tick_locator = ticker.MaxNLocator(nbins=nbins)
    cb.locator = tick_locator
    cb.update_ticks()

# ==========================================================================================================================

def create_nested_dict_ND(n, type):
    '''
    Jan 27, 2021:
    create a nested dictionary with number of dimensions and type.
    inputs:
    n -- number of dimension for new dictionary
    type -- type of each element, like float, str
    '''
    from collections import defaultdict
    if n == 1:
        return defaultdict(type)
    else:
        return defaultdict(lambda: create_nested_dict_ND(n-1, type))

# ==========================================================================================================================

def create_nested_dict():
    ''' 
    Jan 27, 2021:
    create a nested dictionary with no limit on dimensions and type.
    ref: https://stackoverflow.com/questions/5369723/multi-level-defaultdict-with-variable-depth/8702435#8702435
    '''
    from collections import defaultdict
    return defaultdict(lambda: create_nested_dict())


#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
def betai(a, b, x):
    import scipy.special as special
    x = np.asarray(x)
    x = np.where(x < 1.0, x, 1.0)  # if x > 1 then return 1.0
    return special.betainc(a, b, x)

def _chk_asarray(a, axis):
    if axis is None:
        a = np.ravel(a)
        outaxis = 0
    else:
        a = np.asarray(a)
        outaxis = axis
    return a, outaxis

def ss(a, axis=0):
    """
    Squares each element of the input array, and returns the sum(s) of that.
    ref: https://github.com/scipy/scipy/blob/f2ec91c4908f9d67b5445fbfacce7f47518b35d1/scipy/stats/stats.py#L4274
    Parameters
    ----------
    """
    a, axis = _chk_asarray(a, axis)
    return np.sum(a*a, axis)

def pearsonr_nd(x, y):
    """
    Parameters
    ----------
    x : (N,,,) array_like
        Input
    y : (N,,,) array_like
        Input
    #<qinyi 2021-02-12 #------------------
    Description: revised based on 1D pearsonr function. Please keep the correlation axis at the most left.
    #>qinyi 2021-02-12 #------------------
    Returns
    -------
    (Pearson's correlation coefficient,
     2-tailed p-value)
    References
    ----------
    http://www.statsoft.com/textbook/glosp.html#Pearson%20Correlation
    """
    # x and y should have same length.
    x = np.asarray(x)
    y = np.asarray(y)
    n = x.shape[0]
    mx = np.mean(x,axis=0)
    my = np.mean(y,axis=0)
    
    if len(x.shape) == 1:
        newshape = n
    if len(x.shape) == 2:
        newshape = np.append(n,1)
    elif len(x.shape) == 3:
        newshape = np.append(n,[1,1])
        
    mx_new = np.tile(mx,newshape)
    my_new = np.tile(my,newshape)
    
    xm, ym = x-mx_new, y-my_new    
    r_num = np.add.reduce(xm * ym,axis=0)
    r_den = np.sqrt(ss(xm,axis=0) * ss(ym,axis=0))
    r = r_num / r_den

    # Presumably, if abs(r) > 1, then it is only some small artifact of floating
    # point arithmetic.
    r = np.where(r>1.0, 1.0, r)
    r = np.where(r<-1.0, -1.0, r)

    df = n-2
    t_squared = r*r * (df / ((1.0 - r) * (1.0 + r)))
    prob = betai(0.5*df, 0.5, df / (df + t_squared))    
    prob = np.where(r==1.0,0.0,prob)
    
    return r, prob

# ===================================================================================
def get_scaled_lat(lats, spec_lats):
    '''
    Feb 8, 2021: get scaled latitude values based on number of lats and specified latitudes. 
    '''
    clat = np.cos(np.deg2rad(lats))
    clat1 = clat/MV.sum(clat)
    clat1[0] = 0.
    
    clats = np.zeros(len(clat1))
    
    for ilat in range(len(clat1)):
        clats[ilat] = np.sum(clat1[:ilat+1])
    
    clats[0] = 0.
    
    N = [i for i in range(len(lats)) if lats[i] in spec_lats]
    spec_clats = list(np.array(clats)[N])

    return clats,spec_clats

# =========================================================================================
def remove_F2010(case):
    # eleminate F2010 for output title to shorten the output case name 
    if 'F2010' in case:
        if len(case.split('-')) > 2:
            case_out = '-'.join(case.split('-')[1:])
        else:
            case_out = case.split('-')[1]
    else:
        case_out = case

    return case_out

# =========================================================================================