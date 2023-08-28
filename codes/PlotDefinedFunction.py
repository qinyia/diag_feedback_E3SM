
# Aug 31, 2020: This py file is used to save all functions for plotting. 

import numpy as np
import numpy 
import os
import pandas as pd
import time
from scipy import stats
import numpy.ma as ma
import scipy as sp
import scipy.interpolate
from global_land_mask import globe

from sklearn import datasets, linear_model
from sklearn.metrics import r2_score


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
    import matplotlib.pyplot as plt 

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

    #tick_locator = ticker.MaxNLocator(nbins=nbins)
    #cb.locator = tick_locator
    #cb.update_ticks()

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

    from matplotlib import ticker

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
    clat1 = clat/np.ma.sum(clat)
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

def pattern_cor(x,y,wts,opt):
    '''
    Revised from contributed.ncl by Yi Qin
    This function is used to compute the pattern correlation between two fileds.
    (lat,lon),...,wgt(lat)

    also, calculate normalized root mean square error.
    NRMSE = sqrt ( sum ( (X_obs - X_simulate) ** 2) / sum (X_obs**2) )

    added RMSE output

    note: input x should be the 'truth', for NRMSE calculation.

    '''
    import numpy as np

    x = np.ma.array(x)
    y = np.ma.array(y)

    # if all are missing...
    if np.all(x.mask):
        rFill = np.nan
        return rFill

    if np.all(y.mask):
        rFill = np.nan
        return rFill

    dimx = x.shape
    dimy = y.shape

    if dimx != dimy:
        print("pattern_cor: Fatal: x and y do not have the same dimension sizes")
        print("dimx: ", dimx, "dimy: ", dimy)
        exit()

    dimw = wts.shape
    rankw = len(dimw)

    if rankw > 2:
        print("pattern_cor: Fatal: wts can be a scalar, w[*] or w[*][*]")
        print("rankw: ", rankw)
        exit()

    if rankw == 2 and dimx != dimw:
        print("pattern_cor: Fatal: w[*][*] must have the same dimensions as x[*][*]")
        print("dimx: ", dimx, "    dimw=", dimw)
        exit()

    #print('rankw=',rankw)
    #print('dimw=',dimw[0])
    #print('dimx=',dimx[0])

    if rankw == 1:
        if dimw == 1:
            WGT = np.empty((dimx))
            WGT = 1.0

        if dimx[0] == dimw[0]: # broadcaset into 2d
            WGT = np.transpose(np.tile(wts,(dimx[1],1)),(1,0))
    else:
        WGT = wts

    #print(WGT)

    # if x/y  has _Fillvalue attribute; set WGT=0.0 where x or y = _FillValue
    WGT = np.where(x.mask, 0.0, WGT)
    WGT = np.where(y.mask, 0.0, WGT)
    #print(WGT)

    if opt == 0: # centered correlation
        sumWGT = np.ma.sum(WGT)
        xAvgArea = np.ma.sum(x*WGT)/sumWGT
        yAvgArea = np.ma.sum(y*WGT)/sumWGT

        xAnom = x - xAvgArea
        yAnom = y - yAvgArea

        xyCov = np.ma.sum(WGT*xAnom*yAnom)
        xAnom2 = np.ma.sum(WGT*xAnom**2)
        yAnom2 = np.ma.sum(WGT*yAnom**2)

    else:
        xyCov = np.ma.sum(WGT*x*y)
        xAnom2 = np.ma.sum(WGT*x**2)
        yAnom2 = np.ma.sum(WGT*y**2)

    NRMSE = np.sqrt(np.ma.sum(WGT*(x - y)**2) / np.ma.sum(WGT*x**2))

    RMSE = np.sqrt(np.ma.sum(WGT*(x - y)**2))

    if xAnom2 > 0.0 and yAnom2 > 0.0:
        r = xyCov/(np.sqrt(xAnom2)*np.sqrt(yAnom2))

    else:
        r = np.nan

    return r,NRMSE,RMSE

# =========================================================================================

def map_white (num_half_non_white,perc,cmap):
    '''
    Dec 30, 2020: create a new colormap with white in center with specific values
    ref: https://stackoverflow.com/questions/45196233/python-matplotlib-create-normal-colorbar-with-white-interval-at-specific-values
    https://stackoverflow.com/questions/41806285/how-to-change-colorbars-color-in-some-particular-value-interval

    parameters:
    num_half_no_white: number of colors in the upper/lower part
    perc: the relative percent of white colors w/ num_half_no_white

    example:

    arr = np.linspace(0, 50, 100).reshape((10, 10))
    fig, ax = plt.subplots(ncols=2)

    cmap = plt.get_cmap('jet')
    n = 30
    new_cmap = map_white(n,cmap)
    ax[0].imshow(arr, interpolation='nearest', cmap=cmap)
    ax[1].imshow(arr, interpolation='nearest', cmap=new_cmap)
    plt.show()

    '''
    import matplotlib
    import numpy as np
    import matplotlib.pyplot as plt

    n=num_half_non_white
    x = 0.5
    lower = cmap(np.linspace(0, x, n))
#     white = np.ones((100-2*n,4))
    white = np.ones((int(perc*n),4))
    upper = cmap(np.linspace(1-x, 1, n))
    colors = np.vstack((lower, white, upper))
    print('colors.shape=',colors.shape,'upper.shape=',upper.shape,'white.shape=',white.shape,'lower.shape=',lower.shape)
    tmap = matplotlib.colors.LinearSegmentedColormap.from_list('map_white', colors)

    return tmap

# =========================================================================================
def mask_land(lons,lats,data,land=True):
    lons_here = np.where(lons>180,lons-360,lons)
    lon_grid,lat_grid = np.meshgrid(lons_here,lats)
    globe_land_mask = globe.is_land(lat_grid,lon_grid)
    if len(data.shape) == len(globe_land_mask.shape):
        globe_land_mask_nd = globe_land_mask
    elif len(data.shape) - len(globe_land_mask.shape) == 1: # have time dimension
        globe_land_mask_nd = np.tile(globe_land_mask,(data.shape[0],1,1))
    elif len(data.shape) - len(globe_land_mask.shape) == 2: # have ctp and tau dimension
        globe_land_mask_nd = np.tile(globe_land_mask,(data.shape[0],data.shape[1],1,1))
    elif len(data.shape) - len(globe_land_mask.shape) == 3: # have time, ctp and tau dimension
        globe_land_mask_nd = np.tile(globe_land_mask,(data.shape[0],data.shape[1],data.shape[2],1,1))
    #print(globe_land_mask_nd.shape)

    data_new = np.ma.masked_where(globe_land_mask_nd==land,data)

    return data_new

# ========================================================================
# Multiple dimension linear regression along one dimension.
# following genutil.statistics.linearregression() function
# Link: https://github.com/CDAT/genutil/blob/master/Lib/statistics.py#L86
# ========================================================================

def linearregression_nd(y, x, error=None, probability=None,
                       noslope=None, nointercept=None):
    """
    returns slope/intercept for linear regression of dependant var y and indep var x
    also possibly returns error and P values.
    
    YQIN: keep the dimension to do the regression as the first dimension. For example, the shape of y is (5,3), and the shape of x is (5,1). 
    """
    if (not (noslope is None or noslope == 0)) and (
            not (nointercept is None or nointercept == 0)):
        raise StatisticsError(
            'Error in __linearregression, at least one of the following argument as to be None:' +
            'noslope or nointercept, you are requesting nothing back !')
    if (probability is not None) and (error is None):
        raise StatisticsError(
            'Error in __linearregression, error must not be None if probability is defined, probability is:' +
            str(probability))
    if error is not None:
        if error > 3:
            raise StatisticsError(
                "Error in __linearregression, error must be None (0), 1, ,2 or 3")

    xmean = numpy.ma.average(x, axis=0)
    ymean = numpy.ma.average(y, axis=0)
    x = x - xmean
    y = y - ymean
    xy = numpy.ma.sum(y * x, axis=0)
    xx = numpy.ma.sum(x * x, axis=0)
    slope = xy / xx
    intercept = ymean - slope * xmean
    V = []
    if noslope is None or noslope == 0:
        V.append(slope)
    if nointercept is None or nointercept == 0:
        V.append(intercept)
    if error is None or error == 0:
        return V
    elif error == 1:
        E = []
        n1 = numpy.ma.count(y, axis=0)
        # Unadjusted errors
        res = (y + ymean) - (intercept + (x + xmean) *
                             numpy.ma.resize(slope, numpy.ma.shape(y)))  # x2
        ssd = numpy.ma.sum(res * res, axis=0)
        amsd1 = ssd / (n1 - 2.)  # amsd1=ssd/idfd1
        if noslope is None or noslope == 0:
            E.append(numpy.ma.sqrt(amsd1 / xx))
        if nointercept is None or nointercept == 0:
            s1 = xmean * xmean / xx + 1. / n1
            E.append(numpy.ma.sqrt(amsd1 * s1))
        if probability is None or probability == 0:
            return V, E
        else:
            Pt1 = []
            Pt2 = []
            Pf1 = []
            Pf2 = []
            f = numpy.ma.sum(y * y, axis=0) - ssd  # ssr
            f = f / amsd1
            aa1 = n1 / 2.0
            if noslope is None or noslope == 0:
                tb1 = slope / E[0]
                xx3 = n1 / (n1 + (tb1 * tb1))
                Pt1.append(__betai1(aa1, .5, xx3))
                Pt2.append(None)
                Pf1.append(__probf1(f, 1, n1 - 2, 1))
                Pf2.append(__probf1(f, 1, n1 - 2, 2))
            if nointercept is None or nointercept == 0:
                tb1 = V[-1] / E[-1]
                xx3 = n1 / (n1 + (tb1 * tb1))
                Pt1.append(__betai1(aa1, .5, xx3))
                Pt2.append(None)
                Pf1.append(__probf1(f, 1, n1 - 2, 1))
                Pf2.append(__probf1(f, 1, n1 - 2, 2))
            return V, E, Pt1, Pt2, Pf1, Pf2
    else:
        E = []
        # Adjusted error from residual
        n1 = numpy.ma.count(y, axis=0)
        res = (y + ymean) - (intercept + numpy.ma.resize(slope,
                                                         numpy.ma.shape(y)) * (x + xmean))  # x2
        ssd = numpy.ma.sum(res * res, axis=0)
        if error == 2:
            ac = __autocorrelation(res, 1, centered=1, partial=0)
            rdfd2 = 1.0 + ac
            rdfd2 = (1.0 - ac) / rdfd2
        elif error == 3:
            ac = __autocorrelation(y + ymean, 1, centered=1, partial=0)
            rdfd2 = 1.0 + ac
            rdfd2 = (1.0 - ac) / rdfd2
        rneff = n1 * rdfd2  # rneff
        amsd2 = ssd / (rneff - 2.)   # ssd/rdfd2
        if noslope is None or noslope == 0:
            E.append(numpy.ma.sqrt(amsd2 / xx))
        if nointercept is None or nointercept == 0:
            s1 = xmean * xmean / xx + 1. / n1
            E.append(numpy.ma.sqrt(amsd2 * s1))
        if probability is None or probability == 0:
            return V, E
        else:
            Pt1 = []
            Pt2 = []
            Pf1 = []
            Pf2 = []
            f = numpy.ma.sum(y * y, axis=0) - ssd  # ssr = amsr
            amsd1 = ssd / (n1 - 2.)  # amsd1=ssd/idfd1
            f = f / amsd1  # amsr/ssd
            aa1 = n1 / 2.0
            aa2 = rneff / 2.0
            if noslope is None or noslope == 0:
                tb2 = slope / E[0]
                xx1 = n1 / (n1 + (tb2 * tb2))
                xx2 = rneff / (rneff + (tb2 * tb2))
                Pt1.append(__betai1(aa1, .5, xx1))
                Pt2.append(__betai1(aa2, .5, xx2))
                Pf1.append(__probf1(f, 1, n1 - 2, 1))
                Pf2.append(__probf1(f, 1, n1 - 2, 2))
            if nointercept is None or nointercept == 0:
                tb2 = V[-1] / E[-1]
                xx1 = n1 / (n1 + (tb2 * tb2))
                xx2 = rneff / (rneff + (tb2 * tb2))
                Pt1.append(__betai1(aa1, .5, xx1))
                Pt2.append(__betai1(aa2, .5, xx2))
                Pf1.append(__probf1(f, 1, n1 - 2, 1))
                Pf2.append(__probf1(f, 1, n1 - 2, 2))
            return V, E, Pt1, Pt2, Pf1, Pf2

# ----------------------------------------------------------------
def __laggedcovariance(x, y, lag=1, centered=1, partial=1):
    """
    Function: __laggedcovariance

    Description of function:
        Does the main computation for returning lagged covariance. See
        documentation of laggedcovariance() for details.
    """
    if lag == 0:
        return __covariance(x, y, centered=centered)

    if partial == 1:
        if lag > 0:
            x = x[lag:]
            y = y[:-lag]
        else:
            x = x[:lag]
            y = y[-lag:]

    if centered == 1:
        xmean = numpy.ma.average(x, axis=0)
        ymean = numpy.ma.average(y, axis=0)
    else:
        xmean = 0.
        ymean = 0.
    x = x - xmean
    y = y - ymean
    del(xmean)
    del(ymean)

    if partial == 1:
        tmp = x * y
    else:
        if lag > 0:
            tmp = x[lag:] * y[:-lag]
        else:
            tmp = x[:-lag] * y[lag:]
    return numpy.ma.sum(tmp, axis=0) / numpy.ma.count(x * y, axis=0)

# ----------------------------------------------------------------
def __covariance(x, y, weights=None, centered=1, biased=1):
    """
    Function: __covariance

    Description of function:
        Does the main computation for returning covariance. See documentation
        of covariance() for details.
    """
    if weights is not None and biased != 1:
        raise StatisticsError(
            'Error in covariance, you cannot have weights and unbiased together')

    if centered == 1:
        xmean = numpy.ma.average(x, weights=weights, axis=0)
        ymean = numpy.ma.average(y, weights=weights, axis=0)
        x = x - xmean
        y = y - ymean
        del(xmean)
        del(ymean)
    #
    if weights is None:
        weights = numpy.ma.ones(x.shape, dtype=x.dtype.char)
    if not ((x.mask is None) or (x.mask is numpy.ma.nomask)):
        weights = numpy.ma.masked_where(x.mask, weights)
    if biased == 1:
        cov = numpy.ma.sum(x * y * weights, axis=0) / \
            numpy.ma.sum(weights, axis=0)
    else:
        cov = numpy.ma.sum(x * y, axis=0) / (numpy.ma.count(x * y, axis=0) - 1)

    return cov

# ----------------------------------------------------------------
def __autocovariance(x, lag, centered=1, partial=1):
    """
    Function: __autocovariance

    Description of function:
        Does the main computation for returning autocovariance. See
        documentation of autocovariance() for details.
    """
    return __laggedcovariance(x, x, lag, centered=centered, partial=partial)

# ----------------------------------------------------------------
def __autocorrelation(x, lag, centered=1, partial=1):
    """
    Function: __autocorrelation

    Description of function:
        Does the main computation for returning autocorrelation. See
        documentation of autocorrelation() for details.
    """
    if partial == 1 and centered == 1 and lag != 0:
        mean1 = numpy.ma.average(x[:-lag], axis=0)
        mean2 = numpy.ma.average(x[lag:], axis=0)
        x1 = x[:-lag] - mean1
        x2 = x[lag:] - mean2
        num = numpy.ma.sum(x1 * x2, axis=0)
        den = numpy.ma.sum(numpy.ma.power(x1, 2), axis=0) * \
            numpy.ma.sum(numpy.ma.power(x2, 2), axis=0)
        return num / numpy.ma.sqrt(den)
    else:
        return __autocovariance(x, lag, centered=centered, partial=partial) / \
            __autocovariance(x, 0, centered=centered, partial=partial)

# ----------------------------------------------------------------
def __gammln1(x):
    cof = [76.18009172947146, -86.50532032941677,
           24.01409824083091, -1.231739572450155,
           .1208650973866179E-2, -.5395239384953E-5]
    stp = 2.5066282746310005
    y = x * 1.
    tmp = x + 5.5
    tmp = (x + 0.5) * numpy.ma.log(tmp) - tmp
    ser = 1.000000000190015
    for j in range(6):
        y = y + 1.
        ser = ser + cof[j] / y
    return tmp + numpy.ma.log(stp * ser / x)

# ----------------------------------------------------------------
def __betacf1(a, b, x):
    MAXIT = 100
    EPS = 3.E-7
    FPMIN = 1.E-30
    qab = a + b
    qap = a + 1.
    qam = a - 1.
    c = 1.
    d = 1. - qab * x / qap
    d = numpy.ma.where(numpy.ma.less(numpy.ma.absolute(d), FPMIN), FPMIN, d)
    d = 1. / d
    h = d
    for m in range(1, MAXIT + 1):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1. + aa * d
        d = numpy.ma.where(
            numpy.ma.less(
                numpy.ma.absolute(d),
                FPMIN),
            FPMIN,
            d)
        c = 1. + aa / c
        c = numpy.ma.where(
            numpy.ma.less(
                numpy.ma.absolute(c),
                FPMIN),
            FPMIN,
            c)
        d = 1. / d
        h = h * d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1. + aa * d
        d = numpy.ma.where(
            numpy.ma.less(
                numpy.ma.absolute(d),
                FPMIN),
            FPMIN,
            d)
        c = 1. + aa / c
        c = numpy.ma.where(
            numpy.ma.less(
                numpy.ma.absolute(c),
                FPMIN),
            FPMIN,
            c)
        d = 1. / d
        delet = d * c
        h = h * delet
        if numpy.ma.allclose(delet, numpy.ones(
                delet.shape), atol=EPS, rtol=0.):
            break
    h = numpy.ma.masked_where(
        numpy.ma.greater(
            numpy.ma.absolute(
                delet - 1.), EPS), h)
    return h

# ----------------------------------------------------------------
def __betai1(a, b, x):
    bt = numpy.ma.logical_or(numpy.ma.equal(x, 0.), numpy.ma.equal(x, 1.))
    bt = numpy.ma.where(bt, 0., numpy.ma.exp(
        __gammln1(a + b) - __gammln1(a) - __gammln1(b) +
        a * numpy.ma.log(x) + b * numpy.ma.log(1. - x)
    )
    )
    return numpy.ma.where(numpy.ma.less(x, (a + 1.) / (a + b + 2.)),
                          bt * __betacf1(a, b, x) / a,
                          1. - bt * __betacf1(b, a, 1. - x) / b)


# ----------------------------------------------------------------
def __probnd1(x):
    """
    FUNCTION PROBND1.

    Calculates the area under a normal curve (mean=0.0, variance=1.0)
    to the right of x. The accuracy is better than 7.5 * 10.**-8.

    REFERENCE:

    M. Abramowitz and I.A. Stegun.
    Handbook of Mathematical Functions.
    Dover, 1970, pages 931-932 (26.2.1 and 26.2.17).
"""
    b1 = 0.319381530
    b2 = -0.356563782
    b3 = 1.781477937
    b4 = -1.821255978
    b5 = 1.330274429
    p = 0.2316419
    t = 1.0 / (1.0 + (p * x))
    term1 = ((((b1 * t) + (b2 * (t ** 2))) +
              (b3 * (t ** 3))) + (b4 * (t ** 4))) + \
        (b5 * (t ** 5))
    z = (1.0 / numpy.ma.sqrt(2.0 * numpy.pi)) * numpy.ma.exp(- ((x * x) / 2.0))
    return numpy.ma.where(numpy.ma.greater(x, 7.), 0., z * term1)


# ----------------------------------------------------------------
def __probf1(y, n1, n2, id):
    """
    FUNCTION PROBF1.

    The output is either the one- or two-tailed test area: i.e., the
    area under an F-curve (with N1 and N2 degrees of freedom) to the
    right of X if X exceeds 1.0 (one-tailed test) or twice this area
    (two-tailed test).

    Note: if X is less than 1.0, this function gives the area to the
    right of 1/X with reversed order for the degrees of freedom. This
    ensures the accuracy of the numerical algorithm.

    REFERENCE:

    M. Abramowitz and I.A. Stegun.
    Handbook of Mathematical Functions.
    Dover, 1970, page 947 (26.6.15).

    ** INPUT **
    real y            Calculated F-value
    real x            Inverse of Y if Y is less than 1.0
    integer n1, n2    Degrees of freedom
    integer id        Identifier for one- or two-tailed test

    ** OUTPUT **
    real probf1       Significance level (p-value) for F-value

    EXTERNALS:

    function PROBND1 - Calculates the area under a normal curve.
    """
    ly = numpy.ma.less(y, 1.)
    x = numpy.ma.where(ly, 1. / numpy.ma.array(y), y)
    n = numpy.ma.where(ly, n1, n2)
    n1 = numpy.ma.where(ly, n2, n1)
    n2 = numpy.ma.where(ly, n, n2)
    term1 = 2.0 / (9.0 * n1)
    term2 = 2.0 / (9.0 * n2)
    term3 = ((x ** (1.0 / 3.0)) * (1.0 - term2)) - (1.0 - term1)
    term4 = numpy.ma.sqrt(term1 + ((x ** (2.0 / 3.0)) * term2))
    term5 = term3 / term4
    probf1 = id * __probnd1(term5)

    #     The numerical algorithm can have problems when the F-value is
    #     close to 1.0 and the degrees of freedom are small. Therefore,
    #     insure that the probabilities returned cannot exceed 1.0.

    return numpy.ma.where(numpy.ma.greater(probf1, 1.), 1., probf1)

# ======================================================
def weighted_annual_mean(time, obs):
  """
  weight by days in each month
  """
  import xarray as xr 

  # Determine the month length
  month_length = time.dt.days_in_month

  # Calculate the weights
  wgts = month_length.groupby("time.year") / month_length.groupby("time.year").sum()

  # Make sure the weights in each year add up to 1
  np.testing.assert_allclose(wgts.groupby("time.year").sum(xr.ALL_DIMS), 1.0)

  # Setup our masking for nan values
  cond = obs.isnull()
  ones = xr.where(cond, 0.0, 1.0)

  # Calculate the numerator
  obs_sum = (obs * wgts).resample(time="AS").sum(dim="time")

  # Calculate the denominator
  ones_out = (ones * wgts).resample(time="AS").sum(dim="time")

  # Return the weighted average
  return obs_sum / ones_out

# ======================================================
def area_averager(data_plot_xr):
    '''
    calculate weighted area mean
    input data is xarray DataArray
    '''
    weights = np.cos(np.deg2rad(data_plot_xr.lat))
    weights.name = "weights"
    # available in xarray version 0.15 and later
    data_weighted = data_plot_xr.weighted(weights)

    weighted_mean = data_weighted.mean(("lat", "lon"))

    return weighted_mean

# ======================================================
def save_big_dataset(dic_mod,outfile,comment):
    '''
    create a big dataset based on all variables in a dictionary and save to netcdf file.
    '''
    import xarray as xr 

    datalist = []
    for svar in dic_mod.keys():
        data = xr.DataArray(dic_mod[svar],name=svar)
        datalist.append(data)

    data_big = xr.merge(datalist,compat='override')
    data_big.attrs['comments'] = comment 

    #data_big.to_netcdf(outfile,encoding={'time':{'dtype': 'i4'},'bin':{'dtype':'i4'}})
    data_big.to_netcdf(outfile)


