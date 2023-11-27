import scanpy as sc
import pandas as pd
import numpy as np
import scipy
import matplotlib.pyplot as plt
from sklearn.feature_selection import f_regression, mutual_info_regression, mutual_info_classif
from scipy.stats import gaussian_kde
from scipy.optimize import curve_fit
import matplotlib as mpl
from matplotlib import rc_context
import seaborn as sns
mpl.rc('image', cmap='binary')
mpl.rcParams['pdf.fonttype'] = 42

def dynamic_expression_for_de_genes_ct(adata, ct, days_to_include, ):
    df = adata.to_df()
    df['Order'] = adata.obs['Order']
    adata.obs['Order_int'] = [int(x) for x in adata.obs['Order']]
    df['Order_int'] = adata.obs['Order_int']
    master_dict = {}
    de = sc.get.rank_genes_groups_df(adata, group=ct)
    de = de[de['pvals_adj']<0.05]
    tmp = df[df['Order'].isin(days_to_include)]
    for x in de['names']:
        di = {}
        for y in de['names']:
            if x!=y:
                print('===================================================')
                dremi = scprep.stats.knnDREMI(tmp[x], tmp[y], plot=True)
                plot_1f(adata[adata.obs['Order'].isin(days_to_include)],
                        genes=[x, y], ct=ct)
                gb = tmp.groupby(by='Order_int').mean()
                fig, ax = plt.subplots()
                ax.scatter(x=gb.index, y=gb[x], c='b', label=x)
                ax.scatter(x=gb.index, y=gb[y], c='r', label=y)
                plt.legend(loc='upper left')
                plt.show()
                di[y] = dremi
        master_dict[x] = di
    return master_dict

def dynamic_expression_for_de_genes_ct_noplot(adata, ct, days_to_include, ):
    df = adata.to_df()
    df['Order'] = adata.obs['Order']
    df['Order'] = [str(x) for x in df['Order']]
    adata.obs['Order_int'] = [int(x) for x in adata.obs['Order']]
    df['Order_int'] = adata.obs['Order_int']
    master_dict = {}
    de = sc.get.rank_genes_groups_df(adata, group=ct)
    de = de[de['pvals_adj']<0.05]
    tmp = df[df['Order'].isin(days_to_include)]
    for x in de['names']:
        di = {}
        print(x)
        for y in de['names']:
            if x!=y:
                dremi = scprep.stats.knnDREMI(tmp[x], tmp[y])
                di[y] = dremi
        master_dict[x] = di
    return master_dict

def ct_mi(regs, days, thresh):
    tmp = regs[regs.obs['Order'].isin(days)]
    regs_input = tmp.to_df()
    mi = mutual_info_classif(X=regs_input, y=tmp.obs['Order_int'], discrete_features=True)
    mi_dict = {k: v for k, v in sorted(dict(zip(regs_input.columns, mi)).items(), key=lambda item: item[1], 
                                       reverse=True)}
    regs_input['Order_int'] = tmp.obs['Order_int']
    
    for x in list(mi_dict.keys()):
        if mi_dict[x] > thresh:
            print('========================')
            print('MI is ', str(mi_dict[x]))
            z_total = []
            input_t = tmp.obs.sort_values(by='Order_int')['Order_int']
            input_r = regs_input.sort_values(by='Order_int')[x]
            xy = np.vstack([input_t,input_r])
            z = list(gaussian_kde(xy)(xy))

            z_normed_by_day=[]
            prev_index = 0
            for q in set(tmp.obs['Order_int']):
                tmp_t = input_t[input_t==q]
                norm = [float(i)/max(z[prev_index:prev_index+len(tmp_t)]) for i in z[prev_index:prev_index+len(tmp_t)]]
                prev_index +=len(tmp_t)
                z_normed_by_day+=norm
            plt.scatter(x=input_t, y=input_r, c=z)
            plt.xlabel('Order')
            plt.ylabel(x)
            plt.show()


            plt.scatter(x=input_t, y=input_r, c=z_normed_by_day)
            plt.xlabel('Order')
            plt.ylabel(x)
            plt.show()

            sc.pl.umap(regs, color=x)
    return mi_dict

def ct_mi_no_plot(regs, days, thresh):
    tmp = regs[regs.obs['Order'].isin(days)]
    regs_input = tmp.to_df()
    mi = mutual_info_classif(X=regs_input, y=tmp.obs['Order_int'], discrete_features=True)
    mi_dict = {k: v for k, v in sorted(dict(zip(regs_input.columns, mi)).items(), key=lambda item: item[1], 
                                       reverse=True)}
    return mi_dict

def func(x, a, b, c, m): # Hill sigmoidal equation from zunzun.com
    return  (a * np.power(x, b) / (np.power(c, b) + np.power(x, b))) + m

def fit_reg_to_curve(xData, yData, input_params, reg):
    initialParameters = np.array(input_params)
    fittedParameters, pcov = curve_fit(func, xData, yData, initialParameters, 
                                    #bounds=param_bounds,
                                      maxfev=1000000
                                      )
    modelPredictions = func(xData, *fittedParameters)
    absError = modelPredictions - yData

    SE = np.square(absError) # squared errors
    MSE = np.mean(SE) # mean squared errors
    RMSE = np.sqrt(MSE) # Root Mean Squared Error, RMSE
    Rsquared = 1.0 - (np.var(absError) / np.var(yData))
    
    # Plot
    graphWidth = 800
    graphHeight = 600
    f = plt.figure(figsize=(graphWidth/100.0, graphHeight/100.0), dpi=100)
    axes = f.add_subplot(111)

    # first the raw data as a scatter plot
    z_total = []

    xy = np.vstack([xData,yData])
    z = list(gaussian_kde(xy)(xy))

    z_normed_by_day=[]
    prev_index = 0
    for q in set(xData):
        tmp_t = [i for i in xData if i==q]
        norm = [float(i)/max(z[prev_index:prev_index+len(tmp_t)]) for i in z[prev_index:prev_index+len(tmp_t)]]
        prev_index +=len(tmp_t)
        z_normed_by_day+=norm


    axes.scatter(x=xData, y=yData, c=z_normed_by_day)

    # create data for the fitted equation plot
    xModel = np.linspace(min(xData), max(xData))
    yModel = func(xModel, *fittedParameters)

    # now the model as a line plot
    axes.plot(xModel, yModel)

    axes.set_xlabel('Day', fontsize=20) # X axis data label
    axes.set_ylabel(reg, fontsize=20) # Y axis data label
    plt.xticks(fontsize=12)
    plt.tight_layout()
    plt.savefig('../../data/isscr-figures/reg_time'+reg+'.pdf', dpi=150)
    plt.show()
    return fittedParameters

def compute_hill_for_all_regs(regs, regs_names, days, mi_dict, ct):
    params_dict = {}
    for x in regs_names:
        print(x)
        tmp = regs[(regs.obs['Cell Types']==ct)&(regs.obs['Order'].isin(days))]
        regs_input = tmp.to_df()
        regs_input['Order_int'] = tmp.obs['Order_int']
        input_x = list(tmp.obs.sort_values(by='Order_int')['Order_int'])
        input_y = list(regs_input.sort_values(by='Order_int')[x])
        max_yData = max(input_y)
        input_y = [y/max_yData for y in input_y]
        
        # estimate params
        m = min(input_y)
        a = max(input_y)
        b = 20 # estimation of b shouldn't affect output too much
        days_mean = regs_input.sort_values(by='Order_int').groupby('Order_int').mean()[x]
        target = a - m
        c = days_mean.index[(np.abs(days_mean - target)).argmin()]
        xmax = max(input_x)
        xmin = min(input_x)
        h = abs(a - m)
        
        param_est = [a, b, c, m]

        params = fit_reg_to_curve(input_x, input_y, param_est, reg=x)
        params_dict[x] = params
        
    params_df = pd.DataFrame.from_dict(params_dict, orient='index', columns=['a', 'b', 'c', 'm'])
    params_df['MI'] = [mi_dict[x] for x in params_df.index]
    params_df['a-m'] = params_df['a'] - params_df['m']
    params_df['b*a'] = params_df['a'] * params_df['b']
    
    return params_df

def compute_hill_for_all_regs_v2(regs, regs_names, days):
    params_dict = {}
    for x in regs_names:
        print(x)
        tmp = regs[regs.obs['Order'].isin(days)]
        regs_input = tmp.to_df()
        regs_input['Order_int'] = tmp.obs['Order_int']
        input_x = list(tmp.obs.sort_values(by='Order_int')['Order_int'])
        input_y = list(regs_input.sort_values(by='Order_int')[x])
        
        # estimate params
        m = min(input_y)
        a = max(input_y)
        b = 20 # estimation of b shouldn't affect output too much
        days_mean = regs_input.sort_values(by='Order_int').groupby('Order_int').mean()[x]
        target = a - m
        c = days_mean.index[(np.abs(days_mean - target)).argmin()]
        xmax = max(input_x)
        xmin = min(input_x)
        h = abs(a - m)
        
        param_est = [a, b, c, m]

        params = fit_reg_to_curve(input_x, input_y, param_est, reg=x)
        params_dict[x] = params
        
    params_df = pd.DataFrame.from_dict(params_dict, orient='index', columns=['a', 'b', 'c', 'm'])
    params_df['a-m'] = params_df['a'] - params_df['m']
    params_df['b*a'] = params_df['a'] * params_df['b']
    
    return params_df

def plot_mean_expression_by_day(tc_regs, regs_to_plot, cells_in_group_1, 
                                cells_in_group_2, days, cell_anno1, cell_anno2, color1='red', color2='blue'):
    tmp = tc_regs[tc_regs.obs['Order'].isin(days)]
    tmp1 = tmp[tmp.obs_names.isin(cells_in_group_1)]
    tmp2 = tmp[tmp.obs_names.isin(cells_in_group_2)]
    
    tmp1_df = tmp1.to_df()
    tmp2_df = tmp2.to_df()
    tmp1_df['Order'] = tmp1.obs['Order']
    tmp2_df['Order'] = tmp2.obs['Order']
    
    tmp1_df = tmp1_df.groupby('Order').mean()
    tmp2_df = tmp2_df.groupby('Order').mean()
    
    for x in regs_to_plot:
        fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(12, 5))
        ax.plot(tmp1_df.index, tmp1_df[x], marker='o', color=color1, label=cell_anno1)
        ax.plot(tmp2_df.index, tmp2_df[x], marker='o', color=color2, label=cell_anno2)
        plt.xlabel('Day')
        plt.ylabel('Mean '+x+' Activity')
        ax.legend()
        plt.show()
        
def plot_mean_expression_by_day_v2(tc_regs, regs_to_plot, cells_in_each_group=[],
                                days=[], cell_annos=[], anno_colors=['red', 'blue', 'green'],
                                  min_cells_day=0):
            
 
    tmp = tc_regs[tc_regs.obs['Order'].isin(days)]

    df_list = []

    for i, x in enumerate(cells_in_each_group):
        tmp1 = tmp[tmp.obs_names.isin(cells_in_each_group[i])]
        tmp1_df = tmp1.to_df()
        tmp1_df['Order'] = tmp1.obs['Order']
        ct_days = []
        for y in set(tmp1_df['Order']):
            if len(tmp1_df[tmp1_df['Order']==y].index) >= min_cells_day:
                ct_days.append(y)
        tmp1_df = tmp1_df[tmp1_df['Order'].isin(ct_days)]
        tmp1_df = tmp1_df.groupby('Order').mean()
        df_list.append(tmp1_df)
    
    for x in regs_to_plot:
        fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(12, 5))
        for i, y in enumerate(df_list):
            ax.plot(y.index, y[x], marker='o', color=anno_colors[i], label=cell_annos[i])
            plt.xlabel('Day')
            plt.ylabel('Mean '+x+' Activity')
        
        ax.legend()
        plt.show()
        
def fit_reg_to_curve_v3(xData, yData, input_params, reg):
    initialParameters = np.array(input_params)
    fittedParameters, pcov = curve_fit(func, xData, yData, initialParameters, 
                                    #bounds=param_bounds,
                                      maxfev=1000000
                                      )

    modelPredictions = func(xData, *fittedParameters)
    absError = modelPredictions - yData
    
    stdev = np.sqrt(np.diag(pcov))
    #print('STDEV', stdev)
    
    #print('Fitted var / std dev', fittedParameters / stdev)

    SE = np.square(absError) # squared errors
    MSE = np.mean(SE) # mean squared errors
    RMSE = np.sqrt(MSE) # Root Mean Squared Error, RMSE
    Rsquared = 1.0 - (np.var(absError) / np.var(yData))
    param_variance = np.diag(pcov)
    
    #print('Fitted params', fittedParameters)
    #print('Mean squared errors', MSE)
    print('R squared', Rsquared)
    
    fit_check = Rsquared >= 0.1
    
    #     print('Test', abs(fittedParameters[1]) - abs(3*stdev[1]))
    #     if (abs(fittedParameters[1]) - abs(3*stdev[1])) > 0:
    #         print("Passed filter")
    
    # Plot
    graphWidth = 800
    graphHeight = 600
    f = plt.figure(figsize=(graphWidth/100.0, graphHeight/100.0), dpi=100)
    axes = f.add_subplot(111)

    # first the raw data as a scatter plot
    z_total = []

    xy = np.vstack([xData,yData])
    z = list(gaussian_kde(xy)(xy))

    z_normed_by_day=[]
    prev_index = 0
    for q in set(xData):
        tmp_t = [i for i in xData if i==q]
        norm = [float(i)/max(z[prev_index:prev_index+len(tmp_t)]) for i in z[prev_index:prev_index+len(tmp_t)]]
        prev_index +=len(tmp_t)
        z_normed_by_day+=norm


    axes.scatter(x=xData, y=yData, c=z_normed_by_day)

    # create data for the fitted equation plot
    xModel = np.linspace(min(xData), max(xData))
    yModel = func(xModel, *fittedParameters)

    # now the model as a line plot
    axes.plot(xModel, yModel, color='red')

    axes.set_xlabel('Day', fontsize=20) # X axis data label
    axes.set_ylabel(reg, fontsize=20) # Y axis data label
    plt.xticks(fontsize=12)
    plt.tight_layout()
    plt.show()
    return fittedParameters, fit_check, Rsquared

def compute_hill_for_all_regs_v3(regs, regs_names, days):
    params_dict = {}
    fit_check_ls = []
    rsquared_ls = []
    for x in regs_names:
        print(x)
        tmp = regs[regs.obs['Order'].isin(days)]
        regs_input = tmp.to_df()
        regs_input['Order_int'] = tmp.obs['Order_int']
        input_x = list(tmp.obs.sort_values(by='Order_int')['Order_int'])
        input_y = list(regs_input.sort_values(by='Order_int')[x])
        
        # estimate params
        m = min(input_y)
        a = max(input_y)
        b = 20 # estimation of b shouldn't affect output too much
        days_mean = regs_input.sort_values(by='Order_int').groupby('Order_int').mean()[x]
        target = a - m
        c = days_mean.index[(np.abs(days_mean - target)).argmin()]
        xmax = max(input_x)
        xmin = min(input_x)
        h = abs(a - m)
        
        param_est = [a, b, c, m]

        params, fit_check, rsquared = fit_reg_to_curve_v3(input_x, input_y, param_est, reg=x)
        params_dict[x] = params
        fit_check_ls.append(fit_check)
        rsquared_ls.append(rsquared)
        
    params_df = pd.DataFrame.from_dict(params_dict, orient='index', columns=['a', 'b', 'c', 'm'])
    params_df['a-m'] = params_df['a'] - params_df['m']
    params_df['b*a'] = params_df['a'] * params_df['b']
    params_df['Fit Check'] = fit_check_ls
    params_df['Rsquared'] = rsquared_ls
    
    return params_df

def compute_hill_for_all_regs_v4(regs, regs_names, days):
    params_dict = {}
    fit_check_ls = []
    rsquared_ls = []
    for x in regs_names:
        tmp = regs[regs.obs['Order'].isin(days)]
        regs_input = tmp.to_df()
        regs_input['Order_int'] = tmp.obs['Order_int']
        input_x = list(tmp.obs.sort_values(by='Order_int')['Order_int'])
        input_y = list(regs_input.sort_values(by='Order_int')[x])
        
        # estimate params
        m = min(input_y)
        a = max(input_y)
        b = 20 # estimation of b shouldn't affect output too much
        days_mean = regs_input.sort_values(by='Order_int').groupby('Order_int').mean()[x]
        target = a - m
        c = days_mean.index[(np.abs(days_mean - target)).argmin()]
        xmax = max(input_x)
        xmin = min(input_x)
        h = abs(a - m)
        
        param_est = [a, b, c, m]

        params, fit_check, rsquared, stdev = fit_reg_to_curve_v4(input_x, input_y, param_est, reg=x)
        params_stdev = list(params) + list(stdev)

        params_dict[x] = params_stdev
        fit_check_ls.append(fit_check)
        rsquared_ls.append(rsquared)
    
    params_df = pd.DataFrame.from_dict(params_dict, orient='index',
                                       columns=['a', 'b', 'c', 'm', 'a_stdev', 'b_stdev', 'c_stdev', 'm_stdev'])
    params_df['a-m'] = params_df['a'] - params_df['m']
    params_df['b*a'] = params_df['a'] * params_df['b']
    params_df['Fit Check'] = fit_check_ls
    params_df['Rsquared'] = rsquared_ls

    return params_df

def fit_reg_to_curve_v4(xData, yData, input_params, reg):
    initialParameters = np.array(input_params)
    fittedParameters, pcov = curve_fit(func, xData, yData, initialParameters, 
                                    #bounds=param_bounds,
                                      maxfev=1000000
                                      )

    modelPredictions = func(xData, *fittedParameters)
    absError = modelPredictions - yData
    
    stdev = np.sqrt(np.diag(pcov))

    SE = np.square(absError) # squared errors
    MSE = np.mean(SE) # mean squared errors
    RMSE = np.sqrt(MSE) # Root Mean Squared Error, RMSE
    Rsquared = 1.0 - (np.var(absError) / np.var(yData))
    param_variance = np.diag(pcov)

    fit_check = Rsquared >= 0.1

    return fittedParameters, fit_check, Rsquared, stdev

   
def plot_c_stdev(input_df, title='', path=False, axticks=False, markerfacecolor='red', xrange=False): # takes output of compute_hill_for_all_regs_v3
    #fig.subplots_adjust(left=None, right=None, bottom=None, wspace=None, hspace=0, top=len(input_df.index)/4)
    input_df = input_df.sort_values('c', ascending=False)
    ls=[]
    input_df = input_df[(input_df['c'] < 100) & (input_df['c_stdev'] < 1)]
    height = len(input_df.index)/4
    fig, ax = plt.subplots(figsize=(15, height))
    for x in input_df['c_stdev']:
        ls.append(x)
    input_df['c_stdev'] = ls
    input_y = [x for x in range(len(input_df.index))]
    input_y.reverse()
    ytick_label = list(input_df.index)
    height_per_reg = height / len(ytick_label)
    yticks = [x*height_per_reg for x in range(len(ytick_label))]
    input_x = [x*24 for x in input_df['c']]
    xerr = [x*24 for x in input_df['c_stdev']]
    ax.errorbar(x=input_x, y=yticks, xerr=xerr, fmt='o', 
                ecolor='black', markerfacecolor=markerfacecolor, markeredgecolor=markerfacecolor, markersize=4,
               elinewidth=1)
    ax.set_yticks(yticks)
    ax.set_yticklabels(ytick_label)
    ymin=min(yticks) - height_per_reg/2
    ymax=max(yticks) + height_per_reg/2
    ax.set_ylim([ymin, ymax])
    ax.set_xlabel('Hours in Culture')
    ax.set_ylabel('Regulon')
    ax.set_title(title)
    if xrange:
        ax.set_xlim(xrange)
    if axticks:
        ax.set_xticks(axticks)
    if path:
        plt.savefig(path, dpi=150)
    plt.show()
    return input_x, input_y
    

def fit_reg_to_curve_save(xData, yData, input_params, reg, path):
    initialParameters = np.array(input_params)
    fittedParameters, pcov = curve_fit(func, xData, yData, initialParameters, 
                                    #bounds=param_bounds,
                                      maxfev=1000000
                                      )

    modelPredictions = func(xData, *fittedParameters)
    absError = modelPredictions - yData
    
    stdev = np.sqrt(np.diag(pcov))
    #print('STDEV', stdev)
    
    #print('Fitted var / std dev', fittedParameters / stdev)

    SE = np.square(absError) # squared errors
    MSE = np.mean(SE) # mean squared errors
    RMSE = np.sqrt(MSE) # Root Mean Squared Error, RMSE
    Rsquared = 1.0 - (np.var(absError) / np.var(yData))
    param_variance = np.diag(pcov)
    
    #print('Fitted params', fittedParameters)
    #print('Mean squared errors', MSE)
    print('R squared', Rsquared)
    
    fit_check = Rsquared >= 0.1
    
    #     print('Test', abs(fittedParameters[1]) - abs(3*stdev[1]))
    #     if (abs(fittedParameters[1]) - abs(3*stdev[1])) > 0:
    #         print("Passed filter")
    
    # Plot
    graphWidth = 800
    graphHeight = 600
    with rc_context({'figure.figsize': (6.05, 5)}):
        f = plt.figure(figsize=(graphWidth/100.0, graphHeight/100.0), dpi=150)
        axes = f.add_subplot(111)

        # first the raw data as a scatter plot
        z_total = []

        xy = np.vstack([xData,yData])
        z = list(gaussian_kde(xy)(xy))

        z_normed_by_day=[]
        prev_index = 0
        for q in set(xData):
            tmp_t = [i for i in xData if i==q]
            norm = [float(i)/max(z[prev_index:prev_index+len(tmp_t)]) for i in z[prev_index:prev_index+len(tmp_t)]]
            prev_index +=len(tmp_t)
            z_normed_by_day+=norm


        axes.scatter(x=xData, y=yData, c=z_normed_by_day)

        # create data for the fitted equation plot
        xModel = np.linspace(min(xData), max(xData))
        yModel = func(xModel, *fittedParameters)

        # now the model as a line plot
        axes.plot(xModel, yModel, color='red')

        axes.set_xlabel('Day', fontsize=20) # X axis data label
        axes.set_ylabel(reg, fontsize=20) # Y axis data label

        plt.xticks(fontsize=12)
        plt.tight_layout()
        if path:
            plt.savefig(path+'_'+reg+'.pdf', dpi=200)
        plt.show()
    return fittedParameters, fit_check, Rsquared


def compute_hill_for_all_regs_save(regs, regs_names, days, path=False):
    params_dict = {}
    fit_check_ls = []
    rsquared_ls = []
    for x in regs_names:
        tmp = regs[regs.obs['Order'].isin(days)]
        regs_input = tmp.to_df()
        regs_input['Order_int'] = tmp.obs['Order_int']
        input_x = list(tmp.obs.sort_values(by='Order_int')['Order_int'])
        input_y = list(regs_input.sort_values(by='Order_int')[x])
        
        # estimate params
        m = min(input_y)
        a = max(input_y)
        b = 20 # estimation of b shouldn't affect output too much
        days_mean = regs_input.sort_values(by='Order_int').groupby('Order_int').mean()[x]
        target = a - m
        c = days_mean.index[(np.abs(days_mean - target)).argmin()]
        xmax = max(input_x)
        xmin = min(input_x)
        h = abs(a - m)
        
        param_est = [a, b, c, m]
        fit_reg_to_curve_save(input_x, input_y, param_est, reg=x, path=path)