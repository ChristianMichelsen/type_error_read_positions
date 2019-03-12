#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 18 12:15:39 2018

@author: michelsen
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from tqdm import tqdm

from src.data.prepare_reads import get_filename_and_lenght 
from src.extra_functions import (
                                 is_linux,
                                 get_ML_res,
                                 get_error_rates,
                                 get_read_lengt_count,
                                 get_strand_count,
                                 plot_confusion_matrix,
                                 remove_Ns,
                                 precision_recall_bla_to_df,
                                 )


save_plots = True

do = 'ancient'
do = 'modern'
# do = 'gargamel'


verbose = True
force_rerun = False
do_plotting = True
do_remove_Ns = True

cores = 20


if do == 'ancient':
    print("\nRunning on ancient DNA")
    filename = 'ESW_LRUE_MA2621_error_test.txt' # ancient data
    filename_bam = 'ESW_LRUE_MA2621_L1_S1_S83_L006.sort.rmdup.bam'
    refname = 'hs37d5.fa'
elif do == 'modern':
    print("\nRunning on modern DNA")
    filename = 'NA12400_error_test.txt' # modern dna
    # filename = 'NA12400_small_bam_1_000_000_error_test.txt' # modern dna
    filename_bam = 'NA12400_small_bam_1_000_000.bam'
    refname = 'hs37d5.fa'
elif do == 'gargamel':
    print("\nRunning on gargamel simulated DNA")
    filename = 'gargamel_1_000_000.txt'
    filename_bam = 'gargamel_1_000_000.bam'
    refname = 'horse_chrom31.fa'


print(filename)
plot_prefix = f"../figures/{filename.split('_')[0]}_plot_"
file_processed_in, N_reads = get_filename_and_lenght(filename)


#%% =============================================================================
#  Compare own reading with mapDamage reads from bam-file
# =============================================================================

if False:
    
    from mapdamage_custom import mapDamage_main_test
    mapDamage_main_test(refname, filename_bam, filename)


#%% =============================================================================
# https://www.blopig.com/blog/2016/08/processing-large-files-using-python/
# =============================================================================

df = get_ML_res(file_processed_in, cores, force_rerun, N_reads, N_splits=100)

if do_remove_Ns:
    df = remove_Ns(df)

df_forward = df.loc[((df['strand'] & 0x10) == 0)]
df_reverse = df.loc[((df['strand'] & 0x10) != 0)]

df_error_rate_forward = get_error_rates(df_forward, ['C2T', 'G2A'])
df_error_rate_reverse = get_error_rates(df_reverse, ['C2T', 'G2A'])


d_lengths_count_forward = get_read_lengt_count(df_forward)
d_lengths_count_reverse = get_read_lengt_count(df_reverse)

d_strand_count = get_strand_count(df)


x=x

#%% =============================================================================
# 
# =============================================================================


from extra_functions import base2index, ACGT_names

def get_transmission_emission_matrices(df, max_len=30):
    
    A = {}
    E = {}
    
    ind_names = ['ref_'+s for s in ACGT_names[:-2]]
    col_names = ['obs_'+s for s in ACGT_names[:-2]]
    
    for i in range(1, max_len+1):
        
        a = np.zeros((len(ACGT_names[:-2]), len(ACGT_names[:-2])))
        for k, ref in enumerate(ACGT_names[:-2]):
            for l, ref_next in enumerate(ACGT_names[:-2]):
                mask_pos = (df['position']== i+1)
                mask_l = (df['ref_base'] == base2index[ref_next])
                mask_k = (df['prev_ref_base'] == base2index[ref])
                
                num = df[mask_pos & mask_k & mask_l]['counts'].sum()
                den = df[mask_pos & mask_k]['counts'].sum()
                
                if den == 0:
                    a[k,l] = 0
                else:
                    a[k,l] = num/den
        a = pd.DataFrame(a, index=ind_names, columns=col_names)
        A[i] = a 
            
            
        e = np.zeros((len(ACGT_names[:-2]), len(ACGT_names[:-2])))
        for k, ref in enumerate(ACGT_names[:-2]):
            for b, obs in enumerate(ACGT_names[:-2]):
                
                mask_pos = (df['position']== i)
                mask_k = (df['ref_base'] == base2index[ref])
                mask_b = (df['obs_base'] == base2index[obs])
                
                num = df[mask_pos & mask_k & mask_b]['counts'].sum()
                den = df[mask_pos & mask_k]['counts'].sum()
                
                if den == 0:
                    e[k,b] = 0
                else:
                    e[k,b] = num/den
        e = pd.DataFrame(e, index=ind_names, columns=col_names)
        E[i] = e
    
    return A, E


A_forward, E_forward = get_transmission_emission_matrices(df_forward)
A_reverse, E_reverse = get_transmission_emission_matrices(df_reverse)

x=x

#%% =============================================================================
# Error rate plots
# =============================================================================

if not is_linux() and do_plotting:
    
    names = df_error_rate_forward.columns
    
    fig_error, ax_error = plt.subplots(1, 2, figsize=(10, 6))
    ax_error = ax_error.flatten()
    
    for i, (ax, name) in enumerate(zip(ax_error, names)):
        
        for strand, df_error_rates in zip(['Forward', 'Reverse'], 
                                          [df_error_rate_forward, 
                                           df_error_rate_reverse]):
            x = df_error_rates.index
            y = df_error_rates.loc[:, name]
            ax.plot(x, y, '-', label=strand)
        
        ax.set(xlabel='Read Pos', ylabel='Error Rate', title=name,
               xlim=(0, 25), ylim=(0, 0.3))
        
        ax.legend()
        
    fig_error.tight_layout()
    if save_plots:
        fig_error.savefig(plot_prefix+'error_rates.pdf', dpi=600)
        plt.close('all')


# =============================================================================
# Bar chart of lengths
# =============================================================================
    
    
    len_max = max([max(d.keys()) for d in [d_lengths_count_forward, d_lengths_count_reverse]])
    len_min = min([min(d.keys()) for d in [d_lengths_count_forward, d_lengths_count_reverse]])
    
    fig_read_length, ax_read_length = plt.subplots(figsize=(10, 6))
    histrange = (len_min-1, len_max+1)
    width = 0.4
    keys_forward = np.fromiter(d_lengths_count_forward.keys(), dtype=int)
    keys_reverse = np.fromiter(d_lengths_count_reverse.keys(), dtype=int)
    
    ax_read_length.bar(keys_forward-width/2, d_lengths_count_forward.values(), width, label='Forward')
    ax_read_length.bar(keys_reverse+width/2, d_lengths_count_reverse.values(), width, label='Reverse')
    
    ax_read_length.set(xlabel='Read Lenght', ylabel='Counts', 
                       title='Counts of read lenghts',
                        xlim=histrange,
                       )
    
    ax_read_length.legend()
    fig_read_length.tight_layout()
    if save_plots:
        fig_read_length.savefig(plot_prefix+'read_lengths.pdf', dpi=600)
        plt.close('all')
        

# =============================================================================
# Bar chart of strand flag
# =============================================================================
    
    
    fig_strands, ax_strands = plt.subplots(figsize=(20, 6))
    
    keys = np.fromiter(d_strand_count.keys(), dtype=int)
    x = np.arange(len(d_strand_count))
    ax_strands.bar(x, d_strand_count.values())
    ax_strands.set(xlabel='Strand flag', ylabel='Counts', 
           title='Counts of strand flags')
    ax_strands.set_xticks(x)
    ax_strands.set_xticklabels(keys)
    
    fig_strands.tight_layout()
    if save_plots:
        fig_strands.savefig(plot_prefix+'strand_flags.pdf', dpi=600)
        plt.close('all')
    




#%% =============================================================================
# 
# =============================================================================


from sklearn.metrics import roc_curve, auc
from sklearn.preprocessing import label_binarize
from scipy import interp

def ROC_curve(y_test, y_pred_proba):

    n_classes = y_test.nunique()

    # Compute ROC curve and ROC area for each class
    fpr = {}
    tpr = {}
    roc_auc = {}
    for i in range(n_classes):
        fpr[i], tpr[i], _ = roc_curve(label_binarize(y_test, classes=[0, 1, 2, 3, 4])[:, i], y_pred_proba[:, i])
        roc_auc[i] = auc(fpr[i], tpr[i])
    
    # Compute micro-average ROC curve and ROC area
    fpr["micro"], tpr["micro"], _ = roc_curve(label_binarize(y_test, classes=[0, 1, 2, 3, 4]).ravel(), y_pred_proba.ravel())
    roc_auc["micro"] = auc(fpr["micro"], tpr["micro"])
    
    
    # Compute macro-average ROC curve and ROC area
    
    # First aggregate all false positive rates
    all_fpr = np.unique(np.concatenate([fpr[i] for i in range(n_classes)]))
    
    # Then interpolate all ROC curves at this points
    mean_tpr = np.zeros_like(all_fpr)
    for i in range(n_classes):
        mean_tpr += interp(all_fpr, fpr[i], tpr[i])
    
    # Finally average it and compute AUC
    mean_tpr /= n_classes
    
    fpr["macro"] = all_fpr
    tpr["macro"] = mean_tpr
    roc_auc["macro"] = auc(fpr["macro"], tpr["macro"])
    
    
    lw = 2
    # Plot all ROC curves
    fig, ax = plt.subplots(figsize=(12, 8))
    ax.plot(fpr["micro"], tpr["micro"],
             label=f'micro-average ROC curve (area = {roc_auc["micro"]:0.4f})',
             linestyle=':', linewidth=2)
    
    ax.plot(fpr["macro"], tpr["macro"],
             label=f'macro-average ROC curve (area = {roc_auc["macro"]:0.4f})',
             linestyle=':', linewidth=2)
    
    for i in range(n_classes):
        plt.plot(fpr[i], tpr[i], lw=lw,
                 label=f'ROC curve of class {i} (area = {roc_auc[i]:0.4f})')
    
    ax.plot([0, 1], [0, 1], 'k--', lw=lw)
    ax.set(xlim=(0.0, 1.0), ylim=(0.0, 1.05), xlabel='False Positive Rate', 
           ylabel='True Positive Rate', title='Some extension of Receiver operating characteristic to multi-class')
    ax.legend(loc="lower right")

    return fig, ax



#%% =============================================================================
# 
# =============================================================================


use_rescaled_weights = False


from sklearn.model_selection import train_test_split

y = df.loc[:, 'ref_base']
X = df.drop(columns='ref_base')

n_classes = y.nunique()

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

weights_unscaled = X['counts']

weights_scaled = (df.loc[:, 'obs_base'] == df.loc[:, 'ref_base'])
mean_weight = np.average(weights_scaled, weights=weights_unscaled)

mask = weights_scaled[:]
weights_scaled[mask] = mean_weight
weights_scaled[~mask] = 1/(1-mean_weight)

weights_scaled *= weights_unscaled
weights_scaled /= weights_scaled.sum() 
weights_scaled *= weights_unscaled.sum()


df['weights_unscaled'] = weights_unscaled
df['weights_scaled'] = weights_scaled


weights = weights_scaled if use_rescaled_weights else weights_unscaled
    
weights_train = weights.loc[X_train.index]
weights_test  = weights.loc[X_test.index]



#%% =============================================================================
#  Naive model: true base = observed base
# =============================================================================

from pycm import ConfusionMatrix #, online_help

# online_help("J")

y_pred_naive = X_test['obs_base']

cm = ConfusionMatrix(actual_vector=y_test.values, 
                     predict_vector=y_pred_naive.values, 
                     sample_weight=weights_test.values/1) # Create CM From Data
# cm.relabel(mapping={i: val for i, val in enumerate(ACGT_names[:-1])})

cm_matrix = pd.DataFrame(cm.matrix, dtype=int).values

fig, ax = plot_confusion_matrix(conf_mat=cm_matrix,
                                colorbar=True,
                                show_absolute=False,
                                show_normed=True,
                                label_names=ACGT_names[:-1],
                                str_format='.3f',
                                x_ticks_position='top')

df_eval_naive = precision_recall_bla_to_df(y_test, y_pred_naive, 
                                           sample_weight=weights_test)
    


#%% =============================================================================
# 
# =============================================================================


import lightgbm as lgb


params = {
	# ...
    'objective': 'multiclass',
    'num_class': n_classes,
    'metric':    'multi_logloss',
    'num_threads': cores,
    # "min_data_in_leaf":1,  
}


weight = (df.loc[:, 'obs_base'] == df.loc[:, 'ref_base'])
mean_weight = weight.mean()
mask = weight[:]
weight[mask] = mean_weight
weight[~mask] = 1/(1-mean_weight)
weight_train = weight.loc[X_train.index]

lgb_train = lgb.Dataset(data=X_train, label=y_train, categorical_feature=X.columns.to_list())
lgb_train = lgb.Dataset(data=X_train, label=y_train, categorical_feature=X.columns.to_list(), weight=weight_train)
# lgb_train = lgb.Dataset(data=pd.get_dummies(X_train, columns=X.columns.to_list()), label=y_train)


# cv = lgb.cv(lgb_params, 
#               lgb_train, 
#               num_boost_round=100, 
#               early_stopping_rounds=15,
#               stratified=False, 
#               verbose_eval=50)

model = lgb.train(params, lgb_train, num_boost_round=100, categorical_feature=X.columns.to_list())

y_pred_proba = model.predict(X_test)
y_pred = y_pred_proba.argmax(1) 



#%% =============================================================================
# 
# =============================================================================


from catboost import CatBoostClassifier, Pool


# test_data = catboost_pool = Pool(train_data, train_labels)


# model = CatBoostClassifier(iterations=2, 
#                            depth=2, 
#                            learning_rate=1, 
#                            loss_function='Logloss', 
#                            logging_level='Verbose')
# #train the model
# model.fit(train_data, train_labels)
# # make the prediction using the resulting model
# preds_class = model.predict(test_data)
# preds_proba = model.predict_proba(test_data)
# print("class = ", preds_class)
# print("proba = ", preds_proba)



params = {
    'iterations': 10,
    'learning_rate': 0.1,
    # 'eval_metric': 'Accuracy',
    'depth': 2,
    'loss_function': 'MultiClass',
    'random_seed': 42,
    'logging_level': 'Verbose',
    'cat_features': np.arange(X.shape[1]),
    # 'use_best_model': False
}

train_pool = Pool(X_train, y_train, cat_features=np.arange(X.shape[1]),
                    weight=weight_train,
                  )
# validate_pool = Pool(X_validation, y_validation, cat_features=categorical_features_indices)



# Initialize CatBoostClassifier
model = CatBoostClassifier(**params)
# Fit model
model.fit(train_pool, 
          # plot=True, 
          logging_level='Verbose',  # you can uncomment this for text output)
          )
# Get predicted classes
preds_class = model.predict(X_test)
# Get predicted probabilities for each class
preds_proba = model.predict_proba(X_test)
# Get predicted RawFormulaVal
preds_raw = model.predict(X_test, prediction_type='RawFormulaVal')                            

ROC_curve(y_test, preds_proba)




#%% =============================================================================
# 
# =============================================================================




fig, ax = ROC_curve(y_test, y_pred_proba)





# from sklearn.metrics import multilabel_confusion_matrix
# multilabel_confusion_matrix(y_true, y_pred, labels=range(n_classes))

