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
# import re
# from prepare_reads_v3 import (do_correct_file_if_not_exists, 
#                            get_corrected_filename, 
#                            # file_length,
#                            )

from src.data.prepare_reads import get_filename_and_lenght 
from src.extra_functions import (
                                 get_error_rates_dataframe, 
                                 is_linux,
                                  build_alignment_reference_seq,
                                  correct_reverse_strans,
                                  fill_results,
                                 _read_txtfile,
                                 calc_ACGT_content,
                                 MAX_LENGTH,
                                 mismatch_to_dataframe_collapsed,
                                 )


save_plots = True

do = 'ancient'
# do = 'modern'
do = 'gargamel'


verbose = True
force_rerun = True
do_plotting = True


cores = 1


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




plot_prefix = f"../figures/{filename.split('_')[0]}_plot_"


def append_corrected_to_filename(filename):
    return filename.replace('.txt', '_corrected.txt')

def get_filename_input(filename):
    return f'../data/raw/{filename}'
def get_filename_output(filename):
    return f'../data/processed/{append_corrected_to_filename(filename)}'
def get_filename_processed(filebame):
    return get_filename_output(filename)

file_processed_in = get_filename_processed(filename)
file_bam_in = get_filename_input(filename_bam)
ref_in = get_filename_input(refname)


_, N_reads = get_filename_and_lenght(filename)


#%% =============================================================================
# 
# =============================================================================

if False:
    
    # import pysam
    # from datetime import datetime
    from mapdamage_custom import mapDamage_main_test
    
    mapDamage_main_test(ref_in, file_bam_in, file_processed_in)



#%% =============================================================================
# #
# =============================================================================


#%% =============================================================================
# https://www.blopig.com/blog/2016/08/processing-large-files-using-python/
# =============================================================================


import multiprocessing as mp, os
from functools import partial

def save_d_res(d_res, filename_mismatch):
    max_l = max([max(d_res['lengths'][i].keys()) for i in ['+', '-']])
    for i in ['+', '-']:
        d_res['mismatch'][i] = d_res['mismatch'][i][:max_l, :, :]
    with open(filename_mismatch, 'wb') as file:  
        pickle.dump(d_res, file)
    return None

def load_d_res(filename_mismatch):
    with open(filename_mismatch, 'rb') as file:  
        d_res = pickle.load(file)    
    return d_res


from collections import Counter
def init_d():
    d = {'mismatch': {'+': np.zeros((MAX_LENGTH, 7,7), dtype=int), 
                      '-': np.zeros((MAX_LENGTH, 7,7), dtype=int)},
         'lengths':  {'+': Counter(), 
                      '-': Counter()}, 
         'strand': Counter(),
         }
    return d


def process_line(line_txt, d_res):
    strand, cigar, read, md_tag = line_txt
    seq, ref = build_alignment_reference_seq(read, cigar, md_tag)
    is_reverse, ref, seq = correct_reverse_strans(strand, seq, ref)
    fill_results(is_reverse, ref, seq, strand, d_res, verbose=True)
    return None


# =============================================================================
# 
# =============================================================================


ACGT_names = ['A', 'C', 'G', 'T', 'N', '0']
base2index = {val: i for i, val in enumerate(ACGT_names)}

strand2index = {'+': 1, '-': 0}
is_reverse2index = {True: 0, False: 1}

def ACGTN_correct_string(string):
    return ''.join(char if char in ACGT_names else 'N' for char in string)




def tidy_ML_seq_ref_data(ref, seq, is_reverse, res):
    
    seq = ACGTN_correct_string(seq)
    ref = ACGTN_correct_string(ref)
    
    for i, (s_ref, s_seq) in enumerate(zip(ref, seq)):
        if i == 0:
            prev_ref_base = '0'
        else:
            prev_ref_base = ref[i-1]
            
        # res.append((base2index[s_seq], base2index[prev_ref_base], i+1, 
                    # is_reverse2index[is_reverse], base2index[s_ref]))
        
        # res[((base2index[s_seq], base2index[prev_ref_base], i+1, 
        #             is_reverse2index[is_reverse], base2index[s_ref]))] += 1
        
        obs_base = base2index[s_seq]
        prev_ref_base = base2index[prev_ref_base]
        position = i+1
        strand = is_reverse2index[is_reverse]
        ref_base = base2index[s_ref]
        is_mismatch = (obs_base != ref_base)
        
        res[(obs_base, prev_ref_base, position, strand, ref_base, is_mismatch)] += 1
    
    
    
    return None
    

# def tidy_ML_seq_ref_data(ref, seq, is_reverse):
    
#     seq = ACGTN_correct_string(seq)
#     ref = ACGTN_correct_string(ref)
    
#     for i, (s_ref, s_seq) in enumerate(zip(ref, seq)):
#         if i == 0:
#             prev_ref_base = '0'
#         else:
#             prev_ref_base = ref[i-1]
            
#         yield (base2index[s_seq], base2index[prev_ref_base], i+1, 
#                     is_reverse2index[is_reverse], base2index[s_ref])
        


def process_line_tidy_ML(line_txt, res):
    strand, cigar, read, md_tag = line_txt
    seq, ref = build_alignment_reference_seq(read, cigar, md_tag)
    is_reverse, ref, seq = correct_reverse_strans(strand, seq, ref)
    # strand = '-' if is_reverse else '+'
    tidy_ML_seq_ref_data(ref, seq, is_reverse, res)
    # tidy_ML_seq_ref_data(ref, seq, strand, res)
    # return None


# =============================================================================
# 
# =============================================================================


def process_chunk(chunk_start, chunk_size, is_last):
    d_res = init_d()
    with open(file_processed_in) as f:
        f.seek(chunk_start)
        lines = f.read(chunk_size).splitlines()
        it = enumerate(_read_txtfile(lines))
        if is_last:
        # if True:
            it = tqdm(it, total=len(lines))
        for iline, line_txt in it:
            process_line(line_txt, d_res)
    return d_res


def chunkify(fname, cores):
    file_end = os.path.getsize(fname)
    if cores == 1:
        yield 0, file_end
    else:
        size = file_end // cores
        with open(fname, 'r') as f:
            chunkEnd = f.tell()
            while True:
                chunk_start = chunkEnd
                f.seek(chunk_start+size)
                f.readline()
                chunkEnd = f.tell()
                chunk_size = chunkEnd - chunk_start
                
                is_last = True if chunkEnd > file_end else False
                
                yield chunk_start, chunk_size, is_last
                if chunkEnd > file_end:
                    break


def collect_result(d_res, d_chunk):
    for i in ['+', '-']:
        d_res['mismatch'][i] += d_chunk['mismatch'][i]
        d_res['lengths'][i] += d_chunk['lengths'][i]
    d_res['strand'] += d_chunk['strand']
    return None



import h5py
# with h5py.File('several_datasets.hdf5', 'a') as f_out:
#    dset_int_1 = f_out.create_dataset('integers', (None, 5), dtype='i1')
#    dset.resize(20, axis=0)   # or dset.resize((20,1024))

# # with open('AE_data.dat', 'wb') as f_out:
# #     np.savetxt(f, [], header=header)
# #     for i in range(201):
# #         data = np.column_stack((x[i], y[i]))
# #         np.savetxt(f, data)
# #         f.flush()
# #         sleep(0.1)


import pickle
from pathlib import Path

filename_mismatch = file_processed_in.replace(f'corrected.txt', 
                                              f'mismatch_results_{cores}cores.pkl')
filename_mismatch_ML = filename_mismatch.replace('.pkl', '_ML.pkl')



header = ['obs_base', 'prev_ref_base', 'position', 'strand', 'ref_base', 'is_mismatch', 'counts']

do_mapDamage_analysis = True
do_ML_analysis = True


if do_mapDamage_analysis:
    
    if not Path(filename_mismatch).is_file() or force_rerun:
    
        d_res = init_d()
        
        if cores == 1:
            
            with open(file_processed_in, 'r') as f_processed:
                for iline, line_txt in tqdm(enumerate(_read_txtfile(f_processed)), 
                                            total=N_reads):
                    process_line(line_txt, d_res)
        else:
            
            #init objects
            pool = mp.Pool(cores)
            
            #create jobs
            for chunk_start, chunk_size, is_last in chunkify(file_processed_in, 10*cores):
                # print(chunk_start, chunk_size)
                #http://blog.shenwei.me/python-multiprocessing-pool-difference-between-map-apply-map_async-apply_async/
                pool.apply_async(process_chunk, args=(chunk_start, chunk_size, is_last), 
                                 callback=partial(collect_result, d_res))
            
            pool.close() # Prevents any more tasks from being submitted to the pool
            pool.join() # Wait for the worker processes to exit
            
        print('Finished parsing the MD-tags, now saving the file')
        save_d_res(d_res, filename_mismatch)

    
    else:
        
        print("Loading custom mapDamage results")
        d_res = load_d_res(filename_mismatch)



if do_ML_analysis:

    if not Path(filename_mismatch_ML).is_file() or force_rerun:
        
        print(f"Parsing MD-tags to get mismatch matrix using {cores} cores", flush=True)
        
        
        if cores == 1:
            
            # res = []
            ML_res = Counter()
        
            with open(file_processed_in, 'r') as f_processed:    
                # current_length = 0
                for iline, line_txt in tqdm(enumerate(_read_txtfile(f_processed)), 
                                            total=N_reads):
                    process_line_tidy_ML(line_txt, ML_res)
            
            # converts Counter dict to list of tuples that Pandas can read 
            df = []
            for key, val in dict(ML_res).items():
                df.append((*key, val))
            
            # create dataframe from list of tuples
            df = pd.DataFrame(df, columns=header)
            
            #downcast to unsigned integers (since always positive counts)
            for col in df:
                df.loc[:, col] = pd.to_numeric(df[col], downcast='unsigned')
                
            # save dataframe
            df.to_pickle(filename_mismatch_ML)
            
        else:
            print("More cores in ML analysis not implemented yet")
            assert False
                
                # #init objects
                # pool = mp.Pool(cores)
                
                # #create jobs
                # for chunk_start, chunk_size, is_last in chunkify(file_processed_in, 10*cores):
                #     # print(chunk_start, chunk_size)
                #     #http://blog.shenwei.me/python-multiprocessing-pool-difference-between-map-apply-map_async-apply_async/
                #     pool.apply_async(process_chunk, args=(chunk_start, chunk_size, is_last), 
                #                      callback=partial(collect_result, d_res))
                
                # pool.close() # Prevents any more tasks from being submitted to the pool
                # pool.join() # Wait for the worker processes to exit
                
        # print('Finished parsing the MD-tags, now saving the file')
        # save_d_res(d_res, filename_mismatch)


    # load files
    else:
        print("Loading 1 core ML dataframe")
        df = pd.read_pickle(filename_mismatch_ML)
        

x=x




df.loc[df['strand']==0, []]






#%% =============================================================================
# 
# =============================================================================





import numpy_indexed as npi
with h5py.File('several_datasets.hdf5', 'r') as f_in:
    print(f_in['data'].shape, flush=True)
    d = f_in['data'][:100_000_000]
    
    d2 = np.ones((d.shape[0], d.shape[1]+1), dtype='uint8')
    d2[:, :-1] = d
    
    is_mismatch = (d[:, 0] != d[:, 4])
    d2[is_mismatch, -1] = int(1/(is_mismatch.sum() / len(is_mismatch)))
    
    groupby = npi.group_by(d2)
    
    df = pd.DataFrame(groupby.unique, columns=header+['mismach_weight'], dtype=np.uint8)
    df['count'] = pd.to_numeric(groupby.count, downcast='unsigned')

df

d_res = load_d_res(filename_mismatch)



mismatch_forward = d_res['mismatch']['+']
mismatch_reverse = d_res['mismatch']['-']

df_mismatch_collapsed_forward = mismatch_to_dataframe_collapsed(mismatch_forward)
df_mismatch_collapsed_reverse = mismatch_to_dataframe_collapsed(mismatch_reverse)

df_error_rates_forward = get_error_rates_dataframe(mismatch_forward, max_len=30)
df_error_rates_reverse = get_error_rates_dataframe(mismatch_reverse, max_len=30)

ACGT_content = calc_ACGT_content(mismatch_forward, mismatch_reverse)




def compare_two_d_res_dicts(file1, file2):
    
    d_res1 = load_d_res(file1)
    d_res2 = load_d_res(file2)
    
    assert d_res1['strand'] == d_res2['strand']

    assert d_res1['lengths']['+'] == d_res2['lengths']['+']
    assert d_res1['lengths']['-'] == d_res2['lengths']['-']

    assert np.array_equal(d_res1['mismatch']['+'], d_res2['mismatch']['+'])
    assert np.array_equal(d_res1['mismatch']['-'], d_res2['mismatch']['-'])
    

other_file = (filename_mismatch.replace('1cores', '6cores') if cores == 1 
              else filename_mismatch.replace('6cores', '1cores'))
compare_two_d_res_dicts(filename_mismatch, other_file)
    

#%% =============================================================================
# Error rate plots
# =============================================================================

if not is_linux() and do_plotting:
    
    names = [col for col in df_error_rates_forward.columns if not '_s' in col]
    
    fig_error, ax_error = plt.subplots(1, 2, figsize=(10, 6))
    ax_error = ax_error.flatten()
    
    for i, (ax, name) in enumerate(zip(ax_error, names)):
        
        for strand, df_error_rates in zip(['Forward', 'Reverse'], 
                                          [df_error_rates_forward, 
                                           df_error_rates_reverse]):
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
    
    
    c_forward = d_res['lengths']['+']
    c_reverse = d_res['lengths']['-']
    len_max = max([max(d.keys()) for d in [c_forward, c_reverse]])
    len_min = min([min(d.keys()) for d in [c_forward, c_reverse]])
    
    
    fig_read_length, ax_read_length = plt.subplots(figsize=(10, 6))
    histrange = (len_min-1, len_max+1)
    width = 0.4
    keys_forward = np.fromiter(c_forward.keys(), dtype=int)
    keys_reverse = np.fromiter(c_reverse.keys(), dtype=int)
    
    ax_read_length.bar(keys_forward-width/2, c_forward.values(), width, label='Forward')
    ax_read_length.bar(keys_reverse+width/2, c_reverse.values(), width, label='Reverse')
    
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
    
    
    c_strands = d_res['strand']
    
    fig_strands, ax_strands = plt.subplots(figsize=(20, 6))
    
    keys = np.fromiter(c_strands.keys(), dtype=int)
    x = np.arange(len(c_strands))
    ax_strands.bar(x, c_strands.values())
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

from sklearn.model_selection import train_test_split
import lightgbm as lgb
from sklearn.metrics import roc_curve, auc
from sklearn.preprocessing import label_binarize
from scipy import interp


y = df.loc[:, 'ref_base']
X = df.drop(columns='ref_base')
# y = label_binarize(y, classes=[0, 1, 2, 3, 4])



n_classes = y.nunique()


X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

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


fig, ax = ROC_curve(y_test, y_pred_proba)





# from sklearn.metrics import multilabel_confusion_matrix
# multilabel_confusion_matrix(y_true, y_pred, labels=range(n_classes))

