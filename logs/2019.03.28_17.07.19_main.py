#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 18 12:15:39 2018

@author: michelsen
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from tqdm import tqdm, trange

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
    print("\nRunning on ancient DNA", flush=True)
    filename = 'ESW_LRUE_MA2621_error_test.txt' # ancient data
    filename_bam = 'ESW_LRUE_MA2621_L1_S1_S83_L006.sort.rmdup.bam'
    refname = 'hs37d5.fa'
elif do == 'modern':
    print("\nRunning on modern DNA", flush=True)
    filename = 'NA12400_error_test.txt' # modern dna
    # filename = 'NA12400_small_bam_1_000_000_error_test.txt' # modern dna
    filename_bam = 'NA12400_small_bam_1_000_000.bam'
    refname = 'hs37d5.fa'
elif do == 'gargamel':
    print("\nRunning on gargamel simulated DNA", flush=True)
    filename = 'gargamel_0.03_0.4_0.01_0.5_s.txt'
    filename_bam = 'gargamel_0.03_0.4_0.01_0.5_s.bam'
    refname = 'horse_chrom31.fa'


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


x=x




if do_remove_Ns:
    df = remove_Ns(df)

df_forward = df.loc[((df['strand'] & 0x10) == 0)]
df_reverse = df.loc[((df['strand'] & 0x10) != 0)]

df_error_rate_forward = get_error_rates(df_forward, ['C2T', 'G2A'])
df_error_rate_reverse = get_error_rates(df_reverse, ['C2T', 'G2A'])


d_lengths_count_forward = get_read_lengt_count(df_forward)
d_lengths_count_reverse = get_read_lengt_count(df_reverse)

d_strand_count = get_strand_count(df)


def get_df_mismatches_only(df):
    return df.loc[df['ref_base'] != df['obs_base']]

df_mismatches = get_df_mismatches_only(df)



import itertools
from scipy.sparse import csr_matrix
def df_to_sparse_csv_matrix(df, d_merge_cols, min_counts=None):
    
    """ 
        Input: 
        
            a dataframe and a dict with 
            keys: list of columns in dataframe to use 
            values: iterator containing the unique values of the column
        
        Returns: 
        
            Sparse array with dimensions (N, M) where N is the number
            of total counts in df['counts'] and M is the total number of
            unique columns. 
    """
    
    if min_counts is None:
        mask_min_counts = np.ones_like(df['counts'], dtype=bool)
    elif 0 <= min_counts <= 1:
        quantile = df['counts'].quantile(min_counts)
        print(f"min_counts = {min_counts:.1%} corresponds to min_counts = {quantile}")
        mask_min_counts = (quantile <= df['counts'])
    else:
        mask_min_counts = (min_counts <= df['counts'])
    
    
    d_tuple_to_index = {}
    for index, dimension_tuple in enumerate(itertools.product(*d_merge_cols.values())):
        d_tuple_to_index[dimension_tuple] = index

    def tuple_to_index(dimension_tuple):
        """ 
            Get the index of the dimension tuple descrived by the keys of 
            the input dict d_merge_cols. If not in index, return '-1'.
        """
        try:
            return d_tuple_to_index[tuple(dimension_tuple)]
        except KeyError:
            return -1
    
    # get the indices
    indices = np.apply_along_axis(tuple_to_index, 1, df[d_merge_cols.keys()])
    #mask out any indices equal to -1
    mask_indices = (indices != -1)
    mask = np.logical_and(mask_indices, mask_min_counts)
    # get the index of each column by the mask_indicesed indices times the counts
    col_ind = np.repeat(indices[mask], df['counts'][mask])
    # We have one sample pr. row, thus a range from 0, 1, ..., sum(count)    
    row_ind = np.arange(df['counts'][mask].sum())
    # We have a single data point in each sample (ie. row)
    data = np.ones(df['counts'][mask].sum())
    # Create a sparse matrix of the data, and rows and columns indices
    csr = csr_matrix((data, (row_ind, col_ind)), dtype=int)
    return csr


df_gargamel1 = get_ML_res('../data/processed/gargamel_0.03_0.4_0.01_0.1_s_corrected.txt',         cores, force_rerun, N_reads, N_splits=100)
df_gargamel3 = get_ML_res('../data/processed/gargamel_0.03_0.4_0.01_0.3_s_corrected.txt',         cores, force_rerun, N_reads, N_splits=100)
df_gargamel5 = get_ML_res('../data/processed/gargamel_0.03_0.4_0.01_0.5_s_corrected.txt',         cores, force_rerun, N_reads, N_splits=100)
df_modern =   get_ML_res('../data/processed/NA12400_error_test_corrected.txt',         cores, force_rerun, N_reads, N_splits=100)
df_ancient =  get_ML_res('../data/processed/ESW_LRUE_MA2621_error_test_corrected.txt', cores, force_rerun, N_reads, N_splits=100)




pos_max = 40 # TODO: Fix this maximum
d_merge_cols = {'obs_base': range(4), 
                'ref_base': range(4), 
                'position': range(1, pos_max+1),
                }

# csr_gargamel1 = df_to_sparse_csv_matrix(get_df_mismatches_only(df_gargamel1), d_merge_cols, min_counts=10)
# csr_gargamel3 = df_to_sparse_csv_matrix(get_df_mismatches_only(df_gargamel3), d_merge_cols, min_counts=10)
# csr_gargamel5 = df_to_sparse_csv_matrix(get_df_mismatches_only(df_gargamel5), d_merge_cols, min_counts=10)
# csr_modern = df_to_sparse_csv_matrix(get_df_mismatches_only(df_modern), d_merge_cols, min_counts=10)
# csr_ancient = df_to_sparse_csv_matrix(get_df_mismatches_only(df_ancient), d_merge_cols, min_counts=10)


# import scipy.sparse as sp
# csr = sp.vstack((csr_modern, csr_ancient), format='csr')

# from sklearn.decomposition import TruncatedSVD
# svd = TruncatedSVD(n_components=5, n_iter=7, random_state=42)
# svd.fit(csr)  


def dataframe_to_word_frequency(df, d_merge_cols):
    df2 = pd.DataFrame(df['counts'], index=df.index)
    df2['merged_cols'] = list(zip(*[df[i] for i in d_merge_cols.keys()]))
    s_mismatches = (df2.groupby('merged_cols')['counts']
                       .sum()
                       # .reindex(itertools.product(*d_merge_cols.values()))
                       # .fillna(0)
                       # .astype(int)
                    )
    
    return s_mismatches


# s_mismatches_gargamel = dataframe_to_word_frequency(get_df_mismatches_only(df_gargamel), d_merge_cols)
s_mismatches_gargamel1 = dataframe_to_word_frequency(get_df_mismatches_only(df_gargamel1), d_merge_cols)
s_mismatches_gargamel3 = dataframe_to_word_frequency(get_df_mismatches_only(df_gargamel3), d_merge_cols)
s_mismatches_gargamel5 = dataframe_to_word_frequency(get_df_mismatches_only(df_gargamel5), d_merge_cols)
s_mismatches_modern = dataframe_to_word_frequency(get_df_mismatches_only(df_modern), d_merge_cols)
s_mismatches_ancient = dataframe_to_word_frequency(get_df_mismatches_only(df_ancient), d_merge_cols)

x=x


from collections import Counter
def generate_bootstrap_sample(s_mismatch, N_samples=None):
    
    if N_samples is None:
        N_samples = s_mismatch.sum()
    elif 0 <= N_samples <= 1:
        N_samples = int(N_samples * s_mismatch.sum())
    
    random_indices = np.random.choice(s_mismatch.index, 
                                      size=N_samples, 
                                      replace=True, 
                                      p=s_mismatch/s_mismatch.sum())
    counter = Counter(random_indices)
    s = pd.DataFrame.from_dict(counter, orient='index', columns=['counts'])
    return s
    

N_bootstraps = 10
N_samples = 100_000
N_samples = 0.01
do_normalize = True
K_topics = 3

# s_bootstraps_gargamel = [generate_bootstrap_sample(
#                          s_mismatches_gargamel, N_samples=N_samples) 
#                          for _ in trange(N_bootstraps)]

s_bootstraps_gargamel135 = [s_mismatches_gargamel1, s_mismatches_gargamel3, s_mismatches_gargamel5]

s_bootstraps_modern = [generate_bootstrap_sample(
                       s_mismatches_modern, N_samples=N_samples) 
                       for _ in trange(N_bootstraps)]
s_bootstraps_modern = [s_mismatches_modern]


s_bootstraps_ancient = [generate_bootstrap_sample(
                        s_mismatches_ancient, N_samples=N_samples) 
                        for _ in trange(N_bootstraps)]
s_bootstraps_ancient = [s_mismatches_ancient]


def join_bootstrap_lists(*lists):
    """ Both merges lists and counts the elements in each list """
    counts = ([len(l) for l in lists])
    joined = sum(lists, [])
    return joined, counts


s_bootstraps_all, N_counts = join_bootstrap_lists(s_bootstraps_gargamel135, s_bootstraps_modern, s_bootstraps_ancient)
df_concat = (pd.concat(s_bootstraps_all, axis=1, join='outer', sort=True)
               .fillna(0)
               .T
               .astype(int)
               .reset_index(drop=True))

if do_normalize:
    max_count = df_concat.sum(1).max()
    df_concat = df_concat.multiply(max_count//df_concat.sum(1), axis=0)
    




do_pca_reduction = False
if do_pca_reduction:
    df_concat = df_concat

N_gargamel = len(s_bootstraps_gargamel135)
N_modern = len(s_bootstraps_modern)
N_all = len(s_bootstraps_all)

X_all = df_concat
X = df_concat.iloc[:N_gargamel+N_modern, :]
X_ancient = df_concat.iloc[N_gargamel+N_modern:, :]

y = np.zeros(len(X), dtype=int)
y[N_gargamel:] = 1


y_all = np.zeros(len(X_all), dtype=int)
y_all[N_gargamel:N_gargamel+N_modern] = 1
y_all[N_gargamel+N_modern:] = 2



#%% =============================================================================
# 
# =============================================================================


do_tsne = False

if do_tsne:
    
    # import seaborn as sns; sns.set()
    from fast_tsne import fast_tsne
    
    # # 10 nice colors
    col = np.array(['#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99', '#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a'])
    
    
    X_sne = X
    X_sne = X_all
    
    y_sne = y
    y_sne = y_all
    
    
    perplexity = 50
    
    # if(N - 1 < 3 * perplexity) { printf("Perplexity too large for the number of data points!\n"); exit(1); 
    if (len(X_sne) - 1) < (3 * perplexity):
        print(f"\nperplexity is too large at {perplexity}")
        perplexity = (len(X_sne) - 1) // 3 
        print(f"Using new value of {perplexity}\n")
        
    Z = fast_tsne(X_sne, perplexity=perplexity, seed=42, map_dims=2, max_iter=1000, nbody_algo='FFT', nthreads=None, return_loss=False)
    
    
    plt.figure(figsize=(10,10))
    plt.scatter(Z[:,0], Z[:,1], c=y_sne, s=10)
    plt.tight_layout()



#%% =============================================================================
# 
# =============================================================================


def plot_topic_distribution(prob_topic, topics=None, document_names=None, 
                            figsize=(10, 10), text_pos=-0.05, **kwargs):
    
    N, K = prob_topic.shape
    
    if topics is None:
        topics = [f'Topic {i}' for i in range(K)]

    df_prob_topic = pd.DataFrame(prob_topic, columns=topics)
    
    fig, ax = plt.subplots(figsize=figsize)
    df_prob_topic.plot.barh(stacked=True, width=1, ax=ax)
    ax.set(ylim=(0-0.5, N-0.5))
    
    ax.legend(loc='lower left', bbox_to_anchor= (0.0, -0.075), ncol=2, 
            borderaxespad=0, frameon=False)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    
    ax.get_yaxis().set_ticks([])
    
    ax.invert_yaxis()
    
    
    if document_names is not None:
        for key, val in document_names.items():
            mean = np.mean(val) 
            ax.errorbar(text_pos, mean-0.5, (max(val)-mean)*0.99, fmt='', ecolor='k', capsize=5)
            ax.text(text_pos, mean-0.5, key, fontsize=12, 
                    horizontalalignment='center', verticalalignment='center',
                    bbox=dict(facecolor='white', alpha=1))
    
    ax.set(**kwargs)
    fig.tight_layout()
    
    return fig, ax

    
d = {'Gargamel': [0, np.cumsum(N_counts)[0]],
     'Modern': [np.cumsum(N_counts)[0], np.cumsum(N_counts)[1]],
     'Ancient': [np.cumsum(N_counts)[1], np.cumsum(N_counts)[2]]}




if K_topics == 2:
    topics=['Modern', 'Gargamel']
    topics=['1', '2']
else:
    topics=['Modern', 'Ancient', 'Gargamel']
    topics=['1', '2', '3']




# 
#%% =============================================================================
# 
# =============================================================================


from sklearn.decomposition import LatentDirichletAllocation


no_topics = K_topics

# Run LDA
lda = LatentDirichletAllocation(n_components=no_topics, 
                                max_iter=100, 
                                learning_method='online', 
                                learning_offset=50., 
                                random_state=0,
                                n_jobs=-1,
                                )
lda = lda.fit(X_all)
prob_topic_all_lda = lda.transform(X_all)
lda.components_.round(3)
lda_component_normed = lda.components_ / lda.components_.sum(axis=1)[:, np.newaxis]


# topic_models_components = pd.DataFrame({
#                                         'lda_topic0': lda_component_normed[0, :],
#                                         'lda_topic1': lda_component_normed[1, :],
#                                         }, index=X.columns)

    
plot_topic_distribution(prob_topic_all_lda, 
                        # topics=['Modern', 'Ancient', 'Gargamel'],
                        # topics=['Modern', 'Gargamel'],
                        topics=topics,
                        document_names=d,
                        title=f'Topic Distribution with K={K_topics} from sklearn LDA, do_normalize={do_normalize}')



#%% =============================================================================
# 
# =============================================================================
   

from slda.topic_models import LDA


# Estimate parameters
K = K_topics
alpha = np.ones(K)
beta = np.repeat(1/X.shape[1], X.shape[1])
n_iter = 50

lda2 = LDA(K, alpha, beta, n_iter, seed=42)

lda2.fit(X_all.values)
prob_topic_all_lda2 = lda2.transform(X_all.values)
results2_phi = lda2.phi 
results2_phi_normed = results2_phi / results2_phi.sum(axis=1)[:, np.newaxis]
results2_theta = lda2.theta 
    

plot_topic_distribution(prob_topic_all_lda2, 
                        # topics=['Modern', 'Ancient', 'Gargamel'],
                        # topics=['Modern', 'Gargamel'],
                        topics=topics,
                        document_names=d,
                        title=f'Topic Distribution with K={K_topics} from slda LDA, do_normalize={do_normalize}')
    


#%% =============================================================================
# 
# =============================================================================

x=x

from slda.topic_models import SLDA

K = K_topics
alpha = np.ones(K)
beta = np.repeat(1/X.shape[1], X.shape[1])
mu = 0
nu2 = 5
sigma2 = 1
n_iter = 100

slda = SLDA(K, alpha, beta, mu, nu2, sigma2, n_iter, seed=42)

slda.fit(X.values, y)
prob_topic_all_slda = slda.transform(X_all.values)
results_phi = slda.phi 
results_phi_normed = results_phi / results_phi.sum(axis=1)[:, np.newaxis]
results_theta = slda.theta 
results_eta = slda.eta 


plot_topic_distribution(prob_topic_all_slda, 
                        # topics=['Modern', 'Gargamel', 'Ancient'],
                        topics=topics,
                        document_names=d,
                        title=f'Topic Distribution with K={K_topics} from slda SLDA')
    


#%% =============================================================================
# 
# =============================================================================

x=x


import gensim
from gensim.test.utils import common_texts
from gensim.corpora.dictionary import Dictionary
from gensim.models.ldamulticore import LdaMulticore


class Bow2Corpus(object):
    """
    Treat dense numpy array as a sparse, streamed gensim corpus.

    No data copy is made (changes to the underlying matrix imply changes in the
    corpus).

    This is the mirror function to `corpus2dense`
    
    documents_columns: 
        Documents in dense represented as columns, as opposed to rows.

    """
    def __init__(self, dense, documents_columns=False):
        if documents_columns:
            self.dense = dense.T
        else:
            self.dense = dense


    def __iter__(self):
        for doc in self.dense:
            yield full2sparse(doc.flat)

    def __len__(self):
        return len(self.dense)


def full2sparse(vec, eps=1e-9):
    """
    Convert a dense numpy array into the sparse document format (sequence of 2-tuples).

    Values of magnitude < `eps` are treated as zero (ignored).

    This is the mirror function to `sparse2full`.

    """
    vec = np.asarray(vec, dtype=float)
    nnz = np.nonzero(abs(vec) > eps)[0]
    return list(zip(nnz, vec.take(nnz)))




common_dictionary = Dictionary(common_texts)
common_corpus = [common_dictionary.doc2bow(text) for text in common_texts]
lda_gensim = gensim.models.ldamodel.LdaModel(csr, num_topics=2)

bla = np.array([[1, 5], [0, 2], [3, 1]])

for row in Bow2Corpus(bla):
    print(row)


#%% =============================================================================
# 
# =============================================================================


d = {}
for (i,j), name_type in zip([(0, N_bootstraps), (N_bootstraps, 2*N_bootstraps), (2*N_bootstraps, 3*N_bootstraps)], 
                ['Gargam', 'Modern', 'Ancient']):
    
    for prob_topic, name_model in zip([prob_topic_all_lda, prob_topic_all_lda2, prob_topic_all_slda], 
                          ['Scikit LDA', 'slda LDA', 'slda SLDA']):
        mu = prob_topic[i:j].mean(0)
        sigma = prob_topic[i:j].std(0)
        
        d[(name_type, name_model)] = f"Topic 0: {mu[0]:.2%} +/- {sigma[0]:.2%}, Topic 1: {mu[1]:.2%} +/- {sigma[1]:.2%}"

df_comparison = pd.Series(d).unstack()


x=x

#%% =============================================================================
#     
# =============================================================================

# from sklearn.datasets import load_digits
# d = load_digits()
# with open("digits.npy", "wb") as f:
#     np.save(f, d.data.astype(np.int32).T)
#     np.save(f, d.target.astype(np.int32))

# from subprocess import call
# call('fslda train digits.npy model.npy')


# import numpy as np
# with open("model.npy") as f:
#     alpha = np.load(f)
#     beta = np.load(f)
#     eta = np.load(f)



   

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

