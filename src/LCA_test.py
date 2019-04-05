#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 16:06:04 2019

@author: michelsen
"""

import numpy as np
import pandas as pd
from tqdm import tqdm


file_LCA = '../data/raw/Mex_9_S50_L005_R1_001_holi_LCA.lca'
file_bam_in = '../data/raw/Mex_9_S50_L005_R1_001_holi.sort.bam'


    
def extract_int_name_type_from_phylo_comb(string):
    lst = string.split(':')
    for i, elem in enumerate(lst):
        if i==0:
            if elem.isdigit():
                continue
            else:
                print("first element in phylo_comb is not digit")
                assert False
        else:
            if elem.isdigit():
                p_int = int(lst[0])
                p_name = ':'.join(lst[1:i-1])
                p_type = lst[i-1]
                return p_int, p_name, p_type
    return -1, '-1', '-1'






res = []

with open(file_LCA, 'r') as f:
    header = f.readline()
    for line in f:
        line = line.split()
        name_seq_L_N = line[0]
        phylo_comb = ':'.join(line[1:])
        
        splitted = name_seq_L_N.split(':')
        name = ':'.join(splitted[:-3])
        seq, L, N = splitted[-3:]
        
        p = extract_int_name_type_from_phylo_comb(phylo_comb)        
        phylo_int, phylo_name, phylo_type = p

        res.append({'name': name, 
                    'L': int(L), 
                    'N': int(N),
                    'phylo_int': phylo_int,
                    'phylo_name': phylo_name,
                    'phylo_type': phylo_type,
                    'phylo_comb': phylo_comb, 
                    'seq': seq})
    
    
df = pd.DataFrame.from_dict(res)

# remove bad rows
df = df.loc[df['phylo_type'] != '-1']
# sort by N
df = df.sort_values('N', ascending=False)

N_phylos = 10
# get the N_pylos most observed phylos
phylo_int_top = (df.groupby('phylo_int')['N']
                   .sum()
                   .sort_values(ascending=False)
                   .iloc[:N_phylos]
                 )

# get the reads from the original df that are in the top10 most observed phylo
# and sort the dataframe
df_top = df.loc[df['phylo_int'].isin(phylo_int_top.index)].sort_values('N', ascending=False)

set_names_top = df_top.groupby('phylo_int')['name'].apply(set).to_dict()

name_to_phylo_int = {}
for key, vals in set_names_top.items():
    for val in vals:
        name_to_phylo_int[val] = key

phylo_int_to_str = {x: f'{z}:{y}' for x,y,z in 
                    df_top[['phylo_int', 'phylo_name', 'phylo_type']].values}

phylo_int_unique = df_top['phylo_int'].unique()


qname = 'D00708:299:CB6VWANXX:5:2307:8147:26950'
# qname2 = ':'.join(qname.split(':')[:-1])


_BAM_UNMAPPED  = 0x4
_BAM_SECONDARY = 0x100
_BAM_FAILED_QC = 0x200
_BAM_PCR_DUPE  = 0x400
_BAM_CHIMERIC  = 0x800

filtered_flags = (_BAM_UNMAPPED  | \
                  # _BAM_SECONDARY | \
                  _BAM_FAILED_QC | \
                  _BAM_PCR_DUPE  | \
                  _BAM_CHIMERIC)

def _read_bamfile(bamfile):
    for read in bamfile:
        if not (read.flag & filtered_flags):
            yield read


import pysam

# open SAM/BAM file
in_bam = pysam.AlignmentFile(file_bam_in)

N_reads = 619956074

it = tqdm(_read_bamfile(in_bam), total=N_reads)
            
# main loop

reads_matching = {}
# N_icounter = 10_000_000


from collections import Counter


from extra_functions import build_alignment_reference_seq
from extra_functions import correct_reverse_strans
from extra_functions import tidy_ML_seq_ref_data



use_tuple = True


from contextlib import ExitStack

filenames = {key: f'./tmp/phylo_int_{key}_corrected.txt' for key in phylo_int_unique}


with ExitStack() as stack:
    files = {key: stack.enter_context(open(filename, 'w')) for key, filename in filenames.items()}
    
    for icounter, read_bam in enumerate(it, start=1):
        
        phylo_int = name_to_phylo_int.get(read_bam.qname)
        
        if phylo_int is not None:
            
        
            read = read_bam.query_sequence
            cigar = read_bam.cigarstring
            flag = read_bam.flag
            md_tag = read_bam.get_tag('MD')
            tupl = (str(flag), cigar, read, md_tag)
            
            string_out = '\t'.join(tupl)
            f_out = files[phylo_int]
            print(string_out, file=f_out)
            
            
            # reads_matching.setdefault(phylo_int, []).append(tupl)
            
            # # print(f"Got one at counter={counter}")
            # reads_matching.setdefault(phylo_int, []).append(read_bam)
            
            
            # counter_reads_matching = reads_matching.setdefault(phylo_int, Counter())
            # tidy_ML_read_bam(read_bam, counter_reads_matching)
        

print({key: len(vals) for key, vals in reads_matching.items()})

cores = 6
force_rerun = False
do_remove_Ns = True
save_plots = True

import matplotlib.pyplot as plt
from extra_functions import (
                                 # is_linux,
                                 get_ML_res,
                                 get_error_rates,
                                 get_read_lengt_count,
                                  get_flag_count,
                                 # plot_confusion_matrix,
                                 remove_Ns,
                                 # precision_recall_bla_to_df,
                                 )


for key, filename in filenames.items():
    
    assert False
    
    plot_prefix = f"../figures/{filename}_plot_"
    
    df = get_ML_res(filename, cores, force_rerun, N_reads=None, N_splits=100)
    
    if do_remove_Ns:
        df = remove_Ns(df)
    
    df_forward = df.loc[((df['flag'] & 0x10) == 0)]
    df_reverse = df.loc[((df['flag'] & 0x10) != 0)]
    
    df_error_rate_forward = get_error_rates(df_forward, ['C2T', 'G2A'])
    df_error_rate_reverse = get_error_rates(df_reverse, ['C2T', 'G2A'])
    
    
    d_lengths_count_forward = get_read_lengt_count(df_forward)
    d_lengths_count_reverse = get_read_lengt_count(df_reverse)
    
    d_flag_count = get_flag_count(df)

# =============================================================================
# Bar chart of lengths
# =============================================================================
    
    names = df_error_rate_forward.columns
    
    fig_error, ax_error = plt.subplots(1, 2, figsize=(10, 6))
    ax_error = ax_error.flatten()
    
    for i, (ax, name) in enumerate(zip(ax_error, names)):
        
        for flag, df_error_rates in zip(['Forward', 'Reverse'], 
                                          [df_error_rate_forward, 
                                           df_error_rate_reverse]):
            x = df_error_rates.index
            y = df_error_rates.loc[:, name]
            ax.plot(x, y, '-', label=flag)
        
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
# Bar chart of flag flag
# =============================================================================
    
    
    fig_flags, ax_flags = plt.subplots(figsize=(20, 6))
    
    keys = np.fromiter(d_flag_count.keys(), dtype=int)
    x = np.arange(len(d_flag_count))
    ax_flags.bar(x, d_flag_count.values())
    ax_flags.set(xlabel='flag', ylabel='Counts', 
           title='Counts of flag')
    ax_flags.set_xticks(x)
    ax_flags.set_xticklabels(keys)
    
    fig_flags.tight_layout()
    if save_plots:
        fig_flags.savefig(plot_prefix+'flags.pdf', dpi=600)
        plt.close('all')
    












if use_tuple:

    import pickle
    with open('filename.pickle', 'wb') as handle:
        pickle.dump(reads_matching, handle)

# with open('filename.pickle', 'rb') as handle:
#     b = pickle.load(handle)

x=x



def worker(list_reads):
    
    counter_reads_matching = Counter()
    
    if use_tuple:
    
        for tupl in list_reads:
            flag, cigar, read, md_tag = tupl
            seq, ref = build_alignment_reference_seq(read, cigar, md_tag)
            is_reverse, ref, seq = correct_reverse_strans(flag, seq, ref)
            tidy_ML_seq_ref_data(ref, seq, flag, counter_reads_matching, only_mapDamage_vars=True)    
    
    # else:
    #     for read_bam in list_reads:
    #         read = read_bam.query_sequence
    #         cigar = read_bam.cigarstring
    #         flag = read_bam.flag
    #         md_tag = read_bam.get_tag('MD')
            
    #         seq, ref = build_alignment_reference_seq(read, cigar, md_tag)
    #         is_reverse, ref, seq = correct_reverse_strans(flag, seq, ref)
    #         tidy_ML_seq_ref_data(ref, seq, flag, counter_reads_matching, only_mapDamage_vars=True)

    
    return counter_reads_matching


import multiprocessing as mp
from functools import partial


counter_reads_matching_dict = {}
def combine_counters(counter_reads_matching, key):
    counter_reads_matching_dict[key] = counter_reads_matching


cores = 6

#init pool with cores
pool = mp.Pool(cores)
#create jobs
for key, list_reads in reads_matching.items():
    #http://blog.shenwei.me/python-multiprocessing-pool-difference-between-map-apply-map_async-apply_async/
    pool.apply_async(worker, 
                     args=(list_reads), #  pbar, i
                     callback=partial(combine_counters, key)) # pbar, i
pool.close() # Prevents any more tasks from being submitted to the pool
pool.join() # Wait for the worker processes to exit







