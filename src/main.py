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
do_ancient = True
verbose = True
force_rerun = False
do_plotting = True

cores = 6


if do_ancient:
    print("\nRunning on ancient DNA")
    filename = 'ESW_LRUE_MA2621_error_test.txt' # ancient data
    filename_bam = 'ESW_LRUE_MA2621_L1_S1_S83_L006.sort.rmdup.bam'
else:
    print("\nRunning on modern DNA")
    filename = 'NA12400_error_test.txt' # modern dna
    # filename = 'NA12400_small_bam_1_000_000_error_test.txt' # modern dna
    filename_bam = 'NA12400_small_bam_1_000_000.bam'
refname = 'hs37d5.fa'


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


import pickle
from pathlib import Path

filename_mismatch = file_processed_in.replace(f'corrected.txt', 
                                              f'mismatch_results_{cores}cores.pkl')



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


def process_chunk(chunk_start, chunk_size):
    d_res = init_d()
    with open(file_processed_in) as f:
        f.seek(chunk_start)
        lines = f.read(chunk_size).splitlines()
        it = enumerate(_read_txtfile(lines))
        if chunk_start == 0:
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
                yield chunk_start, chunk_size
                if chunkEnd > file_end:
                    break


def collect_result(d_res, d_chunk):
    for i in ['+', '-']:
        d_res['mismatch'][i] += d_chunk['mismatch'][i]
        d_res['lengths'][i] += d_chunk['lengths'][i]
    d_res['strand'] += d_chunk['strand']
    return None



if not Path(filename_mismatch).is_file() or force_rerun:
    
    print(f"Parsing MD-tags to get mismatch matrix using {cores} cores", flush=True)
    
    #init dict to store results in
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
        for chunk_start, chunk_size in chunkify(file_processed_in, 10*cores):
            # print(chunk_start, chunk_size)
            #http://blog.shenwei.me/python-multiprocessing-pool-difference-between-map-apply-map_async-apply_async/
            pool.apply_async(process_chunk, args=(chunk_start, chunk_size), 
                             callback=partial(collect_result, d_res))
        
        pool.close() # Prevents any more tasks from being submitted to the pool
        pool.join() # Wait for the worker processes to exit
        
    print('Finished parsing the MD-tags, now saving the file')
    save_d_res(d_res, filename_mismatch)


d_res = load_d_res(filename_mismatch)



mismatch_forward = d_res['mismatch']['+']
mismatch_reverse = d_res['mismatch']['-']

df_mismatch_collapsed_forward = mismatch_to_dataframe_collapsed(mismatch_forward)
df_mismatch_collapsed_reverse = mismatch_to_dataframe_collapsed(mismatch_reverse)

df_error_rates_forward = get_error_rates_dataframe(mismatch_forward, max_len=50)
df_error_rates_reverse = get_error_rates_dataframe(mismatch_reverse, max_len=50)

ACGT_content = calc_ACGT_content(mismatch_forward, mismatch_reverse)


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
               xlim=(1, 25), ylim=(0, 0.3))
        
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
    
