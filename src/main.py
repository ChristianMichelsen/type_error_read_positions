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
from src.extra_functions import (get_error_rates_dataframe, 
                                 fill_mismatch_matrix, is_linux)

save_plots = True
do_ancient = True
verbose = True
force_rerun = True
do_plotting = False


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

from extra_functions import (build_alignment_sequence, 
                             build_reference_sequence, 
                             # get_cigar_parts, 
                             comp, 
                             _read_txtfile)

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

filename_mismatch = file_processed_in.replace('corrected.txt', 'mismatch_results.pkl')


def main():
    
    if not Path(filename_mismatch).is_file() or force_rerun:
        
        d_mismatch_forward = {}
        d_mismatch_reverse = {}
        
        lengths_forward = []
        lengths_reverse = []
        
        list_strand = []
        
        print("Parsing MD-tags to get mismatch matrix: ", flush=True)
        
        with open(file_processed_in, 'r') as f_processed:
            for iline, line_txt in tqdm(enumerate(_read_txtfile(f_processed)), total=N_reads):
                
                strand_processed, cigar_processed, read_processed, md_tag_processed = line_txt
                # seq, cigar, md_tag = read_processed, cigar_processed, md_tag_processed
                
                seq_processed = build_alignment_sequence(read_processed, cigar_processed, md_tag_processed)
                ref_processed = build_reference_sequence(read_processed, cigar_processed, md_tag_processed)
                
                is_reverse = ((strand_processed & 0x10) != 0)
                if is_reverse:
                    ref_processed = comp(ref_processed)
                    seq_processed = comp(seq_processed)
                    
                L = len(seq_processed)
                
                list_strand.append(strand_processed)
                
                if not is_reverse:
                    fill_mismatch_matrix(ref_processed, seq_processed, d_mismatch_forward, verbose=verbose)
                    lengths_forward.append(L)
                else:
                    fill_mismatch_matrix(ref_processed, seq_processed, d_mismatch_reverse, verbose=verbose)
                    lengths_reverse.append(L)
        
                
        list_strand = np.array(list_strand)
        lengths_forward = np.array(lengths_forward)        
        lengths_reverse = np.array(lengths_reverse)
    
        print('Finished parsing the MD-tags, now saving the file')
    
        mismatch_results = [d_mismatch_forward, d_mismatch_reverse, 
                            list_strand, lengths_forward, lengths_reverse]
    
        with open(filename_mismatch, 'wb') as file:  
            # store the data as binary data stream
            pickle.dump(mismatch_results, file)
    
    else:
    
        with open(filename_mismatch, 'rb') as file:  
            # read the data as binary data stream
            p = pickle.load(file)
            d_mismatch_forward, d_mismatch_reverse, list_strand, lengths_forward, lengths_reverse = p

    return d_mismatch_forward, d_mismatch_reverse, list_strand, lengths_forward, lengths_reverse
    

d_mismatch_forward, d_mismatch_reverse, list_strand, lengths_forward, lengths_reverse = main()


df_error_rates_forward = get_error_rates_dataframe(d_mismatch_forward)
df_error_rates_reverse = get_error_rates_dataframe(d_mismatch_reverse)


#%% =============================================================================
# Error rate plots
# =============================================================================

if not is_linux() and do_plotting:
    
    names = [col for col in df_error_rates_forward.columns if not '_s' in col]
    
    fig_error, ax_error = plt.subplots(1, 3, figsize=(10, 6))
    ax_error = ax_error.flatten()
    
    for i, (ax, name) in enumerate(zip(ax_error, names)):
        
        for strand, df_error_rates in zip(['Forward', 'Reverse'], 
                                          [df_error_rates_forward, 
                                           df_error_rates_reverse]):
            x = df_error_rates.index
            y = df_error_rates.loc[:, name]
            # sy = df_error_rates.loc[:, name+'_s']
            # ax.errorbar(x, y, sy, fmt='-', label=strand)
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
    
    
    from collections import Counter
    c_forward = Counter(sorted(lengths_forward)[::-1])
    c_reverse = Counter(sorted(lengths_reverse)[::-1])
    
    
    fig_read_length, ax_read_length = plt.subplots(figsize=(10, 6))
    histrange = (lengths_forward.min()-1, lengths_forward.max()+1)
    width = 0.4
    keys = np.fromiter(c_forward.keys(), dtype=float)
    ax_read_length.bar(keys-width/2, c_forward.values(), width, label='Forward')
    ax_read_length.bar(keys+width/2, c_reverse.values(), width, label='Reverse')
    
    ax_read_length.set(xlim=histrange[::-1], xlabel='Read Lenght', ylabel='Counts', 
           title='Counts of read lenghts')
    
    ax_read_length.legend()
    fig_read_length.tight_layout()
    if save_plots:
        fig_read_length.savefig(plot_prefix+'read_lengths.pdf', dpi=600)
        plt.close('all')
        

# =============================================================================
# Bar chart of strand flag
# =============================================================================
    
    
    c_strands = Counter(sorted(list_strand))
    
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
    
