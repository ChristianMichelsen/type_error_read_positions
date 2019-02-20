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
from src.extra_functions import parse_md_tag, get_error_rates_dataframe, make_reference

save_plots = False

filename = 'NA12400_error_test.txt' # modern dana
filename = 'ESW_LRUE_MA2621_error_test.txt' # ancient data


ACGT_names = ['A', 'C', 'G', 'T', 'N']
base2index = {val: i for i, val in enumerate(ACGT_names)}

plot_prefix = f"../figures/{filename.split('_')[0]}_plot_"

filename, file_len = get_filename_and_lenght(filename)


#%% =============================================================================
# #
# =============================================================================


d_mismatch_forward = {}
d_mismatch_reverse = {}

lengths_forward = []
lengths_reverse = []

list_strand = []


all_seqs = []

print("Parsing MD-tags to get mismatch matrix: ", flush=True)

with open(filename, 'r') as f:
    for iline, line in tqdm(enumerate(f), total=file_len):
        
        strand, cigar, read, md_tag = line.split()
        ref, seq = make_reference(read, md_tag, cigar)

        if 'I' in cigar:
            assert False

        if strand == '0':
            L = parse_md_tag(read, md_tag, cigar, strand, d_mismatch_forward)
            lengths_forward.append(L)
             
        elif strand == '16':
            L = parse_md_tag(read, md_tag, cigar, strand, d_mismatch_reverse)
            lengths_reverse.append(L)
        else:
            print(strand)
            assert False
            
            
        all_seqs.append(read)
        list_strand.append(int(strand))
        
        
list_strand = np.array(list_strand)
lengths_forward = np.array(lengths_forward)        
lengths_reverse = np.array(lengths_reverse)        

df_error_rates_forward = get_error_rates_dataframe(d_mismatch_forward)
df_error_rates_reverse = get_error_rates_dataframe(d_mismatch_reverse)


#%% =============================================================================
# Error rate plots
# =============================================================================

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
    
    ax.set(xlabel='Read Pos', ylabel='Error Rate', title=name)
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





