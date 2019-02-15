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
import re
from prepare_reads_v3 import (do_correct_file_if_not_exists, 
                           get_corrected_filename, 
                           # file_length,
                           )


save_plots = True

filename = 'ESW_LRUE_MA2621_error_test.txt' # ancient data
filename = 'NA12400_error_test.txt' # modern dana

file_len = do_correct_file_if_not_exists(filename)
if not "corrected" in filename:
    filename = get_corrected_filename(filename)

ACGT_names = ['A', 'C', 'G', 'T', 'N']
base2index = {val: i for i, val in enumerate(ACGT_names)}


plot_prefix = filename.split('_')[0] + '_plot_'


#%%


def seq_to_match(seq, read_pos, N_soft_clipped_beginning, d_mismatch):
    
    mask = any([s in ACGT_names for s in seq])
    if not mask:
        print(seq)
        assert False

    counter = read_pos + N_soft_clipped_beginning
    for base in seq:
        index = base2index[base]
        dict_to_mismatch(d_mismatch, counter)[index, index] += 1
        counter += 1
        
    return None
    

def seq_to_mismatch(seq, char, read_pos, N_soft_clipped_beginning, d_mismatch):
    
    assert len(char) == 1
    index_ref = base2index[char]
    index_seq = base2index[seq[read_pos]]
    dict_to_mismatch(d_mismatch, read_pos+N_soft_clipped_beginning)[index_ref, index_seq] += 1
    return None

cigar_implemented = list('MISD')
cigar_all_implemented = list('MIDNSHP=X')
cigar_not_implemented = [c for c in cigar_all_implemented if not c in cigar_implemented]


def find_all(a_str, sub):
    start = 0
    while True:
        start = a_str.find(sub, start)
        if start == -1: 
            return
        yield start
        start += len(sub) # use start += 1 to find overlapping matches



def parse_md_tag(seq, md_tag, cigar, strand, d_mismatch):
    
    L = len(seq)
    
    read_pos = 0
    N_inserts = 0
    N_soft_clipped_beginning = 0
    N_soft_clipped_end = 0
    
    
    if 'I' in cigar:
        cigar_split = list(filter(None, re.split(r'(\d+)', cigar)))
        for i, split in enumerate(cigar_split):
            if split == 'I':
                N_inserts += int(cigar_split[i-1])
        
    if 'S' in cigar:
        cigar_split = list(filter(None, re.split(r'(\d+)', cigar)))
        
        if cigar_split[1] == 'S' or cigar_split[-1] == 'S':
            if cigar_split[1] == 'S':
                N_soft_clipped_beginning = int(cigar_split[0])
                seq = seq[N_soft_clipped_beginning:]
            if cigar_split[-1] == 'S':
                N_soft_clipped_end = int(cigar_split[-2])
                seq = seq[:-N_soft_clipped_end]
        else:
            print("Soft clipping in the middle not implemented yet")
            assert False 

    if any([c in cigar_not_implemented for c in cigar]):
        print(f'Cigar character not implemented yet, {cigar}')
        assert False
        
    
    for md_split in filter(None, re.split(r'(\d+)', md_tag)):
        if md_split == '0':
            continue
        elif md_split.isdigit():
            i = int(md_split)
            seq_to_match(seq[read_pos:read_pos+i], read_pos, N_soft_clipped_beginning, d_mismatch)
            read_pos += i

        # ignore inserts for now            
        elif md_split[0] == '^':
            continue

        else:
            for char in md_split:
                if char.upper() in ACGT_names:
                    seq_to_mismatch(seq, char, read_pos, N_soft_clipped_beginning, d_mismatch)
                    read_pos += 1
                else: 
                    print(f"The characters {char} is not recognized")
                    assert False
    
    assert read_pos == len(seq) - N_inserts
    
    return L 


def gen_zero_matrix(ACGT_names):
    return np.zeros((len(ACGT_names), len(ACGT_names)), dtype=int)


def dict_to_mismatch(d_mismatch, read_pos):
    read_pos_1_indexed = read_pos + 1
    if not read_pos_1_indexed in d_mismatch:
        d_mismatch[read_pos_1_indexed] = gen_zero_matrix(ACGT_names)
    return d_mismatch[read_pos_1_indexed]



#%% =============================================================================
# #
# =============================================================================


d_mismatch_forward = {}
d_mismatch_reverse = {}

lengths_forward = []
lengths_reverse = []

list_strand = []

with open(filename, 'r') as f:
    for iline, line in tqdm(enumerate(f), total=file_len):
        parts = line.split()
        
        strand, cigar, seq, md_tag = parts
        
        list_strand.append(int(strand))
        
        if strand == '0':
            L = parse_md_tag(seq, md_tag, cigar, strand, d_mismatch_forward)
            lengths_forward.append(L)
             
        elif strand == '16':
            L = parse_md_tag(seq, md_tag, cigar, strand, d_mismatch_reverse)
            lengths_reverse.append(L)
        else:
            continue
        
list_strand = np.array(list_strand)
lengths_forward = np.array(lengths_forward)        
lengths_reverse = np.array(lengths_reverse)        




#%% =============================================================================
# 
# =============================================================================



def mismatch_to_error_rate(mismatch):
    diag_sum = np.diag(mismatch[:4, :4]).sum()
    sum_all = mismatch[:4, :4].sum()
    
    e_all = 1 - diag_sum / sum_all
    se_all = np.sqrt(e_all*(1-e_all)/sum_all)
    
    e_C2T = mismatch[base2index['C'], base2index['T']] / sum_all
    se_C2T = np.sqrt(e_C2T*(1-e_C2T)/sum_all)
    
    e_G2A = mismatch[base2index['G'], base2index['A']] / sum_all
    se_G2A = np.sqrt(e_G2A*(1-e_G2A)/sum_all)
    
    return (e_all, se_all, e_C2T, se_C2T, e_G2A, se_G2A)


def get_error_rates_dataframe(d_mismatch):
    
    d_error_rates = {}
    for key in d_mismatch.keys():
        
        p = mismatch_to_error_rate(d_mismatch[key])
        (e_all, se_all, e_C2T, se_C2T, e_G2A, se_G2A) = p
        
        d_error_rates[key] = {'all':    e_all, 
                              'all_s': se_all,
                              'C2T':    e_C2T, 
                              'C2T_s': se_C2T,
                              'G2A':    e_G2A, 
                              'G2A_s': se_G2A,
                              }
        
    df_error_rates = pd.DataFrame(d_error_rates).T
    return df_error_rates.sort_index()

df_error_rates_forward = get_error_rates_dataframe(d_mismatch_forward)
df_error_rates_reverse = get_error_rates_dataframe(d_mismatch_reverse)


# =============================================================================
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
        sy = df_error_rates.loc[:, name+'_s']
        ax.errorbar(x, y, sy, fmt='-', label=strand)
    
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





