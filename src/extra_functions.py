#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 15 10:51:05 2019

@author: michelsen
"""

import numpy as np
import re
import pandas as pd


ACGT_names = ['A', 'C', 'G', 'T', 'N']
base2index = {val: i for i, val in enumerate(ACGT_names)}


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
cigar_all_options = list('MIDNSHP=X')
cigar_not_implemented = [c for c in cigar_all_options if not c in cigar_implemented]


def find_all(a_str, sub):
    start = 0
    while True:
        start = a_str.find(sub, start)
        if start == -1: 
            return
        yield start
        start += len(sub) # use start += 1 to find overlapping matches



def gen_zero_matrix(ACGT_names):
    return np.zeros((len(ACGT_names), len(ACGT_names)), dtype=int)


def dict_to_mismatch(d_mismatch, read_pos):
    read_pos_1_indexed = read_pos + 1
    if not read_pos_1_indexed in d_mismatch:
        d_mismatch[read_pos_1_indexed] = gen_zero_matrix(ACGT_names)
    return d_mismatch[read_pos_1_indexed]



def parse_md_tag(seq, md_tag, cigar, strand, d_mismatch):
    
    L = len(seq)
    
    read_pos = 0
    N_inserts = 0
    N_soft_clipped_beginning = 0
    N_soft_clipped_end = 0
    N_deletions = 0
    
    
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
            # print(f"Ignoring deletions for now: {seq}, {md_tag}, {cigar}")
            N_deletions += len(md_split[1:])

        else:
            for char in md_split:
                if char.upper() in ACGT_names:
                    seq_to_mismatch(seq, char, read_pos, N_soft_clipped_beginning, d_mismatch)
                    read_pos += 1
                else: 
                    print(f"The characters {char} is not recognized")
                    assert False
    
    assert read_pos == len(seq) - N_inserts
    
    return L + N_deletions





#%% =============================================================================
# 
# =============================================================================



def mismatch_to_error_rate(mismatch):
    diag_sum = np.diag(mismatch[:4, :4]).sum()
    sum_all = mismatch[:4, :4].sum()
    
    e_all = 1 - diag_sum / sum_all
    # se_all = np.sqrt(e_all*(1-e_all)/sum_all)
    
    e_C2T = mismatch[base2index['C'], base2index['T']] / mismatch[base2index['C'], :].sum()
    # se_C2T = np.sqrt(e_C2T*(1-e_C2T)/sum_all)
    
    e_G2A = mismatch[base2index['G'], base2index['A']] / mismatch[base2index['G'], :].sum()
    # se_G2A = np.sqrt(e_G2A*(1-e_G2A)/sum_all)
    
    # return (e_all, se_all, e_C2T, se_C2T, e_G2A, se_G2A)
    return (e_all, e_C2T, e_G2A)


def get_error_rates_dataframe(d_mismatch):
    
    d_error_rates = {}
    for key in d_mismatch.keys():
        
        p = mismatch_to_error_rate(d_mismatch[key])
        # (e_all, se_all, e_C2T, se_C2T, e_G2A, se_G2A) = p
        (e_all, e_C2T, e_G2A) = p
        
        d_error_rates[key] = {'all':    e_all, 
                              # 'all_s': se_all,
                              'C2T':    e_C2T, 
                              # 'C2T_s': se_C2T,
                              'G2A':    e_G2A, 
                              # 'G2A_s': se_G2A,
                              }
        
    df_error_rates = pd.DataFrame(d_error_rates).T
    return df_error_rates.sort_index()


