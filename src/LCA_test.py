#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 16:06:04 2019

@author: michelsen
"""

import numpy as np
import pandas as pd


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



qname = 'D00708:299:CB6VWANXX:5:2307:8147:26950'


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

it = _read_bamfile(in_bam)
            
# main loop
for counter, read_bam in enumerate(it, start=1):
    # counter += 1
    
    if read_bam.qname == qname:
        assert False
        
        
        


























