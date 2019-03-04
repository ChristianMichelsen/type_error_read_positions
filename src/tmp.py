#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 27 11:58:27 2019

@author: michelsen
"""

import numpy as np

#pythran export dprod(int list, int list)
def dprod(l0,l1):
    """WoW, generator expression, zip and sum."""
    return sum(x * y for x, y in zip(l0, l1))



ACGT_names = ['A', 'C', 'G', 'T', 'N', '-', 'Other']
base2index = {val: i for i, val in enumerate(ACGT_names)}

# http://www.dnabaser.com/articles/IUPAC%20ambiguity%20codes.html
IGNORE_LETTERS = 'YRWSKMDVHBXN'





def init_zero_matrix():
    return np.zeros((7,7), dtype=int)
    # return np.zeros((len(ACGT_names), len(ACGT_names)), dtype=int)



# def fill_mismatch_matrix(ref, seq, d_mismatch):
#     """ seq1 = ref, seq2 = seq """
    
#     for i, (s_ref, s_seq) in enumerate(zip(ref, seq)):
#         i += 1
#         if not i in d_mismatch:
#             d_mismatch[i] = init_zero_matrix()
        
#         # try:
#         d_mismatch[i][base2index[s_ref], base2index[s_seq]] += 1
        
#         # except KeyError as e:
#         #     if verbose:
#         #         print(f'Found {e} in sequence. Ignoring it and using "Other" instead')
#         #     if s_ref in ACGT_names:
#         #         d_mismatch[i][base2index[s_ref], base2index['Other']] += 1
#         #     else:
#         #         d_mismatch[i][base2index['Other'], base2index[s_seq]] += 1
    
#     return d_mismatch




# pythran export fill_mismatch_matrix(str, str, int[:][7][7], bool)
def fill_mismatch_matrix(ref, seq, mismatch, verbose=True):
    """ seq1 = ref, seq2 = seq """
    
    for i, (s_ref, s_seq) in enumerate(zip(ref, seq)):
        try:
            mismatch[i, base2index[s_ref], base2index[s_seq]] += 1
        except KeyError: # as e:
            if s_ref in ACGT_names:
                if verbose:
                    print(s_seq + ' in sequence. Ignoring it and using "Other" instead')
                mismatch[i, base2index[s_ref], base2index['Other']] += 1
            else:
                if verbose:
                    print(s_seq + ' in reference. Ignoring it and using "Other" instead')
                mismatch[i, base2index['Other'], base2index[s_seq]] += 1
        i += 1
    
    return None


# pythran export fill_mismatch_forward_reverse(bool, str, str, int[:][7][7], int[:][7][7], bool)
def fill_mismatch_forward_reverse(is_reverse, ref, seq, 
                                  mismatch_forward, mismatch_reverse, 
                                  verbose=True):
    if not is_reverse:
        fill_mismatch_matrix(ref, seq, mismatch_forward, verbose=verbose)
    else:
        fill_mismatch_matrix(ref, seq, mismatch_reverse, verbose=verbose)

    return None
