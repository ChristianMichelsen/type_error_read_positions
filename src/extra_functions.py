#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 15 10:51:05 2019

@author: michelsen
"""

import numpy as np
import re
import pandas as pd












# =============================================================================
# 
# =============================================================================


ACGT_names = ['A', 'C', 'G', 'T', 'N', '-']
base2index = {val: i for i, val in enumerate(ACGT_names)}


# from Martin Kircher, to complement DNA
TABLE = str.maketrans('TGCAMRWSYKVHDBtgcamrwsykvhdb', \
                      'ACGTKYWSRMBDHVacgtkywsrmbdhv')

def comp(seq):
    """ return reverse complemented string """
    return seq.translate(TABLE)

def seq_to_match(read, read_pos, N_soft_clipped_beginning, d_mismatch):
    
    mask = any([s in ACGT_names for s in read])
    if not mask:
        print(read)
        assert False

    counter = read_pos + N_soft_clipped_beginning
    for base in read:
        index = base2index[base]
        dict_to_mismatch(d_mismatch, counter)[index, index] += 1
        counter += 1
        
    return None
    

def seq_to_mismatch(read, char, read_pos, N_soft_clipped_beginning, d_mismatch):
    
    assert len(char) == 1
    index_ref = base2index[char]
    index_seq = base2index[read[read_pos]]
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



def parse_md_tag(read, md_tag, cigar, strand, d_mismatch):
    
    L = len(read)
    
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
                read = read[N_soft_clipped_beginning:]
            if cigar_split[-1] == 'S':
                N_soft_clipped_end = int(cigar_split[-2])
                read = read[:-N_soft_clipped_end]
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
            seq_to_match(read[read_pos:read_pos+i], read_pos, N_soft_clipped_beginning, d_mismatch)
            read_pos += i

        # ignore inserts for now            
        elif md_split[0] == '^':
            # print(f"Ignoring deletions for now: {read}, {md_tag}, {cigar}")
            N_deletions += len(md_split[1:])

        else:
            for char in md_split:
                if char.upper() in ACGT_names:
                    seq_to_mismatch(read, char, read_pos, N_soft_clipped_beginning, d_mismatch)
                    read_pos += 1
                else: 
                    print(f"The characters {char} is not recognized")
                    assert False
    
    assert read_pos == len(read) - N_inserts
    
    return L + N_deletions

#%% =============================================================================
# 
# =============================================================================
import re
# import string

# from Martin Kircher, to complement DNA
TABLE = str.maketrans('TGCAMRWSYKVHDBtgcamrwsykvhdb', \
                         'ACGTKYWSRMBDHVacgtkywsrmbdhv')

def revcomp(seq):
    """ return reverse complemented string """
    return seq.translate(TABLE)[::-1]




def init_ref_from_cigar(cigar):
    ref = ''
    inserts = []
    
    if 'I' in cigar:
        cigar_split = list(filter(None, re.split(r'(\d+)', cigar)))
        # counter = 0
        for i, split in enumerate(cigar_split):
            if split == 'M':
                ref += int(cigar_split[i-1]) * ' '
            if split == 'I':
                ref += int(cigar_split[i-1]) * '-'
                inserts.append(len(ref)-1)
    else:
        ref += int(cigar[:-1]) * ' '
    return ref, inserts


def init_seq_from_md_tag(md_tag):
    seq = ''
    
    md_tags_splitted = list(filter(None, re.split(r'(\d+)', md_tag)))
    for md_split in md_tags_splitted: # TODO: remove 'list' when done
       
        if md_split == '0':
            continue
        
        elif md_split.isdigit():
            seq += int(md_split) * ' '
            
        # ignore inserts for now            
        elif md_split[0] == '^':
            assert False
            seq += len(md_split[1:]) * '-' 

        else:
            for char in md_split:
                if char.upper() in ACGT_names:
                    seq += char
                else: 
                    print(f"The characters {char.upper()} is not recognized")
                    assert False
                    
    return seq



def get_cumsum_md_tag(x):
    res = []
    counter = 0
    for val in x:
        if val.isdigit():
            counter += int(val)
        else:
            if val[0] == '^':
                counter += len(val) - 1
            else:
                counter += len(val)
        res.append(counter)
    return np.array(res)

def get_cumsum_cigar_and_number_of_deletions(x):
    res = []
    counter = 0
    deletions = []
    
    for i,val in enumerate(x):
        if val.upper() == 'M':
            counter += int(x[i-1])
            # res.append(counter)
        elif val.upper() == 'I':
            res.append(counter)
            counter += int(x[i-1])
            deletions.append(int(x[i-1]))
            
    return res, deletions


def cigar_both_ins_and_dels(cigar_split):
    for i, val in enumerate(cigar_split):
        if val == 'D':
            index = i
    counter = int(cigar_split[index-3]) + int(cigar_split[index-1]) + int(cigar_split[index+1])
    cigar = ''.join(cigar_split[:index-3]) + str(counter) + 'M' + ''.join(cigar_split[index+3:])
    return list(filter(None, re.split(r'(\d+)', cigar)))


def cigar_merge_operation(cigar_split, operation='I'):
    counter = 0
    for i,val in enumerate(cigar_split):
        if val.upper() == 'M':
            counter += int(cigar_split[i-1])
        elif val.upper() == operation.upper():
            counter += int(cigar_split[i-1])
            try:
                counter += int(cigar_split[i+1])
                return str(counter) + 'M' + ''.join(cigar_split[i+3:])
            except IndexError:
                return str(counter) + 'M'


def cigar_correct_for_deletion(cigar_split):
    cigar = cigar_merge_operation(cigar_split, operation='D')
    cigar_split = list(filter(None, re.split(r'(\d+)', cigar)))
    return cigar_split

def make_single_round(md_tag, cigar, insertion_symbol):
    
    md_tag_split = list(filter(None, re.split(r'(\d+)', md_tag)))
    md_tag_cumsum = get_cumsum_md_tag(md_tag_split)
    
    cigar_split = list(filter(None, re.split(r'(\d+)', cigar)))
    if 'D' in cigar.upper():
        if 'I' in cigar.upper():
            cigar_split = cigar_both_ins_and_dels(cigar_split)
        else:
            cigar_split = cigar_correct_for_deletion(cigar_split)
    cigar_cumsum, cigar_deletions = get_cumsum_cigar_and_number_of_deletions(cigar_split)
    
    index_pos = np.searchsorted(md_tag_cumsum, cigar_cumsum)
    res = md_tag_split[:]
    md_tag_cumsum = np.r_[0, md_tag_cumsum]
    

    res[index_pos[0]:index_pos[0]+1] = [str(cigar_cumsum[0]-md_tag_cumsum[index_pos[0]]), 
                                            cigar_deletions[0] * insertion_symbol, 
                                            str(md_tag_cumsum[index_pos[0]+1]-cigar_cumsum[0])]
    
    md_tag = ''.join(res).upper()
    cigar = cigar_merge_operation(cigar_split, operation='I')

    return md_tag, cigar, len(index_pos)


def make_md_tag_cigar_fixed(md_tag, cigar, insertion_symbol='-'):
    
    if 'S' in cigar:
        md_tag = correct_for_softclipping(md_tag, cigar)
    
    if 'I' not in cigar:
        return md_tag
    
    md_tag, cigar, l_pos = make_single_round(md_tag, cigar, insertion_symbol)
    while l_pos > 1:
        md_tag, cigar, l_pos = make_single_round(md_tag, cigar, insertion_symbol)
    return md_tag
    


def correct_for_softclipping(md_tag, cigar):
    cigar_split = list(filter(None, re.split(r'(\d+)', cigar)))
    
    if cigar_split[1] == 'S':
        return f"{'P'*int(cigar_split[0])}{md_tag}"
    elif cigar_split[-1]: 
        return f"{md_tag}{'S'*int(cigar_split[-2])}"
    else:
        print("Have to implement soft clipping for other positions than start or end")
        assert False


def test_make_md_tag_cigar_fixed():
    
    # test when first number of cigar is bigger than md_tag
    assert make_md_tag_cigar_fixed('4A5T2', '7M1I6M') == '4A2-3T2'
    
    # test when first number of md_tag is bigger than cigar
    assert make_md_tag_cigar_fixed('6T2', '2M1I7M') == '2-4T2'
    
    # test for 2 insertions
    assert make_md_tag_cigar_fixed('7', '2M1I2M2I3M') == '2-2--3'

    # test for more insertions
    assert make_md_tag_cigar_fixed('9', '2M1I2M1I3M2I2M') == '2-2-3--2'

    return None



    


def make_reference(read, md_tag, strand, cigar=None, is_md_cigar=False):
    """ Takes in a read, md_tag and cigar string and recreates the reference 
        and sequence """
    
    is_reverse = (int(strand)>=16)
    
    if cigar is not None and is_md_cigar:
        print("Ignoring cigar input since is_md_cigar == True")
    
    if is_md_cigar:
        md_cigar = md_tag
    else:
        md_cigar = make_md_tag_cigar_fixed(md_tag, cigar, insertion_symbol='-')
    
    
    seq = ''
    ref = '' # init_ref_from_cigar(cigar)
    counter = 0
    soft_clip_counter_end = 0
    
    
    md_tags_splitted = list(filter(None, re.split(r'(\d+)', md_cigar)))
    
    # handle soft clipping (in the beginning only)
    if 'P' in md_tags_splitted[0]:
        counter += len(md_tags_splitted[0])
    
    
    for md_split in md_tags_splitted: # TODO: remove 'list' when done
       
        # handle 0's
        if md_split == '0':
            continue
        
        # handle normal numbers
        elif md_split.isdigit():
            i = int(md_split)
            ref += read[counter:counter+i]
            seq += read[counter:counter+i]
            counter += i
            
        # Handle inserts            
        elif md_split[0] == '^':
            ref += md_split[1:]
            i = len(md_split[1:])
            seq += i*'-' 
        
        # handle soft clipping (in the end only)
        elif 'S' in md_split:
            soft_clip_counter_end += len(md_split)
        
        # handle soft clipping (in the beginning only)
        elif 'P' in md_split:
            continue
        
        # handle everything else
        else:
            for char in md_split:
                if char.upper() in ACGT_names:
                    ref += char.upper()
                    seq += read[counter]
                else: 
                    print(f"The characters '{char.upper()}' is not recognized")
                    assert False
                counter += 1
        
    
    seq += read[counter:-soft_clip_counter_end]
    ref += read[counter:-soft_clip_counter_end]

    assert len(ref) == len(seq)
    
    if is_reverse:
        seq = comp(seq)
        ref = comp(ref)
    

    return ref, seq


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


# =============================================================================
# 
# =============================================================================
