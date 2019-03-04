#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 15 10:51:05 2019

@author: michelsen
"""

import numpy as np
import re
import pandas as pd


_BAM_UNMAPPED  = 0x4
_BAM_SECONDARY = 0x100
_BAM_FAILED_QC = 0x200
_BAM_PCR_DUPE  = 0x400
_BAM_CHIMERIC  = 0x800

filtered_flags = _BAM_UNMAPPED | \
                 _BAM_SECONDARY | \
                 _BAM_FAILED_QC | \
                 _BAM_PCR_DUPE | \
                 _BAM_CHIMERIC
                 
BAM_CMATCH      = 0
BAM_CINS        = 1
BAM_CDEL        = 2
BAM_CREF_SKIP   = 3
BAM_CSOFT_CLIP  = 4
BAM_CHARD_CLIP  = 5
BAM_CPAD        = 6
BAM_CEQUAL      = 7
BAM_CDIFF       = 8
BAM_CBACK       = 9


# from Martin Kircher, to complement DNA
TABLE = str.maketrans('TGCAMRWSYKVHDBtgcamrwsykvhdb', \
                      'ACGTKYWSRMBDHVacgtkywsrmbdhv')



def _read_txtfile(txtfile):
    """
    Reads the txt file, but skip failed reads (unmapped, secondary, 
    failed qc, duplicates..) based on their strand information. 
    Iterator.
    """
    for line in txtfile:
        strand, cigar, read, md_tag = line.split()
        strand = int(strand)
        
        if not (strand & filtered_flags):
            yield (strand, cigar, read, md_tag)

def comp(seq):
    """ return reverse complemented string """
    return seq.translate(TABLE)


CIGAR_REGEX = re.compile("(\d+)([MIDNSHP=XB])")
CODE2CIGAR= "MIDNSHP=XB"

CIGAR2CODE = dict([y, x] for x, y in enumerate(CODE2CIGAR))
maketrans = TABLE

def get_cigar_parts(cigar):
    parts = CIGAR_REGEX.findall(cigar)
    cigar_parts = [(y, int(x)) for (x,y) in parts]
    return cigar_parts
    
def parse_cigar_parts(cigar_parts, ope):
    """ for a specific operation (mismach, match, insertion, deletion... see above)
    return occurences and index in the alignment """
    tlength = 0
    coordinate = []
    # count matches, indels and mismatches
    # oplist = (0, 1, 2, 7, 8)
    oplist = ['M', 'I', 'D', '=', 'X']
    for operation, length in cigar_parts:
        if operation == ope:
                coordinate.append([length, tlength])
        if operation in oplist: 
                tlength += length
    return coordinate


def align_ref(cigar_parts, ref):
    """ insert gaps according to the cigar string 
    insertions: gaps to be inserted into reference sequence """    
    lref = list(ref)
    # for nbr, idx in parse_cigar(cigarlist, 1):
    for nbr, idx in parse_cigar_parts(cigar_parts, 'I'):
        lref[idx:idx] = ["-"] * nbr
    return "".join(lref)


def get_alignment_length(cigar_parts):
    l = 0
    cigar_parts 
    for op, length in cigar_parts:
        if op == 'S' or op == 'H':
            continue
        l += length
    return l


def get_md_reference_length(md_tag):
    md_split = list(filter(None, re.split(r'(\d+)', md_tag)))
    counter = 0
    for md in md_split:
        if md.isdigit():
            counter += int(md)
        elif md[0] == '^':
            counter += len(md[1:])
        else:
            counter += len(md)
    return counter


def build_alignment_sequence(seq, cigar, md_tag):
    """return expanded sequence from CIGAR and checks length with MD-tag.
    """

    cigar_parts = get_cigar_parts(cigar)
    # cigartuples = cigar_string_to_tuples(cigar)
    
    start = sum([y for (x,y) in cigar_parts if x == 'S'])
    
    # get read sequence, taking into account soft-clipping
    read_sequence = seq[start:]
    
    max_len = get_alignment_length(cigar_parts)
    if max_len == 0:
        raise ValueError("could not determine alignment length")

    r_idx = 0
    s = ''
    for op, l in cigar_parts:
        if op == 'M' or op == '=' or op == 'X':
            s += read_sequence[r_idx:r_idx+l]
            r_idx += l
        elif op == 'D':
            s += l*'-'
        elif op == 'N':
            pass
        # encode insertions into reference as lowercase
        elif op == 'I':
            s += read_sequence[r_idx:r_idx+l].lower()
            r_idx += l
        elif op == 'S':
            pass
        # advances neither
        elif op == 'H':
            pass 
        elif op == 'P':
            raise NotImplementedError(
                "Padding (BAM_CPAD, 6) is currently not supported. "
                "Please implement. Sorry about that.")

    # Check if MD tag is valid by matching CIGAR length to MD tag defined length
    # Insertions would be in addition to what is described by MD, so we calculate
    # the number of insertions seperately.
    insertions = sum([y for (x,y) in cigar_parts if x=='I'])

    md_len = get_md_reference_length(md_tag)
    
    if md_len + insertions > max_len:
        raise AssertionError(f"Invalid MD tag: MD length {md_len} mismatch with CIGAR length {max_len} and {insertions} insertions")
    
    return s[:max_len]


def insert_substring_in_string_at_pos(s, s_insert, index):
    return s[:index] + s_insert + s[index+1:]


def md_correct_reference_sequence(ref_seq_md_corrected, md_tag):
    
    ref_seq = ref_seq_md_corrected[:]
    
    md_tag_lst = list(filter(None, re.split(r'(\d+)', md_tag)))
    
    nmatches = 0
    s_idx = 0
    
    for md in md_tag_lst:
        # c is numerical, 0 to 9
        if md.isdigit():
            nmatches += int(md)
            continue
        else:
            # save matches up to this point, skipping insertions
            for x in range(nmatches-1, -1, -1):
                while ref_seq[s_idx] >= 'a':
                    s_idx += 1
                s_idx += 1
            while ref_seq[s_idx] >= 'a':
                s_idx += 1

            nmatches = 0
            if md[0] == '^':
                for s_insert in md[1:]:
                    ref_seq = insert_substring_in_string_at_pos(ref_seq, s_insert, s_idx)
                    s_idx += 1
                
            else:
                # save mismatch, enforce lower case
                s_insert = md.lower()
                ref_seq = insert_substring_in_string_at_pos(ref_seq, s_insert, s_idx)
                s_idx += 1
    
    return ref_seq


def cigar_correct_reference_sequence(ref_seq, cigar_parts):
    
    s = ''
    
    cref_seq = ref_seq[:]
    r_idx = 0
    for op, l in cigar_parts:
        if op == 'M' or op == '=' or op == 'X':
            s += cref_seq[r_idx:r_idx+l]
            r_idx += l
        elif op == 'D':
            s += cref_seq[r_idx:r_idx+l]
            r_idx += l
        elif op == 'N':
            pass
        elif op == 'I':
            r_idx += l
        elif op == 'S':
            pass
        elif op == 'H':
            pass # advances neither
        elif op == 'P':
            raise NotImplementedError(
                "Padding (BAM_CPAD, 6) is currently not supported. "
                "Please implement. Sorry about that.")
            
            
    """ insert gaps according to the cigar string 
    insertions: gaps to be inserted into reference sequence """    
    lref = list(s)
    for nbr, idx in parse_cigar_parts(cigar_parts, 'I'):
        lref[idx:idx] = ["-"] * nbr
    return "".join(lref)


def build_reference_sequence(alignment_seq, cigar, md_tag):
    """return the reference sequence in the region that is covered by the
    alignment of the read to the reference.
    This method requires the MD tag to be set.
    """
    
    cigar_parts = get_cigar_parts(cigar)
    # alignment_seq = build_alignment_sequence(seq, cigar, md_tag)
    
    ref_seq_md_corrected = md_correct_reference_sequence(alignment_seq, md_tag)
    ref_seq = cigar_correct_reference_sequence(ref_seq_md_corrected, cigar_parts)
            
    return ref_seq.upper()


def build_alignment_reference_seq(read, cigar, md_tag):
    
    seq = build_alignment_sequence(read, cigar, md_tag)
    ref = build_reference_sequence(seq, cigar, md_tag)
    
    return seq.upper(), ref.upper()


def correct_reverse_strans(strand, seq, ref):
    is_reverse = ((strand & 0x10) != 0)
    if is_reverse:
        ref = comp(ref)
        seq = comp(seq)
    return is_reverse, ref, seq




#%% =============================================================================
# 
# =============================================================================


def is_linux():
    import platform
    return platform.system() == 'Linux'


ACGT_names = ['A', 'C', 'G', 'T', 'N', '-', 'Other']
base2index = {val: i for i, val in enumerate(ACGT_names)}

# http://www.dnabaser.com/articles/IUPAC%20ambiguity%20codes.html
IGNORE_LETTERS = 'YRWSKMDVHBXN'


MAX_LENGTH = 1000

def init_zero_matrix(ACGT_names):
    return np.zeros((MAX_LENGTH, 7,7), dtype=int)


def fill_mismatch_matrix(ref, seq, mismatch, verbose=True):
    """ seq1 = ref, seq2 = seq """
    
    for i, (s_ref, s_seq) in enumerate(zip(ref, seq)):
        
        try:
            mismatch[i, base2index[s_ref], base2index[s_seq]] += 1
        
        except KeyError as e:
            if s_ref in ACGT_names:
                if verbose:
                    print(f'Found {e} in sequence. Ignoring it and using "Other" instead')
                mismatch[i, base2index[s_ref], base2index['Other']] += 1
            else:
                if verbose:
                    print(f'Found {e} in reference. Ignoring it and using "Other" instead')
                mismatch[i, base2index['Other'], base2index[s_seq]] += 1
                
        i += 1
    
    return mismatch


# def fill_mismatch_forward_reverse(is_reverse, ref, seq, 
#                                   mismatch_forward, mismatch_reverse, 
#                                   verbose=True):
#     if not is_reverse:
#         fill_mismatch_matrix(ref, seq, mismatch_forward, verbose=verbose)
#     else:
#         fill_mismatch_matrix(ref, seq, mismatch_reverse, verbose=verbose)



# def dict_count(d, key):
#     if key in d:
#         d[key] += 1
#     else:
#         d[key] = 1
#     return None


def fill_results(is_reverse, ref, seq, strand, d_res, verbose=True):
    
    L = len(seq)
    d_res['strand'][strand] += 1
    
    # forward
    if not is_reverse:
        fill_mismatch_matrix(ref, seq, d_res['mismatch']['+'], verbose=verbose)
        d_res['lengths']['+'][L] += 1
    
    # reverse
    else:
        fill_mismatch_matrix(ref, seq, d_res['mismatch']['-'], verbose=verbose)
        d_res['lengths']['-'][L] += 1

    return None


# =============================================================================
# 
# =============================================================================


def calc_ACGT_content(mismatch1, mismatch2=None):
    mismatch = mismatch1[:, :4, :4]
    if isinstance(mismatch2, np.ndarray):
        mismatch += mismatch2[:, :4, :4]
    total_sum = mismatch.flatten().sum()
    return mismatch.sum(0).sum(1)/total_sum


def mismatch_to_dataframe_collapsed(mismatch):
    
    ref_names = ['Ref '+s for s in ACGT_names]
    read_names = ['Read '+s for s in ACGT_names]
    
    df = pd.DataFrame(mismatch.sum(0), index=ref_names, columns=read_names)
    return df


def mismatch_to_error_rate(mismatch):
    diag_sum = np.diag(mismatch[:4, :4]).sum()
    sum_all = mismatch[:4, :4].sum()
    
    e_all = 1 - diag_sum / sum_all
    e_C2T = mismatch[base2index['C'], base2index['T']] / mismatch[base2index['C'], :].sum()
    e_G2A = mismatch[base2index['G'], base2index['A']] / mismatch[base2index['G'], :].sum()
    return (e_all, e_C2T, e_G2A)


def get_error_rates_dataframe(mismatch, max_len=None):
    
    if max_len is None:
        max_len = mismatch.shape[0]
    
    d_error_rates = {}
    for i in range(max_len):
        
        p = mismatch_to_error_rate(mismatch[i, :, :])
        (e_all, e_C2T, e_G2A) = p
        
        d_error_rates[i] = {
                              # 'all':    e_all, 
                              'C2T':    e_C2T, 
                              'G2A':    e_G2A, 
                              }
        
    df_error_rates = pd.DataFrame(d_error_rates).T
    return df_error_rates.sort_index()




def phred_symbol_to_Q_score(s):
    Q = ord(s)-33
    return Q

def Q_score_to_probability(Q):
    P = 10**-(Q/10.0)
    return P