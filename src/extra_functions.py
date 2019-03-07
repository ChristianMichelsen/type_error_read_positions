#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 15 10:51:05 2019

@author: michelsen
"""

import numpy as np
import re
import pandas as pd




#%% =============================================================================
# Read txt/bam file and process reads and references
# =============================================================================

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
    
    # TODO: check this thorughly!!!
    # start = sum([y for (x,y) in cigar_parts if x == 'S'])
    start = cigar_parts[0][1] if cigar_parts[0][0] == 'S' else 0 
    
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

import multiprocessing as mp, os
from functools import partial
from collections import Counter
from pathlib import Path
from tqdm import tqdm


ACGT_names = ['A', 'C', 'G', 'T', 'N', '0']
base2index = {val: i for i, val in enumerate(ACGT_names)}

# strand2index = {'+': 1, '-': 0}
# is_reverse2index = {True: 0, False: 1}

header = ['obs_base', 'prev_ref_base', 'position', 'strand', 'ref_base', 'L', 'counts']


def ACGTN_correct_string(string):
    """
    Corrects a string such that it only contains characters from ACGT_names
    else replace the character with an 'N'.     
    """
    return ''.join(char if char in ACGT_names else 'N' for char in string)


# def tidy_ML_seq_ref_data(ref, seq, is_reverse, res):
def tidy_ML_seq_ref_data(ref, seq, strand, res):
    
    seq = ACGTN_correct_string(seq)
    ref = ACGTN_correct_string(ref)
    
    for i, (s_ref, s_seq) in enumerate(zip(ref, seq)):
        if i == 0:
            prev_ref_base = '0'
        else:
            prev_ref_base = ref[i-1]
        
        obs_base = base2index[s_seq]
        prev_ref_base = base2index[prev_ref_base]
        position = i+1
        # strand = is_reverse2index[is_reverse]
        length = len(seq)
        ref_base = base2index[s_ref]
        # is_mismatch = (obs_base != ref_base)
        
        res[(obs_base, prev_ref_base, position, strand, ref_base, length)] += 1
    
    return None
    

def process_line_tidy_ML(line_txt, ML_res):
    strand, cigar, read, md_tag = line_txt
    seq, ref = build_alignment_reference_seq(read, cigar, md_tag)
    is_reverse, ref, seq = correct_reverse_strans(strand, seq, ref)
    tidy_ML_seq_ref_data(ref, seq, strand, ML_res)
    return None    
    

def process_chunk_ML(file_processed_in, chunk_start, chunk_size, is_last):
    ML_res = Counter()
    with open(file_processed_in) as f:
        f.seek(chunk_start)
        lines = f.read(chunk_size).splitlines()
        it = enumerate(_read_txtfile(lines))
        if is_last:
        # if True:
            it = tqdm(it, total=len(lines))
        for iline, line_txt in it:
            process_line_tidy_ML(line_txt, ML_res)
    return ML_res


def collect_result_ML(ML_res, ML_chunk):
    ML_res += ML_chunk
    return None


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
                
                is_last = True if chunkEnd > file_end else False
                
                yield chunk_start, chunk_size, is_last
                if chunkEnd > file_end:
                    break
                

def save_ML_res(ML_res, filename_mismatch_ML):
    # converts Counter dict to list of tuples that Pandas can read 
    df = []
    for key, val in dict(ML_res).items():
        df.append((*key, val))
    
    # create dataframe from list of tuples
    df = pd.DataFrame(df, columns=header)
    
    #downcast to unsigned integers (since always positive counts)
    for col in df:
        df.loc[:, col] = pd.to_numeric(df[col], downcast='unsigned')
        
    # save dataframe
    df.to_pickle(filename_mismatch_ML)
    return df

def load_ML_res(filename_mismatch_ML):
    return pd.read_pickle(filename_mismatch_ML)





def get_ML_res(file_processed_in, cores, force_rerun, N_reads):
    
    filename_mismatch_ML = file_processed_in.replace(f'corrected.txt', 
                                    f'mismatch_results_{cores}cores_ML.pkl')
    
    if not Path(filename_mismatch_ML).is_file() or force_rerun:
        
        print(f"Parsing MD-tags to get mismatch matrix using {cores} cores", flush=True)
        
        
        if cores == 1:
            
            # res = []
            ML_res = Counter()
        
            with open(file_processed_in, 'r') as f_processed:    
                # current_length = 0
                for iline, line_txt in tqdm(enumerate(_read_txtfile(f_processed)), 
                                            total=N_reads):
                    process_line_tidy_ML(line_txt, ML_res)
            
        else:
            
            ML_res = Counter()
            
            #init objects
            pool = mp.Pool(cores)
            
            #create jobs
            for chunk_start, chunk_size, is_last in chunkify(file_processed_in, 10*cores):
                # print(chunk_start, chunk_size, is_last)
                #http://blog.shenwei.me/python-multiprocessing-pool-difference-between-map-apply-map_async-apply_async/
                pool.apply_async(process_chunk_ML, 
                                 args=(file_processed_in, chunk_start, chunk_size, is_last), 
                                 callback=partial(collect_result_ML, ML_res))
            
            pool.close() # Prevents any more tasks from being submitted to the pool
            pool.join() # Wait for the worker processes to exit
        
        
        
        print('Finished parsing the MD-tags, now saving the file')
        df = save_ML_res(ML_res, filename_mismatch_ML)
    
    
    # load files
    else:
        print("Loading 1 core ML dataframe")
        df = load_ML_res(filename_mismatch_ML)

    return df



# file1 = filename_mismatch_ML
# file2 = file1.replace('6cores', '1cores') if cores==6 else file1.replace('1cores', '6cores')
def compare_multi_core_results(file1, file2, header):
    df1 = load_ML_res(file1).sort_values(by=header , ascending=False).reset_index(drop=True)
    df2 = load_ML_res(file2).sort_values(by=header, ascending=False).reset_index(drop=True)
    pd.testing.assert_frame_equal(df1, df2)
    return None



#%% =============================================================================
#  Process dataframe
# =============================================================================


def get_error_rates(df, string_list):
    
    res = []
    
    for string in string_list:
        
        series = {}
        n_len = int(df['position'].max())
        
        for i in range(n_len):
            ref = string[0]
            obs = string[2]
            i += 1
        
            mask_pos = (df['position']== i)
            mask_ref = (df['ref_base'] == base2index[ref])
            mask_obs = (df['obs_base'] == base2index[obs])
            
            num = df[mask_pos & mask_ref & mask_obs]['counts'].sum()
            den = df[mask_pos & mask_ref]['counts'].sum()
            
            if den == 0:
                series[i] = 0
            else:
                series[i] = num/den
            
        res.append(pd.Series(series, name=string))
        
    return pd.concat(res, axis=1)



def get_X_count(df, X):
    # get only one data point from each read by only taking the first position
    df_single_read = df[df['position']==1][[X, 'counts']]
    return df_single_read.groupby(X)['counts'].sum().to_dict()

def get_read_lengt_count(df):
    return get_X_count(df, 'L')

def get_strand_count(df):
    return get_X_count(df, 'strand')


#%% =============================================================================
# 
# =============================================================================



def is_linux():
    import platform
    return platform.system() == 'Linux'

def phred_symbol_to_Q_score(s):
    Q = ord(s)-33
    return Q

def Q_score_to_probability(Q):
    P = 10**-(Q/10.0)
    return P

#% =============================================================================
# 
# =============================================================================



    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    