#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 18 10:37:51 2018

@author: michelsen
"""

from pathlib import Path
from tqdm import tqdm


# filename = 'NA12400_error_test.txt'
# filename = 'ESW_LRUE_MA2621_error_test.txt'


def append_corrected_to_filename(filename):
    return filename.replace('.txt', '_corrected.txt')


def get_file_length(filename):
    with open(filename) as f:
        for i, l in enumerate(f):
            pass
    return i + 1


def preprocess_line(line):
    parts = line.split()
    string_out = ''
    for part in parts[3:]:
        if 'MD:Z' in part.upper():
            # print(*parts[:3], part.strip('MD:Z'), file=f_out)
            string_out = parts[:3] + [part.strip('MD:Z')]
    return '\t'.join(string_out)
    

def preprocess_file(file_in, file_out, use_tqdm=True):
    file_len = get_file_length(file_in)
    
    with open(file_in, 'r') as f, open(file_out, 'w') as f_out:
        it = tqdm(f, total=file_len) if use_tqdm else f
        for line in it:
            string_out = preprocess_line(line)
            print(string_out, file=f_out)
        
        # for line in tqdm(f, total=file_len):
        #     parts = line.split()
        #     for part in parts[3:]:
        #         if 'MD:Z' in part.upper():
        #             print(*parts[:3], part.strip('MD:Z'), file=f_out)
    return file_len


def preprocess_file_parallel(file_in):
    file_out = file_in + '_corrected'
    preprocess_file(file_in, file_out, use_tqdm=False)
    return None



def do_correct_file_if_not_exists(filename, run_parallel=True):
    
    if 'corrected' in filename:
        filename_corrected = filename
    else:
        filename_corrected = append_corrected_to_filename(filename)
    
    filename_corrected = f'../data/processed/{filename_corrected}'
    
    if Path(filename_corrected).is_file():
        print(f'\nDid nothing, correct file already exists, \n"{filename_corrected}"\n\n')
        file_len = get_file_length(filename_corrected)
        
    else:
        print(f'\nPreprocessing "{filename_corrected}"', flush=True)
        file_in = f'../data/raw/{filename}'
        file_out = f'../data/processed/{append_corrected_to_filename(filename)}'
        
        if run_parallel:
            n_cores = multiprocessing.cpu_count()-1 or 1
            file_len = parallel_preprocessing(file_in, file_out, n_cores)  
        else:
            file_len = preprocess_file(file_in, file_out)
        
        print(f"\nFinished Preprocessing corrected file\n\n")
    
    return file_len

# import sys
# filename = sys.argv[1]
# do_correct_file_if_not_exists(filename)
    

def get_filename_and_lenght(filename, run_parallel=True): 
    
    file_len = do_correct_file_if_not_exists(filename, run_parallel)
    if not "corrected" in filename:
        filename = append_corrected_to_filename(filename)
    filename = f'../data/processed/{filename}'
    
    return filename, file_len


#%% =============================================================================
# 
# =============================================================================

filename = 'test.txt'
n_splits = 3


import multiprocessing


def parallel_preprocessing(file_in, file_out, n_splits, out_suffix='_tmp_'):
    
    file_len = get_file_length(file_in)
    chunk_size = file_len//n_splits
    
    import os
    import glob
    # os.system(f'split -l --numeric-suffixes {chunk_size} {filename} {out_suffix}')
    os.system(f'split -l {chunk_size} {file_in} {out_suffix}')
    
    files = glob.glob(f"{out_suffix}*")
    
    pool = multiprocessing.Pool(n_splits)
    pool.map(preprocess_file_parallel, files)
    pool.close()
    
    os.system(f"cat *_corrected > {file_out}")
    os.system(f"rm {out_suffix}*")

    return file_len 


# =============================================================================
# 
# =============================================================================


# from datetime import datetime
# import platform
# if platform.system() == 'Linux':

#     filename = 'NA12400_error_test.txt'
#     filename = 'ESW_LRUE_MA2621_error_test.txt'
    
#     startTime = datetime.now()
#     preprocess_file(filename, append_corrected_to_filename(filename))
#     print(datetime.now() - startTime)
    
#     startTime = datetime.now()
#     n_cores = multiprocessing.cpu_count()-1 or 1
#     parallel_preprocessing(filename, append_corrected_to_filename(filename), n_cores)  
#     print(datetime.now() - startTime)
    
    
