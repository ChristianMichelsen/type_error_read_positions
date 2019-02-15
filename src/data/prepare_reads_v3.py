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


def get_corrected_filename(filename):
    return filename.replace('.txt', '_corrected.txt')


def file_length(filename):
    with open(filename) as f:
        for i, l in enumerate(f):
            pass
    return i + 1


def correctify_file(filename):
    
    file_len = file_length(filename)
    filename_output = get_corrected_filename(filename)
    
    with open(filename, 'r') as f, open(filename_output, 'w') as f_out:
        for line in tqdm(f, total=file_len):
            parts = line.split()
            for part in parts[3:]:
                if 'MD:Z' in part.upper():
                    print(*parts[:3], part.strip('MD:Z'), file=f_out)
    return file_len


def do_correct_file_if_not_exists(filename):
    
    if 'corrected' in filename:
        filename_corrected = filename
    else:
        filename_corrected = get_corrected_filename(filename)
    
    if Path(filename_corrected).is_file():
        print(f'\nDid nothing, correct file already exists, "{filename_corrected}"\n\n')
        file_len = file_length(filename)
        
    else:
        print(f'\nGenerating "{filename_corrected}"', flush=True)
        file_len = correctify_file(filename)
        print(f"\nFinished generating corrected file\n\n")
    
    return file_len

# import sys
# filename = sys.argv[1]
# do_correct_file_if_not_exists(filename)