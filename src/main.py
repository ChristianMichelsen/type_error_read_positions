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
from src.extra_functions import (get_error_rates_dataframe, make_reference,
                                 fill_mismatch_matrix, is_linux)

save_plots = True
do_ancient = False

if do_ancient:
    print("\nRunning on ancient DNA")
    filename = 'ESW_LRUE_MA2621_error_test.txt' # ancient data
    filename_bam = 'ESW_LRUE_MA2621_L1_S1_S83_L006.sort.rmdup.bam'
else:
    print("\nRunning on modern DNA")
    filename = 'NA12400_error_test.txt' # modern dna
    filename = 'NA12400_small_bam_1_000_000_error_test.txt' # modern dna
    filename_bam = 'NA12400_small_bam_1_000_000.bam'
refname = 'hs37d5.fa'


plot_prefix = f"../figures/{filename.split('_')[0]}_plot_"


def append_corrected_to_filename(filename):
    return filename.replace('.txt', '_corrected.txt')

def get_filename_input(filename):
    return f'../data/raw/{filename}'
def get_filename_output(filename):
    return f'../data/processed/{append_corrected_to_filename(filename)}'
def get_filename_processed(filebame):
    return get_filename_output(filename)


file_processed_in = get_filename_processed(filename)
file_bam_in = get_filename_input(filename_bam)
ref_in = get_filename_input(refname)


_, N_reads = get_filename_and_lenght(filename)



#%% =============================================================================
# 
# =============================================================================

import pysam
from datetime import datetime

from mapdamage_custom import (read_fasta_index, initialize_mut, 
                              initialize_comp, initialize_lg, 
                              _read_bamfile, get_coordinates, record_lg,
                              get_around, options, align, revcomp,
                              record_soft_clipping, compare_sequence_dicts,  
                              get_mis, count_read_comp, count_ref_comp,
                              compare_str, compare_str_print)


from extra_functions import (build_alignment_sequence, build_reference_sequence, 
                             align_ref, cigar_string_to_tuples, comp, _read_txtfile)


x=x

def main():

    start_time = datetime.now()

    # fetch all references and associated lengths in nucleotides
    ref_bam = pysam.FastaFile(ref_in) 
    
    # open SAM/BAM file
    in_bam = pysam.AlignmentFile(file_bam_in)
    
    N_reads = int(pysam.view('-c', '-F 4', f'{file_bam_in}')) # approximate
    
    reflengths = dict(zip(in_bam.references, in_bam.lengths))
    # check if references in SAM/BAM are the same in the fasta reference file
    fai_lengths = read_fasta_index(ref_in + ".fai")
    
    
    # #TODO
    # if not fai_lengths:
    #     return 1
    # elif not compare_sequence_dicts(fai_lengths, reflengths):
    #     return 1
    
    refnames = in_bam.references
    
    
    
    # for misincorporation patterns, record mismatches
    misincorp = initialize_mut(refnames, options.length)
    # for fragmentation patterns, record base compositions
    dnacomp =  initialize_comp(refnames, options.around, options.length)
    # for length distributions
    lgdistrib =  initialize_lg()
    
    
    counter = 0
    
    do_tqdm = True
    
    
    with open(file_processed_in, 'r') as f_processed:
        
        it = tqdm(zip(_read_bamfile(in_bam), _read_txtfile(f_processed)), total=N_reads) if do_tqdm else zip(_read_bamfile(in_bam), _read_txtfile(f_processed))
            
        # main loop
        for read_bam, line_txt in it:
            counter += 1
            
            # if read_bam.is_reverse:
            #     continue
            
            strand_processed, cigar_processed, read_processed, md_tag_processed = line_txt
            
            # md_tag, cigar = md_tag_processed, cigar_processed
            # ref_processed, seq_processed = make_reference(read_processed, md_tag_processed, strand_processed, cigar_processed, is_md_cigar=False)
            # ref_processed, seq_processed = ref, seq
            
            cigarlist = cigar_string_to_tuples(cigar_processed)
            
            s_align_seq = build_alignment_sequence(read_processed, cigar_processed, md_tag_processed)
            s_ref_seq_no_gaps = build_reference_sequence(read_processed, cigar_processed, md_tag_processed)
            s_ref_seq = align_ref(cigarlist, s_ref_seq_no_gaps)
            
            
            ref_processed, seq_processed = s_ref_seq.upper(), s_align_seq.upper()
            
            
            is_reverse = (strand_processed & 0x10) != 0
            
            assert read_bam.is_reverse == is_reverse
            
            if is_reverse:
                ref_processed = comp(ref_processed)
                seq_processed = comp(seq_processed)
                
            
            
            # print(compare_str_print(seq_bam, s_align_seq.upper())) # GOOD
            # print(compare_str_print(refseq_bam, s_ref_seq.upper())) # good
            
            
            
            
    
            # external coordinates 5' and 3' , 3' is 1-based offset
            coordinate = get_coordinates(read_bam)
            # record aligned length for single-end read_bams
            lgdistrib = record_lg(read_bam, coordinate, lgdistrib)
            # fetch reference name, chromosome or contig names
            chrom = in_bam.get_reference_name(read_bam.reference_id)
    
            
    
            (before, after) = get_around(coordinate, chrom, reflengths, options.around, ref_bam)
            refseq_bam_no_gaps = ref_bam.fetch(chrom, min(coordinate), max(coordinate)).upper()
            # read_bam.query_alignment_sequence contains aligned sequences while read_bam.seq is the read_bam itself
            seq_bam_no_gaps = read_bam.query_alignment_sequence
    
            # add gaps according to the cigar string
            (seq_bam, refseq_bam) = align(read_bam.cigar, seq_bam_no_gaps, refseq_bam_no_gaps)
    
        
            # reverse complement read_bam and reference when mapped reverse strand
            if read_bam.is_reverse:
                refseq_bam = revcomp(refseq_bam)
                seq_bam = revcomp(seq_bam)
                beforerev = revcomp(after)
                after = revcomp(before)
                before = beforerev
                
                seq_bam = seq_bam[::-1]
                refseq_bam = refseq_bam[::-1]
            
        
        
            if (compare_str(seq_bam, seq_processed) != 1).sum() != 0:
                print("\nSeq")
                print(counter)
                print(seq_bam)
                print(seq_processed)
                print(compare_str_print(seq_bam, seq_processed))
                assert False
                
            if (compare_str(refseq_bam, ref_processed, ignore_letters=True) != 1).sum() != 0:
                print("\nRef")
                print(counter)
                print(refseq_bam)
                print(ref_processed)
                print(compare_str_print(refseq_bam, ref_processed, ignore_letters=True))
                assert False
            
                
            
            
            # record soft clipping when present
            record_soft_clipping(read_bam, misincorp[chrom], options.length)
    
    
            # count misincorparations by comparing read_bam and reference base by base
            get_mis(read_bam, seq_bam, refseq_bam, chrom, options.length, misincorp, '5p')
            # do the same with sequences align to 3'-ends
            get_mis(read_bam, seq_bam, refseq_bam, chrom, options.length, misincorp, '3p')
            # compute base composition for read_bams
            count_read_comp(read_bam, chrom, options.length, dnacomp)
    
            # compute base composition for genomic regions
            count_ref_comp(read_bam, chrom, before, after, dnacomp)
    
    
    # close file handlers
    in_bam.close()
    ref_bam.close()

    print(f'Run completed in {datetime.now() - start_time:f} seconds')

    return 0



#%% =============================================================================
# #
# =============================================================================


import pickle
from pathlib import Path

filename_mismatch = file_processed_in.replace('corrected.txt', 'mismatch_results.pkl')


if not Path(filename_mismatch).is_file():
    
    d_mismatch_forward = {}
    d_mismatch_reverse = {}
    
    lengths_forward = []
    lengths_reverse = []
    
    list_strand = np.zeros(N_reads)
    
    print("Parsing MD-tags to get mismatch matrix: ", flush=True)
    
    with open(file_processed_in, 'r') as f_processed:
        for iline, line_txt in tqdm(enumerate(_read_txtfile(f_processed)), total=N_reads):
            
            strand_processed, cigar_processed, read_processed, md_tag_processed = line_txt
            
            cigarlist = cigar_string_to_tuples(cigar_processed)
            
            s_align_seq = build_alignment_sequence(read_processed, cigar_processed, md_tag_processed)
            s_ref_seq_no_gaps = build_reference_sequence(read_processed, cigar_processed, md_tag_processed)
            s_ref_seq = align_ref(cigarlist, s_ref_seq_no_gaps)
            
            ref_processed, seq_processed = s_ref_seq.upper(), s_align_seq.upper()
            
            is_reverse = (strand_processed & 0x10) != 0
            if is_reverse:
                ref_processed = comp(ref_processed)
                seq_processed = comp(seq_processed)
                



            
            L = len(seq_processed)
            
            is_reverse = (strand_processed>=16)
            list_strand[iline] = strand_processed
            
            if not is_reverse:
                fill_mismatch_matrix(ref_processed, seq_processed, d_mismatch_forward)
                lengths_forward.append(L)
            else:
                fill_mismatch_matrix(ref_processed, seq_processed, d_mismatch_reverse)
                lengths_reverse.append(L)
    
            
    # list_strand = list_strand
    lengths_forward = np.array(lengths_forward)        
    lengths_reverse = np.array(lengths_reverse)

    print('Finished parsing the MD-tags, now saving the file')

    mismatch_results = [d_mismatch_forward, d_mismatch_reverse, 
                        list_strand, lengths_forward, lengths_reverse]

    with open(filename_mismatch, 'wb') as file:  
        # store the data as binary data stream
        pickle.dump(mismatch_results, file)

else:

    with open(filename_mismatch, 'rb') as file:  
        # read the data as binary data stream
        p = pickle.load(file)
        d_mismatch_forward, d_mismatch_reverse, list_strand, lengths_forward, lengths_reverse = p


df_error_rates_forward = get_error_rates_dataframe(d_mismatch_forward)
df_error_rates_reverse = get_error_rates_dataframe(d_mismatch_reverse)


#%% =============================================================================
# Error rate plots
# =============================================================================

if not is_linux():
    
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
        
        
        ax.set(xlabel='Read Pos', ylabel='Error Rate', title=name,
               xlim=(1, 25), ylim=(0, 0.3))
        
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
    
