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
from src.extra_functions import parse_md_tag, get_error_rates_dataframe, make_reference

save_plots = False

filename = 'NA12400_error_test.txt' # modern dana
filename = 'ESW_LRUE_MA2621_error_test.txt' # ancient data
filename_bam = 'ESW_LRUE_MA2621_L1_S1_S83_L006.sort.rmdup.bam'
refname = 'hs37d5.fa'


ACGT_names = ['A', 'C', 'G', 'T', 'N', '-']
base2index = {val: i for i, val in enumerate(ACGT_names)}

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

get_filename_and_lenght(filename)



#%% =============================================================================
# 
# =============================================================================

import pysam
from datetime import datetime
import sys
import collections
# import itertools


# from Martin Kircher, to complement DNA
TABLE = str.maketrans('TGCAMRWSYKVHDBtgcamrwsykvhdb', \
                      'ACGTKYWSRMBDHVacgtkywsrmbdhv')

LETTERS = ("A", "C", "G", "T")
MUTATIONS = ('G>A', 'C>T', 'A>G', 'T>C', 'A>C', 'A>T', 'C>G', 'C>A', 'T>G',
             'T>A', 'G>C', 'G>T', 'A>-', 'T>-', 'C>-', 'G>-', '->A', '->T', 
             '->C', '->G', 'S')
HEADER = LETTERS + ("Total", ) + MUTATIONS


_BAM_UNMAPPED  = 0x4
_BAM_SECONDARY = 0x100
_BAM_FAILED_QC = 0x200
_BAM_PCR_DUPE  = 0x400
_BAM_CHIMERIC  = 0x800


def _filter_reads(bamfile):
    filtered_flags = _BAM_UNMAPPED | \
                     _BAM_SECONDARY | \
                     _BAM_FAILED_QC | \
                     _BAM_PCR_DUPE | \
                     _BAM_CHIMERIC

    for read in bamfile:
        if not (read.flag & filtered_flags):
            yield read

def _read_bamfile(bamfile):
    """
    Takes a subset of the bamfile. Can use a approximate fraction of the
    hits or specific number of reads using reservoir sampling. Returns
    a list in the last case otherwise a iterator.
    """
    
    return _filter_reads(bamfile)


def read_fasta_index(filename):
    """ from a fasta index file, fai, return dictionary of references:lengths """

    def print_err(msg, filename, line):
        sys.stderr.write(f"Error: {msg:s}\n")
        sys.stderr.write(f"       Filename: {filename:s}\n")
        sys.stderr.write(f"       Line:     {repr(line):s}\n")
  
    fai = {}
    with open(filename, 'r') as handle:
        for line in handle:
            ref = line.split("\t")
            if len(ref) != 5:
                print_err(f"Line in fasta index contains wrong number of fields, found {len(ref):i}, expected 5:", filename, line)
            # return None
    
            try:
                fai[ref[0]] = int(ref[1])
            except ValueError:
                print_err("Column 2 in FASTA index did not contain a number, found '{ref[1]:s}':", filename, line)
                # return None
    
    if not fai:
        sys.stderr.write("Error: Index for {filename} does contain any sequences.\n")
        sys.stderr.write("       Please ensure that FASTA file is valid, and\n")
        sys.stderr.write("       re-index file using 'samtool faidx'.\n")
        return None
    
    return fai


def compare_sequence_dicts(fasta_dict, bam_dict):
    """Compares a FASTA and BAM sequence dictionary, and prints any differences.
    Returns true if all required sequences are found, false otherwise."""
    
    if fasta_dict == bam_dict:
        return True

    sys.stderr.write("Sequence dictionaries in FASTA/BAM files differ:\n")
    common = set(fasta_dict) & set(bam_dict)
    if not common:
        sys.stderr.write("FATAL ERROR: No sequences in common!\n")
        return False

    # Check that the lengths of required sequences match (fatal error)
    different = []
    for key in sorted(common):
        if fasta_dict[key] != bam_dict[key]:
            different.append((key, fasta_dict[key], bam_dict[key]))

    if different:
        sys.stderr.write("FATAL ERROR: Length of required sequences differ:\n")
        for values in different:
            sys.stderr.write("    - %s: %i bp vs %i bp\n" % values)

    # Check for sequences only found in the BAM file (fatal errors)
    bam_only = set(bam_dict) - common
    if bam_only:
        sys.stderr.write("FATAL ERROR: Sequences missing from FASTA dictionary:\n")
        for key in bam_only:
            sys.stderr.write("    - %s = %i bp\n" % (key, bam_dict[key]))

    # Check for sequences only found in the BAM file (fatal errors)
    fasta_only = set(fasta_dict) - common
    if fasta_only:
        sys.stderr.write("WARNING: FASTA file contains extra sequences:\n")
        for key in fasta_only:
            sys.stderr.write("    - %s = %i bp\n" % (key, fasta_dict[key]))
        
    sys.stderr.write("\n")

    return not (different or bam_only)


def initialize_mut(ref, length):    
    tab = {}
    for contig in ref:
        tab_contig = tab[contig] = {}
        for end in ('5p','3p'):
            tab_end = tab_contig[end] = {}
            for std in ('+','-'):
                tab_std = tab_end[std] = {}
                for mut in HEADER:
                    tab_std[mut] = dict.fromkeys(range(length), 0)
    return tab


def initialize_comp(ref, around, length):
    keys = {"3p" : list(range(-length, 0)) + list(range(1, around + 1)),
            "5p" : list(range(-around, 0)) + list(range(1, length + 1))}
    tab = {}
    for contig in ref:
        tab_contig = tab[contig] = {}
        for end in ('5p','3p'):
            tab_end = tab_contig[end] = {}
            for std in ('+','-'):
                tab_std = tab_end[std] = {}
                for letters in LETTERS:
                    tab_std[letters] = dict.fromkeys(keys[end], 0)
    return tab


def initialize_lg():
    tab = {}
    for std in ('+','-'):
        tab[std] = collections.defaultdict(int)
    return tab 


class Options():
    def __init__(self):
        self.length = 70
        self.around = 10

options = Options()


def get_coordinates(read):    
    """ return external coordinates of aligned read bases """
    fivep = read.aend if read.is_reverse else read.pos
    threep = read.pos if read.is_reverse else read.aend
    
    return fivep, threep


def record_lg(read, coordinate, tab):
    """ record global length distribution
    don't record paired reads as they are normally not used for aDNA """
    std = '-' if read.is_reverse else '+'
    
    length = (max(coordinate) - min(coordinate))
    if not read.is_paired:
        tab[std][length] = tab[std][length] + 1
        
    return tab


def get_around(coord, chrom, reflengths, length, ref):
    """ return reference sequences before and after the read
    check for extremities and return what is available """
    coord_min = min(coord)
    coord_max = max(coord)

    pos_before = max(0, coord_min - length)
    pos_after    = min(reflengths[chrom], coord_max + length)

    # Uppercased, to be sure that we don't compare A,C,G,T and a,c,g,t
    before = ref.fetch(chrom, pos_before, coord_min).upper()
    after = ref.fetch(chrom, coord_max, pos_after).upper()

    return before, after


def parse_cigar(cigarlist, ope):
    """ for a specific operation (mismach, match, insertion, deletion... see above)
    return occurences and index in the alignment """
    tlength = 0
    coordinate = []
    # count matches, indels and mismatches
    oplist = (0, 1, 2, 7, 8)
    for operation, length in cigarlist:
        if operation == ope:
                coordinate.append([length, tlength])
        if operation in oplist: 
                tlength += length
    return coordinate


def align(cigarlist, seq, ref):
    """ insert gaps according to the cigar string 
    deletion: gaps to be inserted into read sequences, 
    insertions: gaps to be inserted into reference sequence """    
    lref = list(ref)
    for nbr, idx in parse_cigar(cigarlist, 1):
        lref[idx:idx] = ["-"] * nbr

    lread = list(seq)
    for nbr, idx in parse_cigar(cigarlist, 2):
        lread[idx:idx] = ["-"] * nbr

    return "".join(lread), "".join(lref)


def revcomp(seq):
    """ return reverse complemented string """
    return comp(seq)[::-1]

def comp(seq):
    """ return reverse complemented string """
    return seq.translate(TABLE)

def record_soft_clipping(read, tab, length):
    """ record soft clipped bases at extremities """
    
    def update_table(end, std, bases):
        for i in range(0, min(bases, length)):
            tab[end][std]['S'][i] += 1

    strand = '-' if read.is_reverse else '+'
    for (nbases, idx) in parse_cigar(read.cigar, 4):
        if idx == 0:
            # Soft-clipping at the left side of the alignment
            end = '3p' if read.is_reverse else '5p'
        else:
            # Soft-clipping at the right side of the alignment
            end = '5p' if read.is_reverse else '3p'

        update_table(end, strand, nbases)
        
    return None



def get_mis(read, seq, refseq, ref, length, tab, end):
    """ count mismatches using aligned reference and read,
    must be redone since N in reference were randomly replaced by any bases """
    std = '-' if read.is_reverse else '+'
    subtable = tab[ref][end][std]

    for (i, nt_seq, nt_ref) in zip(range(length), seq, refseq):
        if (nt_seq in "ACGT-") and (nt_ref in "ACGT-"):
            if nt_ref != "-":
                # record base composition in the reference, only A, C, G, T
                subtable[nt_ref][i] += 1

            # Most ref/seq pairs will be identical
            if (nt_ref != nt_seq):
                mut = "%s>%s" % (nt_ref, nt_seq)
                subtable[mut][i] += 1


def _update_table(table, sequence, indices):
    for (index, nt) in zip(indices, sequence):
        if nt in table:
            table[nt][index] += 1


def count_ref_comp(read, chrom, before, after, comp):
    """ record basae composition in external genomic regions """
    std = '-' if read.is_reverse else '+'

    _update_table(comp[chrom]['5p'][std], before, range(-len(before), 0))
    _update_table(comp[chrom]['3p'][std], after,    range(1, len(after) + 1))


def count_read_comp(read, chrom, length, comp):
    """ record base composition of read, discard marked nucleotides """
    std, seq = '+', read.query
    if read.is_reverse:
        std, seq = '-', revcomp(seq)

    _update_table(comp[chrom]['5p'][std], seq, range(1, length + 1))
    _update_table(comp[chrom]['3p'][std], reversed(seq), range(-1, - length - 1, -1))


# http://www.dnabaser.com/articles/IUPAC%20ambiguity%20codes.html
IGNORE_LETTERS = 'YRWSKMDVHBXN'

def compare_str(a, b, ignore_letters=False):
    Nmax = max(len(a), len(b))
    Nmin = min(len(a), len(b))
    
    arr = np.zeros(Nmax, dtype=int)
    for i in range(Nmin):
        if a[i] == b[i] :
            arr[i] = 1 # when equal
        elif ignore_letters and (a[i] in IGNORE_LETTERS or b[i] in IGNORE_LETTERS):
            arr[i] = 1 # when one is an N and is inored
        else:
            arr[i] = 0 # when differing
    arr[Nmin:] = 2 # for values in the end for different read lengths
    return arr

def compare_str_print(a, b, ignore_letters=False):
    arr = compare_str(a, b, ignore_letters)
    res = ''
    for val in arr:
        if val == 1: # if eqal
            res += ' '
        elif val == 0: # when differing
            res += '-'
        elif val == 2: # when too many values
            res += '+'
        else:
            res += str(val)
    return res





# def main():

    # start_time = datetime.now()




# fetch all references and associated lengths in nucleotides
ref_bam = pysam.FastaFile(ref_in) 

# open SAM/BAM file
in_bam = pysam.AlignmentFile(file_bam_in)

N_reads = int(pysam.view('-c', '-F 4', f'{file_bam_in}')) # approximate

reflengths = dict(zip(in_bam.references, in_bam.lengths))
# check if references in SAM/BAM are the same in the fasta reference file
fai_lengths = read_fasta_index(ref_in + ".fai")


# TODO
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
    
    it = tqdm(zip(_read_bamfile(in_bam), f_processed), total=N_reads) if do_tqdm else zip(_read_bamfile(in_bam), f_processed)
        
    # main loop
    for read_bam, line_processed in it:
        counter += 1
        
        # if not read_bam.is_reverse:
        #     continue

        strand_processed, cigar_processed, read_processed, md_tag_processed = line_processed.split()
        
        # md_tag, cigar = md_tag_processed, cigar_processed
        ref_processed, seq_processed = make_reference(read_processed, md_tag_processed, strand_processed, cigar_processed, is_md_cigar=False)
        # ref_processed, seq_processed = ref, seq
        

        # external coordinates 5' and 3' , 3' is 1-based offset
        coordinate = get_coordinates(read_bam)
        # record aligned length for single-end read_bams
        lgdistrib = record_lg(read_bam, coordinate, lgdistrib)
        # fetch reference name, chromosome or contig names
        chrom = in_bam.get_reference_name(read_bam.reference_id)

        

        (before, after) = get_around(coordinate, chrom, reflengths, options.around, ref_bam)
        refseq_bam = ref_bam.fetch(chrom, min(coordinate), max(coordinate)).upper()
        # read_bam.query_alignment_sequence contains aligned sequences while read_bam.seq is the read_bam itself
        seq_bam = read_bam.query_alignment_sequence

        # add gaps according to the cigar string
        (seq_bam, refseq_bam) = align(read_bam.cigar, seq_bam, refseq_bam)

    
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
            print("Seq")
            print(counter)
            print(seq_bam)
            print(seq_processed)
            print(compare_str_print(seq_bam, seq_processed))
            assert False
            
        if (compare_str(refseq_bam, ref_processed, ignore_letters=True) != 1).sum() != 0:
            print("Ref")
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

# print(f'Run completed in {datetime.now() - start_time:f} seconds')

    # return 0



x=x

#%% =============================================================================
# 
# =============================================================================



#%% =============================================================================
# #
# =============================================================================


d_mismatch_forward = {}
d_mismatch_reverse = {}

lengths_forward = []
lengths_reverse = []

list_strand = []


all_seqs = []

print("Parsing MD-tags to get mismatch matrix: ", flush=True)

with open(filename, 'r') as f:
    for iline, line in tqdm(enumerate(f), total=file_len):
        
        strand, cigar, read, md_tag = line.split()
        ref, seq = make_reference(read, md_tag, cigar)

        if 'I' in cigar:
            assert False

        if strand == '0':
            L = parse_md_tag(read, md_tag, cigar, strand, d_mismatch_forward)
            lengths_forward.append(L)
             
        elif strand == '16':
            L = parse_md_tag(read, md_tag, cigar, strand, d_mismatch_reverse)
            lengths_reverse.append(L)
        else:
            print(strand)
            assert False
            
            
        all_seqs.append(read)
        list_strand.append(int(strand))
        
        
list_strand = np.array(list_strand)
lengths_forward = np.array(lengths_forward)        
lengths_reverse = np.array(lengths_reverse)        

df_error_rates_forward = get_error_rates_dataframe(d_mismatch_forward)
df_error_rates_reverse = get_error_rates_dataframe(d_mismatch_reverse)


#%% =============================================================================
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
        # sy = df_error_rates.loc[:, name+'_s']
        # ax.errorbar(x, y, sy, fmt='-', label=strand)
        ax.plot(x, y, '-', label=strand)
    
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





