#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 22 13:41:22 2019

@author: michelsen
"""

import sys
import collections
import numpy as np

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
