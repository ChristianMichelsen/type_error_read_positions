#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 22 13:41:22 2019

@author: michelsen
"""

#%%

import sys
import collections
import numpy as np
import re

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


filtered_flags = (_BAM_UNMAPPED  | \
                  # _BAM_SECONDARY | \
                  _BAM_FAILED_QC | \
                  _BAM_PCR_DUPE  | \
                  _BAM_CHIMERIC)


def _filter_reads(bamfile):
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



def _read_txtfile(txtfile):
    """
    Takes a subset of the bamfile. Can use a approximate fraction of the
    hits or specific number of reads using reservoir sampling. Returns
    a list in the last case otherwise a iterator.
    
    Skips failed reads (unmapped, secondary, failed qc, duplicates..)
    
    """
    
    for line in txtfile:
        
        flag, cigar, read, md_tag = line.split()
        flag = int(flag)
        
        if not (flag & filtered_flags):
            yield (flag, cigar, read, md_tag)

# from Martin Kircher, to complement DNA
TABLE = str.maketrans('TGCAMRWSYKVHDBtgcamrwsykvhdb', \
                      'ACGTKYWSRMBDHVacgtkywsrmbdhv')


CIGAR_REGEX = re.compile("(\d+)([MIDNSHP=XB])")
MDZ_REGEX = re.compile("(\d+)")
CODE2CIGAR= "MIDNSHP=XB"

CIGAR2CODE = dict([y, x] for x, y in enumerate(CODE2CIGAR))
maketrans = TABLE



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

BAM_CIGAR_STR   = "MIDNSHP=XB"
BAM_CIGAR_SHIFT = 4
BAM_CIGAR_MASK  = 0xf
BAM_CIGAR_TYPE  = 0x3C1A7

def cigar_tuples_to_string(cigartuples):
    return "".join([ "%i%c" % (y,CODE2CIGAR[x]) for x,y in cigartuples])

def get_cigar_parts(cigar):
    parts = CIGAR_REGEX.findall(cigar)
    return parts
    
def cigar_string_to_tuples(cigar):
    parts = get_cigar_parts(cigar)
    # reverse order
    cigartuples = [(CIGAR2CODE[y], int(x)) for x,y in parts]
    return cigartuples


def getCigarTuplesLength(cigartuples):
    # return sum([y for (x,y) in cigartuples])
    return get_alignment_length(cigartuples)


def getQueryStart(cigartuples):

    start_offset = 0
    L_max = getCigarTuplesLength(cigartuples)
    
    for op, length in cigartuples:
        if op == BAM_CHARD_CLIP:
            if start_offset != 0 and start_offset != L_max:
                raise ValueError('Invalid clipping in CIGAR string')
        elif op == BAM_CSOFT_CLIP:
            start_offset += length
        else:
            break

    return start_offset


def getQueryEnd(cigartuples):
    
    start_offset = getQueryStart(cigartuples)
    
    
    end_offset = getCigarTuplesLength(cigartuples)
    L_max = end_offset

    # if there is no sequence, compute length from cigar string
    if end_offset == 0:
        for op, length in cigartuples:
            if op == BAM_CMATCH or \
               op == BAM_CINS or \
               op == BAM_CEQUAL or \
               op == BAM_CDIFF or \
              (op == BAM_CSOFT_CLIP and end_offset == 0):
                end_offset += length
    else:
        # walk backwards in cigar string
        for op, length in cigartuples:
            if op == BAM_CHARD_CLIP:
                if end_offset != L_max:
                    raise ValueError('Invalid clipping in CIGAR string')
            elif op == BAM_CSOFT_CLIP:
                end_offset -= length
            else:
                break

    return end_offset + start_offset + 1



def getSequenceInRange(seq, start, end):
    return seq[start:end] 


def get_alignment_length(cigartuples):
    l = 0
    for op, length in cigartuples:
        if op == BAM_CSOFT_CLIP or op == BAM_CHARD_CLIP:
            continue
        l += length
    return l


def get_md_reference_length(md_tag):
    l = 0
    md_idx = 0
    nmatches = 0
    
    md_tag_ord = [ord(s) for s in md_tag]
    md_tag_ord.append(0)

    while md_tag_ord[md_idx] != 0:
        # number 0:9
        if md_tag_ord[md_idx] >= 48 and md_tag_ord[md_idx] <= 57:
            nmatches *= 10
            nmatches += md_tag_ord[md_idx] - 48
            md_idx += 1
            continue
        else:
            l += nmatches
            nmatches = 0
            if md_tag_ord[md_idx] == ord('^'):
                md_idx += 1
                # A to Z
                while md_tag_ord[md_idx] >= 65 and md_tag_ord[md_idx] <= 90:
                    md_idx += 1
                    l += 1
            else:
                md_idx += 1
                l += 1

    l += nmatches
    return l


def get_md_reference_length2(md_tag):
    # is just slightly slower than get_md_reference_length ut easier to read
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
    """return expanded sequence from MD tag.
    The sequence includes substitutions and both insertions in the
    reference as well as deletions to the reference sequence. Combine
    with the cigar string to reconstitute the query or the reference
    sequence.
    Positions corresponding to `N` (skipped region from the reference)
    in the CIGAR string will not appear in the returned sequence. The
    MD should correspondingly not contain these. Thus proper tags are::
       Deletion from the reference:   cigar=5M1D5M    MD=5^C5
       Skipped region from reference: cigar=5M1N5M    MD=10
    Returns
    -------
    None, if no MD tag is present.
    """
    
    parts = get_cigar_parts(cigar)
    cigartuples = cigar_string_to_tuples(cigar)

    start = getQueryStart(cigartuples)
    end = getQueryEnd(cigartuples)
    
    # get read sequence, taking into account soft-clipping
    read_sequence = getSequenceInRange(seq, start, end)
    read_sequence = seq[start:]
    
    # s_idx = 0

    max_len = get_alignment_length(cigartuples)
    if max_len == 0:
        raise ValueError("could not determine alignment length")

    r_idx = 0
    s = ''
    for op, l in cigartuples:
        
        # op = "M", "=", "X"
        if op == BAM_CMATCH or op == BAM_CEQUAL or op == BAM_CDIFF:
            s += read_sequence[r_idx:r_idx+l]
            r_idx += l
        
        # op = "D" 
        elif op == BAM_CDEL:
            s += l*'-'
        
        # op = "N"
        elif op == BAM_CREF_SKIP:
            pass
        
        # op = "I"
        elif op == BAM_CINS:
            # encode insertions into reference as lowercase
            s += read_sequence[r_idx:r_idx+l].lower()
            r_idx += l
        
        # op = "S"
        elif op == BAM_CSOFT_CLIP:
            pass
        
        # op = "H"
        elif op == BAM_CHARD_CLIP:
            pass # advances neither
        
        # op = "P"
        elif op == BAM_CPAD:
            raise NotImplementedError(
                "Padding (BAM_CPAD, 6) is currently not supported. "
                "Please implement. Sorry about that.")

    # Check if MD tag is valid by matching CIGAR length to MD tag defined length
    # Insertions would be in addition to what is described by MD, so we calculate
    # the number of insertions seperately.
    insertions = sum([int(x) for (x,y) in parts if y=='I'])
    

    md_len = get_md_reference_length(md_tag)
    
    # TODO: Remove this check if never called
    assert get_md_reference_length(md_tag) == get_md_reference_length2(md_tag) 
    
    if md_len + insertions > max_len:
        raise AssertionError(f"Invalid MD tag: MD length {md_len} mismatch with CIGAR length {max_len} and {insertions} insertions")

    return s[:max_len]



def insert_substring_in_string_at_pos(s, s_insert, index):
    return s[:index] + s_insert + s[index+1:]


def build_reference_sequence(seq, cigar, md_tag):
    """return the reference sequence in the region that is covered by the
    alignment of the read to the reference.
    This method requires the MD tag to be set.
    """
    
    cigartuples = cigar_string_to_tuples(cigar)
    
    # s_idx = 0
    ref_seq = build_alignment_sequence(seq, cigar, md_tag)

    nmatches = 0
    md_idx = 0
    s_idx = 0

    md_tag_ord = [ord(ss) for ss in md_tag]
    md_tag_ord.append(0)

    while md_tag_ord[md_idx] != 0:
        # c is numerical, 0 to 9
        if 48 <= md_tag_ord[md_idx] <= 57:
            nmatches *= 10
            nmatches += md_tag_ord[md_idx] - 48
            md_idx += 1
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
            if md_tag_ord[md_idx] == ord('^'):
                md_idx += 1
                # A to Z
                while 65 <= md_tag_ord[md_idx] <= 90:
                    s_insert = chr(md_tag_ord[md_idx])
                    ref_seq = insert_substring_in_string_at_pos(ref_seq, s_insert, s_idx)
                    s_idx += 1
                    md_idx += 1
            else:
                # save mismatch, enforce lower case
                s_insert = chr(md_tag_ord[md_idx]).lower()
                ref_seq = insert_substring_in_string_at_pos(ref_seq, s_insert, s_idx)
                s_idx += 1
                md_idx += 1

    s = ''
    
    cref_seq = ref_seq[:]
    r_idx = 0
    for op, l in cigartuples:
        
        # M, =, X
        if op == BAM_CMATCH or op == BAM_CEQUAL or op == BAM_CDIFF:
            s += cref_seq[r_idx:r_idx+l]
            r_idx += l
        # D
        elif op == BAM_CDEL:
            s += cref_seq[r_idx:r_idx+l]
            r_idx += l
        
        # N
        elif op == BAM_CREF_SKIP:
            pass
        
        # I
        elif op == BAM_CINS:
            r_idx += l
        
        # S
        elif op == BAM_CSOFT_CLIP:
            pass
        
        #H
        elif op == BAM_CHARD_CLIP:
            pass # advances neither
        
        # P
        elif op == BAM_CPAD:
            raise NotImplementedError(
                "Padding (BAM_CPAD, 6) is currently not supported. "
                "Please implement. Sorry about that.")

    return s




def align_ref(cigarlist, ref):
    """ insert gaps according to the cigar string 
    deletion: gaps to be inserted into read sequences, 
    insertions: gaps to be inserted into reference sequence """    
    lref = list(ref)
    for nbr, idx in parse_cigar(cigarlist, 1):
        lref[idx:idx] = ["-"] * nbr

    return "".join(lref)



from extra_functions import build_alignment_reference_seq, correct_reverse_strans


def append_corrected_to_filename(filename):
        return filename.replace('.txt', '_corrected.txt')
    
def get_filename_input(filename):
    return f'../data/raw/{filename}'

def get_filename_output(filename):
    return f'../data/processed/{append_corrected_to_filename(filename)}'

def get_filename_processed(filename):
    return get_filename_output(filename)
    

def mapDamage_main_test(refname, filename_bam, filename):

    file_processed_in = get_filename_processed(filename)
    file_bam_in = get_filename_input(filename_bam)
    ref_in = get_filename_input(refname)
    
    
    # from datetime import datetime
    import pysam
    from tqdm import tqdm

    
    # start_time = datetime.now()

    # fetch all references and associated lengths in nucleotides
    ref_bam = pysam.FastaFile(ref_in) 
    
    # open SAM/BAM file
    in_bam = pysam.AlignmentFile(file_bam_in)
    
    
    
    # reflengths = dict(zip(in_bam.references, in_bam.lengths))
    # check if references in SAM/BAM are the same in the fasta reference file
    # fai_lengths = read_fasta_index(ref_in + ".fai")
    
    
    # if not fai_lengths:
    #     return 1
    # elif not compare_sequence_dicts(fai_lengths, reflengths):
    #     return 1
    
    # refnames = in_bam.references
    
    
    
    # for misincorporation patterns, record mismatches
    # misincorp = initialize_mut(refnames, options.length)
    # for fragmentation patterns, record base compositions
    # dnacomp =  initialize_comp(refnames, options.around, options.length)
    # for length distributions
    # lgdistrib =  initialize_lg()
    
    
    # counter = 0
    
    do_tqdm = True
    if do_tqdm:
        N_reads = int(pysam.view('-c', '-F 4', f'{file_bam_in}')) # approximate


    with open(file_processed_in, 'r') as f_processed:
        
        if do_tqdm:
            it = tqdm(zip(_read_bamfile(in_bam), _read_txtfile(f_processed)), total=N_reads) 
        else:
            it = zip(_read_bamfile(in_bam), _read_txtfile(f_processed))
            
        # main loop
        for counter, (read_bam, line_txt) in enumerate(it, start=1):
            # counter += 1
            
            flag, cigar, read, md_tag = line_txt
            # read is equal to read_bam.query_sequence
            # read = read_bam.query_sequence
            
            seq, ref = build_alignment_reference_seq(read, cigar, md_tag)
            is_reverse, ref_processed, seq_processed = correct_reverse_strans(flag, seq, ref)
            
            # external coordinates 5' and 3' , 3' is 1-based offset
            coordinate = get_coordinates(read_bam)
            # record aligned length for single-end read_bams
            # lgdistrib = record_lg(read_bam, coordinate, lgdistrib)
            # fetch reference name, chromosome or contig names
            chrom = in_bam.get_reference_name(read_bam.reference_id)
            
    
            # (before, after) = get_around(coordinate, chrom, reflengths, options.around, ref_bam)
            refseq_bam_no_gaps = ref_bam.fetch(chrom, min(coordinate), max(coordinate)).upper()
            # read_bam.query_alignment_sequence contains aligned sequences while read_bam.seq is the read_bam itself
            seq_bam_no_gaps = read_bam.query_alignment_sequence
    
            # add gaps according to the cigar string
            (seq_bam, refseq_bam) = align(read_bam.cigar, seq_bam_no_gaps, refseq_bam_no_gaps)
    
        
            # reverse complement read_bam and reference when mapped reverse strand
            if read_bam.is_reverse:
                refseq_bam = comp(refseq_bam)
                seq_bam = comp(seq_bam)
        
            
            if not is_reverse:
                if read != seq_bam_no_gaps:
                    
                    # if counter != 3241205:
                
                    print("flag \t\t  cigar \tmdtag")
                    print(f"{flag} \t\t{cigar} \t{md_tag}")
                    
                    print("\nDirect read vs query_alignment_sequence")
                    print(counter)
                    print(read)
                    print(read_bam.query_alignment_sequence)
                    print(compare_str_print(read, read_bam.query_alignment_sequence))
                
                    print("\nProcessed read vs seq_bam (processed)")
                    print(counter)
                    print(seq_processed)
                    print(seq_bam)
                    print(compare_str_print(seq_processed, seq_bam))
                    
                    print("\nReference homemade vs reference fasta")
                    print(counter)
                    print(refseq_bam)
                    print(ref_processed)
                    print(compare_str_print(refseq_bam, ref_processed, ignore_letters=True))
                
                
        
        
        
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
