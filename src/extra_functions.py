#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 15 10:51:05 2019

@author: michelsen
"""

import numpy as np
import re
import pandas as pd


ACGT_names = ['A', 'C', 'G', 'T', 'N', '-']
base2index = {val: i for i, val in enumerate(ACGT_names)}

def is_linux():
    import platform
    return platform.system() == 'Linux'


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
    for md_split in md_tags_splitted: 
       
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
    
    is_reverse = (int(strand) & 0x10) != 0
    
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

import re

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


def _read_txtfile(txtfile):
    """
    Takes a subset of the bamfile. Can use a approximate fraction of the
    hits or specific number of reads using reservoir sampling. Returns
    a list in the last case otherwise a iterator.
    
    Skips failed reads (unmapped, secondary, failed qc, duplicates..)
    
    """
    
    for line in txtfile:
        
        strand, cigar, read, md_tag = line.split()
        strand = int(strand)
        
        if not (strand & filtered_flags):
            yield (strand, cigar, read, md_tag)








# from Martin Kircher, to complement DNA
TABLE = str.maketrans('TGCAMRWSYKVHDBtgcamrwsykvhdb', \
                      'ACGTKYWSRMBDHVacgtkywsrmbdhv')

def comp(seq):
    """ return reverse complemented string """
    return seq.translate(TABLE)


CIGAR_REGEX = re.compile("(\d+)([MIDNSHP=XB])")
MDZ_REGEX = re.compile("(\d+)")
CODE2CIGAR= "MIDNSHP=XB"

CIGAR2CODE = dict([y, x] for x, y in enumerate(CODE2CIGAR))
maketrans = TABLE


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

def bam_cigar_op(c):
    return c & BAM_CIGAR_MASK
def bam_cigar_oplen(c):
    return c >> BAM_CIGAR_SHIFT

# // Note that BAM_CIGAR_STR is padded to length 16 bytes below so that
# // the array look-up will not fall off the end.  '?' is chosen as the
# // padding character so it's easy to spot if one is emitted, and will
# // result in a parsing failure (in sam_parse1(), at least) if read.
# def bam_cigar_opchr(c) (BAM_CIGAR_STR "??????" [bam_cigar_op(c)])
# bam_cigar_gen(l, o) ((l)<<BAM_CIGAR_SHIFT|(o))


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
    # if not src:
        # return None

    # md_tag_ptr = bam_aux_get(src, "MD")
    
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
            # for i in range(l-1, -1, -1):
            # # for i from 0 <= i < l:
            #     s[s_idx] = read_sequence[r_idx]
            #     r_idx += 1
            #     s_idx += 1
        
        # op = "D" 
        elif op == BAM_CDEL:
            s += l*'-'
            # s_idx += l
            # for i in range(l-1, -1, -1):
            # # for i from 0 <= i < l:
            #     s[s_idx] = '-'
            #     s_idx += 1
        
        # op = "N"
        elif op == BAM_CREF_SKIP:
            pass
        
        # op = "I"
        elif op == BAM_CINS:
            # encode insertions into reference as lowercase
            s += read_sequence[r_idx:r_idx+l].lower()
            r_idx += l
            # for i in range(l-1, -1, -1):
            # # for i from 0 <= i < l:
            #     # encode insertions into reference as lowercase
            #     s[s_idx] = read_sequence[r_idx].lower()
            #     r_idx += 1
            #     s_idx += 1
        
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
    # return s[:s_idx], s_old[:s_idx]



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
        # c is numerical
        # 0 to 9
        if 48 <= md_tag_ord[md_idx] <= 57:
            nmatches *= 10
            nmatches += md_tag_ord[md_idx] - 48
            md_idx += 1
            continue
        else:
            # save matches up to this point, skipping insertions
             # for x from 0 <= x < nmatches:
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
                    # assert s[s_idx] == '-'
                    # s[s_idx] = chr(md_tag_ord[md_idx])
                    s_insert = chr(md_tag_ord[md_idx])
                    ref_seq = insert_substring_in_string_at_pos(ref_seq, s_insert, s_idx)
                    s_idx += 1
                    md_idx += 1
            else:
                # save mismatch
                # enforce lower case
                # c = md_tag_ord[md_idx]
                # if c <= 90:
                    # c += 32
                # s[s_idx] = chr(md_tag_ord[md_idx]).lower()
                # s = s[:s_idx] + chr(md_tag_ord[md_idx]).lower() + s[s_idx+1:]
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
            # for i in range(l, -1, -1):
            # # for i from 0 <= i < l:
            #     s[s_idx] = cref_seq[r_idx]
            #     r_idx += 1
            #     s_idx += 1
        # D
        elif op == BAM_CDEL:
            # for i from 0 <= i < l:
            s += cref_seq[r_idx:r_idx+l]
            r_idx += l
            # for i in range(l, -1, -1):
            #     s[s_idx] = cref_seq[r_idx]
            #     r_idx += 1
            #     s_idx += 1
        
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




    # def get_reference_sequence(self):
    #     """return the reference sequence in the region that is covered by the
    #     alignment of the read to the reference.
    #     This method requires the MD tag to be set.
    #     """
    #     return force_str(build_reference_sequence(self._delegate))




    # property query_alignment_sequence:
    #     """aligned portion of the read.
    #     This is a substring of :attr:`seq` that excludes flanking
    #     bases that were :term:`soft clipped` (None if not present). It
    #     is equal to ``seq[qstart:qend]``.
    #     SAM/BAM files may include extra flanking bases that are not
    #     part of the alignment.  These bases may be the result of the
    #     Smith-Waterman or other algorithms, which may not require
    #     alignments that begin at the first residue or end at the last.
    #     In addition, extra sequencing adapters, multiplex identifiers,
    #     and low-quality bases that were not considered for alignment
    #     may have been retained.
    #     """

    #     def __get__(self):
    #         if self.cache_query_alignment_sequence:
    #             return self.cache_query_alignment_sequence

    #         cdef bam1_t * src
    #         cdef uint32_t start, end

    #         src = self._delegate

    #         if src.core.l_qseq == 0:
    #             return None

    #         start = getQueryStart(src)
    #         end   = getQueryEnd(src)

    #         self.cache_query_alignment_sequence = force_str(
    #             getSequenceInRange(src, start, end))
    #         return self.cache_query_alignment_sequence



    # property query_sequence:
    #     """read sequence bases, including :term:`soft clipped` bases
    #     (None if not present).
    #     Note that assigning to seq will invalidate any quality scores.
    #     Thus, to in-place edit the sequence and quality scores, copies of
    #     the quality scores need to be taken. Consider trimming for example::
    #        q = read.query_qualities
    #        read.query_squence = read.query_sequence[5:10]
    #        read.query_qualities = q[5:10]
    #     The sequence is returned as it is stored in the BAM file. Some mappers
    #     might have stored a reverse complement of the original read
    #     sequence.
    #     """
    #     def __get__(self):
    #         if self.cache_query_sequence:
    #             return self.cache_query_sequence

    #         cdef bam1_t * src
    #         cdef char * s
    #         src = self._delegate

    #         if src.core.l_qseq == 0:
    #             return None

    #         self.cache_query_sequence = force_str(getSequenceInRange(
    #             src, 0, src.core.l_qseq))
    #         return self.cache_query_sequence







#%% =============================================================================
# 
# =============================================================================


ACGT_names = ['A', 'C', 'G', 'T', 'N', '-', 'Other']
base2index = {val: i for i, val in enumerate(ACGT_names)}


def init_zero_matrix(ACGT_names):
    return np.zeros((len(ACGT_names), len(ACGT_names)), dtype=int)


def fill_mismatch_matrix(seq1, seq2, d_mismatch):
    """ seq1 = ref, seq2 = seq """
    
    for i, (s1, s2) in enumerate(zip(seq1, seq2)):
        i += 1
        if not i in d_mismatch:
            d_mismatch[i] = init_zero_matrix(ACGT_names)
        
        d_mismatch[i][base2index[s1], base2index[s2]] += 1
    
    return d_mismatch




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





#%% =============================================================================
# 
# =============================================================================


# def seq_to_match(read, read_pos, N_soft_clipped_beginning, d_mismatch):
    
#     mask = any([s in ACGT_names for s in read])
#     if not mask:
#         print(read)
#         assert False

#     counter = read_pos + N_soft_clipped_beginning
#     for base in read:
#         index = base2index[base]
#         dict_to_mismatch(d_mismatch, counter)[index, index] += 1
#         counter += 1
        
#     return None
    

# def seq_to_mismatch(read, char, read_pos, N_soft_clipped_beginning, d_mismatch):
    
#     assert len(char) == 1
#     index_ref = base2index[char]
#     index_seq = base2index[read[read_pos]]
#     dict_to_mismatch(d_mismatch, read_pos+N_soft_clipped_beginning)[index_ref, index_seq] += 1
#     return None

# cigar_implemented = list('MISD')
# cigar_all_options = list('MIDNSHP=X')
# cigar_not_implemented = [c for c in cigar_all_options if not c in cigar_implemented]


# def find_all(a_str, sub):
#     start = 0
#     while True:
#         start = a_str.find(sub, start)
#         if start == -1: 
#             return
#         yield start
#         start += len(sub) # use start += 1 to find overlapping matches



# def gen_zero_matrix(ACGT_names):
#     return np.zeros((len(ACGT_names), len(ACGT_names)), dtype=int)


# def dict_to_mismatch(d_mismatch, read_pos):
#     read_pos_1_indexed = read_pos + 1
#     if not read_pos_1_indexed in d_mismatch:
#         d_mismatch[read_pos_1_indexed] = gen_zero_matrix(ACGT_names)
#     return d_mismatch[read_pos_1_indexed]



# def parse_md_tag(read, md_tag, cigar, strand, d_mismatch):
    
#     L = len(read)
    
#     read_pos = 0
#     N_inserts = 0
#     N_soft_clipped_beginning = 0
#     N_soft_clipped_end = 0
#     N_deletions = 0
    
    
#     if 'I' in cigar:
#         cigar_split = list(filter(None, re.split(r'(\d+)', cigar)))
#         for i, split in enumerate(cigar_split):
#             if split == 'I':
#                 N_inserts += int(cigar_split[i-1])
        
#     if 'S' in cigar:
#         cigar_split = list(filter(None, re.split(r'(\d+)', cigar)))
        
#         if cigar_split[1] == 'S' or cigar_split[-1] == 'S':
#             if cigar_split[1] == 'S':
#                 N_soft_clipped_beginning = int(cigar_split[0])
#                 read = read[N_soft_clipped_beginning:]
#             if cigar_split[-1] == 'S':
#                 N_soft_clipped_end = int(cigar_split[-2])
#                 read = read[:-N_soft_clipped_end]
#         else:
#             print("Soft clipping in the middle not implemented yet")
#             assert False 

#     if any([c in cigar_not_implemented for c in cigar]):
#         print(f'Cigar character not implemented yet, {cigar}')
#         assert False
        
    
#     for md_split in filter(None, re.split(r'(\d+)', md_tag)):
#         if md_split == '0':
#             continue
#         elif md_split.isdigit():
#             i = int(md_split)
#             seq_to_match(read[read_pos:read_pos+i], read_pos, N_soft_clipped_beginning, d_mismatch)
#             read_pos += i

#         # ignore inserts for now            
#         elif md_split[0] == '^':
#             # print(f"Ignoring deletions for now: {read}, {md_tag}, {cigar}")
#             N_deletions += len(md_split[1:])

#         else:
#             for char in md_split:
#                 if char.upper() in ACGT_names:
#                     seq_to_mismatch(read, char, read_pos, N_soft_clipped_beginning, d_mismatch)
#                     read_pos += 1
#                 else: 
#                     print(f"The characters {char} is not recognized")
#                     assert False
    
#     assert read_pos == len(read) - N_inserts
    
#     return L + N_deletions

