#!/usr/bin/env python3
'''
Functions for determining concensus. Includes seperate functions for zero-repeat reads with some overalp and more reads with repeats
'''
## Built-in imports
import sys
from pathlib import Path
import os
import subprocess
import re

## Installed imports
import mappy as mm
import numpy as np
import pyabpoa

## BB-8 imports
from lib.consensus import pairwise_consensus


def determine_consensus(args, read, subreads, sub_qual,
                        tmp_dir):
    """
    Generate consensus sequence for a single read using subreads. One reads with one or more repeat get sent here

    Args:
        args: Namespace with parameters.
        read: tuple (name, seq, qual) of the raw read.
        subreads: list of subread sequences.
        sub_qual: list of subread qualities.
        dangling_subreads: list of partial subreads.
        qual_dangling_subreads: qualities for dangling subreads.
        tmp_dir: Path to temporary directory.

    Returns:
        final_cons: polished consensus sequence (str).
        repeats: number of subreads used.
        subs: list of tuples (subread_name, seq, qual) for all contributing subreads.
    """
    name, seq, qual = read  ## Main read with no splitting
    repeats = len(subreads) ## Number of repeats
    subs = []
    tmp_dir = Path(tmp_dir)
    tmp_files = []

    final_cons = ''
    abpoa_cons = ''

    # -----------------------------
    # Abpoa concensus calling
    # -----------------------------
    ## If one repeat, just take the first subread as concensus
    if repeats == 1:
        abpoa_cons = subreads[0]

    ## if more than 1 repeat, generate a concensus with pyabpoa
    ab = pyabpoa.msa_aligner()
    result = ab.msa(subreads, out_cons=True, out_msa=True, max_n_cons=1)

    ## If two repeats, use inbuilt pair-wise concensus to generate repeat from abpoa alignment output
    if repeats == 2:
        abpoa_cons = pairwise_consensus(result.msa_seq[:-1], subreads, sub_qual, name, 'pair')

    ## If more than two repeats, use the abpoa concensus
    else:
        abpoa_cons = result.cons_seq[0] ## Returned in a list from pyABpoa; want the string
    subs = [(f'{name}_{i+1}', s, q) for i, (s, q) in enumerate(zip(subreads, sub_qual))]
    
    return abpoa_cons, repeats, subs
 

def zero_repeat_cons(args, read, dangling_subreads, qual_dangling_subreads, tmp_dir, read_adapter):
    name, seq, qual = read  ## Main read with no splitting
    subs = []
    tmp_dir = Path(tmp_dir)
    tmp_files = []

    final_cons = ''
    abpoa_cons = ''
    if len(dangling_subreads) == 2:
        tmp_final_cons, tmp_subs, tmp_files = zero_repeats(
            name, seq, qual, dangling_subreads, qual_dangling_subreads, [], tmp_dir, [] 
        )
        if tmp_final_cons and len(tmp_final_cons) >= args.mdistcutoff:
            abpoa_cons = tmp_final_cons
            subs = tmp_subs
    
    ## Prepare final read for FASTA
    avg_qual = round(np.mean([ord(Q)-33 for Q in qual]), 1)
    seq_len = len(seq)
    cons_len = len(abpoa_cons)
    final_cons = f'>{name}_{avg_qual}_{seq_len}_0_{cons_len}\n{abpoa_cons}\n'
    return [(final_cons, read_adapter, subs, 1)] ## '1' is the number of peaks; returned by analyze reads, which called this function   #TODO check this return and see what structure is....


def zero_repeats(name, seq, sub_qual, subreads, sub_qual_list, subs, tmp_dir, tmp_files):
    """
    Handle reads with zero repeats using a pairwise alignment of two dangling subreads.

    Args:
        name: read name
        seq: original read sequence
        subreads: list of two dangling subreads
        sub_qual_list: qualities for dangling subreads
        subs: list of subreads collected so far (will be updated)
        tmp_dir: Path to temporary directory
        tmp_files: list of temporary files (will be updated)

    Returns:
        corrected_cons: consensus sequence (str) or empty if failed
        subs: updated subread tuples
        tmp_files: updated list of temporary files
    """

    ## Add dangling subreads to subs list
    for i, (s, q) in enumerate(zip(subreads, sub_qual_list)):
        subs.append((f'{name}_{i}', s, q))

    ## Pairwise alignment using mappy - This is used to find overlapping seqeunces, but the alignment isnt as good as abpoa?
    mm_aligner = mm.Aligner(seq=subreads[0], preset='map-ont', scoring=(20, 7, 10, 5))
    mappy_res = next(mm_aligner.map(subreads[1]), None)

    if not mappy_res:
        return '', subs, tmp_files

    ## Extract overlapping sequences and flanking regions
    left = subreads[1][:mappy_res.q_st]
    right = subreads[0][mappy_res.r_en:]
    overlap_seq1 = subreads[0][mappy_res.r_st:mappy_res.r_en]
    overlap_qual1 = sub_qual_list[0][mappy_res.r_st:mappy_res.r_en]
    overlap_seq2 = subreads[1][mappy_res.q_st:mappy_res.q_en]
    overlap_qual2 = sub_qual_list[1][mappy_res.q_st:mappy_res.q_en]

    ## get pyabpoa alignment for dangling subreads 
    ab = pyabpoa.msa_aligner(match=5)
    result = ab.msa(
        [overlap_seq1, overlap_seq2],
        out_cons=False,
        out_msa=True,
        max_n_cons=1
    )
    ## extract alignemnt information
    msa_seqs = getattr(result, "msa_seq", None)
    if not msa_seqs:
        return '', subs, tmp_files

    ## Generate concensus using built-in pairwise concensus
    abpoa_cons = pairwise_consensus(msa_seqs, [overlap_seq1, overlap_seq2],
                                    [overlap_qual1, overlap_qual2], name, 'zero')

    ## Combine flanking sequences with consensus
    corrected_cons = left + abpoa_cons + right

    return corrected_cons, subs, tmp_files