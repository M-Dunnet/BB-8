#!/usr/bin/env python3
'''
Set of functions for running built-in pairwise concensus
'''
## Installed imports
import mappy

def consensus(sequences, quality_dict, name, seq_type):
    """
    Generate a pairwise consensus sequence from two aligned sequences using quality scores.

    Args:
        sequences: list of two aligned sequences [seqA, seqB]
        quality_dict: dict mapping unaligned sequences to quality strings
        name: read name (for logging)
        seq_type: descriptive type (for logging)

    Returns:
        str: consensus sequence with gaps removed
    """
    seqA, seqB = sequences
    seqA_key, seqB_key = seqA.replace('-', ''), seqB.replace('-', '')

    # Check if sequences exist in quality dict
    missing = [s for s in [seqA_key, seqB_key] if s not in quality_dict]
    if missing:
        print(f"{seq_type} - Missing sequences in quality dict for {name}: {missing}")

    # Normalize quality strings to match aligned sequences
    seqA_qual = normalizeLen(seqA, quality_dict.get(seqA_key, ''))
    seqB_qual = normalizeLen(seqB, quality_dict.get(seqB_key, ''))

    consensus_seq = []
    i = 0
    while i < len(seqA):
        a, b = seqA[i], seqB[i]

        # Case 1: bases match
        if a == b:
            consensus_seq.append(a)

        # Case 2: mismatch, both non-gap
        elif a != '-' and b != '-':
            consensus_seq.append(a if ord(seqA_qual[i]) > ord(seqB_qual[i]) else b)

        # Case 3: gap in one of the sequences
        else:
            gap_seq = seqA if a == '-' else seqB
            gap_len = 1
            while i + gap_len < len(gap_seq) and gap_seq[i + gap_len] == '-':
                gap_len += 1

            # Compare average quality of the gap regions
            avg_a = avgQual(seqA_qual, i, gap_len)
            avg_b = avgQual(seqB_qual, i, gap_len)
            if avg_a > avg_b:
                consensus_seq.append(seqA[i:i + gap_len])
            else:
                consensus_seq.append(seqB[i:i + gap_len])

            i += gap_len
            continue

        i += 1

    return ''.join(consensus_seq).replace('-', '')

def avgQual(qual, i, gapLen):
    '''Returns average quality of a segment.'''
    return sum(ord(x) for x in list(qual[i:i + gapLen])) / gapLen


def normalizeLen(seq, quality):
    '''
    Inserts avg quality scores based on surrounding quality scores where there are gaps in the sequence.
    Returns a new quality string that's the same len as the sequence.
    '''
    seqIndex, qualIndex = 0, 0
    newQuality = ''
    while qualIndex < len(quality):
        if seq[seqIndex] != '-':
            newQuality += quality[qualIndex]
            qualIndex += 1
            seqIndex += 1
        elif seq[seqIndex] == '-' and qualIndex == 0:
            newQuality += quality[qualIndex]
            seqIndex += 1
        else:
            newQuality += chr(int((ord(quality[qualIndex-1]) + ord(quality[qualIndex]))/2))
            seqIndex += 1
    if len(seq) != len(newQuality):
        gapLen = 0
        while seq[-1 - gapLen] == '-':
            newQuality += newQuality[-1]
            gapLen += 1
    return newQuality


def pairwise_consensus(poa_subreads, subreads, sub_quals,name,type1):
    '''
    Runs pairwise concensus calculations. `concensus` function does the heavy lifting
    '''
    seqDict = {}
    for i in range(len(subreads)):
        seqDict[subreads[i]] = sub_quals[i]
    pw_cons = consensus(poa_subreads, seqDict,name,type1)
    return pw_cons
