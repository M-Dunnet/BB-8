#!/usr/bin/env python3
'''
Preprocesses reads for peak-finding. Maps splint sequences to the target read for `call_peaks`
'''
## In-built imports
import os
from pathlib import Path
import subprocess

## Installed imports
import mappy


def preprocess(blat, args, tmp_dir, fastq_file):
    """
    Align splints to reads using BLAT, parse results, and return adapter info.
    Returns:
        adapter_dict: {read_name: [adapter, strand]}
        adapter_set: set of all adapters found
        no_splint_reads: number of reads with no matching splint
    """
    ## Make temp files and folders
    tmp_dir = Path(tmp_dir)
    tmp_dir.mkdir(exist_ok=True)
    tmp_fasta = tmp_dir / 'R2C2_temp_for_BLAT.fasta'
    align_psl = tmp_dir / 'splint_to_read_alignments.psl'

    ## Step 1: Write reads to a temporary FASTA and run BLAT
    process(args, fastq_file, blat, tmp_dir)

    ## Step 2: Parse BLAT PSL output
    tmp_adapter_dict = {}
    adapter_set = set()
    with open(align_psl) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            fields = line.split('\t')
            read_name, adapter, strand = fields[9], fields[13], fields[8]
            gaps, score = float(fields[5]), float(fields[0])

            ## Filter poor alignments
            if gaps < 50 and score > 50:
                if read_name not in tmp_adapter_dict:
                    tmp_adapter_dict[read_name] = [[None, 1, None]]  # placeholder for scoring
                tmp_adapter_dict[read_name].append([adapter, score, strand])
                adapter_set.add(adapter)

    ## Step 3: Choose best adapter per read
    adapter_dict = {}
    no_splint_reads = 0
    for read_name, alignments in tmp_adapter_dict.items():
        best = max(alignments, key=lambda x: x[1])
        if not best[0]:  # no valid adapter
            no_splint_reads += 1
            continue
        adapter_dict[read_name] = [best[0], best[2]]
        adapter_set.add(best[0])

    return adapter_dict, adapter_set, no_splint_reads


def process(args, fastq_file, blat, tmp_dir):
    """
    Convert FASTQ/FASTA to temporary FASTA for BLAT, then run BLAT.
    """
    ## Make temp files and folders
    tmp_dir = Path(tmp_dir)
    tmp_dir.mkdir(exist_ok=True)
    tmp_fa = tmp_dir / 'tmp_for_blat.fasta'
    align_psl = tmp_dir / 'splint_to_read_alignments.psl'
    blat_msgs = tmp_dir / 'blat_messages.log'

    # Write reads to temporary FASTA
    with open(tmp_fa, 'w') as out_fh:
        for header, seq, _q in mappy.fastx_read(fastq_file):
            out_fh.write(f'>{header}\n{seq}\n')

    # Run BLAT
    os.system(
        f"{blat} -noHead -stepSize=1 -t=DNA -q=DNA -minScore=15 "
        f"-minIdentity=10 {args.splint_file} {tmp_fa} {align_psl} > {blat_msgs}"
    )

    # Cleanup
    tmp_fa.unlink()