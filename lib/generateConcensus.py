#!/usr/bin/env python3
'''
Set of functions to generate concensus called reads from `call_peaks` derived subreads.

generate_concensus
│
├── preprocess (adapter + splint assignment)    // External import
│
├── multiprocessing pool
│   └── analyze_reads (per read)                // Internal Function
│       │
│       ├── conk scoring                        // External import
│       ├── call_peaks (repeat boundaries)      // External import
│       ├── extract subreads
│       ├── cluster_splint_umis                 // Internal Function
│       │   └── find_variable_region            // Internal Function
│       └── determine_consensus                 // External import
│
└── write outputs (consensus + subreads + logs)
'''

## Built-in imports
import os
import sys
import argparse
import multiprocessing as mp
import gc
import gzip
import time
from pathlib import Path
import warnings

## Installed imports
import mappy as mm
from conk import conk
import numpy as np
import edlib
with warnings.catch_warnings():
    warnings.filterwarnings("ignore", message="networkx backend defined more than once: nx-loopback")
    import networkx as nx

## BB-8 module imports
PATH = '/'.join(os.path.realpath(__file__).split('/')[:-1]) + '/bin/'
sys.path.append(os.path.abspath(PATH))
from lib.preprocess import preprocess
from lib.call_peaks import call_peaks
from lib.determine_consensus import determine_consensus, zero_repeat_cons

## Define paths for BB-8.py and blat
bb8Path = '/'.join(os.path.realpath(__file__).split('/')[:-1]) + '/'
blat=bb8Path+'/blat/blat'


def rounding(x, base):
    '''Rounds to the nearest base, we use 50'''
    return int(base * round(float(x) / base))


def find_variable_region(seq, max_len=20):
    '''
    Identify a short variable region (splint UMI) between predefined anchor sequences and infer strand orientation.
    '''
    # Anchors
    left_anchor_R  = "TAAGTAGTGAT"
    right_anchor_R = "TATGGAACTCAT"

    left_anchor_F  = "AAATGCACTA"
    right_anchor_F = "GCGATCGAAAAT"

    # --- Reverse strand (R/R) ---
    left_res = edlib.align(left_anchor_R, seq, task="locations", mode="HW")
    right_res = edlib.align(right_anchor_R, seq, task="locations", mode="HW")

    if left_res["locations"] and right_res["locations"]:
        l_start, l_end = left_res["locations"][0]
        r_start, r_end = right_res["locations"][0]

        if l_end < r_start:
            var = seq[l_end:r_start]
            if len(var) <= max_len:
                return var, "R"

    # --- Forward strand (F/F) ---
    left_res = edlib.align(left_anchor_F, seq, task="locations", mode="HW")
    right_res = edlib.align(right_anchor_F, seq, task="locations", mode="HW")

    if left_res["locations"] and right_res["locations"]:
        l_start, l_end = left_res["locations"][0]
        r_start, r_end = right_res["locations"][0]

        if l_end < r_start:
            var = seq[l_end:r_start]
            if len(var) <= max_len:
                return var, "F"

    return None, None


def cluster_splint_umis(subreads):
    '''
    Cluster subreads into groups based on similarity of extracted UMI-like variable regions using edit distance.
    '''

    g = nx.Graph()
    for i, subread in enumerate(subreads, start = 0):
        region = subread[50:100]  # your window
        var_region, strand = find_variable_region(region)
        
        if var_region:
            g.add_node(i, umi=var_region)   ## Add the index. this lets us pull out subreads and their qual from two lists by index.
    
    nodes = list(g.nodes(data=True))
    for i in range(len(nodes)):
        node_i, data_i = nodes[i]
        umi_i = data_i['umi']

        for j in range(i + 1, len(nodes)):
            node_j, data_j = nodes[j]
            umi_j = data_j['umi']

            result = edlib.align(umi_i, umi_j, mode="NW", task="distance")
            dist = result['editDistance']

            if dist <= 2:
                g.add_edge(node_i, node_j, weight=dist)

    clusters = list(nx.connected_components(g))

    return clusters


def analyze_reads(args, read, splint, read_adapter, tmp_dir):
    '''
    Process a single R2C2 read to generate consensus sequences.
    '''

    ## Unpack parameters
    name, seq, qual = read
    seq_len = len(seq)
    use = False
    read_consensus = ''
    subs = []
    
    ## Peak calling parameters
    penalty, iters, window, order = map(int, args.peakFinderSettings.split(','))

    ## Generate scores and peaks
    scores = conk.conk(splint, seq, penalty)
    peaks = call_peaks(scores, args.mdistcutoff, iters, window, order)

    ## If any peaks were found:
    if list(peaks):
        ## shift peaks by half splint length, means subread split happens in the middle of the splint
        peaks = np.array(peaks) + len(splint)//2
        peaks = peaks[peaks < seq_len]  ## remove peaks beyond sequence length

        if len(peaks) > 0:
                use = True

    subreads, qual_subreads, dangling_subreads, qual_dangling_subreads, cons_lst = [], [], [], [], []
    ## Return empty list if there are no peaks
    if use == False:
        return cons_lst ## TODO Possible format error here
    
    ## Split reads if use is True
    if len(peaks) > 1:
        ## Fill subread lists up; default processing if more than one peak
        subread_lens = np.diff(peaks)
        subread_lens_rounded = [rounding(x, 50) for x in subread_lens]
        median_len = np.median(subread_lens_rounded)

        for i, (start, end) in enumerate(zip(peaks[:-1], peaks[1:])):
            if median_len*0.8 <= subread_lens_rounded[i] <= median_len*1.2:
                subreads.append(seq[start:end])
                qual_subreads.append(qual[start:end])

    else:
        # Single peak => everything is dangling
        dangling_subreads = [seq[:peaks[0]], seq[peaks[0]:]]
        qual_dangling_subreads = [qual[:peaks[0]], qual[peaks[0]:]]
        return zero_repeat_cons(args, read, dangling_subreads, qual_dangling_subreads, tmp_dir, read_adapter)   ## TODO Possible format error here

    ## Cluster subreads based on splint UMI 
    splint_clusters = cluster_splint_umis(subreads)
    for i, umi_group in enumerate(splint_clusters):
        if len(umi_group)/len(subreads) <= 0.25: ## Skip clusters that represent less or equal than 1/4 of all subreads. TODO need better solution...
            continue

        ## Group subreads by cluster
        subread_group = [subreads[idx] for idx in umi_group]
        qual_subread_group = [qual_subreads[idx] for idx in umi_group]

        ## Generate consensus per cluster
        consensus, repeats, subs = determine_consensus(
            args, read, subread_group, qual_subread_group,
            tmp_dir
        )

        ## Compile concensus data and append to cons_lst; which is returned. 
        ## Appending to a list is used, becasue based on the number of clusteres, there can be more than one concensus read per orginal ONT read
        if consensus:
            avg_qual = round(np.mean([ord(Q)-33 for Q in qual]), 1) ## TODO qual needs to be sorted out
            cons_len = len(consensus)
            read_consensus = f'>{name}-{i}_{avg_qual}_{seq_len}_{repeats}_{cons_len}\n{consensus}\n'
            cons_lst.append((read_consensus, read_adapter, subs, peaks))    ## Appended as a tuple; unpacked by `generate_concensus`
    
    ## Return consensus list
    return cons_lst


def create_files(adapter,args,outDict,outSubDict,outCountDict):
    '''
    Creates output files per splint that is found and writes the concensus results within those files
    '''
    ## Create splint folder for current splint
    outCountDict[adapter]=0
    if not os.path.exists(args.out_path + adapter):
        os.mkdir(args.out_path + adapter)
    outCountDict[adapter]=0
    
    ## Create output files; Compress output if set
    if args.compress_output:
        writeMode='ab+'
        outDict[adapter]=gzip.open(args.out_path + adapter +'/R2C2_Consensus.fasta.gz',writeMode)
        outSubDict[adapter]=gzip.open(args.out_path + adapter +'/R2C2_Subreads.fastq.gz',writeMode)
    
    ## Create output files uncompressed
    else:
        writeMode='a'
        outDict[adapter]=open(args.out_path + adapter +'/R2C2_Consensus.fasta',writeMode)
        outSubDict[adapter]=open(args.out_path + adapter +'/R2C2_Subreads.fastq',writeMode)
    
    return outDict,outSubDict,outCountDict


def generate_concensus(args, tmp_file_path):
    """
    Generate consensus sequences from long, concatamer R2C2 reads
    Args:
        args: Namespace with parameters (out_path, splint_file, numThreads, lencutoff, compress_output, etc.)
        tmp_file_path: path to input FASTA/FASTQ file
    """
    start_time = time.time()

    ## Make temp folder for BLAT outputs (mapping splint to full length reads)
    out_path = Path(args.out_path)
    tmp_dir = out_path / 'tmp'
    tmp_dir.mkdir(exist_ok=True)  ## Ensure tmp dir exists

    log_path = out_path / 'bb8.log'
    processed_path = tmp_dir / f'{Path(tmp_file_path).name}_processed'

    ## Prepare splints
    splint_dict = {}
    for name, seq, _q in mm.fastx_read(args.splint_file, read_comment=False):
        splint_dict[name] = [seq, mm.revcomp(seq)]

    ## Preprocess adapters and create output files; preprocess maps splint to full length reads
    adapter_dict, adapter_set, _ = preprocess(blat, args, str(tmp_dir), tmp_file_path)
    outDict, outSubDict, outCountDict = {}, {}, {}
    for adapter in adapter_set:
        outDict, outSubDict, outCountDict = create_files(adapter, args, outDict, outSubDict, outCountDict)

    # Counters
    total_reads = short_reads = no_splint_reads = consNumber = 0
    results = {}

    ## Open log and processed files
    with open(log_path, 'a+') as log_file, open(processed_path, 'w') as processed_file, mp.Pool(args.numThreads) as pool:

        ## Submit all reads to pool
        for name, seq, q in mm.fastx_read(tmp_file_path, read_comment=False):
            total_reads += 1
            
            ## Skip if reads are too short
            if len(seq) < args.lencutoff:
                short_reads += 1
                processed_file.write(f'{name}\n')
                continue

            ## Skip if no mapping to splint
            if name not in adapter_dict:
                no_splint_reads += 1
                processed_file.write(f'{name}\n')
                continue

            ## Define the splint for this read
            adapter_name, strand = adapter_dict[name][:2]
            splint = splint_dict[adapter_name][1 if strand == '-' else 0]
            
            ## Run analyze reads to generate concensus
            results[name] = pool.apply_async(
                analyze_reads,
                [args, [name, seq, q], splint, adapter_name, str(tmp_dir)]
            )
            # results[name] = analyze_reads(
            #     args, [name, seq, q], splint, adapter_name, str(tmp_dir)
            # )

        pool.close()
        pool.join()
        gc.collect()


        ## Process and write to disk concensus results of current bacth
        for name, async_result in results.items():
            for cons_idx in async_result.get():   ## Add .get() when using multiprocessing
                consensus, adapter, subs, _peaks = cons_idx
                if consensus:
                    if args.compress_output:
                        consensus = consensus.encode()
                    outDict[adapter].write(consensus)
                    outCountDict[adapter] += 1
                    consNumber += 1

                    for subname, subseq, subq in subs:
                        entry = f'@{subname}\n{subseq}\n+\n{subq}\n'
                        if args.compress_output:
                            entry = entry.encode()
                        outSubDict[adapter].write(entry)

            processed_file.write(f'{name}\n')

        ## Write summary to log
        log_file.write(f'Too short reads: {short_reads} ({short_reads/total_reads*100:.2f}%)\n')
        log_file.write(f'No splint reads: {no_splint_reads} ({no_splint_reads/total_reads*100:.2f}%)\n')
        log_file.write(f'Successful consensus reads: {consNumber} ({consNumber/total_reads*100:.2f}%)\n')

        ## Close adapter files and log per-adapter counts
        for adapter in adapter_set:
            outDict[adapter].close()
            outSubDict[adapter].close()
            with open(log_path, 'a+') as log_file:
                log_file.write(f'\t{outCountDict[adapter]} consensus reads generated for {adapter}\n')

    ## Print time for this batch
    duration = time.time() - start_time
    print(f'Finished generating {consNumber:,} consensus sequences from {total_reads:,} raw reads '
            f'({round(consNumber/total_reads*100)}%) in {round(duration/60,5)} minutes.')

