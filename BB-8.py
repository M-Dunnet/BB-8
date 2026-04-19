#!/usr/bin/env python3
"""
Main script for running BB-8
"""

## Built-in imports
import os
import sys
import argparse
import gc
import gzip
import time
import shutil
import multiprocessing as mp
from pathlib import Path

## Installed imports
import numpy as np
import mappy as mm
from conk import conk

## BB-8 module imports
from lib import generateConcensus as gc_md

# -----------------------------
# Globals
# -----------------------------
PATH = '/'.join(os.path.realpath(__file__).split('/')[:-1]) + '/lib/'
sys.path.append(os.path.abspath(PATH))
bb8Path = '/'.join(os.path.realpath(__file__).split('/')[:-1]) + '/'
blat=bb8Path+'/blat/blat'
VERSION = "v1.0"
# -----------------------------

def parse_args():
    '''
    Parses arguments.
    '''
    parser = argparse.ArgumentParser(description='Makes consensus sequences from R2C2 reads.',
                                     add_help=True,
                                     prefix_chars='-')
    parser.add_argument('--reads', '-r', type=str, action='store',
                          help='FASTQ file that contains the long R2C2 reads or a folder containing multiple of these FASTQ files.')

    parser.add_argument('--splint_file', '-s', type=str, action='store',
                          help='Path to the splint FASTA file.')

    parser.add_argument('--out_path', '-o', type=str, action='store', default=os.getcwd(),
                        help='Directory where all the files will end up. Defaults to your current directory.')

    parser.add_argument('--lencutoff', '-l', type=int, action='store', default=1000,
                        help='Sets the length cutoff for your raw sequences. Anything shorter than the cutoff will be excluded. Defaults to 1000.')

    parser.add_argument('--mdistcutoff', '-d', type=int, action='store', default=500,
                        help='Sets the median distance cutoff for consensus sequences. Anything shorter will be excluded. Defaults to 500.')

    parser.add_argument('--numThreads', '-n', type=int, default=1,
                        help='Number of threads to use during multiprocessing. Defaults to 1.')

    parser.add_argument('--compress_output', '-co', action='store_true', default=False,
                        help='Use to compress (gzip) both the consensus fasta and subread fastq output files.')

    parser.add_argument('--peakFinderSettings', '-p', action='store', default='20,3,41,2',
                        help='Defaults to "20,3,41,2". Try "30,3,15,2" for a short splint. No promises though. We only tested BB-8 for splints >100nt')

    parser.add_argument('--resume', '-u', action='store_true', default=False,
                        help='''If set, BB-8 will look for processed.log file in output directory. If processed.log exists, reads marked as 
                                processed in the input will be skipped. Output will be appended to existing output files.''')

    parser.add_argument('--version', '-v', action='version', version=VERSION, help='Prints the BB-8 version.')

    ## Print help text if no arguments are supplied.
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    return parser.parse_args()


def getFileList(query_path, done):
    '''
    Takes a path and returns a list of FASTQ files.
    '''
    file_list=[]
    for file in os.listdir(query_path):
        if 'fastq' in file and 'temp' not in file:  ## Skip temp files produced during batching
            if os.path.abspath(query_path+file) not in done:
                file_list.append(os.path.abspath(query_path+file))

    return file_list


def main(args):
    '''
    Runs the BB-8 pipeline.
    '''
    print(f'BB-8 {VERSION} \nGenerating consensus sequences from R2C2 read data')
    done = set()
    processed_reads_set = set()

    # -----------------------------
    # Handle resume / previous log
    # -----------------------------
    write_mode = 'w'
    processed_log_path = os.path.join(args.out_path, 'processed_reads.log')
    processed_file_log_path = os.path.join(args.out_path, 'processed_files.log')
    ## Read processed_log file and fills processed read set with them to skip later
    if args.resume:
        print(f'\n--resume option is True: Looking for existing log file in {args.out_path}')
        ## Check processed reads
        if os.path.isfile(processed_log_path):
            print('Processed reads log found')
            write_mode = 'a'
            with open(processed_log_path) as f:
                processed_reads_set.update(line.strip() for line in f)
        else:
            print('No processed reads log found')
        
        ## Check processed files
        if os.path.isfile(processed_file_log_path):
            print('Processed files log found')
            write_mode = 'a'
            with open(processed_file_log_path) as f:
                done.update(line.strip() for line in f)
        else:
            print('No processed files log found')

        print(f'{len(processed_reads_set)} processed reads found in log file. They will be skipped\n')
    
    ## Clear up old UMI_Splint folders if not resuming
    else:
        print('\nRemoving old results\n')
        if os.path.isfile(args.splint_file):
            for adapter, seq, q in mm.fastx_read(args.splint_file, read_comment=False):
                adapter_dir = os.path.join(args.out_path, adapter)
                if os.path.isdir(adapter_dir):
                    os.system(f'rm -r {adapter_dir}')

    # -----------------------------
    # Setup log and tmp directories
    # -----------------------------
    log_file_path = os.path.join(args.out_path, 'bb8.log')
    tmp_dir = os.path.join(args.out_path, 'tmp')
    shutil.rmtree(tmp_dir, ignore_errors=True)  ## remove if exists, then recreate empty
    os.makedirs(tmp_dir)

    ## Open the both log files to keep track of progress
    print('Starting consensus calling iteration')
    with open(processed_log_path, write_mode) as processed_reads, \
         open(log_file_path, 'a+') as log_file, \
         open(processed_file_log_path, 'a+') as processed_files:

        log_file.write(f'BB-8 version: {VERSION}\nnew iteration\n')

        input_path = args.reads
        file_list = []
        
        ## Check if input is a folder: 
        if os.path.isdir(input_path):
            print('\tRead input is a folder')
            if not input_path.endswith('/'):
                input_path += '/'
            file_list = getFileList(input_path, done)
        
        ## Check if input is a file:
        elif os.path.isfile(input_path):
            print('\tRead input is a file')
            file_list = [os.path.abspath(input_path)]
        
        ## Explicty say no files were provided
        else:
            print('\tNo file provided')

        print(f'\t{len(file_list)} file(s) provided')
        log_file.write(f'Total files to process: {len(file_list)}\n')

        ## Open a temp folder and file for read bacthing.
        tmp_file_path = os.path.join(tmp_dir, 'tmp_file')
        

        # -----------------------------
        # Process each file iteratively 
        # -----------------------------
        for file in file_list:
            tmp_file = open(tmp_file_path, 'w')
            total_reads = 0
            print(f'\tProcessing file {file}')
            log_file.write(f'Processing file {file}\n')

            ## Iterate through each read in the file; skip if it has been processed.
            for name, seq, q in mm.fastx_read(file, read_comment=False):
                if name in processed_reads_set:
                    continue

                ## Write the read to the tmp file for batching
                total_reads += 1
                tmp_file.write(f'@{name}\n{seq}\n+\n{q}\n')
                
                ## Batch processing every 10k reads
                if total_reads % 10000 == 0:
                    tmp_file.close()

                    processed_tmp = f'{tmp_file_path}_processed'
                    if os.path.exists(processed_tmp):
                        os.remove(processed_tmp)

                    gc_md.generate_concensus(args, tmp_file_path)   ## gc_md is the import name for generate concensus.
                    
                    if not os.path.exists(processed_tmp):
                        raise RuntimeError(f"Consensus step failed: {processed_tmp} not created")
                    
                    ## Write each processed read to the processed file log
                    with open(processed_tmp) as pf:
                        for line in pf:
                            line = line.strip()
                            processed_reads.write(line + '\n')
                            processed_reads_set.add(line)

                    ## Clear temp file
                    tmp_file = open(tmp_file_path, 'w')

            ## Process remaining reads if there is less than 10k left
            tmp_file.close()
            if os.path.getsize(tmp_file_path) > 0:
                processed_tmp = f'{tmp_file_path}_processed'
                if os.path.exists(processed_tmp):
                    os.remove(processed_tmp)

                gc_md.generate_concensus(args, tmp_file_path)
                if not os.path.exists(processed_tmp):
                    raise RuntimeError(f"Consensus step failed: {processed_tmp} not created")
                with open(processed_tmp) as pf:
                    for line in pf:
                        line = line.strip()
                        processed_reads.write(line + '\n')
                        processed_reads_set.add(line)

            ## Add file to set of procssed files so we don't re-run it.
            print("REACHED DONE ADD:", file)
            done.add(file)
            print("REACHED FILE WRITE:", file)
            processed_files.write(file+'\n')
            print("REACHED END OF LOOP:", file)
            print('\n')

        tmp_file.close()
    print('\n')


if __name__ == '__main__':
    ## Spawn workers for multiprocessing and parse arguments
    mp.set_start_method("spawn")
    args = parse_args()
    
    ## Ensure the output location is a folder and make sure it exists
    if not args.out_path.endswith('/'):
        args.out_path += '/'
    if not os.path.exists(args.out_path):
        os.mkdir(args.out_path)
    
    ## Run script.
    main(args)

