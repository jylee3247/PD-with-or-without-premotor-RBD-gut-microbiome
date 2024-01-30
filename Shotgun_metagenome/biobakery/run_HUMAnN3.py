#!/usr/bin/env python
# coding: utf-8

# Manuscript: Distinct characteristics of gut microbiome in Parkinson’s disease patients with premotor REM sleep behavior disorders: exploring alpha-synuclein pathways
# Authors:
# Code written by Jae-Yun Lee
# 2024-01-30

# Job description
## Functional profiling using HUMAnN3

# Program information
## HUMAnN version 3.6
## Beghini et al., Elife, 2021. doi:10.7554/eLife.65088

import argparse
import os
import subprocess

def run_humann(input_fastq_path, outdir, out_prefix, n_threads, mpa_db, verbose, remove_temp):

    print("=====Run HUMAnN3=====")
    os.makedirs(outdir, exist_ok=True)
    
    # MetaPhlAn option
    metaphlan_options = f'--metaphlan-options=-t rel_ab --bowtie2db {mpa_db}'

    # Write command
    cmd = ["humann", "--input", input_fastq_path, "--output", outdir, "--threads", str(n_threads), metaphlan_options]

    if verbose:
        cmd += ["--verbose"]
        
    if remove_temp:
        cmd += ["--remove-temp-output"]

    # Print current trial info
    print(f"\tCurrent file is {out_prefix}")
    print('\tCommand:', end=' ')
    print(' '.join(cmd), end='\n\n')

    # Run command
    try:
        subprocess.run(cmd, check = True)
    except subprocess.CalledProcessError as e:
        raise Exception(f"Error occurred while running MetaPhlAn.\n\t\t{e.output}")

if __name__ == "__main__":
    # Arguments
    parser = argparse.ArgumentParser(description="Run HUMAnN3.")
    parser.add_argument("-i", "--input_path", type=str, help="Path to input file.")
    parser.add_argument("-o", "--outdir", type=str, help="Output directory.")
    parser.add_argument("--out_prefix", type=str, help="Output file prefix.")
    parser.add_argument("-t", "--n_threads", type=int, default=4, help="Number of processing threads.")
    parser.add_argument("-db", "--mpa_db", type=str, default="DB/mpa_vJan21_CHOCOPhlAnSGB/", help="Path to MetaPhlAn bowtie2 database.")
    parser.add_argument("--verbose", action='store_true', help="Enable verbose output.")
    parser.add_argument("--remove_temp", action='store_true', help="Remove temporary output files.")
    
    args = parser.parse_args()
    
    # Run HUMAnN3
    run_humann(args.input_path, args.outdir, args.out_prefix, args.n_threads, args.mpa_db, args.verbose, args.remove_temp)