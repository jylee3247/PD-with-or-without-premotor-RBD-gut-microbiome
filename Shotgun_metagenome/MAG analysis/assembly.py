#!/usr/bin/env python
# coding: utf-8

# Manuscript: Distinct gut microbiome characteristics in patients with Parkinson’s disease based on the presence of premotor rapid-eye movement sleep behavior disorders
# Code written by Jae-Yun Lee
# 2024-01-30

# Job description
## Assemble metagenome to contigs.

# Program information
## SPAdes v3.15.5
## Nurk et al., Genome research, 2017. doi:10.1101/gr.213959.116

import argparse
import os
import subprocess
import datetime

def run_metaspades(forward_fastq_path, reverse_fastq_path, out_dir, n_threads = 4):
    print("=====Run SPAdes assembler=====")
    print("Mode: '--meta'")
    os.makedirs(out_dir, exist_ok=True)

    print(f"Current inputs are {forward_fastq_path}, {reverse_fastq_path}")
    print(f"Outputs will be saved in {out_dir}")
    print("\tStart running...")
    start_time = datetime.datetime.now()
    cmd = ["spades.py", "--meta", "-1", forward_fastq_path, "-2", reverse_fastq_path, 
               "-t", str(n_threads), "-o", out_dir]

    try:
        subprocess.run(cmd, check = True)
    except subprocess.CalledProcessError as e:
        raise Exception(e.stderr)
    end_time = datetime.datetime.now()
    runtime = end_time - start_time
    print(f"\tSPAdes finished. RUNTIME: {runtime}-HH:MM:SS.ms\n\n")

if __name__ == "__main__":
    # Arguments
    parser = argparse.ArgumentParser(description='Run SPAdes assembler for metagenomic data.')
    parser.add_argument("--1", '--forward_fastq', type=str, help='Path to forward fastq file')
    parser.add_argument("--2", '--reverse_fastq', type=str, help='Path to reverse fastq file')
    parser.add_argument("--o", 'output_dir', type=str, help='Output directory')
    parser.add_argument('-t', '--threads', type=int, default=4, help='Number of threads')
    args = parser.parse_args()
    
    # Run SPAdes
    run_metaspades(args.forward_fastq, args.reverse_fastq, args.output_dir, args.threads)

