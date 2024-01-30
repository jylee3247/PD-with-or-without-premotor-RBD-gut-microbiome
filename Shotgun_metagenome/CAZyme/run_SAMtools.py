#!/usr/bin/env python
# coding: utf-8

# Manuscript: Distinct characteristics of gut microbiome in Parkinson’s disease patients with premotor REM sleep behavior disorders: exploring alpha-synuclein pathways
# Authors:
# Code written by Jae-Yun Lee
# 2024-01-30

# Job description
## Filtering SAM file and convert to fastq file

# Program information
## SAMtools version 1.16.1
## Petr Danecek et al., GigaScience, 2021. doi: 10.1093/gigascience/giab008

import argparse
import os
import subprocess
import re
from pathlib import Path

def run_command(cmd, error_head):    
    try:
        subprocess.run(cmd, check = True)
    except subprocess.CalledProcessError as e:
        raise Exception(f"{error_head}: {e.stderr}")

def run_samtools(input_path, outdir, n_threads, mode):
    print("Run SAMtools.")
    print(f"Mode: {mode}")
    os.makedirs(outdir, exist_ok=True)

    if mode == "remove_host_read":
        # Convert SAM to BAM and filtering out mapped reads
        out_bam_path = os.path.join(outdir, re.sub("\.sam$", "\.bam$", os.path.basename(input_path)))
        cmd = ["samtools",
                   "view", "-bS", "-f", "12", "-F", "256", "--threads", str(n_threads), "-o", out_bam_path, input_path]
        run_command(cmd, "SAMtools SAM to BAM conversion error")

        # SAMtools sorting BAM
        sort_bam_path = os.path.join(outdir, re.sub("\.bam$", "\.sorted.bam$", os.path.basename(out_bam_path)))    
        cmd = ["samtools", "sort", "-n", "--threads", str(n_threads), "-o", sort_bam_path, out_bam_path]
        run_command(cmd, "SAMtools BAM sorting error")

        # SAMtools convert BAM to fastq
        fastq_prefix = Path(input_path).stem
        f_fastq_path = os.path.join(outdir, fastq_prefix + "_1.fastq.gz")
        r_fastq_path = os.path.join(outdir, fastq_prefix + "_2.fastq.gz")
    
        cmd = ["samtools", "fastq", "--threads", str(n_threads), "-1", f_fastq_path, "-2", r_fastq_path, sort_bam_path]
        run_command(cmd, "SAMtools BAM to FASTQ conversion error")
            
    elif mode == "contig_mapping":
        # Convert SAM to BAM
        out_bam_path = os.path.join(outdir, re.sub("sam$", "bam$", os.path.basename(input_path)))
        cmd = ["samtools", "view", "-b", "--threads", str(n_threads), "-o", out_bam_path, input_path]
        run_command(cmd, "SAMtools SAM to BAM conversion error")
    
        # SAMtools sorting BAM
        sort_bam_path = os.path.join(outdir, re.sub("bam$", "sorted.bam$", os.path.basename(out_bam_path)))
        cmd = ["samtools", "sort", "--threads", str(n_threads), "-o", sort_bam_path, out_bam_path]
        run_command(cmd, "SAMtools BAM sorting error")
    else:
        raise Exception("Mode is not selected. Mode options: 'remove_host_read', 'contig_mapping'")

if __name__ == "__main__":
    # Arguments
    parser = argparse.ArgumentParser(description="Run SAMtools")
    parser.add_argument("-i", "--input_path", type=str, required=True, help="Input SAM file path.")
    parser.add_argument("-o", "--outdir", type=str, required=True, help="Output directory.")
    parser.add_argument("-t", "--n_threads", type=int, default=4, help="Number of threads.")
    parser.add_argument("-m", "--mode", type=str, choices=["remove_host_read", "contig_mapping"], help="Mode of operation for SAMtools")
    args = parser.parse_args()
    
    # Run SAMtools
    run_samtools(args.input_path, args.outdir, args.n_threads, args.mode)