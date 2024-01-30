#!/usr/bin/env python
# coding: utf-8

# Manuscript: Distinct characteristics of gut microbiome in Parkinson’s disease patients with premotor REM sleep behavior disorders: exploring alpha-synuclein pathways
# Authors:
# Code written by Jae-Yun Lee
# 2024-01-30

# Job description
## Mapping metagenomes to contigs using Bowtie2

# Program information
## Bowtie2 version 2.5.1
## Langmead, B., Salzberg, S., Nature Methods, 2012. doi:10.1038/nmeth.1923

import argparse
import os
import subprocess


def run_bowtie2(f_fastq_path, r_fastq_path, bt2db, outpath, n_threads, mode, delete_intermediates):
    print("=====Run Bowtie2=====")
    print(f"\tInput: {" ".join([f_fastq_path, r_fastq_path])}\n\tDB: {bt2db}\n\tMode: {mode}\n\n")

    cmd = ["bowtie2", "-1", f_fastq_path, "-2", r_fastq_path, "-x", bt2db, "-S", outpath, "-p", str(n_threads), mode, "--no-unal"]

    # Run bowtie2
    try:
        subprocess.run(cmd, check = True)
    except subprocess.CalledProcessError as e:
        raise Exception(f"Bowtie2 run error: {e.output}")

# Arguments
parser = argparse.ArgumentParser(description="Run Bowtie2")
parser.add_argument("-1", "--forward_fastq_path", type=str, required = True, help="Input forward fastq path.")
parser.add_argument("-2", "--reverse_fastq_path", type=str, required = True, help="Input reverse fastq path.")
parser.add_argument("-x", "--bt2db", type=str, required = True, help="Path to Bowtie2 database.")
parser.add_argument("-o", "--outpath", type=str, required = True, help="Output file path.")
parser.add_argument("-t", "--n_threads", type=int, default=4, help="Number of threads.")
parser.add_argument("-m", "--mode", type=str, default="--very-sensitive", help="Mode for Bowtie2.")
parser.add_argument("-d", "--delete_intermediates", action='store_true', help="Delete intermediate files.")
args = parser.parse_args()

# Run bowtie2
run_bt2(args.forward_fastq_path, args.reverse_fastq_path, args.bt2db, args.outpath, args.n_threads, args.mode, args.delete_intermediates)