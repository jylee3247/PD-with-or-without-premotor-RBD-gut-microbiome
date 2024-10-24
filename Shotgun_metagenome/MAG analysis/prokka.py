#!/usr/bin/env python
# coding: utf-8

# Manuscript: Distinct gut microbiome characteristics in patients with Parkinson’s disease based on the presence of premotor rapid-eye movement sleep behavior disorders
# Code written by Jae-Yun Lee
# 2024-08-30

# Job description
## Functional annotation of metagenome-assembled genomes reconstructed from current study using Prokka

# Program information
## Prokka version 1.14.6
## Torsten Seemann, Bioinformatics, 2014. doi:10.1093/bioinformatics/btu153

import argparse
import os
import subprocess

def run_prokka(
    input_mag_path,
    out_dir,
    prefix,
    n_threads = 4
):
    print(f"Run prokka.\n Input: {input_mag_path}, Output: {out_dir}/{prefix}")
    cmd = ["prokka", "--outdir", out_dir, "--prefix", prefix, "--locustag", prefix,
           "--cpus", str(n_threads), input_mag_path]

    # Run command
    try:
        subprocess.run(cmd, check = True)
    except subprocess.CalledProcessError as e:
        raise Exception(f"Prokka error.\n\t\t{e.output}")

if __name__ == "__main__":
    # Parameters
    parser.add_argument("-i", "--inputs", type=str, required=True, help="Input assembly file")
    parser.add_argument("-o", "--output_dir", type=str, required=True, help="Output directory")
    parser.add_argument("--out_prefix", type=str, required=True, help="Prefix for the output files")
    parser.add_argument("-t", "--threads", type=int, default=4, help="Number of threads (default: 4)")
    args = parser.parse_args()
    
    # Run Prokka
    run_prokka(args.inputs, args.output_dir, args.out_prefix, args.threads)