#!/usr/bin/env python
# coding: utf-8

# Manuscript: Distinct gut microbiome characteristics in patients with Parkinson’s disease based on the presence of premotor rapid-eye movement sleep behavior disorders
# Code written by Jae-Yun Lee
# 2024-01-30

# Job description
## Check the quality of the metagenome-assembled genomes.

# Program information
## CheckM v1.2.2
## Parks et al., Genome Research, 2015. doi:10.1101/gr.186072.114

import argparse
import subprocess

if __name__ == "__main__":
    # Arguments
    parser = argparse.ArgumentParser(description="Run CheckM lineage-specific workflow")
    parser.add_argument("--bin_input", help="Directory containing bins", required=True)
    parser.add_argument("--out_table", help="Path to output table", required=True)
    parser.add_argument("--output_dir", help="Output directory", required=True)
    parser.add_argument("--n_threads", type=int, help="Number of threads", default=4)
    parser.add_argument("--extension", help="File extension", default="fa")
    
    args = parser.parse_args()

    # Make command
    cmd = ["checkm", "lineage_wf",
           "--file", args.out_table,
           "--extension", args.extension,
           "-t", str(args.n_threads),
           "--pplacer_threads", str(int(args.n_threads/2)),
           "--tab_table",
           args.bin_input, args.output_dir]
    
    try:
        subprocess.run(cmd, check = True)
    except subprocess.CalledProcessError as e:
        raise Exception(f"CheckM error: {e.stderr.decode()}")

