#!/usr/bin/env python
# coding: utf-8

# Manuscript: Distinct gut microbiome characteristics in patients with Parkinson’s disease based on the presence of premotor rapid-eye movement sleep behavior disorders
# Code written by Jae-Yun Lee
# 2024-01-30

# Job description
## Refining bins from metagenome

# Program information
## DAS tool v1.1.6
## Sieber et al., Nature microbiology, 2018. doi:10.1038/s41564-018-0171-1.

import argparse
import os
import subprocess
import re
from glob import glob
from pathlib import Path

def find_binner_and_extension(in_dir, valid_binners):
    path_list = sorted(glob(os.path.join(in_dir, "*")))
    binner_found = [binner for binner in valid_binners if any(binner in path for path in path_list)]
    if len(binner_found) > 1:
        raise Exception(f"find_binner error. Find more than one mathced binner in path_list: {binner_found} are detected.\n\tPlease check name of bins in {os.path.dirname(path_list[0])}")
    binner = binner_found[0]
    extension = path_list[0].split(".")[-1]
    
    return binner, extension

def dastool_fasta_to_contig2bin(out_dir, in_dir):
    valid_binners = ["metabat2", "maxbin2", "concoct"]
    binner, extension = find_binner_and_extension(in_dir, valid_binners)
    cmd = ["Fasta_to_Contig2Bin.sh", "-e", extension, "-i", in_dir, ">", os.path.join(out_dir, f"{binner}_contigs2bin.tsv")]
    print(f"\tCommand: {' '.join(cmd)}")
    try:
        subprocess.run(' '.join(cmd), shell = True, check = True)
    except subprocess.CalledProcessError as e:
        raise Exception(f"Fasta_to_Contig2Bin.sh error occurred: {e.stderr}")
        
    return binner, extension

def run_dastool(contigs2bin, input_contig, out_path, n_threads):

    labels = [re.sub("_contigs2bin\.tsv$", "", x) for x in contigs2bin]
    
    cmd = ["DAS_Tool", 
           "-i", ",".join(contigs2bin), 
           "-c", input_contig,
           "-l", ",".join(labels),
           "-o", out_path,
           "-t", str(n_threads),
           "--write_bin_evals",
           "--write_bins",
           "--write_unbinned"
          ]
    
    print("Run DAS_Tool")
    print(f"\tInputs:\n\t\tBins: {','.join(contigs2bin)}\n\t\tContig: {input_contig}\n\t\tOutput: {out_path}")
    try:
        subprocess.run(cmd, check = True)
    except subprocess.CalledProcessError as e:
        raise Exception(f"DAS Tool error occurred: {e.stderr.decode()}")

if __name__ == "__main__":
    # Arguments
    parser = argparse.ArgumentParser(description='Run DAS Tool for bin refinement')
    parser.add_argument("-b", '--bin_dir', type=str, required=True, help='Directories containing bins to run DAS Tool')
    parser.add_argument("-c", '--input_contig', type=str, required=True, help='Path to input contig file')
    parser.add_argument("-o", '--out_dir', type=str, required=True, help='Output directory')
    parser.add_argument('-t', '--threads', type=int, default=4, help='Number of threads')
    args = parser.parse_args()

    # Generate contigs2bin mappings
    contigs2bin_files = []
    bin_dir_list = glob(f"{args.bind_dir}/*")
    for bin_dir in bin_dir_list:
        binner, ext = dastool_fasta_to_contig2bin(args.out_dir, bin_dir)
        contigs2bin_files.append(os.path.join(args.out_dir, f"{binner}_contigs2bin.tsv"))
        
    # Run DAS_tool
    out_path = os.path.join(args.out_dir, "refined_bins")
    os.makedirs(out_path, exist_ok=True)
    run_dastool(contigs2bin_files, args.input_contig, out_path, args.threads)
