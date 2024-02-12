#!/usr/bin/env python
# coding: utf-8

# Manuscript: Distinct gut microbiome characteristics in patients with Parkinson’s disease based on the presence of premotor rapid-eye movement sleep behavior disorders: exploring the alpha-synuclein pathways 
# Code written by Jae-Yun Lee
# 2024-01-30

# Job description
## Binning assembled contigs

# Program information
## MetaBAT2 v2.15
## Kang et al., PeerJ, 2019. doi:10.7717/peerj.7359
## MaxBin2 v2.2.7
## Wu et al., Bioinformatics, 2015. doi:10.1093/bioinformatics/btv638
## CONCOCT v1.1.0
## Alneberg et al., Nature Methods, 2014. doi:10.1038/nmeth.3103

import argparse
import os
import subprocess
from glob import glob
from pathlib import Path

def run_command(cmd, error_head):    
    try:
        subprocess.run(cmd, check = True)
    except subprocess.CalledProcessError as e:
        raise Exception(f"{error_head}: {e.stderr}")

def run_Metabat2(input_contig, depth_info, out_dir, min_contig_size = 1000, saveCls = True, n_threads = 4):
    out_file = os.path.join(out_dir, Path(input_contig).stem)
    cmd = ["metabat2", "-i", input_contig, "-a", depth_info, "-m", str(min_contig_size), "-t", str(n_threads), "--unbinned"]
    if saveCls:
        cmd.extend(["--saveCls"])
    cmd.extend(["-o", out_file])
    
    print(f"METABAT2 binning Run command:\n\t{' '.join(cmd)}")
    run_command(cmd, "METABAT2 binning Error occurred")


def run_Maxbin2(input_contig, abd_info, out_dir, out_prefix, min_contig_size = 1000, n_threads = 4):
    out_file = os.path.join(out_dir, Path(input_contig).stem)
    cmd = ["run_MaxBin.pl", "-contig", input_contig, "-out", out_file, "-abund_list", abd_info, "-min_contig_length", str(min_contig_size), "-thread", str(n_threads)]
    print(f"MAXBIN2 binning Run command:\n\t{' '.join(cmd)}")
    run_command(cmd, "MAXBIN2 binning Error occurred")


def run_CONCOCT(input_contig, depth_info, out_dir, n_threads = 4):
    out_file = os.path.join(out_dir, Path(input_contig).stem)
    cmd = ["concoct", "--threads", str(n_threads), "--coverage_file", depth_info, "--composition_file", input_contig, "-b", out_file]
    print(f"CONCOCT 'concoct' Run command:\n\t{' '.join(cmd)}")
    run_command(cmd, "CONCOCT 'concoct' Error occurred")

def concoct_merge_cutup_clustering(input_dir):
    f_list = glob(os.path.join(input_dir, "*"))
    cluster = [x for x in f_list if re.search("clustering_gt1000.csv$", x)][0]
    out_path = re.sub("gt1000", "merged", cluster)
    cmd = ["merge_cutup_clustering.py", cluster, ">", out_path]
    print(f"CONCOCT 'merge_cutup_clustering.py' Run command:\n\t{' '.join(cmd)}")
    try:
        subprocess.run(" ".join(cmd), shell = True, check = True, capture_output = True)
    except subprocess.CalledProcessError as e:
          raise Exception(f"CONCOCT 'merge_cutup_clustering.py' Error occurred: {e.stderr.decode()}")    

def concoct_extract_fasta_bins(input_contig, cluster, out_dir):
    cmd = ["extract_fasta_bins.py", input_contig, cluster, "--output_path", out_dir]
    print(f"CONCOCT 'extract_fasta_bins.py' Run command:\n\t{' '.join(cmd)}")
    run_command(cmd, "CONCOCT 'extract_fasta_bins.py' Error occurred")

if __name__ == "__main__":
    if program == "metabat":
        run_Metabat2()
        
    elif program == "maxbin":
        run_Maxbin2()
        
    elif program == "concoct":
        run_CONCOCT()
        concoct_merge_cutup_clustering()
        concoct_extract_fasta_bins()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Binning contigs')
    parser.add_argument("-p", '--program', type=str, choices=['metabat', 'maxbin', 'concoct'], help='Binning program to use')
    parser.add_argument("-c", '--input_contig', type=str, help='Path to input contig file')
    parser.add_argument("-d", '--depth_info', type=str, help='Path to depth information file')
    parser.add_argument("-o", '--out_dir', type=str, help='Output directory')
    parser.add_argument('-m', '--min_contig_size', type=int, default=1000, help='Minimum contig size. default: 1000')
    parser.add_argument('-t', '--threads', type=int, default=4, help='Number of threads')
    parser.add_argument('--saveCls', action='store_true', help='Save class information for Metabat2')
    args = parser.parse_args()

    if args.program == "metabat":
        run_Metabat2(args.input_contig, args.depth_info, args.out_dir, args.min_contig_size, args.saveCls, args.threads)
        
    elif args.program == "maxbin":
        run_Maxbin2(args.input_contig, args.depth_info, args.out_dir, args.out_dir, args.min_contig_size, args.threads)
        
    elif args.program == "concoct":
        run_CONCOCT(args.input_contig, args.depth_info, args.out_dir, args.threads)
        concoct_merge_cutup_clustering(args.out_dir)
        concoct_extract_fasta_bins(args.input_contig, os.path.join(args.out_dir, "clustering_merged.csv"), args.out_dir)

