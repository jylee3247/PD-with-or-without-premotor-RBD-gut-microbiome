#!/usr/bin/env python
# coding: utf-8

# Manuscript: Distinct gut microbiome characteristics in patients with Parkinson’s disease based on the presence of premotor rapid-eye movement sleep behavior disorders: exploring the alpha-synuclein pathways 
# Code written by Jae-Yun Lee
# 2024-01-30

# Job description
## Taxonomy classification of metagenome-assembled genomes.

# Program information
## GTDB-Tk v2.2.6
## Chaumeil et al., Bioinformatics, 2022. doi:10.1093/bioinformatics/btac672

import argparse
import subprocess

if __name__ == "__main__":
    # Arguments
    parser = argparse.ArgumentParser(description="Run GTDB-Tk classify workflow")
    parser.add_argument("--genome_dir", help="Directory containing genomes", required=True)
    parser.add_argument("--out_dir", help="Output directory", required=True)
    parser.add_argument("--mash_db", help="Path to mash database", required=True)
    parser.add_argument("--extension", help="File extension for genome files", default="fa")
    parser.add_argument("--n_threads", type=int, help="Number of threads", default=4)
    args = parser.parse_args()
    
    # Make command
    cmd = ["gtdbtk", "classify_wf",
           "--genome_dir", args.genome_dir,
           "--out_dir", args.out_dir,
           "--mash_db", args.mash_db,
           "--extension", args.extension,
           "--cpus", str(args.n_threads),
           "--pplacer_cpus", str(int(args.n_threads/2)),
           "--keep_intermediates",
           "--write_single_copy_genes"]

    # Run command
    try:
        subprocess.run(cmd, check = True)
    except subprocess.CalledProcessError as e:
        raise Exception(f"GTDBtk Error: {e.stderr.decode()}")

