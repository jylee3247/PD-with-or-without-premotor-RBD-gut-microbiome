#!/usr/bin/env python
# coding: utf-8

# Manuscript: Distinct gut microbiome characteristics in patients with Parkinson’s disease based on the presence of premotor rapid-eye movement sleep behavior disorders
# Code written by Jae-Yun Lee
# 2024-01-30

# Job description
## Count read coverages of genomic features

# Program information
## featureCounts (v2.0.4) from the Subread package
## featureCounts: Liao et al., Bioinformatics, 2014. doi:10.1093/bioinformatics/btt656
## Subread package: https://subread.sourceforge.net

import argparse
import subprocess
from glob import glob

if __name__ == "__main__":
    # Arguments
    parser = argparse.ArgumentParser(description="Run featureCounts from Subread package")
    parser.add_argument("-g", "--input_gff", required=True, help="Path to the GFF file containing genomic features.")
    parser.add_argument("-b", "--bam_dir", required=True, help="Directory containing BAM files.")
    parser.add_argument("-o", "--out_path", required=True, help="Path to the output file")
    parser.add_argument("-f", "--feature_type", default="CDS", help="Specify feature type.")
    parser.add_argument("-a", "--attribute_type", default="ID", help="Specify attribute type.")
    parser.add_argument("-m", "--min_mapping_quality", type=int, default=0, help="Minimum mapping quality of reads to be counted.")
    parser.add_argument("-t", "--n_threads", type=int, default=4, help="Number of threads.")
    parser.add_argument("--primary", action='store_true', help="Count only primary alignments.")
    parser.add_argument("--paired_end_mode", action='store_true', help="Run in paired-end mode.")
    
    args = parser.parse_args()

    # Make command
    cmd = [
        "featureCounts",
        "-a", args.input_gff,
        "-o", args.out_path,
        "-t", args.feature_type,
        "-g", args.attribute_type,
        "-T", str(args.n_threads),
        "-Q", str(args.min_mapping_quality)
    ]

    if args.paired_end_mode:
        cmd.extend(["-p"])

    if args.primary:
        cmd.extend(["--primary"])

    bam_list = glob(f"{args.bam_dir}/*.bam")
    cmd.extend(bam_list)
    
    print(f"Run command:\n\t{' '.join(cmd)}")

    try:
        subprocess.run(" ".join(cmd), shell = True, check = True)
    except subprocess.CalledProcessError as e:
        raise Exception(f"featureCounts error: {e.stderr.decode()}")

