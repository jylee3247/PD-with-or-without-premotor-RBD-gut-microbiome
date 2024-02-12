#!/usr/bin/env python
# coding: utf-8

# Manuscript: Distinct gut microbiome characteristics in patients with Parkinson’s disease based on the presence of premotor rapid-eye movement sleep behavior disorders: exploring the alpha-synuclein pathways 
# Code written by Jae-Yun Lee
# 2024-01-30

# Job description
## Trimming PCR primer sequences from the 16S rRNA gene amplicon raw data.

# Program information
## Cutadapt 4.0
## Cutadapt: Marcel Martin, EMBnet.journal, 2011. doi:10.14806/ej.17.1.200 (https://cutadapt.readthedocs.io/en/stable/index.html)

import argparse
import os
import re
import subprocess
from pathlib import Path

def run_cutadapt(forward_read, reverse_read, outdir, f_primer_seq, r_primer_seq):

    acc_num = re.sub("_1\.fastq\.gz$", "", Path(forward_read).name)
    print(f'Accession number: {acc_num}')
    print(f'\tForward Primer sequence: {f_primer_seq}')
    print(f'\tReverse Primer sequence: {r_primer_seq}')

    cmd = [
        'cutadapt',
        '--cores', '0',
        '--error-rate', '0.1',
        '--times', '1',
        '--overlap', '3',
        '--minimum-length', '1',
        '--json', os.path.join(outdir, f"{acc_num}.cutadapt.json"),
        '-o', os.path.join(outdir, Path(forward_read).name),
        "-p", os.path.join(outdir, Path(reverse_read).name),
        '-g', f_primer_seq,
        '-G', r_primer_seq,
        "--match-read-wildcards",
        "--discard-untrimmed",
        forward_read,
        reverse_read
    ]

    print('Command:', ' '.join(cmd), '\n')

    try:
        out = subprocess.run(cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        raise Exception(f"Error occurred: {e.stderr}")
        
    log_name = os.path.join(outdir, f"{acc_num}.log")
    with open(log_name, "wt") as f:
        f.write(out.stdout)

if __name__ == "__main__":
    # Arguments
    parser = argparse.ArgumentParser(description="Run cutadapt for primer trimming")
    parser.add_argument("-1", "--forward_fastq", type=str, required=True, help="Forward fastq file path")
    parser.add_argument("-2", "--reverse_fastq", type=str, required=True, help="Reverse fastq file path")
    parser.add_argument("-o", "--output_dir", type=str, required=True, help="Output directory path")
    parser.add_argument("-f", "--forward_primer", type=str, required=True, default="CCTACGGGNGGCWGCAG", help="Forward primer sequence")
    parser.add_argument("-r", "--reverse_primer", type=str, required=True, default="GACTACHVGGGTATCTAATCC", help="Reverse primer sequence")
    args = parser.parse_args()

    # Run cutadapt function
    run_cutadapt(args.forward_fastq, args.reverse_fastq, args.outdir, args.forward_primer, args.reverse_primer)

    # Output files had been imported into the Qiime2 platform (https://qiime2.org/) and underwent further analysis as described in Methods section of our manuscript.