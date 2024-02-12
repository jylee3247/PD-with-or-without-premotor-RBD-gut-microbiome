#!/usr/bin/env python
# coding: utf-8

# Manuscript: Distinct gut microbiome characteristics in patients with Parkinson’s disease based on the presence of premotor rapid-eye movement sleep behavior disorders: exploring the alpha-synuclein pathways 
# Code written by Jae-Yun Lee
# 2024-01-30

# Job description
## Adapter trimming of raw reads before aseembly-based analysis.

# Program information
## BBTools v39.01
## https://jgi.doe.gov/data-and-tools/software-tools/bbtools/

import argparse
import os
import subprocess
from pathlib import Path

def bbduk_adapter_trimming(f_fastq_path, r_fastq_path, out_dir, bbduk_path, ktrim, k, mink, hdist, tbo, tpe, ordered):
    
    print("=====Run BBduk from BBtools=====")
    print("\tMode:\tAdapter trimming")
    f_out_path = os.path.join(out_dir, Path(f_fastq_path).name)
    r_out_path = os.path.join(out_dir, Path(r_fastq_path).name)
    log_path = os.path.join(out_dir, Path(f_fastq_path).stem + ".log")
    in1 = f"in1={f_fastq_path}"
    in2 = f"in2={r_fastq_path}"
    out1 = f"out1={f_out_path}"
    out2 = f"out2={r_out_path}"
    min_len = "minlen=70"
    ref = "ref=adapters"
    ktrim = f"ktrim={r}"
    k = f"k=23"
    mink = f"mink={mink}"
    hdist = f"hdist={hdist}"

    cmd = [bbduk_path, in1, in2, out1, out2, min_len, ref, ktrim, k, mink, hdist, tbo, tpe, ordered]

    if tbo:
        cmd.extend(["tbo"])
    if tpe:
        cmd.extend(["tpe"])
    if ordered:
        cmd.extend(["ordered"])
    
    try: 
        res = subprocess.run(cmd, capture_output=True, check=True, text=True)
        with open(log_path, "wt") as f:
            f.write(res.stderr)
    except subprocess.CalledProcessError as e:
        raise Exception(e.stderr)

if __name__ == "__main__":
    # Arguments
    parser = argparse.ArgumentParser(description="Adapter trimming using bbduk")
    parser.add_argument("-1", "--forward_fastq_path", help="Path to forward fastq file")
    parser.add_argument("-2", "--reverse_fastq_path", help="Path to reverse fastq file")
    parser.add_argument("-o", "--out_dir", help="Output directory")
    parser.add_argument("--bbduk_path", default="bbduk.sh", help="Path to bbduk.sh")
    parser.add_argument("--min_read_length", type=int, default=70, help="minimum read length")
    parser.add_argument("--ktrim", default="r", help="direction of k-mer trimming")
    parser.add_argument("--k", type=int, default=23, help="k-mer size")
    parser.add_argument("--mink", type=int, default=11, help="minimum k-mer number")
    parser.add_argument("--hdist", type=int, default=1, help="hamming distance")
    parser.add_argument("--tbo", action='store_true', help="trim adapters based on pair overlap detection")
    parser.add_argument("--tpe", action='store_true', help="trim both reads to the same length")
    parser.add_argument("--ordered", action='store_true', help="retain order of output reads")

    args = parser.parse_args()

    # Run bbduk
    bbduk_adapter_trimming(
        args.forward_fastq_path,
        args.reverse_fastq_path,
        args.out_dir,
        args.bbduk_path,
        args.min_read_length,
        args.ktrim,
        args.k,
        args.mink,
        args.hdist,
        args.tbo,
        args.tpe,
        args.ordered
    )

