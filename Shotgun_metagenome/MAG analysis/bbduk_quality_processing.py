#!/usr/bin/env python
# coding: utf-8

# Manuscript: Distinct gut microbiome characteristics in patients with Parkinson’s disease based on the presence of premotor rapid-eye movement sleep behavior disorders
# Code written by Jae-Yun Lee
# 2024-01-30

# Job description
## Quality trimming of raw reads before aseembly-based analysis.

# Program information
## BBTools v39.01
## https://jgi.doe.gov/data-and-tools/software-tools/bbtools/

import argparse
import os
import subprocess
from pathlib import Path

def bbduk_quality_trimming(f_fastq_path, r_fastq_path, out_dir, bbduk_path, qtrim_q, min_read_length, max_n, read_average_q, entropy_cut, ordered):
    
    print("=====Run BBduk from BBtools=====")
    print("\tMode:\tQuality trimming")
    
    f_out_path = os.path.join(out_dir, Path(f_fastq_path).name)
    r_out_path = os.path.join(out_dir, Path(r_fastq_path).name)
    log_path = os.path.join(out_dir, Path(f_fastq_path).stem + ".log")

    in1 = f"in1={f_fastq_path}"
    in2 = f"in2={r_fastq_path}"
    out1 = f"out1={f_out_path}"
    out2 = f"out2={r_out_path}"
    qtrim = "qtrim=r"
    trimq = f"trimq={qtrim_q}"
    min_len = f"minlen={min_read_length}"
    ordered = "ordered"
    maxns = f"maxns={max_n}"
    maq = f"maq={read_average_q}"
    entropy = f"entropy={entropy_cut}"

    cmd = [bbduk_path, in1, in2, out1, out2, qtrim, trimq, min_len, ordered, maxns, maq, entropy]
    if ordered:
        cmd.extend(["ordered"])

    try: 
        res = subprocess.run(cmd, capture_output = True, check = True, text = True)
        with open(log_path, "wt") as f:
            f.write(res.stderr)
    except subprocess.CalledProcessError as e:
        raise Exception(e.stderr)

if __name__ == "__main__":
    # Arguments
    parser = argparse.ArgumentParser(description="Quality trimming using bbduk")
    parser.add_argument("-1", "--forward_fastq_path", help="Path to forward fastq file")
    parser.add_argument("-2", "--reverse_fastq_path", help="Path to reverse fastq file")
    parser.add_argument("-o", "--out_dir", help="Output directory")
    parser.add_argument("--bbduk_path", default="bbduk.sh", help="Path to bbduk.sh")
    parser.add_argument("--qtrim_q", type=int, default=10, help="quality trimming threshold")
    parser.add_argument("--min_read_length", type=int, default=70, help="minimum read length")
    parser.add_argument("--max_n", type=int, default=0, help="maximum number of N bases allowed")
    parser.add_argument("--read_average_q", type=int, default=8, help="minimum average quality")
    parser.add_argument("--entropy_cut", type=float, default=0.90, help="entropy cutoff for trimming")
    parser.add_argument("--ordered", action='store_true', help="retain order of output reads")

    args = parser.parse_args()

    # bbduk
    bbduk_quality_trimming(
        args.foward_fastq_path,
        args.reverse_fastq_path,
        args.out_dir,
        args.bbduk_path,
        args.qtrim_q,
        args.min_read_length,
        args.max_n,
        args.read_average_q,
        args.entropy_cut,
        args.ordered
    )

