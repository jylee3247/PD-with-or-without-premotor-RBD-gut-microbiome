#!/usr/bin/env python
# coding: utf-8

# Manuscript: Distinct characteristics of gut microbiome in Parkinson’s disease patients with premotor REM sleep behavior disorders: exploring alpha-synuclein pathways
# Authors:
# Code written by Jae-Yun Lee
# 2024-01-30

# Job description
## Artifact filtering of raw reads before aseembly-based analysis.

# Program information
## BBTools v39.01
## https://jgi.doe.gov/data-and-tools/software-tools/bbtools/

import argparse
import os
import subprocess
from pathlib import Path

def bbduk_artifacts_filtering(f_fastq_path, r_fastq_path, out_dir, bbduk_path):
    
    print("=====Run BBduk from BBtools=====")
    print("\tMode:\tArtifacts,PhiX filtering")
    f_out_path = os.path.join(out_dir, Path(f_fastq_path).name)
    r_out_path = os.path.join(out_dir, Path(r_fastq_path).name)
    log_path = os.path.join(out_dir, Path(f_fastq_path).stem + ".log")

    in1 = f"in1={f_fastq_path}"
    in2 = f"in2={r_fastq_path}"
    out1 = f"out1={f_out_path}"
    out2 = f"out2={r_out_path}"
    k = "k=31"
    ref = "ref=artifacts,phix"
    ordered = "ordered"
    cardinality = "cardinality"

    cmd = [bbduk_path, in1, in2, out1, out2, ref, k, ordered, cardinality]
    try: 
        res = subprocess.run(cmd, capture_output=True, check=True, text=True)
        with open(log_path, "wt") as f:
            f.write(res.stderr)
    except subprocess.CalledProcessError as e:
        raise Exception(e.stderr)


if __name__ == "__main__":
    # Arguments
    parser = argparse.ArgumentParser(description="Artifact filtering using bbduk")
    parser.add_argument("-1", "--forward_fastq_path", help="Path to forward fastq file")
    parser.add_argument("-2", "--reverse_fastq_path", help="Path to reverse fastq file")
    parser.add_argument("-o", "--out_dir", help="Output directory")
    parser.add_argument("--bbduk_path", default="bbduk.sh", help="Path to bbduk.sh")

    args = parser.parse_args()

    # Run bbduk
    bbduk_artifacts_filtering(
        args.foward_fastq_path,
        args.reverse_fastq_path,
        args.out_dir,
        args.bbduk_path
    )

