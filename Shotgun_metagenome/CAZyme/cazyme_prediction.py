#!/usr/bin/env python
# coding: utf-8

# Manuscript: Distinct gut microbiome characteristics in patients with Parkinson’s disease based on the presence of premotor rapid-eye movement sleep behavior disorders
# Code written by Jae-Yun Lee
# 2024-01-30

# Job description
## Predict carbohydrate-active enzymes of metagenome-assembled genomes.

# Program information
## run_dbCAN (Standalone version of dbCAN3)
## run_dbCAN: https://github.com/linnabrown/run_dbcan
## dbCAN3: Zheng et al., Nucleic acids Research, 2023. doi:10.1093/nar/gkad328

import argparse
import subprocess
import datetime

def run_dbcan(input_path, out_dir, out_prefix, db_dir, input_type, cgc_finder, n_threads):
    
    cmd = ["run_dbcan",
           input_path, input_type,
           "--out_dir", out_dir,
           "--out_pre", f"{out_prefix}.",
           "--db_dir", db_dir,
           "--dia_cpu", str(n_threads),
           "--hmm_cpu", str(n_threads),
           "--dbcan_thread", str(n_threads),
           "--tf_cpu", str(n_threads),
           "--stp_cpu", str(n_threads),
          ]

    if cgc_finder:
        cmd.extend(["-c", "cluster"])

    print(f"Start 'run_dbcan':")
    print(f"\tInput fasta: {input_path}")
    print(f"\tOutput files: {out_dir}/{out_prefix}")
    print(f"\tInput type: {input_type}")
    print(f"\tRun command:\n\t\t{' '.join(cmd)}")

    start_time = datetime.datetime.now()
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        error_message = e.stderr.decode()
        error_log = f"{out_prefix}.error.log"
        with open(error_log, "w") as f:
            f.writelines(error_message)
            
        raise Exception(f"Error occurred: {error_message}")

    end_time = datetime.datetime.now()
    runtime = end_time - start_time
    print(f"Run time: {runtime}-HH:MM:SS.ms\n")

if __name__ == "__main__":
    # Arguments
    parser = argparse.ArgumentParser(description="Run run_dbcan: Standalone version of dbCAN3")
    parser.add_argument("-i", "--input_path", help="Path to the input file")
    parser.add_argument("-e", "--input_type", help="Type of the input file (prok/faa)", default="prok")
    parser.add_argument("-d", "--db_dir", help="Path to the dbcan database directory")
    parser.add_argument("-o", "--out_dir", help="Output directory")
    parser.add_argument("-p", "--out_prefix", help="Output prefix for generated files")
    parser.add_argument("-t", "--n_threads", help="Number of threads", type=int, default=4)
    parser.add_argument("--cgc_finder", help="Enable CGC finder", action='store_true')

    args = parser.parse_args()

    # Run function
    run_dbcan(args.input_path, 
              args.out_dir, 
              args.out_prefix, 
              args.db_dir, 
              args.input_type, 
              args.cgc_finder,
              args.n_threads)

