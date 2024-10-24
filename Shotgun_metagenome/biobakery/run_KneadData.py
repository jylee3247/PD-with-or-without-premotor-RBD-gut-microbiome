#!/usr/bin/env python
# coding: utf-8

# Manuscript: Distinct gut microbiome characteristics in patients with Parkinson’s disease based on the presence of premotor rapid-eye movement sleep behavior disorders
# Code written by Jae-Yun Lee
# 2024-01-30

# Job description
## Read pre-processing before using biobkaery programs.

# Program information
## KneadData version 0.12.0
## https://huttenhower.sph.harvard.edu/kneaddata/

import argparse
import subprocess

# KneadData function
def run_kneaddata(input1_path, input2_path, 
                  output_path, bt2_db_path, 
                  trimmomatic_path, threads, remove_intermediate):

    print("=====Run KneadData=====")
    process = int(threads/2) # Number of core

    # Default options
    default_options = {
        "input1": None,
        "input2": None,
        "output": None,
        "db": None,
        "run_fastqc_start": True,
        "run_fastqc_end": True,
        "verbose": True,
        "trimmomatic": None,
        "threads": None,
        "process": None,
        "remove_intermediate": None
    }

    # Set parameters
    default_options["input1"] = input1_path
    default_options["input2"] = input2_path
    default_options["output"] = output_path
    default_options["db"] = bt2_db_path
    default_options["trimmomatic"] = trimmomatic_path
    default_options["threads"] = threads
    default_options["process"] = process
    default_options["remove_intermediate"] = remove_intermediate
    input1 = default_options["input1"]
    input2 = default_options["input2"]
    output = default_options["output"]
    db = default_options["db"]
    run_fastqc_start = default_options["run_fastqc_start"]
    run_fastqc_end = default_options["run_fastqc_end"]
    verbose = default_options["verbose"]
    trimmomatic = default_options["trimmomatic"]
    threads = str(default_options["threads"])
    process = str(default_options["process"])

    # KneadData command
    cmd = ["kneaddata", "-db", db, "--output", output, "--trimmomatic", trimmomatic, "-t", threads, "-p", process,]

    if input2 is not None:
        cmd += ["--input1", input1, "--input2", input2]
    else:
        cmd += ["--unpaired", input1]   

    if verbose:
        cmd += ["--verbose"]

    if run_fastqc_start:
        cmd += ["--run-fastqc-start"]

    if run_fastqc_end:
        cmd += ["--run-fastqc-end"]
        
    if remove_intermediate:
        cmd += ["--remove-intermediate-output"]

    print('Command:', end=' ')
    print(' '.join(cmd), end='\n\n')

    # Run command
    try:
        subprocess.run(cmd, check = True)
    except subprocess.CalledProcessError as e:
        raise Exception(f"Error occurred while running kneaddata.\n\t\t{e.output}")

if __name__ == "__main__":
    # Arguments
    parser = argparse.ArgumentParser(description="Run kneadData")
    parser.add_argument("-i1", "--input1_path", type=str, help="Path to forward fastq file.")
    parser.add_argument("-i2", "--input2_path", type=str, help="Path to reverse fastq file.")
    parser.add_argument("-o", "--output_path", type=str, help="Output directory path")
    parser.add_argument("-db", "--bt2_db_path", type=str, default="DB/GRCh38", help="Bowtie2 database path.")
    parser.add_argument("-t", "--threads", type=int, default=4, help="Number of threads")
    parser.add_argument("--trimmomatic_path", type=str, default="Trimmomatic-0.39", help="Trimmomatic path.")
    parser.add_argument("--remove_intermediate", action='store_true', help="Remove intermediate files.")
    
    args = parser.parse_args()
    
    # Run KneadData
    run_kneaddata(args.input1_path, args.input2_path, args.output_path, args.bt2_db_path, args.trimmomatic_path, args.threads, args.remove_intermediate)
