import argparse
import os
import sys
import subprocess
import yaml
import re
import json
import requests
import random
import pandas as pd
from pathlib import Path
import pathlib
from collections import defaultdict
from holosimulator.utils import *

#####
# HoloSimulator installation path
#####

PACKAGE_DIR = Path(__file__).parent
CONFIG_PATH = PACKAGE_DIR / "workflow" / "config.yaml"

def load_config():
    """Load fixed variables from config.yaml."""
    if CONFIG_PATH.exists():
        with open(CONFIG_PATH, "r") as f:
            return yaml.safe_load(f)
    return {}

config_vars = load_config()

###
# Define text colors
###

HEADER1 = "\033[1;95m"
ERROR = "\033[1;31m"
INFO = "\033[1;34m"
RESET = "\033[0m"

#####
# Function definitions
#####

def run_unlock(module, output_dir, threads):

    unlock_command = [
        "/bin/bash", "-c",  # Ensures the module system works properly
        "snakemake "
        f"-s {PACKAGE_DIR / 'bin' / 'unlock.smk'} "
        f"--directory {output_dir} "
        f"--configfile {CONFIG_PATH} "
        f"--config cores={threads}"
        f"--unlock "
    ]

    subprocess.run(unlock_command, shell=False, check=True)
    print(f"The output directory {output_dir} has been succesfully unlocked")

def run_staging(module, output_dir, threads, input, host, microbiome):
    snakemake_command = [
        "/bin/bash", "-c",
        "snakemake "
        f"-s {PACKAGE_DIR / 'workflow' / 'staging.smk'} "
        f"--directory {output_dir} "
        f"--cores {threads} "
        f"--configfile {CONFIG_PATH} "
        f"--config package_dir={PACKAGE_DIR} module={module} output_dir={output_dir} input={input} host={host}"
    ]
    subprocess.run(snakemake_command, shell=False, check=True)

def run_genomics(module, output_dir, threads, input, sequencing_model, seed):
    snakemake_command = [
        "/bin/bash", "-c",
        "snakemake "
        f"-s {PACKAGE_DIR / 'workflow' / 'genomics.smk'} "
        f"--directory {output_dir} "
        f"--cores {threads} "
        f"--configfile {CONFIG_PATH} "
        f"--config package_dir={PACKAGE_DIR} module={module} output_dir={output_dir} input={input} sequencing_model={sequencing_model} seed={seed}"
    ]
    subprocess.run(snakemake_command, shell=False, check=True)

def run_transcriptomics(module, output_dir, threads, input, sequencing_model, seed):
    snakemake_command = [
        "/bin/bash", "-c",
        "snakemake "
        f"-s {PACKAGE_DIR / 'workflow' / 'transcriptomics.smk'} "
        f"--directory {output_dir} "
        f"--cores {threads} "
        f"--configfile {CONFIG_PATH} "
        f"--config package_dir={PACKAGE_DIR} module={module} output_dir={output_dir} input={input} sequencing_model={sequencing_model} seed={seed}"
    ]
    subprocess.run(snakemake_command, shell=False, check=True)

#####
# HoloSimulator execution
#####

def main():
    parser = argparse.ArgumentParser(
        description="HoloSimulator: flexible holo-omic dataset simulator",
        formatter_class=argparse.RawTextHelpFormatter
    )
    subparsers = parser.add_subparsers(dest="module", help="Available workflows")

    # Arguments for Hologenomics module
    DEFAULT_GENOMICS = PACKAGE_DIR / "bin" / "default_genomics.csv"
    subparser_genomics = subparsers.add_parser("genomics", help="Simulate hologenomics reads")
    subparser_genomics.add_argument("-i", "--input", default=DEFAULT_GENOMICS, required=False, help="Input file")
    subparser_genomics.add_argument("-o", "--output", required=False, default=os.getcwd(), type=pathlib.Path, help="Working directory. Default is the directory from which HoloSimulator is called.")
    subparser_genomics.add_argument("-e", "--host", required=False, help="Host genome(s)")
    subparser_genomics.add_argument("-m", "--microbiome", required=False, help="Microbial genome(s)")
    subparser_genomics.add_argument("-n", "--sample-size", dest="sample_size", default=1, required=False, help="Number of simulated samples")
    subparser_genomics.add_argument("-d", "--sequencing-depth", dest="sequencing_depth", default=3000000, required=False, help="Average sequencing depth per sample (Default: 3000000")
    subparser_genomics.add_argument("-z", "--sequencing-depth-variance", dest="sequencing_depth_variance", default=2, required=False, help="Variance of sequencing depth per sample (Default: 2)")
    subparser_genomics.add_argument("-f", "--host-fraction", dest="host_fraction", default=0.2, required=False, help="Fraction of host DNA (Default: 0.2)")
    subparser_genomics.add_argument("-w", "--host-fraction-variance", dest="host_fraction_variance", default=2, required=False, help="Across-sample host-microbiome ratio variance (Default: 5)")
    subparser_genomics.add_argument("-v", "--microbiome-variance", dest="microbiome_variance", default=2, required=False, help="Across-sample microbiome variance percentage (Default: 5)")
    subparser_genomics.add_argument("-q", "--sequencing-model", dest="sequencing_model", default="hiseq", required=False, help="Sequencing model for ISS (Default: HiSeq)")
    subparser_genomics.add_argument("-s", "--seed", required=False, type=int, default=random.randint(0, 9999), help="Random seed for reproducibility. If not set, results will vary across runs ")   
    subparser_genomics.add_argument("-t", "--threads", default=1, help="Number of threads to use (Default: 1)")   

    # Arguments for Holotranscriptomics module
    DEFAULT_TRANSCRIPTOMICS = PACKAGE_DIR / "bin" / "default_transcriptomics.csv"
    subparser_transcriptomics = subparsers.add_parser("transcriptomics", help="Simulate holotranscriptomics reads")
    subparser_transcriptomics.add_argument("-i", "--input", default=DEFAULT_TRANSCRIPTOMICS, required=False, help="Input file")
    subparser_transcriptomics.add_argument("-o", "--output", required=False, default=os.getcwd(), type=pathlib.Path, help="Working directory. Default is the directory from which HoloSimulator is called.")
    subparser_transcriptomics.add_argument("-e", "--host", required=False, help="Host genome(s)")
    subparser_transcriptomics.add_argument("-m", "--microbiome", required=False, help="Microbial genome(s)")
    subparser_transcriptomics.add_argument("-n", "--sample-size", dest="sample_size", required=False, help="Number of simulated samples")
    subparser_transcriptomics.add_argument("-d", "--sequencing-depth", dest="sequencing_depth", default=3000000, required=False, help="Average sequencing depth per sample (Default: 3000000")
    subparser_transcriptomics.add_argument("-z", "--sequencing-depth-variance", dest="sequencing_depth_variance", default=2, required=False, help="Variance of sequencing depth per sample (Default: 2)")
    subparser_transcriptomics.add_argument("-f", "--host-fraction", dest="host_fraction", default=0.2, required=False, help="Fraction of host DNA (Default: 0.2)")
    subparser_transcriptomics.add_argument("-w", "--host-fraction-variance", dest="host_fraction_variance", default=2, required=False, help="Across-sample host-microbiome ratio variance (Default: 5)")
    subparser_transcriptomics.add_argument("-v", "--microbiome-variance", dest="microbiome_variance", required=False, default=2, help="Across-sample microbiome variance percentage (Default: 5)")
    subparser_transcriptomics.add_argument("-q", "--sequencing-model", dest="sequencing_model", default="hiseq", required=False, help="Sequencing model for ISS (Default: HiSeq)")
    subparser_transcriptomics.add_argument("-s", "--seed", required=False, type=int, default=random.randint(0, 9999), help="Random seed for reproducibility. If not set, results will vary across runs")   
    subparser_transcriptomics.add_argument("-t", "--threads", default=1, help="Number of threads to use (Default: 1)")   

    # Arguments for unlock
    subparser_unlock = subparsers.add_parser("unlock", help="Unlock output directory")
    subparser_unlock.add_argument("-o", "--output", required=False, type=pathlib.Path, default=os.getcwd(), help="Output directory. Default is the directory from which drakkar is called.")

    # Arguments for update
    subparser_update = subparsers.add_parser("update", help="Reinstall HoloSimulator from the Git repo (forces reinstall in this environment)")

    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    ###
    # Unlock or update
    ###

    if args.module == "unlock":
        print(f"{HEADER1}Unlocking HoloSimulator Directory...{RESET}", flush=True)
        print(f"", flush=True)
        run_unlock(args.module, args.output, args.threads)

    elif args.module == "update":
        pip_cmd = [
                    sys.executable, "-m", "pip", "install",
                    "--upgrade", "--force-reinstall", "--no-deps",
                    "git+https://github.com/alberdilab/holosimulator.git"
                ]
        print(f"{HEADER1}Updating HoloSimulator...{RESET}", flush=True)
        print(f"", flush=True)
        try:
            update_code = subprocess.run(pip_cmd)
            print(f"", flush=True)
            print(f"{HEADER1}HoloSimulator was succesfully updated!{RESET}", flush=True)
        except Exception as e:
            print(f"Update failed: {e}", file=sys.stderr, flush=True)
            sys.exit(1)

    else:            
        project_name = os.path.basename(os.path.normpath(args.output))

        # Check if directory is locked
        if is_snakemake_locked(args.output):
            print(f"{ERROR}Error: directory is locked by Snakemake{RESET}", file=sys.stderr)
            print(f"This usually happens when a workflow is interrupted suddenly,")
            print(f"often while Slurm jobs are still running. Make sure all pending")
            print(f"jobs are finished or killed before unlocking the working directory.")
            print(f"{INFO}Run 'holosimulator unlock' to unlock the working directory{RESET}")
            sys.exit(2)

    ###
    # Run pipeline
    ###

    if not args.module in ("unlock", "update"):

        GENOMES_JSON = "genomes.json"
        GENES_JSON = "genes.json"
        GENOMES_JSON.parent.mkdir(parents=True, exist_ok=True) 

        if (args.host and str(args.host).strip()) or (args.microbiome and str(args.microbiome).strip()):
            args_to_genomics_json(
                host=args.host,
                microbiome=args.microbiome,
                sample_size=args.sample_size,
                sequencing_depth=args.sequencing_depth,
                sequencing_depth_variance=args.sequencing_depth_variance,
                host_fraction=args.host_fraction,
                host_fraction_variance=args.host_fraction_variance,
                microbiome_variance=args.microbiome_variance,
                seed=args.seed,
                output_json=GENOMES_JSON)
        else:
            csv_to_inputs_json(args.input,GENOMES_JSON)

        print(f"{HEADER1}Staging reference genomes...{RESET}", flush=True)

        # Check genomes and yield errors if necessary
        bad = {p: ok for p, ok in check_genomics_paths(GENOMES_JSON).items() if not ok}
        if bad:
            print(f"{ERROR}Some genome paths are invalid:{RESET}", flush=True)
            for p in bad:
                print(" -", p)
            raise SystemExit(1)

        run_staging(
            args.module, 
            Path(args.output).resolve(), 
            args.threads, 
            GENOMES_JSON, 
            args.host, 
            args.microbiome)

        if args.module == "genomics":            
            print(f"{HEADER1}Simulating reads...{RESET}", flush=True)
            run_genomics(
                args.module, 
                Path(args.output).resolve(), 
                args.threads, 
                GENOMES_JSON, 
                args.sequencing_model,
                args.seed)

        if args.module == "transcriptomics":
            print(f"{HEADER1}Calculating read allocation...{RESET}", flush=True)
            genomics_to_transcriptomics_json(
                GENOMES_JSON, 
                GENES_JSON, 
                Path(args.output).resolve(),
                allocation="tpm",
                allocation_params={
                    "read_length": 150,
                    "ln_mu": 0.0, "ln_sigma": 1.2,
                    "gc_bias_strength": 2.0, "gc_bias_linear": 0.0,
                    "overdispersion": 50.0,      # smaller => noisier
                    "dropout_rate": 0.01,
                    "min_reads": 0,
                },
                pairs_to_reads=2)

            print(f"{HEADER1}Simulating reads...{RESET}", flush=True)
            run_transcriptomics(
                args.module, 
                Path(args.output).resolve(), 
                args.threads, 
                GENES_JSON, 
                args.sequencing_model,
                args.seed)
        ###
        # Run snakemake
        ###



if __name__ == "__main__":
    main()
