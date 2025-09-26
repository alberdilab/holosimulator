import argparse
import os
import sys
import subprocess
import yaml
import re
import json
import requests
import pandas as pd
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

def unlock_snakemake(output_dir, profile):
    unlock_command = [
        "/bin/bash", "-c",  # Ensures the module system works properly
        f"module load {config_vars['SNAKEMAKE_MODULE']} && "
        "snakemake "
        f"-s {PACKAGE_DIR / 'workflow' / 'Snakefile'} "
        f"--directory {output_dir} "
        f"--configfile {CONFIG_PATH} "
        f"--workflow-profile {PACKAGE_DIR / 'profile' / profile} "
        f"--unlock"
    ]

    subprocess.run(unlock_command, shell=False, check=True)
    print(f"The output directory {output_dir} has been succesfully unlocked.")

def run_snakemake(workflow, output_dir, profile):
    snakemake_command = [
        "/bin/bash", "-c",
        f"module load {config_vars['SNAKEMAKE_MODULE']} && "
        "snakemake "
        f"-s {PACKAGE_DIR / 'workflow' / 'Snakefile'} "
        f"--directory {output_dir} "
        f"--workflow-profile {PACKAGE_DIR / 'profile' / profile} "
        f"--configfile {CONFIG_PATH} "
        f"--config package_dir={PACKAGE_DIR} workflow={workflow} output_dir={output_dir}"
    ]
    subprocess.run(snakemake_command, shell=False, check=True)

#####
# HoloSimulator execution
#####

def main():
    parser = argparse.ArgumentParser(
        description="HoloSimulator: Python software for simulating holo-omic datasets",
        formatter_class=argparse.RawTextHelpFormatter
    )
    subparsers = parser.add_subparsers(dest="command", help="Available workflows")

    # Arguments for Metagenomics module
    subparser_macro = subparsers.add_parser("metagenomics", help="Upload macro-scale nucleotide data to ENA")
    subparser_macro.add_argument("-i", "--input", required=True, help="Input file")
    subparser_macro.add_argument("-o", "--output", required=True, type=pathlib.Path, help="Output directory")
    subparser_macro.add_argument("-h", "--host", required=True, help="Host genome(s)")
    subparser_macro.add_argument("-m", "--microbiome", required=True, help="Microbial genome(s)")
    subparser_macro.add_argument("-n", "--sample-size", required=True, help="Number of simulated samples")
    subparser_macro.add_argument("-d", "--sequencing-depth", required=True, help="Average sequencing depth per sample")
    subparser_macro.add_argument("-z", "--sequencing-depth-variance", required=True, help="Variance of sequencing depth per sample (Default: 5)")
    subparser_macro.add_argument("-r", "--ratio", required=True, help="Host-microbiome ratio (Default: 0.25 - 20% host, 80% microbiome)")
    subparser_macro.add_argument("-w", "--ratio-variance", required=True, help="Across-sample host-microbiome ratio variance (Default: 5)")
    subparser_macro.add_argument("-v", "--microbiome-variance", required=True, help="Across-sample microbiome variance percentage (Default: 5)")

    # Arguments for Metatranscriptomics module
    subparser_micro = subparsers.add_parser("metatranscriptomics", help="Upload micro-scale nucleotide data to ENA")
    subparser_macro.add_argument("-i", "--input", required=True, help="Input file")
    subparser_macro.add_argument("-o", "--output", required=True, type=pathlib.Path, help="Output directory")
    subparser_macro.add_argument("-h", "--host", required=True, help="Host genome(s)")
    subparser_macro.add_argument("-m", "--microbiome", required=True, help="Microbial genome(s)")
    subparser_macro.add_argument("-n", "--sample-size", required=True, help="Number of simulated samples")
    subparser_macro.add_argument("-r", "--ratio", required=True, help="Host-microbiome ratio (Default: 0.25 - 20% host, 80% microbiome)")
    subparser_macro.add_argument("-w", "--ratio-variance", required=True, help="Across-sample host-microbiome ratio variance (Default: 5)")
    subparser_macro.add_argument("-v", "--microbiome-variance", required=True, help="Across-sample microbiome variance percentage (Default: 5)")

    # Arguments for unlock
    subparser_unlock = subparsers.add_parser("unlock", help="Unlock output directory")
    subparser_unlock.add_argument("-o", "--output", required=False, type=pathlib.Path, default=os.getcwd(), help="Output directory. Default is the directory from which drakkar is called.")

    # Arguments for update
    subparser_update = subparsers.add_parser("update", help="Reinstall HoloSimulator from the Git repo (forces reinstall in this environment)")

    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    if args.command == "macrosample":
        input_dir = args.output / "input"
        input_dir.mkdir(parents=True, exist_ok=True)
        create_secret(args.username, args.password, str(Path(args.output).resolve() / 'input' / '.secret.yml'))
        create_data_dict(args.metadata, args.data, str(Path(args.output).resolve() / 'input' / 'input.json'))
        create_run_checklists(args.metadata, str(Path(args.output).resolve() / 'checklists' / 'run'))
        create_experiment_checklists(args.metadata, str(Path(args.output).resolve() / 'checklists' / 'experiment'))
        create_sample_checklists(args.metadata, str(Path(args.output).resolve() / 'checklists' / 'sample'))
        run_snakemake(args.command, Path(args.output).resolve(), 'slurm')

    if args.command == "microsample":
        input_dir = args.output / "input"
        input_dir.mkdir(parents=True, exist_ok=True)
        create_secret(args.username, args.password, str(Path(args.output).resolve() / 'input' / '.secret.yml'))
        create_data_dict(args.metadata, args.data, str(Path(args.output).resolve() / 'input' / 'input.json'))
        create_run_checklists(args.metadata, str(Path(args.output).resolve() / 'checklists' / 'run'))
        create_experiment_checklists(args.metadata, str(Path(args.output).resolve() / 'checklists' / 'experiment'))
        create_microsample_checklists(args.metadata, str(Path(args.output).resolve() / 'checklists' / 'sample'))
        run_snakemake(args.command, Path(args.output).resolve(), 'slurm')


    ###
    # Unlock or update
    ###

    if args.command == "unlock":
        print(f"{HEADER1}UNLOCKING HOLOSIMULATOR DIRECTORY...{RESET}", flush=True)
        print(f"", flush=True)
        run_unlock(args.command, args.output, args.profile)

    elif args.command == "update":
        pip_cmd = [
                    sys.executable, "-m", "pip", "install",
                    "--upgrade", "--force-reinstall", "--no-deps",
                    "git+https://github.com/alberdilab/holosimulator.git"
                ]
        print(f"{HEADER1}UPDATING HOLOSIMULATOR...{RESET}", flush=True)
        print(f"", flush=True)
        try:
            update_code = subprocess.run(pip_cmd)
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

if __name__ == "__main__":
    main()
