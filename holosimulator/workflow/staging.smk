configfile: "config.yaml"

import os, json, re, shutil, gzip
from glob import glob

# --- Config ---
INPUT_JSON = config["input"]
OUTDIR     = config["output_dir"]

# --- Parse JSON at parse time ---
with open(INPUT_JSON) as fh:
    META = json.load(fh)

GENOMES = META["genomes"]
IDS        = [g["id"] for g in GENOMES]              # genome IDs like G0001
ID_TO_PATH = {g["id"]: g["path"] for g in GENOMES}

# Paths
GENOME_IN      = os.path.join(OUTDIR, "genomes", gid, "{gid}.fa.gz")
GENOME_OUT      = os.path.join(OUTDIR, "genomes", "{gid}.fa")

rule all:
    input:
        [GENOME_OUT.format(gid=gid) for gid in IDS]

rule stage_genome:
    input:
        INPUT_JSON
    output:
        temp(GENOME_IN)
    threads: 1
    retries: 3
    params:
        src = lambda wc: ID_TO_PATH[wc.gid]
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output})"
        wget --tries=3 --waitretry=5 -O "{output}" "{params.src}"
        """

rule decompress:
    input:
        GENOME_IN
    output:
        GENOME_OUT
    threads: 1
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output})"
        echo "[`date '+%Y-%m-%d %H:%M:%S'`] [Decompress genome] Decompressing {input}"
        gzip -cd "{input}" > "{output}"
        """
