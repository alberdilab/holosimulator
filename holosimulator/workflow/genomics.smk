configfile: "config.yaml"

import os, json, re, shutil, gzip
from glob import glob

# --- Config ---
INPUT_JSON = config["input"]
OUTDIR     = config["output_dir"]
THREADS      = config["threads"]
SEQUENCING_MODEL = config["sequencing_model"]
SEED = config["seed"]
GENE_ALLOC = config["gene_allocation"]

# --- Parse JSON at parse time ---
with open(INPUT_JSON) as fh:
    META = json.load(fh)

SAMPLES = META["samples"]                            # e.g. ["Sample1", ...]
GENOMES = META["genomes"]                            # list of dicts with id, organism, path, abundances

# Maps & lists
TAXA       = [g["organism"] for g in GENOMES]        # organism names (not used in paths)
IDS        = [g["id"] for g in GENOMES]              # genome IDs like G0001
ORG_TO_ID  = {g["organism"]: g["id"] for g in GENOMES}
ID_TO_ORG  = {g["id"]: g["organism"] for g in GENOMES}
ID_TO_PATH = {g["id"]: g["path"] for g in GENOMES}

# Abundances as PAIRS per (organism, sample); we convert to reads later
ABUND = {g["organism"]: {s: int(g["abundances"][s]) for s in SAMPLES} for g in GENOMES}

# Paths
GENOME_OUT      = os.path.join(OUTDIR, "genomes", "{gid}.fa.gz")
ALL_GENOME_FILES = [GENOME_OUT.format(gid=gid) for gid in IDS]

rule all:
    input:
        INPUT_JSON,
        [os.path.join(OUTDIR, f"{s}_1.fq.gz") for s in SAMPLES],
        [os.path.join(OUTDIR, f"{s}_2.fq.gz") for s in SAMPLES]

rule simulate:
    input:
        os.path.join(OUTDIR, "genomes", "{gid}.fa")
    output:
        r1 = temp(os.path.join(OUTDIR, "iss/{sample}/{gid}/reads_R1.fastq")),
        r2 = temp(os.path.join(OUTDIR, "iss/{sample}/{gid}/reads_R2.fastq"))
    params:
        model    = SEQUENCING_MODEL,
        seed     = SEED,
        organism = lambda w: ID_TO_ORG[w.gid],
        # reads = 2 × pairs from JSON
        nreads   = lambda w: 2 * int(ABUND[ID_TO_ORG[w.gid]][w.sample])
    threads: 1
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.r1})"

        if [ "{params.nreads}" -eq 0 ]; then
            : > {output.r1}
            : > {output.r2}
        else
            iss generate \
                --genomes {input} \
                --n_reads {params.nreads} \
                --model {params.model} \
                --cpus {threads} \
                --seed {params.seed} \
                --output iss/{wildcards.sample}/{wildcards.gid}/reads

            rm -f iss/{wildcards.sample}/{wildcards.gid}/reads_abundance.txt
            rm -f iss/{wildcards.sample}/{wildcards.gid}/reads.iss*.vcf
        fi
        """

rule compress:
    input:
        r1 = os.path.join(OUTDIR, "iss/{sample}/{gid}/reads_R1.fastq"),
        r2 = os.path.join(OUTDIR, "iss/{sample}/{gid}/reads_R2.fastq")
    output:
        r1 = temp(os.path.join(OUTDIR, "iss/{sample}/{gid}/reads_R1.fq.gz")),
        r2 = temp(os.path.join(OUTDIR, "iss/{sample}/{gid}/reads_R2.fq.gz"))
    threads: 1
    shell:
        r"""
        set -euo pipefail
        gzip -c "{input.r1}" > "{output.r1}"
        gzip -c "{input.r2}" > "{output.r2}"
        """

rule merge_sample:
    input:
        r1 = lambda w: [os.path.join(OUTDIR, "iss", w.sample, gid, "reads_R1.fq.gz") for gid in IDS],
        r2 = lambda w: [os.path.join(OUTDIR, "iss", w.sample, gid, "reads_R2.fq.gz") for gid in IDS]
    output:
        r1 = os.path.join(OUTDIR, "{sample}_1.fq.gz"),
        r2 = os.path.join(OUTDIR, "{sample}_2.fq.gz")
    wildcard_constraints:
        sample = "|".join(re.escape(s) for s in SAMPLES)
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.r1})"
        cat {input.r1} > {output.r1}
        cat {input.r2} > {output.r2}
        """
