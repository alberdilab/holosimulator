configfile: "config.yaml"

import os, json, re, shutil, gzip
from snakemake.io import expand
from glob import glob

# --- Config ---
INPUT_JSON = config["input"]
OUTDIR     = config["output_dir"]
SEQUENCING_MODEL = config["sequencing_model"]
SEED = config["seed"]

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
        [os.path.join(OUTDIR, "reads", f"{s}_1.fq.gz") for s in SAMPLES],
        [os.path.join(OUTDIR, "reads", f"{s}_2.fq.gz") for s in SAMPLES]

rule simulate:
    input:
        os.path.join(OUTDIR, "genomes", "{gid}.fa")
    output:
        temp(os.path.join(OUTDIR, "simulation/{sample}/{gid}.fastq"))
    log:
        os.path.join(OUTDIR, "log", "art", "{sample}", "{gid}.log")
    params:
        nreads   = lambda w: int(ABUND[ID_TO_ORG[w.gid]][w.sample]),
        art_logdir = lambda w: os.path.join(OUTDIR, "log", "art", w.sample, w.gid)
    threads: 1
    shell:
        r"""
        set -euo pipefail
        mkdir -p "simulation"
        mkdir -p "$(dirname {output})"
        echo "[`date '+%Y-%m-%d %H:%M:%S'`] [Simulate reads] Simulating reads from genome {wildcards.gid} for sample {wildcards.sample}"

        if [ {params.nreads} -eq 0 ]; then
            : > {output}

        else

            # Calculate fold coverage
            GENOME="{input}"
            GENOME_SIZE=$(cat "$GENOME" | awk 'BEGIN{{s=0}} /^>/{{next}} {{s+=length($0)}} END{{printf "%.0f\n", s}}')
            READ_LEN=150
            TOTAL_BASES=$(( {params.nreads} * 2 * READ_LEN ))
            FCOV=$(awk -v b="$TOTAL_BASES" -v g="$GENOME_SIZE" 'BEGIN{{ if (g <= 0) {{ printf "0.00000000"; exit }} printf "%.8f", b/g }}')

            export ART_LOG_DIR="{params.art_logdir}"
            art_modern \
                --mode wgs \
                --lc pe \
                --i-file {input} \
                --o-fastq {output} \
                --parallel {threads} \
                --builtin_qual_file HiSeq2500_150bp \
		        --i-fcov "$FCOV" \
		        --pe_frag_dist_mean 400 \
		        --pe_frag_dist_std_dev 30 \
                --read_len 150

        fi
        """

rule compress:
    input:
        os.path.join(OUTDIR, "simulation/{sample}/{gid}.fastq")
    output:
        r1 = temp(os.path.join(OUTDIR, "simulation/{sample}/{gid}_1.fq.gz")),
        r2 = temp(os.path.join(OUTDIR, "simulation/{sample}/{gid}_2.fq.gz"))
    params:
        nreads   = lambda w: 2 * int(ABUND[ID_TO_ORG[w.gid]][w.sample])
    threads: 1
    shell:
        r"""
        set -euo pipefail

        if [ {params.nreads} -eq 0 ]; then
            : > {output.r1}
            : > {output.r2}

        else

            zcat -f {input} | \
            awk -v R1="{output.r1}" -v R2="{output.r2}" -v T={threads} '
                BEGIN{{ 
                    if (T > 1) {{
                        cmd1 = "pigz -c -p " T " > " R1;
                        cmd2 = "pigz -c -p " T " > " R2;
                    }} else {{
                        cmd1 = "gzip -c > " R1;
                        cmd2 = "gzip -c > " R2;
                    }}
                }}
                {{
                    # read one FASTQ record (4 lines)
                    h=$0; getline s; getline p; getline q;

                    # decide mate by header suffix /1 or /2 (before end or whitespace)
                    mate = 0;
                    if (h ~ /\/1(\s*$)/)      mate = 1;
                    else if (h ~ /\/2(\s*$)/) mate = 2;

                    if (mate==1) {{
                        printf "%s\n%s\n%s\n%s\n", h,s,p,q | cmd1;
                    }} else if (mate==2) {{
                        printf "%s\n%s\n%s\n%s\n", h,s,p,q | cmd2;
                    }}
                    # else: drop non-/1,/2 headers (or add a 3rd output)
                }}
                END{{ close(cmd1); close(cmd2); }}'

        fi
        """

rule merge_sample:
    input:
        r1 = lambda w: expand(os.path.join(OUTDIR, "simulation", "{sample}", "{gid}_1.fq.gz"), sample=w.sample, gid=IDS),
        r2 = lambda w: expand(os.path.join(OUTDIR, "simulation", "{sample}", "{gid}_2.fq.gz"), sample=w.sample, gid=IDS),
    output:
        r1 = os.path.join(OUTDIR, "reads", "{sample}_1.fq.gz"),
        r2 = os.path.join(OUTDIR, "reads", "{sample}_2.fq.gz"),
    threads: 4
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.r1})"
        echo "[`date '+%Y-%m-%d %H:%M:%S'`] [Merge reads] Merging reads for sample {wildcards.sample}"
        cat {input.r1} > {output.r1}
        cat {input.r2} > {output.r2}
        """