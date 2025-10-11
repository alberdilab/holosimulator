configfile: "config.yaml"

import os, json, re, shutil, gzip
from glob import glob

# ------------------
# Config (from YAML)
# ------------------
INPUT_JSON        = config["input"]
OUTDIR            = config["output_dir"]

# -------------------------
# Parse transcriptomics JSON
# -------------------------
with open(INPUT_JSON) as fh:
    META = json.load(fh)

SAMPLES = META["samples"]
TRANSCRIPTOMES = META["transcriptomes"] 

# Maps & lists
TAXA       = [t["organism"] for t in TRANSCRIPTOMES]
IDS        = [t["id"] for t in TRANSCRIPTOMES]
ORG_TO_ID  = {t["organism"]: t["id"] for t in TRANSCRIPTOMES}
ID_TO_ORG  = {t["id"]: t["organism"] for t in TRANSCRIPTOMES}
ID_TO_PATH = {t["id"]: t["path"] for t in TRANSCRIPTOMES}

# Abundances as PAIRS per (organism, sample); we convert to reads later
ABUND = {t["organism"]: {s: int(t["abundances"][s]) for s in SAMPLES} for t in TRANSCRIPTOMES}

# Paths
GENOME_OUT      = os.path.join(OUTDIR, "genomes", "{gid}.fa.gz")
ALL_GENOME_FILES = [GENOME_OUT.format(gid=gid) for gid in IDS]

rule all:
    input:
        INPUT_JSON,
        [os.path.join(OUTDIR, "reads", f"{s}_1.fq.gz") for s in SAMPLES],
        [os.path.join(OUTDIR, "reads", f"{s}_2.fq.gz") for s in SAMPLES]

rule calculate_coverage:
    input:
        os.path.join(OUTDIR, "genomes", "{gid}.fa")
    output:
        os.path.join(OUTDIR, "simulation", "{sample}/{gid}.tsv")
    threads: 1
    params:
        total_reads=10_000_000,
        read_length=150,
        lognorm_mu=-2.0,
        lognorm_sigma=1.2,
        seed=42       
    run:
        import math, numpy as np, random
        print(f"[{ts()}] [Simulate reads] Calculating gene coverages for genome {wildcards.gid} in sample {wildcards.sample}", file=sys.stderr, flush=True)
        # --- read FASTA quickly without external deps ---
        ids = []
        lens = []
        with open(input.fasta, "r") as fh:
            cur_id = None
            cur_len = 0
            for line in fh:
                if line.startswith(">"):
                    # flush previous
                    if cur_id is not None:
                        ids.append(cur_id)
                        lens.append(cur_len)
                    # get new id up to first whitespace
                    cur_id = line[1:].strip().split()[0]
                    cur_len = 0
                else:
                    cur_len += len(line.strip())
            # flush last
            if cur_id is not None:
                ids.append(cur_id)
                lens.append(cur_len)

        if not ids:
            raise ValueError("No sequences found in FASTA.")

        lens = np.asarray(lens, dtype=np.int64)

        # --- sample realistic expression & allocate reads ---
        rng = np.random.default_rng(params.seed)
        # log-normal expression weights (skewed like real transcriptomes)
        weights = rng.lognormal(mean=params.lognorm_mu,
                                sigma=params.lognorm_sigma,
                                size=len(ids)).astype(float)

        # normalize to probabilities; guard against numeric issues
        weights_sum = float(weights.sum())
        if weights_sum == 0.0 or not np.isfinite(weights_sum):
            # fallback to uniform if something odd happens
            probs = np.ones_like(weights) / len(weights)
        else:
            probs = weights / weights_sum

        total_reads = int(params.total_reads)
        if total_reads <= 0:
            raise ValueError("params.total_reads must be > 0")

        # multinomial allocation to ensure read counts sum exactly
        read_counts = rng.multinomial(total_reads, probs)

        # --- compute coverage per transcript ---
        read_len = float(params.read_length)
        lens_safe = lens.astype(float)
        # coverage = (reads * read_length) / transcript_length
        coverage = (read_counts * read_len) / lens_safe

        # --- write TSV (no header): <id>\t<coverage> ---
        with open(output.tsv, "w") as out:
            for tid, cov in zip(ids, coverage):
                # format with reasonable precision
                out.write(f"{tid}\t{cov:.6f}\n")

rule simulate_transcriptome:
    input:
        fa=os.path.join(OUTDIR, "genomes", "{gid}.fa"),
        cov=os.path.join(OUTDIR, "simulation/{sample}/{gid}.tsv")
    output:
        temp(os.path.join(OUTDIR, "simulation/{sample}/{gid}.fastq"))
    params:
        nreads   = lambda w: int(ABUND[ID_TO_ORG[w.gid]][w.sample]),
        art_logdir = lambda w: os.path.join(OUTDIR, "log", "simulate", w.sample, w.gid)
    threads: 1
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.r1})"
        echo "[`date '+%Y-%m-%d %H:%M:%S'`] [Simulate reads] Simulating transcriptomic reads from genome {wildcards.gid} for sample {wildcards.sample}"

        if [ {params.nreads} -eq 0 ]; then
            : > {output}

        else

            art_modern \
            --mode trans \
            --lc pe \
            --i-file {input.fa} \
            --o-fastq {output}
            --builtin_qual_file HiSeq2500_150bp \
            --i-fcov {input.cov} \
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