configfile: "config.yaml"

import os, json, re, shutil, gzip
from glob import glob

# ------------------
# Config (from YAML)
# ------------------
INPUT_JSON        = config["input"]                  # augmented JSON (with contigs)
OUTDIR            = config["output_dir"]
SEQUENCING_MODEL  = config["sequencing_model"]
SEED              = int(config["seed"])

# -------------------------
# Parse transcriptomics JSON
# -------------------------
with open(INPUT_JSON) as fh:
    META = json.load(fh)

SAMPLES = META["samples"]
GENOMES = META["genomes"]

TAXA       = [g["organism"] for g in GENOMES]
IDS        = [g["id"] for g in GENOMES]
ORG_TO_ID  = {g["organism"]: g["id"] for g in GENOMES}
ID_TO_ORG  = {g["id"]: g["organism"] for g in GENOMES}
ID_TO_PATH = {g["id"]: g["path"] for g in GENOMES}

CONTIGS_BY_GID = { g["id"]: g.get("contigs", []) for g in GENOMES }
SAMPLE_INDEX   = { s:i for i,s in enumerate(SAMPLES) }

ABUND = {g["organism"]: {s: int(g["abundances"][s]) for s in SAMPLES} for g in GENOMES}

# From transcriptomics JSON:
CONTIGS = { g["id"]: [c["id"] for c in g.get("contigs", [])] for g in GENOMES }
CONTIG_LEN = { g["id"]: {c["id"]: int(c["length"]) for c in g.get("contigs", [])} for g in GENOMES }

# Wildcard constraints so {sample} and {gid} don't swallow slashes
import re
WCS_SAMPLE = "|".join(re.escape(s) for s in SAMPLES)
WCS_GID    = "|".join(re.escape(gid) for gid in IDS)

# --- paths ---
GENOME_IN_FA   = os.path.join(OUTDIR, "genomes", "{gid}.fa")
GENOME_IN_FAGZ = os.path.join(OUTDIR, "genomes", "{gid}.fa.gz")

ALLOC_TSV       = os.path.join(OUTDIR, "iss", "{sample}", "{gid}", "alloc.tsv")
GENOME_MERGED_R1 = os.path.join(OUTDIR, "iss/{sample}/{gid}/reads_R1.fq.gz")
GENOME_MERGED_R2 = os.path.join(OUTDIR, "iss/{sample}/{gid}/reads_R2.fq.gz")
SAMPLE_R1        = os.path.join(OUTDIR, "{sample}_1.fq.gz")
SAMPLE_R2        = os.path.join(OUTDIR, "{sample}_2.fq.gz")

rule all:
    input:
        INPUT_JSON,
        [SAMPLE_R1.format(sample=s) for s in SAMPLES],
        [SAMPLE_R2.format(sample=s) for s in SAMPLES]

# Helper: prefer .fa then .fa.gz, else common alternates

FA_EXTS = (".fa", ".fna", ".fasta")
def _find_local_fasta(gid: str):
    base = os.path.join(OUTDIR, "genomes", gid)
    for ext in (".fa", ".fa.gz", ".fna", ".fna.gz", ".fasta", ".fasta.gz"):
        p = base + (ext if ext.startswith(".") else "." + ext)
        if os.path.exists(p): return p
    # final fallback
    if os.path.exists(base + ".fa"): return base + ".fa"
    return base + ".fa"  # may fail later if truly missing

def _iter_fasta_ids(path):
    import gzip
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                yield line[1:].strip().split()[0]

# ---------- NEW: allocate per-contig read counts (no splitting needed) ----------

rule allocate_gene_reads:
    input:
        fa = lambda w: _find_local_fasta(w.gid)
    output:
        alloc = ALLOC_TSV
    run:
        import os
        os.makedirs(os.path.dirname(output.alloc), exist_ok=True)

        contigs = CONTIGS_BY_GID.get(wildcards.gid, [])
        si = SAMPLE_INDEX[wildcards.sample]
        # JSON map: contig_id -> reads for this sample
        json_counts = { c["id"]: (int(c.get("reads", [0]*len(SAMPLES))[si]) if c.get("reads") else 0)
                        for c in contigs }

        # Write one line **for every FASTA record**, in the exact ID spelling
        with open(output.alloc, "w") as out:
            for rid in _iter_fasta_ids(input.fa):
                n = json_counts.get(rid, 0)  # 0 if JSON didn't have this id
                out.write(f"{rid}\t{n}\n")

# ---------- NEW: single ISS call per genome+sample ----------
# Uses --draft + --readcount_file so ISS distributes reads per contig internally
rule simulate_genome:
    input:
        fa    = lambda w: _find_local_fasta(w.gid),  # multi-FASTA with all contigs
        alloc = ALLOC_TSV
    output:
        r1 = GENOME_MERGED_R1,
        r2 = GENOME_MERGED_R2
    params:
        model = SEQUENCING_MODEL,
        seed  = SEED
    threads: 1
    wildcard_constraints:
        sample = WCS_SAMPLE,
        gid    = WCS_GID
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.r1})"

        # total reads from alloc (no header, 2 cols). If file missing/empty → zero
        total=$(awk '{{s+=$2}} END{{print (s>0)?s:0}}' "{input.alloc}" 2>/dev/null || echo 0)
        if [ "$total" -eq 0 ]; then
            : | gzip -c > "{output.r1}"
            : | gzip -c > "{output.r2}"
            exit 0
        fi

        iss generate \
          --genomes "{input.fa}" \
          --readcount_file "{input.alloc}" \
          --model "{params.model}" \
          --cpus {threads} \
          --seed {params.seed} \
          --compress \
          --output "iss/{wildcards.sample}/{wildcards.gid}/reads"

        mv "iss/{wildcards.sample}/{wildcards.gid}/reads_R1.fastq.gz" "{output.r1}"
        mv "iss/{wildcards.sample}/{wildcards.gid}/reads_R2.fastq.gz" "{output.r2}"

        rm -f "iss/{wildcards.sample}/{wildcards.gid}/reads_abundance.txt" \
              "iss/{wildcards.sample}/{wildcards.gid}/reads.iss"*.vcf || true
        """

# ---------- final sample merge (unchanged) ----------
rule merge_sample:
    input:
        r1 = lambda w: [GENOME_MERGED_R1.format(sample=w.sample, gid=gid) for gid in IDS],
        r2 = lambda w: [GENOME_MERGED_R2.format(sample=w.sample, gid=gid) for gid in IDS]
    output:
        r1 = SAMPLE_R1,
        r2 = SAMPLE_R2
    wildcard_constraints:
        sample = WCS_SAMPLE
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.r1})"
        cat {input.r1} > "{output.r1}"
        cat {input.r2} > "{output.r2}"
        """