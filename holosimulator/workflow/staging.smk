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
GENOME_IN      = os.path.join(OUTDIR, "genomes", "{gid}.fa.gz")
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
    run:
        import urllib.parse, gzip, shutil, os
        from datetime import datetime

        def ts():
            return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        os.makedirs(os.path.dirname(output[0]), exist_ok=True)

        gid  = wildcards.gid
        src  = ID_TO_PATH[gid]

        def is_url(s):
            try:
                from urllib.parse import urlparse
                p = urlparse(s)
                return p.scheme in ("http","https","ftp")
            except Exception:
                return False

        if is_url(src):
            print(f"[{ts()}] [stage_genome] Downloading genome {gid} from {src}", flush=True)
            shell(
                r'''
                set -euo pipefail
                wget \
                  --tries=3 \
                  --waitretry=5 \
                  -O "{output}" "{src}"
                '''
            )
        else:
            print(f"[{ts()}] [stage_genome] Copying genome {gid} from {src}", flush=True)
            shutil.copy(src, output[0])

        # Compress only if needed
        if not output[0].endswith(".gz") and not src.endswith(".gz"):
            print(f"[{ts()}] [stage_genome] Compressing {output[0]}", flush=True)
            tmp_unzipped = output[0] + ".tmp"
            shutil.move(output[0], tmp_unzipped)
            with open(tmp_unzipped, "rb") as fin, gzip.open(output[0], "wb") as fout:
                shutil.copyfileobj(fin, fout)
            os.remove(tmp_unzipped)

        print(f"[{ts()}] [stage_genome] Genome {gid} staged -> {output[0]}", flush=True)

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
