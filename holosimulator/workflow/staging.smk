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
    run:
        import urllib.parse, tempfile, shutil, gzip, os
        from datetime import datetime
        def ts():
            return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        os.makedirs(os.path.dirname(output[0]), exist_ok=True)

        gid  = wildcards.gid
        src  = ID_TO_PATH[gid]

        def is_url(s):
            try:
                p = urllib.parse.urlparse(s)
                return p.scheme in ("http","https","ftp")
            except Exception:
                return False

        tf = tempfile.NamedTemporaryFile(delete=False); tf.close()
        
        if is_url(src):
            print(f"[{ts()}] [Download genome] Downloading genome {gid} from {src}", flush=True)
            shell(f'wget -q -O "{tf.name}" "{src}"')
            print(f"[{ts()}] [Stage genome] Genome {gid} staged", flush=True)
        else:
            print(f"[{ts()}] [Copy genome] Copying genome {gid} from {src}", flush=True)
            shutil.copy(src, tf.name)
            print(f"[{ts()}] [Stage genome] Genome {gid} staged", flush=True)

        if src.endswith(".gz") or tf.name.endswith(".gz"):
            shutil.move(tf.name, output[0])
        else:
            with open(tf.name, "rb") as fin, gzip.open(output[0], "wb") as fout:
                shutil.copyfileobj(fin, fout)
            os.remove(tf.name)

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
