from __future__ import annotations
import os
import io
import sys
import yaml
import re
import json
import math
import gzip
import tempfile
import requests
import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime
from collections import defaultdict
from urllib.parse import urlparse
from typing import Iterable, List, Dict, Any, Optional

def ts():
            return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

# Check if snakemake directory is locked

def is_snakemake_locked(workdir: str) -> bool:
    locks_dir = os.path.join(workdir, ".snakemake", "locks")
    return os.path.isdir(locks_dir) and len(os.listdir(locks_dir)) > 0

# CSV inputs

def csv_to_json_genomes(input_csv: str, output_json: str) -> None:
    os.makedirs(os.path.dirname(output_json), exist_ok=True)
    df = pd.read_csv(input_csv, sep=None, engine="python")
    df.columns = df.columns.str.strip().str.replace(r"\ufeff", "", regex=True)
    if not {"Organism", "Path"}.issubset(df.columns):
        missing = {"Organism", "Path"} - set(df.columns)
        raise ValueError(f"Input file missing required columns: {missing}")
    df["Organism"] = df["Organism"].astype(str)

    sample_cols = [c for c in df.columns if str(c).startswith("Sample")]
    if not sample_cols:
        raise ValueError("No sample columns found (expected columns starting with 'Sample').")

    # deterministic ID assignment in row order
    genomes = []
    for i, r in enumerate(df.itertuples(index=False), start=1):
        gid = f"G{i:04d}"
        org = getattr(r, "Organism")
        path = getattr(r, "Path")
        cat = getattr(r, "Category", None) if "Category" in df.columns else None
        abund = {s: float(getattr(r, s)) for s in sample_cols}
        genomes.append({
            "id": gid,
            "organism": org,
            "category": None if pd.isna(cat) else str(cat),
            "path": str(path),
            "abundances": abund,
        })

    data = {"samples": sample_cols, "genomes": genomes}
    with open(output_json, "w") as fh:
        json.dump(data, fh, indent=2)

def csv_to_json_transcriptomes(input_csv: str, output_json: str) -> None:
    os.makedirs(os.path.dirname(output_json), exist_ok=True)
    df = pd.read_csv(input_csv, sep=None, engine="python")
    df.columns = df.columns.str.strip().str.replace(r"\ufeff", "", regex=True)
    if not {"Organism", "Path"}.issubset(df.columns):
        missing = {"Organism", "Path"} - set(df.columns)
        raise ValueError(f"Input file missing required columns: {missing}")
    df["Organism"] = df["Organism"].astype(str)

    sample_cols = [c for c in df.columns if str(c).startswith("Sample")]
    if not sample_cols:
        raise ValueError("No sample columns found (expected columns starting with 'Sample').")

    # deterministic ID assignment in row order
    transcriptomes = []
    for i, r in enumerate(df.itertuples(index=False), start=1):
        gid = f"G{i:04d}"
        org = getattr(r, "Organism")
        path = getattr(r, "Path")
        cat = getattr(r, "Category", None) if "Category" in df.columns else None
        abund = {s: float(getattr(r, s)) for s in sample_cols}
        transcriptomes.append({
            "id": gid,
            "organism": org,
            "category": None if pd.isna(cat) else str(cat),
            "path": str(path),
            "abundances": abund,
        })

    data = {"samples": sample_cols, "transcriptomes": transcriptomes}
    with open(output_json, "w") as fh:
        json.dump(data, fh, indent=2)

# Argument input

def args_to_genomics_json(
    host, microbiome, sample_size, sequencing_depth, sequencing_depth_variance,
    host_fraction, host_fraction_variance, microbiome_variance, seed,
    output_json):

    # ---------- small helpers (scoped here so this function is self-contained) ----------
    def _parse_csv_list(s):
        if not s:
            return []
        return [x.strip() for x in str(s).split(",") if x.strip()]

    def as_fraction(x: float) -> float:
        """Accept 0..1 or 0..100 (percent) and return fraction 0..1."""
        x = float(x)
        return x / 100.0 if x > 1.0 and x <= 100.0 else x

    def safe_normal(mean, sd):
        if sd <= 0:
            return max(0.0, float(mean))
        val = np.random.normal(loc=mean, scale=sd)
        return max(0.0, float(val))

    def dirichlet_with_variance(k: int, var_pct: float):
        """
        Map a user variance % into a Dirichlet concentration.
        v=0 -> low variance (alpha~10), v=100 -> high variance (alpha~0.5).
        """
        var_pct = min(max(float(var_pct), 0.0), 100.0)
        alpha_scale = 10.0 - 9.5 * (var_pct / 100.0)
        alpha = np.full(k, alpha_scale, dtype=float)
        return np.random.dirichlet(alpha)

    def basename_no_ext(p: str) -> str:
        name = os.path.basename(p)
        for ext in [".fa.gz", ".fna.gz", ".fasta.gz", ".fa", ".fna", ".fasta", ".gz"]:
            if name.endswith(ext):
                return name[: -len(ext)]
        return name

    def _to_float_pairs(genomes, samples):
        """Serialize abundances as floats (e.g., 20300.0) to match your example JSON."""
        for g in genomes:
            for s in samples:
                g["abundances"][s] = float(g["abundances"][s])
        return genomes

    def _csv_to_inputs_json(input_csv: str, output_json_path: Path):
        import pandas as pd  # only needed for CSV fallback
        output_json_path.parent.mkdir(parents=True, exist_ok=True)
        df = pd.read_csv(input_csv, sep=None, engine="python")
        df.columns = df.columns.str.strip().str.replace(r"\ufeff", "", regex=True)
        if not {"Organism", "Path"}.issubset(df.columns):
            missing = {"Organism", "Path"} - set(df.columns)
            raise ValueError(f"Input file missing required columns: {missing}")
        df["Organism"] = df["Organism"].astype(str)
        sample_cols = [c for c in df.columns if str(c).startswith("Sample")]
        if not sample_cols:
            raise ValueError("No sample columns found (expected 'Sample*' columns).")

        genomes = []
        for i, r in enumerate(df.itertuples(index=False), start=1):
            gid = f"G{i:04d}"
            org = getattr(r, "Organism")
            path = getattr(r, "Path")
            cat  = getattr(r, "Category", None) if "Category" in df.columns else None
            abund = {s: float(getattr(r, s)) for s in sample_cols}
            genomes.append({
                "id": gid,
                "organism": org,
                "category": None if pd.isna(cat) else str(cat),
                "path": str(path),
                "abundances": abund,
            })
        data = {"samples": sample_cols, "genomes": genomes}
        output_json_path.write_text(json.dumps(data, indent=2))
        return output_json_path

def args_to_transcriptomics_json(
    host, microbiome, sample_size, sequencing_depth, sequencing_depth_variance,
    host_fraction, host_fraction_variance, microbiome_variance, seed,
    output_json):

    # ---------- small helpers (scoped here so this function is self-contained) ----------
    def _parse_csv_list(s):
        if not s:
            return []
        return [x.strip() for x in str(s).split(",") if x.strip()]

    def as_fraction(x: float) -> float:
        """Accept 0..1 or 0..100 (percent) and return fraction 0..1."""
        x = float(x)
        return x / 100.0 if x > 1.0 and x <= 100.0 else x

    def safe_normal(mean, sd):
        if sd <= 0:
            return max(0.0, float(mean))
        val = np.random.normal(loc=mean, scale=sd)
        return max(0.0, float(val))

    def dirichlet_with_variance(k: int, var_pct: float):
        """
        Map a user variance % into a Dirichlet concentration.
        v=0 -> low variance (alpha~10), v=100 -> high variance (alpha~0.5).
        """
        var_pct = min(max(float(var_pct), 0.0), 100.0)
        alpha_scale = 10.0 - 9.5 * (var_pct / 100.0)
        alpha = np.full(k, alpha_scale, dtype=float)
        return np.random.dirichlet(alpha)

    def basename_no_ext(p: str) -> str:
        name = os.path.basename(p)
        for ext in [".fa.gz", ".fna.gz", ".fasta.gz", ".fa", ".fna", ".fasta", ".gz"]:
            if name.endswith(ext):
                return name[: -len(ext)]
        return name

    def _to_float_pairs(genomes, samples):
        """Serialize abundances as floats (e.g., 20300.0) to match your example JSON."""
        for g in genomes:
            for s in samples:
                g["abundances"][s] = float(g["abundances"][s])
        return genomes

    def _csv_to_inputs_json(input_csv: str, output_json_path: Path):
        import pandas as pd  # only needed for CSV fallback
        output_json_path.parent.mkdir(parents=True, exist_ok=True)
        df = pd.read_csv(input_csv, sep=None, engine="python")
        df.columns = df.columns.str.strip().str.replace(r"\ufeff", "", regex=True)
        if not {"Organism", "Path"}.issubset(df.columns):
            missing = {"Organism", "Path"} - set(df.columns)
            raise ValueError(f"Input file missing required columns: {missing}")
        df["Organism"] = df["Organism"].astype(str)
        sample_cols = [c for c in df.columns if str(c).startswith("Sample")]
        if not sample_cols:
            raise ValueError("No sample columns found (expected 'Sample*' columns).")

        transcriptomes = []
        for i, r in enumerate(df.itertuples(index=False), start=1):
            gid = f"G{i:04d}"
            org = getattr(r, "Organism")
            path = getattr(r, "Path")
            cat  = getattr(r, "Category", None) if "Category" in df.columns else None
            abund = {s: float(getattr(r, s)) for s in sample_cols}
            transcriptomes.append({
                "id": gid,
                "organism": org,
                "category": None if pd.isna(cat) else str(cat),
                "path": str(path),
                "abundances": abund,
            })
        data = {"samples": sample_cols, "transcriptomes": transcriptomes}
        output_json_path.write_text(json.dumps(data, indent=2))
        return output_json_path
    # ---------------------------------------------------------------------

    # Prepare output path
    out = Path(output_json)
    out.parent.mkdir(parents=True, exist_ok=True)

    # Normalize inputs
    host_paths  = _parse_csv_list(host)
    micro_paths = _parse_csv_list(microbiome)

    # LIST MODE (preferred)
    if host_paths or micro_paths:
        if seed is not None:
            np.random.seed(int(seed))

        # numeric params / defaults
        n_samples          = int(sample_size) if sample_size is not None else 1
        depth_mean_pairs   = float(sequencing_depth) if sequencing_depth is not None else 3000000
        depth_var_pct      = float(sequencing_depth_variance) if sequencing_depth_variance is not None else 2
        host_frac_mean     = as_fraction(host_fraction if host_fraction is not None else 0.2)
        host_frac_var_pct  = float(host_fraction_variance) if host_fraction_variance is not None else 2
        micro_var_pct      = float(microbiome_variance) if microbiome_variance is not None else 2

        samples = [f"Sample{i}" for i in range(1, n_samples + 1)]
        depth_sd = (depth_var_pct / 100.0) * depth_mean_pairs
        host_frac_sd = as_fraction(host_frac_var_pct) * host_frac_mean  # proportional SD

        genomes = []
        gid_counter = 1

        # Register hosts
        for p in host_paths:
            gid = f"G{gid_counter:04d}"; gid_counter += 1
            genomes.append({
                "id": gid,
                "organism": basename_no_ext(p) or f"host_{gid}",
                "category": "Host",
                "path": p,
                "abundances": {s: 0 for s in samples}
            })

        # Register microbes
        for p in micro_paths:
            gid = f"G{gid_counter:04d}"; gid_counter += 1
            genomes.append({
                "id": gid,
                "organism": basename_no_ext(p) or f"microbe_{gid}",
                "category": "Microbiome",
                "path": p,
                "abundances": {s: 0 for s in samples}
            })

        hosts    = [g for g in genomes if g["category"] == "Host"]
        microbes = [g for g in genomes if g["category"] == "Microbiome"]
        k_host   = max(1, len(hosts))
        k_micro  = max(1, len(microbes))

        # Fill per-sample abundances (pairs)
        for s in samples:
            total_pairs = safe_normal(depth_mean_pairs, depth_sd)
            host_frac   = min(max(safe_normal(host_frac_mean, host_frac_sd), 0.0), 1.0)
            host_pairs  = total_pairs * host_frac
            micro_pairs = total_pairs - host_pairs

            if k_host > 0:
                per_host = host_pairs / k_host
                for g in hosts:
                    g["abundances"][s] = int(round(per_host))
            if k_micro > 0:
                weights = dirichlet_with_variance(k_micro, micro_var_pct)
                for g, w in zip(microbes, weights):
                    g["abundances"][s] = int(round(micro_pairs * w))

            # rounding drift correction
            current_sum = sum(g["abundances"][s] for g in genomes)
            drift = int(round(total_pairs)) - current_sum
            if drift != 0:
                pool = microbes if k_micro > 0 else hosts
                if pool:
                    j = max(range(len(pool)), key=lambda idx: pool[idx]["abundances"][s])
                    pool[j]["abundances"][s] = max(0, pool[j]["abundances"][s] + drift)

        data = {
            "samples": samples,
            "genomes": _to_float_pairs(genomes, samples),  # float serialize
            "by_id": {g["id"]: {"organism": g["organism"], "path": g["path"], "category": g["category"]} for g in genomes},
            "by_organism": {g["organism"]: g["id"] for g in genomes}
        }

        out.write_text(json.dumps(data, indent=2))
        return out

# Check existence of genomes
def check_genomics_paths(in_json: str | Path, timeout: int = 10) -> dict[str, bool]:

    in_path = Path(in_json)
    cfg = json.loads(in_path.read_text())
    results: dict[str, bool] = {}

    for g in cfg.get("genomes", []):
        p = g.get("path", "")
        if not p:
            results[p] = False
            continue
        parsed = urlparse(p)
        if parsed.scheme in ("http", "https", "ftp"):
            try:
                r = requests.head(p, allow_redirects=True, timeout=timeout)
                if r.status_code >= 400:
                    # Some servers don't support HEAD; try GET with stream
                    r = requests.get(p, stream=True, timeout=timeout)
                results[p] = (r.status_code < 400)
            except Exception:
                results[p] = False
        else:
            results[p] = os.path.exists(p)
    return results

# Genomics to transcriptomics simulation

FA_EXTS: tuple[str, ...] = (".fa", ".fna", ".fasta")


def _open_text(fp: Path):
    """Open a local file (optionally .gz) as text (utf-8)."""
    p = str(fp)
    if p.endswith(".gz"):
        return io.TextIOWrapper(gzip.GzipFile(filename=p, mode="rb"))
    return open(p, "rt", encoding="utf-8")


def _candidate_paths(genomes_dir: Path, gid: str, source_path: str) -> Iterable[Path]:
    """
    Yield plausible local filenames for a genome FASTA, in priority order:
      1) genomes/{gid}.{fa|fna|fasta}[.gz]
      2) genomes/<basename from original path or URL> (as staged)
    """
    # 1) Canonical names using gid
    for ext in FA_EXTS:
        yield genomes_dir / f"{gid}{ext}"
        yield genomes_dir / f"{gid}{ext}.gz"

    # 2) Original filename (if staging kept it)
    parsed = urlparse(source_path)
    base = os.path.basename(parsed.path) if parsed.scheme else os.path.basename(source_path)
    if base:
        p = genomes_dir / base
        yield p
        if not base.endswith(".gz"):
            yield genomes_dir / (base + ".gz")


def _find_local_fasta(genomes_dir: Path, gid: str, source_path: str) -> Optional[Path]:
    """Return first existing candidate path, or None if not found."""
    for cand in _candidate_paths(genomes_dir, gid, source_path):
        if cand.exists():
            return cand
    return None


def _scan_contigs(local_fasta: Path, max_contigs: Optional[int] = None) -> List[Dict[str, Any]]:
    """
    Stream a FASTA and return list of dicts:
      [{"id": <contig_id>, "length": <bp>, "gc": <0..1>}, ...]
    - Contig id = first token after '>' (header up to first whitespace).
    - Length is the sum of non-whitespace characters across sequence lines.
    - GC is fraction across A/C/G/T (Ns ignored); falls back to 0.5 if unknown.
    """
    out: List[Dict[str, Any]] = []
    with _open_text(local_fasta) as fh:
        hdr: Optional[str] = None
        L = 0
        gc_n = 0
        atgc_n = 0

        def _flush():
            nonlocal hdr, L, gc_n, atgc_n
            if hdr is None:
                return
            gc_frac = (gc_n / atgc_n) if atgc_n > 0 else 0.5
            out.append({"id": hdr, "length": L, "gc": gc_frac})

        for raw in fh:
            if raw.startswith(">"):
                _flush()
                hdr = raw[1:].strip().split()[0]
                L = 0
                gc_n = 0
                atgc_n = 0
            else:
                s = raw.strip().upper()
                L += len(s)
                # quick GC/AT count without allocating another string
                for ch in s:
                    if ch in ("A", "C", "G", "T"):
                        atgc_n += 1
                        if ch in ("G", "C"):
                            gc_n += 1
        _flush()

    if max_contigs is not None:
        out = out[:max_contigs]
    return out


def iter_fasta_ids(path: Path):
    """Yield exact record IDs from a FASTA (first token after '>')."""
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                yield line[1:].strip().split()[0]


def _eligible_mask(contig_ids, contig_len, read_length: int, require_effective: bool = True):
    """
    Return mask of contigs that can be simulated by ISS for the given read_length.
    require_effective=True filters contigs with effective length <= 0 (L - R + 1 <= 0).
    Otherwise we filter by raw length (L < R).
    """
    mask = []
    for c in contig_ids:
        L = max(0, int(contig_len.get(c, 0)))
        if require_effective:
            ok = (L - read_length + 1) > 0
        else:
            ok = L >= read_length
        mask.append(ok)
    return mask

# ------------------------------------------------------------
# Allocation strategies (registry-based)
# ------------------------------------------------------------

def _round_with_drift(target_total: int, weights: List[float]) -> List[int]:
    """Scale weights to the exact target_total via rounding + drift correction."""
    if target_total <= 0 or not weights:
        return [0] * len(weights)
    s = sum(weights)
    if s <= 0:
        raw = [0] * len(weights)
    else:
        raw = [round(target_total * (w / s)) for w in weights]
    drift = target_total - sum(raw)
    i = 0
    n = len(raw)
    while drift != 0 and n > 0:
        raw[i] += 1 if drift > 0 else -1
        drift += -1 if drift > 0 else 1
        i = (i + 1) % n
    return [int(max(0, r)) for r in raw]


def _apply_min_reads(counts: List[int], min_reads: int) -> List[int]:
    if min_reads <= 0:
        return counts
    return [r if r == 0 or r >= min_reads else min_reads for r in counts]


def alloc_uniform(*, contig_ids, contig_len, total_reads, rng, min_reads=0,
                  read_length: int = 150, drop_short: bool = True, require_effective: bool = True, **kwargs):
    elig = [True]*len(contig_ids) if not drop_short else _eligible_mask(contig_ids, contig_len, read_length, require_effective)
    weights = [1.0 if e else 0.0 for e in elig]
    counts = _round_with_drift(total_reads, weights)
    counts = _apply_min_reads(counts, min_reads)
    # zero out ineligible explicitly (in case min_reads bumped them)
    counts = [0 if not e else c for e, c in zip(elig, counts)]
    return dict(zip(contig_ids, counts))

def alloc_length(*, contig_ids, contig_len, total_reads, rng, min_reads=0,
                 read_length: int = 150, drop_short: bool = True, require_effective: bool = True, **kwargs):
    elig = [True]*len(contig_ids) if not drop_short else _eligible_mask(contig_ids, contig_len, read_length, require_effective)
    weights = [max(1, int(contig_len[c])) if e else 0.0 for c, e in zip(contig_ids, elig)]
    counts = _round_with_drift(total_reads, weights)
    counts = _apply_min_reads(counts, min_reads)
    counts = [0 if not e else c for e, c in zip(elig, counts)]
    return dict(zip(contig_ids, counts))

def alloc_zipf(*, contig_ids, contig_len, total_reads, rng, s: float = 1.1, min_reads=0,
               read_length: int = 150, drop_short: bool = True, require_effective: bool = True, **kwargs):
    elig = {c: True for c in contig_ids} if not drop_short else {
        c: ok for c, ok in zip(contig_ids, _eligible_mask(contig_ids, contig_len, read_length, require_effective))
    }
    ranked = [c for c in sorted(contig_ids, key=lambda x: contig_len.get(x, 1), reverse=True) if elig[c]]
    z = [1.0 / ((i+1) ** max(0.01, float(s))) for i in range(len(ranked))]
    counts_ranked = _round_with_drift(total_reads, z)
    counts_ranked = _apply_min_reads(counts_ranked, min_reads)
    out = {c: 0 for c in contig_ids}
    for c, k in zip(ranked, counts_ranked): out[c] = k
    return out

def alloc_tpm(*, contig_ids, contig_len, total_reads, rng,
              read_length: int = 150, ln_mu: float = 0.0, ln_sigma: float = 1.0,
              gc=None, gc_bias_strength: float = 0.0, gc_bias_linear: float = 0.0,
              overdispersion: float = 50.0, dropout_rate: float = 0.0, min_reads: int = 0,
              drop_short: bool = True, require_effective: bool = True, **kwargs):
    n = len(contig_ids)
    if n == 0 or total_reads <= 0:
        return {c: 0 for c in contig_ids}

    elig = [True]*n if not drop_short else _eligible_mask(contig_ids, contig_len, read_length, require_effective)
    # latent expression
    expr = rng.lognormal(mean=ln_mu, sigma=ln_sigma, size=n).tolist()
    # GC bias
    if gc_bias_strength != 0.0 or gc_bias_linear != 0.0:
        adj = []
        for c, w in zip(contig_ids, expr):
            g = (gc or {}).get(c, 0.5); d = g - 0.5
            adj.append(w * math.exp(gc_bias_strength*(d*d) + gc_bias_linear*d))
        expr = adj
    # effective length
    effL = [max(1, int(contig_len[c]) - read_length + 1) for c in contig_ids]
    weights = [(w / el) if e else 0.0 for w, el, e in zip(expr, effL, elig)]

    s = sum(weights)
    probs = [w / s for w in weights] if s > 0 else [0.0]*n
    # Dirichlet-like gamma draw
    alpha = max(1e-3, float(overdispersion))
    gam = rng.gamma(shape=alpha * np.maximum(probs, 1e-16), scale=1.0)
    sg = float(np.sum(gam)) or 1.0
    theta = (gam / sg).tolist()

    counts = _round_with_drift(total_reads, theta)

    # dropouts & min_reads, and zero ineligible
    if dropout_rate > 0.0:
        mask = rng.random(n) < dropout_rate
        counts = [0 if m else c for m, c in zip(mask, counts)]
    counts = _apply_min_reads(counts, min_reads)
    counts = [0 if not e else c for e, c in zip(elig, counts)]
    return dict(zip(contig_ids, counts))


# Strategy registry
ALLOC_STRATEGIES: Dict[str, callable] = {
    "uniform": alloc_uniform,
    "length":  alloc_length,
    "zipf":    alloc_zipf,
    "tpm":     alloc_tpm,
}

# ------------------------------------------------------------
# Path validation helper (optional preflight)
# ------------------------------------------------------------

def check_genomics_paths(in_json: str | Path, timeout: int = 10) -> Dict[str, bool]:
    """
    Validate that all genome 'path' entries in a genomics JSON exist.

    - If 'path' is local -> os.path.exists()
    - If 'path' is http(s)/ftp URL -> try HEAD (falls back to GET if HEAD not allowed)

    Returns: Dict {path: True/False}
    """
    in_path = Path(in_json)
    cfg = json.loads(in_path.read_text())
    results: Dict[str, bool] = {}

    for g in cfg.get("genomes", []):
        p = g.get("path", "")
        if not p:
            results[p] = False
            continue
        parsed = urlparse(p)
        if parsed.scheme in ("http", "https", "ftp"):
            try:
                r = requests.head(p, allow_redirects=True, timeout=timeout)
                if r.status_code >= 400:
                    r = requests.get(p, stream=True, timeout=timeout)
                results[p] = (r.status_code < 400)
            except Exception:
                results[p] = False
        else:
            results[p] = os.path.exists(p)
    return results

# ------------------------------------------------------------
# Mutations utils
# ------------------------------------------------------------

def is_url(s: str) -> bool:
    try:
        u = urlparse(s)
        return u.scheme in ("http", "https")
    except Exception:
        return False

def download_to_temp(url: str, *, suffix: str | None = None, chunk_bytes: int = 1 << 20, timeout: int = 60) -> str:
    """
    Stream a remote file to a local temporary file and return its path.
    Suffix helps keep .gz extension so gzip-aware openers behave.
    """
    if suffix is None:
        # try to inherit suffix from URL path (.fna.gz, .fa.gz, etc)
        suffix = os.path.splitext(urlparse(url).path)[1]  # ".gz"
        # if .gz, try to keep the pre-suffix too
        if suffix == ".gz":
            base = os.path.basename(urlparse(url).path)  # genome.fna.gz
            if base.count(".") >= 2:
                suffix = "." + ".".join(base.split(".")[-2:])  # ".fna.gz"

    with requests.get(url, stream=True, timeout=timeout) as r:
        r.raise_for_status()
        fd, tmp = tempfile.mkstemp(prefix="holosim_dl_", suffix=suffix or "")
        with os.fdopen(fd, "wb") as f:
            for chunk in r.iter_content(chunk_size=chunk_bytes):
                if chunk:
                    f.write(chunk)
    return tmp
