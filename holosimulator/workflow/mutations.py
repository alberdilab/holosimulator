# mutations.py
from __future__ import annotations
import os, time, random, gzip, tempfile
from typing import IO, Optional
from holosimulator.utils import is_url, download_to_temp, ts

# Optional acceleration (exact hypergeometric per chunk)
try:
    import numpy as _np  # noqa: F401
    _HAVE_NUMPY = True
except Exception:
    _np = None
    _HAVE_NUMPY = False

###
# Define text colors
###
HEADER1 = "\033[1;95m"
ERROR   = "\033[1;31m"
INFO    = "\033[1;34m"
RESET   = "\033[0m"
END     = "\033[1;92m"

def _open_maybe_gzip(path: str, mode: str) -> IO:
    return gzip.open(path, mode) if str(path).endswith((".gz", ".bgz", ".bgzip")) else open(path, mode)

def _is_mutable(b: str) -> bool:
    return b in "ACGTacgt"

def _mutate_base(b: str, titv: float, rng: random.Random) -> str:
    b_u = b.upper()
    if b_u not in "ACGT":
        return b
    transitions = {"A":"G","G":"A","C":"T","T":"C"}
    tv = {"A":["C","T"],"G":["C","T"],"C":["A","G"],"T":["A","G"]}
    total = titv + 2.0
    r = rng.random() * total
    alt_u = transitions[b_u] if r < titv else (tv[b_u][0] if rng.random() < 0.5 else tv[b_u][1])
    return alt_u if b.isupper() else alt_u.lower()

def _parse_ani(val: str | float) -> float:
    s = str(val).strip().replace("%",""); x = float(s)
    if x > 1.0: x /= 100.0
    if not (0.0 <= x <= 1.0): raise ValueError("ANI must be in [0,1] or [0–100%]")
    return x

# ---------- progress helpers ----------

class _Progress:
    def __init__(self, total: int, enabled: bool, desc: str = "Mutating"):
        self.total = max(1, total)
        self.enabled = enabled
        self.desc = desc
        self.done = 0
        self._last_print_t = 0.0

    def update(self, n: int = 1):
        if not self.enabled: return
        self.done += n
        now = time.time()
        # Print at most ~2 times per second
        if now - self._last_print_t >= 0.5:
            pct = 100.0 * min(self.done, self.total) / self.total
            barw = 30
            filled = int(barw * pct / 100.0)
            bar = "#" * filled + "-" * (barw - filled)
            print(f"\r{self.desc} [{bar}] {pct:6.2f}% ({self.done:,}/{self.total:,})", end="", flush=True)
            self._last_print_t = now

    def close(self):
        if not self.enabled: return
        pct = 100.0 * min(self.done, self.total) / self.total
        barw = 30
        bar = "#" * barw
        print(f"\r{self.desc} [{bar}] {pct:6.2f}% ({self.done:,}/{self.total:,})", flush=True)

def _count_mutable_total(in_path: str) -> int:
    """Fast single pass to count all A/C/G/T bases (ignores headers)."""
    total = 0
    with _open_maybe_gzip(in_path, "rt") as fh:
        for line in fh:
            if not line or line[0] == ">":  # header
                continue
            s = line.strip()
            for ch in s:
                total += 1 if _is_mutable(ch) else 0
    return total

# ---------- hypergeometric draw per chunk ----------

def _draw_k_in_block(R: int, K: int, m: int, np_rng: "_np.random.Generator|None", py_rng: random.Random) -> int:
    """
    Decide how many SNPs (X) to place in a block of size m, given R remaining mutable and K remaining SNPs.
    - If NumPy is available, use exact Hypergeometric.
    - Else, fallback to Binomial(m, K/R) per-trial (accurate for small m), clamped.
    """
    if m <= 0 or R <= 0 or K <= 0:
        return 0
    if np_rng is not None:
        # numpy hypergeometric: number of 'good' draws among m samples without replacement
        X = int(np_rng.hypergeometric(ngood=K, nbad=R - K, nsample=m))
        if X > K: X = K
        return X
    # Fallback: per-trial Bernoulli (Binomial approx)
    p = K / R
    # quick skip for very small p and large m
    if p < 1e-6 and m >= 1024 and py_rng.random() > (p * m):
        return 0
    succ = 0
    for _ in range(m):
        if py_rng.random() < p:
            succ += 1
    if succ > K: succ = K
    return succ

# ---------- main streaming function with chunking ----------

def mutate_fasta_by_ani_streaming(
    *,
    fasta_in: str,          # local path OR http(s) URL
    fasta_out: str,
    vcf_out: Optional[str],
    target_ani: float | str,
    titv: float = 2.0,
    seed: Optional[int] = None,
    wrap: int = 80,
    progress: bool = True,
    chunk_size: int = 65536,        # NEW: process per-chunk (fast when X==0)
) -> dict:
    rng = random.Random(seed)
    np_rng = _np.random.default_rng(seed) if _HAVE_NUMPY else None
    ani = _parse_ani(target_ani)
    divergence = 1.0 - ani

    # Resolve URL → temp file if needed
    print(f"[{ts()}] Retrieving genome", flush=True)
    tmp_path = None
    in_path = fasta_in
    if is_url(fasta_in):
        tmp_path = download_to_temp(fasta_in)
        in_path = tmp_path

    try:
        # PRE-PASS for progress: count total mutable positions
        print(f"[{ts()}] Calculating mutable positions in genome", flush=True)
        total_mutable_all = _count_mutable_total(in_path) if progress else 0

        if total_mutable_all > 0:
            expected_snps = int(round((1.0 - ani) * total_mutable_all))
            print(f"    {INFO}Mutable positions: {total_mutable_all:,}{RESET}", flush=True)
            print(f"    {INFO}Expected SNPs (target ANI={ani:.4f}): {expected_snps:,}{RESET}", flush=True)

        print(f"[{ts()}] Adding SNP variants to genome ", flush=True)
        prog = _Progress(total_mutable_all, enabled=progress, desc=f"   Mutating (ANI={ani:.4f})")

        out_fa = _open_maybe_gzip(fasta_out, "wt")
        out_vcf = _open_maybe_gzip(vcf_out, "wt") if vcf_out else None
        if out_vcf:
            out_vcf.write("##fileformat=VCFv4.2\n")
            out_vcf.write("##source=holosimulator.workflow.mutations_streaming\n")
            out_vcf.write("##INFO=<ID=ANI,Number=1,Type=Float,Description=\"Target ANI used for mutation\">\n")
            out_vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

        # helper for FASTA wrapping across chunks/contigs
        fa_col = 0
        def _fa_write(s: str):
            nonlocal fa_col
            if wrap <= 0:
                out_fa.write(s)
                return
            i = 0
            n = len(s)
            while i < n:
                take = min(wrap - fa_col, n - i)
                out_fa.write(s[i:i + take])
                fa_col += take
                i += take
                if fa_col >= wrap:
                    out_fa.write("\n")
                    fa_col = 0

        total_mutable = 0
        total_snps = 0

        # --- Stream input, split into per-contig temp files (no RAM blow-up) ---
        print(f"[{ts()}] Parsing contigs", flush=True)
        with _open_maybe_gzip(in_path, "rt") as fh:
            header = None
            tmpf = None
            tmp_contig_path = None
            contig_mutable = 0
            contig_len = 0

            def _finalize_contig():
                nonlocal header, tmpf, tmp_contig_path, contig_mutable, contig_len
                if header is None or tmpf is None:
                    return
                tmpf.close()
                # per-contig K
                K_target = int(round(divergence * contig_mutable))
                # write FASTA header
                out_fa.write(f">{header} | mutated_snp={K_target};target_ani={ani:.6f}\n")
                # mutate from temp file in chunks
                _mutate_contig_from_tmp(
                    header=header,
                    tmp_path=tmp_contig_path,
                    N_mutable=contig_mutable,
                    K_target=K_target,
                    titv=titv,
                    rng=rng,
                    np_rng=np_rng,
                    out_fa_write=_fa_write,
                    out_vcf=out_vcf,
                    ani=ani,
                    wrap=wrap,
                    prog=prog,
                    chunk_size=chunk_size
                )
                # finalize line break if needed
                if wrap > 0 and fa_col != 0:
                    out_fa.write("\n")
                    fa_col = 0
                # update totals
                nonlocal total_mutable, total_snps
                total_mutable += contig_mutable
                total_snps += K_target
                # cleanup
                try: os.remove(tmp_contig_path)
                except OSError: pass
                # reset
                header, tmpf, tmp_contig_path, contig_mutable, contig_len = None, None, None, 0, 0

            for line in fh:
                if line.startswith(">"):
                    _finalize_contig()
                    header = line[1:].strip()
                    fd, tmp_contig_path = tempfile.mkstemp(prefix="holosim_contig_", suffix=".seq")
                    tmpf = os.fdopen(fd, "wt")
                    contig_mutable = 0
                    contig_len = 0
                else:
                    s = line.strip()
                    if not s: continue
                    tmpf.write(s)
                    contig_len += len(s)
                    # count mutable in this chunk
                    for ch in s:
                        if _is_mutable(ch):
                            contig_mutable += 1

            # finalize last contig
            _finalize_contig()

        out_fa.close()
        if out_vcf:
            out_vcf.close()
        prog.close()

        achieved_ani = 1.0 - (total_snps / total_mutable if total_mutable else 0.0)
        return {
            "total_mutable": total_mutable,
            "snp_count": total_snps,
            "achieved_ani": achieved_ani,
            "target_ani": ani,
        }
    finally:
        if tmp_path and os.path.exists(tmp_path):
            try: os.remove(tmp_path)
            except OSError: pass

def _mutate_contig_from_tmp(
    *,
    header: str,
    tmp_path: str,
    N_mutable: int,
    K_target: int,
    titv: float,
    rng: random.Random,
    np_rng: "_np.random.Generator|None",
    out_fa_write,               # callable that handles wrap-aware writing
    out_vcf: Optional[IO],
    ani: float,
    wrap: int,
    prog: Optional[_Progress],
    chunk_size: int,
) -> None:
    """
    Read the contig sequence from a temp file in fixed-size chunks,
    allocate SNPs per-chunk via Hypergeometric (exact with NumPy; otherwise binomial approx),
    mutate only those positions, and stream out FASTA/VCF.
    """
    k_remaining = K_target
    R_remaining = N_mutable
    pos1 = 0  # 1-based within contig (over all bases)

    with open(tmp_path, "rt") as cf:
        while True:
            block = cf.read(chunk_size)
            if not block:
                break
            # indices of mutable bases within this block
            m_idx = [i for i, ch in enumerate(block) if _is_mutable(ch)]
            m = len(m_idx)

            X = _draw_k_in_block(R_remaining, k_remaining, m, np_rng, rng) if m else 0

            if X == 0:
                # fast path: bulk copy
                out_fa_write(block)
                if m and prog: prog.update(m)
            else:
                out_chars = list(block)
                choose = set(rng.sample(m_idx, X))
                for j in choose:
                    ref = out_chars[j]
                    alt = _mutate_base(ref, titv, rng)
                    out_chars[j] = alt
                    if out_vcf:
                        out_vcf.write(f"{header}\t{pos1 + j + 1}\t.\t{ref.upper()}\t{alt.upper()}\t.\tPASS\tANI={ani:.6f}\n")
                out_fa_write("".join(out_chars))
                if prog: prog.update(m)

            if m:
                R_remaining -= m
            pos1 += len(block)

# (legacy helper retained for reference; no longer used)
# def _stream_mutate_and_write(...): pass
