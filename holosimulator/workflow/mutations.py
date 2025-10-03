# mutations_streaming.py
from __future__ import annotations
import os, time, random, gzip
from typing import IO, Optional
from holosimulator.utils import is_url, download_to_temp

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
        # Print at most ~2 times per second, and only if progress advanced
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
            # manual count is faster than sum with function call per char
            # but this is clear and still fast enough
            for ch in s:
                total += 1 if _is_mutable(ch) else 0
    return total

# ---------- main streaming function with progress ----------

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
) -> dict:
    rng = random.Random(seed)
    ani = _parse_ani(target_ani)
    divergence = 1.0 - ani

    # Resolve URL → temp file if needed
    tmp_path = None
    in_path = fasta_in
    if is_url(fasta_in):
        tmp_path = download_to_temp(fasta_in)
        in_path = tmp_path

    try:
        # PRE-PASS for progress: count total mutable positions
        total_mutable_all = _count_mutable_total(in_path) if progress else 0
        prog = _Progress(total_mutable_all, enabled=progress, desc=f"Mutating (ANI={ani:.4f})")

        out_fa = _open_maybe_gzip(fasta_out, "wt")
        out_vcf = _open_maybe_gzip(vcf_out, "wt") if vcf_out else None
        if out_vcf:
            out_vcf.write("##fileformat=VCFv4.2\n")
            out_vcf.write("##source=holosimulator.workflow.mutations_streaming\n")
            out_vcf.write("##INFO=<ID=ANI,Number=1,Type=Float,Description=\"Target ANI used for mutation\">\n")
            out_vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

        total_mutable = 0
        total_snps = 0

        with _open_maybe_gzip(in_path, "rt") as fh:
            header = None
            seq_chunks = []
            for line in fh:
                if line.startswith(">"):
                    if header is not None:
                        seq_str = "".join(seq_chunks)
                        N = sum(1 for ch in seq_str if _is_mutable(ch))
                        k = int(round(divergence * N))
                        _stream_mutate_and_write(header, seq_str, N, k, titv, rng, out_fa, out_vcf, ani, wrap, prog)
                        total_mutable += N
                        total_snps += k
                    header = line[1:].strip()
                    seq_chunks = []
                else:
                    seq_chunks.append(line.strip())

            if header is not None:
                seq_str = "".join(seq_chunks)
                N = sum(1 for ch in seq_str if _is_mutable(ch))
                k = int(round(divergence * N))
                _stream_mutate_and_write(header, seq_str, N, k, titv, rng, out_fa, out_vcf, ani, wrap, prog)
                total_mutable += N
                total_snps += k

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

def _stream_mutate_and_write(
    header: str,
    seq: str,
    N_mutable: int,
    K_target: int,
    titv: float,
    rng: random.Random,
    out_fa: IO,
    out_vcf: Optional[IO],
    ani: float,
    wrap: int,
    prog: Optional[_Progress],
) -> None:
    annotated = f"{header} | mutated_snp={K_target};target_ani={ani:.6f}"
    out_fa.write(f">{annotated}\n")

    k_remaining = K_target
    r_remaining = N_mutable
    pos1 = 0

    out_buf = []
    out_buf_len = 0

    for ch in seq:
        pos1 += 1
        if _is_mutable(ch):
            # sequential selection without replacement
            if k_remaining > 0 and rng.random() < (k_remaining / r_remaining):
                ref_base = ch.upper()
                alt_base = _mutate_base(ch, titv, rng)
                out_buf.append(alt_base)
                k_remaining -= 1
                if out_vcf:
                    out_vcf.write(f"{header}\t{pos1}\t.\t{ref_base}\t{alt_base.upper()}\t.\tPASS\tANI={ani:.6f}\n")
            else:
                out_buf.append(ch)
            r_remaining -= 1
            if prog: prog.update(1)  # advance by one mutable base processed
        else:
            out_buf.append(ch)

        if wrap and len(out_buf) - out_buf_len >= wrap:
            out_fa.write("".join(out_buf[out_buf_len: out_buf_len + wrap]) + "\n")
            out_buf_len += wrap

    if out_buf_len < len(out_buf):
        out_fa.write("".join(out_buf[out_buf_len:]) + "\n")
