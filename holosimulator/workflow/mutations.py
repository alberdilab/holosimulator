# holosimulator/workflow/mutations_streaming.py
from __future__ import annotations
import random, gzip
from typing import Tuple, IO, Iterable, Optional
from holosimulator.utils import is_url, download_to_temp

def _open_maybe_gzip(path: str, mode: str) -> IO:
    if str(path).endswith((".gz", ".bgz", ".bgzip")):
        # text mode with gzip, line buffering handled by Python
        return gzip.open(path, mode if "b" in mode else mode.replace("t", "") + "t")
    return open(path, mode)

def _is_mutable(b: str) -> bool:
    return b in "ACGTacgt"

def _mutate_base(b: str, titv: float, rng: random.Random) -> str:
    b_u = b.upper()
    if b_u not in "ACGT":
        return b
    transitions = {"A":"G", "G":"A", "C":"T", "T":"C"}
    tv = {"A":["C","T"], "G":["C","T"], "C":["A","G"], "T":["A","G"]}
    total = titv + 2.0
    r = rng.random() * total
    alt_u = transitions[b_u] if r < titv else (tv[b_u][0] if rng.random() < 0.5 else tv[b_u][1])
    # preserve original case
    return alt_u if b.isupper() else alt_u.lower()

def _parse_ani(val: str | float) -> float:
    s = str(val).strip().replace("%","")
    x = float(s)
    if x > 1.0: x /= 100.0
    if not (0.0 <= x <= 1.0):
        raise ValueError("ANI must be in [0,1] or [0–100%]")
    return x

def mutate_fasta_by_ani_streaming(
    *,
    fasta_in: str,          # local path OR http(s) URL
    fasta_out: str,
    vcf_out: Optional[str],
    target_ani: float | str,
    titv: float = 2.0,
    seed: Optional[int] = None,
    wrap: int = 80,
) -> dict:
    rng = random.Random(seed)
    ani = _parse_ani(target_ani)
    divergence = 1.0 - ani

    # --- NEW: resolve URL → temp file; remember for cleanup
    tmp_path = None
    in_path = fasta_in
    if is_url(fasta_in):
        tmp_path = download_to_temp(fasta_in)  # preserves ".fna.gz" etc
        in_path = tmp_path

    try:
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
                        _stream_mutate_and_write(header, seq_str, N, k, titv, rng, out_fa, out_vcf, ani, wrap)
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
                _stream_mutate_and_write(header, seq_str, N, k, titv, rng, out_fa, out_vcf, ani, wrap)
                total_mutable += N
                total_snps += k

        out_fa.close()
        if out_vcf:
            out_vcf.close()

        achieved_ani = 1.0 - (total_snps / total_mutable if total_mutable else 0.0)
        return {
            "total_mutable": total_mutable,
            "snp_count": total_snps,
            "achieved_ani": achieved_ani,
            "target_ani": ani,
        }
    finally:
        if tmp_path and os.path.exists(tmp_path):
            try:
                os.remove(tmp_path)
            except OSError:
                pass

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
) -> None:
    """
    Given a contig sequence (string), mutate exactly K_target of the N_mutable
    A/C/G/T positions using streaming selection without replacement, and write
    FASTA (wrapped) and VCF rows as we go. Memory: O(1).
    """
    # FASTA header annotated (no need to know per-contig K beforehand—but we do)
    annotated = f"{header} | mutated_snp={K_target};target_ani={ani:.6f}"
    out_fa.write(f">{annotated}\n")

    k_remaining = K_target
    r_remaining = N_mutable
    pos1 = 0  # 1-based genomic position within contig

    # We also need the REF base for VCF; since we have seq, we can take it
    out_buf = []
    out_buf_len = 0

    for ch in seq:
        pos1 += 1
        if _is_mutable(ch):
            # accept with probability k_remaining / r_remaining
            if k_remaining > 0 and rng.random() < (k_remaining / r_remaining):
                ref_base = ch.upper()
                alt_base = _mutate_base(ch, titv, rng)
                out_buf.append(alt_base)
                k_remaining -= 1
                if out_vcf:
                    out_vcf.write(
                        f"{header}\t{pos1}\t.\t{ref_base}\t{alt_base.upper()}\t.\tPASS\tANI={ani:.6f}\n"
                    )
            else:
                out_buf.append(ch)
            r_remaining -= 1
        else:
            out_buf.append(ch)

        # Wrap output lines as we stream
        if wrap and len(out_buf) - out_buf_len >= wrap:
            out_fa.write("".join(out_buf[out_buf_len: out_buf_len + wrap]) + "\n")
            out_buf_len += wrap

    # flush any residual chars
    if out_buf_len < len(out_buf):
        out_fa.write("".join(out_buf[out_buf_len:]) + "\n")
