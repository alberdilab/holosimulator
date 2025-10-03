from __future__ import annotations
import random
from dataclasses import dataclass
from typing import Iterable, List, Tuple, Dict

@dataclass
class MutationResult:
    mutated_records: List[Tuple[str, str]]  # (header, seq)
    snp_count: int
    total_mutable: int
    achieved_ani: float
    chosen_sites: int  # same as snp_count, for clarity

def _parse_fasta(path: str) -> Iterable[Tuple[str, str]]:
    header, chunks = None, []
    with open(path, "r") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header = line[1:].strip()
                chunks = []
            else:
                chunks.append(line.strip())
        if header is not None:
            yield header, "".join(chunks)

def _write_fasta(records: Iterable[Tuple[str, str]], path: str) -> None:
    with open(path, "w") as out:
        for h, s in records:
            out.write(f">{h}\n")
            for i in range(0, len(s), 80):
                out.write(s[i:i+80] + "\n")

def _parse_ani(val: str | float) -> float:
    s = str(val).strip().replace("%", "")
    x = float(s)
    if x > 1.0:
        x /= 100.0
    if not (0.0 <= x <= 1.0):
        raise ValueError("ANI must be in [0,1] or [0,100]")
    return x

def _alt_base(ref: str, titv_ratio: float, rng: random.Random) -> str:
    ref = ref.upper()
    if ref not in "ACGT":
        return ref
    transitions = {"A": "G", "G": "A", "C": "T", "T": "C"}
    transversions = {
        "A": ["C", "T"], "G": ["C", "T"],
        "C": ["A", "G"], "T": ["A", "G"]
    }
    total = titv_ratio + 2.0
    r = rng.random() * total
    if r < titv_ratio:
        return transitions[ref]
    return transversions[ref][0] if rng.random() < 0.5 else transversions[ref][1]

def mutate_fasta_by_ani(
    fasta_in: str,
    target_ani: float | str,
    *,
    titv: float = 2.0,
    seed: int | None = None,
) -> MutationResult:
    """Return mutated sequences hitting the requested ANI via random SNPs."""
    rng = random.Random(seed)
    ani = _parse_ani(target_ani)

    records = list(_parse_fasta(fasta_in))
    if not records:
        raise ValueError("No sequences found in FASTA")

    contigs: List[Tuple[str, List[str]]] = []
    mutable_positions: List[Tuple[int, int]] = []
    for idx, (hdr, seq) in enumerate(records):
        s = seq.upper()
        contigs.append((hdr, list(s)))
        for pos, b in enumerate(s):
            if b in "ACGT":
                mutable_positions.append((idx, pos))

    total = len(mutable_positions)
    if total == 0:
        raise ValueError("No A/C/G/T positions to mutate")

    divergence = 1.0 - ani
    target_snps = round(divergence * total)
    target_snps = max(0, min(target_snps, total))

    chosen_idx = set(rng.sample(range(total), target_snps))
    per_contig: Dict[int, int] = {}

    for i, (cidx, pos) in enumerate(mutable_positions):
        if i not in chosen_idx:
            continue
        ref = contigs[cidx][1][pos]
        alt = _alt_base(ref, titv, rng)
        if alt == ref:
            # very unlikely, but ensure a change
            for base in "ACGT":
                if base != ref:
                    alt = base
                    break
        contigs[cidx][1][pos] = alt
        per_contig[cidx] = per_contig.get(cidx, 0) + 1

    mutated_records: List[Tuple[str, str]] = []
    for cidx, (hdr, seq_list) in enumerate(contigs):
        mutcount = per_contig.get(cidx, 0)
        new_hdr = f"{hdr} | mutated_snp={mutcount};target_ani={ani:.6f}"
        if seed is not None:
            new_hdr += f";seed={seed}"
        mutated_records.append((new_hdr, "".join(seq_list)))

    achieved_ani = 1.0 - (len(chosen_idx) / total) if total else 1.0

    return MutationResult(
        mutated_records=mutated_records,
        snp_count=len(chosen_idx),
        total_mutable=total,
        achieved_ani=achieved_ani,
        chosen_sites=len(chosen_idx),
    )

def write_outputs(
    result: MutationResult,
    *,
    fasta_in: str,
    fasta_out: str,
    vcf_out: str | None,
    target_ani: float | str,
) -> None:
    """Write the mutated FASTA and (optionally) a VCF by re-diffing input vs mutated."""
    # Write FASTA
    _write_fasta(result.mutated_records, fasta_out)

    if not vcf_out:
        return

    # Load originals to recover REF (we only stored the mutated sequences in result)
    originals = {h: s for h, s in _parse_fasta(fasta_in)}

    # Normalize ANI string → float for INFO field
    ani = _parse_ani(target_ani)

    with open(vcf_out, "w") as v:
        v.write("##fileformat=VCFv4.2\n")
        v.write("##source=holosimulator.workflow.mutations\n")
        v.write("##INFO=<ID=ANI,Number=1,Type=Float,Description=\"Target ANI used for mutation\">\n")
        v.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

        for h_mut, s_mut in result.mutated_records:
            # Strip our annotation suffix so we can look up the original header key
            h_key = h_mut.split("|", 1)[0].strip()
            if h_key not in originals:
                # Header mismatch; skip safely
                continue
            s_ref = originals[h_key].upper()
            s_alt = s_mut.upper()

            # SNP-only; lengths should match. If not, diff over the shared span.
            L = min(len(s_ref), len(s_alt))
            for pos1, (r, a) in enumerate(zip(s_ref[:L], s_alt[:L]), start=1):
                if r in "ACGT" and a in "ACGT" and r != a:
                    v.write(f"{h_key}\t{pos1}\t.\t{r}\t{a}\t.\tPASS\tANI={ani:.6f}\n")
