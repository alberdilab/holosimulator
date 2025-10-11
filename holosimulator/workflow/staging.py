#!/usr/bin/env python3
"""Simple sequential downloader: JSON -> outdir/genomes/{gid}.fa (decompressed)."""

from __future__ import annotations
import json
import os
import shutil
import sys
import time
import gzip
from pathlib import Path
from urllib.parse import urlparse
from urllib.request import Request, urlopen

def _ts() -> str:
    return time.strftime("%Y-%m-%d %H:%M:%S")

def _log(msg: str) -> None:
    print(f"[{_ts()}] {msg}", file=sys.stderr, flush=True)

def _is_url(s: str) -> bool:
    try:
        return urlparse(s).scheme in ("http", "https", "ftp")
    except Exception:
        return False

def _ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)

def _stream_copy(fin, fout, chunk: int = 1024 * 1024) -> None:
    shutil.copyfileobj(fin, fout, length=chunk)

def _download_to_fa(src: str, fa_path: Path, *, retries: int = 3, timeout: float = 60.0) -> None:
    """Download src (URL or local file). If it looks gzipped, decompress on the fly -> fa_path.part, then atomic rename."""
    part = fa_path.with_suffix(fa_path.suffix + ".part")

    for attempt in range(1, retries + 1):
        try:
            if _is_url(src):
                req = Request(src, headers={"User-Agent": "holosim-minidownloader/1.0"})
                with urlopen(req, timeout=timeout) as resp:
                    # Decide whether to gunzip while streaming
                    ct = getattr(resp, "headers", {}).get("Content-Type", "") or ""
                    ce = getattr(resp, "headers", {}).get("Content-Encoding", "") or ""
                    looks_gz = src.endswith(".gz") or "gzip" in ct or "gzip" in ce

                    with open(part, "wb") as fout:
                        if looks_gz:
                            with gzip.GzipFile(fileobj=resp) as fin:
                                _stream_copy(fin, fout)
                        else:
                            _stream_copy(resp, fout)
            else:
                # Local file path
                src_path = Path(src).expanduser()
                if not src_path.exists():
                    raise FileNotFoundError(src)
                with open(part, "wb") as fout:
                    if src_path.suffix == ".gz":
                        with gzip.open(src_path, "rb") as fin:
                            _stream_copy(fin, fout)
                    else:
                        with open(src_path, "rb") as fin:
                            _stream_copy(fin, fout)

            if part.stat().st_size == 0:
                raise IOError("wrote 0 bytes")
            os.replace(part, fa_path)  # atomic publish in same dir
            return
        except Exception as e:
            if part.exists():
                try:
                    part.unlink()
                except Exception:
                    pass
            if attempt < retries:
                _log(f"[retry {attempt}/{retries}] {src}: {e}")
                time.sleep(2 * attempt)
                continue
            raise

def staging_genomes(json_path: str | Path, outdir: str | Path, *, overwrite: bool = False) -> list[Path]:
    """
    Read genomes from JSON and write decompressed FASTA files to <outdir>/genomes/{gid}.fa
    Returns the list of produced/kept .fa paths. Sequential; no parallelism.
    """
    json_path = Path(json_path)
    outdir = Path(outdir)
    genomes_dir = outdir / "genomes"
    _ensure_dir(genomes_dir)

    with open(json_path) as fh:
        meta = json.load(fh)

    genomes = meta.get("genomes", [])
    if not genomes:
        _log("No genomes found in JSON (.genomes[]).")
        return []

    produced: list[Path] = []
    for g in genomes:
        gid = g["id"]
        src = g["path"]
        fa = genomes_dir / f"{gid}.fa"

        if fa.exists() and fa.stat().st_size > 0 and not overwrite:
            _log(f"[Skip download] {gid} already present -> {fa}")
            produced.append(fa)
            continue

        _log(f"[Download genome] {gid} from {src} -> {fa}")
        _ensure_dir(fa.parent)
        _download_to_fa(src, fa)
        produced.append(fa)

    return produced

def staging_transcriptomes(json_path: str | Path, outdir: str | Path, *, overwrite: bool = False) -> list[Path]:
    """
    Read genomes from JSON and write decompressed FASTA files to <outdir>/genomes/{gid}.fa
    Returns the list of produced/kept .fa paths. Sequential; no parallelism.
    """
    json_path = Path(json_path)
    outdir = Path(outdir)
    genomes_dir = outdir / "transcriptomes"
    _ensure_dir(genomes_dir)

    with open(json_path) as fh:
        meta = json.load(fh)

    genomes = meta.get("transcriptomes", [])
    if not genomes:
        _log("No transcriptomes found in JSON (.genomes[]).")
        return []

    produced: list[Path] = []
    for g in genomes:
        gid = g["id"]
        src = g["path"]
        fa = genomes_dir / f"{gid}.fa"

        if fa.exists() and fa.stat().st_size > 0 and not overwrite:
            _log(f"[Skip download] {gid} already present -> {fa}")
            produced.append(fa)
            continue

        _log(f"[Download transcriptome] {gid} from {src} -> {fa}")
        _ensure_dir(fa.parent)
        _download_to_fa(src, fa)
        produced.append(fa)

    return produced