#!/usr/bin/env python3
"""
genome_downloader.py — Sequential, robust downloader for genomes.json.

Public API:
    download_genomes(json_path, outdir, *, layout="by-id", decompress=False,
                     retries=3, connect_timeout=30.0, read_timeout=900.0,
                     speed_limit_kbps=0, speed_time=30, resume=True, log=print)

Call from cli.py:
    from genome_downloader import download_genomes
    download_genomes("holosimulator_test/genomes.json", "holosimulator_test", decompress=True)
"""

from __future__ import annotations
import gzip
import io
import json
import os
import shutil
import socket
import sys
import time
from pathlib import Path
from typing import Callable, Iterable, List, Tuple
from urllib.parse import urlparse
from urllib.request import Request, build_opener, urlopen

LogFn = Callable[[str], None]

# ---------------- utilities ----------------

def _ts() -> str:
    return time.strftime("%Y-%m-%d %H:%M:%S")

def _default_log(msg: str) -> None:
    print(f"[{_ts()}] {msg}", file=sys.stderr, flush=True)

def _is_url(s: str) -> bool:
    try:
        return urlparse(s).scheme in ("http", "https", "ftp")
    except Exception:
        return False

def _ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)

# ---------------- core ops ----------------

def _http_download(
    url: str,
    part_path: Path,
    *,
    retries: int,
    connect_timeout: float,
    read_timeout: float,
    speed_limit_kbps: int,
    speed_time_window: int,
    resume: bool,
    log: LogFn,
) -> None:
    """
    Stream an HTTP(S) URL to part_path. Supports resume via Range if enabled.
    Aborts after bounded retries; includes a simple speed watchdog.
    """
    socket.setdefaulttimeout(connect_timeout)
    opener = build_opener()

    attempt = 0
    while True:
        attempt += 1
        try:
            headers = {"User-Agent": "holosim-downloader/1.0"}
            pos = 0
            if resume and part_path.exists():
                pos = part_path.stat().st_size
                if pos > 0:
                    headers["Range"] = f"bytes={pos}-"

            req = Request(url, headers=headers)
            with opener.open(req, timeout=connect_timeout) as resp:
                code = getattr(resp, "status", getattr(resp, "code", None))
                mode = "ab" if ("Range" in headers and code == 206) else "wb"

                window_start = time.time()
                bytes_in_window = 0
                CHUNK = 1024 * 1024

                with open(part_path, mode) as out:
                    while True:
                        # urllib uses the socket timeout; we also keep a speed watchdog below
                        chunk = resp.read(CHUNK)
                        if not chunk:
                            break
                        out.write(chunk)
                        bytes_in_window += len(chunk)

                        if speed_limit_kbps > 0:
                            now = time.time()
                            if now - window_start >= speed_time_window:
                                kbps = bytes_in_window / 1024.0 / (now - window_start)
                                if kbps < speed_limit_kbps:
                                    raise TimeoutError(
                                        f"Transfer stalled: {kbps:.1f} KB/s < {speed_limit_kbps} KB/s "
                                        f"over {speed_time_window}s"
                                    )
                                window_start = now
                                bytes_in_window = 0

            if part_path.stat().st_size == 0:
                raise IOError("Downloaded size is zero")
            return
        except Exception as e:
            if attempt <= retries:
                log(f"[retry {attempt}/{retries}] {url} -> {part_path.name}: {e}")
                time.sleep(2 * attempt)
                continue
            raise

def _ftp_download(url: str, part_path: Path, *, retries: int, log: LogFn) -> None:
    attempt = 0
    while True:
        attempt += 1
        try:
            with urlopen(url) as r, open(part_path, "wb") as out:
                shutil.copyfileobj(r, out, length=1024 * 1024)
            if part_path.stat().st_size == 0:
                raise IOError("Downloaded size is zero")
            return
        except Exception as e:
            if attempt <= retries:
                log(f"[retry {attempt}/{retries}] {url} -> {part_path.name}: {e}")
                time.sleep(2 * attempt)
                continue
            raise

def _copy_local(src: Path, part_path: Path) -> None:
    with open(src, "rb") as fin, open(part_path, "wb") as fout:
        shutil.copyfileobj(fin, fout, length=1024 * 1024)

def _maybe_decompress(gz_path: Path, fa_path: Path) -> None:
    with gzip.open(gz_path, "rb") as fin, open(fa_path, "wb") as fout:
        shutil.copyfileobj(fin, fout, length=1024 * 1024)

# ---------------- public API ----------------

def download_genomes(
    json_path: str | Path,
    outdir: str | Path,
    *,
    layout: str = "by-id",            # "by-id" -> out/genomes/G0001/G0001.fa.gz, "flat" -> out/genomes/G0001.fa.gz
    decompress: bool = False,         # also write .fa
    retries: int = 3,
    connect_timeout: float = 30.0,
    read_timeout: float = 900.0,      # used by socket default; long to allow big files
    speed_limit_kbps: int = 0,        # 0=disabled; e.g., 4 (KB/s) to detect stalls
    speed_time: int = 30,             # seconds window for speed watchdog
    resume: bool = True,              # resume HTTP(S) partial downloads
    log: LogFn = _default_log,        # pass `log=lambda m: None` to silence
) -> List[Path]:
    """
    Download all genomes defined in .genomes[] from JSON to outdir/genomes.
    Returns a list of produced .fa.gz Paths (and .fa if decompress=True).

    Sequential by design (no parallelism).
    """
    json_path = Path(json_path)
    outdir = Path(outdir)
    _ensure_dir(outdir)

    with open(json_path) as fh:
        meta = json.load(fh)

    genomes = meta.get("genomes", [])
    if not genomes:
        log("No genomes found in JSON (.genomes[]).")
        return []

    def _paths_for(gid: str) -> Tuple[Path, Path, Path]:
        if layout == "by-id":
            dest_dir = outdir / "genomes" / gid
            gz = dest_dir / f"{gid}.fa.gz"
            fa = dest_dir / f"{gid}.fa"
        elif layout == "flat":
            dest_dir = outdir / "genomes"
            gz = dest_dir / f"{gid}.fa.gz"
            fa = dest_dir / f"{gid}.fa"
        else:
            raise ValueError('layout must be "by-id" or "flat"')
        return dest_dir, gz, fa

    produced: List[Path] = []

    for g in genomes:
        gid: str = g["id"]
        src: str = g["path"]
        dest_dir, gz_path, fa_path = _paths_for(gid)
        _ensure_dir(dest_dir)

        if gz_path.exists() and gz_path.stat().st_size > 0:
            log(f"[skip] {gid} already present -> {gz_path}")
        else:
            part = gz_path.with_suffix(gz_path.suffix + ".part")
            log(f"[download] {gid} from {src} -> {part}")
            try:
                if _is_url(src):
                    scheme = urlparse(src).scheme
                    if scheme in ("http", "https"):
                        _http_download(
                            src,
                            part,
                            retries=retries,
                            connect_timeout=connect_timeout,
                            read_timeout=read_timeout,
                            speed_limit_kbps=speed_limit_kbps,
                            speed_time_window=speed_time,
                            resume=resume,
                            log=log,
                        )
                    elif scheme == "ftp":
                        _ftp_download(src, part, retries=retries, log=log)
                    else:
                        raise ValueError(f"Unsupported scheme: {scheme}")
                else:
                    _copy_local(Path(src), part)

                os.replace(part, gz_path)  # atomic publish within same directory/FS
                log(f"[publish]  {gid} -> {gz_path}")
            finally:
                # Clean orphaned zero-length .part if any
                try:
                    if part.exists() and part.stat().st_size == 0:
                        part.unlink(missing_ok=True)
                except Exception:
                    pass

        produced.append(gz_path)

        if decompress:
            if fa_path.exists() and fa_path.stat().st_size > 0:
                log(f"[skip] decompressed present -> {fa_path}")
            else:
                log(f"[decompress] {gid} -> {fa_path}")
                _maybe_decompress(gz_path, fa_path)
                produced.append(fa_path)

    log("All genomes processed.")
    return produced

# Optional: keep a tiny CLI for ad-hoc use
if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="Sequential genome downloader from JSON.")
    ap.add_argument("--json", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--layout", choices=["by-id", "flat"], default="by-id")
    ap.add_argument("--decompress", action="store_true")
    ap.add_argument("--retries", type=int, default=3)
    ap.add_argument("--connect-timeout", type=float, default=30.0)
    ap.add_argument("--read-timeout", type=float, default=900.0)
    ap.add_argument("--speed-limit-kbps", type=int, default=0)
    ap.add_argument("--speed-time", type=int, default=30)
    ap.add_argument("--no-resume", action="store_true")
    args = ap.parse_args()

    download_genomes(
        args.json,
        args.outdir,
        layout=args.layout,
        decompress=args.decompress,
        retries=args.retries,
        connect_timeout=args.connect_timeout,
        read_timeout=args.read_timeout,
        speed_limit_kbps=args.speed_limit_kbps,
        speed_time=args.speed_time,
        resume=not args.no_resume,
        log=_default_log,
    )
