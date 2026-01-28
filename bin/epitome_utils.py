#!/usr/bin/env python3
# epitome_utils.py

from __future__ import annotations

import csv
import gzip
import json
import logging
import os
import re
import sys
from datetime import datetime
from typing import Any, Dict, List, Mapping, Tuple
import statistics
import matplotlib.pyplot as plt

import numpy as np
import screed
import sourmash as sm

# ---------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------


def logging_config(log_level: str = "INFO") -> logging.Logger:
    """Configure root logging once and return a script-named logger."""
    script_name = os.path.basename(sys.argv[0]).replace(".py", "")
    logger = logging.getLogger(script_name)

    # Configure logging (idempotent-ish; basicConfig is no-op if already configured)
    logging.basicConfig(
        level=getattr(logging, log_level.upper(), logging.INFO),
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    return logger


# ---------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------


def sanitize_filename(s: str) -> str:
    """Remove/replace characters unsafe for filesystem paths."""
    return re.sub(r"[^A-Za-z0-9._-]", "_", s)


def normalize_keys(d: Mapping[str, Any]) -> Dict[str, Any]:
    """Convert mapping keys to lowercase snake_case (simple camelCase -> snake)."""
    return {
        re.sub(r"([a-z])([A-Z])", r"\1_\2", str(k)).lower(): v
        for k, v in d.items()
    }


def is_date_string(v: Any) -> bool:
    """True if string parses via datetime.fromisoformat()."""
    if not isinstance(v, str):
        return False
    try:
        datetime.fromisoformat(v)
        return True
    except Exception:
        return False


# ---------------------------------------------------------------------
# JSON / CSV Readers
# ---------------------------------------------------------------------


def _open_maybe_gzip(path: str, mode: str):
    """Open normally or via gzip based on file suffix."""
    if path.endswith(".gz"):
        if "t" in mode:
            return gzip.open(path, mode, encoding="utf-8")
        return gzip.open(path, mode)

    # Plain file
    if "t" in mode:
        return open(path, mode, encoding="utf-8", newline="")
    return open(path, mode)


def read_json(path: str) -> List[Dict[str, Any]]:
    """Load a JSON or JSON.gz file whose root is a list of objects."""
    with _open_maybe_gzip(path, "rt") as f:
        obj = json.load(f)

    if not isinstance(obj, list):
        raise ValueError(f"{path}: JSON root must be a list of objects")

    for i, rec in enumerate(obj):
        if not isinstance(rec, dict):
            raise ValueError(f"{path}: element {i} is not an object")

    return obj


def read_csv(path: str) -> List[Dict[str, Any]]:
    """Load a CSV or CSV.gz into list of dicts; trims strings & converts empty to None."""
    recs: List[Dict[str, Any]] = []

    with _open_maybe_gzip(path, "rt") as f:
        reader = csv.DictReader(f)
        if reader.fieldnames is None:
            raise ValueError(f"{path}: missing CSV header row")

        for row in reader:
            clean: Dict[str, Any] = {}
            for k, v in row.items():
                k2 = k.strip() if isinstance(k, str) else k
                if isinstance(v, str):
                    v2 = v.strip()
                    v2 = v2 if v2 != "" else None
                else:
                    v2 = v
                clean[k2] = v2
            recs.append(clean)

    return recs


def detect_and_read(path: str) -> Tuple[str, List[Dict[str, Any]]]:
    """
    Read JSON/CSV (optionally gzipped), with JSON→CSV fallback.
    Extensions handled:
      - .json, .json.gz
      - .csv, .csv.gz
    """
    base, ext = os.path.splitext(path)
    ext = ext.lower()

    # Handle gz suffix: ext = ".gz" and base endswith .json/.csv
    if ext == ".gz":
        _, ext2 = os.path.splitext(base)
        ext2 = ext2.lower()
        if ext2 == ".json":
            return path, read_json(path)
        if ext2 == ".csv":
            return path, read_csv(path)
        # Otherwise fallback below

    if ext == ".json":
        return path, read_json(path)
    if ext == ".csv":
        return path, read_csv(path)

    # Fallback: try JSON then CSV
    try:
        return path, read_json(path)
    except Exception:
        return path, read_csv(path)


# ---------------------------------------------------------------------
# JSONL Helpers
# ---------------------------------------------------------------------


def read_jsonl(path: str) -> List[Dict[str, Any]]:
    """Load JSON Lines (one JSON object per line). Supports .jsonl and .jsonl.gz."""
    out: List[Dict[str, Any]] = []

    if path.endswith(".gz"):
        fh = gzip.open(path, "rt", encoding="utf-8")
    else:
        fh = open(path, "rt", encoding="utf-8")

    try:
        for line in fh:
            line = line.strip()
            if line:
                out.append(json.loads(line))
    finally:
        fh.close()

    logging.info("Loaded JSONL file: %s", path)
    return out


def load_json(path: str) -> List[Dict[str, Any]]:
    """
    For compatibility with prior usage: wrap a single JSON object into a list.
    If the file already contains a list, it becomes [that_list].
    Supports .json and .json.gz.
    """
    if path.endswith(".gz"):
        fh = gzip.open(path, "rt", encoding="utf-8")
    else:
        fh = open(path, "rt", encoding="utf-8")

    try:
        out = [json.load(fh)]
    finally:
        fh.close()

    logging.info("Loaded JSON file: %s", path)
    return out



def load_json_or_jsonl(path: str) -> List[Dict[str, Any]]:
    """Try JSON first, then JSONL."""
    try:
        return load_json(path)
    except json.JSONDecodeError:
        return read_jsonl(path)


def load_jsonl_gz_by_key(filepath: str, key: str) -> Dict[str, Dict[str, Any]]:
    """Load a .jsonl.gz into a dict keyed by `key`, removing `key` from each value dict."""
    logging.info("Loading file: %s using key: %s", filepath, key)
    data: Dict[str, Dict[str, Any]] = {}
    count = 0

    with gzip.open(filepath, "rt", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            entry = json.loads(line)
            if isinstance(entry, dict) and key in entry:
                data[str(entry[key])] = {k: v for k, v in entry.items() if k != key}
                count += 1

    logging.info("Loaded %d records from %s", count, filepath)
    return data


# ---------------------------------------------------------------------
# Distance cache + MinHash utilities
# ---------------------------------------------------------------------


class DistanceCache:
    """Symmetric key cache: cache[(scope, a, b)] == cache[(scope, b, a)]."""

    def __init__(self) -> None:
        self._d: Dict[Tuple[str, Any, Any], Any] = {}

    def _key(self, scope: str, a: Any, b: Any) -> Tuple[str, Any, Any]:
        return (scope, a, b) if a <= b else (scope, b, a)

    def get(self, scope: str, a: Any, b: Any) -> Any:
        return self._d.get(self._key(scope, a, b))

    def set(self, scope: str, a: Any, b: Any, val: Any) -> None:
        self._d[self._key(scope, a, b)] = val


def build_full_minhash_map(seqs: Dict[str, str], ksize: int, scaled: int) -> Dict[str, sm.MinHash]:
    """Build a MinHash object for each sequence."""
    mh_map: Dict[str, sm.MinHash] = {}
    for sid, s in seqs.items():
        mh = sm.MinHash(n=0, ksize=ksize, scaled=scaled, track_abundance=True)
        mh.add_sequence(s, force=True)
        mh_map[sid] = mh
    return mh_map


def _distance(
    a: str,
    b: str,
    *,
    scope: str,
    mh_map: Dict[str, sm.MinHash],
    cache: DistanceCache,
) -> float:
    """Cached pairwise distance between two IDs for a given scope ('full' or 'win:XYZ')."""
    cached = cache.get(scope, a, b)
    if cached is not None:
        return float(cached)

    d = mh_map[a].containment_ani(mh_map[b]).dist
    cache.set(scope, a, b, float(d))
    return float(d)


def compute_matrix(
    data: Dict[str, dict],
    window_scope: str,
    dist_cache: DistanceCache,
) -> Tuple[np.ndarray, List[str]]:
    """Compute pairwise distance matrix for sequences in a window using cached distances."""
    ids = sorted(data)
    n = len(ids)
    mat = np.zeros((n, n), dtype=float)
    mh_map = {sid: data[sid]["mh"] for sid in ids}

    for i in range(n):
        ai = ids[i]
        for j in range(i + 1, n):
            aj = ids[j]
            d = _distance(ai, aj, scope=window_scope, mh_map=mh_map, cache=dist_cache)
            mat[i, j] = mat[j, i] = d

    np.fill_diagonal(mat, 0.0)
    return mat, ids


# ---------------------------------------------------------------------
# FASTA loading
# ---------------------------------------------------------------------


def load_fastas(files: List[str]) -> Dict[str, str]:
    """Load one or more FASTA/FASTQ files (gz ok via screed) into a {name: sequence} map."""
    fasta: Dict[str, str] = {}
    n_files = 0
    n_seqs = 0

    for f in files:
        n_files += 1
        for rec in screed.open(f):
            name = rec.name.split()[0]
            if name in fasta:
                raise ValueError(f"Multiple sequences supplied for {name}")
            fasta[name] = rec.sequence
            n_seqs += 1

    return fasta

# ---------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------

# Modified Z-Score (Median + MAD)
def modified_zscore(data, z_threshold):
    if len(data) >= 2:
        med_len = statistics.median(data)
        abs_devs = [abs(x - med_len) for x in data]
        mad = statistics.median(abs_devs)
        mad_sigma = 1.4826 * mad

        if mad_sigma > 0:
            lower = med_len - (z_threshold * mad_sigma)
            upper = med_len + (z_threshold * mad_sigma)
        else:
            lower = upper = med_len

    elif len(data) == 1:
        med_len = data[0]
        mad = mad_sigma = 0.0
        lower = upper = med_len
    else:
        med_len = mad = mad_sigma = lower = upper = 0

    return med_len, mad, mad_sigma, lower, upper

# Standard z-score (mean / stddev)
def standard_zscore(data, z_threshold):
    if len(data) >= 2:
        mean = statistics.mean(data)
        stdev = statistics.stdev(data)  # sample standard deviation

        if stdev > 0:
            lower = mean - (z_threshold * stdev)
            upper = mean + (z_threshold * stdev)
        else:
            lower = upper = mean

    elif len(data) == 1:
        mean = data[0]
        stdev = 0.0
        lower = upper = mean
    else:
        mean = stdev = lower = upper = 0.0

    return mean, stdev, lower, upper

def plot_distribution(values, outfile, bins=50, title=None, cutoffs=None, xlab = "Value"):
    """
    Plot a histogram with lower, median, and upper vertical lines.

    Args:
        values (list or array-like): Numeric values to plot
        outfile (str): Path to output image
        bins (int): Number of histogram bins
        title (str, optional): Plot title
    """
    values = np.asarray(values)

    if cutoffs:
        lower, middle, upper = cutoffs

    plt.figure()
    plt.hist(values, bins=bins)
    plt.axvline(lower, color="red", linestyle="--", label=f"Lower: {lower:.2f}")
    plt.axvline(middle, color="black", linestyle="-",  label=f"Middle: {middle:.2f}")
    plt.axvline(upper, color="red", linestyle="--", label=f"Upper: {upper:.2f}")

    if title:
        plt.title(title)

    plt.xlabel(xlab)
    plt.ylabel("Count")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outfile)
    plt.close()