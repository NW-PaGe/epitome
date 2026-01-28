#!/usr/bin/env python3

# epitome_inputs.py
# Author: Jared Johnson, jared.johnson@doh.wa.gov

import argparse
import gzip
import json
import logging
import os
import re
import sys
from typing import Any, Dict, List, Tuple

from epitome_utils import (
    detect_and_read,
    load_fastas,
    logging_config,
    normalize_keys,
)

# -------------------------------
#  GLOBAL CONFIG
# -------------------------------

LOGGER = logging_config()

STRING_KEYS = ["accession", "segment", "tax_id"]
YEAR_KEYS = ["collection_date"]
MISSING_DATA = ["na", "n/a", "null", "none", "nan"]


# -------------------------------
#  FUNCTIONS
# -------------------------------

def _sanitize_str(s: str) -> str:
    """Keep letters, numbers, underscores, dashes, and dots; replace others with '_'."""
    s0 = s
    s = s.strip().replace(" ", "_")
    s = re.sub(r"[^A-Za-z0-9._-]", "_", s)
    s = re.sub(r"__+", "_", s)
    s = s.lower()
    if s != s0:
        LOGGER.debug("Sanitized string '%s' -> '%s'", s0, s)
    return s


def _extract_year(rec: Dict[str, Any]) -> None:
    """Extract a 4-digit year from YEAR_KEYS fields and replace the value with an int if found."""
    year_re = re.compile(r"^\d{4}")
    for k in YEAR_KEYS:
        v = rec.get(k)
        if v:
            m = year_re.search(str(v))
            if m:
                old = rec[k]
                rec[k] = int(m.group())
                if rec[k] != old:
                    LOGGER.debug("Extracted year for key '%s': %s -> %s", k, old, rec[k])


def _coerce_to_str(rec: Dict[str, Any]) -> None:
    """Convert certain keys in a record to strings."""
    for k in STRING_KEYS:
        if k in rec and rec[k] is not None and rec[k] != "" and not isinstance(rec[k], str):
            old = rec[k]
            rec[k] = str(rec[k])
            LOGGER.debug("Coerced key '%s' to str: %r -> %r", k, old, rec[k])


def _update_missing(rec: Dict[str, Any]) -> None:
    """Replace certain missing value markers with None."""
    for k, v in list(rec.items()):
        if isinstance(v, str) and v.lower() in MISSING_DATA:
            LOGGER.debug("Setting missing value for key '%s': %r -> None", k, v)
            rec[k] = None


def merge_records(
    metadata: List[Tuple[str, List[Dict[str, Any]]]],
    fasta: Dict[str, str],
    segmented: bool
) -> List[Dict[str, Any]]:
    """Merge normalized metadata records with sequences from provided FASTA map."""
    LOGGER.info("Merging metadata with FASTA sequences")
    merged: List[Dict[str, Any]] = []
    required_keys = ["accession"]
    accs: set[str] = set()

    for src_path, recs in metadata:
        LOGGER.info("Source %s: %d metadata records", src_path, len(recs))
        for rec in recs:
            rec = normalize_keys(rec)
            _coerce_to_str(rec)
            _extract_year(rec)
            _update_missing(rec)

            if not all(k in rec for k in required_keys):
                raise ValueError(
                    f"One or more required key ({required_keys}) is missing from record {rec}"
                )

            acc = rec["accession"]
            if acc not in fasta:
                LOGGER.warning("%s not found in supplied assembly files", acc)
                continue

            if acc in accs:
                raise ValueError(f"Multiple metadata entries for accession {acc}")
            accs.add(acc)

            rec["sequence"] = fasta[acc]
            
            seg = rec.get("segment") if segmented else 'wg'
            rec["segment"] = _sanitize_str(seg) if seg else None

            merged.append(rec)

    LOGGER.info("Merged %d records with sequences (from %d sources)", len(merged), len(metadata))
    return merged


def drop_na(data: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Remove keys that are None for all records (only considering keys that exist in records)."""
    if not data:
        return data

    n_rec = len(data)
    na_counts: Dict[str, int] = {}
    rm_set: set[str] = set()

    for rec in data:
        for k, v in rec.items():
            if v is None:
                na_counts[k] = na_counts.get(k, 0) + 1
                if na_counts[k] == n_rec:
                    rm_set.add(k)

    if rm_set:
        LOGGER.info("Dropping all-NA keys across %d records: %s", n_rec, ", ".join(sorted(rm_set)))
        for rec in data:
            for k in rm_set:
                rec.pop(k, None)

    return data


# -------------------------------
#  MAIN
# -------------------------------

def main() -> None:
    """CLI entry point for merging JSON/CSV metadata into JSONL.GZ."""
    version = "1.0"

    parser = argparse.ArgumentParser(
        description="Merge JSON/CSV by 'accession' into flat JSONL.GZ metadata."
    )
    parser.add_argument("--taxon", default="null", help="Taxon name.")
    parser.add_argument("--assembly", nargs="+", help="Assembly FASTA(s).")
    parser.add_argument("--metadata", nargs="+", help="Metadata file(s).")
    parser.add_argument("--segmented",action="store_true",help="Sequences are from a segmented virus",)
    parser.add_argument("--outdir", default="./", help="Output directory.")
    parser.add_argument("--version", action="version", version=version, help="Show script version and exit.")
    args = parser.parse_args()

    LOGGER.info("%s v%s", os.path.basename(__file__).replace(".py", ""), version)
    LOGGER.info("Author: Jared Johnson")
    LOGGER.info(
        "Args: taxon=%s, segmented=%s, outdir=%s",
        args.taxon,
        args.segmented,
        args.outdir,
    )

    os.makedirs(args.outdir, exist_ok=True)

    # Load all inputs via detector
    metadata: List[Tuple[str, List[Dict[str, Any]]]] = []
    for p in (args.metadata or []):
        try:
            src_path, recs = detect_and_read(p)
            LOGGER.info("Read %d metadata records from %s", len(recs), src_path)
            metadata.append((src_path, recs))
        except Exception:
            LOGGER.exception("Error reading %s", p)
            sys.exit(2)

    fastas = load_fastas(args.assembly)

    # Merge and write
    merged = merge_records(metadata, fastas, args.segmented)
    no_na = drop_na(merged)

    outfile = os.path.join(args.outdir, f"{args.taxon}.jsonl.gz")
    n_rows = 0

    with gzip.open(outfile, "wt", encoding="utf-8") as mf:
        buf = []
        buf_bytes = 0
        flush_at = 4 * 1024 * 1024  # ~4MB of text

        for row in no_na:
            if not isinstance(row, dict):
                continue
            s = json.dumps(row, separators=(",", ":")) + "\n"
            buf.append(s)
            n_rows += 1
            buf_bytes += len(s)

            if buf_bytes >= flush_at:
                mf.writelines(buf)
                buf.clear()
                buf_bytes = 0

        if buf:
            mf.writelines(buf)

    LOGGER.info("Wrote %d records to %s", n_rows, outfile)



if __name__ == "__main__":
    main()
