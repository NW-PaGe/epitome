#!/usr/bin/env python3

# epitome_ncbi.py
# Author: Jared Johnson, jared.johnson@doh.wa.gov

import argparse
import logging
import sys
from pathlib import Path
import screed
import re
import os
from typing import Optional, Iterable, Mapping, Sequence

from epitome_utils import load_json_or_jsonl, logging_config

# -------------------------------
#  GLOBAL CONFIG
# -------------------------------

LOGGER = logging_config()

# -------------------------------
#  FUNCTIONS
# -------------------------------

def load_fasta_map(path: str) -> dict[str, str]:
    """Load FASTA(.gz ok via screed) and map both full and base accessions to sequences."""
    LOGGER.info("Loading FASTA sequences: %s", path)
    seq_by_acc: dict[str, str] = {}
    count = 0
    for rec in screed.open(path):
        header = rec.name
        acc_full = header.split()[0]
        seq = rec.sequence
        if acc_full not in seq_by_acc:
            seq_by_acc[acc_full] = seq
        acc_base = acc_full.split(".", 1)[0]
        if acc_base not in seq_by_acc:
            seq_by_acc[acc_base] = seq
        count += 1
        if count % 10000 == 0:
            LOGGER.debug("Processed %d FASTA records...", count)
    LOGGER.info(f"Loaded {len(seq_by_acc)} sequences into FASTA map (keys include full and base accessions)"
    )
    return seq_by_acc


SUBTYPE_TARGETS = {
    "subtype",
    "serotype",
    "genotype",
    "serogroup",
    "variety",
    "biotype",
    "lineage",
    "clade",
    "subclade",
    "pathotype",
    "group",
    "subgroup",
}

def parse_edirect_json_docsum(edirect_json_obj: Iterable[Mapping]) -> dict[str, dict[str, str]]:
    """Extract subtype-like attributes from EDirect esummary/efetch docsum JSON."""
    LOGGER.info("Parsing EDirect JSON docsum for subtype-like attributes")
    data: dict[str, dict[str, str]] = {}
    rows = 0
    for row in edirect_json_obj:
        rows += 1
        result = row.get("result", row)
        if not isinstance(result, dict):
            LOGGER.debug("Skipping non-dict result at row %d", rows)
            continue
        uids = result.get("uids")
        records = [result.get(u, {}) for u in uids] if isinstance(uids, list) else [result]
        for value in records:
            if not isinstance(value, dict):
                continue
            accession = value.get("accessionversion") or value.get("AccessionVersion")
            if not accession:
                continue
            subtype_keys = value.get("SubType") or value.get("subtype") or []
            subtype_vals = value.get("SubName") or value.get("subname") or []
            if isinstance(subtype_keys, str):
                subtype_keys = subtype_keys.split("|")
            if isinstance(subtype_vals, str):
                subtype_vals = subtype_vals.split("|")
            m: dict[str, str] = {}
            for i, k in enumerate(subtype_keys):
                k_norm = (k or "").strip().lower()
                if k_norm in SUBTYPE_TARGETS and i < len(subtype_vals):
                    v = str(subtype_vals[i]).strip()
                    if v:
                        m[k_norm] = v
            if m:
                data[accession] = m
    LOGGER.info(f"Extracted subtype data for {len(data)} accessions (from {rows} rows)")
    return data


def extract_taxids(ds_taxa_json: list[dict]) -> dict[str, str]:
    """Extract a mapping from taxid to species name from NCBI datasets taxonomy JSON."""
    LOGGER.info("Extracting taxid -> species mapping")
    out: dict[str, str] = {}

    if not ds_taxa_json or not isinstance(ds_taxa_json[0], dict):
        LOGGER.warning("Unexpected taxonomy JSON structure; returning empty map")
        return out

    reports = ds_taxa_json[0].get("reports", [])
    for idx, row in enumerate(reports, 1):
        taxonomy = row.get("taxonomy", {})
        name = taxonomy.get("current_scientific_name", {}).get("name")
        primary_tid = taxonomy.get("tax_id")
        secondary_tids = taxonomy.get("secondary_tax_ids", [])

        if not name or not primary_tid:
            LOGGER.debug("No species info in report #%d", idx)
            continue

        tids = [primary_tid] + (secondary_tids if isinstance(secondary_tids, list) else [])
        for taxid in tids:
            if str(taxid) in out:
                LOGGER.debug("Overwriting existing mapping for taxid %s", taxid)
            out[str(taxid)] = name

    LOGGER.info("Extracted %d taxids", len(out))
    return out


def lineage_to_species_and_taxid(
    lineage_list: Optional[list], taxid_name_map: dict[object, str]
) -> tuple[Optional[str], Optional[object]]:
    """Map a lineage list to (species_name, taxid) using a taxid->name map."""
    if not isinstance(lineage_list, list):
        return None, None
    for item in lineage_list:
        taxid = item.get("taxId")
        if taxid is None:
            continue
        taxid = str(taxid)
        if taxid in taxid_name_map:
            return taxid_name_map[taxid], taxid
    return None, None


def write_taxon_files(
    taxon: str,
    rows: Sequence[dict],
    fasta_map: dict[str, str],
    out_dir: Path,
) -> None:
    """Write a single JSON and FASTA for the taxon, preserving raw segment info."""
    LOGGER.info("Writing consolidated outputs for taxon '%s' (n=%d)", taxon, len(rows))
    json_path = out_dir / f"{taxon}.json"
    fasta_path = out_dir / f"{taxon}.fa.gz"

    # JSON
    import json
    with json_path.open("w", encoding="utf-8") as f:
        json.dump(rows, f, ensure_ascii=False, indent=2)
    LOGGER.debug("Wrote JSON: %s", json_path)

    # FASTA (dedupe on accession)
    import gzip
    n_fasta = 0
    seen: set[str] = set()
    with gzip.open(f"{fasta_path}", "wt", encoding="utf-8") as f:
        for r in rows:
            acc = r.get("accession")
            if not acc:
                continue
            if acc in seen:
                continue
            seq = fasta_map.get(acc) or fasta_map.get(acc.split(".", 1)[0])
            if not seq:
                LOGGER.warning("No sequence found for accession %s in FASTA map", acc)
                continue
            seen.add(acc)
            f.write(f">{acc}\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i : i + 60] + "\n")
            n_fasta += 1
    LOGGER.info("Wrote FASTA: %s (records=%d)", fasta_path, n_fasta)


# -------------------------------
#  MAIN
# -------------------------------

def main() -> None:
    version = "2.1"

    """CLI entry point: merge NCBI datasets & EDirect docsum; write single taxon outputs."""
    parser = argparse.ArgumentParser(
        description="Merge NCBI datasets & EDirect subtype data, then write a single JSON/FASTA for the taxon."
    )
    parser.add_argument("--taxon", required=True, help="Taxon name.")
    parser.add_argument("--datasets_genome_fasta", required=True, help="Path to genomic.fna from 'datasets download virus genome taxon'.")
    parser.add_argument("--datasets_genome_json", required=True, help="Path to 'data_report.jsonl' from 'datasets download virus genome taxon'.")
    parser.add_argument("--datasets_taxonomy_json", required=True, help="Path to JSON from 'datasets summary taxonomy taxon'.")
    parser.add_argument("--edirect_json",default=None, help="Optional: Path to JSON or JSONL from EDirect esummary/efetch docsum. If omitted, no subtype fields are added.")
    parser.add_argument("--min-length",default=800, type=int, help="Minimum sequence length to be included.")
    parser.add_argument("--out_dir", default=".", help="Output directory.")
    parser.add_argument("--version", action="version", version=version, help="Show script version and exit.")
    args = parser.parse_args()

    LOGGER.info(f"{os.path.basename(__file__).replace('.py', '')} v{version}")
    LOGGER.info("Author: Jared Johnson")

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    LOGGER.debug(f"Output directory ensured: {out_dir}")

    try:
        fasta_map = load_fasta_map(args.datasets_genome_fasta)
    except Exception as e:
        LOGGER.exception(f"Failed to load FASTA: {e}")
        sys.exit(2)

    try:
        ds_genome_list = load_json_or_jsonl(args.datasets_genome_json)
        LOGGER.info(f"Loaded genome report ({args.datasets_genome_json})")
    except Exception as e:
        LOGGER.exception(f"Failed to load genome JSON/JSONL: {e}")
        sys.exit(2)

    try:
        ds_taxa_json = load_json_or_jsonl(args.datasets_taxonomy_json)
        LOGGER.info(f"Loaded taxonomy summary ({args.datasets_taxonomy_json})")
    except Exception as e:
        LOGGER.exception(f"Failed to load taxonomy JSON/JSONL: {e}")
        sys.exit(2)

    taxid_name_map = extract_taxids(ds_taxa_json)
    LOGGER.debug(f"taxid->species entries: {len(taxid_name_map)}")

    subtypes_map: dict[str, dict[str, str]] = {}
    if args.edirect_json:
        try:
            edirect_json_obj = load_json_or_jsonl(args.edirect_json)
            LOGGER.info("Loaded EDirect docsum (%s)", args.edirect_json)
        except Exception as e:
            LOGGER.exception("Failed to load EDirect JSON/JSONL: %s", e)
            sys.exit(2)

        subtypes_map = parse_edirect_json_docsum(edirect_json_obj)
        LOGGER.debug("subtype entries: %d", len(subtypes_map))
    else:
        LOGGER.info("No --edirect_json provided; proceeding without subtype annotations.")

    rows: list[dict] = []
    missing_accession = 0
    short_sequences   = 0
    for i, rec in enumerate(ds_genome_list, 1):
        accession = rec.get("accession")
        if not accession:
            missing_accession += 1
            if missing_accession <= 5:
                LOGGER.warning(f"Record #{i} missing 'accession'; skipping")
            continue

        length = rec.get('length')
        if not length or int(length) < int(args.min_length):
            short_sequences += 1
            continue

        # Preserve raw segment value exactly; if absent, None.
        segment = rec.get("segment")
        if segment is None:
            segment = None

        virus_lineage = (rec.get("virus") or {}).get("lineage") or []
        species_name, species_taxid = lineage_to_species_and_taxid(virus_lineage, taxid_name_map)

        host_lineage = rec.get("host", {}).get("lineage") or []
        host = host_lineage[-1].get("name") if host_lineage else None
        collection_date = rec.get("isolate", {}).get("collectionDate")
        geographic_region = rec.get("location", {}).get("geographicRegion")

        base = {
            "accession": accession,
            "taxon": args.taxon,
            "segment": segment,
            "species": species_name,
            "taxId": species_taxid,
            "host": host,
            "collectionDate": collection_date,
            "geographicRegion": geographic_region,
            "canonical": ((rec.get("sourceDatabase", '').lower() == "refseq" ) and (rec.get("completeness", '').lower() == "complete")),
        }


        sub = (subtypes_map.get(accession) or {})
        rows.append({**base, **sub})

        if i % 1000 == 0:
            LOGGER.debug(f"Processed {i} genome records...")

    if missing_accession:
        LOGGER.warning(f"Skipped {missing_accession} records without accession")

    if short_sequences:
        LOGGER.warning(f"Skipped {short_sequences} records with short sequences (<{args.min_length})")

    # Simple visibility on what we saw, without transforming names
    unique_segments = sorted({str(r.get("segment")) for r in rows})
    LOGGER.info(f"Observed {len(unique_segments)} unique segment values: {unique_segments}", )

    # Single consolidated outputs by taxon
    write_taxon_files(args.taxon, rows, fasta_map, out_dir)

    LOGGER.info(
        f"Processing complete for taxon '{args.taxon}' (rows={len(rows)}, outdir={out_dir})"
    )


if __name__ == "__main__":
    main()
