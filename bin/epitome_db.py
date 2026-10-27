#!/usr/bin/env python3
"""
epitome_db.py - incremental NCBI virus genome store (DuckDB) for EPITOME.

  sync  : pull new/updated records for a taxon from NCBI datasets into the DB
  build : write {taxon}.json / {taxon}.fa.gz from the DB, in the same shape
          epitome_ncbi.py produced, so downstream steps can stay unchanged

Replaces epitome_ncbi.py. The three pure helpers from that script
(extract_taxids, parse_edirect_json_docsum, lineage_to_species_and_taxid)
are carried over unchanged in behavior.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import logging
import shlex
import subprocess
import sys
import tempfile
import zipfile
from collections import defaultdict
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Iterator, Optional

import duckdb

LOG = logging.getLogger("epitome_db")

DATA_REPORT = "ncbi_dataset/data/data_report.jsonl"
GENOMIC_FNA = "ncbi_dataset/data/genomic.fna"

SUBTYPE_TARGETS = {
    "subtype", "serotype", "genotype", "serogroup", "variety", "biotype",
    "lineage", "clade", "subclade", "pathotype", "group", "subgroup",
}

SCHEMA = """
CREATE TABLE IF NOT EXISTS records (
    accession         VARCHAR PRIMARY KEY,   -- versioned, e.g. MN908947.3
    accession_base    VARCHAR NOT NULL,
    version           INTEGER,
    length            INTEGER,
    segment           VARCHAR,
    source_database   VARCHAR,
    completeness      VARCHAR,
    host              VARCHAR,
    collection_date   VARCHAR,
    geographic_region VARCHAR,
    report            JSON,                  -- raw datasets record, never re-fetch for new fields
    report_sha        VARCHAR,
    status            VARCHAR DEFAULT 'current',  -- current | superseded | removed
    subtypes_fetched  BOOLEAN DEFAULT FALSE,
    first_seen        TIMESTAMP,
    last_seen         TIMESTAMP
);
CREATE TABLE IF NOT EXISTS sequences (
    accession VARCHAR PRIMARY KEY,           -- versioned accession => sequence is immutable
    sequence  VARCHAR NOT NULL
);
CREATE TABLE IF NOT EXISTS subtypes (
    accession VARCHAR, key VARCHAR, value VARCHAR,
    PRIMARY KEY (accession, key)
);
CREATE TABLE IF NOT EXISTS membership (      -- which taxon query returned which record
    taxon VARCHAR, accession VARCHAR, last_seen TIMESTAMP,
    PRIMARY KEY (taxon, accession)
);
CREATE TABLE IF NOT EXISTS taxa (taxid VARCHAR PRIMARY KEY, name VARCHAR);
CREATE TABLE IF NOT EXISTS sync_log (
    run_id BIGINT, taxon VARCHAR, mode VARCHAR, scope VARCHAR,
    started TIMESTAMP, finished TIMESTAMP, updated_since VARCHAR,
    n_in_pull INTEGER, n_new INTEGER, n_changed INTEGER, n_removed INTEGER,
    datasets_version VARCHAR, success BOOLEAN, error VARCHAR
);
"""

STAGE_COLS = [
    "accession", "accession_base", "version", "length", "segment",
    "source_database", "completeness", "host", "collection_date",
    "geographic_region", "report", "report_sha",
]
STAGE_DDL = """accession VARCHAR, accession_base VARCHAR, version INTEGER,
    length INTEGER, segment VARCHAR, source_database VARCHAR,
    completeness VARCHAR, host VARCHAR, collection_date VARCHAR,
    geographic_region VARCHAR, report VARCHAR, report_sha VARCHAR"""

UPSERT_SQL = f"""
INSERT INTO records ({", ".join(STAGE_COLS)}, status, first_seen, last_seen)
SELECT {", ".join(STAGE_COLS)}, 'current', ?, ? FROM stage
ON CONFLICT (accession) DO UPDATE SET
    {", ".join(f"{c} = excluded.{c}" for c in STAGE_COLS if c != "accession")},
    status = 'current', last_seen = excluded.last_seen
"""


# ---------------------------------------------------------------------------
# Helpers carried over from epitome_ncbi.py
# ---------------------------------------------------------------------------

def extract_taxids(ds_taxa_json: list[dict]) -> dict[str, str]:
    """taxid -> species name from `datasets summary taxonomy` JSON."""
    out: dict[str, str] = {}
    if not ds_taxa_json or not isinstance(ds_taxa_json[0], dict):
        LOG.warning("Unexpected taxonomy JSON structure; returning empty map")
        return out
    for row in ds_taxa_json[0].get("reports", []):
        taxonomy = row.get("taxonomy", {})
        name = taxonomy.get("current_scientific_name", {}).get("name")
        primary = taxonomy.get("tax_id")
        secondary = taxonomy.get("secondary_tax_ids", [])
        if not name or not primary:
            continue
        for tid in [primary] + (secondary if isinstance(secondary, list) else []):
            out[str(tid)] = name
    return out


def parse_edirect_json_docsum(objs) -> dict[str, dict[str, str]]:
    """accession -> {subtype-like key: value} from EDirect esummary JSON."""
    data: dict[str, dict[str, str]] = {}
    for row in objs:
        result = row.get("result", row)
        if not isinstance(result, dict):
            continue
        uids = result.get("uids")
        records = [result.get(u, {}) for u in uids] if isinstance(uids, list) else [result]
        for value in records:
            if not isinstance(value, dict):
                continue
            accession = value.get("accessionversion") or value.get("AccessionVersion")
            if not accession:
                continue
            keys = value.get("SubType") or value.get("subtype") or []
            vals = value.get("SubName") or value.get("subname") or []
            if isinstance(keys, str):
                keys = keys.split("|")
            if isinstance(vals, str):
                vals = vals.split("|")
            m = {}
            for i, k in enumerate(keys):
                k_norm = (k or "").strip().lower()
                if k_norm in SUBTYPE_TARGETS and i < len(vals):
                    v = str(vals[i]).strip()
                    if v:
                        m[k_norm] = v
            if m:
                data[accession] = m
    return data


def lineage_to_species_and_taxid(lineage_list, taxid_name_map):
    """First lineage entry whose taxId is in the taxid->name map."""
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


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

def _split_version(acc: str) -> tuple[str, Optional[int]]:
    base, _, ver = acc.partition(".")
    return base, int(ver) if ver.isdigit() else None


def report_to_row(rec: dict) -> Optional[dict]:
    """Flatten one datasets virus report. Tolerates null sub-objects."""
    acc = rec.get("accession")
    if not acc:
        return None
    base, ver = _split_version(acc)
    host_lineage = (rec.get("host") or {}).get("lineage") or []
    host = host_lineage[-1].get("name") if host_lineage else None
    raw = json.dumps(rec, sort_keys=True, separators=(",", ":"))
    length = rec.get("length")
    return {
        "accession": acc,
        "accession_base": base,
        "version": ver,
        "length": int(length) if length else None,
        "segment": rec.get("segment"),
        "source_database": rec.get("sourceDatabase") or "",
        "completeness": rec.get("completeness") or "",
        "host": host,
        "collection_date": (rec.get("isolate") or {}).get("collectionDate"),
        "geographic_region": (rec.get("location") or {}).get("geographicRegion"),
        "report": raw,
        "report_sha": hashlib.sha256(raw.encode()).hexdigest(),
    }


def iter_fasta(fh) -> Iterator[tuple[str, str]]:
    """Stream (accession, sequence); accession = first token of the header."""
    acc, parts = None, []
    for line in fh:
        line = line.rstrip()
        if line.startswith(">"):
            if acc is not None:
                yield acc, "".join(parts)
            tokens = line[1:].split()
            acc, parts = (tokens[0] if tokens else None), []
        elif acc is not None:
            parts.append(line)
    if acc is not None:
        yield acc, "".join(parts)


# ---------------------------------------------------------------------------
# NCBI calls (verify flags against `datasets ... --help` for your CLI version)
# ---------------------------------------------------------------------------

def datasets_version() -> str:
    try:
        out = subprocess.run(["datasets", "--version"], capture_output=True, text=True)
        return (out.stdout or out.stderr).strip()
    except FileNotFoundError:
        return "unknown"


def run_datasets_download(taxon: str, zip_path: Path, updated_since: Optional[str],
                          complete_only: bool) -> None:
    cmd = ["datasets", "download", "virus", "genome", "taxon", taxon,
           "--filename", str(zip_path), "--no-progressbar"]
    if complete_only:
        cmd.append("--complete-only")
    if updated_since:
        cmd += ["--updated-since", updated_since]
    LOG.info("Running: %s", " ".join(cmd))
    # NOTE: check how your datasets version behaves when nothing matches
    # --updated-since (empty zip vs non-zero exit) and handle it here if needed.
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        raise RuntimeError(f"datasets download failed ({res.returncode}): {res.stderr.strip()}")


def fetch_taxonomy(taxon: str, extra_args: str) -> list[dict]:
    """Use the same flags your current Nextflow module passes (via --taxonomy-args)."""
    cmd = ["datasets", "summary", "taxonomy", "taxon", taxon] + shlex.split(extra_args)
    res = subprocess.run(cmd, capture_output=True, text=True, check=True)
    return [json.loads(res.stdout)]


def fetch_subtypes(con, chunk: int = 200) -> int:
    """EDirect esummary for records that have not been looked up yet (new ones only)."""
    todo = [r[0] for r in con.execute(
        "SELECT accession FROM records WHERE NOT subtypes_fetched AND status = 'current'"
    ).fetchall()]
    done = 0
    for i in range(0, len(todo), chunk):
        ids = todo[i:i + chunk]
        try:
            res = subprocess.run(
                ["esummary", "-db", "nuccore", "-id", ",".join(ids), "-mode", "json"],
                capture_output=True, text=True, check=True)
            found = parse_edirect_json_docsum([json.loads(res.stdout)])
        except (subprocess.CalledProcessError, json.JSONDecodeError, FileNotFoundError) as e:
            LOG.warning("esummary failed for chunk starting %d (%s); will retry next sync", i, e)
            break  # leave subtypes_fetched = FALSE so the next run retries
        rows = [(a, k, v) for a, kv in found.items() for k, v in kv.items()]
        if rows:
            con.executemany(
                "INSERT INTO subtypes VALUES (?, ?, ?) ON CONFLICT (accession, key) "
                "DO UPDATE SET value = excluded.value", rows)
        con.execute("UPDATE records SET subtypes_fetched = TRUE WHERE list_contains(?, accession)", [ids])
        done += len(ids)
    LOG.info("Subtype lookup complete for %d/%d records", done, len(todo))
    return done


# ---------------------------------------------------------------------------
# Ingest
# ---------------------------------------------------------------------------

def ingest_zip(con, taxon: str, zip_path: Path, taxonomy: list[dict],
               run_started: datetime) -> dict:
    """Load one datasets virus package into the DB. Caller owns the transaction."""
    with zipfile.ZipFile(zip_path) as zf:
        names = set(zf.namelist())
        if DATA_REPORT not in names:
            LOG.warning("No %s in package; nothing to ingest", DATA_REPORT)
            return {"n_pull": 0, "n_new": 0, "n_changed": 0}

        con.execute(f"CREATE OR REPLACE TEMP TABLE stage ({STAGE_DDL})")
        placeholders = ", ".join("?" for _ in STAGE_COLS)
        seen: set[str] = set()
        batch: list[list] = []

        def flush():
            if batch:
                con.executemany(f"INSERT INTO stage VALUES ({placeholders})", batch)
                batch.clear()

        with io.TextIOWrapper(zf.open(DATA_REPORT), encoding="utf-8") as fh:
            for line in fh:
                if not line.strip():
                    continue
                row = report_to_row(json.loads(line))
                if not row or row["accession"] in seen:
                    continue
                seen.add(row["accession"])
                batch.append([row[c] for c in STAGE_COLS])
                if len(batch) >= 5000:
                    flush()
        flush()

        n_new = con.execute(
            "SELECT count(*) FROM stage s WHERE NOT EXISTS "
            "(SELECT 1 FROM records r WHERE r.accession = s.accession)").fetchone()[0]
        n_changed = con.execute(
            "SELECT count(*) FROM stage s JOIN records r ON r.accession = s.accession "
            "WHERE r.report_sha <> s.report_sha").fetchone()[0]

        con.execute(UPSERT_SQL, [run_started, run_started])
        con.execute(
            "INSERT INTO membership SELECT ?, accession, ? FROM stage "
            "ON CONFLICT (taxon, accession) DO UPDATE SET last_seen = excluded.last_seen",
            [taxon, run_started])

        tax_rows = list(extract_taxids(taxonomy).items())
        if tax_rows:
            con.executemany(
                "INSERT INTO taxa VALUES (?, ?) ON CONFLICT (taxid) "
                "DO UPDATE SET name = excluded.name", tax_rows)

        # Sequences: versioned accessions are immutable, so only load ones we lack.
        need = {a for (a,) in con.execute(
            "SELECT accession FROM stage WHERE accession NOT IN (SELECT accession FROM sequences)"
        ).fetchall()}
        if need and GENOMIC_FNA in names:
            seq_batch: list[tuple[str, str]] = []
            with io.TextIOWrapper(zf.open(GENOMIC_FNA), encoding="utf-8") as fh:
                for acc, seq in iter_fasta(fh):
                    if acc in need:
                        need.discard(acc)
                        seq_batch.append((acc, seq))
                        if len(seq_batch) >= 2000:
                            con.executemany("INSERT INTO sequences VALUES (?, ?)", seq_batch)
                            seq_batch.clear()
            if seq_batch:
                con.executemany("INSERT INTO sequences VALUES (?, ?)", seq_batch)
            if need:
                LOG.warning("%d accessions had no FASTA record (e.g. %s)", len(need), sorted(need)[:3])

    return {"n_pull": len(seen), "n_new": n_new, "n_changed": n_changed}


def maintain(con, taxon: str, run_started: datetime, full: bool) -> int:
    """Supersession everywhere; removal detection only after a full pull. Returns n_removed."""
    con.execute("""
        UPDATE records SET status = 'superseded'
        WHERE status = 'current' AND EXISTS (
            SELECT 1 FROM records n
            WHERE n.accession_base = records.accession_base AND n.version > records.version)""")
    n_removed = 0
    if full:
        con.execute("DELETE FROM membership WHERE taxon = ? AND last_seen < ?", [taxon, run_started])
        cond = ("status = 'current' AND NOT EXISTS "
                "(SELECT 1 FROM membership m WHERE m.accession = records.accession)")
        n_removed = con.execute(f"SELECT count(*) FROM records WHERE {cond}").fetchone()[0]
        if n_removed:
            con.execute(f"UPDATE records SET status = 'removed' WHERE {cond}")
    return n_removed


# ---------------------------------------------------------------------------
# Sync
# ---------------------------------------------------------------------------

def _log_run(con, taxon, mode, scope, started, since, counts, ver, success, error=None):
    con.execute(
        "INSERT INTO sync_log SELECT (SELECT coalesce(max(run_id), 0) + 1 FROM sync_log), "
        "?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?",
        [taxon, mode, scope, started, datetime.now(timezone.utc).replace(tzinfo=None), since,
         counts.get("n_pull"), counts.get("n_new"), counts.get("n_changed"),
         counts.get("n_removed"), ver, success, error])


def sync(con, taxon: str, force_full: bool = False, overlap_days: int = 3,
         complete_only: bool = False, subtypes: bool = True, taxonomy_args: str = "") -> dict:
    run_started = datetime.now(timezone.utc).replace(tzinfo=None)
    scope = f"complete_only={complete_only}"
    ver = datasets_version()

    since = None
    if not force_full:
        last = con.execute(
            "SELECT max(started) FROM sync_log WHERE taxon = ? AND success AND scope = ?",
            [taxon, scope]).fetchone()[0]
        if last:
            since = (last - timedelta(days=overlap_days)).strftime("%Y-%m-%d")
    full = since is None
    mode = "full" if full else "incremental"
    LOG.info("Sync '%s': mode=%s updated_since=%s", taxon, mode, since)

    counts: dict = {}
    try:
        with tempfile.TemporaryDirectory() as tmp:
            zip_path = Path(tmp) / "virus.zip"
            run_datasets_download(taxon, zip_path, since, complete_only)
            taxonomy = fetch_taxonomy(taxon, taxonomy_args)
            con.begin()
            try:
                counts = ingest_zip(con, taxon, zip_path, taxonomy, run_started)
                con.commit()
            except Exception:
                con.rollback()
                raise
        con.begin()
        try:
            counts["n_removed"] = maintain(con, taxon, run_started, full)
            con.commit()
        except Exception:
            con.rollback()
            raise
    except Exception as e:
        _log_run(con, taxon, mode, scope, run_started, since, counts, ver, False, str(e))
        raise

    _log_run(con, taxon, mode, scope, run_started, since, counts, ver, True)
    if subtypes:
        fetch_subtypes(con)  # best effort; never fails the sync
    LOG.info("Sync done: %s", counts)
    return counts


# ---------------------------------------------------------------------------
# Build (DB -> same per-taxon files epitome_ncbi.py wrote)
# ---------------------------------------------------------------------------

def build(con, taxon: str, out_dir: Path, min_length: int = 800) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    taxid_name_map = dict(con.execute("SELECT taxid, name FROM taxa").fetchall())

    subtype_map: dict[str, dict[str, str]] = defaultdict(dict)
    for acc, k, v in con.execute(
            "SELECT s.accession, s.key, s.value FROM subtypes s "
            "JOIN membership m ON m.accession = s.accession WHERE m.taxon = ?", [taxon]).fetchall():
        subtype_map[acc][k] = v

    query = """
        SELECT r.accession, r.segment, r.host, r.collection_date, r.geographic_region,
               r.source_database, r.completeness, CAST(r.report AS VARCHAR), q.sequence
        FROM records r
        JOIN membership m ON m.accession = r.accession
        LEFT JOIN sequences q ON q.accession = r.accession
        WHERE m.taxon = ? AND r.status = 'current' AND r.length >= ?
        ORDER BY r.accession"""

    rows: list[dict] = []
    seqs: list[tuple[str, str]] = []
    for acc, seg, host, cdate, region, src, comp, report, seq in con.execute(query, [taxon, min_length]).fetchall():
        lineage = ((json.loads(report).get("virus") or {}).get("lineage")) or []
        species, taxid = lineage_to_species_and_taxid(lineage, taxid_name_map)
        base = {
            "accession": acc, "taxon": taxon, "segment": seg, "species": species,
            "taxId": taxid, "host": host, "collectionDate": cdate,
            "geographicRegion": region,
            "canonical": src.lower() == "refseq",
        }
        if not seq:
            LOG.warning("No sequence stored for %s; skipping", acc)
            continue
        rows.append({**base, **subtype_map.get(acc, {})})
        seqs.append((acc, seq))

    json_path = out_dir / f"{taxon}.json"
    fasta_path = out_dir / f"{taxon}.fa.gz"
    with json_path.open("w", encoding="utf-8") as f:
        json.dump(rows, f, ensure_ascii=False, indent=2)
    with gzip.open(fasta_path, "wt", encoding="utf-8") as f:
        for acc, seq in seqs:
            f.write(f">{acc}\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i:i + 60] + "\n")

    acc_hash = hashlib.sha256("\n".join(a for a, _ in seqs).encode()).hexdigest()
    last_sync = con.execute("SELECT max(finished) FROM sync_log WHERE taxon = ? AND success", [taxon]).fetchone()[0]
    manifest = {
        "taxon": taxon, "n_records": len(rows), "min_length": min_length,
        "accession_list_sha256": acc_hash,
        "last_successful_sync": last_sync.isoformat() if last_sync else None,
        "built_at": datetime.now(timezone.utc).isoformat(),
    }
    (out_dir / f"{taxon}.manifest.json").write_text(json.dumps(manifest, indent=2))
    LOG.info("Built %s: %d records, accession hash %s", taxon, len(rows), acc_hash[:12])
    return manifest


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def connect(path: str):
    con = duckdb.connect(path)
    con.execute(SCHEMA)
    return con


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("sync", help="Pull new/updated records for a taxon into the DB")
    s.add_argument("--db", required=True)
    s.add_argument("--taxon", required=True)
    s.add_argument("--full", action="store_true", help="Ignore the watermark; also detects removals")
    s.add_argument("--overlap-days", type=int, default=3)
    s.add_argument("--complete-only", action="store_true",
                   help="Pass --complete-only to datasets (off by default: NCBI 'complete' drops good assemblies)")
    s.add_argument("--no-subtypes", action="store_true")
    s.add_argument("--taxonomy-args", default="", help="Extra flags for `datasets summary taxonomy taxon`")

    b = sub.add_parser("build", help="Write {taxon}.json / .fa.gz from the DB")
    b.add_argument("--db", required=True)
    b.add_argument("--taxon", required=True)
    b.add_argument("--out_dir", default=".")
    b.add_argument("--min-length", type=int, default=800)

    args = p.parse_args()
    con = connect(args.db)
    try:
        if args.cmd == "sync":
            sync(con, args.taxon, force_full=args.full, overlap_days=args.overlap_days,
                 complete_only=args.complete_only, subtypes=not args.no_subtypes,
                 taxonomy_args=args.taxonomy_args)
        else:
            build(con, args.taxon, Path(args.out_dir), args.min_length)
    except Exception as e:
        LOG.error("%s failed: %s", args.cmd, e)
        sys.exit(1)
    finally:
        con.close()


if __name__ == "__main__":
    main()