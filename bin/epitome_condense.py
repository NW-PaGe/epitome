#!/usr/bin/env python3

# epitome_condense.py
# Author: Jared Johnson, jared.johnson@doh.wa.gov
version = 2.0

import sourmash
import screed
import numpy as np
from sklearn.cluster import DBSCAN
import random
import argparse
import json
import os
from collections import defaultdict
import logging
import gzip
from typing import Dict, Tuple, List, Any
import math
import statistics

from epitome_utils import sanitize_filename, DistanceCache, logging_config

# -------------------------------
#  GLOBAL CONFIG
# -------------------------------

LOGGER = logging_config()

# -------------------------------
#  FUNCTIONS
# -------------------------------

def compute_matrix(
    data: Dict[str, dict],
    scope: str,
    dist_cache: DistanceCache
) -> Tuple[np.ndarray, List[str]]:
    """
    Compute pairwise distance matrix for sequences.

    Args:
        data: Mapping of sequence IDs to dict containing 'mh' key.
        scope: Cache key prefix for this comparison scope.
        dist_cache: Distance cache object.

    Returns:
        Tuple of (distance matrix, sorted list of sequence IDs).
    """
    ids = sorted(data.keys())
    n = len(ids)
    mat = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            id1, id2 = ids[i], ids[j]
            cached = dist_cache.get(scope, id1, id2)
            if cached is None:
                d = data[id1]['mh'].containment_ani(data[id2]['mh']).dist
                dist_cache.set(scope, id1, id2, d)
            else:
                d = cached
            mat[i, j] = mat[j, i] = d
    return mat, ids


def create_clusters(
    data: Dict[str, dict],
    threshold: float,
    start: int,
    stage: str,
    outdir: str,
    scope: str,
    dist_cache: DistanceCache
) -> Dict[str, int]:
    """
    Cluster sequences using DBSCAN.

    Args:
        data: Mapping of sequence IDs to dict containing 'mh'.
        threshold: DBSCAN epsilon distance.
        start: Starting cluster number (unused here for offsetting).
        stage: Stage name for logging.
        outdir: Output directory path.
        scope: Cache scope string.
        dist_cache: Distance cache object.

    Returns:
        Dict mapping sequence ID to cluster number.
    """
    mat, ids = compute_matrix(data, scope, dist_cache)
    db = DBSCAN(eps=threshold, min_samples=1, metric='precomputed').fit(mat)
    clusters = [int(i) + 1 for i in db.labels_]
    result = {}
    for id_, cluster in zip(ids, clusters):
        result[id_] = cluster
    return result


def select_best(data: Dict[str, dict]) -> Dict[str, dict]:
    """
    Select the best representative sequence for each merged cluster.
    """
    LOGGER.info('Selecting best sequence for each merged cluster')
    result = {}

    # Group seq ids by their merged_cluster key
    grouped = defaultdict(list)
    for sid, info in data.items():
        mcid = info.get('merged_cluster')
        if mcid:
            grouped[mcid].append(sid)

    # For each merged cluster, pick a winner by n, then seqLen
    for mcid, ids in grouped.items():
        # highest n wins
        max_n = max(int(data[i]['n']) for i in ids)
        candidates = [i for i in ids if int(data[i]['n']) == max_n]

        # tie-breaker: longest sequence
        if len(candidates) > 1:
            max_len = max(int(data[i]['seqLen']) for i in candidates)
            candidates = [i for i in candidates if int(data[i]['seqLen']) == max_len]

        winner = candidates[0]
        # copy winner and aggregate members from losers
        result[winner] = data[winner]
        for i in ids:
            if i != winner:
                LOGGER.info(f"Condensing {i} -> {winner}")
                result[winner]['members'].extend(data[i]['members'])

    return result


def condense_seqs(
    windows: Dict[str, Dict[str, dict]],
    data: Dict[str, dict],
    threshold: float,
    stage: int,
    outdir: str,
    dist_cache: DistanceCache
) -> Tuple[Dict[str, dict], str]:
    """
    Condense sequences by merging similar ones based on windowed MinHash clustering.

    Args:
        windows: Mapping of window keys to sequence ID -> {'mh': MinHash}.
        data: Mapping of sequence IDs to sequence metadata.
        threshold: Distance threshold for DBSCAN.
        stage: Current condensation stage number.
        outdir: Output directory.
        dist_cache: Distance cache object.

    Returns:
        Tuple of (updated data dict, status string: 'done' or 'not done').
    """
    LOGGER.info(f"Stage {stage} (n = {len(data)})")
    if len(data) == 1:
        LOGGER.info("Only one sequence - done.")
        return data, 'done'

    for win_key, win_data in windows.items():
        mh_subset = {k: win_data[k] for k in data if k in win_data}
        if len(mh_subset) == len(data):
            scope = f"win:{win_key}"
            window_clusters = create_clusters(mh_subset, threshold, 0, f'{win_key}_{stage}', outdir, scope, dist_cache)
            for cid, mcid in window_clusters.items():
                data[cid]['merged_cluster'] = str(data[cid].get('merged_cluster', '')) + str(mcid)

    merged_clusters = [v['merged_cluster'] for v in data.values() if 'merged_cluster' in v]
    LOGGER.info(f'Created {len(set(merged_clusters))} merged clusters')

    if len(data) == len(set(merged_clusters)):
        LOGGER.info("Nothing to condense")
        result, status = data, 'done'
    else:
        result = select_best(data)

        # stop if nothing reduced to avoid infinite loop
        if len(result) == len(data):
            LOGGER.warning("No reduction after selection; stopping to avoid infinite loop.")
            status = 'done'
        else:
            status = 'not done'

    for v in result.values():
        v.pop('merged_cluster', None)
    return result, status


def save_output(
    data: List[dict],
    taxon: str,
    segment: str,
    outdir: str
) -> None:
    """
    Save condensed clustering results to compressed JSONL.

    Args:
        data: List of condensed cluster dictionaries.
        taxon: Taxon name.
        segment: Segment name.
        outdir: Output directory.

    Returns:
        None
    """
    safe_taxon = sanitize_filename(taxon)
    safe_segment = sanitize_filename(segment)
    prefix = os.path.join(outdir, f"{safe_taxon}-{safe_segment}")
    outfile = prefix + ".condensed.jsonl.gz"
    with gzip.open(outfile, "wt", encoding="utf-8") as f:
        for row in data:
            json.dump(row, f, sort_keys=True)
            f.write("\n")
    LOGGER.info(f"Saved metrics to {outfile}")

# -------------------------------
#  MAIN
# -------------------------------

def main() -> None:
    """
    CLI entry point for condensing clustered sequences.

    Loads sequences and cluster metadata, computes MinHash sketches,
    performs iterative condensation until no further merging is possible,
    and saves condensed output.
    """
    version = "2.0"

    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", nargs='+', required=True, help="Path(s) to FASTA files.")
    parser.add_argument("--clusters", required=True, help="Path to cluster metadata JSONL.gz.")
    parser.add_argument("--taxon", default="null", help="Taxon name.")
    parser.add_argument("--segment", default="null", help="Segment name.")
    parser.add_argument("--dist", type=float, default=0.02, help="Distance threshold (1 - ANI/100).")
    parser.add_argument("--scaled", type=int, default=1000, help="MinHash scaled factor.")
    parser.add_argument("--ksize", type=int, default=31, help="MinHash k-mer size.")
    parser.add_argument("--window_size", type=int, default=8000, help="Window size for MinHash splitting.")
    parser.add_argument("--seed", type=int, default=11, help="Random seed.")
    parser.add_argument("--outdir", default="./", help="Output directory.")
    parser.add_argument("--version", action="version", version=version, help="Show script version and exit.")
    args = parser.parse_args()

    LOGGER.info(f"{os.path.basename(__file__).replace('.py', '')} v{version}")
    LOGGER.info(f"Author: Jared Johnson")

    os.makedirs(args.outdir, exist_ok=True)
    random.seed(args.seed)

    norm_taxon = args.taxon.replace(' ', '_')
    norm_segment = args.segment.replace(' ', '_')
    expected_prefix = f"{norm_taxon}-{norm_segment}-"

    data: Dict[str, dict] = {}
    window_size = args.window_size
    for fasta in args.fasta:
        for rec in screed.open(fasta):
            sid, seq = rec.name, rec.sequence

            if not sid.startswith(expected_prefix):
                raise ValueError(
                    f"{sid} is malformed for provided taxon/segment. "
                    f"Expected it to start with '{expected_prefix}'"
                )
            cid = sid[len(expected_prefix):].strip()
            data[str(cid)] = {
                "sequence": seq,
                "seqLen": len(seq)
            }
            window_size = min(window_size, len(seq))
    
    LOGGER.info(f"Loaded {len(data)} sequences.")

    windows: Dict[str, Dict[str, dict]] = {}
    for cid, info in data.items():
        j = 0
        for i in range(window_size - 1, info["seqLen"], window_size):
            mh = sourmash.MinHash(n=0, ksize=args.ksize, scaled=args.scaled)
            mh.add_sequence(info["sequence"][j:(i - 1)], True)
            win_key = f"{j}-{i}"
            windows.setdefault(win_key, {})[cid] = {"mh": mh}
            j = i

    start_len = len(data)

    with gzip.open(args.clusters, "rt", encoding="utf-8") as f:
        for line in f:
            row = json.loads(line.strip())
            if "cluster" in row and "members" in row:
                cid = str(row["cluster"])
                if cid in data:
                    data[str(cid)]["n"] = len(row["members"])
                    data[str(cid)]["members"] = row["members"]

    for key, value in data.items():
        req_subkeys = ['sequence', 'seqLen', 'members', 'n']
        for subkey in req_subkeys:
            if subkey not in value:
                raise ValueError(f"Problem finding {subkey} for sequence {key}")

    dist_cache = DistanceCache()

    stage = 0
    status = "not done"
    while status != "done":
        stage += 1
        data, status = condense_seqs(windows, data, args.dist, stage, args.outdir, dist_cache)

    condensed = [
        {
            "cluster": cid,
            "members": v.get("members"),
            "sequence": v.get("sequence")
        } for cid, v in data.items()
    ]

    save_output(condensed, args.taxon, args.segment, args.outdir)

    print(f"""
    Raw: {start_len} seqs
    Condensed: {len(condensed)} seqs

    Results saved to {args.outdir}
    """, flush=True)


if __name__ == "__main__":
    main()
