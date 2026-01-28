#!/usr/bin/env python3
# epitome_inputs.py
# Author: Jared Johnson, jared.johnson@doh.wa.gov

import argparse
import json
import os
import sys
import gzip
from typing import Dict, List, Tuple, Set, Optional, Any

import sourmash as sm

from epitome_utils import (
    build_full_minhash_map,
    logging_config,
    read_jsonl,
    sanitize_filename,
)

LOGGER = logging_config()


# ============================================================================
# MINHASH OPERATIONS
# ============================================================================

def safe_contained_by(query: sm.MinHash, subject: sm.MinHash) -> float:
    """Calculate weighted containment, returning 0.0 for empty sketches."""
    if not getattr(query, "hashes", None) or not getattr(subject, "hashes", None):
        return 0.0
    try:
        return float(query.contained_by_weighted(subject))
    except ZeroDivisionError:
        return 0.0


def build_segment_globals(
    seg_mhs: Dict[str, List[sm.MinHash]],
    ksize: int,
    scaled: int,
) -> Dict[str, sm.MinHash]:
    """Build segment-level MinHash by summing abundances from member sketches."""
    segment_globals = {}
    
    for seg, members in seg_mhs.items():
        # Aggregate all hashes from member sketches
        combined_hashes: Dict[int, int] = {}
        for mh in members:
            for h, cnt in mh.hashes.items():
                combined_hashes[h] = combined_hashes.get(h, 0) + cnt
        
        # Create segment-level MinHash
        mh = sm.MinHash(n=0, scaled=scaled, ksize=ksize, track_abundance=True)
        mh.set_abundances(combined_hashes)
        segment_globals[seg] = mh
    
    return segment_globals


def compute_overlaps(
    globals_by_seg: Dict[str, sm.MinHash],
    threshold: float,
) -> Tuple[Dict[str, Set[str]], List[Tuple[str, str, float, float]]]:
    """
    Compute directed containment overlaps between segment globals.
    
    Returns:
        overlaps: {segment: set of overlapping segments}
        stats: [(seg1, seg2, seg1→seg2, seg2→seg1), ...]
    """
    segments = sorted(globals_by_seg.keys())
    overlaps = {s: set() for s in segments}
    stats = []
    
    for i, s1 in enumerate(segments):
        for s2 in segments[i + 1:]:
            # Bidirectional containment
            s1_in_s2 = safe_contained_by(globals_by_seg[s1], globals_by_seg[s2])
            s2_in_s1 = safe_contained_by(globals_by_seg[s2], globals_by_seg[s1])
            
            stats.append((s1, s2, s1_in_s2, s2_in_s1))
            
            if s1_in_s2 > threshold:
                overlaps[s1].add(s2)
            if s2_in_s1 > threshold:
                overlaps[s2].add(s1)
    
    return overlaps, stats


def log_overlaps(
    overlaps: Dict[str, Set[str]],
    stats: List[Tuple[str, str, float, float]],
    threshold: float,
    iteration: Optional[int] = None,
) -> None:
    """Log overlap statistics."""
    header = f"Iteration {iteration}" if iteration else "Final"
    LOGGER.info("=" * 70)
    LOGGER.info("%s overlaps (threshold=%.3f)", header, threshold)
    LOGGER.info("-" * 70)
    
    # Pairwise details
    for s1, s2, d12, d21 in stats:
        flag12 = "✓" if d12 > threshold else " "
        flag21 = "✓" if d21 > threshold else " "
        LOGGER.info(
            "  %s %-8s → %-8s: %.4f  |  %s %-8s → %-8s: %.4f",
            flag12, s1, s2, d12,
            flag21, s2, s1, d21,
        )
    
    # Summary
    LOGGER.info("-" * 70)
    for seg in sorted(overlaps.keys()):
        targets = sorted(overlaps[seg])
        LOGGER.info("  %-10s → %s", seg, ", ".join(targets) if targets else "(isolated)")

# ============================================================================
# FINAL MERGE LOGIC (POST-REASSIGNMENT)
# ============================================================================

def merge_overlapping_segments(
    seg_to_acc: Dict[str, List[str]],
    seg_to_mh: Dict[str, List[sm.MinHash]],
    ksize: int,
    scaled: int,
    merge_dist: float,
) -> Tuple[Dict[str, sm.MinHash], int]:
    """
    Evaluate final overlaps and merge segments when containment overlap exceeds merge_dist.
    The smaller segment (by member count) is merged into the larger segment.
    """
    merges = 0

    while True:
        globals_by_seg = build_segment_globals(seg_to_mh, ksize, scaled)
        overlaps, stats = compute_overlaps(globals_by_seg, merge_dist)

        # Find a merge candidate (deterministic): pick first pair in stats meeting threshold
        merge_pair = None
        for s1, s2, d12, d21 in stats:
            if d12 > merge_dist or d21 > merge_dist:
                merge_pair = (s1, s2, d12, d21)
                break

        if not merge_pair:
            return globals_by_seg, merges

        s1, s2, d12, d21 = merge_pair

        # Decide winner (keep) vs loser (merge into winner) by member count
        n1 = len(seg_to_acc.get(s1, []))
        n2 = len(seg_to_acc.get(s2, []))

        if n1 > n2:
            keep, drop = s1, s2
        elif n2 > n1:
            keep, drop = s2, s1
        else:
            # Tie-break deterministically
            keep, drop = (s1, s2) if s1 <= s2 else (s2, s1)

        LOGGER.info(
            "Merging segments (merge_dist=%.3f): %s(%d) + %s(%d) [%.4f/%.4f] → %s",
            merge_dist,
            keep, len(seg_to_acc.get(keep, [])),
            drop, len(seg_to_acc.get(drop, [])),
            d12, d21,
            keep,
        )

        # Move members + sketches
        seg_to_acc.setdefault(keep, []).extend(seg_to_acc.get(drop, []))
        seg_to_mh.setdefault(keep, []).extend(seg_to_mh.get(drop, []))

        # Remove dropped segment
        seg_to_acc.pop(drop, None)
        seg_to_mh.pop(drop, None)

        merges += 1


# ============================================================================
# DATA LOADING AND GROUPING
# ============================================================================

def get_segment(record: Dict[str, Any]) -> Optional[str]:
    """Extract segment from record (top-level or metadata)."""
    if seg := record.get("segment"):
        return seg
    if metadata := record.get("metadata"):
        if isinstance(metadata, dict):
            return metadata.get("segment")
    return None


def set_segment(record: Dict[str, Any], segment: Optional[str]) -> None:
    """Update segment in record (top-level and metadata if present)."""
    record["segment"] = segment
    if metadata := record.get("metadata"):
        if isinstance(metadata, dict):
            metadata["segment"] = segment


def load_records(
    input_paths: List[str]
) -> Tuple[Dict[str, str], Dict[str, str], Dict[str, str], Dict[str, Optional[str]], Dict[str, Dict]]:
    """
    Load and parse JSONL input files.
    
    Returns:
        canon_seq: {accession: sequence} for canonical
        canon_seg: {accession: segment} for canonical
        noncanon_seq: {accession: sequence} for non-canonical
        noncanon_seg: {accession: segment/None} for non-canonical
        all_records: {accession: record_dict} for all records
    """
    canon_seq, canon_seg = {}, {}
    noncanon_seq, noncanon_seg = {}, {}
    all_records = {}
    taxon = set()
    
    for path in input_paths:
        try:
            count = 0
            for i, rec in enumerate(read_jsonl(path), start=1):
                acc = rec.get("accession")
                if not acc:
                    raise ValueError(f"Missing accession in record {i}")
                
                seq = rec.get("sequence")
                if not seq:
                    raise ValueError(f"Missing sequence for {acc} in record {i}")

                _taxon = rec.get("taxon")
                if _taxon:
                    taxon.add(_taxon)
                
                seg = get_segment(rec)
                all_records[acc] = rec
                count += 1
                
                if seg and rec.get("canonical", False):
                    canon_seq[acc] = seq
                    canon_seg[acc] = seg
                else:
                    noncanon_seq[acc] = seq
                    noncanon_seg[acc] = seg

            if len(taxon) > 1:
                raise ValueError("Multiple taxon values in input!")
            else:
                taxon = list(taxon)[0]
            
            LOGGER.info("Loaded %d records from %s", count, path)
        
        except Exception as e:
            LOGGER.error("Failed to read %s: %s", path, e)
            sys.exit(2)
    
    if not canon_seq:
        raise ValueError("No canonical sequences found - at least one required")
    
    return canon_seq, canon_seg, noncanon_seq, noncanon_seg, all_records, taxon


def group_by_segment(
    canon_mh: Dict[str, sm.MinHash],
    canon_seg: Dict[str, str],
) -> Tuple[Dict[str, List[str]], Dict[str, List[sm.MinHash]]]:
    """Group canonical accessions and sketches by segment."""
    seg_to_acc = {}
    seg_to_mh = {}
    
    for acc, mh in canon_mh.items():
        seg = canon_seg[acc]
        seg_to_acc.setdefault(seg, []).append(acc)
        seg_to_mh.setdefault(seg, []).append(mh)
    
    return seg_to_acc, seg_to_mh


# ============================================================================
# REASSIGNMENT LOGIC
# ============================================================================

def plan_reassignments(
    overlaps: Dict[str, Set[str]],
    seg_to_acc: Dict[str, List[str]],
    canon_mh: Dict[str, sm.MinHash],
    globals_by_seg: Dict[str, sm.MinHash],
) -> List[Tuple[str, str, str]]:
    """
    Identify sequences that should be reassigned to different segments.
    
    For each segment with overlaps, check if member sequences fit better
    in overlapping segments.
    
    Returns: [(accession, old_segment, new_segment), ...]
    """
    reassignments = []
    
    for seg, overlap_targets in overlaps.items():
        if not overlap_targets:
            continue
        
        # Test against original segment + all overlapping segments
        candidates = [seg] + sorted(overlap_targets)
        
        for acc in seg_to_acc.get(seg, []):
            acc_mh = canon_mh[acc]
            
            # Find best-fitting segment
            best_seg = seg
            best_score = safe_contained_by(acc_mh, globals_by_seg[seg])
            
            for cand in overlap_targets:
                score = safe_contained_by(acc_mh, globals_by_seg[cand])
                if score > best_score:
                    best_score = score
                    best_seg = cand
            
            # Reassign if improvement found
            if best_seg != seg:
                reassignments.append((acc, seg, best_seg))
    
    return reassignments


def apply_reassignments(
    reassignments: List[Tuple[str, str, str]],
    seg_to_acc: Dict[str, List[str]],
    seg_to_mh: Dict[str, List[sm.MinHash]],
    canon_mh: Dict[str, sm.MinHash],
) -> None:
    """Apply reassignments in-place to segment groupings."""
    for acc, old_seg, new_seg in reassignments:
        mh = canon_mh[acc]
        
        seg_to_acc[old_seg].remove(acc)
        seg_to_mh[old_seg].remove(mh)
        
        seg_to_acc.setdefault(new_seg, []).append(acc)
        seg_to_mh.setdefault(new_seg, []).append(mh)


def iterative_reassignment(
    seg_to_acc: Dict[str, List[str]],
    seg_to_mh: Dict[str, List[sm.MinHash]],
    canon_mh: Dict[str, sm.MinHash],
    ksize: int,
    scaled: int,
    threshold: float,
    max_iterations: int = 50,
    verbose: bool = True,
) -> Tuple[Dict[str, sm.MinHash], int]:
    """
    Iteratively reassign sequences until stable or max iterations reached.
    
    Returns: (final_segment_globals, total_reassignment_count)
    """
    total_moved = 0
    
    for iteration in range(1, max_iterations + 1):
        # Build current segment globals
        globals_by_seg = build_segment_globals(seg_to_mh, ksize, scaled)
        
        # Check for empty segments
        empty = [s for s, mh in globals_by_seg.items() if not mh.hashes]
        if empty:
            LOGGER.warning("Iteration %d: Empty segments: %s", iteration, ", ".join(empty))
        
        # Compute overlaps
        overlaps, stats = compute_overlaps(globals_by_seg, threshold)
        
        if verbose:
            log_overlaps(overlaps, stats, threshold, iteration)
        
        # Plan reassignments
        reassignments = plan_reassignments(overlaps, seg_to_acc, canon_mh, globals_by_seg)
        
        if not reassignments:
            LOGGER.info("Iteration %d: Stable (no reassignments)", iteration)
            return globals_by_seg, total_moved
        
        # Log reassignments
        LOGGER.info("Iteration %d: %d reassignments planned", iteration, len(reassignments))
        for acc, old, new in reassignments[:20]:
            LOGGER.info("  %s: %s → %s", acc, old, new)
        if len(reassignments) > 20:
            LOGGER.info("  ... +%d more", len(reassignments) - 20)
        
        # Apply reassignments
        apply_reassignments(reassignments, seg_to_acc, seg_to_mh, canon_mh)
        total_moved += len(reassignments)
    
    LOGGER.warning("Max iterations (%d) reached without stabilizing", max_iterations)
    return build_segment_globals(seg_to_mh, ksize, scaled), total_moved


# ============================================================================
# ASSIGNMENT OF NON-CANONICAL SEQUENCES
# ============================================================================

def assign_noncanonical(
    noncanon_seq: Dict[str, str],
    noncanon_seg: Dict[str, Optional[str]],
    final_globals: Dict[str, sm.MinHash],
    seg_to_acc: Dict[str, List[str]],
    ksize: int,
    scaled: int,
    threshold: float,
    batch_size: int = 10000,
) -> Tuple[Dict[str, str], List[str]]:
    """
    Assign non-canonical sequences to segments.

    Returns: (assigned_map, unassigned_list)
        assigned_map: {accession: segment}
        unassigned_list: [accession, ...]
    """
    assigned: Dict[str, str] = {}
    unassigned: List[str] = []

    items = list(noncanon_seq.items())
    total = len(items)

    LOGGER.info(f"Assigning non-canonical in batches of {batch_size}")

    for start in range(0, total, batch_size):
        end = min(start + batch_size, total)
        batch = dict(items[start:end])

        # Build sketches for this batch only (limits memory)
        noncanon_mh = build_full_minhash_map(batch, ksize, scaled)

        batch_assigned = 0
        batch_unassigned = 0

        for acc, acc_mh in noncanon_mh.items():
            # Try original segment first if available
            orig_seg = noncanon_seg.get(acc)
            if orig_seg and orig_seg in final_globals:
                score = safe_contained_by(acc_mh, final_globals[orig_seg])
                if score >= threshold:
                    seg_to_acc.setdefault(orig_seg, []).append(acc)
                    assigned[acc] = orig_seg
                    batch_assigned += 1
                    continue

            # Find best-fitting segment
            best_seg = None
            best_score = -1.0

            for seg, global_mh in final_globals.items():
                score = safe_contained_by(acc_mh, global_mh)
                if score > best_score:
                    best_score = score
                    best_seg = seg

            if best_seg and best_score >= threshold:
                seg_to_acc.setdefault(best_seg, []).append(acc)
                assigned[acc] = best_seg
                batch_assigned += 1
            else:
                unassigned.append(acc)
                batch_unassigned += 1

        LOGGER.info(
            "Non-canonical batch %d–%d / %d: %d assigned, %d unassigned",
            start + 1,
            end,
            total,
            batch_assigned,
            batch_unassigned,
        )

    return assigned, unassigned




# ============================================================================
# OUTPUT GENERATION
# ============================================================================

def write_jsonl(path: str, records: List[Dict[str, Any]]) -> None:
    """Write records to JSONL file."""
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8") as f:
        for rec in records:
            f.write(json.dumps(rec, ensure_ascii=False) + "\n")


def prepare_output_record(record: Dict[str, Any], new_segment: Optional[str]) -> Dict[str, Any]:
    """Prepare record for output with original and updated segment info."""
    out = dict(record)
    out["segment_raw"] = get_segment(record)
    set_segment(out, new_segment)
    return out


def write_outputs(
    seg_to_acc: Dict[str, List[str]],
    unassigned: List[str],
    all_records: Dict[str, Dict[str, Any]],
    taxon: str,
    outdir: str,
) -> int:
    """Write segment-specific and unassigned JSONL files."""
    total_written = 0
    
    # Write per-segment files
    for seg in sorted(seg_to_acc.keys()):
        records = []
        for acc in seg_to_acc[seg]:
            if rec := all_records.get(acc):
                records.append(prepare_output_record(rec, seg))
        
        if records:
            filename = sanitize_filename(f"{taxon}-{seg}.jsonl.gz")
            path = os.path.join(outdir, filename)
            write_jsonl(path, records)
            total_written += len(records)
            LOGGER.info("Wrote %d records → %s", len(records), filename)
    
    # Write unassigned file
    unassigned_records = []
    for acc in unassigned:
        if rec := all_records.get(acc):
            unassigned_records.append(prepare_output_record(rec, None))
    
    if unassigned_records:
        path = os.path.join(outdir, "unassigned.jsonl")
        write_jsonl(path, unassigned_records)
        LOGGER.info("Wrote %d records → unassigned.jsonl", len(unassigned_records))
    
    return total_written


# ============================================================================
# MAIN
# ============================================================================

def main() -> None:
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Compute and iteratively refine segment assignments using MinHash containment"
    )
    parser.add_argument("inputs", nargs="+", help="Input JSONL/JSONL.GZ files")
    parser.add_argument("--ksize", type=int, default=16, help="k-mer size")
    parser.add_argument("--scaled", type=int, default=10, help="MinHash scaled factor")
    parser.add_argument("--containment", type=float, default=0.1, help="Containment threshold")
    parser.add_argument("--merge-dist", type=float, default=0.25, help="Merge threshold for final segment overlaps")
    parser.add_argument("--outdir", default="./", help="Output directory")
    parser.add_argument("--max-iters", type=int, default=50, help="Maximum iterations")
    parser.add_argument("--quiet-iters", action="store_true", help="Suppress per-iteration logging")
    parser.add_argument("--version", action="version", version="1.0")
    args = parser.parse_args()

    # Log configuration
    LOGGER.info("=" * 70)
    LOGGER.info("epitome_inputs v1.0")
    LOGGER.info("=" * 70)
    LOGGER.info("Configuration:")
    LOGGER.info("  Inputs: %d files", len(args.inputs))
    LOGGER.info("  k-mer size: %d", args.ksize)
    LOGGER.info("  Scaled: %d", args.scaled)
    LOGGER.info("  Containment threshold: %.3f", args.containment)
    LOGGER.info("  Merge threshold: %.3f", args.merge_dist)
    LOGGER.info("  Output directory: %s", args.outdir)
    LOGGER.info("  Max iterations: %d", args.max_iters)

    os.makedirs(args.outdir, exist_ok=True)

    # Load data
    LOGGER.info("-" * 70)
    canon_seq, canon_seg, noncanon_seq, noncanon_seg, all_records, taxon = load_records(args.inputs)
    LOGGER.info("Loaded: %d canonical, %d non-canonical sequences",
                len(canon_seq), len(noncanon_seq))

    # Build canonical sketches and group by segment
    canon_mh = build_full_minhash_map(canon_seq, args.ksize, args.scaled)
    seg_to_acc, seg_to_mh = group_by_segment(canon_mh, canon_seg)

    segment_summary = ", ".join(f"{s}({len(a)})" for s, a in sorted(seg_to_acc.items()))
    LOGGER.info("Initial segments: %s", segment_summary)

    # Iteratively reassign
    LOGGER.info("-" * 70)
    final_globals, total_moves = iterative_reassignment(
        seg_to_acc, seg_to_mh, canon_mh,
        args.ksize, args.scaled, args.containment,
        args.max_iters, not args.quiet_iters
    )
    LOGGER.info("Total reassignments: %d", total_moves)

    # Log final overlaps (containment threshold)
    overlaps, stats = compute_overlaps(final_globals, args.containment)
    log_overlaps(overlaps, stats, args.containment)

    n_overlapping = sum(1 for ov in overlaps.values() if ov)
    n_isolated = len(final_globals) - n_overlapping
    LOGGER.info("Final: %d overlapping segments, %d isolated", n_overlapping, n_isolated)

    # Merge overlapping segments (merge_dist), smaller into larger by member count
    LOGGER.info("-" * 70)
    final_globals, n_merges = merge_overlapping_segments(
        seg_to_acc, seg_to_mh,
        args.ksize, args.scaled,
        args.merge_dist,
    )
    LOGGER.info("Final merges applied: %d", n_merges)

    segment_summary = ", ".join(f"{s}({len(a)})" for s, a in sorted(seg_to_acc.items()))
    LOGGER.info("Post-merge segments: %s", segment_summary)

    # Assign non-canonical sequences
    LOGGER.info("-" * 70)
    assigned_noncanon, unassigned = assign_noncanonical(
        noncanon_seq, noncanon_seg, final_globals, seg_to_acc,
        args.ksize, args.scaled, args.containment
    )
    LOGGER.info("Non-canonical: %d assigned, %d unassigned",
                len(assigned_noncanon), len(unassigned))

    # Write outputs
    LOGGER.info("-" * 70)
    total_written = write_outputs(seg_to_acc, unassigned, all_records, taxon, args.outdir)
    LOGGER.info("Total assigned records written: %d", total_written)
    LOGGER.info("=" * 70)


if __name__ == "__main__":
    main()