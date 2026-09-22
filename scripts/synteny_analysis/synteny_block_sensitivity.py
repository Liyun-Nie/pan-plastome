#!/usr/bin/env python3
"""
Block parameter sensitivity from cached PAF (--from-cache).

Recomputes dual-axis merge + inter-block events across a parameter grid
without re-running minimap2.
"""

import argparse
import itertools
from pathlib import Path

import pandas as pd

from collinear_block_merge import merge_alignment_segments
from synteny_block_analysis import detect_inter_block_events, paf_to_primary_records
from synteny_shared import (
    DEFAULT_CONFIG,
    DATA_ROOT,
    PAN_ROOT,
    cache_dir_from_cfg,
    load_config,
    ref_fasta_path,
    resolve_sample_fastas,
)
from window_rearrangement_test import SPECIES_LIST

try:
    from Bio import SeqIO
except ImportError:
    raise SystemExit("ERROR: biopython required")


def species_ref_lengths(ref_dir, species_list, cache_dir):
    """Return {species: ref_len_bp} from cached GenBank exports."""
    lengths = {}
    cache_dir = Path(cache_dir)
    for species in species_list:
        ref_fa = ref_fasta_path(ref_dir, species, cache_dir)
        lengths[species] = len(str(next(SeqIO.parse(ref_fa, "fasta")).seq))
    return lengths


def load_paf_records(paf_path, min_mapq, min_alen_floor=300):
    paf_text = Path(paf_path).read_text(encoding="utf-8")
    q_len = None
    for line in paf_text.splitlines():
        if line.strip() and not line.startswith("#"):
            q_len = int(line.split("\t")[1])
            break
    records = paf_to_primary_records(paf_text, min_mapq=min_mapq, min_alen=min_alen_floor)
    return q_len, records


def summarize_from_records(records, min_alen, max_gap, overlap_tol, translocation_jump, q_len):
    filtered = [r for r in records if r["alen"] >= min_alen]
    blocks = merge_alignment_segments(filtered, max_gap=max_gap, overlap_tol=overlap_tol)
    events = detect_inter_block_events(blocks, translocation_jump=translocation_jump)
    aligned_bp = sum(b["block_len"] for b in blocks)
    aligned_fraction = round(aligned_bp / q_len, 6) if q_len else 0.0
    return {
        "n_blocks": len(blocks),
        "n_inversions": sum(1 for e in events if e["event_type"] == "inversion"),
        "n_translocations": sum(1 for e in events if e["event_type"] == "translocation"),
        "has_synteny_event": len(events) > 0,
        "aligned_fraction": aligned_fraction,
        "n_inter_block_events": len(events),
    }


def analyze_paf_path(
    paf_path,
    min_mapq,
    min_alen,
    max_gap,
    overlap_tol,
    translocation_jump,
):
    q_len, records = load_paf_records(paf_path, min_mapq, min_alen_floor=min_alen)
    return summarize_from_records(
        records, min_alen, max_gap, overlap_tol, translocation_jump, q_len
    )


def run_sensitivity_grid(cfg, paf_root, ref_dir, species_list, sample_dir):
    mm = cfg.get("minimap2", {})
    blk = cfg.get("block", {})
    sens = cfg.get("block_sensitivity", {})
    min_mapq = int(mm.get("min_mapq", 1))
    max_gap = int(sens.get("max_gap", blk.get("max_gap", 200)))
    overlap_tol = int(blk.get("overlap_tol", 50))
    min_lens = [int(x) for x in sens.get("min_len", [300, 500, 1000])]
    jump_fracs = [float(x) for x in sens.get("translocation_jump_frac", [0.05, 0.10])]
    min_alen_floor = min(min_lens)

    ref_cache = Path(PAN_ROOT / "results/synteny_analysis/ref_cache")
    if not ref_cache.exists():
        ref_cache = Path(PAN_ROOT / "results/synteny_analysis/_align_work" / "ref_cache")
    ref_lens = species_ref_lengths(ref_dir, species_list, ref_cache)

    grid = []
    for min_alen, jump_frac in itertools.product(min_lens, jump_fracs):
        grid.append({
            "min_alen": int(min_alen),
            "translocation_jump_frac": float(jump_frac),
        })

    # Per grid cell accumulators
    stats = {
        (g["min_alen"], g["translocation_jump_frac"]): {
            "n_samples": 0,
            "n_inter_block_total": 0,
            "n_samples_with_events": 0,
            "n_blocks_sum": 0,
            "n_inversions_total": 0,
            "n_translocations_total": 0,
        }
        for g in grid
    }
    missing_paf_samples = 0

    n_total = 0
    for species in species_list:
        fastas = resolve_sample_fastas(species, sample_dir)
        if not fastas:
            continue
        ref_len = ref_lens[species]
        sp_dir = Path(paf_root) / species
        for sid in fastas:
            n_total += 1
            paf_path = sp_dir / f"{sid}.paf"
            if not paf_path.exists():
                missing_paf_samples += 1
                continue
            q_len, records = load_paf_records(paf_path, min_mapq, min_alen_floor)
            for g in grid:
                key = (g["min_alen"], g["translocation_jump_frac"])
                trans_jump = max(1, int(round(ref_len * g["translocation_jump_frac"])))
                row = summarize_from_records(
                    records,
                    min_alen=g["min_alen"],
                    max_gap=max_gap,
                    overlap_tol=overlap_tol,
                    translocation_jump=trans_jump,
                    q_len=q_len,
                )
                acc = stats[key]
                acc["n_samples"] += 1
                acc["n_inter_block_total"] += row["n_inter_block_events"]
                acc["n_samples_with_events"] += int(row["has_synteny_event"])
                acc["n_blocks_sum"] += row["n_blocks"]
                acc["n_inversions_total"] += row["n_inversions"]
                acc["n_translocations_total"] += row["n_translocations"]

    rows = []
    for g in grid:
        key = (g["min_alen"], g["translocation_jump_frac"])
        acc = stats[key]
        mean_ref = sum(ref_lens.values()) / len(ref_lens) if ref_lens else 0
        jump_bp_typical = int(round(mean_ref * g["translocation_jump_frac"]))
        rows.append({
            "aligner": "minimap2",
            "min_alen": g["min_alen"],
            "max_gap": max_gap,
            "overlap_tol": overlap_tol,
            "translocation_jump_frac": g["translocation_jump_frac"],
            "translocation_jump_bp_typical": jump_bp_typical,
            "n_samples": acc["n_samples"],
            "missing_paf_samples": missing_paf_samples,
            "n_inter_block_total": acc["n_inter_block_total"],
            "n_inversions_total": acc["n_inversions_total"],
            "n_translocations_total": acc["n_translocations_total"],
            "n_samples_with_events": acc["n_samples_with_events"],
            "mean_n_blocks": round(acc["n_blocks_sum"] / acc["n_samples"], 6)
            if acc["n_samples"] else 0.0,
        })
    return pd.DataFrame(rows), n_total, missing_paf_samples


def main():
    parser = argparse.ArgumentParser(description="Block parameter sensitivity from cached PAF")
    parser.add_argument("--config", default=str(DEFAULT_CONFIG))
    parser.add_argument("--from-cache", action="store_true", required=True,
                        help="Recompute blocks from PAF cache (required)")
    parser.add_argument("--paf-dir", default=None, help="Override PAF cache directory")
    parser.add_argument("--ref-dir", default=str(DATA_ROOT / "ref/single_IR"))
    parser.add_argument("--sample-dir", default=str(DATA_ROOT / "00-rawfa/single_IR"))
    parser.add_argument("--output", default=str(
        PAN_ROOT / "results/synteny_analysis/sensitivity/block_param_sensitivity.tsv"
    ))
    parser.add_argument("--species", nargs="*", default=None)
    args = parser.parse_args()

    cfg = load_config(args.config)
    if "block_sensitivity" not in cfg:
        cfg["block_sensitivity"] = {
            "min_len": [300, 500, 1000],
            "translocation_jump_frac": [0.05, 0.10],
            "max_gap": 200,
        }

    paf_root = Path(args.paf_dir) if args.paf_dir else cache_dir_from_cfg(
        cfg, "paf", "results/synteny_analysis/paf"
    )
    species_list = args.species if args.species else SPECIES_LIST

    print(f"PAF cache: {paf_root}", flush=True)
    print(f"Grid: min_len={cfg['block_sensitivity']['min_len']}, "
          f"jump_frac={cfg['block_sensitivity']['translocation_jump_frac']}", flush=True)

    df, n_expected, missing = run_sensitivity_grid(
        cfg, paf_root, args.ref_dir, species_list, args.sample_dir
    )

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, sep="\t", index=False)

    print(f"Wrote {out}", flush=True)
    print(df.to_string(index=False), flush=True)
    if missing > 0:
        print(f"WARNING: missing PAF for {missing} samples", flush=True)
    print("Done.", flush=True)


if __name__ == "__main__":
    main()
