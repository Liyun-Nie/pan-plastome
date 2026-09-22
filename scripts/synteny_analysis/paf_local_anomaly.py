#!/usr/bin/env python3
"""
Scan cached PAF for within-chain local collinearity anomalies.

Flags adjacent segments (query gap <= max_query_gap) with:
  - strand_flip
  - ref_backtrack (ref coordinates non-monotonic on same strand)
"""

import argparse
from pathlib import Path

import pandas as pd

from paf_viz_shared import load_paf_segments, qc_pass_mask, species_paf_dir
from synteny_shared import DEFAULT_CONFIG, DATA_ROOT, PAN_ROOT, load_config, resolve_sample_fastas
from window_rearrangement_test import SPECIES_LIST

ANOMALY_COLS = [
    "species", "sample_id", "anomaly_type", "seg_i", "seg_j", "query_gap", "ref_gap",
    "qstart_i", "qend_i", "qstart_j", "qend_j", "tstart_i", "tend_i", "tstart_j", "tend_j",
    "strand_i", "strand_j",
]


def empty_anomaly_frame():
    return pd.DataFrame(columns=ANOMALY_COLS)


def detect_local_anomalies(
    segments,
    max_query_gap=200,
    overlap_tol=50,
    min_segment_len=500,
    report_ref_decrease=True,
):
    segs = [s for s in segments if s["alen"] >= min_segment_len]
    segs.sort(key=lambda s: (s["qstart"], s["tstart"]))
    hits = []
    for i in range(len(segs) - 1):
        a, b = segs[i], segs[i + 1]
        qgap = b["qstart"] - a["qend"]
        if qgap > max_query_gap:
            continue

        if a["strand"] != b["strand"]:
            hits.append({
                "anomaly_type": "strand_flip",
                "seg_i": i,
                "seg_j": i + 1,
                "query_gap": qgap,
                "ref_gap": None,
                "qstart_i": a["qstart"],
                "qend_i": a["qend"],
                "qstart_j": b["qstart"],
                "qend_j": b["qend"],
                "tstart_i": a["tstart"],
                "tend_i": a["tend"],
                "tstart_j": b["tstart"],
                "tend_j": b["tend"],
                "strand_i": a["strand"],
                "strand_j": b["strand"],
            })
            continue

        if not report_ref_decrease:
            continue

        if a["strand"] == "+":
            ref_gap = b["tstart"] - a["tend"]
            if ref_gap < -overlap_tol:
                hits.append({
                    "anomaly_type": "ref_backtrack",
                    "seg_i": i,
                    "seg_j": i + 1,
                    "query_gap": qgap,
                    "ref_gap": ref_gap,
                    "qstart_i": a["qstart"],
                    "qend_i": a["qend"],
                    "qstart_j": b["qstart"],
                    "qend_j": b["qend"],
                    "tstart_i": a["tstart"],
                    "tend_i": a["tend"],
                    "tstart_j": b["tstart"],
                    "tend_j": b["tend"],
                    "strand_i": a["strand"],
                    "strand_j": b["strand"],
                })
        else:
            ref_gap = a["tstart"] - b["tend"]
            if ref_gap < -overlap_tol:
                hits.append({
                    "anomaly_type": "ref_backtrack",
                    "seg_i": i,
                    "seg_j": i + 1,
                    "query_gap": qgap,
                    "ref_gap": ref_gap,
                    "qstart_i": a["qstart"],
                    "qend_i": a["qend"],
                    "qstart_j": b["qstart"],
                    "qend_j": b["qend"],
                    "tstart_i": a["tstart"],
                    "tend_i": a["tend"],
                    "tstart_j": b["tstart"],
                    "tend_j": b["tend"],
                    "strand_i": a["strand"],
                    "strand_j": b["strand"],
                })
    return hits


def process_species(species, cfg, paf_root, sample_dir, min_mapq, min_alen, anomaly_cfg, qc_only=True):
    summary_path = PAN_ROOT / "results/synteny_analysis/blocks" / f"{species}_synteny_summary.tsv"
    summary = pd.read_csv(summary_path, sep="\t") if summary_path.exists() else pd.DataFrame()
    qc_ids = set()
    if qc_only and not summary.empty:
        qc_ids = set(summary[qc_pass_mask(summary, cfg)]["sample_id"].astype(str))

    fastas = resolve_sample_fastas(species, sample_dir)
    rows = []
    sp_paf = Path(paf_root) / species
    for sid in sorted(fastas.keys()):
        if qc_ids and sid not in qc_ids:
            continue
        paf = sp_paf / f"{sid}.paf"
        if not paf.exists():
            continue
        segs, _ = load_paf_segments(
            paf, min_mapq=min_mapq, min_alen=int(anomaly_cfg.get("min_segment_len", 500))
        )
        hits = detect_local_anomalies(
            segs,
            max_query_gap=int(anomaly_cfg.get("max_query_gap", cfg.get("block", {}).get("max_gap", 200))),
            overlap_tol=int(cfg.get("block", {}).get("overlap_tol", 50)),
            min_segment_len=int(anomaly_cfg.get("min_segment_len", 500)),
            report_ref_decrease=bool(anomaly_cfg.get("report_ref_decrease", True)),
        )
        for h in hits:
            rows.append({"species": species, "sample_id": sid, **h})
    return rows


def main():
    parser = argparse.ArgumentParser(description="PAF local collinearity anomaly scan")
    parser.add_argument("--config", default=str(DEFAULT_CONFIG))
    parser.add_argument("--species", nargs="*", default=None)
    parser.add_argument("--all-species", action="store_true")
    parser.add_argument("--include-fail-qc", action="store_true")
    parser.add_argument("--output-dir", default=str(PAN_ROOT / "results/synteny_analysis/local_anomalies"))
    args = parser.parse_args()

    cfg = load_config(args.config)
    mm = cfg.get("minimap2", {})
    anomaly_cfg = cfg.get("paf_anomaly", {})
    min_mapq = int(mm.get("min_mapq", 1))
    min_alen = int(anomaly_cfg.get("min_segment_len", 500))
    paf_root = species_paf_dir(cfg)
    sample_dir = str(DATA_ROOT / "00-rawfa/single_IR")
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    species_list = SPECIES_LIST if args.all_species or not args.species else args.species
    all_rows = []
    for species in species_list:
        print(f"{species}...", flush=True)
        rows = process_species(
            species, cfg, paf_root, sample_dir, min_mapq, min_alen, anomaly_cfg,
            qc_only=not args.include_fail_qc,
        )
        df = pd.DataFrame(rows)
        if df.empty:
            df = empty_anomaly_frame()
        df.to_csv(out_dir / f"{species}_local_anomalies.tsv", sep="\t", index=False)
        print(f"  {len(df)} anomalies", flush=True)
        all_rows.extend(rows)

    all_df = pd.DataFrame(all_rows)
    if all_df.empty:
        all_df = empty_anomaly_frame()
    all_path = out_dir / "all_local_anomalies.tsv"
    all_df.to_csv(all_path, sep="\t", index=False)
    print(f"Total anomalies: {len(all_df)} -> {all_path}", flush=True)


if __name__ == "__main__":
    main()
