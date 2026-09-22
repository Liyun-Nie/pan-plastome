#!/usr/bin/env python3
"""Reproducible audit of results/synteny_analysis/ for review reports."""

import argparse
import json
from pathlib import Path

import pandas as pd

from paf_viz_shared import resolve_manifest_path

SCRIPT_DIR = Path(__file__).resolve().parent
PAN_ROOT = SCRIPT_DIR.parent.parent
SA = PAN_ROOT / "results" / "synteny_analysis"
ARCH = PAN_ROOT / "results" / "_archive"


def read_tsv(p):
    if not p.exists() or p.stat().st_size <= 1:
        return None
    try:
        return pd.read_csv(p, sep="\t")
    except pd.errors.EmptyDataError:
        return None


def run_audit():
    out = {"issues": [], "warnings": [], "blocks": {}, "local_anomalies": {}}

    species_summaries = sorted((SA / "blocks").glob("*_synteny_summary.tsv"))
    all_summ = [read_tsv(p) for p in species_summaries]
    all_summ = [d for d in all_summ if d is not None]
    if all_summ:
        syn = pd.concat(all_summ, ignore_index=True)
        out["blocks"] = {
            "n_samples": int(len(syn)),
            "n_species": int(syn["species"].nunique()),
            "has_synteny_event_sum": int(syn["has_synteny_event"].sum()),
            "mean_aligned_fraction": float(syn["aligned_fraction"].mean()),
            "min_aligned_fraction": float(syn["aligned_fraction"].min()),
            "n_blocks_gt1": int((syn["n_blocks"] > 1).sum()),
            "multi_block_samples": syn[syn["n_blocks"] > 1][
                ["species", "sample_id", "n_blocks", "aligned_fraction", "has_synteny_event"]
            ].to_dict(orient="records"),
        }

    as_df = read_tsv(SA / "all_samples_synteny_summary.tsv")
    if as_df is not None and out.get("blocks"):
        if len(as_df) != out["blocks"]["n_samples"]:
            out["issues"].append(
                f"blocks sum ({out['blocks']['n_samples']}) != all_samples ({len(as_df)})"
            )

    out["paf_files"] = len(list((SA / "paf").rglob("*.paf")))
    if out.get("blocks") and out["paf_files"] != out["blocks"]["n_samples"]:
        out["warnings"].append(
            f"PAF count ({out['paf_files']}) != sample count ({out['blocks']['n_samples']})"
        )

    la = read_tsv(SA / "local_anomalies" / "all_local_anomalies.tsv")
    out["local_anomalies"]["rows"] = int(len(la)) if la is not None else 0
    if (SA / "local_anomalies" / "all_local_anomalies.tsv").exists():
        la_path = SA / "local_anomalies" / "all_local_anomalies.tsv"
        if la_path.stat().st_size <= 1:
            out["issues"].append("all_local_anomalies.tsv is empty (no header)")

    fm = read_tsv(SA / "figures" / "phase1b_figure_manifest.tsv")
    if fm is not None:
        missing = [
            p for p in fm["path"]
            if not resolve_manifest_path(p, PAN_ROOT).exists()
        ]
        if missing:
            out["issues"].append(f"figure manifest missing {len(missing)} files")
        out["figure_manifest"] = {
            "total": int(len(fm)),
            "present": int(len(fm) - len(missing)),
        }

    return out


def main():
    parser = argparse.ArgumentParser(description="Audit synteny_analysis results")
    parser.add_argument("--json", action="store_true", help="Print JSON to stdout")
    args = parser.parse_args()
    result = run_audit()
    if args.json:
        print(json.dumps(result, indent=2, ensure_ascii=False, default=str))
    else:
        print(f"Samples: {result.get('blocks', {}).get('n_samples', 'NA')}")
        print(f"Inter-block events (primary): {result.get('blocks', {}).get('has_synteny_event_sum', 'NA')}")
        print(f"Issues: {result.get('issues')}")
        print(f"Warnings: {result.get('warnings')}")
        fm = result.get("figure_manifest")
        if fm:
            print(f"Figure manifest: {fm.get('present')}/{fm.get('total')} present")


if __name__ == "__main__":
    main()
