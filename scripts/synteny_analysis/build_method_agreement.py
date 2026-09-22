#!/usr/bin/env python3
"""Build representative sample list and minimap2 vs blastn agreement table (Phase 2a-7)."""

import argparse
from pathlib import Path

import pandas as pd

from synteny_shared import DEFAULT_CONFIG, PAN_ROOT, load_config
from window_rearrangement_test import SPECIES_LIST

VALIDATION_DIR = PAN_ROOT / "results/synteny_analysis/validation"
REP_LIST = PAN_ROOT / "config/representative_samples.tsv"
SUMMARY_ALL = PAN_ROOT / "results/synteny_analysis/all_samples_synteny_summary.tsv"

# C2: ensure minus-strand HSP coverage in cross-check panel
FORCED_REPRESENTATIVES = {
    "Oryza_sativa": "I11311",
}


def pick_representatives(summary_path, out_path):
    syn = pd.read_csv(summary_path, sep="\t")
    rows = []
    for species in SPECIES_LIST:
        sub = syn[syn["species"] == species].copy()
        if sub.empty:
            continue
        if species in FORCED_REPRESENTATIVES:
            sid = FORCED_REPRESENTATIVES[species]
            pick = sub[sub["sample_id"].astype(str) == sid]
            if pick.empty:
                pick = sub.iloc[[(sub["aligned_fraction"] - sub["aligned_fraction"].median()).abs().argmin()]]
                selection = "median_fallback"
            else:
                selection = "forced_c2_minus_strand"
        else:
            med = sub["aligned_fraction"].median()
            idx = (sub["aligned_fraction"] - med).abs().idxmin()
            pick = sub.loc[[idx]]
            sid = str(pick.iloc[0]["sample_id"])
            selection = "median_aligned_fraction"
        row = pick.iloc[0]
        rows.append({
            "species": species,
            "sample_id": str(sid),
            "selection": selection,
            "aligned_fraction_minimap2": row["aligned_fraction"],
            "n_blocks_minimap2": int(row["n_blocks"]),
            "n_inversions_minimap2": int(row["n_inversions"]),
            "n_translocations_minimap2": int(row["n_translocations"]),
            "has_synteny_event_minimap2": bool(row["has_synteny_event"]),
        })
    rep = pd.DataFrame(rows)
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    rep[["species", "sample_id", "selection"]].to_csv(out_path, sep="\t", index=False)
    rep.to_csv(out_path.with_suffix(".full.tsv"), sep="\t", index=False)
    print(f"Wrote {len(rep)} representatives -> {out_path}")
    return rep


def build_agreement(blast_summary_dir, rep_full_path, out_path):
    rep = pd.read_csv(rep_full_path, sep="\t")
    rep["sample_id"] = rep["sample_id"].astype(str)
    blast_rows = []
    for species in rep["species"]:
        spath = Path(blast_summary_dir) / f"{species}_synteny_summary.tsv"
        if not spath.exists():
            raise FileNotFoundError(f"Missing BLAST summary: {spath}")
        sub = pd.read_csv(spath, sep="\t")
        blast_rows.append(sub)
    blast = pd.concat(blast_rows, ignore_index=True)
    blast["sample_id"] = blast["sample_id"].astype(str)
    blast = blast.rename(columns={
        "aligned_fraction": "aligned_fraction_blastn",
        "n_blocks": "n_blocks_blastn",
        "n_inversions": "n_inversions_blastn",
        "n_translocations": "n_translocations_blastn",
        "has_synteny_event": "has_synteny_event_blastn",
    })

    merged = rep.merge(
        blast[["species", "sample_id", "aligned_fraction_blastn", "n_blocks_blastn",
               "n_inversions_blastn", "n_translocations_blastn", "has_synteny_event_blastn"]],
        on=["species", "sample_id"],
        how="left",
    )
    # Normalize sample_id types (BLAST TSV may read numeric IDs as int)
    merged["sample_id"] = merged["sample_id"].astype(str)

    def agree(row):
        keys = (
            row["n_blocks_minimap2"] == row["n_blocks_blastn"],
            row["has_synteny_event_minimap2"] == row["has_synteny_event_blastn"],
            row["n_inversions_minimap2"] == row["n_inversions_blastn"],
            row["n_translocations_minimap2"] == row["n_translocations_blastn"],
        )
        return "concordant" if all(keys) else "discordant"

    merged["agreement"] = merged.apply(agree, axis=1)
    merged["af_delta"] = (merged["aligned_fraction_blastn"] - merged["aligned_fraction_minimap2"]).round(6)

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(out_path, sep="\t", index=False)

    n_con = (merged["agreement"] == "concordant").sum()
    print(f"Agreement: {n_con}/{len(merged)} concordant")
    if (merged["agreement"] == "discordant").any():
        print("Discordant samples:")
        d = merged[merged["agreement"] == "discordant"]
        print(d[["species", "sample_id", "n_blocks_minimap2", "n_blocks_blastn",
                 "has_synteny_event_minimap2", "has_synteny_event_blastn"]].to_string(index=False))
    print(f"Wrote {out_path}")
    return merged


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pick-only", action="store_true", help="Only build representative_samples.tsv")
    parser.add_argument("--agreement-only", action="store_true", help="Only build agreement from existing BLAST")
    parser.add_argument("--summary", default=str(SUMMARY_ALL))
    parser.add_argument("--rep-list", default=str(REP_LIST))
    parser.add_argument("--blast-blocks-dir", default=str(VALIDATION_DIR / "blast_blocks"))
    parser.add_argument("--out", default=str(VALIDATION_DIR / "method_agreement.tsv"))
    args = parser.parse_args()

    rep_full = Path(args.rep_list).with_suffix(".full.tsv")
    if not args.agreement_only:
        pick_representatives(args.summary, args.rep_list)
    if not args.pick_only:
        build_agreement(args.blast_blocks_dir, rep_full, args.out)


if __name__ == "__main__":
    main()
