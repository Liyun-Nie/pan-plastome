#!/usr/bin/env python3
"""Aggregate L1 block summaries and cross-check L2 window flags."""

import argparse
from pathlib import Path

import pandas as pd

from synteny_shared import DEFAULT_CONFIG, PAN_ROOT, canonical_sample_id, load_config


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default=str(DEFAULT_CONFIG))
    parser.add_argument("--blocks-dir", default=str(PAN_ROOT / "results/synteny_analysis/blocks"))
    parser.add_argument("--output-dir", default=str(PAN_ROOT / "results/synteny_analysis"))
    args = parser.parse_args()

    cfg = load_config(args.config)
    out_dir = Path(args.output_dir)
    blocks_dir = Path(args.blocks_dir)

    summaries = []
    events = []
    for spath in sorted(blocks_dir.glob("*_synteny_summary.tsv")):
        summaries.append(pd.read_csv(spath, sep="\t"))
    for epath in sorted(blocks_dir.glob("*_inter_block_events.tsv")):
        df = pd.read_csv(epath, sep="\t")
        if len(df) > 0:
            events.append(df)

    syn = pd.concat(summaries, ignore_index=True)
    syn.to_csv(out_dir / "all_samples_synteny_summary.tsv", sep="\t", index=False)

    species_agg = syn.groupby("species").agg(
        n_samples=("sample_id", "count"),
        mean_aligned_fraction=("aligned_fraction", "mean"),
        median_aligned_fraction=("aligned_fraction", "median"),
        n_synteny_events=("has_synteny_event", "sum"),
        mean_n_blocks=("n_blocks", "mean"),
    ).reset_index()
    species_agg.to_csv(out_dir / "species_synteny_aggregate.tsv", sep="\t", index=False)

    if events:
        ev = pd.concat(events, ignore_index=True)
        ev.to_csv(out_dir / "all_inter_block_events.tsv", sep="\t", index=False)
    else:
        pd.DataFrame(columns=["species", "sample_id", "event_type"]).to_csv(
            out_dir / "all_inter_block_events.tsv", sep="\t", index=False
        )

    # Cross-check L2 window monotonicity flags (rearrangement_details.tsv per species)
    win200_dir = PAN_ROOT / cfg["results_dirs"]["win200"]
    flag_frames = []
    for fpath in sorted(win200_dir.glob("*_rearrangement_details.tsv")):
        df = pd.read_csv(fpath, sep="\t")
        if len(df) == 0:
            continue
        species = fpath.name.replace("_rearrangement_details.tsv", "")
        df["species"] = species
        flag_frames.append(df)
    if flag_frames:
        flags = pd.concat(flag_frames, ignore_index=True)
        flags["sample_key"] = flags.apply(
            lambda r: (r["species"], canonical_sample_id(str(r["sample_id"]), r["species"])),
            axis=1,
        )
        syn["sample_key"] = syn.apply(
            lambda r: (r["species"], canonical_sample_id(str(r["sample_id"]), r["species"])),
            axis=1,
        )
        merged = flags.merge(
            syn[["sample_key", "aligned_fraction", "n_blocks", "has_synteny_event", "n_inversions", "n_translocations"]],
            on="sample_key",
            how="left",
        )
        bridge_dir = PAN_ROOT / "results/_archive/l2_bridge"
        bridge_dir.mkdir(parents=True, exist_ok=True)
        merged.to_csv(
            bridge_dir / "l2_flag_l1_crosscheck.tsv",
            sep="\t", index=False,
        )
        print(f"Win200 flagged samples: {len(flags)}")

    print(f"Samples: {len(syn)}")
    print(f"L1 synteny events: {syn['has_synteny_event'].sum()}")
    print(f"Species aggregate: {out_dir / 'species_synteny_aggregate.tsv'}")


if __name__ == "__main__":
    main()
