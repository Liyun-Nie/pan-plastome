#!/usr/bin/env python3
"""Species-level PAF raster: all QC-pass samples vs reference bp (n <= 500)."""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from paf_viz_shared import (
    DEFAULT_CONFIG,
    RASTER_MAX_SAMPLES,
    choose_viz_mode,
    figures_dir,
    load_config,
    load_paf_segments,
    load_species_summary,
    qc_pass_mask,
    species_paf_dir,
)
from window_rearrangement_test import SPECIES_LIST


def plot_species_raster(species, sample_ids, paf_root, cfg, out_path, min_mapq=1, min_alen=500):
    ref_len = None
    rows = []
    for sid in sample_ids:
        paf = Path(paf_root) / species / f"{sid}.paf"
        if not paf.exists():
            continue
        segs, rl = load_paf_segments(paf, min_mapq, min_alen)
        if ref_len is None:
            ref_len = rl
        rows.append((sid, segs))

    if not rows or not ref_len:
        print(f"SKIP {species}: no raster data")
        return False

    n = len(rows)
    fig_h = max(4, min(20, 0.012 * n + 2))
    fig, ax = plt.subplots(figsize=(12, fig_h))

    for yi, (sid, segs) in enumerate(rows):
        for seg in segs:
            if seg["strand"] == "+":
                t0, t1 = seg["tstart"], seg["tend"]
            else:
                t0, t1 = seg["tend"], seg["tstart"]
            ax.plot([t0 / 1000, t1 / 1000], [yi, yi], color="#2166ac", lw=0.8, alpha=0.7)

    ax.set_xlabel("Reference position (kb)")
    ax.set_ylabel(f"Sample index (n={n} QC-pass)")
    ax.set_title(f"{species.replace('_', ' ')} — PAF collinearity raster", fontsize=11)
    ax.set_xlim(0, ref_len / 1000.0)
    ax.set_ylim(-1, n)
    fig.tight_layout()
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150)
    fig.savefig(out_path.with_suffix(".png"), dpi=150)
    plt.close(fig)
    print(f"Saved {out_path} ({n} samples)")
    return True


def main():
    parser = argparse.ArgumentParser(description="PAF raster plot per species")
    parser.add_argument("--config", default=str(DEFAULT_CONFIG))
    parser.add_argument("--species", nargs="*", default=None)
    parser.add_argument("--all-species", action="store_true")
    parser.add_argument("--include-fail-qc", action="store_true")
    args = parser.parse_args()

    cfg = load_config(args.config)
    mm = cfg.get("minimap2", {})
    blk = cfg.get("block", {})
    min_mapq = int(mm.get("min_mapq", 1))
    min_alen = int(blk.get("min_len", 500))
    paf_root = species_paf_dir(cfg)

    species_list = args.species if args.species else SPECIES_LIST
    if args.all_species:
        species_list = SPECIES_LIST

    for species in species_list:
        summary = load_species_summary(species)
        if not args.include_fail_qc:
            summary = summary[qc_pass_mask(summary, cfg)]
        n = len(summary)
        if choose_viz_mode(n) != "raster":
            print(f"SKIP {species}: n={n} > {RASTER_MAX_SAMPLES} (use envelope)")
            continue
        sids = summary.sort_values("sample_id")["sample_id"].astype(str).tolist()
        out = figures_dir(cfg, species) / f"{species}_all_samples_raster.pdf"
        plot_species_raster(species, sids, paf_root, cfg, out, min_mapq, min_alen)


if __name__ == "__main__":
    main()
