#!/usr/bin/env python3
"""Species-level collinearity envelope from cached PAF (n > 500 samples)."""

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


def build_coverage_envelope(sample_ids, species, paf_root, ref_len, bin_bp, min_mapq, min_alen):
    if bin_bp <= 0:
        raise ValueError(f"bin_bp must be positive, got {bin_bp}")

    covered = None
    n_bins = None
    n_used = 0
    ref_len_out = ref_len

    for sid in sample_ids:
        paf = Path(paf_root) / species / f"{sid}.paf"
        if not paf.exists():
            continue
        segs, rl = load_paf_segments(paf, min_mapq, min_alen)
        if not segs:
            continue
        if ref_len_out <= 0:
            ref_len_out = rl
            n_bins = int(np.ceil(ref_len_out / bin_bp))
            covered = np.zeros(n_bins, dtype=np.int32)
        n_used += 1
        hit = np.zeros(n_bins, dtype=bool)
        for seg in segs:
            t0 = min(seg["tstart"], seg["tend"])
            t1 = max(seg["tstart"], seg["tend"])
            b0 = max(0, int(t0 // bin_bp))
            b1 = min(n_bins, int(np.ceil(t1 / bin_bp)))
            hit[b0:b1] = True
        covered += hit.astype(np.int32)

    if covered is None or n_used == 0 or n_bins is None:
        return np.array([]), np.array([]), 0, ref_len_out if ref_len_out > 0 else 0

    frac = covered / max(1, n_used)
    x_kb = (np.arange(n_bins) + 0.5) * bin_bp / 1000.0
    return x_kb, frac, n_used, ref_len_out


def plot_envelope(species, x_kb, frac, n_samples, ref_len, out_path):
    fig, ax = plt.subplots(figsize=(12, 4))
    ax.fill_between(x_kb, 0, frac, color="#2166ac", alpha=0.35, step="post")
    ax.plot(x_kb, frac, color="#2166ac", lw=0.8)
    ax.set_xlabel("Reference position (kb)")
    ax.set_ylabel("Fraction of samples with alignment")
    ax.set_ylim(0, 1.05)
    ax.set_xlim(0, ref_len / 1000.0)
    ax.set_title(
        f"{species.replace('_', ' ')} — collinearity envelope (n={n_samples} QC-pass)",
        fontsize=11,
    )
    fig.tight_layout()
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150)
    fig.savefig(out_path.with_suffix(".png"), dpi=150)
    plt.close(fig)
    print(f"Saved {out_path} ({n_samples} samples)")


def main():
    parser = argparse.ArgumentParser(description="PAF collinearity envelope per species")
    parser.add_argument("--config", default=str(DEFAULT_CONFIG))
    parser.add_argument("--species", nargs="*", default=None)
    parser.add_argument("--all-species", action="store_true")
    parser.add_argument("--bin-bp", type=int, default=None)
    parser.add_argument("--include-fail-qc", action="store_true")
    args = parser.parse_args()

    cfg = load_config(args.config)
    mm = cfg.get("minimap2", {})
    blk = cfg.get("block", {})
    dot = cfg.get("dotplot", {})
    min_mapq = int(mm.get("min_mapq", 1))
    min_alen = int(blk.get("min_len", 500))
    bin_bp = int(args.bin_bp or dot.get("raster_ref_bin_bp", 500))
    paf_root = species_paf_dir(cfg)

    species_list = SPECIES_LIST if args.all_species or not args.species else args.species

    for species in species_list:
        summary = load_species_summary(species)
        if not args.include_fail_qc:
            summary = summary[qc_pass_mask(summary, cfg)]
        n = len(summary)
        if choose_viz_mode(n) != "envelope":
            print(f"SKIP {species}: n={n} <= {RASTER_MAX_SAMPLES} (use raster)")
            continue
        sids = summary["sample_id"].astype(str).tolist()
        x_kb, frac, n_used, ref_len = build_coverage_envelope(
            sids, species, paf_root, 0, bin_bp, min_mapq, min_alen
        )
        if n_used == 0 or len(frac) == 0:
            print(f"SKIP {species}: no envelope data")
            continue
        out = figures_dir(cfg, species) / f"{species}_collinearity_envelope.pdf"
        plot_envelope(species, x_kb, frac, n_used, ref_len, out)


if __name__ == "__main__":
    main()
