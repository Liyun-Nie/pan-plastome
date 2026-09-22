#!/usr/bin/env python3
"""Bar chart: QC-pass vs total sample counts per species (all 17 species)."""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from paf_viz_shared import (
    DEFAULT_CONFIG,
    PAN_ROOT,
    figures_dir,
    load_config,
    load_species_summary,
    manifest_relpath,
    qc_pass_mask,
)
from window_rearrangement_test import SPECIES_LIST


def build_qc_summary(cfg, species_list=None):
    species_list = species_list or SPECIES_LIST
    qc_cfg = cfg.get("qc", {})
    min_af = float(qc_cfg.get("min_aligned_fraction", 0.95))
    ratio = qc_cfg.get("query_ref_len_ratio", [0.95, 1.05])
    lo, hi = float(ratio[0]), float(ratio[1])

    rows = []
    for species in species_list:
        summary = load_species_summary(species)
        mask = qc_pass_mask(summary, cfg)
        n_total = len(summary)
        n_pass = int(mask.sum())
        rows.append({
            "species": species,
            "species_label": species.replace("_", " "),
            "n_total": n_total,
            "n_qc_pass": n_pass,
            "n_qc_fail": n_total - n_pass,
            "qc_pass_fraction": n_pass / n_total if n_total else 0.0,
        })

    df = pd.DataFrame(rows).sort_values("n_total", ascending=True)
    return df, min_af, lo, hi


def plot_qc_bars(df, min_af, len_lo, len_hi, out_path, show_pct=True):
    labels = df["species_label"].tolist()
    n_total = df["n_total"].to_numpy()
    n_pass = df["n_qc_pass"].to_numpy()
    n_fail = df["n_qc_fail"].to_numpy()
    y = np.arange(len(labels))

    fig_h = max(5, 0.35 * len(labels) + 1.5)
    fig, ax = plt.subplots(figsize=(10, fig_h))

    bar_h = 0.72
    ax.barh(y, n_total, height=bar_h, color="#d9d9d9", edgecolor="#888888", lw=0.6,
            label="Total samples", zorder=1)
    ax.barh(y, n_pass, height=bar_h * 0.55, color="#2166ac", edgecolor="none",
            label="QC-pass samples", zorder=2)

    xmax = max(n_total.max() * 1.18, 100)
    for yi, (tot, pas, fail) in enumerate(zip(n_total, n_pass, n_fail)):
        pct = 100.0 * pas / tot if tot else 0.0
        text = f"{pas:,} / {tot:,} ({pct:.1f}%)" if show_pct else f"{pas:,} / {tot:,}"
        ax.text(tot + xmax * 0.01, yi, text, va="center", ha="left", fontsize=8)

    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=9)
    ax.set_xlabel("Sample count")
    ax.set_xlim(0, xmax)
    ax.set_title(
        "Collinearity QC: sample counts by species\n"
        f"(QC-pass: aligned_fraction ≥ {min_af:.0%}, "
        f"query/ref length ∈ [{len_lo:.2f}, {len_hi:.2f}])",
        fontsize=11,
    )
    ax.legend(loc="lower right", framealpha=0.9)
    ax.grid(axis="x", alpha=0.25, lw=0.6)
    ax.set_axisbelow(True)

    tot_all = int(n_total.sum())
    pass_all = int(n_pass.sum())
    fig.text(
        0.01, 0.01,
        f"Panel total: {pass_all:,} / {tot_all:,} QC-pass ({100 * pass_all / tot_all:.2f}%)",
        fontsize=9, color="#333333",
    )

    fig.tight_layout(rect=[0, 0.03, 1, 1])
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150)
    fig.savefig(out_path.with_suffix(".png"), dpi=150)
    plt.close(fig)
    print(f"Saved {out_path}")


def main():
    parser = argparse.ArgumentParser(description="QC-pass vs total sample bar chart")
    parser.add_argument("--config", default=str(DEFAULT_CONFIG))
    parser.add_argument("--output", default=None)
    parser.add_argument("--tsv", default=None)
    args = parser.parse_args()

    cfg = load_config(args.config)
    df, min_af, lo, hi = build_qc_summary(cfg)

    fig_base = figures_dir(cfg) / "qc_sample_counts_by_species"
    out = Path(args.output) if args.output else fig_base.with_suffix(".pdf")
    tsv_out = Path(args.tsv) if args.tsv else (
        PAN_ROOT / "results/synteny_analysis/qc_sample_counts_by_species.tsv"
    )

    tsv_out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(tsv_out, sep="\t", index=False)
    print(f"Wrote {tsv_out}")

    plot_qc_bars(df, min_af, lo, hi, out)
    print(f"Manifest path: {manifest_relpath(out)}")


if __name__ == "__main__":
    main()
