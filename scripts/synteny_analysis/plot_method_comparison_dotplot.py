#!/usr/bin/env python3
"""Side-by-side minimap2 (PAF) vs BLAST (merged blocks) dotplots for discordant representatives."""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from paf_viz_shared import (
    DEFAULT_CONFIG,
    PAN_ROOT,
    figures_dir,
    load_config,
    load_paf_segments,
    species_paf_dir,
)


def _draw_segments(ax, segments, title, kb=True):
    scale = 1000.0 if kb else 1.0
    for seg in segments:
        q0, q1 = seg["qstart"] / scale, seg["qend"] / scale
        t0, t1 = seg["tstart"] / scale, seg["tend"] / scale
        color = "#2166ac" if seg.get("strand", "+") == "+" else "#b2182b"
        ax.plot([q0, q1], [t0, t1], color=color, alpha=0.5, lw=1.2)
        ax.scatter([(q0 + q1) / 2], [(t0 + t1) / 2], c=color, s=8, alpha=0.8, edgecolors="none")
    if segments:
        lim = max(
            max(s["qend"] for s in segments),
            max(max(s["tstart"], s["tend"]) for s in segments),
        ) / scale * 1.02
        ax.plot([0, lim], [0, lim], "--", color="gray", lw=0.8, alpha=0.5)
    ax.set_xlabel("Query position (kb)")
    ax.set_ylabel("Reference position (kb)")
    ax.set_title(title, fontsize=9)
    ax.set_aspect("equal", adjustable="box")


def load_blocks_tsv(path, sample_id):
    df = pd.read_csv(path, sep="\t")
    sub = df[df["sample_id"].astype(str) == str(sample_id)]
    rows = []
    for _, r in sub.iterrows():
        rows.append({
            "qstart": int(r["qstart"]),
            "qend": int(r["qend"]),
            "tstart": int(r["tstart"]),
            "tend": int(r["tend"]),
            "strand": str(r["strand"]),
        })
    return rows


def plot_comparison(species, sample_id, cfg, blast_blocks_dir, out_path):
    mm = cfg.get("minimap2", {})
    blk = cfg.get("block", {})
    min_mapq = int(mm.get("min_mapq", 1))
    min_alen = int(blk.get("min_len", 500))

    paf = species_paf_dir(cfg) / species / f"{sample_id}.paf"
    paf_segs, _ = load_paf_segments(paf, min_mapq, min_alen)
    mm_blocks = load_blocks_tsv(
        PAN_ROOT / "results/synteny_analysis/blocks" / f"{species}_synteny_blocks.tsv",
        sample_id,
    )
    blast_blocks = load_blocks_tsv(
        Path(blast_blocks_dir) / f"{species}_synteny_blocks.tsv",
        sample_id,
    )

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    _draw_segments(
        axes[0], paf_segs,
        f"minimap2 PAF segments (n={len(paf_segs)})",
    )
    _draw_segments(
        axes[1], mm_blocks,
        f"minimap2 merged blocks (n={len(mm_blocks)})",
    )
    _draw_segments(
        axes[2], blast_blocks,
        f"BLAST merged blocks (n={len(blast_blocks)})",
    )
    fig.suptitle(
        f"{species.replace('_', ' ')} — {sample_id}\n"
        f"Method discordance: minimap2 vs BLAST cross-validation",
        fontsize=11,
    )
    fig.tight_layout()
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150)
    fig.savefig(out_path.with_suffix(".png"), dpi=150)
    plt.close(fig)
    print(f"Saved {out_path}")


def discordant_pairs(method_agreement_path):
    df = pd.read_csv(method_agreement_path, sep="\t")
    disc = df[df["agreement"] == "discordant"]
    return list(zip(disc["species"], disc["sample_id"].astype(str)))


def main():
    parser = argparse.ArgumentParser(description="minimap2 vs BLAST comparison dotplots")
    parser.add_argument("--config", default=str(DEFAULT_CONFIG))
    parser.add_argument("--species", default=None)
    parser.add_argument("--sample-id", default=None)
    parser.add_argument("--discordant-all", action="store_true")
    parser.add_argument(
        "--blast-blocks-dir",
        default=str(PAN_ROOT / "results/synteny_analysis/validation/blast_blocks"),
    )
    parser.add_argument(
        "--method-agreement",
        default=str(PAN_ROOT / "results/synteny_analysis/validation/method_agreement.tsv"),
    )
    args = parser.parse_args()
    cfg = load_config(args.config)

    if args.discordant_all:
        pairs = discordant_pairs(args.method_agreement)
    elif args.species and args.sample_id:
        pairs = [(args.species, args.sample_id)]
    else:
        pairs = discordant_pairs(args.method_agreement)

    for species, sid in pairs:
        out = figures_dir(cfg, species) / f"method_comparison_dotplot_{sid}.pdf"
        plot_comparison(species, sid, cfg, args.blast_blocks_dir, out)


if __name__ == "__main__":
    main()
