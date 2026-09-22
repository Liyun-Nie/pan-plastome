#!/usr/bin/env python3
"""PAF bp dotplot for a single sample (query vs reference coordinates)."""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from paf_viz_shared import (
    DEFAULT_CONFIG,
    figures_dir,
    load_config,
    load_paf_segments,
    species_paf_dir,
)


def plot_paf_dotplot(segments, species, sample_id, out_path, ref_len=None):
    if not segments:
        raise ValueError(f"No PAF segments for {species}/{sample_id}")

    fig, ax = plt.subplots(figsize=(7, 6))
    for seg in segments:
        qx = (seg["qstart"] + seg["qend"]) / 2 / 1000.0
        ry = (seg["tstart"] + seg["tend"]) / 2 / 1000.0
        color = "#2166ac" if seg["strand"] == "+" else "#b2182b"
        ax.plot(
            [seg["qstart"] / 1000.0, seg["qend"] / 1000.0],
            [seg["tstart"] / 1000.0, seg["tend"] / 1000.0],
            color=color, alpha=0.35, lw=1.0,
        )
        ax.scatter([qx], [ry], c=color, s=6, alpha=0.7, edgecolors="none")

    lim = max(
        max(s["qend"] for s in segments),
        max(max(s["tstart"], s["tend"]) for s in segments),
    ) / 1000.0 * 1.02
    ax.plot([0, lim], [0, lim], "--", color="gray", lw=0.8, alpha=0.5)
    ax.set_xlabel("Query position (kb)")
    ax.set_ylabel("Reference position (kb)")
    ax.set_title(
        f"{species.replace('_', ' ')} — {sample_id}\n"
        f"PAF segments: {len(segments)}",
        fontsize=10,
    )
    ax.set_aspect("equal", adjustable="box")
    fig.tight_layout()
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150)
    fig.savefig(out_path.with_suffix(".png"), dpi=150)
    plt.close(fig)
    print(f"Saved {out_path}")


def main():
    parser = argparse.ArgumentParser(description="PAF bp dotplot for one sample")
    parser.add_argument("--config", default=str(DEFAULT_CONFIG))
    parser.add_argument("--species", required=True)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--output", default=None)
    parser.add_argument("--representatives", action="store_true",
                        help="Plot all entries in config/representative_samples.tsv")
    args = parser.parse_args()

    cfg = load_config(args.config)
    mm = cfg.get("minimap2", {})
    blk = cfg.get("block", {})
    min_mapq = int(mm.get("min_mapq", 1))
    min_alen = int(blk.get("min_len", 500))
    paf_root = species_paf_dir(cfg)

    if args.representatives:
        from paf_viz_shared import load_representative_list
        reps = load_representative_list()
        for _, row in reps.iterrows():
            sp = row["species"]
            sid = str(row["sample_id"])
            paf = paf_root / sp / f"{sid}.paf"
            if not paf.exists():
                print(f"SKIP missing PAF: {sp}/{sid}")
                continue
            segs, _ = load_paf_segments(paf, min_mapq, min_alen)
            out = figures_dir(cfg, sp) / f"paf_dotplot_{sid}.pdf"
            plot_paf_dotplot(segs, sp, sid, out)
        return

    paf = paf_root / args.species / f"{args.sample_id}.paf"
    segs, _ = load_paf_segments(paf, min_mapq, min_alen)
    out = Path(args.output) if args.output else (
        figures_dir(cfg, args.species) / f"paf_dotplot_{args.sample_id}.pdf"
    )
    plot_paf_dotplot(segs, args.species, args.sample_id, out)


if __name__ == "__main__":
    main()
