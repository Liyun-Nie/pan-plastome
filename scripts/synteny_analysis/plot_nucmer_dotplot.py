#!/usr/bin/env python3
"""Dotplot from archived nucmer delta (show-coords) for validation panel samples."""

import argparse
import subprocess
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from paf_viz_shared import DEFAULT_CONFIG, PAN_ROOT, figures_dir, load_config
from syri_representative_screen import parse_show_coords


def clusters_to_segments(clusters):
    segs = []
    for c in clusters:
        d = c["direction"].upper()
        strand = "-" if d.startswith(("R", "-", "N")) or "REVERSE" in d else "+"
        segs.append({
            "qstart": min(c["qry_start"], c["qry_end"]),
            "qend": max(c["qry_start"], c["qry_end"]),
            "tstart": min(c["ref_start"], c["ref_end"]),
            "tend": max(c["ref_start"], c["ref_end"]),
            "strand": strand,
        })
    return segs


def plot_nucmer_dotplot(segments, species, sample_id, out_path, verdict=""):
    fig, ax = plt.subplots(figsize=(7, 6))
    for seg in segments:
        q0, q1 = seg["qstart"] / 1000.0, seg["qend"] / 1000.0
        t0, t1 = seg["tstart"] / 1000.0, seg["tend"] / 1000.0
        color = "#2166ac" if seg["strand"] == "+" else "#b2182b"
        ax.plot([q0, q1], [t0, t1], color=color, alpha=0.5, lw=1.2)
        ax.scatter([(q0 + q1) / 2], [(t0 + t1) / 2], c=color, s=8, alpha=0.8, edgecolors="none")
    if segments:
        lim = max(
            max(s["qend"] for s in segments),
            max(max(s["tstart"], s["tend"]) for s in segments),
        ) / 1000.0 * 1.02
        ax.plot([0, lim], [0, lim], "--", color="gray", lw=0.8, alpha=0.5)
    ax.set_xlabel("Query position (kb)")
    ax.set_ylabel("Reference position (kb)")
    title = f"{species.replace('_', ' ')} — {sample_id}\nnucmer clusters: {len(segments)}"
    if verdict:
        title += f" ({verdict})"
    ax.set_title(title, fontsize=10)
    ax.set_aspect("equal", adjustable="box")
    fig.tight_layout()
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150)
    fig.savefig(out_path.with_suffix(".png"), dpi=150)
    plt.close(fig)
    print(f"Saved {out_path}")


def main():
    parser = argparse.ArgumentParser(description="nucmer delta dotplot from archived work dir")
    parser.add_argument("--config", default=str(DEFAULT_CONFIG))
    parser.add_argument("--delta", required=True, help="Path to out.delta")
    parser.add_argument("--species", required=True)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--verdict", default="")
    parser.add_argument("--output", default=None)
    args = parser.parse_args()
    cfg = load_config(args.config)

    delta = Path(args.delta)
    if not delta.exists():
        raise FileNotFoundError(delta)
    clusters = parse_show_coords(delta)
    segs = clusters_to_segments(clusters)
    out = Path(args.output) if args.output else (
        figures_dir(cfg, args.species) / f"nucmer_dotplot_{args.sample_id}.pdf"
    )
    plot_nucmer_dotplot(segs, args.species, args.sample_id, out, args.verdict)


if __name__ == "__main__":
    main()
