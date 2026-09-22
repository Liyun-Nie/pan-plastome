#!/usr/bin/env python3
"""PAF dotplots at two min_alen thresholds (e.g. 500 bp primary vs 300 bp sensitivity)."""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from collinear_block_merge import merge_alignment_segments
from paf_viz_shared import (
    DEFAULT_CONFIG,
    figures_dir,
    load_config,
    load_paf_segments,
    manifest_relpath,
    species_paf_dir,
)
from synteny_block_analysis import detect_inter_block_events

DEFAULT_TRIO = ("hh311", "hh361", "hh391")
SPECIES = "Hemerocallis_citrina"


def _summarize_blocks(segments, max_gap, overlap_tol, trans_jump, q_len):
    blocks = merge_alignment_segments(segments, max_gap=max_gap, overlap_tol=overlap_tol)
    events = detect_inter_block_events(blocks, translocation_jump=trans_jump)
    n_inv = sum(1 for e in events if e["event_type"] == "inversion")
    n_trans = sum(1 for e in events if e["event_type"] == "translocation")
    aligned_bp = sum(b["block_len"] for b in blocks)
    af = round(aligned_bp / q_len, 6) if q_len else 0.0
    return {
        "n_blocks": len(blocks),
        "n_inversions": n_inv,
        "n_translocations": n_trans,
        "has_event": len(events) > 0,
        "aligned_fraction": af,
        "blocks": blocks,
    }


def _draw_panel(ax, segments, title_suffix):
    if not segments:
        ax.text(0.5, 0.5, "No segments", ha="center", va="center", transform=ax.transAxes)
        ax.set_title(title_suffix, fontsize=9)
        return

    for seg in segments:
        q0 = seg["qstart"] / 1000.0
        q1 = seg["qend"] / 1000.0
        t0 = seg["tstart"] / 1000.0
        t1 = seg["tend"] / 1000.0
        color = "#2166ac" if seg["strand"] == "+" else "#b2182b"
        ax.plot([q0, q1], [t0, t1], color=color, alpha=0.45, lw=1.2)
        ax.scatter([(q0 + q1) / 2], [(t0 + t1) / 2], c=color, s=10, alpha=0.85, edgecolors="none")

    lim = max(
        max(s["qend"] for s in segments),
        max(max(s["tstart"], s["tend"]) for s in segments),
    ) / 1000.0 * 1.02
    ax.plot([0, lim], [0, lim], "--", color="gray", lw=0.8, alpha=0.5)
    ax.set_xlabel("Query (kb)")
    ax.set_ylabel("Reference (kb)")
    ax.set_title(title_suffix, fontsize=9)
    ax.set_aspect("equal", adjustable="box")


def plot_single_threshold(segments, species, sample_id, min_alen, stats, out_path):
    event_txt = "event" if stats["has_event"] else "no event"
    if stats["n_inversions"]:
        event_txt = f"{stats['n_inversions']} inversion"
    fig, ax = plt.subplots(figsize=(7, 6))
    title = (
        f"{species.replace('_', ' ')} — {sample_id}\n"
        f"min_alen={min_alen} bp | segments={len(segments)} | "
        f"blocks={stats['n_blocks']} | {event_txt} | AF={stats['aligned_fraction']:.4f}"
    )
    _draw_panel(ax, segments, title)
    fig.tight_layout()
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150)
    fig.savefig(out_path.with_suffix(".png"), dpi=150)
    plt.close(fig)
    print(f"Saved {out_path}")


def plot_comparison(segments_lo, segments_hi, species, sample_id, min_lo, min_hi,
                    stats_lo, stats_hi, out_path):
    def _subtitle(min_alen, segs, st):
        ev = "no event"
        if st["n_inversions"]:
            ev = f"{st['n_inversions']} inversion"
        elif st["has_event"]:
            ev = "event"
        return (
            f"min_alen={min_alen} bp\n"
            f"segments={len(segs)}, blocks={st['n_blocks']}, {ev}\n"
            f"AF={st['aligned_fraction']:.4f}"
        )

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    _draw_panel(axes[0], segments_hi, _subtitle(min_hi, segments_hi, stats_hi))
    _draw_panel(axes[1], segments_lo, _subtitle(min_lo, segments_lo, stats_lo))
    fig.suptitle(
        f"{species.replace('_', ' ')} — {sample_id}\n"
        f"Primary ({min_hi} bp) vs sensitivity ({min_lo} bp) PAF dotplot",
        fontsize=11,
    )
    fig.tight_layout()
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150)
    fig.savefig(out_path.with_suffix(".png"), dpi=150)
    plt.close(fig)
    print(f"Saved {out_path}")


def process_sample(species, sample_id, cfg, min_lo, min_hi, paf_root, out_dir):
    mm = cfg.get("minimap2", {})
    blk = cfg.get("block", {})
    min_mapq = int(mm.get("min_mapq", 1))
    max_gap = int(blk.get("max_gap", 200))
    overlap_tol = int(blk.get("overlap_tol", 50))
    trans_jump = int(blk.get("translocation_jump", 10000))

    paf = Path(paf_root) / species / f"{sample_id}.paf"
    if not paf.exists():
        print(f"SKIP missing PAF: {paf}")
        return []

    segs_hi, _ = load_paf_segments(paf, min_mapq, min_hi)
    segs_lo, _ = load_paf_segments(paf, min_mapq, min_lo)
    q_len = max((s["qend"] for s in segs_lo), default=0)
    if not q_len and segs_hi:
        q_len = max(s["qend"] for s in segs_hi)

    stats_hi = _summarize_blocks(segs_hi, max_gap, overlap_tol, trans_jump, q_len)
    stats_lo = _summarize_blocks(segs_lo, max_gap, overlap_tol, trans_jump, q_len)

    paths = []
    cmp_out = out_dir / f"sensitivity_dotplot_{sample_id}_min{min_hi}_vs_{min_lo}.pdf"
    plot_comparison(segs_lo, segs_hi, species, sample_id, min_lo, min_hi,
                    stats_lo, stats_hi, cmp_out)
    paths.append(cmp_out)

    for min_alen, segs, stats in (
        (min_hi, segs_hi, stats_hi),
        (min_lo, segs_lo, stats_lo),
    ):
        single = out_dir / f"paf_dotplot_{sample_id}_minalen{min_alen}.pdf"
        plot_single_threshold(segs, species, sample_id, min_alen, stats, single)
        paths.append(single)

    return paths


def main():
    parser = argparse.ArgumentParser(description="PAF dotplots at 500 vs 300 bp min_alen")
    parser.add_argument("--config", default=str(DEFAULT_CONFIG))
    parser.add_argument("--species", default=SPECIES)
    parser.add_argument("--sample-id", action="append", default=None)
    parser.add_argument("--hemerocallis-trio", action="store_true",
                        help="Plot hh311, hh361, hh391 (default when no --sample-id)")
    parser.add_argument("--min-alen-lo", type=int, default=300)
    parser.add_argument("--min-alen-hi", type=int, default=500)
    args = parser.parse_args()

    cfg = load_config(args.config)
    paf_root = species_paf_dir(cfg)
    out_dir = figures_dir(cfg, args.species)

    if args.sample_id:
        sample_ids = args.sample_id
    else:
        sample_ids = list(DEFAULT_TRIO)

    all_paths = []
    for sid in sample_ids:
        all_paths.extend(process_sample(
            args.species, sid, cfg, args.min_alen_lo, args.min_alen_hi, paf_root, out_dir
        ))

    for p in all_paths:
        print(f"  manifest: {manifest_relpath(p)}")


if __name__ == "__main__":
    main()
