#!/usr/bin/env python3
"""Orchestrate Phase 1b PAF visualization for all 17 species."""

import argparse
from pathlib import Path

import pandas as pd

from paf_viz_shared import (
    DEFAULT_CONFIG,
    RASTER_MAX_SAMPLES,
    choose_viz_mode,
    figures_dir,
    load_config,
    load_representative_list,
    load_species_summary,
    manifest_relpath,
    qc_pass_mask,
)
from plot_paf_dotplot import plot_paf_dotplot
from plot_synteny_envelope import build_coverage_envelope, plot_envelope
from plot_synteny_raster import plot_species_raster
from paf_viz_shared import load_paf_segments, species_paf_dir
from window_rearrangement_test import SPECIES_LIST


def main():
    parser = argparse.ArgumentParser(description="Run Phase 1b viz batch")
    parser.add_argument("--config", default=str(DEFAULT_CONFIG))
    parser.add_argument("--species", nargs="*", default=None)
    parser.add_argument("--skip-dotplots", action="store_true")
    args = parser.parse_args()

    cfg = load_config(args.config)
    mm = cfg.get("minimap2", {})
    blk = cfg.get("block", {})
    dot = cfg.get("dotplot", {})
    min_mapq = int(mm.get("min_mapq", 1))
    min_alen = int(blk.get("min_len", 500))
    bin_bp = int(dot.get("raster_ref_bin_bp", 500))
    paf_root = species_paf_dir(cfg)

    species_list = args.species if args.species else SPECIES_LIST
    manifest = []

    if not args.skip_dotplots:
        reps = load_representative_list()
        for _, row in reps.iterrows():
            sp, sid = row["species"], str(row["sample_id"])
            paf = paf_root / sp / f"{sid}.paf"
            if not paf.exists():
                continue
            segs, _ = load_paf_segments(paf, min_mapq, min_alen)
            out = figures_dir(cfg, sp) / f"paf_dotplot_{sid}.pdf"
            plot_paf_dotplot(segs, sp, sid, out)
            manifest.append({
                "species": sp,
                "figure_type": "paf_dotplot",
                "path": manifest_relpath(out),
            })

    for species in species_list:
        summary = load_species_summary(species)
        qc = summary[qc_pass_mask(summary, cfg)]
        n = len(qc)
        mode = choose_viz_mode(n)
        sids = qc.sort_values("sample_id")["sample_id"].astype(str).tolist()

        if mode == "raster":
            out = figures_dir(cfg, species) / f"{species}_all_samples_raster.pdf"
            ok = plot_species_raster(
                species, sids, paf_root, cfg, out, min_mapq, min_alen
            )
            if ok:
                manifest.append({
                    "species": species,
                    "figure_type": "raster",
                    "path": manifest_relpath(out),
                    "n_qc_pass": n,
                })
        else:
            x_kb, frac, n_used, ref_len = build_coverage_envelope(
                sids, species, paf_root, 0, bin_bp, min_mapq, min_alen
            )
            out = figures_dir(cfg, species) / f"{species}_collinearity_envelope.pdf"
            plot_envelope(species, x_kb, frac, n_used, ref_len, out)
            manifest.append({
                "species": species,
                "figure_type": "envelope",
                "path": manifest_relpath(out),
                "n_qc_pass": n_used,
            })

    out_manifest = figures_dir(cfg) / "phase1b_figure_manifest.tsv"
    pd.DataFrame(manifest).to_csv(out_manifest, sep="\t", index=False)
    print(f"Manifest: {out_manifest} ({len(manifest)} figures)")


if __name__ == "__main__":
    main()
