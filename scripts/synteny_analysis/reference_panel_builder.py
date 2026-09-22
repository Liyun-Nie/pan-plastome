#!/usr/bin/env python3
"""
Build multi-reference panel and assign samples to best reference (Oryza/Ginkgo pilot).

Uses minimap2 aligned_fraction vs R0 (GenBank) + R1/R2 (medoid samples).
"""

import argparse
from pathlib import Path

import pandas as pd

from synteny_block_analysis import analyze_sample
from synteny_shared import (
    DEFAULT_CONFIG,
    DATA_ROOT,
    PAN_ROOT,
    canonical_sample_id,
    load_config,
    ref_fasta_path,
    resolve_sample_fastas,
    write_fasta,
)
from window_rearrangement_test import extract_ref_sequence


def pick_medoids(species, sample_dir, skip_summary_path, n_medoid=3, min_mf=0.80):
    df = pd.read_csv(skip_summary_path, sep="\t")
    sub = df[df["species"] == species].copy()
    sub["mapped_fraction"] = sub["mapped_windows"] / sub["total_windows"]
    sub = sub[sub["mapped_fraction"] >= min_mf].sort_values("mapped_fraction", ascending=False)
    fastas = resolve_sample_fastas(species, sample_dir)
    medoids = []
    for _, row in sub.head(n_medoid * 3).iterrows():
        sid = canonical_sample_id(row["sample_id"], species)
        if sid in fastas and sid not in [m[0] for m in medoids]:
            medoids.append((sid, fastas[sid], row["mapped_fraction"]))
        if len(medoids) >= n_medoid:
            break
    return medoids


def build_ref_panel(species, ref_dir, medoids, cache_dir):
    panel = {"R0": ref_fasta_path(ref_dir, species, cache_dir)}
    for i, (sid, fa, _) in enumerate(medoids, start=1):
        seq = str(next(__import__("Bio.SeqIO", fromlist=["SeqIO"]).parse(fa, "fasta")).seq).upper()
        out = Path(cache_dir) / f"{species}_R{i}_{sid}.fa"
        write_fasta(out, f"{species}_R{i}", seq)
        panel[f"R{i}"] = out
    return panel


def assign_samples(species, sample_dir, ref_panel, cfg, work_dir):
    fastas = resolve_sample_fastas(species, sample_dir)
    rows = []
    for sid, fa in sorted(fastas.items()):
        best_ref, best_af = None, -1.0
        scores = {}
        for ref_id, ref_fa in ref_panel.items():
            try:
                summary, _, _ = analyze_sample(sid, fa, ref_fa, cfg, work_dir / species / ref_id)
                scores[ref_id] = summary["aligned_fraction"]
                if summary["aligned_fraction"] > best_af:
                    best_af = summary["aligned_fraction"]
                    best_ref = ref_id
            except Exception as ex:
                scores[ref_id] = None
                print(f"  WARN {sid} vs {ref_id}: {ex}")
        rows.append({
            "species": species,
            "sample_id": sid,
            "assigned_ref": best_ref,
            "aligned_fraction_best": best_af,
            **{f"af_{k}": v for k, v in scores.items()},
        })
    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser(description="Multi-reference panel builder (pilot)")
    parser.add_argument("--config", default=str(DEFAULT_CONFIG))
    parser.add_argument("--sample-dir", default=str(DATA_ROOT / "00-rawfa/single_IR"))
    parser.add_argument("--ref-dir", default=str(DATA_ROOT / "ref/single_IR"))
    parser.add_argument("--skip-summary",
                        default=str(PAN_ROOT / "results/_archive/window_rearrangement_results"
                                    / "window_rearrangement_results_allSample_win200"
                                    / "window_skip_summary_by_sample.tsv"))
    parser.add_argument("--output-dir", default=str(PAN_ROOT / "results/synteny_analysis/ref_panel"))
    parser.add_argument("--work-dir", default=str(PAN_ROOT / "results/synteny_analysis/_ref_panel_work"))
    parser.add_argument("--species", nargs="*", default=None)
    args = parser.parse_args()

    cfg = load_config(args.config)
    pilot = cfg.get("multi_reference_pilot", {})
    species_list = args.species or pilot.get("species", ["Oryza_sativa", "Ginkgo_biloba"])
    n_medoid = int(pilot.get("medoid_candidates", 3))
    min_mf = float(pilot.get("min_mapped_fraction_for_medoid", 0.80))

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    work_dir = Path(args.work_dir)
    cache_dir = work_dir / "ref_cache"

    all_assign = []
    medoid_rows = []

    for species in species_list:
        print(f"{species}: picking medoids...", flush=True)
        medoids = pick_medoids(species, args.sample_dir, args.skip_summary, n_medoid, min_mf)
        if not medoids:
            print(f"  SKIP: no medoid candidates for {species}")
            continue
        for i, (sid, fa, mf) in enumerate(medoids, start=1):
            medoid_rows.append({
                "species": species,
                "ref_label": f"R{i}",
                "sample_id": sid,
                "mapped_fraction_win200": round(mf, 4),
            })
        ref_panel = build_ref_panel(species, args.ref_dir, medoids, cache_dir)
        panel_manifest = pd.DataFrame([
            {"species": species, "ref_id": k, "fasta_path": str(v)} for k, v in ref_panel.items()
        ])
        panel_manifest.to_csv(out_dir / f"{species}_ref_panel_manifest.tsv", sep="\t", index=False)

        print(f"  Assigning {len(resolve_sample_fastas(species, args.sample_dir))} samples...", flush=True)
        assign_df = assign_samples(species, args.sample_dir, ref_panel, cfg, work_dir)
        assign_df.to_csv(out_dir / f"{species}_sample_ref_assignment.tsv", sep="\t", index=False)
        all_assign.append(assign_df)

        # compare mean aligned_fraction R0 vs best
        for ref_col in [c for c in assign_df.columns if c.startswith("af_R")]:
            pass
        r0_mean = assign_df["af_R0"].mean() if "af_R0" in assign_df.columns else None
        best_mean = assign_df["aligned_fraction_best"].mean()
        print(f"  {species}: mean af_R0={r0_mean:.4f}, mean best={best_mean:.4f}", flush=True)

    if medoid_rows:
        pd.DataFrame(medoid_rows).to_csv(out_dir / "medoid_candidates.tsv", sep="\t", index=False)
    if all_assign:
        pd.concat(all_assign, ignore_index=True).to_csv(
            out_dir / "all_sample_ref_assignment.tsv", sep="\t", index=False
        )
    print("Done.", flush=True)


if __name__ == "__main__":
    main()
