#!/usr/bin/env python3
"""Build STRUCTURAL_SCREENING_SUMMARY.tsv and Structural_Screening_Results.md."""

import argparse
from datetime import date
from pathlib import Path

import pandas as pd

from paf_viz_shared import figures_dir, load_representative_list, qc_pass_mask
from synteny_shared import DEFAULT_CONFIG, PAN_ROOT, load_config
from window_rearrangement_test import SPECIES_LIST


def load_optional(path):
    path = Path(path)
    if not path.exists() or path.stat().st_size <= 1:
        return pd.DataFrame()
    try:
        return pd.read_csv(path, sep="\t")
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def species_figure_path(cfg, species):
    manifest = figures_dir(cfg) / "phase1b_figure_manifest.tsv"
    if not manifest.exists():
        return ""
    df = pd.read_csv(manifest, sep="\t")
    sub = df[(df["species"] == species) & (df["figure_type"].isin(["raster", "envelope"]))]
    if sub.empty:
        return ""
    return str(sub.iloc[0]["path"])


def build_summary(cfg):
    sensitivity = load_optional(PAN_ROOT / "results/synteny_analysis/sensitivity/block_param_sensitivity.tsv")
    sens_max_events = int(sensitivity["n_inter_block_total"].max()) if len(sensitivity) else 0
    agreement = load_optional(PAN_ROOT / "results/synteny_analysis/validation/method_agreement.tsv")
    n_concordant = int((agreement["agreement"] == "concordant").sum()) if len(agreement) else 0
    n_agreement = len(agreement)
    anomalies = load_optional(PAN_ROOT / "results/synteny_analysis/local_anomalies/all_local_anomalies.tsv")

    rows = []
    for species in SPECIES_LIST:
        sp_sum = pd.read_csv(
            PAN_ROOT / "results/synteny_analysis/blocks" / f"{species}_synteny_summary.tsv",
            sep="\t",
        )
        qc = sp_sum[qc_pass_mask(sp_sum, cfg)]
        sp_anom = anomalies[anomalies["species"] == species] if len(anomalies) else pd.DataFrame()
        n_anom_samples = sp_anom["sample_id"].nunique() if len(sp_anom) else 0
        n_anom_events = len(sp_anom)
        n_inter_primary = int(sp_sum["has_synteny_event"].sum())

        rows.append({
            "species": species,
            "n_samples": len(sp_sum),
            "n_qc_pass": len(qc),
            "mean_aligned_fraction": round(sp_sum["aligned_fraction"].mean(), 6),
            "median_aligned_fraction": round(sp_sum["aligned_fraction"].median(), 6),
            "n_inter_block_primary": n_inter_primary,
            "n_inter_block_sensitivity_max": sens_max_events,
            "n_local_anomaly_samples": n_anom_samples,
            "n_local_anomaly_events": n_anom_events,
            "primary_figure": species_figure_path(cfg, species),
        })

    summary = pd.DataFrame(rows)
    summary["method_agreement_concordant"] = f"{n_concordant}/{n_agreement}" if n_agreement else "NA"
    summary["pggb_note"] = summary["species"].map(
        lambda s: "Glycine macro graph: no SV (external)" if s.startswith("Glycine") else ""
    )
    return summary, agreement, anomalies, sensitivity


def df_to_md_table(df):
    if df.empty:
        return "_None_"
    cols = list(df.columns)
    lines = ["| " + " | ".join(cols) + " |", "| " + " | ".join("---" for _ in cols) + " |"]
    for _, row in df.iterrows():
        lines.append("| " + " | ".join(str(row[c]) for c in cols) + " |")
    return "\n".join(lines)


def write_markdown(summary, agreement, anomalies, sensitivity, syri, out_path):
    total_samples = int(summary["n_samples"].sum())
    total_inter = int(summary["n_inter_block_primary"].sum())
    n_anom = len(anomalies)
    discordant = agreement[agreement["agreement"] == "discordant"] if len(agreement) else pd.DataFrame()
    syri_counts = syri["verdict"].value_counts().to_dict() if len(syri) and "verdict" in syri.columns else {}

    lines = [
        "# Structural Screening Results (Alignment-first v4.1)",
        "",
        f"**Generated**: {date.today().isoformat()}",
        f"**Panel**: {total_samples} single-IR plastome assemblies, 17 species",
        "",
        "## Executive summary",
        "",
        f"- **Primary block screen** (minimap2, dual-axis merge, 500/200/10 kb): "
        f"**{total_inter}** inter-block inversion/translocation events across the panel.",
        f"- **Parameter sensitivity** (min_alen=300): max **{int(sensitivity['n_inter_block_total'].max()) if len(sensitivity) else 0}** "
        f"inter-block events (Hemerocallis_citrina only; see sensitivity table).",
        f"- **PAF local anomalies** (within-chain): **{n_anom}** events in "
        f"**{anomalies['sample_id'].nunique() if len(anomalies) else 0}** samples.",
        f"- **Nucmer/SyRI validation** ({len(syri) if len(syri) else 0} representatives; SyRI unavailable → nucmer): "
        f"{syri_counts}.",
        f"- **BLAST cross-check** (17 representatives): "
        f"**{(agreement['agreement'] == 'concordant').sum() if len(agreement) else 'NA'}/{len(agreement) if len(agreement) else 'NA'}** concordant.",
        "",
        "## Scope statement (for co-authors)",
        "",
        "Nucleotide landscape in the main text = Snippy SNV/InDel (max INDEL 57 bp). "
        "This package = supplementary **macro-synteny screening**; does not claim zero structural variation.",
        "",
        "## Species summary",
        "",
        df_to_md_table(summary),
        "",
        "## Discordant aligner cross-check (SyRI/nucmer priority)",
        "",
    ]
    if len(discordant):
        lines.append(df_to_md_table(discordant[["species", "sample_id", "n_blocks_minimap2", "n_blocks_blastn"]]))
    else:
        lines.append("_None_")

    if len(syri):
        sub = syri[syri["verdict"] != "no_macro_sv"][
            ["species", "sample_id", "verdict", "detail", "selection_reason"]
        ]
        lines.extend([
            "",
            "## Nucmer non-collinear verdicts (subset)",
            "",
            df_to_md_table(sub) if len(sub) else "_All representatives: no_macro_sv (contiguous collinear nucmer alignment)_",
        ])

    lines.extend([
        "",
        "## Key paths",
        "",
        "- `results/synteny_analysis/blocks/` — post-fix block summaries",
        "- `results/synteny_analysis/sensitivity/block_param_sensitivity.tsv`",
        "- `results/synteny_analysis/validation/method_agreement.tsv`",
        "- `results/synteny_analysis/local_anomalies/`",
        "- `results/synteny_analysis/figures/phase1b_figure_manifest.tsv`",
        "- `results/synteny_analysis/syri/syri_representative_verdict.tsv`",
        "",
    ])
    Path(out_path).write_text("\n".join(lines), encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description="Structural screening report orchestrator")
    parser.add_argument("--config", default=str(DEFAULT_CONFIG))
    parser.add_argument("--output-dir", default=str(PAN_ROOT / "results/synteny_analysis"))
    args = parser.parse_args()

    cfg = load_config(args.config)
    out_dir = Path(args.output_dir)
    summary, agreement, anomalies, sensitivity = build_summary(cfg)
    syri = load_optional(out_dir / "syri/syri_representative_verdict.tsv")

    tsv_path = out_dir / "STRUCTURAL_SCREENING_SUMMARY.tsv"
    summary.to_csv(tsv_path, sep="\t", index=False)
    md_path = out_dir / "Structural_Screening_Results.md"
    write_markdown(summary, agreement, anomalies, sensitivity, syri, md_path)
    print(f"Wrote {tsv_path}")
    print(f"Wrote {md_path}")


if __name__ == "__main__":
    main()
