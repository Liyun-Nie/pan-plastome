#!/usr/bin/env python3
"""
Representative assembly validation via SyRI (if available) or nucmer fallback.

Phase 4: ~20-30 priority samples for macro-SV validation attempt.
"""

import argparse
import re
import subprocess
import sys
from pathlib import Path

import pandas as pd

from synteny_shared import (
    DEFAULT_CONFIG,
    DATA_ROOT,
    PAN_ROOT,
    load_config,
    ref_fasta_path,
    resolve_sample_fastas,
)

try:
    from Bio import SeqIO
except ImportError:
    SeqIO = None


def build_priority_panel(cfg):
    """Assemble ~25 samples: anomaly > discordant > sensitivity flags > representatives."""
    rows = []
    seen = set()

    def add(species, sample_id, priority, reason):
        key = (species, str(sample_id))
        if key in seen:
            return
        seen.add(key)
        rows.append({
            "species": species,
            "sample_id": str(sample_id),
            "priority": priority,
            "selection_reason": reason,
        })

    anom_path = PAN_ROOT / "results/synteny_analysis/local_anomalies/all_local_anomalies.tsv"
    if anom_path.exists() and anom_path.stat().st_size > 1:
        try:
            anom = pd.read_csv(anom_path, sep="\t")
        except pd.errors.EmptyDataError:
            anom = pd.DataFrame()
        if len(anom) and "sample_id" in anom.columns:
            for sid in anom["sample_id"].astype(str).unique():
                sp = anom.loc[anom["sample_id"].astype(str) == sid, "species"].iloc[0]
                add(sp, sid, 1, "local_anomaly")

    agree_path = PAN_ROOT / "results/synteny_analysis/validation/method_agreement.tsv"
    if agree_path.exists():
        agree = pd.read_csv(agree_path, sep="\t")
        disc = agree[agree["agreement"] == "discordant"]
        for _, r in disc.iterrows():
            add(r["species"], r["sample_id"], 2, "blast_minimap2_discordant")

    for sid in ("hh311", "hh361", "hh391"):
        add("Hemerocallis_citrina", sid, 3, "sensitivity_min_alen300_inversion")

    rep_path = PAN_ROOT / "config/representative_samples.tsv"
    if rep_path.exists():
        reps = pd.read_csv(rep_path, sep="\t")
        for _, r in reps.iterrows():
            add(r["species"], r["sample_id"], 4, "species_representative")

    for item in cfg.get("priority_dotplot_panel", []):
        add(item["species"], item["sample_id"], 5, f"priority_panel_{item.get('role', '')}")

    df = pd.DataFrame(rows).sort_values(["priority", "species", "sample_id"])
    return df.head(30)


def run_nucmer(ref_fa, query_fa, work_dir, prefix="out"):
    work_dir = Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)
    delta = work_dir / f"{prefix}.delta"
    for cmd in (
        ["nucmer", "--maxmatch", "-l", "40", "-c", "100", "-p", str(work_dir / prefix),
         str(ref_fa), str(query_fa)],
    ):
        proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if proc.returncode != 0:
            raise RuntimeError(proc.stderr[:800])
    return delta


def parse_show_coords(delta_path):
    proc = subprocess.run(
        ["show-coords", "-r", "-l", "-d", str(delta_path)],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
    )
    if proc.returncode != 0:
        raise RuntimeError(proc.stderr[:500])
    clusters = []
    for line in proc.stdout.splitlines():
        line = line.strip()
        if not line or line.startswith("=") or line.startswith("[") or "NUCMmer" in line:
            continue
        if line.endswith(".fa") or line.endswith(".fasta") or "\t" in line and "|" not in line:
            continue
        parts = [p.strip() for p in line.split("|")]
        if len(parts) < 3:
            continue
        try:
            ref_nums = parts[0].split()
            qry_nums = parts[1].split()
            len_nums = parts[2].split()
            if len(ref_nums) < 2 or len(qry_nums) < 2:
                continue
            ref_start, ref_end = int(ref_nums[0]), int(ref_nums[1])
            qry_start, qry_end = int(qry_nums[0]), int(qry_nums[1])
            alen = int(len_nums[0]) if len_nums else abs(qry_end - qry_start)
            direction = parts[4].strip() if len(parts) > 4 else "F"
            clusters.append({
                "ref_start": ref_start,
                "ref_end": ref_end,
                "qry_start": qry_start,
                "qry_end": qry_end,
                "align_len": alen,
                "direction": direction,
            })
        except (ValueError, IndexError):
            continue
    return clusters


def classify_nucmer_clusters(clusters, ref_len, qry_len):
    if not clusters:
        return "inconclusive", "no nucmer alignments", 0, 0

    def is_reverse(c):
        d = c["direction"].upper()
        return d.startswith(("R", "-", "N")) or "REVERSE" in d

    n_rev = sum(1 for c in clusters if is_reverse(c))
    n_fwd = len(clusters) - n_rev
    q_cov = sum(c.get("align_len", abs(c["qry_end"] - c["qry_start"])) for c in clusters)
    q_frac = q_cov / max(1, qry_len)

    if n_rev >= 1 and n_fwd >= 1:
        return "inversion_candidate", f"{n_rev} reverse + {n_fwd} forward blocks", n_rev, n_fwd
    if len(clusters) >= 3 and q_frac > 0.9:
        return "rearrangement_candidate", f"{len(clusters)} alignment clusters", n_rev, n_fwd
    if len(clusters) <= 2 and q_frac >= 0.95 and n_rev == 0:
        return "no_macro_sv", "contiguous collinear alignment", n_rev, n_fwd
    if q_frac >= 0.90 and n_rev == 0:
        return "no_macro_sv", f"high coverage ({q_frac:.3f}), no reverse blocks", n_rev, n_fwd
    return "inconclusive", f"q_cov={q_frac:.3f}, clusters={len(clusters)}", n_rev, n_fwd


def run_syri_if_available(ref_fa, query_fa, work_dir):
    proc = subprocess.run(["which", "syri"], stdout=subprocess.PIPE, text=True)
    if proc.returncode != 0:
        return None
    work_dir = Path(work_dir)
    out_prefix = work_dir / "syri_out"
    cmd = ["syri", "-c", str(query_fa), "-r", str(ref_fa), "--prefix", str(out_prefix)]
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if proc.returncode != 0:
        return {"tool": "syri", "status": "failed", "detail": proc.stderr[:300]}
    return {"tool": "syri", "status": "completed", "detail": str(out_prefix)}


def screen_sample(species, sample_id, cfg, sample_dir, ref_dir, work_root):
    fastas = resolve_sample_fastas(species, sample_dir)
    query_fa = fastas.get(str(sample_id))
    if not query_fa:
        return {"verdict": "inconclusive", "detail": "sample fasta not found", "tool": "none"}

    ref_cache = work_root / "ref_cache"
    ref_fa = ref_fasta_path(ref_dir, species, ref_cache)
    if SeqIO is None:
        raise RuntimeError("biopython required")
    q_len = len(str(next(SeqIO.parse(query_fa, "fasta")).seq))
    r_len = len(str(next(SeqIO.parse(ref_fa, "fasta")).seq))

    sample_work = work_root / species / str(sample_id)
    tool = "nucmer"
    detail = ""
    try:
        syri_res = run_syri_if_available(ref_fa, query_fa, sample_work)
        if syri_res and syri_res.get("status") == "completed":
            tool = "syri"
            detail = syri_res["detail"]
            verdict = "syri_completed"
            return {
                "tool": tool,
                "verdict": verdict,
                "detail": detail,
                "n_reverse_blocks": "",
                "n_forward_blocks": "",
                "query_len": q_len,
                "ref_len": r_len,
            }

        delta = run_nucmer(ref_fa, query_fa, sample_work)
        clusters = parse_show_coords(delta)
        verdict, detail, n_rev, n_fwd = classify_nucmer_clusters(clusters, r_len, q_len)
        return {
            "tool": tool,
            "verdict": verdict,
            "detail": detail,
            "n_reverse_blocks": n_rev,
            "n_forward_blocks": n_fwd,
            "n_alignment_clusters": len(clusters),
            "query_len": q_len,
            "ref_len": r_len,
        }
    except Exception as ex:
        return {
            "tool": tool,
            "verdict": "inconclusive",
            "detail": str(ex)[:300],
            "n_reverse_blocks": "",
            "n_forward_blocks": "",
            "query_len": q_len,
            "ref_len": r_len,
        }


def main():
    parser = argparse.ArgumentParser(description="SyRI/nucmer representative validation")
    parser.add_argument("--config", default=str(DEFAULT_CONFIG))
    parser.add_argument("--sample-dir", default=str(DATA_ROOT / "00-rawfa/single_IR"))
    parser.add_argument("--ref-dir", default=str(DATA_ROOT / "ref/single_IR"))
    parser.add_argument("--panel-tsv", default=None, help="Override priority panel TSV")
    parser.add_argument("--output", default=str(PAN_ROOT / "results/synteny_analysis/syri/syri_representative_verdict.tsv"))
    parser.add_argument("--work-dir", default=str(PAN_ROOT / "results/synteny_analysis/syri/_work"))
    args = parser.parse_args()

    cfg = load_config(args.config)
    if args.panel_tsv:
        panel = pd.read_csv(args.panel_tsv, sep="\t")
    else:
        panel = build_priority_panel(cfg)
        panel_path = Path(args.output).parent / "syri_priority_panel.tsv"
        panel_path.parent.mkdir(parents=True, exist_ok=True)
        panel.to_csv(panel_path, sep="\t", index=False)
        print(f"Priority panel ({len(panel)}): {panel_path}", flush=True)

    results = []
    for i, row in panel.iterrows():
        sp, sid = row["species"], str(row["sample_id"])
        print(f"[{i+1}/{len(panel)}] {sp}/{sid}...", flush=True)
        res = screen_sample(sp, sid, cfg, args.sample_dir, args.ref_dir, Path(args.work_dir))
        results.append({**row.to_dict(), **res})

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(results).to_csv(out, sep="\t", index=False)
    print(f"Wrote {out} ({len(results)} samples)", flush=True)


if __name__ == "__main__":
    main()
