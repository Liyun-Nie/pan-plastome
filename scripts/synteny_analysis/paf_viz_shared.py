"""Shared helpers for Phase 1b PAF visualization."""

from pathlib import Path

import numpy as np
import pandas as pd

from synteny_block_analysis import parse_paf_line
from synteny_shared import (
    DEFAULT_CONFIG,
    PAN_ROOT,
    cache_dir_from_cfg,
    load_config,
)
from window_rearrangement_test import SPECIES_LIST

RASTER_MAX_SAMPLES = 500


def manifest_relpath(path, root=PAN_ROOT):
    """Store figure paths relative to project root (cross-platform manifest)."""
    p = Path(path)
    if not p.is_absolute():
        return p.as_posix()
    try:
        return p.resolve().relative_to(Path(root).resolve()).as_posix()
    except ValueError:
        return p.as_posix()


def resolve_manifest_path(path, root=PAN_ROOT):
    """Resolve manifest path (relative to PAN_ROOT or absolute legacy entry)."""
    p = Path(path)
    if p.is_absolute():
        return p
    return (Path(root) / p).resolve()


def figures_dir(cfg, species=None):
    base = PAN_ROOT / cfg.get("results_dirs", {}).get(
        "synteny_output", "results/synteny_analysis"
    ) / "figures"
    return base / species if species else base


def load_species_summary(species):
    path = PAN_ROOT / "results/synteny_analysis/blocks" / f"{species}_synteny_summary.tsv"
    if not path.exists():
        raise FileNotFoundError(path)
    return pd.read_csv(path, sep="\t")


def qc_pass_mask(df, cfg):
    qc = cfg.get("qc", {})
    min_af = float(qc.get("min_aligned_fraction", 0.95))
    ratio = qc.get("query_ref_len_ratio", [0.95, 1.05])
    lo, hi = float(ratio[0]), float(ratio[1])
    length_ratio = df["query_len"] / df["ref_len"].replace(0, np.nan)
    return (
        (df["aligned_fraction"] >= min_af)
        & (length_ratio >= lo)
        & (length_ratio <= hi)
    )


def load_paf_segments(paf_path, min_mapq=1, min_alen=300):
    segments = []
    ref_len = None
    with open(paf_path, encoding="utf-8") as fh:
        for line in fh:
            rec = parse_paf_line(line)
            if rec is None:
                continue
            if ref_len is None:
                ref_len = rec["tlen"]
            if rec["mapq"] < min_mapq or rec["alen"] < min_alen:
                continue
            segments.append(rec)
    return segments, ref_len or 0


def species_paf_dir(cfg):
    return cache_dir_from_cfg(cfg, "paf", "results/synteny_analysis/paf")


def choose_viz_mode(n_samples, max_raster=RASTER_MAX_SAMPLES):
    return "raster" if n_samples <= max_raster else "envelope"


def load_representative_list(path=None):
    path = Path(path or PAN_ROOT / "config/representative_samples.tsv")
    return pd.read_csv(path, sep="\t")
