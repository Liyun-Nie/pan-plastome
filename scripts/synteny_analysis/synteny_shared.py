"""Shared helpers for synteny block analysis modules."""

from pathlib import Path

import yaml
from Bio import SeqIO

from audit_window_skips import canonical_sample_id, resolve_sample_fastas
from window_rearrangement_test import SPECIES_LIST, extract_ref_sequence

SCRIPT_DIR = Path(__file__).resolve().parent
SCRIPTS_ROOT = SCRIPT_DIR.parent
ANALYSIS_ROOT = SCRIPTS_ROOT.parent
PAN_ROOT = ANALYSIS_ROOT
DATA_ROOT = ANALYSIS_ROOT / "data"
DEFAULT_CONFIG = PAN_ROOT / "config" / "synteny_screening.yml"


def load_config(path=None):
    path = Path(path or DEFAULT_CONFIG)
    with open(path, encoding="utf-8") as f:
        return yaml.safe_load(f)


def cache_dir_from_cfg(cfg, section, default_rel, root=None):
    """Resolve cache_dir from config section relative to PAN_ROOT."""
    root = Path(root or PAN_ROOT)
    rel = cfg.get(section, {}).get("cache_dir", default_rel)
    path = Path(rel)
    return path if path.is_absolute() else root / path


def write_fasta(path, seq_id, sequence):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(f">{seq_id}\n")
        for i in range(0, len(sequence), 80):
            fh.write(sequence[i:i + 80] + "\n")


def ref_fasta_path(ref_dir, species, cache_dir):
    """Export species reference from GenBank to cached FASTA."""
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    out = cache_dir / f"{species}.ref.fa"
    if out.exists() and out.stat().st_size > 0:
        return out
    gb = Path(ref_dir) / f"{species}.gb"
    seq = extract_ref_sequence(gb)
    if not seq:
        raise ValueError(f"Cannot extract reference sequence from {gb}")
    write_fasta(out, species, seq.upper())
    return out
