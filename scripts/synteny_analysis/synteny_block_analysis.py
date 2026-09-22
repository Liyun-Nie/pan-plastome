#!/usr/bin/env python3
"""
Level-1 synteny block analysis via minimap2 or blastn.

Per sample vs species reference:
  - aligner: minimap2 (default) or blastn (cross-validation subset)
  - dual-axis collinear block merge (collinear_block_merge.py)
  - detect inter-block inversion / translocation events
"""

import argparse
import subprocess
import sys
from pathlib import Path

import pandas as pd

from collinear_block_merge import merge_alignment_segments, normalize_blast_hsp
from synteny_shared import (
    DEFAULT_CONFIG,
    DATA_ROOT,
    PAN_ROOT,
    cache_dir_from_cfg,
    load_config,
    ref_fasta_path,
    resolve_sample_fastas,
)
from window_rearrangement_test import SPECIES_LIST

try:
    from Bio import SeqIO
except ImportError:
    sys.exit("ERROR: biopython required")


def parse_paf_line(line):
    if not line.strip() or line.startswith("#"):
        return None
    p = line.rstrip("\n").split("\t")
    if len(p) < 12:
        return None
    return {
        "qname": p[0],
        "qlen": int(p[1]),
        "qstart": int(p[2]),
        "qend": int(p[3]),
        "strand": p[4],
        "tname": p[5],
        "tlen": int(p[6]),
        "tstart": int(p[7]),
        "tend": int(p[8]),
        "nmatch": int(p[9]),
        "alen": int(p[10]),
        "mapq": int(p[11]),
    }


def run_minimap2(query_fa, ref_fa, threads=4, preset="asm5"):
    cmd = [
        "minimap2", f"-x{preset}", "-c", "-t", str(threads),
        str(ref_fa), str(query_fa),
    ]
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if proc.returncode != 0:
        raise RuntimeError(proc.stderr[:500])
    return proc.stdout


def paf_to_primary_records(paf_text, min_mapq=1, min_alen=500):
    rows = []
    for line in paf_text.splitlines():
        rec = parse_paf_line(line)
        if rec is None:
            continue
        if rec["mapq"] < min_mapq or rec["alen"] < min_alen:
            continue
        rows.append(rec)
    if not rows:
        return []
    rows.sort(key=lambda r: (r["qstart"], r["tstart"]))
    return rows


BLAST_OUTFMT = (
    "6 qseqid sseqid pident length mismatch gapopen "
    "qstart qend sstart send evalue bitscore sstrand"
)


def _blast_str_opt(value, default="no"):
    """Coerce YAML booleans to BLAST yes/no strings."""
    if value is None:
        return default
    if value is True:
        return "yes"
    if value is False:
        return "no"
    return str(value)


def run_blastn(query_fa, ref_fa, work_dir, threads=4, blast_cfg=None):
    """Run blastn (dc-megablast) query=sample vs subject=reference."""
    blast_cfg = blast_cfg or {}
    work_dir = Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)
    db_path = work_dir / "ref_db"
    out_path = work_dir / "blast.out"

    for cmd in (
        ["makeblastdb", "-in", str(ref_fa), "-dbtype", "nucl", "-out", str(db_path)],
        [
            "blastn",
            "-task", str(blast_cfg.get("task", "dc-megablast")),
            "-query", str(query_fa),
            "-db", str(db_path),
            "-out", str(out_path),
            "-outfmt", BLAST_OUTFMT,
            "-num_threads", str(threads),
            "-evalue", str(blast_cfg.get("evalue", "1e-20")),
            "-dust", _blast_str_opt(blast_cfg.get("dust"), "no"),
            "-word_size", str(blast_cfg.get("word_size", 11)),
        ],
    ):
        proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if proc.returncode != 0:
            raise RuntimeError(proc.stderr[:500])
    return out_path.read_text(encoding="utf-8")


def parse_blast_records(blast_text, min_identity=90, min_alen=500):
    """Parse BLAST outfmt6 into segment dicts with normalized coordinates (C2)."""
    rows = []
    for line in blast_text.splitlines():
        if not line.strip():
            continue
        p = line.split("\t")
        if len(p) < 13:
            continue
        pident = float(p[2])
        alen = int(p[3])
        if pident < min_identity or alen < min_alen:
            continue
        norm = normalize_blast_hsp(int(p[6]), int(p[7]), int(p[8]), int(p[9]), p[12])
        norm["nmatch"] = int(round(alen * pident / 100.0))
        norm["alen"] = alen
        rows.append(norm)
    if not rows:
        return []
    rows.sort(key=lambda r: (r["qstart"], r["tstart"]))
    return rows


def detect_inter_block_events(blocks, translocation_jump=10000):
    events = []
    if len(blocks) < 2:
        return events

    ordered = sorted(blocks, key=lambda b: b["qstart"])
    for i in range(1, len(ordered)):
        prev = ordered[i - 1]
        curr = ordered[i]
        event = None

        if prev["strand"] != curr["strand"]:
            event = "inversion"
        elif prev["strand"] == "+" and curr["tstart"] < prev["tend"]:
            event = "inversion"
        elif prev["strand"] == "-" and curr["tend"] > prev["tstart"]:
            event = "inversion"
        else:
            if prev["strand"] == "+":
                jump = curr["tstart"] - prev["tend"]
            else:
                jump = prev["tstart"] - curr["tend"]
            if abs(jump) > translocation_jump:
                event = "translocation"

        if event:
            events.append({
                "block_i": prev["block_id"],
                "block_j": curr["block_id"],
                "event_type": event,
                "prev_qend": prev["qend"],
                "curr_qstart": curr["qstart"],
                "prev_tend": prev["tend"],
                "curr_tstart": curr["tstart"],
                "prev_strand": prev["strand"],
                "curr_strand": curr["strand"],
            })
    return events


def analyze_sample(
    sample_id,
    sample_fa,
    ref_fa,
    cfg,
    work_dir,
    aligner="minimap2",
    save_paf=False,
    save_blast=False,
    paf_save_dir=None,
    blast_save_dir=None,
):
    mm = cfg.get("minimap2", {})
    blast_cfg = cfg.get("blastn", {})
    blk = cfg.get("block", {})
    threads = int(mm.get("threads", blast_cfg.get("threads", 4)))
    preset = mm.get("preset", "asm5")
    min_mapq = int(mm.get("min_mapq", 1))
    min_alen = int(blk.get("min_len", 500))
    max_gap = int(blk.get("max_gap", 200))
    overlap_tol = int(blk.get("overlap_tol", 50))
    trans_jump = int(blk.get("translocation_jump", 10000))
    min_identity = float(blast_cfg.get("min_identity", 90))

    q_len = len(str(next(SeqIO.parse(sample_fa, "fasta")).seq))
    r_len = len(str(next(SeqIO.parse(ref_fa, "fasta")).seq))

    work_dir = Path(work_dir)
    if aligner == "blastn":
        blast_text = run_blastn(
            sample_fa, ref_fa, work_dir / "_blast_db", threads=threads, blast_cfg=blast_cfg
        )
        if save_blast:
            out_dir = Path(blast_save_dir or work_dir)
            out_dir.mkdir(parents=True, exist_ok=True)
            (out_dir / f"{sample_id}.blast6").write_text(blast_text, encoding="utf-8")
        records = parse_blast_records(
            blast_text, min_identity=min_identity, min_alen=min_alen
        )
        n_segments_key = "n_blast_hsps"
    else:
        cached_paf = None
        if paf_save_dir is not None:
            cached_paf = Path(paf_save_dir) / f"{sample_id}.paf"
        if cached_paf is not None and cached_paf.exists() and cached_paf.stat().st_size > 0:
            paf_text = cached_paf.read_text(encoding="utf-8")
        else:
            paf_text = run_minimap2(sample_fa, ref_fa, threads=threads, preset=preset)
            if save_paf:
                out_dir = Path(paf_save_dir or work_dir)
                out_dir.mkdir(parents=True, exist_ok=True)
                (out_dir / f"{sample_id}.paf").write_text(paf_text, encoding="utf-8")
        records = paf_to_primary_records(paf_text, min_mapq=min_mapq, min_alen=min_alen)
        n_segments_key = "n_paf_segments"

    blocks = merge_alignment_segments(
        records, max_gap=max_gap, overlap_tol=overlap_tol
    )
    events = detect_inter_block_events(blocks, translocation_jump=trans_jump)

    aligned_bp = sum(b["block_len"] for b in blocks)
    aligned_fraction = round(aligned_bp / q_len, 6) if q_len else 0.0

    summary = {
        "sample_id": sample_id,
        "aligner": aligner,
        "query_len": q_len,
        "ref_len": r_len,
        "n_blocks": len(blocks),
        n_segments_key: len(records),
        "total_aligned_bp": aligned_bp,
        "aligned_fraction": aligned_fraction,
        "n_inversions": sum(1 for e in events if e["event_type"] == "inversion"),
        "n_translocations": sum(1 for e in events if e["event_type"] == "translocation"),
        "has_synteny_event": len(events) > 0,
    }

    block_rows = []
    for b in blocks:
        block_rows.append({
            "sample_id": sample_id,
            **b,
        })

    event_rows = [{"sample_id": sample_id, **e} for e in events]
    return summary, block_rows, event_rows


def process_species(
    species,
    sample_dir,
    ref_dir,
    out_dir,
    cfg,
    work_dir,
    aligner="minimap2",
    sample_filter=None,
    save_paf=False,
    save_blast=False,
    paf_cache_dir=None,
    blast_cache_dir=None,
    max_samples=None,
):
    sample_fastas = resolve_sample_fastas(species, sample_dir)
    if not sample_fastas:
        print(f"  SKIP {species}: no FASTAs")
        return

    ref_fa = ref_fasta_path(ref_dir, species, Path(work_dir) / "ref_cache")
    species_work = Path(work_dir) / species
    species_paf_dir = Path(paf_cache_dir) / species if paf_cache_dir else None
    species_blast_dir = Path(blast_cache_dir) / species if blast_cache_dir else None

    items = sorted(sample_fastas.items())
    if sample_filter is not None:
        items = [(sid, fa) for sid, fa in items if sid in sample_filter]
    if max_samples:
        items = items[:max_samples]

    print(f"  {species}: {len(items)} samples", flush=True)

    summaries, all_blocks, all_events = [], [], []
    for i, (sid, fa) in enumerate(items):
        if (i + 1) % 50 == 0 or i == 0:
            print(f"    [{i+1}/{len(items)}] {sid}", flush=True)
        try:
            summary, blocks, events = analyze_sample(
                sid,
                fa,
                ref_fa,
                cfg,
                species_work,
                aligner=aligner,
                save_paf=save_paf,
                save_blast=save_blast,
                paf_save_dir=species_paf_dir,
                blast_save_dir=species_blast_dir,
            )
            summary["species"] = species
            summaries.append(summary)
            for b in blocks:
                b["species"] = species
            all_blocks.extend(blocks)
            for e in events:
                e["species"] = species
            all_events.extend(events)
        except Exception as ex:
            print(f"    WARN {sid}: {ex}", flush=True)

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    if species_paf_dir:
        species_paf_dir.mkdir(parents=True, exist_ok=True)
    if species_blast_dir:
        species_blast_dir.mkdir(parents=True, exist_ok=True)

    if summaries:
        pd.DataFrame(summaries).to_csv(out_dir / f"{species}_synteny_summary.tsv", sep="\t", index=False)
    if all_blocks:
        pd.DataFrame(all_blocks).to_csv(out_dir / f"{species}_synteny_blocks.tsv", sep="\t", index=False)
    if all_events:
        pd.DataFrame(all_events).to_csv(out_dir / f"{species}_inter_block_events.tsv", sep="\t", index=False)
    else:
        pd.DataFrame(columns=[
            "species", "sample_id", "block_i", "block_j", "event_type",
        ]).to_csv(out_dir / f"{species}_inter_block_events.tsv", sep="\t", index=False)

    print(f"  Saved {species}: {len(summaries)} summaries, {len(all_events)} events", flush=True)


def load_sample_list(path):
    """Load sample filter TSV with columns species, sample_id."""
    df = pd.read_csv(path, sep="\t")
    by_species = {}
    for _, row in df.iterrows():
        by_species.setdefault(str(row["species"]), set()).add(str(row["sample_id"]))
    return by_species


def main():
    parser = argparse.ArgumentParser(description="Synteny block analysis (minimap2 or blastn)")
    parser.add_argument("--config", default=str(DEFAULT_CONFIG))
    parser.add_argument("--aligner", choices=["minimap2", "blastn"], default=None,
                        help="Alignment engine (default from config or minimap2)")
    parser.add_argument("--sample-dir", default=str(DATA_ROOT / "00-rawfa/single_IR"))
    parser.add_argument("--ref-dir", default=str(DATA_ROOT / "ref/single_IR"))
    parser.add_argument("--output-dir", default=str(PAN_ROOT / "results/synteny_analysis/blocks"))
    parser.add_argument("--work-dir", default=str(PAN_ROOT / "results/synteny_analysis/_align_work"))
    parser.add_argument("--species", nargs="*", default=None)
    parser.add_argument("--sample-list", default=None,
                        help="TSV (species, sample_id) to restrict samples")
    parser.add_argument("--max-samples", type=int, default=None)
    parser.add_argument("--save-paf-all", action="store_true", help="Save PAF for every sample")
    parser.add_argument("--save-blast", action="store_true", help="Save BLAST outfmt6 for every sample")
    args = parser.parse_args()

    cfg = load_config(args.config)
    aligner = args.aligner or cfg.get("aligner", "minimap2")
    if "block" not in cfg:
        cfg["block"] = {
            "min_len": 500,
            "max_gap": 200,
            "overlap_tol": 50,
            "translocation_jump": 10000,
        }
    if "minimap2" not in cfg:
        cfg["minimap2"] = {"preset": "asm5", "threads": 4, "min_mapq": 1}
    if "blastn" not in cfg:
        cfg["blastn"] = {
            "task": "dc-megablast",
            "threads": 4,
            "evalue": "1e-20",
            "dust": "no",
            "word_size": 11,
            "min_identity": 90,
        }

    species_list = args.species if args.species else SPECIES_LIST
    sample_by_species = load_sample_list(args.sample_list) if args.sample_list else {}

    paf_cache_dir = None
    blast_cache_dir = None
    if args.save_paf_all and aligner == "minimap2":
        paf_cache_dir = cache_dir_from_cfg(cfg, "paf", "results/synteny_analysis/paf")
        print(f"PAF cache: {paf_cache_dir}", flush=True)
    if args.save_blast and aligner == "blastn":
        blast_cache_dir = cache_dir_from_cfg(cfg, "blast_cache", "results/synteny_analysis/blast")
        print(f"BLAST cache: {blast_cache_dir}", flush=True)

    print(f"Aligner: {aligner}", flush=True)
    print(f"Species: {len(species_list)}", flush=True)
    for species in species_list:
        print(f"{species}...", flush=True)
        sample_filter = sample_by_species.get(species)
        if args.sample_list and sample_filter is None:
            continue
        process_species(
            species,
            args.sample_dir,
            args.ref_dir,
            args.output_dir,
            cfg,
            Path(args.work_dir),
            aligner=aligner,
            sample_filter=sample_filter,
            save_paf=args.save_paf_all and aligner == "minimap2",
            save_blast=args.save_blast and aligner == "blastn",
            paf_cache_dir=paf_cache_dir,
            blast_cache_dir=blast_cache_dir,
            max_samples=args.max_samples,
        )
    print("Done.", flush=True)


if __name__ == "__main__":
    main()
