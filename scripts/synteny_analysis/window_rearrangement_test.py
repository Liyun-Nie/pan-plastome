#!/usr/bin/env python3
"""
DEPRECATED (2026-05-27): Window-based structural rearrangement test.

Superseded by alignment-first screening: synteny_block_analysis.py (v4.1).
See docs/.archive/2026-06/WINDOW_PIPELINE_DEPRECATED.md. Do not use for new manuscript evidence.

Window-based Structural Rearrangement Test for Plastomes
========================================================
Splits sample and reference plastome sequences into non-overlapping windows,
aligns each sample window to the reference window database via BLAST, and
checks whether the mapping order is monotonically increasing (= no rearrangement).

Responds to Reviewer #1 Comment 2: limitations of short-read SV detection.

Dependencies: biopython, pandas, numpy
External tools: makeblastdb, blastn (NCBI BLAST+)

Usage:
    python window_rearrangement_test.py \
        --sample-dir <00-rawfa/single_IR/> \
        --ref-dir <ref/single_IR/> \
        --output-dir <window_rearrangement_results/> \
        --window-size 500
"""

import argparse
import os
import sys
import subprocess
import tempfile
import shutil
import warnings
from pathlib import Path
from collections import defaultdict

import pandas as pd
import numpy as np

try:
    from Bio import SeqIO
    from Bio.Seq import Seq
except ImportError:
    sys.exit("ERROR: BioPython is required. Install with: pip install biopython")

# ---------------------------------------------------------------------------
# Species list
# ---------------------------------------------------------------------------

SPECIES_LIST = [
    'Citrullus_lanatus', 'Fagopyrum_tataricum', 'Ginkgo_biloba',
    'Glycine_max', 'Glycine_soja',
    'Gossypium_barbadense', 'Gossypium_hirsutum',
    'Hemerocallis_citrina', 'Nelumbo_nucifera',
    'Oryza_rufipogon', 'Oryza_sativa',
    'Pisum_sativum',
    'Prunus_mume',
    'Solanum_brevicaule', 'Solanum_candolleanum',
    'Zea_mays', 'Ziziphus_jujuba',
]


# ---------------------------------------------------------------------------
# Utility functions
# ---------------------------------------------------------------------------

def split_into_windows(sequence, window_size=500):
    """Split a sequence into non-overlapping windows."""
    windows = []
    seq_len = len(sequence)
    for i in range(0, seq_len, window_size):
        win_seq = sequence[i:i + window_size]
        if len(win_seq) >= window_size * 0.5:
            windows.append((i, win_seq))
    return windows


def write_windows_fasta(windows, prefix, out_file):
    """Write windows to a FASTA file."""
    records = []
    for idx, (start, seq) in enumerate(windows):
        from Bio.SeqRecord import SeqRecord
        rec = SeqRecord(Seq(seq), id=f"{prefix}_win{idx}", description=f"start={start}")
        records.append(rec)
    SeqIO.write(records, out_file, 'fasta')
    return len(records)


def extract_ref_sequence(gb_file):
    """Extract nucleotide sequence from GenBank file."""
    for record in SeqIO.parse(gb_file, 'genbank'):
        return str(record.seq)
    return None


def run_blast(query_fasta, db_path, output_file, evalue=1e-10):
    """Run blastn and return output file path."""
    cmd = [
        'blastn',
        '-query', str(query_fasta),
        '-db', str(db_path),
        '-out', str(output_file),
        '-outfmt', '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore',
        '-evalue', str(evalue),
        '-max_target_seqs', '5',
        '-num_threads', '1',
    ]
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        return True
    except subprocess.CalledProcessError as e:
        print(f"  WARNING: BLAST failed: {e.stderr[:200]}")
        return False
    except FileNotFoundError:
        return False


def parse_blast_results(blast_output, min_identity=90, min_coverage=0.8, window_size=500):
    """Parse BLAST tabular output and extract best hit per query window."""
    if not os.path.exists(blast_output) or os.path.getsize(blast_output) == 0:
        return {}

    df = pd.read_csv(blast_output, sep='\t', header=None,
                     names=['qseqid', 'sseqid', 'pident', 'length', 'mismatch',
                            'gapopen', 'qstart', 'qend', 'sstart', 'send',
                            'evalue', 'bitscore'])

    df = df[df['pident'] >= min_identity]
    df = df[df['length'] >= window_size * min_coverage]

    best_hits = {}
    for qid, group in df.groupby('qseqid'):
        best = group.loc[group['bitscore'].idxmax()]
        query_win_num = int(qid.split('_win')[-1])
        subj_win_num = int(best['sseqid'].split('_win')[-1])
        best_hits[query_win_num] = {
            'ref_win': subj_win_num,
            'pident': best['pident'],
            'length': best['length'],
            'bitscore': best['bitscore'],
        }

    return best_hits


def check_monotonicity(mapping, n_windows, edge_tolerance=1):
    """
    Check if the window mapping is monotonically increasing.
    Allow edge_tolerance windows at start/end to deviate (circular genome origin).
    """
    if not mapping:
        return True, []

    sorted_queries = sorted(mapping.keys())
    inner_queries = [q for q in sorted_queries
                     if q >= edge_tolerance and q < n_windows - edge_tolerance]

    if len(inner_queries) < 2:
        return True, []

    violations = []
    for i in range(1, len(inner_queries)):
        q_prev = inner_queries[i - 1]
        q_curr = inner_queries[i]
        ref_prev = mapping[q_prev]['ref_win']
        ref_curr = mapping[q_curr]['ref_win']
        if ref_curr <= ref_prev:
            violations.append({
                'query_win': q_curr,
                'ref_win_expected_gt': ref_prev,
                'ref_win_actual': ref_curr,
            })

    return len(violations) == 0, violations


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

def process_species(species, sample_dir, ref_dir, out_dir, window_size,
                    min_identity, min_coverage, edge_tolerance, max_samples):
    """Process all samples for one species."""
    ref_gb = ref_dir / f"{species}.gb"
    if not ref_gb.exists():
        print(f"  WARNING: {ref_gb} not found, skipping {species}")
        return None

    ref_seq = extract_ref_sequence(ref_gb)
    if not ref_seq:
        print(f"  WARNING: Could not extract sequence from {ref_gb}")
        return None

    ref_windows = split_into_windows(ref_seq, window_size)
    n_ref_windows = len(ref_windows)

    sample_species_dir = sample_dir / species
    if species == 'Glycine_max' or species == 'Glycine_soja':
        pass

    if not sample_species_dir.is_dir():
        glycine_dir = sample_dir / 'Glycine'
        if glycine_dir.is_dir() and species in ('Glycine_max', 'Glycine_soja'):
            sample_txt = sample_dir / species / 'sample.txt'
            if sample_txt.exists():
                with open(sample_txt) as f:
                    valid_samples = {line.strip() for line in f if line.strip()}
            else:
                valid_samples = set()
            sample_species_dir = glycine_dir
        else:
            print(f"  WARNING: Sample directory {sample_species_dir} not found")
            return None
    else:
        valid_samples = None

    sample_fastas = sorted(sample_species_dir.glob('*.fasta')) + \
                    sorted(sample_species_dir.glob('*.fa'))

    if valid_samples is not None:
        sample_fastas = [f for f in sample_fastas
                         if f.stem.replace('.single_IR', '') in valid_samples
                         or f.stem in valid_samples]

    if max_samples and len(sample_fastas) > max_samples:
        rng = np.random.RandomState(42)
        indices = rng.choice(len(sample_fastas), max_samples, replace=False)
        sample_fastas = [sample_fastas[i] for i in sorted(indices)]

    if not sample_fastas:
        print(f"  WARNING: No sample FASTA files found for {species}")
        return None

    print(f"  Reference: {n_ref_windows} windows, Samples: {len(sample_fastas)}")

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)

        ref_fasta = tmpdir / 'ref_windows.fa'
        write_windows_fasta(ref_windows, 'ref', ref_fasta)

        db_path = tmpdir / 'ref_db'
        try:
            subprocess.run(
                ['makeblastdb', '-in', str(ref_fasta), '-dbtype', 'nucl',
                 '-out', str(db_path)],
                check=True, capture_output=True, text=True
            )
        except (subprocess.CalledProcessError, FileNotFoundError) as e:
            print(f"  WARNING: makeblastdb failed for {species}: {e}")
            return None

        all_mappings = []
        rearranged_samples = []

        for sample_fa in sample_fastas:
            sample_id = sample_fa.stem.replace('.single_IR', '')

            try:
                sample_seq = str(next(SeqIO.parse(sample_fa, 'fasta')).seq)
            except Exception:
                continue

            sample_windows = split_into_windows(sample_seq, window_size)
            n_sample_windows = len(sample_windows)

            query_fasta = tmpdir / 'query.fa'
            write_windows_fasta(sample_windows, sample_id, query_fasta)

            blast_out = tmpdir / 'blast_result.txt'
            success = run_blast(query_fasta, db_path, blast_out)
            if not success:
                continue

            mapping = parse_blast_results(
                blast_out, min_identity=min_identity,
                min_coverage=min_coverage, window_size=window_size
            )

            is_monotonic, violations = check_monotonicity(
                mapping, n_sample_windows, edge_tolerance=edge_tolerance
            )

            mapped_wins = len(mapping)
            for q_win, hit in mapping.items():
                all_mappings.append({
                    'sample_id': sample_id,
                    'query_win': q_win,
                    'ref_win': hit['ref_win'],
                    'pident': hit['pident'],
                    'bitscore': hit['bitscore'],
                })

            if not is_monotonic:
                rearranged_samples.append({
                    'sample_id': sample_id,
                    'n_violations': len(violations),
                    'violation_positions': ';'.join(
                        f"q{v['query_win']}->r{v['ref_win_actual']}"
                        for v in violations
                    ),
                    'mapped_windows': mapped_wins,
                    'total_windows': n_sample_windows,
                })

    mapping_df = pd.DataFrame(all_mappings)
    mapping_file = out_dir / f'{species}_window_mapping.tsv'
    mapping_df.to_csv(mapping_file, sep='\t', index=False)

    summary = {
        'species': species,
        'n_samples': len(sample_fastas),
        'n_ref_windows': n_ref_windows,
        'n_rearranged': len(rearranged_samples),
        'rearrangement_rate': len(rearranged_samples) / len(sample_fastas) if sample_fastas else 0,
    }

    if rearranged_samples:
        rearr_df = pd.DataFrame(rearranged_samples)
        rearr_file = out_dir / f'{species}_rearrangement_details.tsv'
        rearr_df.to_csv(rearr_file, sep='\t', index=False)

    return summary


def main():
    parser = argparse.ArgumentParser(
        description="Window-based structural rearrangement test for plastomes")
    parser.add_argument('--sample-dir', required=True,
                        help='Directory containing {Species}/ subdirectories with FASTA files')
    parser.add_argument('--ref-dir', required=True,
                        help='Directory containing {Species}.gb GenBank reference files')
    parser.add_argument('--output-dir', default='window_rearrangement_results',
                        help='Output directory')
    parser.add_argument('--window-size', type=int, default=500,
                        help='Window size in bp (default: 500)')
    parser.add_argument('--min-identity', type=float, default=90,
                        help='Minimum BLAST identity percent (default: 90)')
    parser.add_argument('--min-coverage', type=float, default=0.8,
                        help='Minimum alignment coverage fraction (default: 0.8)')
    parser.add_argument('--edge-tolerance', type=int, default=1,
                        help='Number of edge windows to tolerate (default: 1)')
    parser.add_argument('--max-samples', type=int, default=None,
                        help='Max samples per species (for testing; default: all)')
    parser.add_argument('--species', nargs='*', default=None,
                        help='Specific species to process (default: all)')
    args = parser.parse_args()

    warnings.warn(
        "window_rearrangement_test.py is DEPRECATED (v4.1). "
        "Use synteny_block_analysis.py instead. "
        "See docs/.archive/2026-06/WINDOW_PIPELINE_DEPRECATED.md.",
        DeprecationWarning,
        stacklevel=1,
    )

    sample_dir = Path(args.sample_dir)
    ref_dir = Path(args.ref_dir)
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    for tool in ['makeblastdb', 'blastn']:
        if shutil.which(tool) is None:
            sys.exit(f"ERROR: {tool} not found in PATH. Install NCBI BLAST+.")

    species_list = args.species if args.species else SPECIES_LIST

    print("=" * 60)
    print("Window-based Structural Rearrangement Test")
    print(f"Window size: {args.window_size} bp")
    print(f"Species: {len(species_list)}")
    print("=" * 60)

    all_summaries = []
    for species in species_list:
        print(f"\nProcessing {species}...")
        summary = process_species(
            species, sample_dir, ref_dir, out_dir,
            window_size=args.window_size,
            min_identity=args.min_identity,
            min_coverage=args.min_coverage,
            edge_tolerance=args.edge_tolerance,
            max_samples=args.max_samples,
        )
        if summary:
            all_summaries.append(summary)
            status = "REARRANGEMENT DETECTED" if summary['n_rearranged'] > 0 else "No rearrangement"
            print(f"  Result: {status} ({summary['n_rearranged']}/{summary['n_samples']})")

    if all_summaries:
        summary_df = pd.DataFrame(all_summaries)
        summary_file = out_dir / 'all_species_summary.tsv'
        summary_df.to_csv(summary_file, sep='\t', index=False)
        print(f"\n{'=' * 60}")
        print(f"Summary saved to {summary_file}")
        print(f"Total species processed: {len(all_summaries)}")
        total_rearr = sum(s['n_rearranged'] for s in all_summaries)
        total_samples = sum(s['n_samples'] for s in all_summaries)
        print(f"Total rearranged samples: {total_rearr}/{total_samples}")

    print("\nDone.")


if __name__ == '__main__':
    main()
