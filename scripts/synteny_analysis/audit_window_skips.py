#!/usr/bin/env python3
"""
DEPRECATED for screening (2026-05-27): Window skip audit from L2 BLAST results.

Superseded by aligned_fraction QC + paf_local_anomaly.py (v4.1).
Shared helpers (canonical_sample_id, resolve_sample_fastas) remain imported by active scripts.
See docs/.archive/2026-06/WINDOW_PIPELINE_DEPRECATED.md.

Audit skipped BLAST windows from window rearrangement test results.

Computes per-sample/per-species skip counts and correlates skipped windows
with variants in all_combined_data.csv (reference coordinates).
"""

import argparse
import sys
import warnings
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

try:
    from Bio import SeqIO
except ImportError:
    sys.exit("ERROR: BioPython required. pip install biopython")

from window_rearrangement_test import SPECIES_LIST, split_into_windows


def canonical_sample_id(sample_id, species):
    """Align sample IDs across FASTA, mapping TSV, and variant table."""
    sid = str(sample_id)
    if species in ('Glycine_max', 'Glycine_soja'):
        return sid.lstrip('0') or '0'
    return sid


def resolve_sample_fastas(species, sample_dir):
    """Return {sample_id: fasta_path} matching window_rearrangement_test logic."""
    sample_dir = Path(sample_dir)
    sample_species_dir = sample_dir / species
    valid_samples = None

    if not sample_species_dir.is_dir():
        glycine_dir = sample_dir / 'Glycine'
        if glycine_dir.is_dir() and species in ('Glycine_max', 'Glycine_soja'):
            sample_txt = sample_dir / species / 'sample.txt'
            if sample_txt.exists():
                with open(sample_txt) as f:
                    valid_samples = {line.strip() for line in f if line.strip()}
            sample_species_dir = glycine_dir
        else:
            return {}

    fastas = sorted(sample_species_dir.glob('*.fasta')) + sorted(sample_species_dir.glob('*.fa'))
    result = {}
    for fa in fastas:
        sid = canonical_sample_id(fa.stem.replace('.single_IR', ''), species)
        if valid_samples is not None:
            valid = {canonical_sample_id(v, species) for v in valid_samples}
            if sid not in valid:
                continue
        result[sid] = fa
    return result


def count_windows_for_sample(fasta_path, window_size):
    seq = str(next(SeqIO.parse(fasta_path, 'fasta')).seq)
    return len(split_into_windows(seq, window_size))


def load_mapping_tsv(path):
    df = pd.read_csv(path, sep='\t')
    df['query_win'] = df['query_win'].astype(int)
    df['ref_win'] = df['ref_win'].astype(int)
    return df


def build_consensus_ref_win(mapping_df, exclude_samples=None):
    """Median ref_win per query_win, optionally excluding outlier samples."""
    exclude = set(exclude_samples or [])
    sub = mapping_df[~mapping_df['sample_id'].isin(exclude)]
    if sub.empty:
        sub = mapping_df
    consensus = sub.groupby('query_win')['ref_win'].median().astype(int)
    confidence = {}
    for q, grp in sub.groupby('query_win'):
        vals = grp['ref_win'].values
        confidence[q] = 'high' if len(np.unique(vals)) == 1 else 'moderate'
    return consensus, confidence


def load_rearrangement_samples(results_dir, species):
    path = Path(results_dir) / f'{species}_rearrangement_details.tsv'
    if not path.exists():
        return set()
    df = pd.read_csv(path, sep='\t')
    return {canonical_sample_id(s, species) for s in df['sample_id'].astype(str)}


def load_variant_data(path):
    df = pd.read_csv(path, dtype=str)
    col_map = {
        df.columns[0]: 'sample_id',
        'V2': 'species', 'V3': 'position', 'V4': 'var_type',
    }
    df = df.rename(columns=col_map)
    df['position'] = pd.to_numeric(df['position'], errors='coerce')
    df = df[['sample_id', 'species', 'position', 'var_type']].dropna(subset=['position'])
    df['var_type'] = df['var_type'].str.lower()
    return df


def variant_lookup_key(sample_id, species):
    return (canonical_sample_id(sample_id, species), species)


class VariantIndex:
    """Fast interval queries on per-sample variant positions."""

    EMPTY = {
        'n_snp': 0, 'n_ins': 0, 'n_del': 0, 'n_complex': 0, 'n_mnp': 0,
        'n_indel': 0, 'has_indel': False, 'has_any_variant': False,
    }

    def __init__(self, variant_df):
        self._data = {}
        for (sid, sp), grp in variant_df.groupby(['sample_id', 'species'], sort=False):
            key = variant_lookup_key(sid, sp)
            order = np.argsort(grp['position'].values)
            self._data[key] = {
                'pos': grp['position'].values[order].astype(np.int64),
                'types': grp['var_type'].values[order],
            }

    def count_interval(self, sample_id, species, ref_start, ref_end):
        entry = self._data.get(variant_lookup_key(sample_id, species))
        if entry is None:
            return dict(self.EMPTY)
        pos = entry['pos']
        types = entry['types']
        i0 = int(np.searchsorted(pos, ref_start, side='left'))
        i1 = int(np.searchsorted(pos, ref_end, side='left'))
        sub = types[i0:i1]
        counts = {
            'n_snp': int((sub == 'snp').sum()),
            'n_ins': int((sub == 'ins').sum()),
            'n_del': int((sub == 'del').sum()),
            'n_complex': int((sub == 'complex').sum()),
            'n_mnp': int((sub == 'mnp').sum()),
        }
        counts['n_indel'] = counts['n_ins'] + counts['n_del']
        counts['has_indel'] = counts['n_indel'] > 0
        counts['has_any_variant'] = sum(counts[k] for k in counts if k.startswith('n_')) > 0
        return counts


def infer_skip_category(row, indel_p75_mapped):
    if row.get('is_gap_context'):
        return 'mapping_gap_context'
    if row['n_indel'] >= 1 or row['n_indel'] > indel_p75_mapped:
        return 'variant_indel_rich'
    if row['has_any_variant'] and not row['has_indel']:
        return 'variant_snp_only'
    return 'no_called_variant'


def mark_gap_context(skipped_wins):
    """Mark query_win in runs of >=3 consecutive skipped windows."""
    if not skipped_wins:
        return set()
    wins = sorted(skipped_wins)
    gap_set = set()
    run_start = 0
    for i in range(1, len(wins) + 1):
        if i == len(wins) or wins[i] != wins[i - 1] + 1:
            run_len = i - run_start
            if run_len >= 3:
                gap_set.update(wins[run_start:i])
            run_start = i
    return gap_set


def process_species(species, results_dir, sample_dir, variant_index, window_size, output_dir):
    mapping_path = Path(results_dir) / f'{species}_window_mapping.tsv'
    if not mapping_path.exists():
        print(f"  SKIP {species}: no mapping file")
        return None, None, None, None

    mapping_df = load_mapping_tsv(mapping_path)
    mapping_df['sample_id'] = mapping_df['sample_id'].map(
        lambda s: canonical_sample_id(s, species)
    )
    rearranged = load_rearrangement_samples(results_dir, species)
    sample_fastas = resolve_sample_fastas(species, sample_dir)
    if not sample_fastas:
        print(f"  SKIP {species}: no sample FASTAs")
        return None, None, None, None

    mapped_by_sample = (
        mapping_df.groupby('sample_id')['query_win']
        .apply(lambda s: set(s.unique()))
        .to_dict()
    )

    sample_map_lookup = {}
    for sid, grp in mapping_df.groupby('sample_id'):
        sample_map_lookup[sid] = dict(zip(grp['query_win'], grp['ref_win']))

    consensus, confidence = build_consensus_ref_win(mapping_df, exclude_samples=rearranged)

    sample_rows = []
    attribution_rows = []
    mapped_variant_rows = []

    for sample_id, fasta_path in sorted(sample_fastas.items()):
        total = count_windows_for_sample(fasta_path, window_size)
        mapped_set = mapped_by_sample.get(sample_id, set())
        mapped_n = len(mapped_set)
        skipped_set = set(range(total)) - mapped_set
        skipped_n = len(skipped_set)
        gap_ctx = mark_gap_context(skipped_set)

        sample_rows.append({
            'species': species,
            'sample_id': sample_id,
            'total_windows': total,
            'mapped_windows': mapped_n,
            'skipped_windows': skipped_n,
            'skip_rate': round(skipped_n / total, 6) if total else np.nan,
            'collinearity_qc_failed': mark_collinearity_qc(mapped_n),
            'collinearity_qc_status': (
                'collinearity_QC_failed' if mapped_n == 0 else 'pass'
            ),
            'is_rearrangement_flagged': sample_id in rearranged,
        })

        for q_win in sorted(skipped_set):
            ref_win = int(consensus.get(q_win, q_win))
            conf = confidence.get(q_win, 'low' if q_win not in consensus.index else 'fallback_query')
            ref_start = ref_win * window_size
            ref_end = (ref_win + 1) * window_size
            vc = variant_index.count_interval(sample_id, species, ref_start, ref_end)
            attribution_rows.append({
                'species': species,
                'sample_id': sample_id,
                'query_win': q_win,
                'consensus_ref_win': ref_win,
                'ref_start': ref_start,
                'ref_end': ref_end,
                'mapping_confidence': conf if q_win in consensus.index else 'low',
                'window_status': 'skipped',
                'collinearity_qc_failed': mapped_n == 0,
                'is_rearrangement_flagged': sample_id in rearranged,
                'is_gap_context': q_win in gap_ctx,
                **vc,
            })

        for q_win in sorted(mapped_set):
            ref_win = int(sample_map_lookup[sample_id][q_win])
            ref_start = ref_win * window_size
            ref_end = (ref_win + 1) * window_size
            vc = variant_index.count_interval(sample_id, species, ref_start, ref_end)
            mapped_variant_rows.append({
                'species': species,
                'sample_id': sample_id,
                'query_win': q_win,
                'ref_win': ref_win,
                'ref_start': ref_start,
                'ref_end': ref_end,
                'window_status': 'mapped',
                **vc,
            })

    return sample_rows, attribution_rows, mapped_variant_rows, rearranged


def mark_collinearity_qc(mapped_windows):
    """Samples with zero mapped windows cannot support macro-rearrangement inference."""
    return mapped_windows == 0


def aggregate_species_summary(sample_df):
    """Build species table with all-sample and QC-pass-only aggregates."""
    qc_pass = sample_df[~sample_df['collinearity_qc_failed']]

    species_all = sample_df.groupby('species').agg(
        n_samples_all=('sample_id', 'count'),
        n_qc_failed=('collinearity_qc_failed', 'sum'),
        mean_total=('total_windows', 'mean'),
        mean_mapped=('mapped_windows', 'mean'),
        mean_skipped=('skipped_windows', 'mean'),
        median_skipped=('skipped_windows', 'median'),
        std_skipped=('skipped_windows', 'std'),
        mean_skip_rate=('skip_rate', 'mean'),
        max_skipped=('skipped_windows', 'max'),
    ).reset_index()

    max_samples = sample_df.loc[sample_df.groupby('species')['skipped_windows'].idxmax()]
    species_all = species_all.merge(
        max_samples[['species', 'sample_id']].rename(columns={'sample_id': 'max_skipped_sample'}),
        on='species', how='left',
    )

    species_pass = qc_pass.groupby('species').agg(
        n_samples_qc_pass=('sample_id', 'count'),
        mean_skipped_qc_pass=('skipped_windows', 'mean'),
        median_skipped_qc_pass=('skipped_windows', 'median'),
        mean_skip_rate_qc_pass=('skip_rate', 'mean'),
        n_rearrangement_flagged_qc_pass=('is_rearrangement_flagged', 'sum'),
    ).reset_index()

    species_df = species_all.merge(species_pass, on='species', how='left')
    round_cols = [
        'mean_total', 'mean_mapped', 'mean_skipped', 'median_skipped', 'std_skipped',
        'mean_skip_rate', 'mean_skipped_qc_pass', 'median_skipped_qc_pass', 'mean_skip_rate_qc_pass',
    ]
    for col in round_cols:
        if col in species_df.columns:
            species_df[col] = species_df[col].round(4)
    species_df['n_qc_failed'] = species_df['n_qc_failed'].fillna(0).astype(int)
    species_df['n_samples_qc_pass'] = species_df['n_samples_qc_pass'].fillna(0).astype(int)
    species_df['n_rearrangement_flagged_qc_pass'] = (
        species_df['n_rearrangement_flagged_qc_pass'].fillna(0).astype(int)
    )
    return species_df, qc_pass


def main():
    parser = argparse.ArgumentParser(description='Audit skipped windows in rearrangement test results')
    parser.add_argument('--results-dir', required=True)
    parser.add_argument('--sample-dir', required=True)
    parser.add_argument('--variant-data', required=True)
    parser.add_argument('--window-size', type=int, default=100)
    parser.add_argument('--output-dir', required=True,
                        help='Output directory (e.g. win100 root or skip_audit/)')
    parser.add_argument('--species', nargs='*', default=None)
    args = parser.parse_args()

    warnings.warn(
        "audit_window_skips.py CLI is DEPRECATED for screening (v4.1). "
        "Use aligned_fraction QC + paf_local_anomaly.py instead. "
        "See docs/.archive/2026-06/WINDOW_PIPELINE_DEPRECATED.md.",
        DeprecationWarning,
        stacklevel=1,
    )

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print('Loading variant data...', flush=True)
    variant_df = load_variant_data(args.variant_data)
    print(f'  {len(variant_df)} variant records', flush=True)
    variant_index = VariantIndex(variant_df)
    del variant_df

    species_list = args.species if args.species else SPECIES_LIST
    all_sample_rows = []
    all_attribution = []
    all_mapped_var = []

    print(f'Processing {len(species_list)} species (window_size={args.window_size})...', flush=True)
    for species in species_list:
        print(f'  {species}...', flush=True)
        sample_rows, attr_rows, mapped_rows, _ = process_species(
            species, args.results_dir, args.sample_dir, variant_index, args.window_size, out_dir
        )
        if not sample_rows:
            continue
        all_sample_rows.extend(sample_rows)
        all_attribution.extend(attr_rows)
        all_mapped_var.extend(mapped_rows)

    sample_df = pd.DataFrame(all_sample_rows)
    sample_path = out_dir / 'window_skip_summary_by_sample.tsv'
    sample_df.to_csv(sample_path, sep='\t', index=False)
    print(f'Saved {sample_path} ({len(sample_df)} rows)')

    qc_failed_df = sample_df[sample_df['collinearity_qc_failed']].copy()
    qc_failed_path = out_dir / 'collinearity_qc_failed_samples.tsv'
    qc_failed_df.to_csv(qc_failed_path, sep='\t', index=False)
    print(f'Saved {qc_failed_path} ({len(qc_failed_df)} QC-failed samples)')

    species_df, qc_pass_df = aggregate_species_summary(sample_df)
    species_path = out_dir / 'window_skip_summary_by_species.tsv'
    species_df.to_csv(species_path, sep='\t', index=False)
    print(f'Saved {species_path}')

    attr_df = pd.DataFrame(all_attribution)
    if not attr_df.empty:
        indel_p75 = {}
        mapped_df = pd.DataFrame(all_mapped_var)
        qc_pass_ids = set(zip(qc_pass_df['species'], qc_pass_df['sample_id']))
        mapped_df['collinearity_qc_failed'] = mapped_df.apply(
            lambda r: (r['species'], r['sample_id']) not in qc_pass_ids, axis=1
        )
        attr_df['collinearity_qc_failed'] = attr_df['collinearity_qc_failed'].astype(bool)

        for sp in mapped_df['species'].unique():
            mp_sp = mapped_df[
                (~mapped_df['collinearity_qc_failed']) & (mapped_df['species'] == sp)
            ]
            if not mp_sp.empty:
                indel_p75[sp] = mp_sp['n_indel'].quantile(0.75)
        attr_df['inferred_skip_category'] = attr_df.apply(
            lambda r: infer_skip_category(r, indel_p75.get(r['species'], 0)), axis=1
        )
        attr_path = out_dir / 'window_skip_variant_attribution.tsv'
        attr_df.to_csv(attr_path, sep='\t', index=False)
        print(f'Saved {attr_path} ({len(attr_df)} rows)')

        comparison_rows = []
        for sp in sample_df['species'].unique():
            sk = attr_df[(attr_df['species'] == sp) & (~attr_df['collinearity_qc_failed'])]
            mp = mapped_df[(mapped_df['species'] == sp) & (~mapped_df['collinearity_qc_failed'])]
            if sk.empty or mp.empty:
                continue
            comparison_rows.append({
                'species': sp,
                'n_samples_qc_pass': int(qc_pass_df.loc[qc_pass_df['species'] == sp, 'sample_id'].count()),
                'n_qc_failed_excluded': int(qc_failed_df.loc[qc_failed_df['species'] == sp, 'sample_id'].count()),
                'n_skipped_windows': len(sk),
                'n_mapped_windows': len(mp),
                'skipped_mean_n_indel': round(sk['n_indel'].mean(), 4),
                'mapped_mean_n_indel': round(mp['n_indel'].mean(), 4),
                'skipped_mean_n_snp': round(sk['n_snp'].mean(), 4),
                'mapped_mean_n_snp': round(mp['n_snp'].mean(), 4),
                'skipped_pct_has_indel': round(100 * sk['has_indel'].mean(), 2),
                'mapped_pct_has_indel': round(100 * mp['has_indel'].mean(), 2),
                'skipped_pct_has_any_variant': round(100 * sk['has_any_variant'].mean(), 2),
                'mapped_pct_has_any_variant': round(100 * mp['has_any_variant'].mean(), 2),
                'skipped_pct_gap_context': round(100 * sk['is_gap_context'].mean(), 2),
                'skip_category_counts': sk['inferred_skip_category'].value_counts().to_dict(),
            })
        cmp_df = pd.DataFrame(comparison_rows)
        cmp_path = out_dir / 'window_skip_vs_mapped_comparison.tsv'
        cmp_df.to_csv(cmp_path, sep='\t', index=False)
        print(f'Saved {cmp_path}')
    else:
        print('No skipped windows found.')

    print('Done.')


if __name__ == '__main__':
    main()
