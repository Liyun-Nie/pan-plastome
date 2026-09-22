#!/usr/bin/env python3
"""
Tajima's D Neutrality Test for Plastome CDS Genes
===================================================
Calculates per-gene Tajima's D from curated variant data (all_combined_data.csv).
Designed for haploid chloroplast genome data: n = number of accessions (not 2*N).

Based on: Tajima F (1989) Genetics 123:585-595
Reference implementation: SNPGenie Tajima_D.R (Chase W. Nelson)

Usage:
    python selection_test_tajima.py \
        --variant-data <all_combined_data.csv> \
        --sample-dir <00-rawfa/single_IR/> \
        --cds-lengths <output_cds_lengths.csv> \
        --output-dir <selection_test_results/tajima_d/>
"""

import argparse
import os
import sys
import math
import warnings
from pathlib import Path
from collections import defaultdict

import pandas as pd
import numpy as np

from gene_function_utils import (
    DEFAULT_GENE_FUNCTION_MAP,
    get_functional_category,
    load_gene_function_map,
)

# ---------------------------------------------------------------------------
# Tajima's D core computation (ported from SNPGenie Tajima_D.R)
# ---------------------------------------------------------------------------

def compute_tajima_d(n, S, pi):
    """
    Compute Tajima's D statistic and significance level.

    Parameters
    ----------
    n : int   - number of sequences (= number of accessions for haploid plastome)
    S : int   - number of segregating sites
    pi : float - nucleotide diversity (mean pairwise differences)

    Returns
    -------
    D : float        - Tajima's D value
    sig_level : str  - significance level ('n.s.', 0.1, 0.05, 0.01, 0.001)
    """
    if n < 4 or S < 1:
        return np.nan, 'NA'

    a1 = sum(1.0 / i for i in range(1, n))
    a2 = sum(1.0 / (i ** 2) for i in range(1, n))

    b1 = (n + 1) / (3.0 * (n - 1))
    b2 = 2.0 * (n ** 2 + n + 3) / (9.0 * n * (n - 1))

    c1 = b1 - 1.0 / a1
    c2 = b2 - (n + 2) / (a1 * n) + a2 / (a1 ** 2)

    e1 = c1 / a1
    e2 = c2 / (a1 ** 2 + a2)

    V = e1 * S + e2 * S * (S - 1)
    if V <= 0:
        return np.nan, 'NA'

    D = (pi - S / a1) / math.sqrt(V)

    sig_level = _get_significance(n, D)
    return D, sig_level


# Critical values from Tajima 1989 Table 2 (subset for common sample sizes)
_CRITICAL_VALUES = [
    (1000, -2.369, 3.722, -2.062, 2.887, -1.715, 2.150, -1.505, 1.772),
    (100,  -2.495, 3.336, -2.160, 2.704, -1.781, 2.073, -1.555, 1.735),
    (50,   -2.505, 3.212, -2.178, 2.627, -1.800, 2.044, -1.570, 1.723),
    (30,   -2.478, 3.073, -2.173, 2.559, -1.807, 2.020, -1.580, 1.714),
    (20,   -2.414, 2.939, -2.146, 2.496, -1.803, 2.001, -1.584, 1.710),
    (15,   -2.329, 2.811, -2.103, 2.436, -1.791, 1.984, -1.584, 1.708),
    (10,   -2.105, 2.640, -1.967, 2.362, -1.733, 1.975, -1.559, 1.719),
    (8,    -1.909, 2.524, -1.830, 2.313, -1.663, 1.975, -1.522, 1.736),
    (6,    -1.556, 2.373, -1.540, 2.255, -1.478, 1.999, -1.405, 1.786),
    (5,    -1.276, 1.913, -1.275, 1.901, -1.269, 1.834, -1.255, 1.737),
    (4,    -0.876, 2.336, -0.876, 2.324, -0.876, 2.232, -0.876, 2.081),
]


def _get_significance(n, D):
    """Look up significance from Tajima 1989 critical value table."""
    for min_n, lo001, hi001, lo01, hi01, lo05, hi05, lo10, hi10 in _CRITICAL_VALUES:
        if n >= min_n:
            if D <= lo001 or D >= hi001:
                return 0.001
            elif D <= lo01 or D >= hi01:
                return 0.01
            elif D <= lo05 or D >= hi05:
                return 0.05
            elif D <= lo10 or D >= hi10:
                return 0.1
            else:
                return 'n.s.'
    return 'n.s.'


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_variant_data(filepath):
    """Load all_combined_data.csv and standardize column names."""
    df = pd.read_csv(filepath, dtype=str)
    col_map = {
        df.columns[0]: 'sample_id',
        'V2': 'species', 'V3': 'position', 'V4': 'var_type',
        'V5': 'ref_allele', 'V6': 'alt_allele', 'V7': 'evidence',
        'V8': 'region_type', 'V9': 'strand', 'V10': 'nt_pos',
        'V11': 'aa_pos', 'V12': 'effect', 'V13': 'locus_tag',
        'V14': 'gene', 'V15': 'product',
    }
    df = df.rename(columns=col_map)
    df['position'] = pd.to_numeric(df['position'], errors='coerce')
    return df


def load_sample_counts(sample_dir):
    """Load sample.txt for each species to get total n."""
    counts = {}
    sample_dir = Path(sample_dir)
    for species_dir in sample_dir.iterdir():
        if not species_dir.is_dir():
            continue
        sample_file = species_dir / 'sample.txt'
        if sample_file.exists():
            with open(sample_file) as f:
                n = sum(1 for line in f if line.strip())
            species_name = species_dir.name
            if species_name == 'Glycine':
                continue
            counts[species_name] = n
    return counts


def load_cds_lengths(filepath):
    """Load CDS length data."""
    df = pd.read_csv(filepath)
    lengths = {}
    for _, row in df.iterrows():
        key = (row['species'], row['Gene'])
        lengths[key] = int(row['CDS_length'])
    return lengths


# ---------------------------------------------------------------------------
# Tajima's D calculation pipeline
# ---------------------------------------------------------------------------

def calculate_tajima_d_per_gene(variant_df, sample_counts, cds_lengths,
                                category_map, min_S=3, min_n=4):
    """
    Calculate Tajima's D for each (species, gene) combination.

    For each segregating site, j = number of samples carrying the derived allele.
    Pi = sum over all sites of [ 2*j*(n-j) / (n*(n-1)) ]
    """
    cds_snps = variant_df[
        (variant_df['var_type'] == 'snp') &
        (variant_df['region_type'].str.upper() == 'CDS') &
        (variant_df['gene'].notna()) &
        (variant_df['gene'] != '')
    ].copy()

    results = []
    skipped = []

    for (species, gene), group in cds_snps.groupby(['species', 'gene']):
        n = sample_counts.get(species)
        if n is None or n < min_n:
            skipped.append({
                'species': species, 'gene': gene,
                'reason': f'n={n} < {min_n}' if n else 'no sample.txt'
            })
            continue

        site_counts = group.groupby('position')['sample_id'].nunique()
        S = len(site_counts)

        if S < min_S:
            skipped.append({
                'species': species, 'gene': gene,
                'reason': f'S={S} < {min_S}'
            })
            continue

        pi = 0.0
        for pos, j in site_counts.items():
            if j > n:
                j = n
            pi += 2.0 * j * (n - j) / (n * (n - 1))

        D, sig = compute_tajima_d(n, S, pi)

        a1 = sum(1.0 / i for i in range(1, n))
        theta_w = S / a1

        gene_length = cds_lengths.get((species, gene), np.nan)
        pi_per_site = pi / gene_length if gene_length and gene_length > 0 else np.nan
        theta_w_per_site = theta_w / gene_length if gene_length and gene_length > 0 else np.nan

        results.append({
            'species': species,
            'gene': gene,
            'n': n,
            'S': S,
            'pi': round(pi, 6),
            'theta_w': round(theta_w, 6),
            'tajima_d': round(D, 6) if not np.isnan(D) else np.nan,
            'sig_level': sig,
            'gene_length': gene_length,
            'pi_per_site': round(pi_per_site, 8) if not np.isnan(pi_per_site) else np.nan,
            'theta_w_per_site': round(theta_w_per_site, 8) if not np.isnan(theta_w_per_site) else np.nan,
            'functional_category': get_functional_category(gene, category_map),
        })

    results_df = pd.DataFrame(results)
    skipped_df = pd.DataFrame(skipped)
    return results_df, skipped_df


# ---------------------------------------------------------------------------
# Heterogeneity analysis
# ---------------------------------------------------------------------------

def analyze_heterogeneity(results_df):
    """Test whether Tajima's D differs across functional categories."""
    from scipy import stats

    valid = results_df.dropna(subset=['tajima_d'])
    if valid.empty:
        return pd.DataFrame()

    groups = valid.groupby('functional_category')['tajima_d']

    summary_rows = []
    for cat, vals in groups:
        summary_rows.append({
            'functional_category': cat,
            'n_genes': len(vals),
            'mean_D': round(vals.mean(), 4),
            'median_D': round(vals.median(), 4),
            'std_D': round(vals.std(), 4),
            'min_D': round(vals.min(), 4),
            'max_D': round(vals.max(), 4),
        })

    cat_groups = [g.values for _, g in groups]
    if len(cat_groups) >= 2:
        try:
            h_stat, p_val = stats.kruskal(*cat_groups)
        except ValueError:
            h_stat, p_val = np.nan, np.nan
        for row in summary_rows:
            row['kruskal_wallis_H'] = round(h_stat, 4) if not np.isnan(h_stat) else np.nan
            row['kruskal_wallis_p'] = p_val

    return pd.DataFrame(summary_rows)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Tajima's D neutrality test for plastome CDS genes")
    parser.add_argument('--variant-data', required=True,
                        help='Path to all_combined_data.csv')
    parser.add_argument('--sample-dir', required=True,
                        help='Directory containing {Species}/sample.txt')
    parser.add_argument('--cds-lengths', required=True,
                        help='Path to output_cds_lengths.csv')
    parser.add_argument('--output-dir', default='selection_test_results/tajima_d',
                        help='Output directory')
    parser.add_argument('--min-S', type=int, default=3,
                        help='Minimum segregating sites per gene (default: 3)')
    parser.add_argument('--min-n', type=int, default=4,
                        help='Minimum sample size per species (default: 4)')
    parser.add_argument('--gene-function-map', default=str(DEFAULT_GENE_FUNCTION_MAP),
                        help='Path to cpopvar gene_function_map.csv')
    args = parser.parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print("Loading variant data...")
    variant_df = load_variant_data(args.variant_data)
    print(f"  Loaded {len(variant_df)} variant records")

    print("Loading sample counts...")
    sample_counts = load_sample_counts(args.sample_dir)
    print(f"  Found {len(sample_counts)} species: " +
          ", ".join(f"{sp}(n={n})" for sp, n in sorted(sample_counts.items())))

    print("Loading CDS lengths...")
    cds_lengths = load_cds_lengths(args.cds_lengths)
    print(f"  Loaded lengths for {len(cds_lengths)} (species, gene) pairs")

    print("Loading gene function map...")
    category_map = load_gene_function_map(args.gene_function_map)
    print(f"  Loaded {len(category_map)} gene annotations")

    print(f"Calculating Tajima's D (min_S={args.min_S}, min_n={args.min_n})...")
    results_df, skipped_df = calculate_tajima_d_per_gene(
        variant_df, sample_counts, cds_lengths, category_map,
        min_S=args.min_S, min_n=args.min_n
    )
    print(f"  Computed D for {len(results_df)} (species, gene) pairs")
    print(f"  Skipped {len(skipped_df)} pairs")

    # Per-species output
    for species in results_df['species'].unique():
        sp_df = results_df[results_df['species'] == species]
        sp_file = out_dir / f"{species}_tajima_d_per_gene.tsv"
        sp_df.to_csv(sp_file, sep='\t', index=False)

    # Summary output
    summary_file = out_dir / 'all_species_tajima_d_summary.tsv'
    results_df.to_csv(summary_file, sep='\t', index=False)
    print(f"  Saved summary to {summary_file}")

    skipped_file = out_dir / 'skipped_genes.tsv'
    skipped_df.to_csv(skipped_file, sep='\t', index=False)

    # Heterogeneity test
    print("Running heterogeneity analysis...")
    het_df = analyze_heterogeneity(results_df)
    if not het_df.empty:
        het_file = out_dir / 'tajima_d_heterogeneity_test.tsv'
        het_df.to_csv(het_file, sep='\t', index=False)
        print(f"  Saved heterogeneity test to {het_file}")
        print("\n  Functional category summary:")
        for _, row in het_df.iterrows():
            print(f"    {row['functional_category']}: "
                  f"mean_D={row['mean_D']}, n={row['n_genes']}")
        if 'kruskal_wallis_p' in het_df.columns:
            p = het_df['kruskal_wallis_p'].iloc[0]
            print(f"  Kruskal-Wallis p-value: {p:.6f}" if not np.isnan(p) else
                  "  Kruskal-Wallis: could not compute")

    # Quick stats
    valid = results_df.dropna(subset=['tajima_d'])
    if not valid.empty:
        print(f"\n  Overall Tajima's D statistics:")
        print(f"    Mean: {valid['tajima_d'].mean():.4f}")
        print(f"    Median: {valid['tajima_d'].median():.4f}")
        print(f"    Range: [{valid['tajima_d'].min():.4f}, {valid['tajima_d'].max():.4f}]")
        sig_count = valid[valid['sig_level'].apply(
            lambda x: isinstance(x, float) and x <= 0.05)].shape[0]
        print(f"    Significant (p<=0.05): {sig_count}/{len(valid)} "
              f"({100*sig_count/len(valid):.1f}%)")

    print("\nDone.")


if __name__ == '__main__':
    main()
