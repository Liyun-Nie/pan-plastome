#!/usr/bin/env python3
"""
Tajima's D Neutrality Test for Plastome IGS (Intergenic Spacer) Regions
========================================================================
Extends the CDS-based `selection_test_tajima.py` analysis to intergenic
spacer (IGS) regions, in order to test whether cultivated species carry a
stronger signature of population expansion following a domestication
bottleneck (more negative Tajima's D) than their wild counterparts.

Rationale
---------
CDS Tajima's D is confounded by purifying/background selection acting on
coding sequence, which biases D negative regardless of demographic history
and can mask a demographic (bottleneck -> expansion) signal. IGS regions
are not fully neutral (some harbour regulatory elements / are subject to
repeat-driven processes), but are the best available lower-selection-
constraint proxy in this dataset for isolating demographic signal from
locus-specific selection, and are already used elsewhere in the manuscript
as the main "less constrained" plastomic feature (Fig. 1g-1j).

Because the plastome is a single, non-recombining locus, the
population-genetically appropriate way to test a genome-wide demographic
signal (bottleneck/expansion) is to pool ALL segregating sites of a given
region type (IGS or CDS) into ONE Tajima's D value per species ("genome-
wide" D), rather than averaging many per-locus D values (which are not
independent replicates of the same coalescent history and would
pseudoreplicate a demographic test). This script therefore reports THREE
tables, with genome-wide D treated as the PRIMARY metric for the wild vs.
cultivated comparison:

  1. Per-(species, IGS region) Tajima's D  -- exploratory / descriptive,
     mirrors the CDS per-gene table for structural symmetry. NOT
     independent replicates (linked, non-recombining genome) - do not use
     this table alone to test the domestication hypothesis.
  2. Genome-wide IGS Tajima's D (one value per species, all intergenic
     SNPs pooled)      -- PRIMARY metric for the domestication test.
  3. Genome-wide CDS Tajima's D (one value per species, all CDS SNPs
     pooled, ignoring gene boundaries) -- computed for a symmetric
     CDS-vs-IGS contrast per species. This is a NEW, different
     aggregation from the already-audited per-gene CDS v2 table
     (`all_species_tajima_d_summary_v2.tsv`) and must not be confused
     with, or used to replace, that table.
  4. TRUE whole-plastome genome-wide Tajima's D (one value per species,
     ALL SNPs pooled regardless of region type: CDS + IGS + intron + RNA),
     normalized by the single-IR reference length in
     `species_genome_regions.csv` (`total_length` column). This is the
     most literal reading of "plastome as a single non-recombining locus"
     and is reported alongside (2) and (3) for a three-way scope
     comparison (Plastome / CDS-only / IGS-only) per species.

Two comparison designs are reported for the wild vs. cultivated
hypothesis, from most to least rigorous:

  A. Congeneric paired contrast (PRIMARY, most rigorous): compares the
     genome-wide D of a cultivated species against its wild congener
     within the SAME genus, which controls for phylogeny, background
     mutation rate, and generation time. Only two such pairs exist in the
     current 17-species panel: Glycine_max (cultivar) vs Glycine_soja
     (wild), and Oryza_sativa (cultivar) vs Oryza_rufipogon (wild).
     (Solanum_brevicaule/candolleanum are BOTH wild, and
     Gossypium_hirsutum/barbadense are BOTH cultivated in this panel, so
     neither pair is usable for this specific contrast; they are reported
     separately for transparency but excluded from the domestication
     test.)
  B. Pooled cross-species contrast (SECONDARY, exploratory): Mann-Whitney
     U test (= two-sample Wilcoxon rank-sum test, consistent with the
     test already used in the main text for the Fig. 1j comparison) of
     genome-wide D across all 13 cultivated vs. 4 wild species. This
     design pools unrelated lineages and is confounded by phylogeny
     (each species is not an independent replicate of "domestication");
     it is reported only as a weaker, illustrative supplement to (A).

Scope / method notes (for Methods drafting)
--------------------------------------------
- Only SNPs (var_type == 'snp') are used, consistent with the CDS Tajima's
  D script and with the infinite-sites assumption underlying the Tajima's
  D formula. This differs from the main-text IGS variant-frequency metric
  (Fig. 1g-1j), which combines SNV + indel + CPX counts; the two analyses
  therefore use different numerators/denominators and should not be
  conflated in text.
- IGS region boundaries and lengths are taken from
  `output_gene_intergenic_intron_pos_length.csv` (Type == 'intergenic'),
  the SAME region_info source used for the main-text IGS length
  normalization, per user instruction.
- The Tajima's D core formula (compute_tajima_d) is intentionally
  DUPLICATED (not imported) from `selection_test_tajima.py`, byte-for-byte
  identical, so that CDS and IGS values are computed with exactly the same
  statistical method while leaving the already-audited CDS v2 pipeline
  file untouched. Any future formula fix must be applied in both files.

Usage
-----
    python selection_test_tajima_igs.py \
        --variant-data <all_combined_data.csv> \
        --sample-dir <00-rawfa/single_IR/> \
        --region-info <output_gene_intergenic_intron_pos_length.csv> \
        --cds-lengths <output_cds_lengths.csv> \
        --group-info <group_info.csv> \
        --output-dir <selection_test_results/tajima_d_igs/>
"""

import argparse
import bisect
import math
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore", category=FutureWarning)

# ---------------------------------------------------------------------------
# Tajima's D core computation (duplicated verbatim from selection_test_tajima.py
# -- see module docstring for rationale. DO NOT let this drift from the
# original; if you fix a bug here, fix it there too.)
# ---------------------------------------------------------------------------


def compute_tajima_d(n, S, pi):
    """Compute Tajima's D statistic and significance level.

    Parameters
    ----------
    n : int   - number of sequences (= number of accessions for haploid plastome)
    S : int   - number of segregating sites
    pi : float - nucleotide diversity (mean pairwise differences)
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
    """Load all_combined_data.csv and standardize column names (same mapping
    as selection_test_tajima.py)."""
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
    """Load sample.txt for each species to get total n (identical logic to
    selection_test_tajima.py)."""
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
    """Load per-(species,gene) CDS length; also return per-species total
    CDS length (sum over genes) for genome-wide CDS pooling."""
    df = pd.read_csv(filepath)
    lengths = {}
    total_by_species = {}
    for _, row in df.iterrows():
        key = (row['species'], row['Gene'])
        length = int(row['CDS_length'])
        lengths[key] = length
        total_by_species[row['species']] = total_by_species.get(row['species'], 0) + length
    return lengths, total_by_species


def load_igs_region_info(filepath):
    """Load output_gene_intergenic_intron_pos_length.csv, keep Type ==
    'intergenic' rows only. Returns:
      regions_by_species: {species: DataFrame[region_name, start, end, length]
                            sorted by start, index reset}
      total_igs_length: {species: total intergenic bp}
    """
    df = pd.read_csv(filepath)
    df = df[df['Type'] == 'intergenic'].copy()
    df = df.rename(columns={
        'Region_Name': 'region_name',
        'Start_Position': 'start',
        'End_Position': 'end',
        'Length': 'length',
    })
    df['start'] = pd.to_numeric(df['start'], errors='coerce')
    df['end'] = pd.to_numeric(df['end'], errors='coerce')
    df['length'] = pd.to_numeric(df['length'], errors='coerce')
    df = df.dropna(subset=['start', 'end', 'length'])

    regions_by_species = {}
    total_igs_length = {}
    for species, grp in df.groupby('species'):
        grp = grp.sort_values('start').reset_index(drop=True)
        regions_by_species[species] = grp
        total_igs_length[species] = grp['length'].sum()
    return regions_by_species, total_igs_length


def load_group_info(filepath):
    """species -> Status ('Wild' / 'Cultivar')."""
    df = pd.read_csv(filepath)
    return dict(zip(df['species'], df['Status']))


def load_genome_total_lengths(filepath):
    """Load species_genome_regions.csv -> {species: total_length}. This is
    the single-IR plastome reference length (one IR copy already removed
    prior to variant calling, per manuscript Methods), i.e. the correct
    denominator for a TRUE whole-plastome genome-wide Tajima's D that pools
    SNPs from ALL region types (CDS + IGS + intron + RNA) together."""
    df = pd.read_csv(filepath)
    df['total_length'] = pd.to_numeric(df['total_length'], errors='coerce')
    return dict(zip(df['species'], df['total_length']))


# ---------------------------------------------------------------------------
# Position -> IGS region assignment (mirrors calculate_igs_frequencies() in
# cpopvar/R/normalize_frequencies.R for methodological consistency with the
# main-text IGS variant-frequency pipeline)
# ---------------------------------------------------------------------------


def assign_igs_regions(igs_snp_df, regions_by_species):
    """For each row in igs_snp_df (species, position, ...), find the
    containing IGS region_name via interval search. Rows whose position
    does not fall inside any known interval (annotation/version mismatch)
    get region_name = NaN and are reported for QC, but ARE STILL RETAINED
    for genome-wide pooling (they are still valid intergenic SNPs per the
    variant caller's own annotation; only the per-region table drops them).
    """
    region_names = [None] * len(igs_snp_df)
    species_arr = igs_snp_df['species'].values
    pos_arr = igs_snp_df['position'].values

    cache = {}
    for i in range(len(igs_snp_df)):
        sp = species_arr[i]
        pos = pos_arr[i]
        if sp not in regions_by_species:
            continue
        if sp not in cache:
            df = regions_by_species[sp]
            cache[sp] = (df, df['start'].values, df['end'].values, df['region_name'].values)
        df, starts, ends, names = cache[sp]
        idx = bisect.bisect_right(starts, pos) - 1
        if idx >= 0 and starts[idx] <= pos <= ends[idx]:
            region_names[i] = names[idx]

    igs_snp_df = igs_snp_df.copy()
    igs_snp_df['igs_region'] = region_names
    return igs_snp_df


# ---------------------------------------------------------------------------
# Tajima's D calculation: per-region (exploratory) and genome-wide (primary)
# ---------------------------------------------------------------------------


def calculate_tajima_d_per_region(igs_snp_df, sample_counts, region_lengths_lookup,
                                   min_S=3, min_n=4):
    """Per-(species, igs_region) Tajima's D -- structural analogue of the
    CDS per-gene table. NOT independent replicates; exploratory only."""
    results = []
    skipped = []

    mapped = igs_snp_df.dropna(subset=['igs_region'])

    for (species, region), group in mapped.groupby(['species', 'igs_region']):
        n = sample_counts.get(species)
        if n is None or n < min_n:
            skipped.append({'species': species, 'igs_region': region,
                             'reason': 'n={} < {}'.format(n, min_n) if n else 'no sample.txt'})
            continue

        site_counts = group.groupby('position')['sample_id'].nunique()
        S = len(site_counts)
        if S < min_S:
            skipped.append({'species': species, 'igs_region': region,
                             'reason': 'S={} < {}'.format(S, min_S)})
            continue

        pi = 0.0
        for pos, j in site_counts.items():
            if j > n:
                j = n
            pi += 2.0 * j * (n - j) / (n * (n - 1))

        D, sig = compute_tajima_d(n, S, pi)
        a1 = sum(1.0 / i for i in range(1, n))
        theta_w = S / a1

        region_length = region_lengths_lookup.get((species, region), np.nan)
        pi_per_site = pi / region_length if region_length and region_length > 0 else np.nan
        theta_w_per_site = theta_w / region_length if region_length and region_length > 0 else np.nan

        results.append({
            'species': species,
            'igs_region': region,
            'n': n,
            'S': S,
            'pi': round(pi, 6),
            'theta_w': round(theta_w, 6),
            'tajima_d': round(D, 6) if not np.isnan(D) else np.nan,
            'sig_level': sig,
            'region_length': region_length,
            'pi_per_site': round(pi_per_site, 8) if not np.isnan(pi_per_site) else np.nan,
            'theta_w_per_site': round(theta_w_per_site, 8) if not np.isnan(theta_w_per_site) else np.nan,
        })

    return pd.DataFrame(results), pd.DataFrame(skipped)


def calculate_tajima_d_genomewide(snp_df, sample_counts, total_length_by_species,
                                   min_S=3, min_n=4, region_label='region'):
    """Pool ALL segregating sites of one region type (IGS or CDS) into a
    SINGLE Tajima's D value per species. This is the population-genetically
    appropriate unit for a non-recombining, single-genealogy genome (the
    whole plastome / whole IGS complement / whole CDS complement behaves as
    one locus), and is the PRIMARY metric for the wild vs. cultivated
    demographic-signal comparison.
    """
    results = []
    skipped = []

    for species, group in snp_df.groupby('species'):
        n = sample_counts.get(species)
        if n is None or n < min_n:
            skipped.append({'species': species, 'region': region_label,
                             'reason': 'n={} < {}'.format(n, min_n) if n else 'no sample.txt'})
            continue

        site_counts = group.groupby('position')['sample_id'].nunique()
        S = len(site_counts)
        if S < min_S:
            skipped.append({'species': species, 'region': region_label,
                             'reason': 'S={} < {}'.format(S, min_S)})
            continue

        pi = 0.0
        for pos, j in site_counts.items():
            if j > n:
                j = n
            pi += 2.0 * j * (n - j) / (n * (n - 1))

        D, sig = compute_tajima_d(n, S, pi)
        a1 = sum(1.0 / i for i in range(1, n))
        theta_w = S / a1

        total_len = total_length_by_species.get(species, np.nan)
        pi_per_site = pi / total_len if total_len and total_len > 0 else np.nan
        theta_w_per_site = theta_w / total_len if total_len and total_len > 0 else np.nan

        results.append({
            'species': species,
            'region': region_label,
            'n': n,
            'S': S,
            'pi': round(pi, 6),
            'theta_w': round(theta_w, 6),
            'tajima_d': round(D, 6) if not np.isnan(D) else np.nan,
            'sig_level': sig,
            'total_length_bp': total_len,
            'pi_per_site': round(pi_per_site, 8) if not np.isnan(pi_per_site) else np.nan,
            'theta_w_per_site': round(theta_w_per_site, 8) if not np.isnan(theta_w_per_site) else np.nan,
        })

    return pd.DataFrame(results), pd.DataFrame(skipped)


# ---------------------------------------------------------------------------
# Wild vs. cultivated comparisons
# ---------------------------------------------------------------------------


def add_status(df, status_map, species_col='species'):
    df = df.copy()
    df['status'] = df[species_col].map(status_map)
    return df


CONGENERIC_DOMESTICATION_PAIRS = [
    ('Glycine_max', 'Glycine_soja', 'Glycine', True),
    ('Oryza_sativa', 'Oryza_rufipogon', 'Oryza', True),
    ('Gossypium_hirsutum', 'Gossypium_barbadense', 'Gossypium', False),   # both cultivated
    ('Solanum_brevicaule', 'Solanum_candolleanum', 'Solanum', False),    # both wild
]


def build_congeneric_pairs_table(plastome_gw_df, igs_gw_df, cds_gw_df):
    """Side-by-side genome-wide D (Plastome / IGS-only / CDS-only) for each
    congeneric pair. `usable_for_domestication_test` flags whether the pair
    is a true wild-vs-cultivated contrast in this panel."""
    plastome_lookup = plastome_gw_df.set_index('species')['tajima_d'].to_dict()
    igs_lookup = igs_gw_df.set_index('species')['tajima_d'].to_dict()
    cds_lookup = cds_gw_df.set_index('species')['tajima_d'].to_dict()
    rows = []
    for sp_a, sp_b, genus, usable in CONGENERIC_DOMESTICATION_PAIRS:
        rows.append({
            'genus': genus,
            'species_a': sp_a, 'status_a': 'Cultivar' if usable else '(see note)',
            'species_b': sp_b, 'status_b': 'Wild' if usable else '(see note)',
            'plastome_D_a': plastome_lookup.get(sp_a), 'plastome_D_b': plastome_lookup.get(sp_b),
            'plastome_D_diff_a_minus_b': (plastome_lookup.get(sp_a) - plastome_lookup.get(sp_b))
                if plastome_lookup.get(sp_a) is not None and plastome_lookup.get(sp_b) is not None else np.nan,
            'igs_D_a': igs_lookup.get(sp_a), 'igs_D_b': igs_lookup.get(sp_b),
            'igs_D_diff_a_minus_b': (igs_lookup.get(sp_a) - igs_lookup.get(sp_b))
                if igs_lookup.get(sp_a) is not None and igs_lookup.get(sp_b) is not None else np.nan,
            'cds_D_a': cds_lookup.get(sp_a), 'cds_D_b': cds_lookup.get(sp_b),
            'cds_D_diff_a_minus_b': (cds_lookup.get(sp_a) - cds_lookup.get(sp_b))
                if cds_lookup.get(sp_a) is not None and cds_lookup.get(sp_b) is not None else np.nan,
            'usable_for_domestication_test': usable,
            'note': '' if usable else 'both species share the same domestication status in this panel; excluded from the wild-vs-cultivated test',
        })
    return pd.DataFrame(rows)


def run_wild_vs_cultivar_tests(plastome_gw_df, igs_gw_df, cds_gw_df, igs_region_df, status_map):
    """Mann-Whitney U (= two-sample Wilcoxon rank-sum test, consistent with
    the main text) comparisons, species-level (primary, three region scopes)
    and per-region pooled (secondary/exploratory)."""
    from scipy import stats

    rows = []

    def _mwu(vals_a, vals_b, label_a, label_b, level):
        if len(vals_a) < 2 or len(vals_b) < 2:
            u, p = np.nan, np.nan
        else:
            try:
                u, p = stats.mannwhitneyu(vals_a, vals_b, alternative='two-sided')
            except ValueError:
                u, p = np.nan, np.nan
        return {
            'comparison_level': level,
            'group_a': label_a, 'n_a': len(vals_a), 'mean_a': np.mean(vals_a) if len(vals_a) else np.nan,
            'median_a': np.median(vals_a) if len(vals_a) else np.nan,
            'group_b': label_b, 'n_b': len(vals_b), 'mean_b': np.mean(vals_b) if len(vals_b) else np.nan,
            'median_b': np.median(vals_b) if len(vals_b) else np.nan,
            'mannwhitney_U': u, 'p_value': p,
            'direction_consistent_with_expansion_hypothesis': bool(
                len(vals_a) and len(vals_b) and np.mean(vals_a) < np.mean(vals_b)
            ),
        }

    plastome_gw = add_status(plastome_gw_df, status_map)
    igs_gw = add_status(igs_gw_df, status_map)
    cds_gw = add_status(cds_gw_df, status_map)

    cult_plastome = plastome_gw.loc[plastome_gw['status'] == 'Cultivar', 'tajima_d'].dropna().values
    wild_plastome = plastome_gw.loc[plastome_gw['status'] == 'Wild', 'tajima_d'].dropna().values
    rows.append(_mwu(cult_plastome, wild_plastome, 'Cultivar', 'Wild',
                      'species-level TRUE whole-plastome genome-wide D (n=17 species)'))

    cult_igs = igs_gw.loc[igs_gw['status'] == 'Cultivar', 'tajima_d'].dropna().values
    wild_igs = igs_gw.loc[igs_gw['status'] == 'Wild', 'tajima_d'].dropna().values
    rows.append(_mwu(cult_igs, wild_igs, 'Cultivar', 'Wild',
                      'species-level genome-wide IGS-only D (n=17 species)'))

    cult_cds = cds_gw.loc[cds_gw['status'] == 'Cultivar', 'tajima_d'].dropna().values
    wild_cds = cds_gw.loc[cds_gw['status'] == 'Wild', 'tajima_d'].dropna().values
    rows.append(_mwu(cult_cds, wild_cds, 'Cultivar', 'Wild',
                      'species-level genome-wide CDS-only D (n=17 species)'))

    if igs_region_df is not None and len(igs_region_df):
        region_with_status = add_status(igs_region_df, status_map)
        cult_region = region_with_status.loc[region_with_status['status'] == 'Cultivar', 'tajima_d'].dropna().values
        wild_region = region_with_status.loc[region_with_status['status'] == 'Wild', 'tajima_d'].dropna().values
        rows.append(_mwu(cult_region, wild_region, 'Cultivar', 'Wild',
                          'EXPLORATORY (pseudoreplicated, non-independent loci): per-IGS-region D pooled records'))

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Folded SFS by domestication status (IGS only)
# ---------------------------------------------------------------------------


def compute_folded_sfs_by_status(igs_snp_df, sample_counts, status_map, n_bins=10):
    """Pool folded site-frequency spectra of intergenic SNPs within each
    domestication-status group (all segregating sites from all species in
    the group, frequency = min(j, n-j)/n). A visual excess of low-frequency
    variants in cultivated vs. wild would support the expansion hypothesis."""
    bin_edges = np.linspace(0, 0.5, n_bins + 1)
    rows = []
    freqs_by_status = {'Wild': [], 'Cultivar': []}

    for species, group in igs_snp_df.groupby('species'):
        status = status_map.get(species)
        if status not in freqs_by_status:
            continue
        n = sample_counts.get(species)
        if not n:
            continue
        site_counts = group.groupby('position')['sample_id'].nunique()
        for j in site_counts.values:
            if j <= 0 or j >= n:
                continue
            folded = min(j, n - j) / float(n)
            freqs_by_status[status].append(folded)

    for status, freqs in freqs_by_status.items():
        freqs = np.array(freqs)
        total = len(freqs)
        for k in range(n_bins):
            lo, hi = bin_edges[k], bin_edges[k + 1]
            if k == n_bins - 1:
                count = int(((freqs >= lo) & (freqs <= hi)).sum())
            else:
                count = int(((freqs >= lo) & (freqs < hi)).sum())
            rows.append({
                'status': status,
                'freq_bin_left': round(lo, 4),
                'freq_bin_right': round(hi, 4),
                'freq_bin_label': '{:.2f}-{:.2f}'.format(lo, hi),
                'count': count,
                'proportion': (count / total) if total else np.nan,
                'total_sites_in_group': total,
            })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Plotting (lightweight; kept separate from visualize_tajima_d.py)
# ---------------------------------------------------------------------------


def make_plots(plastome_gw_df, igs_gw_df, cds_gw_df, sfs_df, pairs_df, status_map, out_dir):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig_dir = Path(out_dir) / 'figures'
    fig_dir.mkdir(parents=True, exist_ok=True)

    status_colors = {'Wild': '#2E7D32', 'Cultivar': '#C62828'}

    # 1. Genome-wide D boxplot: Plastome vs IGS vs CDS, colored by status
    fig, axes = plt.subplots(1, 3, figsize=(13, 5), sharey=True)
    panels = [
        (plastome_gw_df, 'Whole-plastome genome-wide D'),
        (igs_gw_df, 'IGS-only genome-wide D'),
        (cds_gw_df, 'CDS-only genome-wide D'),
    ]
    for ax, (df, title) in zip(axes, panels):
        d = add_status(df, status_map)
        data = [d.loc[d['status'] == s, 'tajima_d'].dropna().values for s in ('Wild', 'Cultivar')]
        bp = ax.boxplot(data, labels=['Wild (n={})'.format(len(data[0])), 'Cultivar (n={})'.format(len(data[1]))],
                         patch_artist=True, widths=0.5)
        for patch, s in zip(bp['boxes'], ('Wild', 'Cultivar')):
            patch.set_facecolor(status_colors[s])
            patch.set_alpha(0.5)
        for i, s in enumerate(('Wild', 'Cultivar')):
            vals = data[i]
            jitter = np.random.uniform(-0.08, 0.08, size=len(vals))
            ax.scatter(np.full(len(vals), i + 1) + jitter, vals, color=status_colors[s],
                       edgecolor='black', linewidth=0.5, zorder=3)
        ax.axhline(0, color='grey', linestyle='--', linewidth=0.8)
        ax.set_title(title)
        ax.set_ylabel("Tajima's D")
    fig.suptitle('Species-level genome-wide Tajima\'s D by domestication status')
    fig.tight_layout()
    fig.savefig(fig_dir / 'tajima_d_genomewide_boxplot_by_status.pdf')
    fig.savefig(fig_dir / 'tajima_d_genomewide_boxplot_by_status.png', dpi=200)
    plt.close(fig)

    # 2. Folded SFS by status
    if sfs_df is not None and len(sfs_df):
        fig, ax = plt.subplots(figsize=(7, 5))
        for status in ('Wild', 'Cultivar'):
            sub = sfs_df[sfs_df['status'] == status].sort_values('freq_bin_left')
            ax.plot(sub['freq_bin_label'], sub['proportion'], marker='o', label=status,
                    color=status_colors[status])
        ax.set_xlabel('Folded derived allele frequency bin')
        ax.set_ylabel('Proportion of intergenic segregating sites')
        ax.set_title('Pooled folded SFS of IGS SNPs by domestication status')
        ax.legend()
        plt.setp(ax.get_xticklabels(), rotation=45, ha='right')
        fig.tight_layout()
        fig.savefig(fig_dir / 'sfs_igs_by_status.pdf')
        fig.savefig(fig_dir / 'sfs_igs_by_status.png', dpi=200)
        plt.close(fig)

    # 3. Congeneric pair dumbbell plots, one per region scope (Plastome / IGS / CDS)
    usable = pairs_df[pairs_df['usable_for_domestication_test']]
    if len(usable):
        for scope, col_a, col_b, fname in [
            ('whole-plastome', 'plastome_D_a', 'plastome_D_b', 'congeneric_pairs_plastome_dumbbell'),
            ('IGS-only', 'igs_D_a', 'igs_D_b', 'congeneric_pairs_igs_dumbbell'),
            ('CDS-only', 'cds_D_a', 'cds_D_b', 'congeneric_pairs_cds_dumbbell'),
        ]:
            fig, ax = plt.subplots(figsize=(6, 3 + len(usable)))
            for i, row in usable.reset_index(drop=True).iterrows():
                y = i
                ax.plot([row[col_b], row[col_a]], [y, y], color='grey', zorder=1)
                ax.scatter(row[col_b], y, color=status_colors['Wild'], s=90, zorder=2,
                           label='Wild' if i == 0 else None)
                ax.scatter(row[col_a], y, color=status_colors['Cultivar'], s=90, zorder=2,
                           label='Cultivar' if i == 0 else None)
                ax.text(max(row[col_a], row[col_b]) + 0.05, y,
                        '{} ({} vs {})'.format(row['genus'], row['species_a'], row['species_b']),
                        va='center', fontsize=9)
            ax.axvline(0, color='grey', linestyle='--', linewidth=0.8)
            ax.set_yticks([])
            ax.set_xlabel("Genome-wide Tajima's D ({})".format(scope))
            ax.set_title('Congeneric wild-vs-cultivated pairs ({} D)'.format(scope))
            ax.legend()
            fig.tight_layout()
            fig.savefig(fig_dir / (fname + '.pdf'))
            fig.savefig(fig_dir / (fname + '.png'), dpi=200)
            plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(
        description="Tajima's D neutrality test for plastome IGS regions, "
                     "with a wild-vs-cultivated domestication-bottleneck test")
    parser.add_argument('--variant-data', required=True, help='Path to all_combined_data.csv')
    parser.add_argument('--sample-dir', required=True, help='Directory containing {Species}/sample.txt')
    parser.add_argument('--region-info', required=True,
                         help='Path to output_gene_intergenic_intron_pos_length.csv')
    parser.add_argument('--cds-lengths', required=True, help='Path to output_cds_lengths.csv')
    parser.add_argument('--group-info', required=True,
                         help='Path to group_info.csv (species, Phylogeny, Life_form, Status)')
    parser.add_argument('--genome-regions', required=True,
                         help='Path to species_genome_regions.csv (used for the total_length column '
                              '= single-IR plastome length, for the TRUE whole-plastome genome-wide D)')
    parser.add_argument('--output-dir', default='selection_test_results/tajima_d_igs',
                         help='Output directory')
    parser.add_argument('--min-S', type=int, default=3, help='Minimum segregating sites (default: 3)')
    parser.add_argument('--min-n', type=int, default=4, help='Minimum sample size per species (default: 4)')
    parser.add_argument('--skip-plots', action='store_true', help='Skip figure generation')
    args = parser.parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print("Loading variant data...")
    variant_df = load_variant_data(args.variant_data)
    print("  Loaded {} variant records".format(len(variant_df)))

    print("Loading sample counts...")
    sample_counts = load_sample_counts(args.sample_dir)
    print("  Found {} species".format(len(sample_counts)))

    print("Loading CDS lengths...")
    cds_lengths, cds_total_by_species = load_cds_lengths(args.cds_lengths)

    print("Loading IGS region info ({})...".format(args.region_info))
    igs_regions_by_species, igs_total_by_species = load_igs_region_info(args.region_info)
    print("  Loaded intergenic regions for {} species, total lengths e.g. {}".format(
        len(igs_regions_by_species),
        {k: igs_total_by_species[k] for k in list(igs_total_by_species)[:3]}))

    print("Loading group info (wild/cultivar status)...")
    status_map = load_group_info(args.group_info)

    print("Loading whole-plastome (single-IR) total lengths ({})...".format(args.genome_regions))
    genome_total_lengths = load_genome_total_lengths(args.genome_regions)

    # --- Extract intergenic SNP records ---
    igs_snp = variant_df[
        (variant_df['var_type'] == 'snp') &
        (variant_df['region_type'].str.lower() == 'intergenic')
    ].copy()
    print("Intergenic SNP records: {}".format(len(igs_snp)))

    print("Assigning intergenic SNPs to IGS regions (interval join)...")
    igs_snp = assign_igs_regions(igs_snp, igs_regions_by_species)
    n_unmapped_records = igs_snp['igs_region'].isna().sum()
    n_unmapped_positions = igs_snp.loc[igs_snp['igs_region'].isna(), ['species', 'position']].drop_duplicates().shape[0]
    print("  Unmapped records: {} ({} distinct positions) -- retained for genome-wide pooling, "
          "excluded from per-region table".format(n_unmapped_records, n_unmapped_positions))

    region_lengths_lookup = {}
    for sp, df in igs_regions_by_species.items():
        for _, row in df.iterrows():
            region_lengths_lookup[(sp, row['region_name'])] = row['length']

    # --- 1. Per-region (exploratory) ---
    print("Calculating per-(species, IGS region) Tajima's D...")
    igs_region_df, igs_region_skipped = calculate_tajima_d_per_region(
        igs_snp, sample_counts, region_lengths_lookup, min_S=args.min_S, min_n=args.min_n)
    print("  Computed D for {} (species, IGS region) pairs; skipped {}".format(
        len(igs_region_df), len(igs_region_skipped)))

    # --- 2. Genome-wide IGS (primary) ---
    print("Calculating genome-wide IGS Tajima's D (one value per species)...")
    igs_gw_df, igs_gw_skipped = calculate_tajima_d_genomewide(
        igs_snp, sample_counts, igs_total_by_species, min_S=args.min_S, min_n=args.min_n,
        region_label='IGS_genomewide')

    # --- 3. Genome-wide CDS (reference / symmetric contrast) ---
    print("Calculating genome-wide CDS Tajima's D (one value per species, for contrast)...")
    cds_snp = variant_df[
        (variant_df['var_type'] == 'snp') &
        (variant_df['region_type'].str.upper() == 'CDS') &
        (variant_df['gene'].notna()) & (variant_df['gene'] != '')
    ].copy()
    cds_gw_df, cds_gw_skipped = calculate_tajima_d_genomewide(
        cds_snp, sample_counts, cds_total_by_species, min_S=args.min_S, min_n=args.min_n,
        region_label='CDS_genomewide')

    # --- 4. TRUE whole-plastome genome-wide D (ALL region types, ALL SNPs) ---
    print("Calculating TRUE whole-plastome genome-wide Tajima's D "
          "(all region types pooled, single-IR reference length)...")
    plastome_snp = variant_df[variant_df['var_type'] == 'snp'].copy()
    plastome_gw_df, plastome_gw_skipped = calculate_tajima_d_genomewide(
        plastome_snp, sample_counts, genome_total_lengths, min_S=args.min_S, min_n=args.min_n,
        region_label='PLASTOME_genomewide')

    # --- Wild vs. cultivated tests ---
    print("Running wild-vs-cultivated comparisons...")
    wvc_df = run_wild_vs_cultivar_tests(plastome_gw_df, igs_gw_df, cds_gw_df, igs_region_df, status_map)

    # --- Congeneric paired comparison ---
    pairs_df = build_congeneric_pairs_table(plastome_gw_df, igs_gw_df, cds_gw_df)

    # --- SFS by status ---
    print("Computing pooled folded SFS by domestication status...")
    sfs_df = compute_folded_sfs_by_status(igs_snp, sample_counts, status_map)

    # --- Write outputs ---
    igs_region_df.to_csv(out_dir / 'all_species_tajima_d_igs_per_region.tsv', sep='\t', index=False)
    igs_region_skipped.to_csv(out_dir / 'skipped_igs_regions.tsv', sep='\t', index=False)
    igs_gw_df.to_csv(out_dir / 'all_species_tajima_d_igs_genomewide.tsv', sep='\t', index=False)
    cds_gw_df.to_csv(out_dir / 'all_species_tajima_d_cds_genomewide.tsv', sep='\t', index=False)
    plastome_gw_df.to_csv(out_dir / 'all_species_tajima_d_plastome_genomewide.tsv', sep='\t', index=False)
    pd.concat([igs_gw_skipped, cds_gw_skipped, plastome_gw_skipped], ignore_index=True).to_csv(
        out_dir / 'skipped_genomewide.tsv', sep='\t', index=False)
    wvc_df.to_csv(out_dir / 'wild_vs_cultivar_comparison.tsv', sep='\t', index=False)
    pairs_df.to_csv(out_dir / 'congeneric_pairs_comparison.tsv', sep='\t', index=False)
    sfs_df.to_csv(out_dir / 'sfs_igs_by_status.tsv', sep='\t', index=False)

    qc_rows = [{
        'n_intergenic_snp_records': len(igs_snp),
        'n_unmapped_records': int(n_unmapped_records),
        'n_unmapped_distinct_positions': int(n_unmapped_positions),
        'pct_unmapped_records': round(100.0 * n_unmapped_records / len(igs_snp), 3) if len(igs_snp) else np.nan,
    }]
    pd.DataFrame(qc_rows).to_csv(out_dir / 'igs_region_assignment_qc.tsv', sep='\t', index=False)

    if not args.skip_plots:
        print("Generating figures...")
        make_plots(plastome_gw_df, igs_gw_df, cds_gw_df, sfs_df, pairs_df, status_map, out_dir)

    # --- Console summary ---
    print("\n=== TRUE whole-plastome genome-wide Tajima's D by species ===")
    plastome_gw_with_status = add_status(plastome_gw_df, status_map).sort_values(['status', 'species'])
    for _, r in plastome_gw_with_status.iterrows():
        print("  {:22s} status={:9s} n={:5d} S={:5d} D={:8.4f} sig={}".format(
            r['species'], str(r['status']), int(r['n']), int(r['S']), r['tajima_d'], r['sig_level']))

    print("\n=== Genome-wide IGS-only Tajima's D by species ===")
    igs_gw_with_status = add_status(igs_gw_df, status_map).sort_values(['status', 'species'])
    for _, r in igs_gw_with_status.iterrows():
        print("  {:22s} status={:9s} n={:5d} S={:4d} D={:8.4f} sig={}".format(
            r['species'], str(r['status']), int(r['n']), int(r['S']), r['tajima_d'], r['sig_level']))

    print("\n=== Genome-wide CDS-only Tajima's D by species ===")
    cds_gw_with_status = add_status(cds_gw_df, status_map).sort_values(['status', 'species'])
    for _, r in cds_gw_with_status.iterrows():
        print("  {:22s} status={:9s} n={:5d} S={:4d} D={:8.4f} sig={}".format(
            r['species'], str(r['status']), int(r['n']), int(r['S']), r['tajima_d'], r['sig_level']))

    print("\n=== Wild vs. Cultivar tests ===")
    for _, r in wvc_df.iterrows():
        print("  [{}] {} (n={}, mean={:.4f}) vs {} (n={}, mean={:.4f}) -> U={}, p={}".format(
            r['comparison_level'], r['group_a'], r['n_a'], r['mean_a'],
            r['group_b'], r['n_b'], r['mean_b'], r['mannwhitney_U'], r['p_value']))

    print("\n=== Congeneric pairs (Plastome / IGS / CDS genome-wide D) ===")
    for _, r in pairs_df.iterrows():
        flag = '' if r['usable_for_domestication_test'] else '  [excluded: same status]'
        print("  {}: plastome {}={} vs {}={} | igs {}={} vs {}={} | cds {}={} vs {}={}{}".format(
            r['genus'],
            r['species_a'], r['plastome_D_a'], r['species_b'], r['plastome_D_b'],
            r['species_a'], r['igs_D_a'], r['species_b'], r['igs_D_b'],
            r['species_a'], r['cds_D_a'], r['species_b'], r['cds_D_b'],
            flag))

    print("\nDone. Outputs written to: {}".format(out_dir))


if __name__ == '__main__':
    main()
