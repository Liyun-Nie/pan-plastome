#!/usr/bin/env python3
"""
Tajima's D & SFS Visualization for Plastome CDS Genes
======================================================
Generates publication-quality figures from Tajima's D results and raw variant data.

Outputs:
  - tajima_d_heatmap.pdf          (Panel A)
  - tajima_d_hotspot_heatmap.pdf  (Core shared hotspot genes)
  - tajima_d_category_boxplot.pdf (Panel B)
  - tajima_d_species_boxplot.pdf  (Panel C)
  - tajima_d_key_genes.pdf        (Panel D)
  - tajima_d_results_composite.pdf (A+B+C+D)
  - sfs_results.pdf               (Folded SFS facet grid)
  - PNG copies at 300 DPI

Style: aligned with cpopvar professional_light theme.
"""

import argparse
import os
import sys
import warnings
from pathlib import Path
from collections import defaultdict

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.lines import Line2D
import seaborn as sns
from scipy import stats

from gene_function_utils import (
    GENE_CATEGORY_COLORS,
    GENE_CATEGORY_ORDER,
    GENE_CATEGORY_SHORT_LABELS,
)

warnings.filterwarnings('ignore', category=UserWarning)

# ---------------------------------------------------------------------------
# cpopvar color constants (from YAML config & ggsci palettes)
# ---------------------------------------------------------------------------

D3_CATEGORY10 = [
    '#1F77B4', '#FF7F0E', '#2CA02C', '#D62728', '#9467BD',
    '#8C564B', '#E377C2', '#7F7F7F', '#BCBD22', '#17BECF',
]

HEATMAP_DIVERGENT = ['#2166AC', '#FFFFFF', '#D62728']

HIST_FILL = '#4DBBD5'
MISSING_FILL = '#F0F0F0'

ALPHA_FILL = 0.7
ALPHA_POINT = 0.4


# ---------------------------------------------------------------------------
# Style setup
# ---------------------------------------------------------------------------

def setup_cpopvar_style():
    """Configure matplotlib rcParams to match cpopvar professional_light theme."""
    font_family = 'Arial'
    try:
        from matplotlib.font_manager import fontManager
        available = {f.name for f in fontManager.ttflist}
        if 'Arial' not in available:
            font_family = 'DejaVu Sans'
    except Exception:
        font_family = 'DejaVu Sans'

    plt.rcParams.update({
        'font.family': font_family,
        'font.size': 10,
        'axes.titlesize': 12,
        'axes.titleweight': 'bold',
        'axes.labelsize': 10,
        'xtick.labelsize': 8,
        'ytick.labelsize': 8,
        'legend.fontsize': 8,
        'legend.title_fontsize': 10,
        'axes.linewidth': 0.5,
        'axes.edgecolor': 'black',
        'axes.facecolor': 'white',
        'figure.facecolor': 'white',
        'axes.grid': False,
        'axes.axisbelow': True,
        'xtick.major.size': 4,
        'ytick.major.size': 4,
        'xtick.major.width': 0.4,
        'ytick.major.width': 0.4,
        'xtick.direction': 'out',
        'ytick.direction': 'out',
        'figure.dpi': 150,
        'savefig.dpi': 300,
        'savefig.bbox': 'tight',
        'savefig.facecolor': 'white',
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
    })


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_tajima_results(path):
    df = pd.read_csv(path, sep='\t')
    return df


def load_species_order(path):
    df = pd.read_csv(path)
    return list(zip(df['species'], df['custom_label']))


def load_sample_counts(sample_dir):
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


def load_variant_data(path):
    df = pd.read_csv(path, dtype=str)
    col_map = {
        df.columns[0]: 'sample_id',
        'V2': 'species', 'V3': 'position', 'V4': 'var_type',
        'V5': 'ref_allele', 'V6': 'alt_allele', 'V7': 'evidence',
        'V8': 'region_type',
    }
    df = df.rename(columns=col_map)
    df['position'] = pd.to_numeric(df['position'], errors='coerce')
    return df


# ---------------------------------------------------------------------------
# Key gene selection
# ---------------------------------------------------------------------------

def load_core_shared_hotspot_genes(path, min_species=3, top_n=20):
    """Load core shared hotspot CDS (n_species >= min_species) from M03 summary."""
    df = pd.read_csv(path)
    df = df[df['n_species'] >= min_species].sort_values(
        ['n_species', 'avg_frequency'], ascending=[False, False]
    ).head(top_n)

    genes = df['gene'].tolist()
    n_shared = dict(zip(df['gene'], df['n_species']))
    membership = {}
    for _, row in df.iterrows():
        species = [s.strip() for s in str(row['species_list']).split(',')]
        membership[row['gene']] = set(species)
    return genes, n_shared, membership


def select_key_genes(df, min_species=3, top_n=4):
    """Select key genes aligned with top shared hotspot CDS from manuscript."""
    gene_stats = df.groupby('gene').agg(
        n_species=('species', 'nunique'),
        mean_S=('S', 'mean'),
    ).reset_index()

    preferred = ['rps15', 'ycf1', 'rps3', 'accD']
    selected = [g for g in preferred if g in gene_stats['gene'].values]

    if len(selected) < top_n:
        candidates = gene_stats[gene_stats['n_species'] >= min_species]
        remaining = candidates[~candidates['gene'].isin(selected)]
        remaining = remaining.nlargest(top_n - len(selected), 'mean_S')
        selected.extend(remaining['gene'].tolist())

    return selected[:top_n]


# ---------------------------------------------------------------------------
# Panel A: Species-Gene Heatmap
# ---------------------------------------------------------------------------

def plot_heatmap(df, species_order, output_path):
    sp_names = [s[0] for s in species_order]
    sp_labels = [s[1] for s in species_order]

    present_sp = df['species'].unique()
    sp_names_filtered = [s for s in sp_names if s in present_sp]
    sp_labels_filtered = [l for s, l in species_order if s in present_sp]

    gene_counts = df.groupby('gene')['species'].nunique()
    genes_in_3plus = gene_counts[gene_counts >= 3].index
    df_filtered = df[df['gene'].isin(genes_in_3plus)]

    gene_mean_d = df_filtered.groupby('gene')['tajima_d'].mean().sort_values()
    gene_order = gene_mean_d.index.tolist()

    pivot = df_filtered.pivot(index='species', columns='gene', values='tajima_d')
    pivot = pivot.reindex(index=sp_names_filtered, columns=gene_order)

    sig_pivot = df_filtered.pivot(index='species', columns='gene', values='sig_level')
    sig_pivot = sig_pivot.reindex(index=sp_names_filtered, columns=gene_order)

    n_genes = len(gene_order)
    fig_width = max(10, n_genes * 0.3 + 3)
    fig_height = max(5, len(sp_names_filtered) * 0.35 + 2)
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    cmap = LinearSegmentedColormap.from_list('divergent', HEATMAP_DIVERGENT, N=256)

    d_values = pivot.values.flatten()
    d_values = d_values[~np.isnan(d_values)]
    vmax = max(abs(d_values.min()), abs(d_values.max())) if len(d_values) > 0 else 2
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)

    mask = pivot.isna()

    ax.set_facecolor(MISSING_FILL)
    im = ax.imshow(
        np.ma.masked_where(mask.values, pivot.values),
        cmap=cmap, norm=norm, aspect='auto', interpolation='nearest'
    )

    for i in range(pivot.shape[0]):
        for j in range(pivot.shape[1]):
            sig_val = sig_pivot.iloc[i, j]
            if pd.notna(sig_val):
                try:
                    sig_num = float(sig_val)
                    if sig_num <= 0.05:
                        text_color = 'white' if abs(pivot.iloc[i, j]) > vmax * 0.6 else 'black'
                        ax.text(j, i, '*', ha='center', va='center',
                                fontsize=9, fontweight='bold', color=text_color)
                except (ValueError, TypeError):
                    pass

    ax.set_xticks(range(n_genes))
    ax.set_xticklabels(gene_order, rotation=90, fontsize=6, fontstyle='italic')
    ax.set_yticks(range(len(sp_names_filtered)))
    ax.set_yticklabels(sp_labels_filtered, fontsize=8, fontstyle='italic')

    cbar = fig.colorbar(im, ax=ax, shrink=0.6, pad=0.02)
    cbar.set_label("Tajima's D", fontsize=9)
    cbar.ax.tick_params(labelsize=7)

    ax.set_title("Tajima's D across species and genes (* p \u2264 0.05)", pad=12)

    fig.tight_layout()
    fig.savefig(output_path)
    _save_png(fig, output_path)
    plt.close(fig)
    print(f"  Saved heatmap: {output_path}")


def plot_hotspot_heatmap(df, species_order, hotspot_genes, hotspot_membership, output_path):
    """Heatmap for core shared hotspot CDS (fixed column order by n_shared)."""
    sp_names = [s[0] for s in species_order]
    sp_labels = [s[1] for s in species_order]

    present_sp = df['species'].unique()
    sp_names_filtered = [s for s in sp_names if s in present_sp]
    sp_labels_filtered = [l for s, l in species_order if s in present_sp]

    gene_order = hotspot_genes

    df_filtered = df[df['gene'].isin(gene_order)]

    pivot = df_filtered.pivot(index='species', columns='gene', values='tajima_d')
    pivot = pivot.reindex(index=sp_names_filtered, columns=gene_order)

    sig_pivot = df_filtered.pivot(index='species', columns='gene', values='sig_level')
    sig_pivot = sig_pivot.reindex(index=sp_names_filtered, columns=gene_order)

    n_genes = len(gene_order)
    fig_width = max(12, n_genes * 0.35 + 3)
    fig_height = max(5, len(sp_names_filtered) * 0.35 + 2)
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    cmap = LinearSegmentedColormap.from_list('divergent', HEATMAP_DIVERGENT, N=256)

    d_values = pivot.values.flatten()
    d_values = d_values[~np.isnan(d_values)]
    vmax = max(abs(d_values.min()), abs(d_values.max())) if len(d_values) > 0 else 2
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)

    mask = pivot.isna()

    ax.set_facecolor(MISSING_FILL)
    im = ax.imshow(
        np.ma.masked_where(mask.values, pivot.values),
        cmap=cmap, norm=norm, aspect='auto', interpolation='nearest'
    )

    for i, species in enumerate(sp_names_filtered):
        for j, gene in enumerate(gene_order):
            sig_val = sig_pivot.iloc[i, j]
            if pd.notna(sig_val):
                try:
                    sig_num = float(sig_val)
                    if sig_num <= 0.05:
                        text_color = 'white' if abs(pivot.iloc[i, j]) > vmax * 0.6 else 'black'
                        ax.text(j, i, '*', ha='center', va='center',
                                fontsize=9, fontweight='bold', color=text_color)
                except (ValueError, TypeError):
                    pass

            if gene in hotspot_membership and species in hotspot_membership[gene]:
                if pd.notna(pivot.iloc[i, j]):
                    ax.add_patch(plt.Rectangle(
                        (j - 0.5, i - 0.5), 1, 1, fill=False,
                        edgecolor='black', linewidth=1.2, zorder=5
                    ))

    ax.set_xticks(range(n_genes))
    ax.set_xticklabels(gene_order, rotation=90, fontsize=6, fontstyle='italic')
    ax.set_yticks(range(len(sp_names_filtered)))
    ax.set_yticklabels(sp_labels_filtered, fontsize=8, fontstyle='italic')

    cbar = fig.colorbar(im, ax=ax, shrink=0.6, pad=0.02)
    cbar.set_label("Tajima's D", fontsize=9)
    cbar.ax.tick_params(labelsize=7)

    ax.set_title(
        "Tajima's D for core shared hotspot CDS (n_shared \u2265 3; * p \u2264 0.05; box = hotspot)",
        pad=12
    )

    fig.tight_layout()
    fig.savefig(output_path)
    _save_png(fig, output_path)
    plt.close(fig)
    print(f"  Saved hotspot heatmap: {output_path}")


# ---------------------------------------------------------------------------
# Panel B: Functional Category Boxplot
# ---------------------------------------------------------------------------

def plot_category_boxplot(df, output_path):
    valid = df.dropna(subset=['tajima_d']).copy()
    cat_order = GENE_CATEGORY_ORDER
    cat_labels = [GENE_CATEGORY_SHORT_LABELS[c] for c in cat_order]
    valid['functional_category'] = pd.Categorical(
        valid['functional_category'], categories=cat_order, ordered=True
    )

    fig, ax = plt.subplots(figsize=(6.5, 5))

    palette = [GENE_CATEGORY_COLORS[c] for c in cat_order]

    bp = sns.boxplot(
        data=valid, x='functional_category', y='tajima_d',
        order=cat_order, palette=palette,
        width=0.6, linewidth=0.8, fliersize=3,
        boxprops=dict(alpha=ALPHA_FILL),
        ax=ax,
    )
    sns.stripplot(
        data=valid, x='functional_category', y='tajima_d',
        order=cat_order, color='black', alpha=ALPHA_POINT,
        size=3, jitter=0.2, ax=ax,
    )

    ax.axhline(y=0, color='black', linestyle='--', linewidth=0.8, alpha=0.6)
    ax.set_xticklabels(cat_labels, fontsize=9)
    ax.set_xlabel('')
    ax.set_ylabel("Tajima's D", fontsize=10)
    ax.set_title("Tajima's D by functional category", pad=12)

    groups = [valid[valid['functional_category'] == c]['tajima_d'].values for c in cat_order]
    groups = [g for g in groups if len(g) > 0]
    if len(groups) >= 2:
        try:
            h_stat, p_val = stats.kruskal(*groups)
            ax.text(0.95, 0.95, f'Kruskal-Wallis p = {p_val:.4f}',
                    transform=ax.transAxes, ha='right', va='top', fontsize=7,
                    bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                              edgecolor='grey', alpha=0.8))
        except Exception:
            pass

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='y', color='#F2F2F2', linewidth=0.25)

    fig.tight_layout()
    fig.savefig(output_path)
    _save_png(fig, output_path)
    plt.close(fig)
    print(f"  Saved category boxplot: {output_path}")


# ---------------------------------------------------------------------------
# Panel C: Per-Species Boxplot
# ---------------------------------------------------------------------------

def plot_species_boxplot(df, species_order, output_path):
    valid = df.dropna(subset=['tajima_d']).copy()
    sp_names = [s[0] for s in species_order]
    sp_labels = [s[1] for s in species_order]

    present_sp = valid['species'].unique()
    sp_names_f = [s for s in sp_names if s in present_sp]
    sp_labels_f = [l for s, l in species_order if s in present_sp]

    valid['species'] = pd.Categorical(valid['species'], categories=sp_names_f, ordered=True)
    valid = valid.sort_values('species')

    fig, ax = plt.subplots(figsize=(12, 5))

    sns.boxplot(
        data=valid, x='species', y='tajima_d',
        order=sp_names_f, color=HIST_FILL,
        width=0.6, linewidth=0.8, fliersize=0,
        boxprops=dict(alpha=ALPHA_FILL),
        ax=ax,
    )
    sns.stripplot(
        data=valid, x='species', y='tajima_d',
        order=sp_names_f, color='#333333', alpha=ALPHA_POINT,
        size=3, jitter=0.15, ax=ax,
    )

    ax.axhline(y=0, color='black', linestyle='--', linewidth=0.8, alpha=0.6)
    ax.set_xticklabels(sp_labels_f, rotation=45, ha='right', fontsize=8, fontstyle='italic')
    ax.set_xlabel('')
    ax.set_ylabel("Tajima's D", fontsize=10)
    ax.set_title("Per-species distribution of Tajima's D", pad=12)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='y', color='#F2F2F2', linewidth=0.25)

    fig.tight_layout()
    fig.savefig(output_path)
    _save_png(fig, output_path)
    plt.close(fig)
    print(f"  Saved species boxplot: {output_path}")


# ---------------------------------------------------------------------------
# Panel D: Key Gene Comparison Line Plot
# ---------------------------------------------------------------------------

def plot_key_gene_lines(df, species_order, genes, output_path):
    valid = df.dropna(subset=['tajima_d']).copy()

    sp_names = [s[0] for s in species_order]
    sp_labels = [s[1] for s in species_order]
    present_sp = valid['species'].unique()
    sp_names_f = [s for s in sp_names if s in present_sp]
    sp_labels_f = [l for s, l in species_order if s in present_sp]

    sp_to_idx = {s: i for i, s in enumerate(sp_names_f)}

    fig, ax = plt.subplots(figsize=(12, 5))

    for gi, gene in enumerate(genes):
        gene_data = valid[valid['gene'] == gene].copy()
        if gene_data.empty:
            continue

        gene_data = gene_data[gene_data['species'].isin(sp_names_f)]
        gene_data['x_idx'] = gene_data['species'].map(sp_to_idx)
        gene_data = gene_data.sort_values('x_idx')

        color = D3_CATEGORY10[gi % len(D3_CATEGORY10)]

        ax.plot(gene_data['x_idx'], gene_data['tajima_d'],
                '-', color=color, linewidth=1.5, label=gene, alpha=0.8, zorder=2)

        for _, row in gene_data.iterrows():
            sig = row['sig_level']
            is_sig = False
            try:
                is_sig = float(sig) <= 0.05
            except (ValueError, TypeError):
                pass

            marker = 'o' if is_sig else 'o'
            fc = color if is_sig else 'white'
            ec = color
            ax.plot(row['x_idx'], row['tajima_d'], marker,
                    markersize=6, markerfacecolor=fc, markeredgecolor=ec,
                    markeredgewidth=1.2, zorder=3)

    ax.axhline(y=0, color='black', linestyle='--', linewidth=0.8, alpha=0.6)
    ax.set_xticks(range(len(sp_names_f)))
    ax.set_xticklabels(sp_labels_f, rotation=45, ha='right', fontsize=8, fontstyle='italic')
    ax.set_xlabel('')
    ax.set_ylabel("Tajima's D", fontsize=10)
    ax.set_title("Tajima's D for key genes across species", pad=12)

    sig_marker = Line2D([0], [0], marker='o', color='grey', markerfacecolor='grey',
                        markersize=6, linewidth=0, label='p \u2264 0.05')
    ns_marker = Line2D([0], [0], marker='o', color='grey', markerfacecolor='white',
                       markeredgecolor='grey', markersize=6, linewidth=0,
                       markeredgewidth=1.2, label='n.s.')
    gene_handles = [Line2D([0], [0], color=D3_CATEGORY10[i], linewidth=1.5, label=g)
                    for i, g in enumerate(genes)]

    ax.legend(handles=gene_handles + [sig_marker, ns_marker],
              loc='lower center', bbox_to_anchor=(0.5, -0.35),
              ncol=len(genes) + 2, frameon=True, edgecolor='grey',
              fontsize=7, handlelength=1.5)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='y', color='#F2F2F2', linewidth=0.25)

    fig.tight_layout()
    fig.savefig(output_path)
    _save_png(fig, output_path)
    plt.close(fig)
    print(f"  Saved key gene lines: {output_path}")


# ---------------------------------------------------------------------------
# Composite figure (A+B+C+D)
# ---------------------------------------------------------------------------

def create_composite_figure(df, species_order, genes, output_path):
    sp_names = [s[0] for s in species_order]
    sp_labels = [s[1] for s in species_order]
    present_sp = df['species'].unique()
    sp_names_f = [s for s in sp_names if s in present_sp]
    sp_labels_f = [l for s, l in species_order if s in present_sp]

    fig = plt.figure(figsize=(16, 18))
    gs = gridspec.GridSpec(3, 2, height_ratios=[1, 0.8, 0.8],
                           hspace=0.35, wspace=0.3)

    # --- Panel A: Heatmap (top, full width) ---
    ax_a = fig.add_subplot(gs[0, :])
    _draw_heatmap_on_ax(ax_a, fig, df, sp_names_f, sp_labels_f)
    ax_a.set_title("(A) Tajima's D across species and genes (* p \u2264 0.05)",
                   pad=12, fontsize=11, fontweight='bold')

    # --- Panel B: Category boxplot (middle-left) ---
    ax_b = fig.add_subplot(gs[1, 0])
    _draw_category_boxplot_on_ax(ax_b, df)
    ax_b.set_title("(B) Tajima's D by functional category",
                   pad=10, fontsize=11, fontweight='bold')

    # --- Panel C: Species boxplot (middle-right) ---
    ax_c = fig.add_subplot(gs[1, 1])
    _draw_species_boxplot_on_ax(ax_c, df, sp_names_f, sp_labels_f)
    ax_c.set_title("(C) Per-species distribution",
                   pad=10, fontsize=11, fontweight='bold')

    # --- Panel D: Key gene lines (bottom, full width) ---
    ax_d = fig.add_subplot(gs[2, :])
    _draw_key_gene_lines_on_ax(ax_d, df, sp_names_f, sp_labels_f, genes)
    ax_d.set_title("(D) Key genes across species",
                   pad=10, fontsize=11, fontweight='bold')

    fig.savefig(output_path)
    _save_png(fig, output_path)
    plt.close(fig)
    print(f"  Saved composite figure: {output_path}")


def _draw_heatmap_on_ax(ax, fig, df, sp_names_f, sp_labels_f):
    gene_counts = df.groupby('gene')['species'].nunique()
    genes_3plus = gene_counts[gene_counts >= 3].index
    df_f = df[df['gene'].isin(genes_3plus)]
    gene_mean = df_f.groupby('gene')['tajima_d'].mean().sort_values()
    gene_order = gene_mean.index.tolist()

    pivot = df_f.pivot(index='species', columns='gene', values='tajima_d')
    pivot = pivot.reindex(index=sp_names_f, columns=gene_order)

    sig_pivot = df_f.pivot(index='species', columns='gene', values='sig_level')
    sig_pivot = sig_pivot.reindex(index=sp_names_f, columns=gene_order)

    cmap = LinearSegmentedColormap.from_list('div', HEATMAP_DIVERGENT, N=256)
    d_vals = pivot.values.flatten()
    d_vals = d_vals[~np.isnan(d_vals)]
    vmax = max(abs(d_vals.min()), abs(d_vals.max())) if len(d_vals) > 0 else 2
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)

    ax.set_facecolor(MISSING_FILL)
    im = ax.imshow(
        np.ma.masked_where(pivot.isna().values, pivot.values),
        cmap=cmap, norm=norm, aspect='auto', interpolation='nearest'
    )

    for i in range(pivot.shape[0]):
        for j in range(pivot.shape[1]):
            sig_val = sig_pivot.iloc[i, j]
            if pd.notna(sig_val):
                try:
                    if float(sig_val) <= 0.05:
                        tc = 'white' if abs(pivot.iloc[i, j]) > vmax * 0.6 else 'black'
                        ax.text(j, i, '*', ha='center', va='center',
                                fontsize=7, fontweight='bold', color=tc)
                except (ValueError, TypeError):
                    pass

    ax.set_xticks(range(len(gene_order)))
    ax.set_xticklabels(gene_order, rotation=90, fontsize=5, fontstyle='italic')
    ax.set_yticks(range(len(sp_names_f)))
    ax.set_yticklabels(sp_labels_f, fontsize=7, fontstyle='italic')

    cbar = fig.colorbar(im, ax=ax, shrink=0.5, pad=0.02)
    cbar.set_label("Tajima's D", fontsize=8)
    cbar.ax.tick_params(labelsize=6)


def _draw_category_boxplot_on_ax(ax, df):
    valid = df.dropna(subset=['tajima_d']).copy()
    cat_order = GENE_CATEGORY_ORDER
    cat_labels = [GENE_CATEGORY_SHORT_LABELS[c][:10] for c in cat_order]
    palette = [GENE_CATEGORY_COLORS[c] for c in cat_order]

    sns.boxplot(data=valid, x='functional_category', y='tajima_d',
                order=cat_order, palette=palette, width=0.6, linewidth=0.8,
                fliersize=3, boxprops=dict(alpha=ALPHA_FILL), ax=ax)
    sns.stripplot(data=valid, x='functional_category', y='tajima_d',
                  order=cat_order, color='black', alpha=ALPHA_POINT,
                  size=2.5, jitter=0.2, ax=ax)

    ax.axhline(y=0, color='black', linestyle='--', linewidth=0.8, alpha=0.6)
    ax.set_xticklabels(cat_labels, fontsize=8)
    ax.set_xlabel('')
    ax.set_ylabel("Tajima's D", fontsize=9)

    groups = [valid[valid['functional_category'] == c]['tajima_d'].values for c in cat_order]
    groups = [g for g in groups if len(g) > 0]
    if len(groups) >= 2:
        try:
            _, p_val = stats.kruskal(*groups)
            ax.text(0.95, 0.95, f'KW p={p_val:.3f}', transform=ax.transAxes,
                    ha='right', va='top', fontsize=6,
                    bbox=dict(boxstyle='round,pad=0.2', fc='white', ec='grey', alpha=0.8))
        except Exception:
            pass

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='y', color='#F2F2F2', linewidth=0.25)


def _draw_species_boxplot_on_ax(ax, df, sp_names_f, sp_labels_f):
    valid = df.dropna(subset=['tajima_d']).copy()
    valid = valid[valid['species'].isin(sp_names_f)]

    sns.boxplot(data=valid, x='species', y='tajima_d',
                order=sp_names_f, color=HIST_FILL, width=0.6, linewidth=0.8,
                fliersize=0, boxprops=dict(alpha=ALPHA_FILL), ax=ax)
    sns.stripplot(data=valid, x='species', y='tajima_d',
                  order=sp_names_f, color='#333333', alpha=ALPHA_POINT,
                  size=2.5, jitter=0.15, ax=ax)

    ax.axhline(y=0, color='black', linestyle='--', linewidth=0.8, alpha=0.6)
    ax.set_xticklabels(sp_labels_f, rotation=45, ha='right', fontsize=6, fontstyle='italic')
    ax.set_xlabel('')
    ax.set_ylabel("Tajima's D", fontsize=9)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='y', color='#F2F2F2', linewidth=0.25)


def _draw_key_gene_lines_on_ax(ax, df, sp_names_f, sp_labels_f, genes):
    valid = df.dropna(subset=['tajima_d']).copy()
    valid = valid[valid['species'].isin(sp_names_f)]
    sp_to_idx = {s: i for i, s in enumerate(sp_names_f)}

    for gi, gene in enumerate(genes):
        gd = valid[valid['gene'] == gene].copy()
        if gd.empty:
            continue
        gd['x_idx'] = gd['species'].map(sp_to_idx)
        gd = gd.dropna(subset=['x_idx']).sort_values('x_idx')
        color = D3_CATEGORY10[gi % len(D3_CATEGORY10)]

        ax.plot(gd['x_idx'], gd['tajima_d'], '-', color=color,
                linewidth=1.2, alpha=0.8, zorder=2)

        for _, row in gd.iterrows():
            is_sig = False
            try:
                is_sig = float(row['sig_level']) <= 0.05
            except (ValueError, TypeError):
                pass
            fc = color if is_sig else 'white'
            ax.plot(row['x_idx'], row['tajima_d'], 'o', markersize=5,
                    markerfacecolor=fc, markeredgecolor=color,
                    markeredgewidth=1.0, zorder=3)

    ax.axhline(y=0, color='black', linestyle='--', linewidth=0.8, alpha=0.6)
    ax.set_xticks(range(len(sp_names_f)))
    ax.set_xticklabels(sp_labels_f, rotation=45, ha='right', fontsize=7, fontstyle='italic')
    ax.set_xlabel('')
    ax.set_ylabel("Tajima's D", fontsize=9)

    gene_handles = [Line2D([0], [0], color=D3_CATEGORY10[i], linewidth=1.2, label=g)
                    for i, g in enumerate(genes)]
    sig_h = Line2D([0], [0], marker='o', color='grey', markerfacecolor='grey',
                   markersize=5, linewidth=0, label='p\u22640.05')
    ns_h = Line2D([0], [0], marker='o', color='grey', markerfacecolor='white',
                  markeredgecolor='grey', markersize=5, linewidth=0,
                  markeredgewidth=1.0, label='n.s.')
    ax.legend(handles=gene_handles + [sig_h, ns_h],
              loc='lower center', bbox_to_anchor=(0.5, -0.30),
              ncol=len(genes) + 2, frameon=True, edgecolor='grey',
              fontsize=6, handlelength=1.5)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='y', color='#F2F2F2', linewidth=0.25)


# ---------------------------------------------------------------------------
# SFS computation
# ---------------------------------------------------------------------------

def compute_folded_sfs(variant_df, sample_counts, n_bins=10):
    """Compute folded SFS for each species from CDS SNPs."""
    cds_snps = variant_df[
        (variant_df['var_type'] == 'snp') &
        (variant_df['region_type'].str.upper() == 'CDS')
    ].copy()

    results = []
    for species, grp in cds_snps.groupby('species'):
        n = sample_counts.get(species)
        if n is None or n < 4:
            continue

        site_j = grp.groupby('position')['sample_id'].nunique()
        total_S = len(site_j)

        freqs = site_j.values / n
        folded = np.minimum(freqs, 1 - freqs)
        folded = folded[folded > 0]

        bin_edges = np.linspace(0, 0.5, n_bins + 1)
        counts, _ = np.histogram(folded, bins=bin_edges)

        for k in range(n_bins):
            results.append({
                'species': species,
                'freq_bin_left': round(bin_edges[k], 4),
                'freq_bin_right': round(bin_edges[k + 1], 4),
                'freq_bin_label': f'{bin_edges[k]:.2f}-{bin_edges[k+1]:.2f}',
                'count': int(counts[k]),
                'n': n,
                'total_S': total_S,
            })

    return pd.DataFrame(results)


# ---------------------------------------------------------------------------
# SFS facet grid
# ---------------------------------------------------------------------------

def plot_sfs_grid(sfs_df, species_order, output_path):
    sp_names = [s[0] for s in species_order]
    sp_labels_map = {s: l for s, l in species_order}
    present_sp = sfs_df['species'].unique()
    sp_ordered = [s for s in sp_names if s in present_sp]

    n_sp = len(sp_ordered)
    ncols = 5
    nrows = (n_sp + ncols - 1) // ncols

    fig, axes = plt.subplots(nrows, ncols, figsize=(14, nrows * 2.8),
                             sharex=True, sharey=False)
    axes = axes.flatten()

    bin_labels = sfs_df['freq_bin_label'].unique()

    for idx, sp in enumerate(sp_ordered):
        ax = axes[idx]
        sp_data = sfs_df[sfs_df['species'] == sp].sort_values('freq_bin_left')

        ax.bar(range(len(sp_data)), sp_data['count'].values,
               color=HIST_FILL, alpha=ALPHA_FILL, edgecolor='#3A9AB5', linewidth=0.5)

        label = sp_labels_map.get(sp, sp)
        n_val = sp_data['n'].iloc[0] if len(sp_data) > 0 else '?'
        s_val = sp_data['total_S'].iloc[0] if len(sp_data) > 0 else '?'

        ax.set_title(f'{label}  (n={n_val}, S={s_val})',
                     fontsize=8, fontstyle='italic', pad=4)

        if idx >= (nrows - 1) * ncols:
            x_labels = sp_data['freq_bin_label'].values
            ax.set_xticks(range(len(x_labels)))
            ax.set_xticklabels(x_labels, rotation=90, fontsize=5)
        else:
            ax.set_xticks(range(len(sp_data)))
            ax.set_xticklabels([])

        ax.tick_params(axis='y', labelsize=6)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(axis='y', color='#F2F2F2', linewidth=0.25)

    for idx in range(n_sp, len(axes)):
        axes[idx].set_visible(False)

    fig.supxlabel('Minor allele frequency', fontsize=10, y=0.02)
    fig.supylabel('Number of segregating sites', fontsize=10, x=0.02)
    fig.suptitle('Folded Site Frequency Spectrum (CDS SNPs)', fontsize=12,
                 fontweight='bold', y=0.98)

    fig.tight_layout(rect=[0.03, 0.04, 1, 0.96])
    fig.savefig(output_path)
    _save_png(fig, output_path)
    plt.close(fig)
    print(f"  Saved SFS grid: {output_path}")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _save_png(fig, pdf_path):
    png_path = str(pdf_path).replace('.pdf', '.png')
    fig.savefig(png_path, dpi=300, facecolor='white')


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Tajima's D & SFS visualization for plastome CDS genes")
    parser.add_argument('--tajima-results', required=True,
                        help='Path to all_species_tajima_d_summary.tsv')
    parser.add_argument('--variant-data', required=True,
                        help='Path to all_combined_data.csv')
    parser.add_argument('--sample-dir', required=True,
                        help='Directory containing {Species}/sample.txt')
    parser.add_argument('--species-order', required=True,
                        help='Path to species_label_order.csv')
    parser.add_argument('--output-dir',
                        default='selection_test_results/tajima_d/figures',
                        help='Output directory for figures')
    parser.add_argument('--hotspot-genes', default=None,
                        help='Path to M03_core_shared_genes_summary.csv for hotspot heatmap')
    args = parser.parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    setup_cpopvar_style()

    print("=" * 60)
    print("Tajima's D & SFS Visualization Pipeline")
    print("=" * 60)

    # 1. Load data
    print("\n[1/6] Loading data...")
    tajima_df = load_tajima_results(args.tajima_results)
    print(f"  Tajima's D results: {len(tajima_df)} records")

    species_order = load_species_order(args.species_order)
    print(f"  Species order: {len(species_order)} species")

    sample_counts = load_sample_counts(args.sample_dir)
    print(f"  Sample counts: {len(sample_counts)} species")

    variant_df = load_variant_data(args.variant_data)
    print(f"  Variant data: {len(variant_df)} records")

    # 2. Select key genes
    print("\n[2/6] Selecting key genes...")
    key_genes = select_key_genes(tajima_df)
    print(f"  Key genes: {key_genes}")

    # 3. Individual panels
    print("\n[3/6] Generating individual Tajima's D panels...")
    plot_heatmap(tajima_df, species_order,
                 out_dir / 'tajima_d_heatmap.pdf')

    if args.hotspot_genes:
        print("\n[3b/6] Generating hotspot-aligned heatmap...")
        hotspot_genes, _, hotspot_membership = load_core_shared_hotspot_genes(
            args.hotspot_genes
        )
        print(f"  Core shared hotspot genes: {len(hotspot_genes)}")
        plot_hotspot_heatmap(
            tajima_df, species_order, hotspot_genes, hotspot_membership,
            out_dir / 'tajima_d_hotspot_heatmap.pdf'
        )

    plot_category_boxplot(tajima_df,
                          out_dir / 'tajima_d_category_boxplot.pdf')
    plot_species_boxplot(tajima_df, species_order,
                         out_dir / 'tajima_d_species_boxplot.pdf')
    plot_key_gene_lines(tajima_df, species_order, key_genes,
                        out_dir / 'tajima_d_key_genes.pdf')

    # 4. Composite
    print("\n[4/6] Generating composite figure...")
    create_composite_figure(tajima_df, species_order, key_genes,
                            out_dir / 'tajima_d_results_composite.pdf')

    # 5. SFS
    print("\n[5/6] Computing and plotting SFS...")
    sfs_df = compute_folded_sfs(variant_df, sample_counts)
    sfs_out = out_dir / 'sfs_data.tsv'
    sfs_df.to_csv(sfs_out, sep='\t', index=False)
    print(f"  SFS data saved: {sfs_out} ({len(sfs_df)} rows)")

    plot_sfs_grid(sfs_df, species_order,
                  out_dir / 'sfs_results.pdf')

    # 6. Summary
    print("\n[6/6] Summary")
    print(f"  Output directory: {out_dir}")
    outputs = list(out_dir.glob('*'))
    for f in sorted(outputs):
        size_kb = f.stat().st_size / 1024
        print(f"    {f.name} ({size_kb:.1f} KB)")

    print("\nDone.")


if __name__ == '__main__':
    main()
