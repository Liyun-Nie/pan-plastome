"""Load gene functional categories from cpopvar gene_function_map.csv."""

from pathlib import Path

import pandas as pd

GENE_CATEGORY_ORDER = [
    'Photosynthesis Related',
    'Self-Replication Related',
    'Other Functional Genes',
    'Unknown Function Genes',
]

GENE_CATEGORY_SHORT_LABELS = {
    'Photosynthesis Related': 'Photosynthesis',
    'Self-Replication Related': 'Self-replication',
    'Other Functional Genes': 'Other functional',
    'Unknown Function Genes': 'Unknown function',
}

GENE_CATEGORY_COLORS = {
    'Photosynthesis Related': '#4CAF50',
    'Self-Replication Related': '#2196F3',
    'Other Functional Genes': '#FFC107',
    'Unknown Function Genes': '#9C27B0',
}

DEFAULT_GENE_FUNCTION_MAP = (
    Path(__file__).resolve().parent / 'gene_function_map.csv'
)


def load_gene_function_map(path=None):
    """Return {gene_id: gene_category} from cpopvar reference file."""
    map_path = Path(path) if path else DEFAULT_GENE_FUNCTION_MAP
    if not map_path.exists():
        raise FileNotFoundError(f'Gene function map not found: {map_path}')

    df = pd.read_csv(map_path, dtype=str)
    df['gene_id'] = df['gene_id'].str.strip()
    df['gene_category'] = df['gene_category'].str.strip()
    return dict(zip(df['gene_id'], df['gene_category']))


def get_functional_category(gene_name, category_map):
    """Look up gene_category; unmapped genes default to Unknown Function Genes."""
    return category_map.get(gene_name, 'Unknown Function Genes')
