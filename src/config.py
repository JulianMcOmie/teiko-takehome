"""Shared configuration: filesystem paths, the cytometry panel, and analysis constants.

Anything more than one pipeline stage has to agree on lives here, so the loader,
the analyses and the dashboard cannot drift apart.
"""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CSV_PATH = ROOT / "cell-count.csv"
DB_PATH = ROOT / "cell_counts.db"
SCHEMA_PATH = ROOT / "src" / "schema.sql"
OUTPUT_DIR = ROOT / "outputs"

# The five measured populations, in the order the CSV presents them. Keys are both
# the CSV column names and the population ids used in the database.
POPULATIONS = {
    "b_cell": "B cell",
    "cd8_t_cell": "CD8+ T cell",
    "cd4_t_cell": "CD4+ T cell",
    "nk_cell": "NK cell",
    "monocyte": "Monocyte",
}

# Significance threshold applied to Benjamini-Hochberg adjusted p-values in Part 3.
ALPHA = 0.05
