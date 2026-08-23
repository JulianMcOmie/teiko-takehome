"""Part 2 - relative frequency of each cell population within each sample.

For every sample the five population counts are summed to a total, and each
population is expressed as a percentage of that total. The arithmetic lives in the
`sample_frequencies` view (see schema.sql) so that this export, the statistics in
Part 3 and the dashboard all read the same definition.
"""

from src.config import OUTPUT_DIR
from src.db import query

SUMMARY_SQL = """
SELECT sample, total_count, population, count, percentage
FROM sample_frequencies
ORDER BY sample, population
"""

OUTPUT_PATH = OUTPUT_DIR / "summary_table.csv"


def summary_table():
    """Return the Part 2 table: one row per sample-population pair."""
    return query(SUMMARY_SQL)


def main():
    table = summary_table()

    # Percentages are rounded only on the way out. Part 3 reads the view directly, so
    # the statistics are computed on unrounded values.
    export = table.assign(percentage=table["percentage"].round(2))

    OUTPUT_DIR.mkdir(exist_ok=True)
    export.to_csv(OUTPUT_PATH, index=False)

    print(f"Part 2 - summary table: {len(export):,} rows "
          f"({export['sample'].nunique():,} samples x {export['population'].nunique()} populations)")
    print(export.head(5).to_string(index=False))
    print(f"Written to {OUTPUT_PATH.relative_to(OUTPUT_DIR.parent)}")


if __name__ == "__main__":
    main()
