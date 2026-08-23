#!/usr/bin/env python3
"""Part 1 - initialize the SQLite database and load cell-count.csv into it.

Usage, from the repository root:

    python load_data.py

The script takes no arguments and writes cell_counts.db next to itself. It is
idempotent: the schema is dropped and rebuilt on every run, so running it twice
yields the same database rather than duplicate-key errors.
"""

from contextlib import closing

import pandas as pd

from src.config import CSV_PATH, DB_PATH, POPULATIONS, SCHEMA_PATH
from src.db import connect

COUNT_COLUMNS = list(POPULATIONS)


def read_source_csv(csv_path=CSV_PATH):
    """Read the CSV and normalize the two values that need it before loading."""
    if not csv_path.exists():
        raise FileNotFoundError(f"Expected the input data at {csv_path}")

    df = pd.read_csv(csv_path)

    missing = set(COUNT_COLUMNS) - set(df.columns)
    if missing:
        raise ValueError(f"Input is missing expected count columns: {sorted(missing)}")

    # Healthy, untreated subjects have no response assessment. pandas reads the blank
    # field as NaN; make that an explicit None so it lands in SQLite as NULL.
    df["response"] = df["response"].where(df["response"].notna(), None)
    return df


def build_schema(conn, schema_path=SCHEMA_PATH):
    """Create tables, indexes and views from schema.sql."""
    conn.executescript(schema_path.read_text())


def subject_frame(df):
    """One row per subject, with the clinical attributes that describe the person.

    The CSV repeats these across a subject's three samples; drop_duplicates collapses
    them, and the verification step confirms the repetition really was consistent.
    """
    columns = ["subject", "project", "condition", "age", "sex", "treatment", "response"]
    subjects = df[columns].drop_duplicates(subset="subject")
    return subjects.rename(columns={"subject": "subject_id", "project": "project_id"})


def long_counts(df):
    """Melt the five wide count columns into one row per (sample, population)."""
    counts = df.melt(
        id_vars="sample",
        value_vars=COUNT_COLUMNS,
        var_name="population",
        value_name="count",
    )
    return counts.rename(columns={"sample": "sample_id"})


def load(conn, df):
    """Insert every table in foreign-key order."""
    projects = pd.DataFrame({"project_id": sorted(df["project"].unique())})
    populations = pd.DataFrame(
        {"population": list(POPULATIONS), "display_name": list(POPULATIONS.values())}
    )
    subjects = subject_frame(df)
    samples = df[
        ["sample", "subject", "sample_type", "time_from_treatment_start"]
    ].rename(columns={"sample": "sample_id", "subject": "subject_id"})

    for name, frame in [
        ("projects", projects),
        ("populations", populations),
        ("subjects", subjects),
        ("samples", samples),
        ("cell_counts", long_counts(df)),
    ]:
        frame.to_sql(name, conn, if_exists="append", index=False)
        print(f"  loaded {len(frame):>6,} rows into {name}")


def verify(conn, df):
    """Check the loaded database against the source file before declaring success.

    A silent partial load is the expensive failure here, because every downstream
    number would still look plausible.
    """
    scalar = lambda sql: conn.execute(sql).fetchone()[0]

    expected_subjects = df["subject"].nunique()
    checks = {
        "samples": (scalar("SELECT COUNT(*) FROM samples"), len(df)),
        "subjects": (scalar("SELECT COUNT(*) FROM subjects"), expected_subjects),
        "cell_counts": (
            scalar("SELECT COUNT(*) FROM cell_counts"),
            len(df) * len(COUNT_COLUMNS),
        ),
    }
    for table, (found, expected) in checks.items():
        if found != expected:
            raise RuntimeError(f"{table}: loaded {found} rows, expected {expected}")

    # The clinical attributes must have been genuinely constant per subject; if they
    # were not, drop_duplicates silently kept an arbitrary one.
    varying = df.groupby("subject")[
        ["project", "condition", "age", "sex", "treatment", "response"]
    ].nunique(dropna=False)
    if (varying > 1).any().any():
        offenders = varying[(varying > 1).any(axis=1)].index.tolist()
        raise RuntimeError(
            f"Subject-level attributes disagree between samples for: {offenders[:5]}"
        )

    # Total cells per sample must survive the wide-to-long reshape untouched.
    loaded_total = scalar("SELECT SUM(count) FROM cell_counts")
    source_total = int(df[COUNT_COLUMNS].to_numpy().sum())
    if loaded_total != source_total:
        raise RuntimeError(
            f"Total cell count mismatch: loaded {loaded_total}, source {source_total}"
        )

    print(
        f"  verified {checks['samples'][0]:,} samples, "
        f"{checks['subjects'][0]:,} subjects, "
        f"{checks['cell_counts'][0]:,} measurements, "
        f"{loaded_total:,} cells total"
    )


def main():
    print(f"Reading {CSV_PATH.name}")
    df = read_source_csv()

    if DB_PATH.exists():
        DB_PATH.unlink()
    print(f"Building {DB_PATH.name}")

    with closing(connect()) as conn:
        build_schema(conn)
        load(conn, df)
        verify(conn, df)
        conn.commit()

    print(f"Done. Database written to {DB_PATH}")


if __name__ == "__main__":
    main()
