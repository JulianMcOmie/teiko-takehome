"""SQLite access helpers shared by the loader, the analyses and the dashboard."""

import sqlite3
from contextlib import closing

import pandas as pd

from src.config import DB_PATH


class DatabaseMissingError(RuntimeError):
    """Raised when an analysis runs before the database has been built."""


def connect(db_path=None):
    """Open a connection with foreign key enforcement turned on.

    SQLite disables foreign keys per connection by default, so the PRAGMA has to be
    reissued every time rather than set once at build time.
    """
    conn = sqlite3.connect(db_path or DB_PATH)
    conn.execute("PRAGMA foreign_keys = ON")
    return conn


def require_database(db_path=None):
    """Fail with an actionable message instead of an empty-table traceback."""
    path = db_path or DB_PATH
    if not path.exists():
        raise DatabaseMissingError(
            f"No database at {path}. Build it first with `python load_data.py` "
            f"(or run the whole pipeline with `make pipeline`)."
        )
    return path


def query(sql, params=(), db_path=None):
    """Run a read-only query and return the result as a DataFrame."""
    require_database(db_path)
    with closing(connect(db_path)) as conn:
        return pd.read_sql_query(sql, conn, params=params)
