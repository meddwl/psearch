#!/usr/bin/env python3
# license         : BSD-3
# ==============================================================================

"""SQLite-backed storage for psearch compound databases.

Each database is a single .db file using only the stdlib sqlite3 module.
The public API is designed so that connections are opened and closed per
operation — RAM usage stays flat regardless of database size.

Public API
----------
init_db(db_path, bin_step, pharm_def=None)
save_model(db_path, mol_name, model)
load_model(db_path, mol_name)
delete_model(db_path, mol_name)
iter_models(db_path)
get_mol_names(db_path)
get_bin_step(db_path)
write_metadata(db_path, key, value)
read_metadata(db_path, key)
"""

import pickle
import sqlite3
from contextlib import contextmanager
from datetime import datetime, timezone
from typing import Generator, NamedTuple

import psearch


# ---------------------------------------------------------------------------
# Schema
# ---------------------------------------------------------------------------

_CREATE_METADATA = """
CREATE TABLE IF NOT EXISTS metadata (
    key   TEXT PRIMARY KEY,
    value TEXT NOT NULL
);
"""

_CREATE_COMPOUNDS = """
CREATE TABLE IF NOT EXISTS compounds (
    mol_name   TEXT    NOT NULL,
    stereo_id  INTEGER NOT NULL,
    mol_data   BLOB    NOT NULL,
    pharm_data BLOB    NOT NULL,
    fp_data    BLOB    NOT NULL,
    PRIMARY KEY (mol_name, stereo_id)
);
"""

_CREATE_INDEX = """
CREATE INDEX IF NOT EXISTS idx_compounds_mol_name
    ON compounds (mol_name);
"""


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

@contextmanager
def _connect(db_path: str) -> Generator[sqlite3.Connection, None, None]:
    """Open a SQLite connection in WAL mode, commit on success, rollback on error.

    Args:
        db_path: Path to the SQLite database file.

    Yields:
        An open sqlite3.Connection.
    """
    conn = sqlite3.connect(db_path)
    conn.execute("PRAGMA journal_mode=WAL")
    try:
        yield conn
        conn.commit()
    except Exception:
        conn.rollback()
        raise
    finally:
        conn.close()


def _ensure_schema(conn: sqlite3.Connection) -> None:
    """Create tables and index if they do not already exist.

    Args:
        conn: An open sqlite3.Connection.
    """
    conn.execute(_CREATE_METADATA)
    conn.execute(_CREATE_COMPOUNDS)
    conn.execute(_CREATE_INDEX)


# ---------------------------------------------------------------------------
# CompoundRecord
# ---------------------------------------------------------------------------

class CompoundRecord(NamedTuple):
    """All database data for one compound across all stereoisomers.

    Attributes:
        mol_dict: Dict mapping stereo_id (int) to an RDKit Mol with embedded conformers.
        pharm_dict: Dict mapping stereo_id (int) to a list of feature-coordinate lists,
                    one list per conformer: [[(label, (x, y, z)), ...], ...].
        fp_dict: Dict mapping stereo_id (int) to a list of pharmacophore fingerprints,
                 one fingerprint (frozenset[str] or dict) per conformer.
    """
    mol_dict: dict
    pharm_dict: dict
    fp_dict: dict


# ---------------------------------------------------------------------------
# Database initialisation
# ---------------------------------------------------------------------------

def init_db(
    db_path: str,
    bin_step: float,
    pharm_def: str | None = None,
) -> None:
    """Create a new psearch SQLite database and write initial metadata.

    Creates the schema (metadata and compounds tables) and records bin_step,
    psearch version, pharmacophore feature definition path, and creation
    timestamp. Raises sqlite3.OperationalError if the tables already exist and
    the caller used a pre-existing file.

    Args:
        db_path: Path to the new .db file to create.
        bin_step: Bin width (Å) used to discretise pharmacophore coordinates.
        pharm_def: Path to a custom pmapper feature-definition file, or None
                   for defaults.
    """
    with _connect(db_path) as conn:
        _ensure_schema(conn)
        conn.executemany(
            "INSERT OR REPLACE INTO metadata (key, value) VALUES (?, ?)",
            [
                ('bin_step', str(bin_step)),
                ('psearch_version', psearch.__version__),
                ('pharm_def', pharm_def if pharm_def is not None else ''),
                ('created_at', datetime.now(timezone.utc).isoformat()),
            ],
        )


# ---------------------------------------------------------------------------
# Metadata helpers
# ---------------------------------------------------------------------------

def get_bin_step(db_path: str) -> float:
    """Return the bin width (Å) stored in the database metadata.

    Args:
        db_path: Path to the psearch SQLite database file.

    Raises:
        KeyError: If 'bin_step' is not present in the metadata table.
    """
    with _connect(db_path) as conn:
        row = conn.execute(
            "SELECT value FROM metadata WHERE key = 'bin_step'"
        ).fetchone()
    if row is None:
        raise KeyError("'bin_step' not found in database metadata")
    return float(row[0])


def write_metadata(db_path: str, key: str, value: str) -> None:
    """Insert or replace a single metadata entry.

    Args:
        db_path: Path to the psearch SQLite database file.
        key: Metadata key string.
        value: Metadata value string.
    """
    with _connect(db_path) as conn:
        conn.execute(
            "INSERT OR REPLACE INTO metadata (key, value) VALUES (?, ?)",
            (key, value),
        )


def read_metadata(db_path: str, key: str) -> str:
    """Return the metadata value for key.

    Args:
        db_path: Path to the psearch SQLite database file.
        key: Metadata key string.

    Raises:
        KeyError: If key is not present in the metadata table.
    """
    with _connect(db_path) as conn:
        row = conn.execute(
            "SELECT value FROM metadata WHERE key = ?", (key,)
        ).fetchone()
    if row is None:
        raise KeyError(f"'{key}' not found in database metadata")
    return row[0]


# ---------------------------------------------------------------------------
# Compound CRUD
# ---------------------------------------------------------------------------

def save_model(db_path: str, mol_name: str, model: CompoundRecord) -> None:
    """Store all stereoisomer data for a compound, replacing any existing rows.

    Serialises mol_data, pharm_data, and fp_data with pickle (protocol 4)
    before storage. One row is written per stereo_id found in model.mol_dict.

    Args:
        db_path: Path to the psearch SQLite database file.
        mol_name: Compound identifier string.
        model: CompoundRecord with mol_dict, pharm_dict, and fp_dict keyed
               by stereo_id (int).
    """
    rows = [
        (
            mol_name,
            stereo_id,
            pickle.dumps(model.mol_dict[stereo_id], protocol=4),
            pickle.dumps(model.pharm_dict[stereo_id], protocol=4),
            pickle.dumps(model.fp_dict[stereo_id], protocol=4),
        )
        for stereo_id in model.mol_dict
    ]
    with _connect(db_path) as conn:
        conn.executemany(
            """INSERT OR REPLACE INTO compounds
               (mol_name, stereo_id, mol_data, pharm_data, fp_data)
               VALUES (?, ?, ?, ?, ?)""",
            rows,
        )


def load_model(db_path: str, mol_name: str) -> CompoundRecord:
    """Load all stereoisomer data for a compound.

    Args:
        db_path: Path to the psearch SQLite database file.
        mol_name: Compound identifier string.

    Returns:
        CompoundRecord with mol_dict, pharm_dict, and fp_dict keyed by
        stereo_id (int).

    Raises:
        KeyError: If mol_name is not present in the database.
    """
    with _connect(db_path) as conn:
        rows = conn.execute(
            "SELECT stereo_id, mol_data, pharm_data, fp_data "
            "FROM compounds WHERE mol_name = ?",
            (mol_name,),
        ).fetchall()
    if not rows:
        raise KeyError(f"Compound '{mol_name}' not found in database")
    mol_dict: dict = {}
    pharm_dict: dict = {}
    fp_dict: dict = {}
    for stereo_id, mol_blob, pharm_blob, fp_blob in rows:
        mol_dict[stereo_id] = pickle.loads(mol_blob)
        pharm_dict[stereo_id] = pickle.loads(pharm_blob)
        fp_dict[stereo_id] = pickle.loads(fp_blob)
    return CompoundRecord(mol_dict, pharm_dict, fp_dict)


def delete_model(db_path: str, mol_name: str) -> None:
    """Remove all rows for a compound from the database.

    No-op if mol_name does not exist.

    Args:
        db_path: Path to the psearch SQLite database file.
        mol_name: Compound identifier string.
    """
    with _connect(db_path) as conn:
        conn.execute(
            "DELETE FROM compounds WHERE mol_name = ?", (mol_name,)
        )


def iter_models(
    db_path: str,
) -> Generator[tuple[str, CompoundRecord], None, None]:
    """Yield (mol_name, CompoundRecord) for every compound in the database.

    Loads one compound at a time so that memory usage stays flat regardless
    of database size.

    Args:
        db_path: Path to the psearch SQLite database file.

    Yields:
        Tuples of (mol_name, CompoundRecord).
    """
    with _connect(db_path) as conn:
        mol_names = [
            row[0]
            for row in conn.execute(
                "SELECT DISTINCT mol_name FROM compounds ORDER BY mol_name"
            ).fetchall()
        ]
    for mol_name in mol_names:
        yield mol_name, load_model(db_path, mol_name)


def get_mol_names(db_path: str) -> tuple[str, ...]:
    """Return a tuple of all compound identifiers stored in the database.

    Args:
        db_path: Path to the psearch SQLite database file.
    """
    with _connect(db_path) as conn:
        rows = conn.execute(
            "SELECT DISTINCT mol_name FROM compounds ORDER BY mol_name"
        ).fetchall()
    return tuple(row[0] for row in rows)
