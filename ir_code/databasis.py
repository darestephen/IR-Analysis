import sqlite3
from pathlib import Path
import pandas as pd


# ============================================================
# CENTRAL DATABASE
# ============================================================

DB_PATH = (
    Path.home()
    / "IR_Spectra"
    / "compound_spectra.db"
)


# ============================================================
# CREATE DATABASE
# ============================================================

def create_database(db_path=DB_PATH):

    db_path = Path(db_path)

    db_path.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    conn = sqlite3.connect(
        str(db_path)
    )

    cursor = conn.cursor()

    cursor.execute(
        """
        CREATE TABLE IF NOT EXISTS compound_spectra
        (
            id INTEGER PRIMARY KEY AUTOINCREMENT,

            compound_name TEXT NOT NULL,

            ir_data TEXT NOT NULL,

            fft_data TEXT NOT NULL,

            created_at TIMESTAMP
                DEFAULT CURRENT_TIMESTAMP
        )
        """
    )

    conn.commit()
    conn.close()

    return db_path


# ============================================================
# SAVE IR + FFT
# ============================================================

def save_spectra_to_db(
    db_path,
    compound_name,
    ir_df,
    fft_df
):

    try:

        db_path = Path(db_path)

        db_path.parent.mkdir(
            parents=True,
            exist_ok=True
        )

        # ----------------------------------------------------
        # Validate
        # ----------------------------------------------------

        if ir_df is None:
            raise ValueError(
                "IR DataFrame is None."
            )

        if fft_df is None:
            raise ValueError(
                "FFT DataFrame is None."
            )

        if len(ir_df) == 0:
            raise ValueError(
                "IR DataFrame is empty."
            )

        if len(fft_df) == 0:
            raise ValueError(
                "FFT DataFrame is empty."
            )

        compound_name = str(
            compound_name
        ).strip()

        if not compound_name:

            raise ValueError(
                "Compound name is empty."
            )

        # ----------------------------------------------------
        # Convert DataFrames to CSV strings
        # ----------------------------------------------------

        ir_csv = ir_df.to_csv(
            index=False
        )

        fft_csv = fft_df.to_csv(
            index=False
        )

        # ----------------------------------------------------
        # Create database
        # ----------------------------------------------------

        create_database(
            db_path
        )

        # ----------------------------------------------------
        # Connect
        # ----------------------------------------------------

        conn = sqlite3.connect(
            str(db_path)
        )

        cursor = conn.cursor()

        # ----------------------------------------------------
        # Insert
        #
        # DO NOT replace existing records.
        # Each spectrum receives its own ID.
        # ----------------------------------------------------

        cursor.execute(
            """
            INSERT INTO compound_spectra
            (
                compound_name,
                ir_data,
                fft_data
            )
            VALUES (?, ?, ?)
            """,
            (
                compound_name,
                ir_csv,
                fft_csv
            )
        )

        record_id = cursor.lastrowid

        conn.commit()
        conn.close()

        # ----------------------------------------------------
        # Verify
        # ----------------------------------------------------

        print()
        print("=" * 80)
        print("DATABASE SAVE SUCCESSFUL")
        print("=" * 80)

        print(
            f"Database : {db_path}"
        )

        print(
            f"ID       : {record_id}"
        )

        print(
            f"Compound : {compound_name}"
        )

        print(
            f"IR points: {len(ir_df)}"
        )

        print(
            f"FFT points: {len(fft_df)}"
        )

        print("=" * 80)

        return record_id

    except Exception as error:

        print()
        print("=" * 80)
        print("❌ DATABASE SAVE FAILED")
        print("=" * 80)

        print(
            f"{type(error).__name__}: {error}"
        )

        print("=" * 80)

        return None


# ============================================================
# LOAD ALL SPECTRA
# ============================================================

def load_database(
    db_path=DB_PATH
):

    db_path = Path(db_path)

    if not db_path.exists():

        print(
            f"❌ Database does not exist:"
        )

        print(
            db_path
        )

        return []

    conn = sqlite3.connect(
        str(db_path)
    )

    cursor = conn.cursor()

    cursor.execute(
        """
        SELECT
            id,
            compound_name,
            ir_data,
            fft_data,
            created_at
        FROM compound_spectra
        ORDER BY id
        """
    )

    rows = cursor.fetchall()

    conn.close()

    return rows


# ============================================================
# LOAD USABLE SPECTRA
# ============================================================

def load_usable_spectra(
    db_path=DB_PATH
):

    db_path = Path(db_path)

    if not db_path.exists():

        return []

    conn = sqlite3.connect(
        str(db_path)
    )

    cursor = conn.cursor()

    cursor.execute(
        """
        SELECT
            id,
            compound_name,
            ir_data,
            fft_data,
            created_at
        FROM compound_spectra
        WHERE ir_data IS NOT NULL
          AND ir_data != ''
          AND fft_data IS NOT NULL
          AND fft_data != ''
        ORDER BY id
        """
    )

    rows = cursor.fetchall()

    conn.close()

    return rows


# ============================================================
# SHOW DATABASE
# ============================================================

def show_database(
    db_path=DB_PATH
):

    db_path = Path(db_path)

    print()
    print("=" * 90)
    print("IR SPECTRA DATABASE")
    print("=" * 90)

    print(
        f"Database:"
    )

    print(
        db_path
    )

    print()

    if not db_path.exists():

        print(
            "❌ Database does not exist."
        )

        print("=" * 90)

        return

    rows = load_database(
        db_path
    )

    if not rows:

        print(
            "❌ Database table is empty."
        )

        print("=" * 90)

        return

    print(
        f"{'ID':<8}"
        f"{'Compound':<40}"
        f"{'IR':<10}"
        f"{'FFT':<10}"
    )

    print("-" * 90)

    for row in rows:

        (
            record_id,
            compound_name,
            ir_data,
            fft_data,
            created_at
        ) = row

        ir_status = (
            "YES"
            if ir_data
            else "NO"
        )

        fft_status = (
            "YES"
            if fft_data
            else "NO"
        )

        print(
            f"{record_id:<8}"
            f"{str(compound_name)[:38]:<40}"
            f"{ir_status:<10}"
            f"{fft_status:<10}"
        )

    print("-" * 90)

    print(
        f"Total records: {len(rows)}"
    )

    print("=" * 90)


# ============================================================
# DATABASE STATISTICS
# ============================================================

def database_statistics(
    db_path=DB_PATH
):

    rows = load_database(
        db_path
    )

    total = len(rows)

    usable = sum(
        1
        for row in rows
        if row[2]
        and row[3]
    )

    compounds = set(
        row[1]
        for row in rows
    )

    print()
    print("=" * 70)
    print("DATABASE STATISTICS")
    print("=" * 70)

    print(
        f"Database         : {db_path}"
    )

    print(
        f"Total records    : {total}"
    )

    print(
        f"Usable IR + FFT  : {usable}"
    )

    print(
        f"Unique compounds : {len(compounds)}"
    )

    print("=" * 70)


# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":

    print()
    print("=" * 80)
    print("IR SPECTRA DATABASE MANAGER")
    print("=" * 80)

    db = create_database()

    print()
    print(
        "Database created/verified:"
    )

    print(
        db
    )

    database_statistics(
        db
    )

    show_database(
        db
    )

    print()
    print(
        "Database manager finished."
    )