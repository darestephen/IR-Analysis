import sqlite3
from pathlib import Path


# ============================================================
# CENTRAL DATABASE LOCATION
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

    # Make sure IR_Spectra directory exists
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
        CREATE TABLE IF NOT EXISTS compound_spectra (

            id INTEGER PRIMARY KEY AUTOINCREMENT,

            compound_name TEXT NOT NULL,

            ir_data TEXT,

            fft_data TEXT,

            created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
        )
        """
    )

    conn.commit()

    conn.close()

    return db_path


# ============================================================
# SAVE SPECTRA
# ============================================================

def save_spectra_to_db(
    db_path,
    compound_name,
    ir_df,
    fft_df
):

    try:

        db_path = Path(db_path)

        # ----------------------------------------------------
        # Make sure database directory exists
        # ----------------------------------------------------

        db_path.parent.mkdir(
            parents=True,
            exist_ok=True
        )

        # ----------------------------------------------------
        # Convert DataFrames to CSV
        # ----------------------------------------------------

        ir_csv_str = ir_df.to_csv(
            index=False
        )

        fft_csv_str = fft_df.to_csv(
            index=False
        )

        # ----------------------------------------------------
        # Connect
        # ----------------------------------------------------

        conn = sqlite3.connect(
            str(db_path)
        )

        cursor = conn.cursor()

        # ----------------------------------------------------
        # Create table
        # ----------------------------------------------------

        cursor.execute(
            """
            CREATE TABLE IF NOT EXISTS compound_spectra (

                id INTEGER PRIMARY KEY AUTOINCREMENT,

                compound_name TEXT NOT NULL,

                ir_data TEXT,

                fft_data TEXT,

                created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            )
            """
        )

        # ----------------------------------------------------
        # INSERT NEW RECORD
        #
        # IMPORTANT:
        # We do NOT use INSERT OR REPLACE.
        #
        # Every spectrum gets its own ID.
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
                ir_csv_str,
                fft_csv_str
            )
        )

        # ----------------------------------------------------
        # Get new database ID
        # ----------------------------------------------------

        record_id = cursor.lastrowid

        # ----------------------------------------------------
        # Commit
        # ----------------------------------------------------

        conn.commit()

        conn.close()

        # ----------------------------------------------------
        # Success message
        # ----------------------------------------------------

        print()
        print("=" * 70)
        print("DATABASE SAVE SUCCESSFUL")
        print("=" * 70)

        print(
            f"Database ID : {record_id}"
        )

        print(
            f"Compound    : {compound_name}"
        )

        print(
            f"Database    : {db_path}"
        )

        print("=" * 70)

        return True

    except Exception as error:

        print()
        print("=" * 70)
        print("DATABASE SAVE FAILED")
        print("=" * 70)

        print(
            f"Error type: {type(error).__name__}"
        )

        print(
            f"Error: {error}"
        )

        print("=" * 70)

        return False


# ============================================================
# LIST DATABASE RECORDS
# ============================================================

def list_database_records(
    db_path=DB_PATH
):

    try:

        db_path = Path(db_path)

        if not db_path.exists():

            print()
            print(
                "❌ Database does not exist:"
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

    except Exception as error:

        print()
        print(
            "❌ Could not read database."
        )

        print(
            f"Error: {error}"
        )

        return []


# ============================================================
# DISPLAY DATABASE
# ============================================================

def show_database(
    db_path=DB_PATH
):

    rows = list_database_records(
        db_path
    )

    print()
    print("=" * 90)

    print(
        "DATABASE CONTENTS"
    )

    print("=" * 90)

    print(
        f"Database:"
    )

    print(
        db_path
    )

    print()

    if not rows:

        print(
            "❌ No records found."
        )

        print("=" * 90)

        return

    print(
        f"{'ID':<6}"
        f"{'Compound':<35}"
        f"{'IR':<10}"
        f"{'FFT':<10}"
    )

    print("-" * 90)

    for (
        record_id,
        compound_name,
        ir_data,
        fft_data,
        created_at
    ) in rows:

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
            f"{record_id:<6}"
            f"{compound_name:<35}"
            f"{ir_status:<10}"
            f"{fft_status:<10}"
        )

    print("-" * 90)

    print(
        f"Total records: {len(rows)}"
    )

    print("=" * 90)


# ============================================================
# TEST DATABASE
# ============================================================

if __name__ == "__main__":

    db = create_database()

    print()
    print(
        "Database location:"
    )

    print(
        db
    )

    show_database(
        db
    )