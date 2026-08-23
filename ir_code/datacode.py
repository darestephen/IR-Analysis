import sqlite3
from pathlib import Path

DB_PATH = Path(__file__).resolve().parent / "compound_spectra.db"

print("=" * 70)
print("DATABASE DIAGNOSTIC")
print("=" * 70)

print(f"\nDatabase location:")
print(DB_PATH)

print(f"\nDatabase exists: {DB_PATH.exists()}")

if not DB_PATH.exists():
    print("\n❌ Database does not exist.")
    raise SystemExit

conn = sqlite3.connect(DB_PATH)

# ------------------------------------------------------------
# SHOW TABLES
# ------------------------------------------------------------

tables = conn.execute("""
    SELECT name
    FROM sqlite_master
    WHERE type='table'
""").fetchall()

print("\nTables in database:")

for table in tables:
    print("  ", table[0])

# ------------------------------------------------------------
# COUNT ALL COMPOUNDS
# ------------------------------------------------------------

try:

    count = conn.execute("""
        SELECT COUNT(*)
        FROM compound_spectra
    """).fetchone()[0]

    print(f"\nTotal database records: {count}")

except Exception as error:

    print("\n❌ Could not count records.")
    print(error)

# ------------------------------------------------------------
# SHOW ALL COMPOUNDS
# ------------------------------------------------------------

try:

    rows = conn.execute("""
        SELECT
            rowid,
            compound_name,
            CASE
                WHEN ir_data IS NULL THEN 'NO'
                ELSE 'YES'
            END AS has_ir,
            CASE
                WHEN fft_data IS NULL THEN 'NO'
                ELSE 'YES'
            END AS has_fft
        FROM compound_spectra
        ORDER BY rowid
    """).fetchall()

    print("\n" + "=" * 70)
    print("ALL DATABASE COMPOUNDS")
    print("=" * 70)

    print(
        f"{'ID':<6}"
        f"{'Compound':<35}"
        f"{'IR':<10}"
        f"{'FFT':<10}"
    )

    print("-" * 70)

    for rowid, name, has_ir, has_fft in rows:

        print(
            f"{rowid:<6}"
            f"{str(name):<35}"
            f"{has_ir:<10}"
            f"{has_fft:<10}"
        )

    # --------------------------------------------------------
    # COUNT USABLE
    # --------------------------------------------------------

    usable = conn.execute("""
        SELECT COUNT(*)
        FROM compound_spectra
        WHERE ir_data IS NOT NULL
        AND fft_data IS NOT NULL
    """).fetchone()[0]

    print("\n" + "=" * 70)

    print(
        f"Total records: {len(rows)}"
    )

    print(
        f"Records with IR + FFT: {usable}"
    )

    print("=" * 70)

except Exception as error:

    print("\n❌ Could not read compound records.")
    print(error)

conn.close()