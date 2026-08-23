import sqlite3
import pandas as pd
import numpy as np
import requests

from pathlib import Path
from io import StringIO, BytesIO
from urllib.parse import quote

import tkinter as tk
from tkinter import filedialog, messagebox

from PIL import Image, ImageTk


# ============================================================
# DATABASE LOCATION
# ============================================================

# Try to use the same DB_PATH used by databasecode.py.
#
# This is important because Fourier.py and Steve.py should
# use the SAME database.

try:

    from databasis import DB_PATH

except ImportError:

    # Fallback location.
    #
    # This matches the location used by the current Fourier.py:
    #
    # ~/project/compound_spectra.db

    DB_PATH = (
        Path.home()
        / "IR_Spectra"
        / "compound_spectra.db"
    )


DB_PATH = Path(DB_PATH).resolve()


APPLICATION_FOLDER = (
    Path.home()
    / "IR_Spectra"
)

RESULTS_FOLDER = (
    APPLICATION_FOLDER
    / "identification_results"
)

try:

    APPLICATION_FOLDER.mkdir(
        parents=True,
        exist_ok=True
    )

    RESULTS_FOLDER.mkdir(
        parents=True,
        exist_ok=True
    )

except Exception as error:

    print(
        f"⚠️ Could not create application folders: {error}"
    )


def select_file(title):

    """
    Open a Tkinter file-selection dialog.

    Returns:
        pathlib.Path
        or None if the user cancels.
    """

    root = None

    try:

        root = tk.Tk()

        root.withdraw()

        try:
            root.attributes(
                "-topmost",
                True
            )
        except Exception:
            pass

        file_path = filedialog.askopenfilename(
            title=title,
            filetypes=[
                (
                    "CSV files",
                    "*.csv"
                ),
                (
                    "All files",
                    "*.*"
                )
            ]
        )

        root.destroy()

        if not file_path:

            return None

        return Path(
            file_path
        ).resolve()

    except Exception as error:

        if root is not None:

            try:
                root.destroy()
            except Exception:
                pass

        print()
        print(
            "❌ Could not open file-selection window."
        )

        print(
            f"Error: {error}"
        )

        return None


# ============================================================
# IDENTIFY CSV TYPE
# ============================================================

def identify_csv_type(file_path):

    """
    Determine whether a CSV file is an IR or FFT file.
    """

    file_path = Path(
        file_path
    )

    try:

        df = pd.read_csv(
            file_path,
            nrows=5
        )

    except Exception as error:

        raise ValueError(
            f"Could not read CSV file: {error}"
        )

    if df.empty:

        return "UNKNOWN"

    columns = [
        str(column).strip().lower()
        for column in df.columns
    ]

    # --------------------------------------------------------
    # FFT indicators
    # --------------------------------------------------------

    fft_indicators = [
        "frequency index",
        "fft magnitude",
        "fft"
    ]

    fft_found = any(
        indicator in column
        for column in columns
        for indicator in fft_indicators
    )

    # --------------------------------------------------------
    # IR indicators
    # --------------------------------------------------------

    ir_indicators = [
        "wavenumber",
        "transmittance",
        "absorbance"
    ]

    ir_found = any(
        indicator in column
        for column in columns
        for indicator in ir_indicators
    )

    if fft_found and not ir_found:

        return "FFT"

    if ir_found and not fft_found:

        return "IR"

    # --------------------------------------------------------
    # Additional checks
    # --------------------------------------------------------

    if (
        "frequency index" in columns
        and "fft magnitude" in columns
    ):

        return "FFT"

    if (
        any(
            "wavenumber" in column
            for column in columns
        )
        and
        any(
            (
                "transmittance" in column
                or
                "absorbance" in column
            )
            for column in columns
        )
    ):

        return "IR"

    return "UNKNOWN"


# ============================================================
# FIND NUMERIC COLUMNS
# ============================================================

def find_numeric_columns(df):

    """
    Return columns containing numerical data.
    """

    numeric_columns = []

    for column in df.columns:

        values = pd.to_numeric(
            df[column],
            errors="coerce"
        )

        if values.notna().sum() > 0:

            numeric_columns.append(
                column
            )

    return numeric_columns


# ============================================================
# LOAD CSV
# ============================================================

def load_csv(
    file_path,
    expected_type=None
):

    """
    Load an IR or FFT CSV file.

    Returns:
        x, y
    """

    file_path = Path(
        file_path
    )

    try:

        df = pd.read_csv(
            file_path
        )

    except Exception as error:

        raise ValueError(
            f"Could not read {file_path.name}: {error}"
        )

    if df.empty:

        raise ValueError(
            f"{file_path.name} is empty."
        )

    print()
    print(
        f"Columns found in {file_path.name}:"
    )

    print(
        list(df.columns)
    )

    detected_type = identify_csv_type(
        file_path
    )

    print(
        f"Detected file type: {detected_type}"
    )

    # --------------------------------------------------------
    # Check expected type
    # --------------------------------------------------------

    if expected_type is not None:

        if detected_type != expected_type:

            raise ValueError(
                f"{file_path.name} appears to be "
                f"{detected_type}, but "
                f"{expected_type} was expected."
            )

    # --------------------------------------------------------
    # Column dictionary
    # --------------------------------------------------------

    lower_columns = {
        str(column).strip().lower(): column
        for column in df.columns
    }

    x_column = None
    y_column = None

    # ========================================================
    # IR
    # ========================================================

    if expected_type == "IR":

        # X = wavenumber

        for (
            column_lower,
            original_column
        ) in lower_columns.items():

            if "wavenumber" in column_lower:

                x_column = original_column

                break

        # Y = transmittance or absorbance

        for (
            column_lower,
            original_column
        ) in lower_columns.items():

            if (
                "transmittance" in column_lower
                or
                "absorbance" in column_lower
            ):

                y_column = original_column

                break

    # ========================================================
    # FFT
    # ========================================================

    elif expected_type == "FFT":

        # X = frequency index

        for (
            column_lower,
            original_column
        ) in lower_columns.items():

            if "frequency index" in column_lower:

                x_column = original_column

                break

        # Y = FFT magnitude

        for (
            column_lower,
            original_column
        ) in lower_columns.items():

            if (
                "fft magnitude" in column_lower
                or
                "magnitude" in column_lower
            ):

                y_column = original_column

                break

    # ========================================================
    # FALLBACK
    # ========================================================

    if (
        x_column is None
        or
        y_column is None
    ):

        numeric_columns = find_numeric_columns(
            df
        )

        if len(numeric_columns) < 2:

            raise ValueError(
                f"{file_path.name} does not contain "
                "two numerical columns."
            )

        x_column = numeric_columns[0]

        y_column = numeric_columns[1]

    print(
        f"Using '{x_column}' as X."
    )

    print(
        f"Using '{y_column}' as Y."
    )

    # ========================================================
    # CONVERT TO NUMERIC
    # ========================================================

    df[x_column] = pd.to_numeric(
        df[x_column],
        errors="coerce"
    )

    df[y_column] = pd.to_numeric(
        df[y_column],
        errors="coerce"
    )

    df = df.replace(
        [np.inf, -np.inf],
        np.nan
    )

    df = df.dropna(
        subset=[
            x_column,
            y_column
        ]
    )

    # ========================================================
    # CHECK DATA SIZE
    # ========================================================

    if len(df) < 10:

        raise ValueError(
            f"{file_path.name} contains too few "
            "valid data points."
        )

    # ========================================================
    # SORT BY X
    # ========================================================

    df = df.sort_values(
        by=x_column
    )

    # Remove duplicate x values.

    df = df.drop_duplicates(
        subset=[x_column]
    )

    # ========================================================
    # ARRAYS
    # ========================================================

    x = df[x_column].to_numpy(
        dtype=float
    )

    y = df[y_column].to_numpy(
        dtype=float
    )

    return x, y


# ============================================================
# LOAD DATABASE CSV
# ============================================================

def load_database_csv(csv_text):

    """
    Convert CSV text stored in SQLite into x/y arrays.
    """

    if csv_text is None:

        raise ValueError(
            "Database spectrum is empty."
        )

    try:

        df = pd.read_csv(
            StringIO(
                str(csv_text)
            )
        )

    except Exception as error:

        raise ValueError(
            f"Could not decode database CSV: {error}"
        )

    if df.empty:

        raise ValueError(
            "Database spectrum contains no data."
        )

    numeric_columns = find_numeric_columns(
        df
    )

    if len(numeric_columns) < 2:

        raise ValueError(
            "Database spectrum does not contain "
            "two numerical columns."
        )

    x_column = numeric_columns[0]

    y_column = numeric_columns[1]

    df[x_column] = pd.to_numeric(
        df[x_column],
        errors="coerce"
    )

    df[y_column] = pd.to_numeric(
        df[y_column],
        errors="coerce"
    )

    df = df.replace(
        [np.inf, -np.inf],
        np.nan
    )

    df = df.dropna(
        subset=[
            x_column,
            y_column
        ]
    )

    if len(df) < 10:

        raise ValueError(
            "Database spectrum contains too few "
            "valid data points."
        )

    df = df.sort_values(
        by=x_column
    )

    df = df.drop_duplicates(
        subset=[x_column]
    )

    x = df[x_column].to_numpy(
        dtype=float
    )

    y = df[y_column].to_numpy(
        dtype=float
    )

    return x, y


# ============================================================
# NORMALIZE ARRAY
# ============================================================

def normalize(values):

    """
    Standard-score normalization.

    This removes differences in absolute intensity
    while preserving the shape of the spectrum.
    """

    values = np.asarray(
        values,
        dtype=float
    )

    if len(values) == 0:

        return values

    mean = np.mean(
        values
    )

    std = np.std(
        values
    )

    centered = (
        values - mean
    )

    if std < 1e-12:

        return centered

    return (
        centered / std
    )


# ============================================================
# RESAMPLE USING ACTUAL X VALUES
# ============================================================

def resample_xy(
    x,
    y,
    new_x
):

    """
    Interpolate y values using the actual x coordinates.
    """

    x = np.asarray(
        x,
        dtype=float
    )

    y = np.asarray(
        y,
        dtype=float
    )

    new_x = np.asarray(
        new_x,
        dtype=float
    )

    # Remove invalid values.

    valid = (
        np.isfinite(x)
        &
        np.isfinite(y)
    )

    x = x[valid]

    y = y[valid]

    if len(x) < 2:

        raise ValueError(
            "Not enough valid points for interpolation."
        )

    # Sort.

    order = np.argsort(
        x
    )

    x = x[order]

    y = y[order]

    # Remove duplicate x values.

    x_unique, unique_indices = np.unique(
        x,
        return_index=True
    )

    y_unique = y[
        unique_indices
    ]

    if len(x_unique) < 2:

        raise ValueError(
            "Not enough unique X values."
        )

    return np.interp(
        new_x,
        x_unique,
        y_unique
    )


# ============================================================
# PREPARE IR SPECTRA FOR COMPARISON
# ============================================================

def prepare_ir_for_comparison(
    unknown_x,
    unknown_y,
    database_x,
    database_y,
    points=2000
):

    """
    Put two IR spectra onto the same wavenumber grid.

    The actual wavenumber values are used.

    This is much better than simply resampling both
    spectra according to array position.
    """

    unknown_x = np.asarray(
        unknown_x,
        dtype=float
    )

    database_x = np.asarray(
        database_x,
        dtype=float
    )

    # --------------------------------------------------------
    # Determine overlapping wavenumber range
    # --------------------------------------------------------

    lower_limit = max(
        np.min(unknown_x),
        np.min(database_x)
    )

    upper_limit = min(
        np.max(unknown_x),
        np.max(database_x)
    )

    if upper_limit <= lower_limit:

        raise ValueError(
            "The unknown and database IR spectra "
            "have no overlapping wavenumber range."
        )

    # --------------------------------------------------------
    # Common wavenumber grid
    # --------------------------------------------------------

    common_x = np.linspace(
        lower_limit,
        upper_limit,
        points
    )

    unknown_interp = resample_xy(
        unknown_x,
        unknown_y,
        common_x
    )

    database_interp = resample_xy(
        database_x,
        database_y,
        common_x
    )

    # --------------------------------------------------------
    # Normalize
    # --------------------------------------------------------

    unknown_norm = normalize(
        unknown_interp
    )

    database_norm = normalize(
        database_interp
    )

    return (
        common_x,
        unknown_norm,
        database_norm
    )


# ============================================================
# PREPARE FFT SPECTRA
# ============================================================

def prepare_fft_for_comparison(
    unknown_x,
    unknown_y,
    database_x,
    database_y,
    points=1000
):

    """
    Put two FFT spectra onto a common normalized
    frequency-index axis.
    """

    unknown_x = np.asarray(
        unknown_x,
        dtype=float
    )

    database_x = np.asarray(
        database_x,
        dtype=float
    )

    # --------------------------------------------------------
    # Normalize frequency axis to 0-1
    # --------------------------------------------------------

    def normalize_x(x):

        minimum = np.min(
            x
        )

        maximum = np.max(
            x
        )

        if maximum == minimum:

            return np.zeros_like(
                x
            )

        return (
            (x - minimum)
            /
            (maximum - minimum)
        )

    unknown_x_norm = normalize_x(
        unknown_x
    )

    database_x_norm = normalize_x(
        database_x
    )

    # --------------------------------------------------------
    # Common frequency grid
    # --------------------------------------------------------

    common_x = np.linspace(
        0,
        1,
        points
    )

    unknown_interp = resample_xy(
        unknown_x_norm,
        unknown_y,
        common_x
    )

    database_interp = resample_xy(
        database_x_norm,
        database_y,
        common_x
    )

    # --------------------------------------------------------
    # Normalize amplitudes
    # --------------------------------------------------------

    unknown_norm = normalize(
        unknown_interp
    )

    database_norm = normalize(
        database_interp
    )

    return (
        common_x,
        unknown_norm,
        database_norm
    )


# ============================================================
# PEARSON CORRELATION
# ============================================================

def calculate_correlation(
    values_a,
    values_b
):

    """
    Calculate Pearson correlation.

    Returns a value from approximately -1 to +1.
    """

    a = np.asarray(
        values_a,
        dtype=float
    )

    b = np.asarray(
        values_b,
        dtype=float
    )

    if len(a) != len(b):

        raise ValueError(
            "Arrays must have equal length."
        )

    if len(a) < 2:

        return 0.0

    std_a = np.std(
        a
    )

    std_b = np.std(
        b
    )

    if (
        std_a < 1e-12
        or
        std_b < 1e-12
    ):

        return 0.0

    correlation = np.corrcoef(
        a,
        b
    )[0, 1]

    if not np.isfinite(
        correlation
    ):

        return 0.0

    return float(
        correlation
    )


# ============================================================
# CORRELATION TO PERCENTAGE
# ============================================================

def correlation_to_score(
    correlation
):

    """
    Convert correlation from -1...+1
    to 0...100%.

    +1 = 100%
     0 = 50%
    -1 = 0%
    """

    correlation = float(
        correlation
    )

    correlation = max(
        -1.0,
        min(
            1.0,
            correlation
        )
    )

    score = (
        (correlation + 1.0)
        / 2.0
        * 100.0
    )

    return score


# ============================================================
# CALCULATE SPECTRAL MATCH
# ============================================================

def calculate_spectral_match(
    unknown_ir_x,
    unknown_ir_y,
    unknown_fft_x,
    unknown_fft_y,
    database_ir_x,
    database_ir_y,
    database_fft_x,
    database_fft_y
):

    """
    Calculate IR similarity, FFT similarity,
    and combined similarity.

    IR receives 70% of the final score.
    FFT receives 30%.

    This weighting is intentional because the IR
    spectrum contains the primary chemical information.
    """

    # ========================================================
    # IR COMPARISON
    # ========================================================

    (
        ir_grid,
        unknown_ir,
        database_ir
    ) = prepare_ir_for_comparison(
        unknown_ir_x,
        unknown_ir_y,
        database_ir_x,
        database_ir_y
    )

    ir_correlation = calculate_correlation(
        unknown_ir,
        database_ir
    )

    ir_score = correlation_to_score(
        ir_correlation
    )

    # ========================================================
    # FFT COMPARISON
    # ========================================================

    (
        fft_grid,
        unknown_fft,
        database_fft
    ) = prepare_fft_for_comparison(
        unknown_fft_x,
        unknown_fft_y,
        database_fft_x,
        database_fft_y
    )

    fft_correlation = calculate_correlation(
        unknown_fft,
        database_fft
    )

    fft_score = correlation_to_score(
        fft_correlation
    )

    # ========================================================
    # COMBINED SCORE
    # ========================================================

    combined_score = (
        0.70 * ir_score
        +
        0.30 * fft_score
    )

    # --------------------------------------------------------
    # Additional RMSE information
    # --------------------------------------------------------

    ir_rmse = float(
        np.sqrt(
            np.mean(
                (
                    unknown_ir
                    -
                    database_ir
                ) ** 2
            )
        )
    )

    fft_rmse = float(
        np.sqrt(
            np.mean(
                (
                    unknown_fft
                    -
                    database_fft
                ) ** 2
            )
        )
    )

    return {
        "ir_correlation": ir_correlation,
        "ir_score": ir_score,
        "fft_correlation": fft_correlation,
        "fft_score": fft_score,
        "combined_score": combined_score,
        "ir_rmse": ir_rmse,
        "fft_rmse": fft_rmse
    }


# ============================================================
# GET PUBCHEM INFORMATION
# ============================================================

def get_pubchem_information(
    compound_name
):

    """
    Search PubChem for the identified compound.

    Returns:
        dictionary
        or None
    """

    print()
    print(
        f"Searching PubChem for: {compound_name}"
    )

    encoded_name = quote(
        str(compound_name),
        safe=""
    )

    properties_url = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/"
        f"compound/name/{encoded_name}/"
        "property/Title,MolecularFormula,CanonicalSMILES/"
        "JSON"
    )

    try:

        response = requests.get(
            properties_url,
            timeout=15
        )

        response.raise_for_status()

        data = response.json()

        properties = (
            data
            .get(
                "PropertyTable",
                {}
            )
            .get(
                "Properties",
                []
            )
        )

        if not properties:

            print(
                "❌ PubChem did not return a compound."
            )

            return None

        compound = properties[0]

        cid = compound.get(
            "CID"
        )

        title = compound.get(
            "Title",
            compound_name
        )

        molecular_formula = compound.get(
            "MolecularFormula",
            "Not available"
        )

        smiles = compound.get(
            "ConnectivitySMILES"
        )

        if not smiles:

            smiles = compound.get(
                "CanonicalSMILES",
                "Not available"
            )


        image_data = None

        if cid is not None:

            image_url = (
                "https://pubchem.ncbi.nlm.nih.gov/rest/pug/"
                f"compound/cid/{cid}/PNG"
                "?image_size=300x300"
            )

            image_response = requests.get(
                image_url,
                timeout=15
            )

            image_response.raise_for_status()

            image_data = (
                image_response.content
            )

        print(
            f"✅ PubChem CID: {cid}"
        )

        return {
            "cid": cid,
            "title": title,
            "formula": molecular_formula,
            "smiles": smiles,
            "image_data": image_data
        }

    except requests.RequestException as error:

        print()
        print(
            "❌ PubChem network error."
        )

        print(
            f"Error: {error}"
        )

        return None

    except Exception as error:

        print()
        print(
            "❌ Could not retrieve PubChem information."
        )

        print(
            f"Error: {error}"
        )

        return None


# ============================================================
# DISPLAY CHEMICAL STRUCTURE
# ============================================================

def show_structure(
    compound_name,
    pubchem_data
):

    """
    Display the chemical structure using Tkinter.
    """

    if pubchem_data is None:

        return

    image_data = pubchem_data.get(
        "image_data"
    )

    if image_data is None:

        print(
            "⚠️ No structure image was returned."
        )

        return

    try:

        image = Image.open(
            BytesIO(
                image_data
            )
        )

    except Exception as error:

        print()
        print(
            "❌ Could not open chemical structure image."
        )

        print(
            f"Error: {error}"
        )

        return

    # ========================================================
    # TKINTER WINDOW
    # ========================================================

    try:

        root = tk.Tk()

        root.title(
            f"Chemical Structure - {compound_name}"
        )

        root.geometry(
            "550x700"
        )

        root.resizable(
            False,
            False
        )

        # ====================================================
        # TITLE
        # ====================================================

        name_label = tk.Label(
            root,
            text=compound_name,
            font=(
                "Arial",
                18,
                "bold"
            )
        )

        name_label.pack(
            pady=10
        )

        # ====================================================
        # STRUCTURE IMAGE
        # ====================================================

        photo = ImageTk.PhotoImage(
            image
        )

        image_label = tk.Label(
            root,
            image=photo
        )

        image_label.image = photo

        image_label.pack(
            pady=10
        )

        # ====================================================
        # CID
        # ====================================================

        cid_label = tk.Label(
            root,
            text=(
                f"PubChem CID: "
                f"{pubchem_data.get('cid', 'N/A')}"
            ),
            font=(
                "Arial",
                12
            )
        )

        cid_label.pack(
            pady=5
        )

        # ====================================================
        # FORMULA
        # ====================================================

        formula_label = tk.Label(
            root,
            text=(
                "Molecular formula: "
                f"{pubchem_data.get('formula', 'N/A')}"
            ),
            font=(
                "Arial",
                12
            )
        )

        formula_label.pack(
            pady=5
        )

        # ====================================================
        # SMILES
        # ====================================================

        smiles_label = tk.Label(
            root,
            text=(
                "SMILES:\n"
                f"{pubchem_data.get('smiles', 'N/A')}"
            ),
            wraplength=500,
            justify="center",
            font=(
                "Arial",
                10
            )
        )

        smiles_label.pack(
            pady=10
        )

        # ====================================================
        # CLOSE BUTTON
        # ====================================================

        close_button = tk.Button(
            root,
            text="Close",
            command=root.destroy,
            width=15
        )

        close_button.pack(
            pady=15
        )

        # ====================================================
        # KEEP WINDOW OPEN
        # ====================================================

        root.mainloop()

    except Exception as error:

        print()
        print(
            "❌ Could not display structure window."
        )

        print(
            f"Error: {error}"
        )


# ============================================================
# CHECK DATABASE
# ============================================================

def check_database():

    """
    Verify that the SQLite database exists and
    contains the expected table.
    """

    print()
    print(
        f"Database location:"
    )

    print(
        DB_PATH
    )

    if not DB_PATH.exists():

        print()
        print(
            "❌ Database file not found."
        )

        print(
            f"Expected:"
        )

        print(
            DB_PATH
        )

        return False

    try:

        conn = sqlite3.connect(
            str(DB_PATH)
        )

        tables = conn.execute(
            """
            SELECT name
            FROM sqlite_master
            WHERE type='table'
            """
        ).fetchall()

        conn.close()

        table_names = [
            row[0]
            for row in tables
        ]

        print()
        print(
            "Database tables:"
        )

        print(
            table_names
        )

        if "compound_spectra" not in table_names:

            print()
            print(
                "❌ The database does not contain "
                "the 'compound_spectra' table."
            )

            return False

        return True

    except Exception as error:

        print()
        print(
            "❌ Could not inspect database."
        )

        print(
            f"Error: {error}"
        )

        return False


# ============================================================
# READ DATABASE SPECTRA
# ============================================================

def read_database_spectra():

    """
    Read all compounds that have both IR and FFT data.
    """

    try:

        conn = sqlite3.connect(
            str(DB_PATH)
        )

        rows = conn.execute(
            """
            SELECT
                compound_name,
                ir_data,
                fft_data
            FROM compound_spectra
            WHERE ir_data IS NOT NULL
            AND fft_data IS NOT NULL
            """
        ).fetchall()

        conn.close()

        return rows

    except Exception as error:

        print()
        print(
            "❌ Failed to read database."
        )

        print(
            f"Error: {error}"
        )

        return []


# ============================================================
# IDENTIFY UNKNOWN COMPOUND
# ============================================================

def identify_compound():

    """
    Main unknown-compound identification routine.
    """

    print()
    print("=" * 80)

    print(
        "        IR + FFT UNKNOWN COMPOUND IDENTIFICATION"
    )

    print("=" * 80)

    # ========================================================
    # DATABASE
    # ========================================================

    if not check_database():

        return None

    # ========================================================
    # SELECT UNKNOWN IR
    # ========================================================

    print()
    print(
        "Select the UNKNOWN IR CSV file."
    )

    ir_file = select_file(
        "Select Unknown IR CSV"
    )

    if ir_file is None:

        print(
            "❌ No IR file selected."
        )

        return None

    print()
    print(
        f"✅ Selected IR file:"
    )

    print(
        ir_file
    )

    # ========================================================
    # VERIFY IR
    # ========================================================

    try:

        detected_type = identify_csv_type(
            ir_file
        )

        if detected_type != "IR":

            print()
            print(
                "❌ Selected file is not an IR CSV."
            )

            print(
                f"Detected type: {detected_type}"
            )

            return None

    except Exception as error:

        print()
        print(
            "❌ Could not identify IR file."
        )

        print(
            f"Error: {error}"
        )

        return None

    # ========================================================
    # SELECT UNKNOWN FFT
    # ========================================================

    print()
    print(
        "Select the corresponding UNKNOWN FFT CSV file."
    )

    fft_file = select_file(
        "Select Unknown FFT CSV"
    )

    if fft_file is None:

        print(
            "❌ No FFT file selected."
        )

        return None

    print()
    print(
        f"✅ Selected FFT file:"
    )

    print(
        fft_file
    )

    # ========================================================
    # VERIFY FFT
    # ========================================================

    try:

        detected_type = identify_csv_type(
            fft_file
        )

        if detected_type != "FFT":

            print()
            print(
                "❌ Selected file is not an FFT CSV."
            )

            print(
                f"Detected type: {detected_type}"
            )

            return None

    except Exception as error:

        print()
        print(
            "❌ Could not identify FFT file."
        )

        print(
            f"Error: {error}"
        )

        return None

    # ========================================================
    # LOAD UNKNOWN IR
    # ========================================================

    print()
    print(
        "Loading unknown IR spectrum..."
    )

    try:

        unknown_ir_x, unknown_ir_y = load_csv(
            ir_file,
            expected_type="IR"
        )

        print(
            f"✅ Loaded {len(unknown_ir_y)} IR points."
        )

    except Exception as error:

        print()
        print(
            "❌ Failed to load unknown IR."
        )

        print(
            f"Error: {error}"
        )

        return None

    # ========================================================
    # LOAD UNKNOWN FFT
    # ========================================================

    print()
    print(
        "Loading unknown FFT spectrum..."
    )

    try:

        unknown_fft_x, unknown_fft_y = load_csv(
            fft_file,
            expected_type="FFT"
        )

        print(
            f"✅ Loaded {len(unknown_fft_y)} FFT points."
        )

    except Exception as error:

        print()
        print(
            "❌ Failed to load unknown FFT."
        )

        print(
            f"Error: {error}"
        )

        return None

    # ========================================================
    # READ DATABASE
    # ========================================================

    rows = read_database_spectra()

    if not rows:

        print()
        print(
            "❌ No compounds containing both "
            "IR and FFT data were found."
        )

        return None

    print()
    print(
        f"Database contains {len(rows)} compounds "
        "with IR and FFT data."
    )

    # ========================================================
    # COMPARE
    # ========================================================

    print()
    print(
        "Comparing unknown spectrum against database..."
    )

    print()

    results = []

    for (
        compound_name,
        ir_data,
        fft_data
    ) in rows:

        try:

            # =================================================
            # LOAD DATABASE IR
            # =================================================

            (
                database_ir_x,
                database_ir_y
            ) = load_database_csv(
                ir_data
            )

            # =================================================
            # LOAD DATABASE FFT
            # =================================================

            (
                database_fft_x,
                database_fft_y
            ) = load_database_csv(
                fft_data
            )

            # =================================================
            # CALCULATE MATCH
            # =================================================

            match = calculate_spectral_match(

                unknown_ir_x,
                unknown_ir_y,

                unknown_fft_x,
                unknown_fft_y,

                database_ir_x,
                database_ir_y,

                database_fft_x,
                database_fft_y
            )

            results.append(
                {
                    "compound": str(
                        compound_name
                    ),
                    **match
                }
            )

        except Exception as error:

            print(
                f"⚠️ Skipping "
                f"{compound_name}: "
                f"{error}"
            )

    # ========================================================
    # CHECK RESULTS
    # ========================================================

    if not results:

        print()
        print(
            "❌ No compounds could be compared."
        )

        return None

    # ========================================================
    # SORT BY COMBINED SCORE
    # ========================================================

    results.sort(
        key=lambda item: item[
            "combined_score"
        ],
        reverse=True
    )

    # ========================================================
    # TOP 10
    # ========================================================

    print()
    print("=" * 100)

    print(
        "                              TOP MATCHES"
    )

    print("=" * 100)

    print(
        f"{'Rank':<6}"
        f"{'Compound':<35}"
        f"{'IR':>10}"
        f"{'FFT':>10}"
        f"{'Combined':>12}"
    )

    print("-" * 100)

    for position, result in enumerate(
        results[:10],
        start=1
    ):

        print(
            f"{position:<6}"
            f"{result['compound'][:34]:<35}"
            f"{result['ir_score']:>9.2f}%"
            f"{result['fft_score']:>9.2f}%"
            f"{result['combined_score']:>11.2f}%"
        )

    print("=" * 100)

    # ========================================================
    # BEST MATCH
    # ========================================================

    best = results[0]

    best_name = best[
        "compound"
    ]

    print()
    print("=" * 80)

    print(
        "                           BEST MATCH"
    )

    print("=" * 80)

    print(
        f"Compound: {best_name}"
    )

    print()

    print(
        f"IR correlation: "
        f"{best['ir_correlation']:.6f}"
    )

    print(
        f"IR similarity: "
        f"{best['ir_score']:.2f}%"
    )

    print()

    print(
        f"FFT correlation: "
        f"{best['fft_correlation']:.6f}"
    )

    print(
        f"FFT similarity: "
        f"{best['fft_score']:.2f}%"
    )

    print()

    print(
        f"COMBINED SIMILARITY: "
        f"{best['combined_score']:.2f}%"
    )

    print()

    print(
        f"IR RMSE: "
        f"{best['ir_rmse']:.6f}"
    )

    print(
        f"FFT RMSE: "
        f"{best['fft_rmse']:.6f}"
    )

    print("=" * 80)

    # ========================================================
    # SAVE RESULTS
    # ========================================================

    try:

        results_dataframe = pd.DataFrame(
            results
        )

        results_dataframe.insert(
            0,
            "Rank",
            range(
                1,
                len(results_dataframe) + 1
            )
        )

        results_path = (
            RESULTS_FOLDER
            / "identification_results.csv"
        )

        results_dataframe.to_csv(
            results_path,
            index=False
        )

        print()
        print(
            "✅ Identification results saved:"
        )

        print(
            results_path
        )

    except Exception as error:

        print()
        print(
            f"⚠️ Could not save results: {error}"
        )

        results_path = None

    # ========================================================
    # PUBCHEM
    # ========================================================

    print()
    print(
        "Retrieving chemical information..."
    )

    pubchem_data = get_pubchem_information(
        best_name
    )

    # ========================================================
    # DISPLAY CHEMICAL INFORMATION
    # ========================================================

    if pubchem_data is not None:

        print()
        print("=" * 80)

        print(
            "                       CHEMICAL INFORMATION"
        )

        print("=" * 80)

        print(
            f"Compound: "
            f"{pubchem_data['title']}"
        )

        print(
            f"PubChem CID: "
            f"{pubchem_data['cid']}"
        )

        print(
            f"Molecular formula: "
            f"{pubchem_data['formula']}"
        )

        print(
            f"SMILES: "
            f"{pubchem_data['smiles']}"
        )

        print("=" * 80)

        print()
        print(
            "Opening chemical structure..."
        )

        show_structure(
            best_name,
            pubchem_data
        )

    else:

        print()
        print(
            "⚠️ Chemical structure could not "
            "be retrieved from PubChem."
        )

    # ========================================================
    # IMPORTANT INFORMATION
    # ========================================================

    print()
    print("=" * 80)

    print(
        "                              IMPORTANT"
    )

    print("=" * 80)

    print(
        "The similarity percentage measures how "
        "closely the unknown spectral shape matches "
        "the reference spectrum."
    )

    print()

    print(
        "It is NOT the probability that the "
        "compound identification is correct."
    )

    print()

    print(
        "IR weighting: 70%"
    )

    print(
        "FFT weighting: 30%"
    )

    print("=" * 80)

    # ========================================================
    # RETURN
    # ========================================================

    return {
        "best_match": best_name,
        "best_score": best[
            "combined_score"
        ],
        "results": results,
        "results_file": (
            str(results_path)
            if results_path is not None
            else None
        ),
        "pubchem": pubchem_data
    }


# ============================================================
# PUBLIC FUNCTION
# ============================================================

def Unknown_compound_identifier():

    """
    Public entry point.

    This name is retained so your existing programs that
    call Unknown_compound_identifier() will continue to work.
    """

    return identify_compound()


# ============================================================
# MAIN PROGRAM
# ============================================================

if __name__ == "__main__":

    Unknown_compound_identifier()