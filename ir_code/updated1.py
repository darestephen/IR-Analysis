#!/usr/bin/env python3

"""
===========================================================================
AUTOMATED IR MIXTURE IDENTIFICATION
===========================================================================

FFT-BASED IDENTIFICATION

Workflow:

1. User selects ONE TSV file using a Tkinter file picker.
2. TSV is loaded.
3. X/Y columns are detected automatically.
4. Absorbance is converted to transmittance when necessary.
5. Processed spectrum is saved.
6. Spectrum is resampled uniformly.
7. FFT is calculated using the project's FFT method.
8. ALL FFT magnitude values are retained.
9. Unknown FFT is compared against ALL database FFT magnitude values.
10. Top matches are displayed.
11. FFT diagrams are generated.
12. PubChem compound information is retrieved for top matches.
13. Mixture names are automatically separated into components.
14. 3D PubChem conformers are downloaded for components when available.
15. Interactive 3D molecular structures are embedded in HTML reports.
16. HTML and CSV reports are created.
17. HTML report and 3D viewer are automatically opened.

IMPORTANT:
No FFT magnitude threshold is used for comparison.

All FFT magnitude values are compared.
"""

import base64
import html
import json
import re
import sqlite3
import subprocess
import sys
import time
import urllib.parse
import urllib.request
import webbrowser
from io import StringIO
from pathlib import Path

import matplotlib

# Non-interactive backend for Chromebook/Linux
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import tkinter as tk
from tkinter import filedialog, messagebox


# ============================================================================
# CONFIGURATION
# ============================================================================

DATABASE_PATH = (
    Path.home()
    / "IR_Spectra"
    / "compound_spectra.db"
)

GENERATED_SPECTRA_DIR = (
    Path.home()
    / "IR_Spectra"
    / "generated_spectra"
)

RESULTS_DIR = (
    Path.home()
    / "IR_Spectra"
    / "identification_results"
)

STRUCTURES_DIR = (
    RESULTS_DIR
    / "top10_3d_structures"
)

FFT_PLOTS_DIR = (
    RESULTS_DIR
    / "fft_plots"
)

TOP_N = 10

# Common comparison length.
#
# IMPORTANT:
# ALL FFT values are retained.
#
# If FFT arrays have different lengths, the COMPLETE
# arrays are interpolated to this common comparison length.
COMPARISON_POINTS = 2048

PUBCHEM_TIMEOUT = 20

USER_AGENT = (
    "Automated-IR-Spectroscopy-Identifier/1.0 "
    "(educational/research application)"
)


# ============================================================================
# PRINTING
# ============================================================================

def print_header():
    print()
    print("=" * 75)
    print("AUTOMATED IR MIXTURE IDENTIFICATION")
    print("=" * 75)
    print()
    print("FFT comparison mode: ALL FFT MAGNITUDE VALUES")
    print("FFT magnitude threshold: NONE")
    print("Every FFT magnitude value will be compared.")
    print()


def print_subheader(title):
    print()
    print("-" * 75)
    print(title)
    print("-" * 75)


# ============================================================================
# TKINTER FILE SELECTION
# ============================================================================

def select_input_file():

    print_subheader("Select Unknown TSV Spectrum")

    print("Opening file-selection window...")

    root = tk.Tk()
    root.withdraw()

    try:
        root.attributes("-topmost", True)
    except Exception:
        pass

    file_path = filedialog.askopenfilename(
        title="Select IR Spectrum TSV File",
        filetypes=[
            ("TSV files", "*.tsv"),
            ("Tab-separated files", "*.txt"),
            ("All files", "*.*"),
        ],
    )

    root.destroy()

    if not file_path:
        print("No file selected.")
        sys.exit(0)

    path = Path(file_path)

    print()
    print("Selected file:")
    print(path)

    return path


# ============================================================================
# COLUMN DETECTION
# ============================================================================

def normalize_column_name(name):

    text = str(name).strip().lower()

    text = (
        text
        .replace("²", "2")
        .replace("−", "-")
        .replace("–", "-")
    )

    text = re.sub(
        r"[^a-z0-9]+",
        "_",
        text
    )

    return text.strip("_")


def detect_xy_columns(df):

    normalized = {
        column: normalize_column_name(column)
        for column in df.columns
    }

    x_candidates = [
        "x",
        "wavenumber",
        "wavenumber_cm_1",
        "wavenumber_cm1",
        "cm_1",
        "cm1",
        "frequency",
    ]

    y_candidates = [
        "y",
        "absorbance",
        "transmittance",
        "transmission",
        "intensity",
        "intensity_percent",
    ]

    x_column = None
    y_column = None

    # Exact matching
    for original, normalized_name in normalized.items():

        if normalized_name in x_candidates:
            x_column = original
            break

    for original, normalized_name in normalized.items():

        if normalized_name in y_candidates:
            y_column = original
            break

    # Name inspection
    if x_column is None:

        for original, normalized_name in normalized.items():

            if (
                "wave" in normalized_name
                or "wavenumber" in normalized_name
                or normalized_name.startswith("x")
            ):
                x_column = original
                break

    if y_column is None:

        for original, normalized_name in normalized.items():

            if (
                "trans" in normalized_name
                or "absorb" in normalized_name
                or "intens" in normalized_name
                or normalized_name.startswith("y")
            ):
                y_column = original
                break

    # Last resort
    if x_column is None or y_column is None:

        numeric_columns = []

        for column in df.columns:

            converted = pd.to_numeric(
                df[column],
                errors="coerce"
            )

            if converted.notna().sum() >= 2:
                numeric_columns.append(column)

        if len(numeric_columns) >= 2:

            if x_column is None:
                x_column = numeric_columns[0]

            if y_column is None:

                for column in numeric_columns:

                    if column != x_column:
                        y_column = column
                        break

    if x_column is None or y_column is None:

        raise ValueError(
            "Could not automatically detect X and Y columns.\n"
            f"Available columns: {list(df.columns)}"
        )

    return x_column, y_column


# ============================================================================
# LOAD TSV
# ============================================================================

def load_tsv_spectrum(file_path):

    print_subheader("Loading TSV Spectrum")

    print(f"Reading:\n{file_path}")

    try:

        df = pd.read_csv(
            file_path,
            sep="\t",
            engine="python"
        )

    except Exception as exc:

        raise RuntimeError(
            f"Could not read TSV file:\n{exc}"
        )

    if df.empty:
        raise ValueError(
            "The selected TSV file is empty."
        )

    print()
    print("Detected columns:")
    print(list(df.columns))

    x_column, y_column = detect_xy_columns(df)

    print()
    print(f"X column: {x_column}")
    print(f"Y column: {y_column}")

    x = pd.to_numeric(
        df[x_column],
        errors="coerce"
    )

    y = pd.to_numeric(
        df[y_column],
        errors="coerce"
    )

    valid = (
        x.notna()
        & y.notna()
        & np.isfinite(x)
        & np.isfinite(y)
    )

    x = x[valid].to_numpy(dtype=float)
    y = y[valid].to_numpy(dtype=float)

    if len(x) < 2:

        raise ValueError(
            "At least two valid X/Y points are required."
        )

    # ------------------------------------------------------------------------
    # Detect absorbance/transmittance
    # ------------------------------------------------------------------------

    y_column_normalized = normalize_column_name(
        y_column
    )

    is_absorbance = (
        "absorb" in y_column_normalized
    )

    is_transmittance = (
        "trans" in y_column_normalized
        or "transmission" in y_column_normalized
    )

    if is_absorbance:

        print()
        print("Y-axis detected as ABSORBANCE.")

        y = 10.0 ** (-y)

        print(
            "Converted absorbance values to "
            "transmittance using T = 10^(-A)."
        )

    elif is_transmittance:

        print()
        print("Y-axis detected as TRANSMITTANCE.")

    else:

        print()
        print(
            "Y-axis type could not be determined "
            "from the column name."
        )

        print(
            "Treating values as transmittance."
        )

    # ------------------------------------------------------------------------
    # Clean
    # ------------------------------------------------------------------------

    y = np.nan_to_num(
        y,
        nan=0.0,
        posinf=0.0,
        neginf=0.0
    )

    y = np.clip(
        y,
        1e-6,
        1.0
    )

    # ------------------------------------------------------------------------
    # Sort
    # ------------------------------------------------------------------------

    order = np.argsort(x)

    x = x[order]
    y = y[order]

    # Remove duplicate X values
    unique_x, unique_indices = np.unique(
        x,
        return_index=True
    )

    x = unique_x
    y = y[unique_indices]

    if len(x) < 2:

        raise ValueError(
            "Not enough unique X values after cleaning."
        )

    print()
    print(
        f"Valid spectrum points: {len(x)}"
    )

    print(
        f"Wavenumber range: "
        f"{x.min():.3f} - {x.max():.3f} cm⁻¹"
    )

    return x, y


# ============================================================================
# SAVE PROCESSED SPECTRUM
# ============================================================================

def save_processed_spectrum(
    x,
    y,
    input_path
):

    safe_name = re.sub(
        r"[^A-Za-z0-9._-]+",
        "_",
        input_path.stem
    )

    output_path = (
        GENERATED_SPECTRA_DIR
        / f"{safe_name}_processed.csv"
    )

    output_df = pd.DataFrame(
        {
            "Wavenumber (cm⁻¹)": x,
            "Transmittance": y,
        }
    )

    output_df.to_csv(
        output_path,
        index=False
    )

    print()
    print(
        "Processed spectrum saved to:"
    )
    print(output_path)

    return output_path


# ============================================================================
# UNIFORM RESAMPLING
# ============================================================================

def resample_uniformly(x, y):

    target_length = len(x)

    uniform_x = np.linspace(
        x.min(),
        x.max(),
        target_length
    )

    uniform_y = np.interp(
        uniform_x,
        x,
        y
    )

    return uniform_x, uniform_y


# ============================================================================
# FFT
# ============================================================================

def create_fft(intensity):

    """
    EXACT FFT METHOD USED BY THE PROJECT.

    DO NOT FILTER FFT MAGNITUDES.
    """

    intensity = np.asarray(
        intensity,
        dtype=float
    )

    if len(intensity) < 2:

        raise ValueError(
            "At least two points are required for FFT."
        )

    intensity = np.nan_to_num(
        intensity,
        nan=0.0,
        posinf=0.0,
        neginf=0.0
    )

    n = len(intensity)

    signal = (
        intensity
        - np.mean(intensity)
    )

    fft_values = np.fft.rfft(
        signal
    )

    fft_magnitude = (
        np.abs(fft_values)
        / n
    )

    if n % 2 == 0:

        if len(fft_magnitude) > 2:

            fft_magnitude[1:-1] *= 2.0

    else:

        if len(fft_magnitude) > 1:

            fft_magnitude[1:] *= 2.0

    fft_magnitude[0] = 0.0

    frequency_index = np.arange(
        len(fft_magnitude),
        dtype=int
    )

    return pd.DataFrame(
        {
            "Frequency Index": frequency_index,
            "FFT Magnitude": fft_magnitude,
        }
    )


# ============================================================================
# GENERATE UNKNOWN FFT
# ============================================================================

def generate_unknown_fft(x, y):

    print_subheader(
        "Generating Unknown FFT"
    )

    uniform_x, uniform_y = (
        resample_uniformly(
            x,
            y
        )
    )

    fft_df = create_fft(
        uniform_y
    )

    print(
        f"Uniform spectrum points: "
        f"{len(uniform_y)}"
    )

    print(
        f"FFT points generated: "
        f"{len(fft_df)}"
    )

    print(
        f"Maximum FFT magnitude: "
        f"{fft_df['FFT Magnitude'].max():.8f}"
    )

    print()
    print(
        "IMPORTANT: ALL FFT magnitude values "
        "will be compared."
    )

    return (
        fft_df,
        uniform_x,
        uniform_y
    )


# ============================================================================
# SAVE FFT
# ============================================================================

def save_fft(
    fft_df,
    input_path
):

    safe_name = re.sub(
        r"[^A-Za-z0-9._-]+",
        "_",
        input_path.stem
    )

    output_path = (
        GENERATED_SPECTRA_DIR
        / f"{safe_name}_fft.csv"
    )

    fft_df.to_csv(
        output_path,
        index=False
    )

    print()
    print(
        "Unknown FFT saved to:"
    )
    print(output_path)

    return output_path


# ============================================================================
# PLOT UNKNOWN FFT
# ============================================================================

def plot_unknown_fft(
    fft_df,
    output_path,
    sample_name
):

    print()
    print(
        "Creating unknown FFT diagram..."
    )

    x = fft_df[
        "Frequency Index"
    ].to_numpy(dtype=float)

    y = fft_df[
        "FFT Magnitude"
    ].to_numpy(dtype=float)

    plt.figure(
        figsize=(12, 6)
    )

    plt.plot(
        x,
        y,
        linewidth=1.0
    )

    plt.xlabel(
        "Frequency Index"
    )

    plt.ylabel(
        "FFT Magnitude"
    )

    plt.title(
        f"Unknown Spectrum FFT - {sample_name}"
    )

    plt.grid(
        True,
        alpha=0.3
    )

    plt.tight_layout()

    plt.savefig(
        output_path,
        dpi=180
    )

    plt.close()

    print(
        "Unknown FFT diagram saved to:"
    )
    print(output_path)

    return output_path


# ============================================================================
# DATABASE FFT PARSER
# ============================================================================

def parse_database_fft(fft_data):

    if fft_data is None:
        return None

    if isinstance(fft_data, bytes):

        fft_data = fft_data.decode(
            "utf-8",
            errors="ignore"
        )

    text = str(
        fft_data
    ).strip()

    if not text:
        return None

    # ------------------------------------------------------------------------
    # CSV
    # ------------------------------------------------------------------------

    try:

        df = pd.read_csv(
            StringIO(text)
        )

        columns = {
            normalize_column_name(column): column
            for column in df.columns
        }

        frequency_column = None
        magnitude_column = None

        for normalized, original in columns.items():

            if (
                "frequency" in normalized
                or normalized == "frequency_index"
            ):
                frequency_column = original

            if (
                "fft_magnitude" in normalized
                or (
                    "fft" in normalized
                    and "magnitude" in normalized
                )
                or normalized == "magnitude"
            ):
                magnitude_column = original

        if (
            frequency_column is not None
            and magnitude_column is not None
        ):

            frequency = pd.to_numeric(
                df[frequency_column],
                errors="coerce"
            )

            magnitude = pd.to_numeric(
                df[magnitude_column],
                errors="coerce"
            )

            valid = (
                frequency.notna()
                & magnitude.notna()
                & np.isfinite(frequency)
                & np.isfinite(magnitude)
            )

            frequency = (
                frequency[valid]
                .to_numpy(dtype=float)
            )

            magnitude = (
                magnitude[valid]
                .to_numpy(dtype=float)
            )

            if len(magnitude) >= 2:

                order = np.argsort(
                    frequency
                )

                frequency = frequency[
                    order
                ]

                magnitude = magnitude[
                    order
                ]

                return pd.DataFrame(
                    {
                        "Frequency Index": frequency,
                        "FFT Magnitude": magnitude,
                    }
                )

    except Exception:
        pass

    # ------------------------------------------------------------------------
    # JSON fallback
    # ------------------------------------------------------------------------

    try:

        parsed = json.loads(
            text
        )

        if isinstance(parsed, list):

            frequencies = []
            magnitudes = []

            for item in parsed:

                if not isinstance(
                    item,
                    dict
                ):
                    continue

                frequency = (
                    item.get(
                        "Frequency Index"
                    )
                    if item.get(
                        "Frequency Index"
                    ) is not None
                    else item.get(
                        "frequency"
                    )
                )

                if frequency is None:

                    frequency = item.get(
                        "frequency_index"
                    )

                magnitude = (
                    item.get(
                        "FFT Magnitude"
                    )
                    if item.get(
                        "FFT Magnitude"
                    ) is not None
                    else item.get(
                        "magnitude"
                    )
                )

                if magnitude is None:

                    magnitude = item.get(
                        "fft_magnitude"
                    )

                if (
                    frequency is not None
                    and magnitude is not None
                ):

                    frequencies.append(
                        float(frequency)
                    )

                    magnitudes.append(
                        float(magnitude)
                    )

            if len(magnitudes) >= 2:

                return pd.DataFrame(
                    {
                        "Frequency Index": frequencies,
                        "FFT Magnitude": magnitudes,
                    }
                )

    except Exception:
        pass

    return None


# ============================================================================
# LOAD DATABASE
# ============================================================================

def load_database_spectra():

    print_subheader(
        "Loading Database FFT Spectra"
    )

    if not DATABASE_PATH.exists():

        raise FileNotFoundError(
            f"Database not found:\n{DATABASE_PATH}"
        )

    connection = sqlite3.connect(
        DATABASE_PATH
    )

    try:

        query = """
        SELECT
            id,
            compound_name,
            fft_data
        FROM compound_spectra
        ORDER BY id
        """

        rows = connection.execute(
            query
        ).fetchall()

    finally:

        connection.close()

    print(
        f"Database records found: {len(rows)}"
    )

    database_spectra = []

    for row_id, compound_name, fft_data in rows:

        parsed_fft = parse_database_fft(
            fft_data
        )

        if parsed_fft is None:

            print(
                f"Skipping database ID {row_id}: "
                f"invalid FFT data."
            )

            continue

        if len(parsed_fft) < 2:

            print(
                f"Skipping database ID {row_id}: "
                f"not enough FFT values."
            )

            continue

        parsed_fft[
            "FFT Magnitude"
        ] = np.nan_to_num(
            parsed_fft[
                "FFT Magnitude"
            ].to_numpy(dtype=float),
            nan=0.0,
            posinf=0.0,
            neginf=0.0,
        )

        database_spectra.append(
            {
                "id": row_id,
                "compound_name": compound_name,
                "fft": parsed_fft,
            }
        )

    print(
        f"Usable database FFT spectra: "
        f"{len(database_spectra)}"
    )

    return database_spectra


# ============================================================================
# FFT RESAMPLING
# ============================================================================

def resample_fft(
    fft_df,
    comparison_points=COMPARISON_POINTS
):

    """
    Resample the COMPLETE FFT magnitude array.

    NO FFT magnitude values are filtered.
    """

    frequency = fft_df[
        "Frequency Index"
    ].to_numpy(dtype=float)

    magnitude = fft_df[
        "FFT Magnitude"
    ].to_numpy(dtype=float)

    if len(frequency) < 2:

        raise ValueError(
            "FFT must contain at least two points."
        )

    order = np.argsort(
        frequency
    )

    frequency = frequency[
        order
    ]

    magnitude = magnitude[
        order
    ]

    unique_frequency, unique_indices = (
        np.unique(
            frequency,
            return_index=True
        )
    )

    frequency = unique_frequency

    magnitude = magnitude[
        unique_indices
    ]

    normalized_frequency = np.linspace(
        0.0,
        1.0,
        len(frequency)
    )

    comparison_axis = np.linspace(
        0.0,
        1.0,
        comparison_points
    )

    resampled_magnitude = np.interp(
        comparison_axis,
        normalized_frequency,
        magnitude
    )

    return (
        comparison_axis,
        resampled_magnitude
    )


# ============================================================================
# FFT COMPARISON
# ============================================================================

def calculate_fft_value_similarity(
    unknown_values,
    database_values
):

    unknown_values = np.asarray(
        unknown_values,
        dtype=float
    )

    database_values = np.asarray(
        database_values,
        dtype=float
    )

    scale = max(
        float(
            np.max(
                np.abs(
                    unknown_values
                )
            )
        ),
        float(
            np.max(
                np.abs(
                    database_values
                )
            )
        ),
        np.finfo(float).eps
    )

    unknown_normalized = (
        unknown_values
        / scale
    )

    database_normalized = (
        database_values
        / scale
    )

    mean_absolute_difference = np.mean(
        np.abs(
            unknown_normalized
            - database_normalized
        )
    )

    similarity = (
        1.0
        - mean_absolute_difference
    ) * 100.0

    similarity = float(
        np.clip(
            similarity,
            0.0,
            100.0
        )
    )

    return similarity


def calculate_fft_correlation(
    unknown_values,
    database_values
):

    unknown_values = np.asarray(
        unknown_values,
        dtype=float
    )

    database_values = np.asarray(
        database_values,
        dtype=float
    )

    if len(unknown_values) != len(
        database_values
    ):

        raise ValueError(
            "FFT arrays must have equal length."
        )

    unknown_std = np.std(
        unknown_values
    )

    database_std = np.std(
        database_values
    )

    if (
        unknown_std <= 1e-15
        or database_std <= 1e-15
    ):

        return 0.0

    correlation_matrix = np.corrcoef(
        unknown_values,
        database_values
    )

    correlation = float(
        correlation_matrix[0, 1]
    )

    if not np.isfinite(
        correlation
    ):

        return 0.0

    correlation_score = (
        (correlation + 1.0)
        / 2.0
        * 100.0
    )

    return float(
        np.clip(
            correlation_score,
            0.0,
            100.0
        )
    )


def compare_fft_spectra(
    unknown_fft,
    database_fft
):

    (
        comparison_axis,
        unknown_values
    ) = resample_fft(
        unknown_fft
    )

    (
        _,
        database_values
    ) = resample_fft(
        database_fft
    )

    value_similarity = (
        calculate_fft_value_similarity(
            unknown_values,
            database_values
        )
    )

    correlation = (
        calculate_fft_correlation(
            unknown_values,
            database_values
        )
    )

    final_score = (
        value_similarity
        + correlation
    ) / 2.0

    return {
        "value_similarity": value_similarity,
        "correlation": correlation,
        "score": final_score,
        "comparison_axis": comparison_axis,
        "unknown_values": unknown_values,
        "database_values": database_values,
    }


# ============================================================================
# BEST MATCH FFT PLOT
# ============================================================================

def plot_best_match_fft(
    comparison,
    compound_name,
    output_path,
    sample_name
):

    print()
    print(
        "Creating best-match FFT comparison diagram..."
    )

    axis = comparison[
        "comparison_axis"
    ]

    unknown_values = comparison[
        "unknown_values"
    ]

    database_values = comparison[
        "database_values"
    ]

    unknown_max = np.max(
        np.abs(
            unknown_values
        )
    )

    database_max = np.max(
        np.abs(
            database_values
        )
    )

    if unknown_max > 0:

        unknown_plot = (
            unknown_values
            / unknown_max
        )

    else:

        unknown_plot = unknown_values

    if database_max > 0:

        database_plot = (
            database_values
            / database_max
        )

    else:

        database_plot = database_values

    plt.figure(
        figsize=(12, 6)
    )

    plt.plot(
        axis,
        unknown_plot,
        label="Unknown FFT",
        linewidth=1.2
    )

    plt.plot(
        axis,
        database_plot,
        label=f"Database: {compound_name}",
        linewidth=1.0,
        alpha=0.8
    )

    plt.xlabel(
        "Normalized FFT Frequency Index"
    )

    plt.ylabel(
        "Normalized FFT Magnitude"
    )

    plt.title(
        f"Unknown vs Best Database FFT\n"
        f"{sample_name} vs {compound_name}"
    )

    plt.legend()

    plt.grid(
        True,
        alpha=0.3
    )

    plt.tight_layout()

    plt.savefig(
        output_path,
        dpi=180
    )

    plt.close()

    print(
        "Best-match FFT diagram saved to:"
    )

    print(output_path)

    return output_path


# ============================================================================
# TOP 10 SCORE PLOT
# ============================================================================

def plot_top10_scores(
    results,
    output_path
):

    print()
    print(
        "Creating top-10 FFT score diagram..."
    )

    if not results:
        return None

    names = [
        str(
            result["compound_name"]
        )
        for result in results
    ]

    scores = [
        float(
            result["score"]
        )
        for result in results
    ]

    display_names = []

    for name in names:

        if len(name) > 30:
            name = name[:27] + "..."

        display_names.append(
            name
        )

    plt.figure(
        figsize=(12, 7)
    )

    positions = np.arange(
        len(display_names)
    )

    plt.barh(
        positions,
        scores
    )

    plt.yticks(
        positions,
        display_names
    )

    plt.xlabel(
        "FFT Match Score (%)"
    )

    plt.title(
        "Top 10 FFT Database Matches"
    )

    plt.xlim(
        0,
        100
    )

    plt.grid(
        axis="x",
        alpha=0.3
    )

    plt.tight_layout()

    plt.savefig(
        output_path,
        dpi=180
    )

    plt.close()

    print(
        "Top-10 FFT score diagram saved to:"
    )

    print(output_path)

    return output_path


# ============================================================================
# PUBCHEM REQUEST
# ============================================================================

def pubchem_request(
    url,
    timeout=PUBCHEM_TIMEOUT
):

    request = urllib.request.Request(
        url,
        headers={
            "User-Agent": USER_AGENT
        }
    )

    with urllib.request.urlopen(
        request,
        timeout=timeout
    ) as response:

        return response.read()


# ============================================================================
# PUBCHEM SINGLE COMPOUND LOOKUP
# ============================================================================

def pubchem_get_compound(
    compound_name
):

    encoded_name = urllib.parse.quote(
        compound_name,
        safe=""
    )

    url = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/"
        "compound/name/"
        f"{encoded_name}/property/"
        "Title,MolecularFormula,MolecularWeight,"
        "CanonicalSMILES,IsomericSMILES,"
        "XLogP,TPSA/JSON"
    )

    try:

        raw = pubchem_request(
            url
        )

        data = json.loads(
            raw.decode("utf-8")
        )

        properties = (
            data
            .get("PropertyTable", {})
            .get("Properties", [])
        )

        if not properties:
            return None

        return properties[0]

    except Exception:

        return None


# ============================================================================
# SPLIT MIXTURE NAME
# ============================================================================

def split_mixture_name(
    compound_name
):
    """
    Split common mixture/database naming formats.

    Examples:

        Benzoic Acid and Maleic anhydride
        Benzoic Acid + Maleic anhydride
        Benzoic Acid & Maleic anhydride

    The original name is NOT modified for FFT matching.
    This function is only used for PubChem/3D lookup.
    """

    name = str(
        compound_name
    ).strip()

    if not name:
        return []

    # First preserve the original exact name.
    separators = [
        r"\s+and\s+",
        r"\s+\+\s+",
        r"\s+&\s+",
        r"\s*;\s+",
    ]

    for separator in separators:

        parts = re.split(
            separator,
            name,
            flags=re.IGNORECASE
        )

        cleaned = [
            part.strip()
            for part in parts
            if part.strip()
        ]

        if len(cleaned) >= 2:

            return cleaned

    return [name]


# ============================================================================
# PUBCHEM LOOKUP INCLUDING MIXTURES
# ============================================================================

def get_pubchem_components(
    compound_name
):

    print()
    print(
        f"PubChem lookup for: {compound_name}"
    )

    # ------------------------------------------------------------------------
    # First try the complete name.
    # ------------------------------------------------------------------------

    exact = pubchem_get_compound(
        compound_name
    )

    if exact:

        return [
            {
                "name": compound_name,
                "properties": exact,
            }
        ]

    # ------------------------------------------------------------------------
    # Exact lookup failed.
    #
    # This is expected for mixture names such as:
    #
    # Benzoic Acid and Maleic anhydride
    #
    # Now split into components.
    # ------------------------------------------------------------------------

    components = split_mixture_name(
        compound_name
    )

    if len(components) <= 1:

        print(
            "No PubChem record found."
        )

        return []

    print()
    print(
        "The database name appears to represent "
        "a mixture."
    )

    print(
        "Searching PubChem for individual components:"
    )

    component_results = []

    for component in components:

        print(
            f"    -> {component}"
        )

        properties = pubchem_get_compound(
            component
        )

        if properties:

            component_results.append(
                {
                    "name": component,
                    "properties": properties,
                }
            )

        else:

            print(
                f"       PubChem record not found "
                f"for {component}"
            )

        time.sleep(0.25)

    return component_results


# ============================================================================
# 3D STRUCTURE DOWNLOAD
# ============================================================================

def save_3d_structure(
    compound_name,
    cid
):

    if cid is None:
        return None

    safe_name = re.sub(
        r"[^A-Za-z0-9._-]+",
        "_",
        str(compound_name)
    )

    output_path = (
        STRUCTURES_DIR
        / f"{safe_name}_{cid}_3d.sdf"
    )

    url = (
        "https://pubchem.ncbi.nlm.nih.gov/"
        "rest/pug/compound/cid/"
        f"{cid}/record/SDF?record_type=3d"
    )

    try:

        raw = pubchem_request(
            url
        )

        if not raw:
            return None

        output_path.write_bytes(
            raw
        )

        return output_path

    except Exception as exc:

        print(
            f"3D structure download failed "
            f"for {compound_name}: {exc}"
        )

        return None


# ============================================================================
# VALIDATE 3D STRUCTURE
# ============================================================================

def validate_3d_structure(
    sdf_path
):

    if sdf_path is None:
        return False

    try:

        text = Path(
            sdf_path
        ).read_text(
            errors="ignore"
        )

        if not text.strip():
            return False

        return (
            "V2000" in text
            or "V3000" in text
        )

    except Exception:
        return False


# ============================================================================
# ENRICH RESULT WITH PUBCHEM INFORMATION
# ============================================================================

def enrich_result_with_pubchem(
    result
):

    compound_name = result[
        "compound_name"
    ]

    component_results = (
        get_pubchem_components(
            compound_name
        )
    )

    result["pubchem_components"] = []

    # ------------------------------------------------------------------------
    # No PubChem information
    # ------------------------------------------------------------------------

    if not component_results:

        result[
            "pubchem_cid"
        ] = ""

        result[
            "molecular_formula"
        ] = ""

        result[
            "molecular_weight"
        ] = ""

        result[
            "canonical_smiles"
        ] = ""

        result[
            "structure_path"
        ] = ""

        return

    # ------------------------------------------------------------------------
    # Store component data
    # ------------------------------------------------------------------------

    for component in component_results:

        component_name = component[
            "name"
        ]

        properties = component[
            "properties"
        ]

        cid = properties.get(
            "CID",
            ""
        )

        formula = properties.get(
            "MolecularFormula",
            ""
        )

        molecular_weight = properties.get(
            "MolecularWeight",
            ""
        )

        smiles = properties.get(
            "ConnectivitySMILES",
            properties.get(
                "CanonicalSMILES",
                ""
            )
        )

        component_record = {
            "name": component_name,
            "cid": cid,
            "formula": formula,
            "molecular_weight": molecular_weight,
            "smiles": smiles,
            "structure_path": "",
        }

        # --------------------------------------------------------------------
        # Download 3D conformer
        # --------------------------------------------------------------------

        if cid:

            structure_path = (
                save_3d_structure(
                    component_name,
                    cid
                )
            )

            if structure_path:

                if validate_3d_structure(
                    structure_path
                ):

                    component_record[
                        "structure_path"
                    ] = str(
                        structure_path
                    )

                    print(
                        f"       3D conformer saved: "
                        f"{structure_path}"
                    )

                else:

                    print(
                        f"       Downloaded structure "
                        f"was not a valid SDF."
                    )

        result[
            "pubchem_components"
        ].append(
            component_record
        )

    # ------------------------------------------------------------------------
    # Backwards-compatible summary fields
    # ------------------------------------------------------------------------

    first = result[
        "pubchem_components"
    ][0]

    if len(
        result["pubchem_components"]
    ) == 1:

        result[
            "pubchem_cid"
        ] = first["cid"]

        result[
            "molecular_formula"
        ] = first["formula"]

        result[
            "molecular_weight"
        ] = first["molecular_weight"]

        result[
            "canonical_smiles"
        ] = first["smiles"]

        result[
            "structure_path"
        ] = first["structure_path"]

    else:

        # Mixture
        result[
            "pubchem_cid"
        ] = "; ".join(
            str(
                component["cid"]
            )
            for component in result[
                "pubchem_components"
            ]
            if component["cid"]
        )

        result[
            "molecular_formula"
        ] = "; ".join(
            str(
                component["formula"]
            )
            for component in result[
                "pubchem_components"
            ]
            if component["formula"]
        )

        result[
            "molecular_weight"
        ] = "; ".join(
            str(
                component["molecular_weight"]
            )
            for component in result[
                "pubchem_components"
            ]
            if component["molecular_weight"]
        )

        result[
            "canonical_smiles"
        ] = "; ".join(
            str(
                component["smiles"]
            )
            for component in result[
                "pubchem_components"
            ]
            if component["smiles"]
        )

        result[
            "structure_path"
        ] = "; ".join(
            str(
                component["structure_path"]
            )
            for component in result[
                "pubchem_components"
            ]
            if component["structure_path"]
        )


# ============================================================================
# IMAGE TO DATA URI
# ============================================================================

def image_to_data_uri(
    image_path
):

    image_path = Path(
        image_path
    )

    if not image_path.exists():
        return ""

    data = image_path.read_bytes()

    encoded = base64.b64encode(
        data
    ).decode("ascii")

    return (
        "data:image/png;base64,"
        + encoded
    )


# ============================================================================
# CSV REPORT
# ============================================================================

def save_results_csv(
    results,
    output_path
):

    rows = []

    for rank, result in enumerate(
        results,
        start=1
    ):

        components = result.get(
            "pubchem_components",
            []
        )

        component_names = "; ".join(
            str(
                component["name"]
            )
            for component in components
        )

        component_cids = "; ".join(
            str(
                component["cid"]
            )
            for component in components
            if component["cid"]
        )

        component_formulas = "; ".join(
            str(
                component["formula"]
            )
            for component in components
            if component["formula"]
        )

        structure_paths = "; ".join(
            str(
                component["structure_path"]
            )
            for component in components
            if component["structure_path"]
        )

        rows.append(
            {
                "Rank": rank,

                "Database ID": (
                    result["id"]
                ),

                "Compound": (
                    result["compound_name"]
                ),

                "FFT Value Similarity (%)": (
                    result["value_similarity"]
                ),

                "FFT Correlation (%)": (
                    result["correlation"]
                ),

                "FFT Final Score (%)": (
                    result["score"]
                ),

                "Unknown FFT Points": (
                    result["unknown_fft_points"]
                ),

                "Database FFT Points": (
                    result["database_fft_points"]
                ),

                "PubChem Components": (
                    component_names
                ),

                "PubChem CIDs": (
                    component_cids
                ),

                "Molecular Formulas": (
                    component_formulas
                ),

                "Molecular Weights": (
                    result.get(
                        "molecular_weight",
                        ""
                    )
                ),

                "Canonical SMILES": (
                    result.get(
                        "canonical_smiles",
                        ""
                    )
                ),

                "3D Structure Files": (
                    structure_paths
                ),
            }
        )

    df = pd.DataFrame(
        rows
    )

    df.to_csv(
        output_path,
        index=False
    )

    print()
    print(
        "Results CSV saved to:"
    )
    print(output_path)


# ============================================================================
# 3D VIEWER JAVASCRIPT
# ============================================================================

def make_3d_viewer_card(
    component,
    viewer_id
):

    structure_path = component.get(
        "structure_path"
    )

    if not structure_path:
        return ""

    structure_path = Path(
        structure_path
    )

    if not validate_3d_structure(
        structure_path
    ):
        return ""

    try:

        sdf_text = structure_path.read_text(
            errors="ignore"
        )

    except Exception:

        return ""

    sdf_json = json.dumps(
        sdf_text
    )

    name = html.escape(
        str(
            component.get(
                "name",
                "Unknown"
            )
        )
    )

    cid = html.escape(
        str(
            component.get(
                "cid",
                ""
            )
        )
    )

    formula = html.escape(
        str(
            component.get(
                "formula",
                ""
            )
        )
    )

    return f"""
    <div class="structure-card">

        <h3>{name}</h3>

        <div class="structure-info">

            <strong>PubChem CID:</strong>
            {cid}

            &nbsp;&nbsp;

            <strong>Formula:</strong>
            {formula}

        </div>

        <p class="small-note">
            PubChem computed 3D conformer.
        </p>

        <div
            id="{viewer_id}"
            class="viewer"
        ></div>

        <script>

        (function() {{

            const sdfData = {sdf_json};

            function loadViewer() {{

                const container =
                    document.getElementById(
                        "{viewer_id}"
                    );

                if (!container) {{
                    return;
                }}

                if (
                    typeof $3Dmol === "undefined"
                ) {{

                    container.innerHTML =
                        "<div class='viewer-error'>" +
                        "3Dmol.js could not be loaded. " +
                        "Please make sure this computer " +
                        "has an internet connection." +
                        "</div>";

                    return;
                }}

                try {{

                    const viewer =
                        $3Dmol.createViewer(
                            container,
                            {{
                                backgroundColor:
                                    "white"
                            }}
                        );

                    viewer.addModel(
                        sdfData,
                        "sdf"
                    );

                    viewer.setStyle(
                        {{}},
                        {{
                            stick: {{}}
                        }}
                    );

                    viewer.zoomTo();

                    viewer.render();

                    // Re-render after the browser finishes
                    // calculating the container dimensions.
                    setTimeout(
                        function() {{
                            viewer.resize();
                            viewer.zoomTo();
                            viewer.render();
                        }},
                        300
                    );

                }} catch (error) {{

                    container.innerHTML =
                        "<div class='viewer-error'>" +
                        "Could not render this 3D structure: " +
                        error.message +
                        "</div>";
                }}
            }}

            if (
                document.readyState === "loading"
            ) {{

                document.addEventListener(
                    "DOMContentLoaded",
                    loadViewer
                );

            }} else {{

                loadViewer();
            }}

        }})();

        </script>

    </div>
    """


# ============================================================================
# CREATE HTML REPORT
# ============================================================================

def create_html_report(
    sample_name,
    selected_file,
    results,
    unknown_fft_plot,
    best_match_plot,
    top10_plot,
    output_path,
    viewer_html_path
):

    print()
    print(
        "Creating HTML results report..."
    )

    unknown_fft_uri = image_to_data_uri(
        unknown_fft_plot
    )

    best_match_uri = image_to_data_uri(
        best_match_plot
    )

    top10_uri = image_to_data_uri(
        top10_plot
    )

    # ------------------------------------------------------------------------
    # Results table
    # ------------------------------------------------------------------------

    rows_html = []

    for rank, result in enumerate(
        results,
        start=1
    ):

        compound = html.escape(
            str(
                result["compound_name"]
            )
        )

        components = result.get(
            "pubchem_components",
            []
        )

        component_text = ""

        for component in components:

            component_text += (
                "<div>"
                "<strong>"
                + html.escape(
                    str(
                        component["name"]
                    )
                )
                + "</strong>"
                + "<br>"
                + "CID: "
                + html.escape(
                    str(
                        component["cid"]
                    )
                )
                + "<br>"
                + "Formula: "
                + html.escape(
                    str(
                        component["formula"]
                    )
                )
                + "</div><br>"
            )

        if not component_text:

            component_text = (
                "No PubChem information found."
            )

        rows_html.append(
            f"""
            <tr>

                <td>{rank}</td>

                <td>{result["id"]}</td>

                <td>
                    <strong>
                        {compound}
                    </strong>
                </td>

                <td>
                    {result["value_similarity"]:.4f}
                </td>

                <td>
                    {result["correlation"]:.4f}
                </td>

                <td class="score">
                    {result["score"]:.4f}
                </td>

                <td>
                    {result["unknown_fft_points"]}
                </td>

                <td>
                    {result["database_fft_points"]}
                </td>

                <td>
                    {component_text}
                </td>

            </tr>
            """
        )

    rows = "\n".join(
        rows_html
    )

    # ------------------------------------------------------------------------
    # 3D structures
    # ------------------------------------------------------------------------

    structure_sections = []

    structure_counter = 0

    for rank, result in enumerate(
        results,
        start=1
    ):

        components = result.get(
            "pubchem_components",
            []
        )

        valid_components = [
            component
            for component in components
            if component.get(
                "structure_path"
            )
            and validate_3d_structure(
                component.get(
                    "structure_path"
                )
            )
        ]

        if not valid_components:
            continue

        cards = []

        for component in valid_components:

            structure_counter += 1

            viewer_id = (
                f"report_viewer_"
                f"{structure_counter}"
            )

            card = make_3d_viewer_card(
                component,
                viewer_id
            )

            if card:
                cards.append(card)

        if cards:

            structure_sections.append(
                f"""
                <section class="match-structure">

                    <h3>
                        #{rank}
                        {html.escape(
                            str(
                                result["compound_name"]
                            )
                        )}
                    </h3>

                    {''.join(cards)}

                </section>
                """
            )

    structures_html = "".join(
        structure_sections
    )

    if not structures_html:

        structures_html = """
        <div class="no-structure">

            No validated PubChem 3D conformers
            were available for the top matches.

        </div>
        """

    generated_time = time.strftime(
        "%Y-%m-%d %H:%M:%S"
    )

    # Relative link
    viewer_filename = (
        Path(
            viewer_html_path
        ).name
    )

    html_content = f"""
<!DOCTYPE html>

<html lang="en">

<head>

<meta charset="UTF-8">

<meta
    name="viewport"
    content="width=device-width,
             initial-scale=1.0"
>

<title>
IR FFT Identification Results -
{html.escape(sample_name)}
</title>

<!--
3Dmol.js is loaded from the official 3Dmol
distribution. An internet connection is required
for the interactive 3D viewer.
-->

<script
src="https://3Dmol.csb.pitt.edu/build/3Dmol-min.js">
</script>

<style>

body {{

    font-family:
        Arial,
        Helvetica,
        sans-serif;

    margin: 0;

    padding: 30px;

    background: #f5f5f5;

    color: #222;

}}

.container {{

    max-width: 1500px;

    margin: auto;

    background: white;

    padding: 30px;

    border-radius: 12px;

    box-shadow:
        0 2px 12px
        rgba(0,0,0,0.08);

}}

h1 {{

    margin-top: 0;

}}

h2 {{

    margin-top: 45px;

    border-bottom:
        2px solid #ddd;

    padding-bottom: 8px;

}}

h3 {{

    margin-top: 0;

}}

.info {{

    background: #eeeeee;

    padding: 15px;

    border-radius: 8px;

    margin-bottom: 25px;

}}

.notice {{

    background: #eef6ff;

    border-left:
        5px solid #4488cc;

    padding: 15px;

    margin:
        20px 0;

}}

.open-viewer {{

    display: inline-block;

    padding:
        12px 18px;

    background: #333;

    color: white;

    text-decoration: none;

    border-radius: 6px;

    margin:
        10px 0 20px 0;

}}

.open-viewer:hover {{

    background: #555;

}}

img {{

    max-width: 100%;

    height: auto;

    border:
        1px solid #cccccc;

    border-radius: 8px;

    margin-top: 10px;

}}

.table-container {{

    overflow-x: auto;

}}

table {{

    border-collapse: collapse;

    width: 100%;

    margin-top: 20px;

}}

th,
td {{

    border:
        1px solid #cccccc;

    padding: 9px;

    text-align: left;

    vertical-align: top;

}}

th {{

    background: #eeeeee;

}}

tr:nth-child(even) {{

    background: #f8f8f8;

}}

.score {{

    font-weight: bold;

}}

.match-structure {{

    margin-top: 30px;

    padding: 20px;

    border:
        1px solid #ddd;

    border-radius: 10px;

    background: #fafafa;

}}

.structure-card {{

    background: white;

    border:
        1px solid #ddd;

    border-radius: 10px;

    padding: 15px;

    margin:
        15px 0;

}}

.structure-info {{

    margin-bottom: 5px;

}}

.small-note {{

    color: #666;

    font-size: 13px;

}}

.viewer {{

    width: 100%;

    height: 500px;

    position: relative;

    border:
        1px solid #ccc;

    border-radius: 8px;

    background: white;

}}

.viewer-error {{

    padding: 30px;

    color: #a00;

    text-align: center;

}}

.no-structure {{

    padding: 20px;

    background: #fff3f3;

    border-radius: 8px;

}}

.footer {{

    margin-top: 40px;

    color: #777;

    font-size: 13px;

}}

</style>

</head>

<body>

<div class="container">

<h1>
Automated IR FFT Identification Results
</h1>

<div class="info">

<p>
<strong>Unknown sample:</strong>
{html.escape(sample_name)}
</p>

<p>
<strong>Input file:</strong>
{html.escape(str(selected_file))}
</p>

<p>
<strong>Comparison method:</strong>
ALL FFT magnitude values
</p>

<p>
<strong>FFT magnitude threshold:</strong>
None
</p>

<p>
<strong>Comparison points:</strong>
{COMPARISON_POINTS}
</p>

<p>
<strong>Generated:</strong>
{generated_time}
</p>

</div>

<div class="notice">

<strong>
Important:
</strong>

The FFT identification score uses
ALL FFT magnitude values.

No FFT magnitude threshold was applied.

The 100% match is therefore preserved
exactly according to the current FFT
comparison method.

</div>

<a
    class="open-viewer"
    href="{html.escape(viewer_filename)}"
    target="_blank"
>
Open separate interactive 3D viewer
</a>


<h2>
Unknown FFT Diagram
</h2>

<p>
This diagram contains the complete FFT
generated from the unknown spectrum.
No FFT magnitude values were removed
by a threshold.
</p>

<img
    src="{unknown_fft_uri}"
    alt="Unknown FFT diagram"
>


<h2>
Top 10 FFT Matches
</h2>

<img
    src="{top10_uri}"
    alt="Top 10 FFT scores"
>


<h2>
Best Match FFT Comparison
</h2>

<img
    src="{best_match_uri}"
    alt="Best match FFT comparison"
>


<h2>
Results
</h2>

<div class="table-container">

<table>

<thead>

<tr>

<th>Rank</th>

<th>DB ID</th>

<th>Compound</th>

<th>FFT Value Similarity (%)</th>

<th>FFT Correlation (%)</th>

<th>Final FFT Score (%)</th>

<th>Unknown FFT Points</th>

<th>Database FFT Points</th>

<th>PubChem Information</th>

</tr>

</thead>

<tbody>

{rows}

</tbody>

</table>

</div>


<h2>
Interactive 3D Molecular Structures
</h2>

<div class="notice">

<strong>
3D structure source:
</strong>

PubChem computed 3D conformers.

For mixture matches, each component
is displayed separately.

The FFT match itself is based on the
database spectrum and is not changed by
the PubChem lookup.

</div>

{structures_html}


<div class="footer">

Generated by Automated IR Mixture Identification.

</div>

</div>

</body>

</html>
"""

    output_path.write_text(
        html_content,
        encoding="utf-8"
    )

    print()
    print(
        "HTML report saved to:"
    )

    print(output_path)

    return output_path


# ============================================================================
# CREATE SEPARATE 3D VIEWER
# ============================================================================

def create_3d_viewer_html(
    results,
    output_path
):

    print()
    print(
        "Creating 3D structure viewer..."
    )

    cards = []

    viewer_counter = 0

    for rank, result in enumerate(
        results,
        start=1
    ):

        components = result.get(
            "pubchem_components",
            []
        )

        valid_components = [
            component
            for component in components
            if component.get(
                "structure_path"
            )
            and validate_3d_structure(
                component.get(
                    "structure_path"
                )
            )
        ]

        if not valid_components:
            continue

        cards.append(
            f"""
            <div class="match-card">

                <h2>
                    #{rank}
                    {html.escape(
                        str(
                            result["compound_name"]
                        )
                    )}
                </h2>

                <p>
                    FFT Score:
                    <strong>
                        {result["score"]:.4f}%
                    </strong>
                </p>
            """
        )

        for component in valid_components:

            viewer_counter += 1

            viewer_id = (
                f"standalone_viewer_"
                f"{viewer_counter}"
            )

            cards.append(
                make_3d_viewer_card(
                    component,
                    viewer_id
                )
            )

        cards.append(
            """
            </div>
            """
        )

    if not cards:

        cards.append(
            """
            <div class="no-structure">

                No validated 3D structures were
                available for the top matches.

            </div>
            """
        )

    html_content = f"""
<!DOCTYPE html>

<html lang="en">

<head>

<meta charset="UTF-8">

<meta
    name="viewport"
    content="width=device-width,
             initial-scale=1.0"
>

<title>
IR Compound 3D Structures
</title>

<script
src="https://3Dmol.csb.pitt.edu/build/3Dmol-min.js">
</script>

<style>

body {{

    font-family:
        Arial,
        Helvetica,
        sans-serif;

    margin: 0;

    padding: 30px;

    background: #f5f5f5;

}}

.container {{

    max-width: 1400px;

    margin: auto;

}}

h1 {{

    background: white;

    padding: 20px;

    border-radius: 10px;

}}

.match-card {{

    background: white;

    padding: 20px;

    margin-bottom: 30px;

    border-radius: 10px;

    box-shadow:
        0 2px 10px
        rgba(0,0,0,0.07);

}}

.structure-card {{

    background: white;

    border:
        1px solid #ddd;

    border-radius: 10px;

    padding: 15px;

    margin:
        20px 0;

}}

.structure-info {{

    margin-bottom: 5px;

}}

.small-note {{

    color: #666;

    font-size: 13px;

}}

.viewer {{

    width: 100%;

    height: 550px;

    position: relative;

    border:
        1px solid #ccc;

    border-radius: 8px;

    background: white;

}}

.viewer-error {{

    padding: 30px;

    color: #a00;

    text-align: center;

}}

.no-structure {{

    background: white;

    padding: 30px;

    border-radius: 10px;

}}

</style>

</head>

<body>

<div class="container">

<h1>
Top IR Compound 3D Structures
</h1>

<p>
Interactive structures are PubChem computed
3D conformers. For mixture matches, the
components are displayed separately.
</p>

{''.join(cards)}

</div>

</body>

</html>
"""

    output_path.write_text(
        html_content,
        encoding="utf-8"
    )

    print()
    print(
        "3D viewer HTML saved to:"
    )

    print(output_path)

    return output_path


# ============================================================================
# OPEN HTML IN BROWSER
# ============================================================================

def open_html_file(
    file_path,
    description
):

    file_path = Path(
        file_path
    ).resolve()

    print()
    print(
        f"Opening {description}:"
    )

    print(file_path)

    try:

        file_url = file_path.as_uri()

        opened = webbrowser.open_new_tab(
            file_url
        )

        if opened:

            print(
                f"{description} opened in "
                "the default browser."
            )

            return True

    except Exception as exc:

        print(
            f"webbrowser failed: {exc}"
        )

    # ------------------------------------------------------------------------
    # Linux fallback
    # ------------------------------------------------------------------------

    try:

        subprocess.Popen(
            [
                "xdg-open",
                str(file_path)
            ],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL
        )

        print(
            f"{description} opened using xdg-open."
        )

        return True

    except Exception as exc:

        print(
            f"Could not automatically open "
            f"{description}: {exc}"
        )

    return False


# ============================================================================
# MAIN
# ============================================================================

def main():

    print_header()

    # ------------------------------------------------------------------------
    # Create directories
    # ------------------------------------------------------------------------

    GENERATED_SPECTRA_DIR.mkdir(
        parents=True,
        exist_ok=True
    )

    RESULTS_DIR.mkdir(
        parents=True,
        exist_ok=True
    )

    STRUCTURES_DIR.mkdir(
        parents=True,
        exist_ok=True
    )

    FFT_PLOTS_DIR.mkdir(
        parents=True,
        exist_ok=True
    )

    # ------------------------------------------------------------------------
    # Select input
    # ------------------------------------------------------------------------

    input_path = select_input_file()

    sample_name = input_path.stem

    safe_name = re.sub(
        r"[^A-Za-z0-9._-]+",
        "_",
        sample_name
    )

    # ------------------------------------------------------------------------
    # Load spectrum
    # ------------------------------------------------------------------------

    x, y = load_tsv_spectrum(
        input_path
    )

    # ------------------------------------------------------------------------
    # Save processed spectrum
    # ------------------------------------------------------------------------

    save_processed_spectrum(
        x,
        y,
        input_path
    )

    # ------------------------------------------------------------------------
    # Generate FFT
    # ------------------------------------------------------------------------

    (
        unknown_fft,
        uniform_x,
        uniform_y
    ) = generate_unknown_fft(
        x,
        y
    )

    # ------------------------------------------------------------------------
    # Save FFT
    # ------------------------------------------------------------------------

    unknown_fft_path = save_fft(
        unknown_fft,
        input_path
    )

    # ------------------------------------------------------------------------
    # Unknown FFT plot
    # ------------------------------------------------------------------------

    unknown_fft_plot = (
        FFT_PLOTS_DIR
        / f"{safe_name}_unknown_fft.png"
    )

    plot_unknown_fft(
        unknown_fft,
        unknown_fft_plot,
        sample_name
    )

    # ------------------------------------------------------------------------
    # Load database
    # ------------------------------------------------------------------------

    database_spectra = (
        load_database_spectra()
    )

    if not database_spectra:

        raise RuntimeError(
            "No usable database FFT spectra were found."
        )

    # ------------------------------------------------------------------------
    # Compare
    # ------------------------------------------------------------------------

    print_subheader(
        "Comparing Unknown Against Database"
    )

    print(
        "Comparison method:"
    )

    print(
        "ALL FFT magnitude values"
    )

    print()

    results = []

    for index, database_record in enumerate(
        database_spectra,
        start=1
    ):

        compound_name = (
            database_record[
                "compound_name"
            ]
        )

        database_fft = (
            database_record[
                "fft"
            ]
        )

        try:

            comparison = compare_fft_spectra(
                unknown_fft,
                database_fft
            )

            result = {

                "id":
                    database_record[
                        "id"
                    ],

                "compound_name":
                    compound_name,

                "value_similarity":
                    comparison[
                        "value_similarity"
                    ],

                "correlation":
                    comparison[
                        "correlation"
                    ],

                "score":
                    comparison[
                        "score"
                    ],

                "unknown_fft_points":
                    len(
                        unknown_fft
                    ),

                "database_fft_points":
                    len(
                        database_fft
                    ),

                "_comparison":
                    comparison,
            }

            results.append(
                result
            )

        except Exception as exc:

            print(
                f"Comparison failed for "
                f"{compound_name}: {exc}"
            )

        if index % 10 == 0:

            print(
                f"Compared {index} / "
                f"{len(database_spectra)} "
                f"database spectra..."
            )

    if not results:

        raise RuntimeError(
            "No database spectra could be compared."
        )

    # ------------------------------------------------------------------------
    # Sort
    # ------------------------------------------------------------------------

    results.sort(
        key=lambda item: item["score"],
        reverse=True
    )

    top_results = results[
        :TOP_N
    ]

    # ------------------------------------------------------------------------
    # Display results
    # ------------------------------------------------------------------------

    print_subheader(
        "TOP FFT MATCHES"
    )

    for rank, result in enumerate(
        top_results,
        start=1
    ):

        print(
            f"{rank:2d}. "
            f"{result['compound_name']}"
        )

        print(
            f"    FFT Value Similarity: "
            f"{result['value_similarity']:.2f}%"
        )

        print(
            f"    FFT Correlation: "
            f"{result['correlation']:.2f}%"
        )

        print(
            f"    FINAL FFT SCORE: "
            f"{result['score']:.2f}%"
        )

        print(
            f"    Unknown FFT points: "
            f"{result['unknown_fft_points']}"
        )

        print(
            f"    Database FFT points: "
            f"{result['database_fft_points']}"
        )

        print()

    # ------------------------------------------------------------------------
    # PubChem
    # ------------------------------------------------------------------------

    print_subheader(
        "PubChem Information"
    )

    for rank, result in enumerate(
        top_results,
        start=1
    ):

        print(
            f"[{rank}/{len(top_results)}] "
            f"{result['compound_name']}"
        )

        enrich_result_with_pubchem(
            result
        )

        components = result.get(
            "pubchem_components",
            []
        )

        if components:

            for component in components:

                print(
                    f"    Component: "
                    f"{component['name']}"
                )

                print(
                    f"        CID: "
                    f"{component['cid']}"
                )

                print(
                    f"        Formula: "
                    f"{component['formula']}"
                )

                if component.get(
                    "structure_path"
                ):

                    print(
                        "        3D conformer: "
                        "AVAILABLE"
                    )

                else:

                    print(
                        "        3D conformer: "
                        "NOT AVAILABLE"
                    )

        else:

            print(
                "    No PubChem information found."
            )

        time.sleep(
            0.25
        )

    # ------------------------------------------------------------------------
    # Best match FFT plot
    # ------------------------------------------------------------------------

    best_result = top_results[0]

    best_match_plot = (
        FFT_PLOTS_DIR
        / f"{safe_name}_best_match_fft.png"
    )

    plot_best_match_fft(
        best_result[
            "_comparison"
        ],
        best_result[
            "compound_name"
        ],
        best_match_plot,
        sample_name
    )

    # ------------------------------------------------------------------------
    # Top 10 plot
    # ------------------------------------------------------------------------

    top10_plot = (
        FFT_PLOTS_DIR
        / f"{safe_name}_top10_fft_scores.png"
    )

    plot_top10_scores(
        top_results,
        top10_plot
    )

    # ------------------------------------------------------------------------
    # Results CSV
    # ------------------------------------------------------------------------

    results_csv = (
        RESULTS_DIR
        / f"{safe_name}_fft_identification_results.csv"
    )

    save_results_csv(
        top_results,
        results_csv
    )

    # ------------------------------------------------------------------------
    # 3D viewer
    #
    # Create this BEFORE HTML report so the report can link to it.
    # ------------------------------------------------------------------------

    viewer_html = (
        RESULTS_DIR
        / f"{safe_name}_top10_3d_viewer.html"
    )

    create_3d_viewer_html(
        top_results,
        viewer_html
    )

    # ------------------------------------------------------------------------
    # HTML report
    # ------------------------------------------------------------------------

    html_report = (
        RESULTS_DIR
        / f"{safe_name}_fft_identification_report.html"
    )

    create_html_report(
        sample_name,
        input_path,
        top_results,
        unknown_fft_plot,
        best_match_plot,
        top10_plot,
        html_report,
        viewer_html
    )

    # ------------------------------------------------------------------------
    # Remove internal comparison data AFTER reports are created.
    # ------------------------------------------------------------------------

    for result in results:

        result.pop(
            "_comparison",
            None
        )

    # ------------------------------------------------------------------------
    # Final output
    # ------------------------------------------------------------------------

    print()
    print("=" * 75)
    print("IDENTIFICATION COMPLETE")
    print("=" * 75)

    print()

    print(
        f"Best FFT match: "
        f"{best_result['compound_name']}"
    )

    print(
        f"FFT Value Similarity: "
        f"{best_result['value_similarity']:.2f}%"
    )

    print(
        f"FFT Correlation: "
        f"{best_result['correlation']:.2f}%"
    )

    print(
        f"Final FFT Score: "
        f"{best_result['score']:.2f}%"
    )

    print()

    print(
        "ALL FFT magnitude values were used "
        "in the comparison."
    )

    print()
    print(
        "Unknown FFT:"
    )
    print(unknown_fft_path)

    print()
    print(
        "Unknown FFT diagram:"
    )
    print(unknown_fft_plot)

    print()
    print(
        "Best-match FFT diagram:"
    )
    print(best_match_plot)

    print()
    print(
        "Top-10 FFT score diagram:"
    )
    print(top10_plot)

    print()
    print(
        "Results CSV:"
    )
    print(results_csv)

    print()
    print(
        "HTML report:"
    )
    print(html_report)

    print()
    print(
        "3D structure viewer:"
    )
    print(viewer_html)

    # ------------------------------------------------------------------------
    # Automatically open both HTML pages.
    # ------------------------------------------------------------------------

    print()
    print(
        "Opening results in your browser..."
    )

    open_html_file(
        html_report,
        "HTML report"
    )

    time.sleep(
        1.0
    )

    open_html_file(
        viewer_html,
        "3D structure viewer"
    )

    print()
    print("=" * 75)


# ============================================================================
# RUN
# ============================================================================

if __name__ == "__main__":

    try:

        main()

    except KeyboardInterrupt:

        print()
        print(
            "Program cancelled by user."
        )

        sys.exit(1)

    except Exception as exc:

        print()
        print("=" * 75)
        print("ERROR")
        print("=" * 75)
        print()
        print(str(exc))
        print()

        try:

            root = tk.Tk()
            root.withdraw()

            messagebox.showerror(
                "IR Identification Error",
                str(exc)
            )

            root.destroy()

        except Exception:
            pass

        sys.exit(1)