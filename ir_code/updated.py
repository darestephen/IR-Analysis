```python
"""
=======================================================================
AUTOMATED IDENTIFICATION OF UNKNOWN COMPOUNDS
VIA INFRARED SPECTROSCOPY AND DATABASE MATCHING

INPUT:
    One TSV file containing X and Y spectral data.

AUTOMATIC WORKFLOW:
    TSV
      |
      +----> CSV
      |
      +----> Fourier Transform
                  |
                  +----> FFT CSV
      |
      +----> IR comparison
      |
      +----> Diagnostic-region comparison
      |
      +----> FFT comparison
      |
      +----> Nitrile analysis
      |
      +----> Mixture identification
      |
      +----> Residual analysis
      |
      +----> PubChem information
      |
      +----> Results CSV

SCORING:
    IR                60%
    Diagnostic        25%
    FFT               15%
    Nitrile bonus     Additional weighted contribution

DATABASE:
    ~/IR_Spectra/compound_spectra.db

TABLE:
    compound_spectra

The program attempts to work with several possible database
column names so that it can accommodate your existing database.
=======================================================================
"""

# =====================================================================
# IMPORTS
# =====================================================================

from pathlib import Path
import sqlite3
import re
import math
import io
import urllib.parse

import numpy as np
import pandas as pd
import requests

import tkinter as tk
from tkinter import filedialog, messagebox

# Optional PubChem image support
try:
    from PIL import Image, ImageTk
    PIL_AVAILABLE = True
except ImportError:
    PIL_AVAILABLE = False


# =====================================================================
# CONFIGURATION
# =====================================================================

DATABASE_PATH = (
    Path.home() /
    "IR_Spectra" /
    "compound_spectra.db"
)

RESULTS_DIRECTORY = (
    Path.home() /
    "IR_Spectra" /
    "identification_results"
)

GENERATED_DIRECTORY = (
    Path.home() /
    "IR_Spectra" /
    "generated_spectra"
)

MINIMUM_POINTS = 10

DEFAULT_MAX_COMPOUNDS = 5

DEFAULT_MINIMUM_SIMILARITY = 60.0

COMPARISON_LENGTH = 4096

FFT_COMPARISON_LENGTH = 2048

DEFAULT_MIN_WAVENUMBER = 400.0

DEFAULT_MAX_WAVENUMBER = 4000.0

EPSILON = 1e-12

REQUEST_TIMEOUT = 20

# Existing weighting system
IR_WEIGHT = 0.60
DIAGNOSTIC_WEIGHT = 0.25
FFT_WEIGHT = 0.15

# Nitrile settings
NITRILE_CENTER = 2250.0
NITRILE_HALF_WIDTH = 55.0
NITRILE_WEIGHT = 4.0


# =====================================================================
# DIAGNOSTIC REGIONS
# =====================================================================

DIAGNOSTIC_REGIONS = [

    {
        "name": "O-H / N-H",
        "min": 3200.0,
        "max": 3650.0,
        "weight": 1.0
    },

    {
        "name": "sp2 C-H",
        "min": 3000.0,
        "max": 3100.0,
        "weight": 1.0
    },

    {
        "name": "sp3 C-H",
        "min": 2800.0,
        "max": 3000.0,
        "weight": 1.0
    },

    {
        "name": "C#N nitrile",
        "min": 2185.0,
        "max": 2305.0,
        "weight": 4.0
    },

    {
        "name": "C=O / C=C",
        "min": 1650.0,
        "max": 1800.0,
        "weight": 1.5
    },

    {
        "name": "C=C / bending",
        "min": 1450.0,
        "max": 1600.0,
        "weight": 1.0
    },

    {
        "name": "C-O / fingerprint",
        "min": 1000.0,
        "max": 1300.0,
        "weight": 1.2
    },

    {
        "name": "Fingerprint",
        "min": 650.0,
        "max": 1000.0,
        "weight": 1.2
    }
]


# =====================================================================
# TKINTER HELPERS
# =====================================================================

def select_unknown_tsv():
    """
    Ask the user to select ONE TSV file.
    """

    root = tk.Tk()
    root.withdraw()

    file_path = filedialog.askopenfilename(
        title="Select Unknown IR TSV File",
        filetypes=[
            ("TSV files", "*.tsv"),
            ("Tab-separated files", "*.txt"),
            ("All files", "*.*")
        ]
    )

    root.destroy()

    if not file_path:
        return None

    return Path(file_path)


# =====================================================================
# NUMERIC CLEANING
# =====================================================================

def clean_numeric(value):
    """
    Convert a value to float.

    Handles:
        commas
        spaces
        percentages
        common non-numeric characters
    """

    if pd.isna(value):
        return np.nan

    if isinstance(value, (int, float, np.integer, np.floating)):
        return float(value)

    text = str(value).strip()

    if not text:
        return np.nan

    text = text.replace(",", "")

    text = text.replace("%", "")

    # Keep numbers, decimal, sign, and scientific notation
    match = re.search(
        r"[-+]?(?:\d*\.\d+|\d+\.?)(?:[eE][-+]?\d+)?",
        text
    )

    if not match:
        return np.nan

    try:
        return float(match.group())
    except Exception:
        return np.nan


# =====================================================================
# COLUMN DETECTION
# =====================================================================

def detect_xy_columns(df):
    """
    Automatically detect X and Y columns.

    Recognized X names include:
        X
        Wavenumber
        Wavenumber (cm-1)
        Frequency
        cm-1

    Recognized Y names include:
        Y
        Intensity
        Absorbance
        Transmittance
        %T
        Signal

    If names are not recognized, the first two numeric columns
    are used.
    """

    columns = list(df.columns)

    normalized = {}

    for col in columns:
        key = str(col).strip().lower()

        key = (
            key
            .replace("−", "-")
            .replace(" ", "")
        )

        normalized[key] = col

    x_names = [
        "x",
        "wavenumber",
        "wavenumber(cm-1)",
        "wavenumber(cm−1)",
        "cm-1",
        "cm−1",
        "frequency"
    ]

    y_names = [
        "y",
        "intensity",
        "absorbance",
        "abs",
        "transmittance",
        "%transmittance",
        "%t",
        "signal"
    ]

    x_column = None
    y_column = None

    for name in x_names:
        key = (
            name
            .strip()
            .lower()
            .replace("−", "-")
            .replace(" ", "")
        )

        if key in normalized:
            x_column = normalized[key]
            break

    for name in y_names:
        key = (
            name
            .strip()
            .lower()
            .replace("−", "-")
            .replace(" ", "")
        )

        if key in normalized:
            y_column = normalized[key]
            break

    # -------------------------------------------------------------
    # Fallback to numeric columns
    # -------------------------------------------------------------

    numeric_columns = []

    for col in columns:

        numeric = pd.to_numeric(
            df[col],
            errors="coerce"
        )

        if numeric.notna().sum() >= MINIMUM_POINTS:
            numeric_columns.append(col)

    if x_column is None:

        if numeric_columns:
            x_column = numeric_columns[0]

    if y_column is None:

        for col in numeric_columns:

            if col != x_column:
                y_column = col
                break

    if x_column is None or y_column is None:

        raise ValueError(
            "Could not determine the X and Y columns.\n\n"
            "The TSV file must contain at least two numeric columns."
        )

    return x_column, y_column


# =====================================================================
# LOAD TSV
# =====================================================================

def load_tsv_spectrum(tsv_path):
    """
    Load X and Y data from TSV.
    """

    # First attempt: true TSV
    try:

        df = pd.read_csv(
            tsv_path,
            sep="\t",
            engine="python"
        )

    except Exception:

        # Fallback for whitespace-separated data
        df = pd.read_csv(
            tsv_path,
            sep=r"\s+",
            engine="python"
        )

    if df.empty:
        raise ValueError(
            "The TSV file contains no data."
        )

    x_column, y_column = detect_xy_columns(df)

    print()
    print("Detected columns:")
    print(f"  X = {x_column}")
    print(f"  Y = {y_column}")

    x = df[x_column].apply(
        clean_numeric
    )

    y = df[y_column].apply(
        clean_numeric
    )

    valid = (
        x.notna() &
        y.notna()
    )

    x = x[valid].to_numpy(
        dtype=float
    )

    y = y[valid].to_numpy(
        dtype=float
    )

    if len(x) < MINIMUM_POINTS:

        raise ValueError(
            f"Only {len(x)} valid points were found."
        )

    # -------------------------------------------------------------
    # Remove duplicate X values
    # -------------------------------------------------------------

    temp = pd.DataFrame({
        "x": x,
        "y": y
    })

    temp = (
        temp
        .groupby("x", as_index=False)["y"]
        .mean()
    )

    x = temp["x"].to_numpy(
        dtype=float
    )

    y = temp["y"].to_numpy(
        dtype=float
    )

    # -------------------------------------------------------------
    # Sort by X
    # -------------------------------------------------------------

    order = np.argsort(x)

    x = x[order]
    y = y[order]

    return x, y


# =====================================================================
# NORMALIZATION
# =====================================================================

def normalize_array(values):
    """
    Normalize array to 0-1.
    """

    values = np.asarray(
        values,
        dtype=float
    )

    minimum = np.nanmin(values)
    maximum = np.nanmax(values)

    difference = maximum - minimum

    if difference < EPSILON:

        return np.zeros_like(values)

    return (
        values - minimum
    ) / difference


# =====================================================================
# BASELINE CORRECTION
# =====================================================================

def baseline_correction(
    x,
    y,
    polynomial_order=3
):
    """
    Polynomial baseline correction.
    """

    x = np.asarray(
        x,
        dtype=float
    )

    y = np.asarray(
        y,
        dtype=float
    )

    if len(x) < polynomial_order + 2:
        return y.copy()

    # Scale x to improve numerical stability
    x_scaled = (
        2 *
        (x - x.min()) /
        (x.max() - x.min() + EPSILON)
        - 1
    )

    try:

        coefficients = np.polyfit(
            x_scaled,
            y,
            polynomial_order
        )

        baseline = np.polyval(
            coefficients,
            x_scaled
        )

        corrected = y - baseline

        return corrected

    except Exception:

        return y.copy()


# =====================================================================
# SMOOTHING
# =====================================================================

def smooth_spectrum(
    y,
    window_length=11
):
    """
    Simple moving-average smoothing.
    """

    y = np.asarray(
        y,
        dtype=float
    )

    if len(y) < 5:
        return y.copy()

    if window_length >= len(y):
        window_length = len(y) - 1

    if window_length < 3:
        return y.copy()

    if window_length % 2 == 0:
        window_length -= 1

    kernel = np.ones(
        window_length
    ) / window_length

    return np.convolve(
        y,
        kernel,
        mode="same"
    )


# =====================================================================
# TRANSMITTANCE / ABSORBANCE
# =====================================================================

def convert_to_absorbance(
    y,
    y_type="auto"
):
    """
    Convert transmittance to absorbance when appropriate.

    If the spectrum appears to be percent transmittance,
    converts:

        A = -log10(%T / 100)

    If it appears to be fractional transmittance,
    converts:

        A = -log10(T)

    Otherwise the values are returned unchanged.
    """

    y = np.asarray(
        y,
        dtype=float
    )

    if y_type == "absorbance":
        return y

    if y_type == "transmittance":

        positive = np.clip(
            y,
            EPSILON,
            None
        )

        if np.nanmax(y) > 1.5:

            return -np.log10(
                np.clip(
                    positive / 100.0,
                    EPSILON,
                    None
                )
            )

        return -np.log10(
            np.clip(
                positive,
                EPSILON,
                None
            )
        )

    # Automatic detection
    maximum = np.nanmax(y)
    minimum = np.nanmin(y)

    # Typical percent transmittance
    if (
        maximum > 1.5
        and maximum <= 110
        and minimum >= 0
    ):

        return -np.log10(
            np.clip(
                y / 100.0,
                EPSILON,
                None
            )
        )

    # Typical fractional transmittance
    if (
        maximum <= 1.5
        and minimum >= 0
    ):

        return -np.log10(
            np.clip(
                y,
                EPSILON,
                None
            )
        )

    # Already appears to be absorbance/intensity
    return y


# =====================================================================
# SAVE CSV
# =====================================================================

def save_spectrum_csv(
    x,
    y,
    output_path
):
    """
    Save processed spectrum as CSV.
    """

    output_path.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    df = pd.DataFrame({
        "Wavenumber": x,
        "Intensity": y
    })

    df.to_csv(
        output_path,
        index=False
    )

    return output_path


# =====================================================================
# INTERPOLATION
# =====================================================================

def resample_spectrum(
    x,
    y,
    target_length=COMPARISON_LENGTH,
    min_x=None,
    max_x=None
):
    """
    Resample spectrum onto an evenly spaced X axis.
    """

    x = np.asarray(
        x,
        dtype=float
    )

    y = np.asarray(
        y,
        dtype=float
    )

    valid = (
        np.isfinite(x) &
        np.isfinite(y)
    )

    x = x[valid]
    y = y[valid]

    if len(x) < 2:

        raise ValueError(
            "Not enough points for resampling."
        )

    order = np.argsort(x)

    x = x[order]
    y = y[order]

    unique_x, indices = np.unique(
        x,
        return_index=True
    )

    x = unique_x
    y = y[indices]

    if min_x is None:
        min_x = max(
            DEFAULT_MIN_WAVENUMBER,
            float(x.min())
        )

    if max_x is None:
        max_x = min(
            DEFAULT_MAX_WAVENUMBER,
            float(x.max())
        )

    if max_x <= min_x:

        min_x = float(x.min())
        max_x = float(x.max())

    target_x = np.linspace(
        min_x,
        max_x,
        target_length
    )

    target_y = np.interp(
        target_x,
        x,
        y
    )

    return target_x, target_y


# =====================================================================
# FFT
# =====================================================================

def calculate_fft(
    x,
    y,
    f
```
