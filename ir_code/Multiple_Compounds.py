from __future__ import annotations

import io
import math
import re
import sqlite3
from pathlib import Path
from urllib.parse import quote

import numpy as np
import pandas as pd
import requests

import tkinter as tk
from tkinter import Tk, filedialog, messagebox
from PIL import Image, ImageTk


from pathlib import Path

DATABASE_PATH = (
    Path.home()
    / "IR_Spectra"
    / "compound_spectra.db"
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


IR_WEIGHT = 0.60
DIAGNOSTIC_WEIGHT = 0.25
FFT_WEIGHT = 0.15



NITRILE_CENTER = 2250.0
NITRILE_HALF_WIDTH = 55.0

NITRILE_WEIGHT = 4.0


DIAGNOSTIC_REGIONS = [

    (
        3650.0,
        3200.0,
        1.0,
        "O-H / N-H"
    ),

    (
        3100.0,
        3000.0,
        1.0,
        "sp2 C-H"
    ),

    (
        3000.0,
        2800.0,
        1.0,
        "sp3 C-H"
    ),

    (
        2305.0,
        2185.0,
        NITRILE_WEIGHT,
        "C#N nitrile"
    ),

    (
        1800.0,
        1650.0,
        1.5,
        "C=O / C=C"
    ),

    (
        1600.0,
        1450.0,
        1.0,
        "C=C / bending"
    ),

    (
        1300.0,
        1000.0,
        1.2,
        "C-O / fingerprint"
    ),

    (
        1000.0,
        650.0,
        1.2,
        "fingerprint"
    ),
]


# ============================================================
# TKINTER HELPERS
# ============================================================

def create_hidden_root():

    root = Tk()

    root.withdraw()

    try:
        root.attributes("-topmost", True)
    except Exception:
        pass

    return root

def select_ir_csv_file():

    root = create_hidden_root()

    path = filedialog.askopenfilename(
        parent=root,
        title="Select UNKNOWN IR CSV file",
        filetypes=[
            ("IR CSV files", "*.csv"),
            ("CSV files", "*.csv"),
            ("All files", "*.*"),
        ],
    )

    root.destroy()

    if not path:
        return None

    return Path(path)


def select_fft_csv_file():

    root = create_hidden_root()

    path = filedialog.askopenfilename(
        parent=root,
        title="Select UNKNOWN FFT CSV file",
        filetypes=[
            ("FFT CSV files", "*.csv"),
            ("CSV files", "*.csv"),
            ("All files", "*.*"),
        ],
    )

    root.destroy()

    if not path:
        return None

    return Path(path)


def select_database_file():
    """
    Return the single, correct IR spectroscopy database.

    Database:
        /home/darestephen/IR_Spectra/compound_spectra.db
    """

    db_path = DATABASE_PATH

    print()
    print("=" * 90)
    print("DATABASE LOCATION")
    print("=" * 90)

    print(
        f"Using database automatically:\n"
        f"{db_path}"
    )

    if not db_path.exists():

        print()
        print("❌ DATABASE NOT FOUND")

        print()
        print(
            f"Expected database:\n"
            f"{db_path}"
        )

        print()
        print(
            "Please make sure the database manager has created:"
        )

        print(
            db_path
        )

        return None

    if not db_path.is_file():

        print()
        print(
            "❌ Database path exists but is not a file:"
        )

        print(
            db_path
        )

        return None

    print()
    print("✅ Database found.")

    return db_path

    # ---------------------------------------------------------
    # CHECK DATABASE IS A FILE
    # ---------------------------------------------------------

    if not db_path.is_file():

        print()
        print(
            "❌ Database path exists but is not a file:"
        )

        print(
            db_path
        )

        return None

    # ---------------------------------------------------------
    # DATABASE FOUND
    # ---------------------------------------------------------

    print()
    print("✅ Database found.")

    print(
        f"Database path:\n"
        f"{db_path}"
    )

    return db_path


# ============================================================
# NUMERIC CLEANING
# ============================================================

def clean_numeric_arrays(x, y):

    x = np.asarray(
        x,
        dtype=float
    )

    y = np.asarray(
        y,
        dtype=float
    )

    mask = (
        np.isfinite(x)
        &
        np.isfinite(y)
    )

    x = x[mask]
    y = y[mask]

    if len(x) < MINIMUM_POINTS:

        raise ValueError(
            "Spectrum contains too few valid "
            "numeric points."
        )

    order = np.argsort(x)

    x = x[order]
    y = y[order]

    unique_x, inverse = np.unique(
        x,
        return_inverse=True
    )

    if len(unique_x) != len(x):

        sums = np.zeros(
            len(unique_x),
            dtype=float
        )

        counts = np.zeros(
            len(unique_x),
            dtype=float
        )

        np.add.at(
            sums,
            inverse,
            y
        )

        np.add.at(
            counts,
            inverse,
            1.0
        )

        y = (
            sums
            /
            np.maximum(
                counts,
                1.0
            )
        )

        x = unique_x

    if len(x) < MINIMUM_POINTS:

        raise ValueError(
            "Spectrum contains too few "
            "unique points."
        )

    return x, y


# ============================================================
# NORMALIZATION
# ============================================================

def normalize_01(y):

    y = np.asarray(
        y,
        dtype=float
    )

    if len(y) == 0:
        return y

    lo = float(
        np.nanmin(y)
    )

    hi = float(
        np.nanmax(y)
    )

    if hi - lo < EPSILON:

        return np.zeros_like(y)

    return (
        (y - lo)
        /
        (hi - lo)
    )


def normalize_vector(v):

    v = np.asarray(
        v,
        dtype=float
    )

    norm = np.linalg.norm(v)

    if norm < EPSILON:

        return np.zeros_like(v)

    return v / norm


# ============================================================
# BASELINE CORRECTION
# ============================================================

def baseline_correct(
    y,
    window_fraction=0.08
):

    y = np.asarray(
        y,
        dtype=float
    ).copy()

    if len(y) < 20:

        return y - np.mean(y)

    window = max(
        3,
        int(
            len(y)
            *
            window_fraction
        )
    )

    if window % 2 == 0:

        window += 1

    maximum_window = (
        len(y)
        if len(y) % 2
        else len(y) - 1
    )

    window = min(
        window,
        maximum_window
    )

    s = pd.Series(y)

    baseline = (
        s.rolling(
            window=window,
            center=True,
            min_periods=1
        )
        .quantile(0.10)
    )

    baseline = (
        baseline
        .interpolate(
            limit_direction="both"
        )
        .to_numpy(
            dtype=float
        )
    )

    return y - baseline


# ============================================================
# SMOOTHING
# ============================================================

def smooth_signal(
    y,
    window=7
):

    y = np.asarray(
        y,
        dtype=float
    )

    if (
        len(y) < window
        or window < 3
    ):

        return y.copy()

    if window % 2 == 0:

        window += 1

    kernel = (
        np.ones(window)
        /
        window
    )

    padded = np.pad(
        y,
        window // 2,
        mode="edge"
    )

    return np.convolve(
        padded,
        kernel,
        mode="valid"
    )


# ============================================================
# TRANSMITTANCE / ABSORBANCE
# ============================================================

def convert_to_absorbance(
    y,
    y_type="transmittance"
):

    y = np.asarray(
        y,
        dtype=float
    )

    kind = str(
        y_type or ""
    ).lower()

    if "absorb" in kind:

        return np.maximum(
            y,
            0.0
        )

    if (
        "percent" in kind
        or "%"
        in kind
    ):

        t = np.clip(
            y / 100.0,
            1e-8,
            1.0
        )

        return -np.log10(t)

    finite = y[
        np.isfinite(y)
    ]

    if (
        len(finite)
        and
        np.nanmax(finite) > 1.5
    ):

        t = np.clip(
            y / 100.0,
            1e-8,
            1.0
        )

        return -np.log10(t)

    t = np.clip(
        y,
        1e-8,
        1.0
    )

    return -np.log10(t)


def absorbance_to_transmittance(
    absorbance
):

    return np.power(
        10.0,
        -np.asarray(
            absorbance,
            dtype=float
        )
    )


# ============================================================
# CSV COLUMN DETECTION
# ============================================================

def _find_column(
    df,
    patterns,
    exclude=None
):

    exclude = set(
        exclude or []
    )

    normalized = {
        str(c): str(c)
        .strip()
        .lower()
        for c in df.columns
    }

    for pattern in patterns:

        regex = re.compile(
            pattern,
            re.I
        )

        for original, lower in normalized.items():

            if original in exclude:

                continue

            if regex.search(lower):

                return original

    return None


# ============================================================
# CSV TYPE
# ============================================================

def identify_csv_type(df):

    columns = [
        str(c)
        .strip()
        .lower()
        for c in df.columns
    ]

    joined = " ".join(
        columns
    )

    if (
        "fft" in joined
        or
        "frequency" in joined
        or
        "magnitude" in joined
    ):

        return "FFT"

    if (
        "wavenumber" in joined
        or
        "cm" in joined
        or
        "transmittance" in joined
        or
        "absorbance" in joined
    ):

        return "IR"

    return "UNKNOWN"


# ============================================================
# LOAD CSV
# ============================================================

def load_csv_file(
    path_or_buffer,
    expected_type=None
):

    try:

        df = pd.read_csv(
            path_or_buffer
        )

    except Exception as error:

        raise ValueError(
            f"Could not read CSV: {error}"
        )

    if df.empty:

        raise ValueError(
            "CSV file is empty."
        )

    df = df.dropna(
        axis=1,
        how="all"
    )

    if len(df.columns) < 2:

        raise ValueError(
            "CSV must contain at least "
            "two columns."
        )

    # --------------------------------------------------------
    # IR
    # --------------------------------------------------------

    if expected_type == "IR":

        x_col = _find_column(
            df,
            [
                r"wavenumber",
                r"wave.?number",
                r"cm.?1",
                r"cm-1",
                r"x$",
                r"^x",
            ],
        )

        y_col = _find_column(
            df,
            [
                r"transmittance",
                r"transmission",
                r"absorbance",
                r"intensity",
                r"y$",
                r"^y",
            ],
            exclude=(
                [x_col]
                if x_col
                else None
            ),
        )

        if y_col and (
            "absorb"
            in str(y_col).lower()
        ):

            y_type = "absorbance"

        elif y_col and (
            "percent"
            in str(y_col).lower()
            or "%"
            in str(y_col)
        ):

            y_type = "percent transmittance"

        else:

            y_type = "transmittance"

    # --------------------------------------------------------
    # FFT
    # --------------------------------------------------------

    elif expected_type == "FFT":

        x_col = _find_column(
            df,
            [
                r"frequency",
                r"freq",
                r"fft.?x",
                r"x$",
                r"^x",
            ],
        )

        y_col = _find_column(
            df,
            [
                r"magnitude",
                r"amplitude",
                r"fft.?y",
                r"fft",
                r"intensity",
                r"y$",
                r"^y",
            ],
            exclude=(
                [x_col]
                if x_col
                else None
            ),
        )

        y_type = "magnitude"

    else:

        x_col = df.columns[0]

        y_col = (
            df.columns[1]
            if len(df.columns) > 1
            else None
        )

        y_type = "unknown"

    # --------------------------------------------------------
    # FALLBACK
    # --------------------------------------------------------

    if (
        x_col is None
        or
        y_col is None
    ):

        numeric_candidates = []

        for column in df.columns:

            numeric = pd.to_numeric(
                df[column],
                errors="coerce"
            )

            if numeric.notna().sum() >= MINIMUM_POINTS:

                numeric_candidates.append(
                    column
                )

        if len(numeric_candidates) >= 2:

            x_col = numeric_candidates[0]
            y_col = numeric_candidates[1]

        else:

            raise ValueError(
                "Could not identify two numeric "
                "CSV columns."
            )

    # --------------------------------------------------------
    # NUMERIC DATA
    # --------------------------------------------------------

    x = pd.to_numeric(
        df[x_col],
        errors="coerce"
    ).to_numpy(
        dtype=float
    )

    y = pd.to_numeric(
        df[y_col],
        errors="coerce"
    ).to_numpy(
        dtype=float
    )

    x, y = clean_numeric_arrays(
        x,
        y
    )

    return x, y, y_type


# ============================================================
# LOAD CSV TEXT STORED IN SQLITE
# ============================================================

def load_database_data(
    value,
    expected_type
):

    if value is None:

        raise ValueError(
            "Database CSV field is empty."
        )

    if isinstance(
        value,
        bytes
    ):

        value = value.decode(
            "utf-8",
            errors="replace"
        )

    text = str(value)

    candidate = Path(
        text
    ).expanduser()

    # If a file path was stored.
    if (
        len(text) < 1000
        and
        candidate.exists()
        and
        candidate.is_file()
    ):

        return load_csv_file(
            candidate,
            expected_type
        )

    # Otherwise it is CSV text.
    return load_csv_file(
        io.StringIO(text),
        expected_type
    )


# ============================================================
# FFT CALCULATION
#
# Used only for residual mixture rounds.
#
# The initial unknown FFT is ALWAYS loaded from the
# user's selected FFT CSV.
# ============================================================

def calculate_fft(
    x,
    y
):

    x, y = clean_numeric_arrays(
        x,
        y
    )

    uniform_x = np.linspace(
        float(x.min()),
        float(x.max()),
        len(x)
    )

    uniform_y = np.interp(
        uniform_x,
        x,
        y
    )

    uniform_y = smooth_signal(
        uniform_y,
        5
    )

    uniform_y -= np.mean(
        uniform_y
    )

    fft_values = np.fft.rfft(
        uniform_y
    )

    magnitude = np.abs(
        fft_values
    )

    frequency = np.fft.rfftfreq(
        len(uniform_y),
        d=(
            uniform_x[1]
            -
            uniform_x[0]
        )
    )

    if len(magnitude) > 1:

        frequency = frequency[1:]
        magnitude = magnitude[1:]

    return frequency, magnitude


# ============================================================
# SPECTRUM RESAMPLING
# ============================================================

def create_common_axis(
    x1,
    x2,
    length=COMPARISON_LENGTH
):

    x1 = np.asarray(
        x1,
        dtype=float
    )

    x2 = np.asarray(
        x2,
        dtype=float
    )

    lo = max(
        float(np.min(x1)),
        float(np.min(x2)),
        DEFAULT_MIN_WAVENUMBER
    )

    hi = min(
        float(np.max(x1)),
        float(np.max(x2)),
        DEFAULT_MAX_WAVENUMBER
    )

    if hi <= lo:

        lo = max(
            float(np.min(x1)),
            float(np.min(x2))
        )

        hi = min(
            float(np.max(x1)),
            float(np.max(x2))
        )

    if hi <= lo:

        raise ValueError(
            "Spectra do not have a common range."
        )

    return np.linspace(
        lo,
        hi,
        length
    )


def resample_spectrum(
    x,
    y,
    axis
):

    x, y = clean_numeric_arrays(
        x,
        y
    )

    return np.interp(
        axis,
        x,
        y
    )


# ============================================================
# SIMILARITY FUNCTIONS
# ============================================================

def calculate_correlation(
    a,
    b
):

    a = np.asarray(
        a,
        dtype=float
    )

    b = np.asarray(
        b,
        dtype=float
    )

    if (
        len(a) < 2
        or
        len(b) < 2
    ):

        return 0.0

    a = a - np.mean(a)
    b = b - np.mean(b)

    da = np.linalg.norm(a)
    db = np.linalg.norm(b)

    if (
        da < EPSILON
        or
        db < EPSILON
    ):

        return 0.0

    return float(
        np.clip(
            np.dot(a, b)
            /
            (
                da * db
            ),
            -1.0,
            1.0
        )
    )


def calculate_cosine_similarity(
    a,
    b
):

    a = np.asarray(
        a,
        dtype=float
    )

    b = np.asarray(
        b,
        dtype=float
    )

    a_norm = np.linalg.norm(a)
    b_norm = np.linalg.norm(b)

    if (
        a_norm < EPSILON
        or
        b_norm < EPSILON
    ):

        return 0.0

    return float(
        np.clip(
            np.dot(a, b)
            /
            (
                a_norm * b_norm
            ),
            -1.0,
            1.0
        )
    )


def calculate_derivative_similarity(
    a,
    b
):

    a = np.asarray(
        a,
        dtype=float
    )

    b = np.asarray(
        b,
        dtype=float
    )

    da = np.gradient(a)
    db = np.gradient(b)

    return max(
        0.0,
        calculate_correlation(
            da,
            db
        )
    )


def calculate_normalized_mse(
    a,
    b
):

    a = np.asarray(
        a,
        dtype=float
    )

    b = np.asarray(
        b,
        dtype=float
    )

    scale = max(
        np.ptp(a),
        np.ptp(b),
        EPSILON
    )

    return float(
        np.mean(
            (
                (a - b)
                /
                scale
            )
            ** 2
        )
    )


def mse_to_score(
    mse
):

    return float(
        np.clip(
            100.0
            *
            np.exp(
                -12.0
                *
                max(
                    0.0,
                    mse
                )
            ),
            0.0,
            100.0
        )
    )


# ============================================================
# DIAGNOSTIC REGION
# ============================================================

def region_similarity(
    unknown_axis,
    unknown_y,
    db_axis,
    db_y,
    lo,
    hi
):

    low = max(
        min(lo, hi),
        float(np.min(unknown_axis)),
        float(np.min(db_axis))
    )

    high = min(
        max(lo, hi),
        float(np.max(unknown_axis)),
        float(np.max(db_axis))
    )

    if high <= low:

        return None

    axis = np.linspace(
        low,
        high,
        256
    )

    u = np.interp(
        axis,
        unknown_axis,
        unknown_y
    )

    d = np.interp(
        axis,
        db_axis,
        db_y
    )

    u_corr = baseline_correct(u)
    d_corr = baseline_correct(d)

    corr = max(
        0.0,
        calculate_correlation(
            u_corr,
            d_corr
        )
    )

    cosine = max(
        0.0,
        calculate_cosine_similarity(
            u_corr,
            d_corr
        )
    )

    u_strength = float(
        np.ptp(u_corr)
    )

    d_strength = float(
        np.ptp(d_corr)
    )

    if (
        u_strength < 0.005
        or
        d_strength < 0.005
    ):

        return None

    ratio = min(
        u_strength
        /
        max(
            d_strength,
            EPSILON
        ),

        d_strength
        /
        max(
            u_strength,
            EPSILON
        )
    )

    return float(
        0.55 * corr
        +
        0.25 * cosine
        +
        0.20 * np.clip(
            ratio,
            0.0,
            1.0
        )
    )


def calculate_diagnostic_similarity(
    unknown_axis,
    unknown_y,
    db_axis,
    db_y
):

    scores = []
    weights = []

    for (
        lo,
        hi,
        weight,
        name
    ) in DIAGNOSTIC_REGIONS:

        score = region_similarity(
            unknown_axis,
            unknown_y,
            db_axis,
            db_y,
            lo,
            hi
        )

        if score is not None:

            scores.append(
                score
            )

            weights.append(
                weight
            )

    if not scores:

        return 0.0

    return float(
        np.average(
            scores,
            weights=weights
        )
    )


def calculate_nitrile_score(
    unknown_axis,
    unknown_y,
    db_axis,
    db_y
):

    lo = (
        NITRILE_CENTER
        -
        NITRILE_HALF_WIDTH
    )

    hi = (
        NITRILE_CENTER
        +
        NITRILE_HALF_WIDTH
    )

    score = region_similarity(
        unknown_axis,
        unknown_y,
        db_axis,
        db_y,
        lo,
        hi
    )

    if score is None:

        return 0.0

    return float(score)


# ============================================================
# IR COMPARISON
# ============================================================

def compare_ir_spectra(
    unknown_x,
    unknown_abs,
    db_x,
    db_abs
):

    axis = create_common_axis(
        unknown_x,
        db_x
    )

    u = resample_spectrum(
        unknown_x,
        unknown_abs,
        axis
    )

    d = resample_spectrum(
        db_x,
        db_abs,
        axis
    )

    u = np.maximum(
        u,
        0.0
    )

    d = np.maximum(
        d,
        0.0
    )

    distance = calculate_normalized_mse(
        u,
        d
    )

    rmse = math.sqrt(
        distance
    )

    correlation = max(
        0.0,
        calculate_correlation(
            u,
            d
        )
    )

    cosine = max(
        0.0,
        calculate_cosine_similarity(
            u,
            d
        )
    )

    derivative = (
        calculate_derivative_similarity(
            u,
            d
        )
    )

    diagnostic = (
        calculate_diagnostic_similarity(
            axis,
            u,
            axis,
            d
        )
    )

    nitrile = (
        calculate_nitrile_score(
            axis,
            u,
            axis,
            d
        )
    )

    mse_score = (
        mse_to_score(
            distance
        )
        /
        100.0
    )

    diagnostic_enhanced = float(
        np.clip(
            0.65 * diagnostic
            +
            0.35 * nitrile,
            0.0,
            1.0
        )
    )

    ir_score = 100.0 * (
        0.25 * correlation
        +
        0.20 * cosine
        +
        0.15 * mse_score
        +
        0.15 * derivative
        +
        0.25 * diagnostic_enhanced
    )

    return {

        "distance":
            distance,

        "rmse":
            rmse,

        "correlation":
            correlation,

        "cosine":
            cosine,

        "derivative":
            derivative,

        "diagnostic":
            diagnostic,

        "nitrile":
            nitrile,

        "ir_score":
            float(
                np.clip(
                    ir_score,
                    0.0,
                    100.0
                )
            ),
    }


# ============================================================
# FFT COMPARISON
#
# IMPORTANT:
# The unknown FFT comes directly from the selected FFT CSV.
#
# We normalize both FFT axes to 0-1 before comparison.
# This prevents otherwise valid FFT spectra from receiving
# 0 simply because their absolute frequency ranges differ.
# ============================================================

def normalize_fft_axis(x):

    x = np.asarray(
        x,
        dtype=float
    )

    lo = float(
        np.min(x)
    )

    hi = float(
        np.max(x)
    )

    if hi - lo < EPSILON:

        return np.linspace(
            0.0,
            1.0,
            len(x)
        )

    return (
        (x - lo)
        /
        (hi - lo)
    )


def compare_fft_spectra(
    unknown_x,
    unknown_y,
    db_x,
    db_y
):

    try:

        unknown_x, unknown_y = (
            clean_numeric_arrays(
                unknown_x,
                unknown_y
            )
        )

        db_x, db_y = (
            clean_numeric_arrays(
                db_x,
                db_y
            )
        )

    except Exception:

        return 0.0

    # --------------------------------------------------------
    # Normalize frequency axes independently.
    # --------------------------------------------------------

    u_axis = normalize_fft_axis(
        unknown_x
    )

    d_axis = normalize_fft_axis(
        db_x
    )

    common_axis = np.linspace(
        0.0,
        1.0,
        FFT_COMPARISON_LENGTH
    )

    u = np.interp(
        common_axis,
        u_axis,
        unknown_y
    )

    d = np.interp(
        common_axis,
        d_axis,
        db_y
    )

    # --------------------------------------------------------
    # Remove invalid values.
    # --------------------------------------------------------

    u = np.nan_to_num(
        u,
        nan=0.0,
        posinf=0.0,
        neginf=0.0
    )

    d = np.nan_to_num(
        d,
        nan=0.0,
        posinf=0.0,
        neginf=0.0
    )

    # --------------------------------------------------------
    # FFT magnitude is normally highly dynamic.
    # Log compression prevents the first few large
    # components from dominating the entire comparison.
    # --------------------------------------------------------

    u = np.log1p(
        np.maximum(
            u,
            0.0
        )
    )

    d = np.log1p(
        np.maximum(
            d,
            0.0
        )
    )

    u = normalize_01(u)
    d = normalize_01(d)

    # --------------------------------------------------------
    # If one spectrum is completely flat, correlation and
    # cosine cannot provide useful information.
    # --------------------------------------------------------

    if (
        np.ptp(u) < EPSILON
        or
        np.ptp(d) < EPSILON
    ):

        return 0.0

    corr = max(
        0.0,
        calculate_correlation(
            u,
            d
        )
    )

    cosine = max(
        0.0,
        calculate_cosine_similarity(
            u,
            d
        )
    )

    mse = calculate_normalized_mse(
        u,
        d
    )

    mse_score = (
        mse_to_score(
            mse
        )
        /
        100.0
    )

    score = 100.0 * (
        0.45 * corr
        +
        0.40 * cosine
        +
        0.15 * mse_score
    )

    return float(
        np.clip(
            score,
            0.0,
            100.0
        )
    )


# ============================================================
# DATABASE COLUMN DETECTION
# ============================================================

def database_column_map(
    columns
):

    lower = {
        str(c)
        .strip()
        .lower(): c
        for c in columns
    }

    def first_match(
        patterns
    ):

        for pattern in patterns:

            regex = re.compile(
                pattern,
                re.I
            )

            for low, original in lower.items():

                if regex.search(low):

                    return original

        return None

    name_col = first_match([
        r"^compound_name$",
        r"^compound$",
        r"name",
        r"title",
    ])

    ir_col = first_match([
        r"^ir_data$",
        r"^ir_csv$",
        r"ir.?csv",
        r"infrared",
        r"ir_spectrum",
    ])

    fft_col = first_match([
        r"^fft_data$",
        r"^fft_csv$",
        r"ir.?fft",
        r"fft.?csv",
        r"fourier",
    ])

    return (
        name_col,
        ir_col,
        fft_col
    )


# ============================================================
# READ DATABASE
# ============================================================

def load_database_records(
    db_path
):

    print()
    print(
        "=" * 90
    )

    print(
        "LOADING DATABASE"
    )

    print(
        "=" * 90
    )

    print(
        f"Database: {db_path}"
    )

    if not db_path.exists():

        print(
            "❌ Database file does not exist."
        )

        return []

    try:

        with sqlite3.connect(
            str(db_path)
        ) as conn:

            tables = conn.execute(
                """
                SELECT name
                FROM sqlite_master
                WHERE type='table'
                """
            ).fetchall()

            table_names = [
                row[0]
                for row in tables
            ]

            if (
                "compound_spectra"
                not in table_names
            ):

                print(
                    "❌ compound_spectra table "
                    "was not found."
                )

                print(
                    "Tables found:",
                    table_names
                )

                return []

            df = pd.read_sql_query(
                "SELECT * FROM compound_spectra",
                conn
            )

    except Exception as error:

        print(
            f"❌ Could not read database: "
            f"{error}"
        )

        return []

    if df.empty:

        print(
            "❌ Database table is empty."
        )

        return []

    (
        name_col,
        ir_col,
        fft_col
    ) = database_column_map(
        df.columns
    )

    print(
        f"Compound column: {name_col}"
    )

    print(
        f"IR CSV column:   {ir_col}"
    )

    print(
        f"FFT CSV column:  {fft_col}"
    )

    if name_col is None:

        print(
            "❌ Compound-name column "
            "not found."
        )

        return []

    if ir_col is None:

        print(
            "❌ IR CSV column "
            "not found."
        )

        return []

    if fft_col is None:

        print(
            "❌ FFT CSV column "
            "not found."
        )

        return []

    records = []

    for index, row in df.iterrows():

        name = str(
            row[name_col]
        ).strip()

        if (
            not name
            or
            name.lower() == "nan"
        ):

            continue

        try:

            ir_x, ir_y, ir_type = (
                load_database_data(
                    row[ir_col],
                    "IR"
                )
            )

            fft_x, fft_y, _ = (
                load_database_data(
                    row[fft_col],
                    "FFT"
                )
            )

            ir_abs = (
                convert_to_absorbance(
                    ir_y,
                    ir_type
                )
            )

            records.append({

                "database_row":
                    index,

                "name":
                    name,

                "ir_x":
                    ir_x,

                "ir_y":
                    ir_y,

                "ir_abs":
                    ir_abs,

                "fft_x":
                    fft_x,

                "fft_y":
                    fft_y,
            })

        except Exception as error:

            print()
            print(
                f"⚠️ Skipping {name}"
            )

            print(
                f"   Reason: {error}"
            )

    print()
    print(
        f"Database records: {len(df)}"
    )

    print(
        f"Usable IR/FFT records: "
        f"{len(records)}"
    )

    return records


# ============================================================
# SELECT DATABASE COMPOUNDS
# ============================================================

def select_database_compounds(
    records
):

    if not records:

        return []

    print()
    print(
        "=" * 90
    )

    print(
        "DATABASE COMPOUNDS"
    )

    print(
        "=" * 90
    )

    print(
        "The following compounds contain "
        "usable IR CSV + FFT CSV data:"
    )

    print()

    for i, record in enumerate(
        records,
        start=1
    ):

        print(
            f"{i:>3}. "
            f"{record['name']}"
        )

    print()
    print(
        "You can enter:"
    )

    print(
        "  ALL       = compare every compound"
    )

    print(
        "  1,2,5     = compare selected compounds"
    )

    print(
        "  1-5       = compare a range"
    )

    print()

    while True:

        selection = input(
            "Select database compounds [ALL]: "
        ).strip()

        if not selection:

            selection = "ALL"

        if selection.upper() == "ALL":

            return records.copy()

        try:

            indexes = set()

            parts = selection.split(",")

            for part in parts:

                part = part.strip()

                if "-" in part:

                    start, end = (
                        part.split(
                            "-",
                            1
                        )
                    )

                    start = int(
                        start.strip()
                    )

                    end = int(
                        end.strip()
                    )

                    if start > end:

                        start, end = (
                            end,
                            start
                        )

                    for n in range(
                        start,
                        end + 1
                    ):

                        indexes.add(n)

                else:

                    indexes.add(
                        int(part)
                    )

            if not indexes:

                raise ValueError

            if any(
                n < 1
                or
                n > len(records)
                for n in indexes
            ):

                raise ValueError

            selected = [
                records[n - 1]
                for n in sorted(
                    indexes
                )
            ]

            print()
            print(
                "Selected database compounds:"
            )

            for record in selected:

                print(
                    f"  ✓ "
                    f"{record['name']}"
                )

            return selected

        except Exception:

            print(
                "❌ Invalid selection."
            )

            print(
                "Example: ALL or "
                "1,2,5 or 1-5"
            )


# ============================================================
# COMPARE ONE DATABASE CANDIDATE
# ============================================================

def compare_candidate(
    unknown_ir_x,
    unknown_ir_abs,
    unknown_fft_x,
    unknown_fft_y,
    record
):

    ir = compare_ir_spectra(
        unknown_ir_x,
        unknown_ir_abs,
        record["ir_x"],
        record["ir_abs"]
    )

    fft_score = compare_fft_spectra(
        unknown_fft_x,
        unknown_fft_y,
        record["fft_x"],
        record["fft_y"]
    )

    combined = (
        IR_WEIGHT
        *
        ir["ir_score"]
        +
        DIAGNOSTIC_WEIGHT
        *
        ir["diagnostic"]
        *
        100.0
        +
        FFT_WEIGHT
        *
        fft_score
    )

    # --------------------------------------------------------
    # NITRILE BONUS
    # --------------------------------------------------------

    if ir["nitrile"] >= 0.75:

        combined += (
            4.0
            *
            ir["nitrile"]
        )

    combined = float(
        np.clip(
            combined,
            0.0,
            100.0
        )
    )

    return {

        "name":
            record["name"],

        "distance":
            ir["distance"],

        "rmse":
            ir["rmse"],

        "correlation":
            ir["correlation"],

        "cosine":
            ir["cosine"],

        "derivative":
            ir["derivative"],

        "diagnostic":
            ir["diagnostic"],

        "nitrile":
            ir["nitrile"],

        "ir_score":
            ir["ir_score"],

        "fft_score":
            fft_score,

        "similarity":
            combined,

        "record":
            record,
    }


# ============================================================
# DATABASE COMPARISON
# ============================================================

def compare_database(
    unknown_ir_x,
    unknown_ir_abs,
    unknown_fft_x,
    unknown_fft_y,
    database
):

    candidates = []

    for record in database:

        try:

            result = compare_candidate(
                unknown_ir_x,
                unknown_ir_abs,
                unknown_fft_x,
                unknown_fft_y,
                record
            )

            candidates.append(
                result
            )

        except Exception as error:

            print(
                f"⚠️ Could not compare "
                f"{record['name']}: "
                f"{error}"
            )

    candidates.sort(
        key=lambda item:
        item["similarity"],
        reverse=True
    )

    return candidates


# ============================================================
# RESIDUAL FITTING
# ============================================================

def fit_candidate(
    residual,
    candidate,
    residual_axis
):

    candidate_axis = np.asarray(
        candidate["record"]["ir_x"],
        dtype=float
    )

    candidate_abs = np.asarray(
        candidate["record"]["ir_abs"],
        dtype=float
    )

    # --------------------------------------------------------
    # Put candidate on residual axis.
    # --------------------------------------------------------

    candidate_on_residual = np.interp(
        residual_axis,
        candidate_axis,
        candidate_abs
    )

    template = baseline_correct(
        candidate_on_residual
    )

    template = np.maximum(
        template,
        0.0
    )

    target = baseline_correct(
        residual
    )

    target = np.maximum(
        target,
        0.0
    )

    denominator = float(
        np.dot(
            template,
            template
        )
    )

    if denominator < EPSILON:

        return (
            0.0,
            residual.copy(),
            0.0
        )

    coefficient = max(
        0.0,
        float(
            np.dot(
                target,
                template
            )
            /
            denominator
        )
    )

    subtraction = (
        coefficient
        *
        template
    )

    old_energy = float(
        np.mean(
            target ** 2
        )
    )

    new_residual = np.maximum(
        residual
        -
        subtraction,
        0.0
    )

    new_target = baseline_correct(
        new_residual
    )

    new_target = np.maximum(
        new_target,
        0.0
    )

    new_energy = float(
        np.mean(
            new_target ** 2
        )
    )

    if old_energy < EPSILON:

        reduction = 0.0

    else:

        reduction = (
            1.0
            -
            new_energy
            /
            old_energy
        )

    reduction = float(
        np.clip(
            reduction,
            0.0,
            1.0
        )
    )

    return (
        coefficient,
        new_residual,
        reduction
    )


# ============================================================
# USER INPUT
# ============================================================

def ask_integer(
    prompt,
    minimum,
    maximum,
    default=None
):

    while True:

        raw = input(
            prompt
        ).strip()

        if (
            raw == ""
            and
            default is not None
        ):

            return default

        try:

            value = int(
                raw
            )

        except ValueError:

            print(
                "Please enter a valid integer."
            )

            continue

        if (
            minimum
            <= value
            <= maximum
        ):

            return value

        print(
            f"Enter a number between "
            f"{minimum} and "
            f"{maximum}."
        )


def ask_float(
    prompt,
    minimum,
    maximum,
    default=None
):

    while True:

        raw = input(
            prompt
        ).strip()

        if (
            raw == ""
            and
            default is not None
        ):

            return default

        try:

            value = float(
                raw
            )

        except ValueError:

            print(
                "Please enter a valid number."
            )

            continue

        if (
            minimum
            <= value
            <= maximum
        ):

            return value

        print(
            f"Enter a value between "
            f"{minimum} and "
            f"{maximum}."
        )


# ============================================================
# PUBCHEM
# ============================================================

def get_pubchem_information(
    compound_name
):

    print(
        f"\nSearching PubChem for: "
        f"{compound_name}"
    )

    encoded = quote(
        str(compound_name),
        safe=""
    )

    url = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/"
        f"compound/name/{encoded}/"
        "property/Title,MolecularFormula,"
        "ConnectivitySMILES,CanonicalSMILES/JSON"
    )

    try:

        response = requests.get(
            url,
            timeout=REQUEST_TIMEOUT
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

            return None

        compound = properties[0]

        cid = compound.get(
            "CID"
        )

        title = compound.get(
            "Title",
            compound_name
        )

        formula = compound.get(
            "MolecularFormula",
            "Not available"
        )

        smiles = compound.get(
            "ConnectivitySMILES",
            compound.get(
                "CanonicalSMILES",
                "Not available"
            )
        )

        image_data = None

        if cid is not None:

            image_url = (
                "https://pubchem.ncbi.nlm.nih.gov/rest/pug/"
                f"compound/cid/{cid}/PNG?image_size=400x400"
            )

            try:

                image_response = requests.get(
                    image_url,
                    timeout=REQUEST_TIMEOUT
                )

                image_response.raise_for_status()

                image_data = (
                    image_response.content
                )

            except Exception as error:

                print(
                    f"⚠️ Structure image unavailable: "
                    f"{error}"
                )

        return {

            "cid":
                cid,

            "title":
                title,

            "formula":
                formula,

            "smiles":
                smiles,

            "image_data":
                image_data,
        }

    except Exception as error:

        print(
            f"⚠️ PubChem lookup failed: "
            f"{error}"
        )

        return None


# ============================================================
# STRUCTURE WINDOW
# ============================================================

def show_structure(
    compounds
):

    if not compounds:

        return

    root = Tk()

    root.title(
        "Identified Mixture Components"
    )

    root.geometry(
        "950x950"
    )

    root.minsize(
        700,
        600
    )

    canvas = tk.Canvas(
        root
    )

    scrollbar = tk.Scrollbar(
        root,
        orient="vertical",
        command=canvas.yview
    )

    scroll_frame = tk.Frame(
        canvas
    )

    canvas_window = (
        canvas.create_window(
            (0, 0),
            window=scroll_frame,
            anchor="nw"
        )
    )

    def update_scroll_region(
        event=None
    ):

        canvas.configure(
            scrollregion=canvas.bbox(
                "all"
            )
        )

    def resize_frame(
        event
    ):

        canvas.itemconfigure(
            canvas_window,
            width=event.width
        )

    scroll_frame.bind(
        "<Configure>",
        update_scroll_region
    )

    canvas.bind(
        "<Configure>",
        resize_frame
    )

    canvas.configure(
        yscrollcommand=scrollbar.set
    )

    canvas.pack(
        side="left",
        fill="both",
        expand=True
    )

    scrollbar.pack(
        side="right",
        fill="y"
    )

    tk.Label(
        scroll_frame,
        text="Identified Mixture Components",
        font=(
            "Arial",
            20,
            "bold"
        )
    ).pack(
        pady=15
    )

    tk.Label(
        scroll_frame,
        text="Structures retrieved from PubChem",
        font=(
            "Arial",
            11
        )
    ).pack(
        pady=(0, 10)
    )

    photos = []

    for position, compound in enumerate(
        compounds,
        start=1
    ):

        pubchem = compound.get(
            "pubchem"
        )

        frame = tk.Frame(
            scroll_frame,
            bd=2,
            relief="groove",
            padx=15,
            pady=15
        )

        frame.pack(
            fill="x",
            padx=20,
            pady=10
        )

        tk.Label(
            frame,
            text=(
                f"{position}. "
                f"{compound['name']}"
            ),
            font=(
                "Arial",
                16,
                "bold"
            )
        ).pack(
            pady=(0, 8)
        )

        tk.Label(
            frame,
            text=(
                f"Combined: "
                f"{compound['similarity']:.2f}%\n"
                f"IR: "
                f"{compound['ir_score']:.2f}%   "
                f"FFT: "
                f"{compound['fft_score']:.2f}%\n"
                f"Nitrile: "
                f"{compound['nitrile'] * 100:.2f}%\n"
                f"Residual reduction: "
                f"{compound['residual_reduction'] * 100:.2f}%"
            ),
            font=(
                "Arial",
                11
            )
        ).pack()

        if pubchem:

            tk.Label(
                frame,
                text=(
                    f"PubChem CID: "
                    f"{pubchem.get('cid', 'N/A')}"
                )
            ).pack()

            tk.Label(
                frame,
                text=(
                    f"Formula: "
                    f"{pubchem.get('formula', 'N/A')}"
                )
            ).pack()

            tk.Label(
                frame,
                text=(
                    f"SMILES: "
                    f"{pubchem.get('smiles', 'N/A')}"
                ),
                wraplength=820,
                justify="left"
            ).pack(
                pady=5
            )

            image_data = pubchem.get(
                "image_data"
            )

            if image_data:

                try:

                    image = Image.open(
                        io.BytesIO(
                            image_data
                        )
                    )

                    image.thumbnail(
                        (350, 350)
                    )

                    photo = ImageTk.PhotoImage(
                        image
                    )

                    photos.append(
                        photo
                    )

                    tk.Label(
                        frame,
                        image=photo
                    ).pack(
                        pady=8
                    )

                except Exception as error:

                    tk.Label(
                        frame,
                        text=(
                            "Could not display "
                            f"structure: {error}"
                        )
                    ).pack()

        else:

            tk.Label(
                frame,
                text=(
                    "PubChem information "
                    "unavailable."
                )
            ).pack()

    root._photos = photos

    root.mainloop()


# ============================================================
# PRINT RESULTS
# ============================================================

def print_candidate_table(
    candidates,
    title
):

    print()
    print(
        "=" * 160
    )

    print(title)

    print(
        "=" * 160
    )

    print(
        f"{'Rank':<6}"
        f"{'Compound':<30}"
        f"{'Combined':<12}"
        f"{'IR':<10}"
        f"{'FFT':<10}"
        f"{'Diagnostic':<12}"
        f"{'Nitrile':<10}"
        f"{'Corr.':<10}"
    )

    print(
        "-" * 160
    )

    for position, item in enumerate(
        candidates[:10],
        start=1
    ):

        print(
            f"{position:<6}"
            f"{item['name'][:28]:<30}"
            f"{item['similarity']:<12.2f}"
            f"{item['ir_score']:<10.2f}"
            f"{item['fft_score']:<10.2f}"
            f"{item['diagnostic'] * 100:<12.2f}"
            f"{item['nitrile'] * 100:<10.2f}"
            f"{item['correlation'] * 100:<10.2f}"
        )


# ============================================================
# MAIN MIXTURE IDENTIFICATION
# ============================================================

def identify_mixture_compounds():

    print()
    print(
        "=" * 90
    )

    print(
        "MULTI-COMPOUND IR MIXTURE IDENTIFICATION"
    )

    print(
        "=" * 90
    )

    print(
        "Unknown: existing IR CSV + existing FFT CSV"
    )

    print(
        "Database: SQLite → IR CSV + FFT CSV"
    )

    print(
        "JCAMP-DX files are NOT required."
    )

    print(
        "Unknown IR and FFT CSV files are selected separately."
    )

    print(
        "Nitrile region around 2250 cm-1 receives extra weighting."
    )

    print(
        "=" * 90
    )

    # ========================================================
    # SELECT UNKNOWN IR CSV
    # ========================================================

    print()
    print(
        "=" * 90
    )

    print(
        "SELECT UNKNOWN IR CSV"
    )

    print(
        "=" * 90
    )

    mixture_ir_file = (
        select_ir_csv_file()
    )

    if mixture_ir_file is None:

        print(
            "❌ No unknown IR CSV selected."
        )

        return None

    print(
        f"Selected unknown IR CSV: "
        f"{mixture_ir_file}"
    )

    # ========================================================
    # SELECT UNKNOWN FFT CSV
    # ========================================================

    print()
    print(
        "=" * 90
    )

    print(
        "SELECT UNKNOWN FFT CSV"
    )

    print(
        "=" * 90
    )

    mixture_fft_file = (
        select_fft_csv_file()
    )

    if mixture_fft_file is None:

        print(
            "❌ No unknown FFT CSV selected."
        )

        return None

    print(
        f"Selected unknown FFT CSV: "
        f"{mixture_fft_file}"
    )

    # ========================================================
    # LOAD UNKNOWN IR
    # ========================================================

    try:

        ir_x, ir_y, ir_type = (
            load_csv_file(
                mixture_ir_file,
                "IR"
            )
        )

    except Exception as error:

        print()
        print(
            "❌ Could not load unknown "
            f"IR CSV: {error}"
        )

        return None

    # ========================================================
    # LOAD UNKNOWN FFT
    # ========================================================

    try:

        fft_x, fft_y, fft_type = (
            load_csv_file(
                mixture_fft_file,
                "FFT"
            )
        )

    except Exception as error:

        print()
        print(
            "❌ Could not load unknown "
            f"FFT CSV: {error}"
        )

        return None

    # ========================================================
    # UNKNOWN FILE INFORMATION
    # ========================================================

    print()
    print(
        "=" * 90
    )

    print(
        "UNKNOWN FILES LOADED"
    )

    print(
        "=" * 90
    )

    print(
        f"IR file:  {mixture_ir_file}"
    )

    print(
        f"FFT file: {mixture_fft_file}"
    )

    print()
    print(
        f"✅ IR points: {len(ir_x)}"
    )

    print(
        f"   IR X range: "
        f"{np.min(ir_x):.2f} - "
        f"{np.max(ir_x):.2f}"
    )

    print(
        f"   IR Y range: "
        f"{np.min(ir_y):.6f} - "
        f"{np.max(ir_y):.6f}"
    )

    print()
    print(
        f"✅ FFT points: {len(fft_x)}"
    )

    print(
        f"   FFT X range: "
        f"{np.min(fft_x):.8f} - "
        f"{np.max(fft_x):.8f}"
    )

    print(
        f"   FFT Y range: "
        f"{np.min(fft_y):.6f} - "
        f"{np.max(fft_y):.6f}"
    )

    sample_name = (
        mixture_ir_file.stem
        .replace(
            "_IR",
            ""
        )
        .replace(
            "_ir",
            ""
        )
    )

    unknown_absorbance = (
        convert_to_absorbance(
            ir_y,
            ir_type
        )
    )

    # ========================================================
    # SELECT DATABASE
    # ========================================================

    # ========================================================
# AUTOMATIC DATABASE
# ========================================================

    # ========================================================
    # AUTOMATIC DATABASE
    # ========================================================

    print()
    print(
        "=" * 90
    )

    print(
        "DATABASE"
    )

    print(
        "=" * 90
    )

    db_path = (
        Path.home()
        / "IR_Spectra"
        / "compound_spectra.db"
    )

    print(
        f"Using database:\n{db_path}"
    )

    if not db_path.exists():

        print()
        print(
            "❌ DATABASE NOT FOUND"
        )

        print(
            f"Expected location:\n{db_path}"
        )

        return None

    if not db_path.is_file():

        print()
        print(
            "❌ Database path is not a file:"
        )

        print(
            db_path
        )

        return None

    print()
    print(
        "✅ Database found."
    )

    # ========================================================
    # LOAD DATABASE
    # ========================================================

    records = load_database_records(
        db_path
    )

    if not records:

        print()
        print(
            "❌ No usable records were loaded "
            "from the database."
        )

        return None

    # ========================================================
    # SELECT DATABASE COMPOUNDS
    # ========================================================

    selected_records = (
        select_database_compounds(
            records
        )
    )

    if not selected_records:

        print(
            "❌ No database compounds selected."
        )

        return None

    # ========================================================
    # MAXIMUM COMPONENTS
    # ========================================================

    maximum_compounds = ask_integer(
        (
            "Maximum number of compounds "
            f"(1-10) "
            f"[default "
            f"{DEFAULT_MAX_COMPOUNDS}]: "
        ),
        1,
        10,
        DEFAULT_MAX_COMPOUNDS
    )

    maximum_compounds = min(
        maximum_compounds,
        len(selected_records)
    )

    # ========================================================
    # MINIMUM SIMILARITY
    # ========================================================

    minimum_similarity = ask_float(
        (
            "Minimum combined similarity "
            f"(0-100) "
            f"[default "
            f"{DEFAULT_MINIMUM_SIMILARITY:.0f}]: "
        ),
        0.0,
        100.0,
        DEFAULT_MINIMUM_SIMILARITY
    )

    # ========================================================
    # INITIAL DATABASE SEARCH
    #
    # IMPORTANT:
    # The user's uploaded FFT CSV is used here.
    # ========================================================

    initial_candidates = (
        compare_database(
            ir_x,
            unknown_absorbance,
            fft_x,
            fft_y,
            selected_records
        )
    )

    if not initial_candidates:

        print(
            "❌ No database spectra "
            "could be compared."
        )

        return None

    print_candidate_table(
        initial_candidates,
        "TOP INITIAL DATABASE MATCHES"
    )

    # ========================================================
    # RESIDUAL AXIS
    # ========================================================

    unknown_min = max(
        float(np.min(ir_x)),
        DEFAULT_MIN_WAVENUMBER
    )

    unknown_max = min(
        float(np.max(ir_x)),
        DEFAULT_MAX_WAVENUMBER
    )

    if unknown_max <= unknown_min:

        unknown_min = float(
            np.min(ir_x)
        )

        unknown_max = float(
            np.max(ir_x)
        )

    mixture_axis = np.linspace(
        unknown_min,
        unknown_max,
        COMPARISON_LENGTH
    )

    mixture_absorbance = (
        resample_spectrum(
            ir_x,
            unknown_absorbance,
            mixture_axis
        )
    )

    residual = np.maximum(
        mixture_absorbance,
        0.0
    )

    identified = []

    used_compounds = set()

    # ========================================================
    # MULTI-COMPONENT SEARCH
    # ========================================================

    for round_number in range(
        maximum_compounds
    ):

        print()
        print(
            "=" * 90
        )

        print(
            f"IDENTIFICATION ROUND "
            f"{round_number + 1}"
        )

        print(
            "=" * 90
        )

        # ----------------------------------------------------
        # For round 1, use the uploaded FFT.
        #
        # For later rounds, calculate the FFT of the residual.
        # ----------------------------------------------------

        if round_number == 0:

            residual_fft_x = fft_x.copy()
            residual_fft_y = fft_y.copy()

            print(
                "FFT source: "
                "uploaded unknown FFT CSV"
            )

        else:

            residual_fft_x, residual_fft_y = (
                calculate_fft(
                    mixture_axis,
                    residual
                )
            )

            print(
                "FFT source: "
                "calculated residual FFT"
            )

        remaining_records = [

            record

            for record in selected_records

            if record["name"]
            not in used_compounds
        ]

        candidates = compare_database(
            mixture_axis,
            residual,
            residual_fft_x,
            residual_fft_y,
            remaining_records
        )

        if not candidates:

            print(
                "No unused database "
                "compounds remain."
            )

            break

        print_candidate_table(
            candidates,
            "TOP RESIDUAL CANDIDATES"
        )

        best = candidates[0]

        # ====================================================
        # BEST NUMERICAL CANDIDATE
        # ====================================================

        print()
        print(
            "BEST NUMERICAL CANDIDATE"
        )

        print(
            "-" * 70
        )

        print(
            f"Compound: "
            f"{best['name']}"
        )

        print(
            f"Combined score: "
            f"{best['similarity']:.2f}%"
        )

        print(
            f"IR score: "
            f"{best['ir_score']:.2f}%"
        )

        print(
            f"FFT score: "
            f"{best['fft_score']:.2f}%"
        )

        print(
            f"Diagnostic score: "
            f"{best['diagnostic'] * 100:.2f}%"
        )

        print(
            f"Nitrile score: "
            f"{best['nitrile'] * 100:.2f}%"
        )

        print(
            f"Correlation: "
            f"{best['correlation'] * 100:.2f}%"
        )

        print(
            f"Cosine: "
            f"{best['cosine'] * 100:.2f}%"
        )

        # ====================================================
        # THRESHOLD TEST
        # ====================================================

        if (
            best["similarity"]
            <
            minimum_similarity
        ):

            print()
            print(
                f"⚠️ Best numerical candidate "
                f"is below the "
                f"{minimum_similarity:.2f}% "
                f"threshold."
            )

            print(
                "❌ It will NOT be identified."
            )

            break

        # ====================================================
        # RESIDUAL FIT
        # ====================================================

        (
            coefficient,
            new_residual,
            reduction
        ) = fit_candidate(
            residual,
            best,
            mixture_axis
        )

        print(
            f"Spectral coefficient: "
            f"{coefficient:.6f}"
        )

        print(
            f"Residual reduction: "
            f"{reduction * 100:.2f}%"
        )

        if coefficient <= EPSILON:

            print(
                "⚠️ Candidate has no positive "
                "spectral contribution."
            )

            break

        if reduction <= 0.0001:

            print(
                "⚠️ Candidate did not improve "
                "the residual enough."
            )

            break

        # ====================================================
        # SAVE IDENTIFIED COMPOUND
        # ====================================================

        identified.append({

            "name":
                best["name"],

            "distance":
                best["distance"],

            "rmse":
                best["rmse"],

            "correlation":
                best["correlation"],

            "cosine":
                best["cosine"],

            "derivative":
                best["derivative"],

            "diagnostic":
                best["diagnostic"],

            "nitrile":
                best["nitrile"],

            "ir_score":
                best["ir_score"],

            "fft_score":
                best["fft_score"],

            "similarity":
                best["similarity"],

            "coefficient":
                coefficient,

            "residual_reduction":
                reduction,
        })

        used_compounds.add(
            best["name"]
        )

        residual = np.maximum(
            new_residual,
            0.0
        )

        residual_energy = float(
            np.mean(
                (
                    residual
                    -
                    np.mean(residual)
                )
                ** 2
            )
        )

        print(
            f"Remaining residual energy: "
            f"{residual_energy:.8f}"
        )

        if residual_energy < 0.0001:

            print(
                "✅ Remaining residual "
                "is very small."
            )

            break

    # ========================================================
    # NO RESULTS
    # ========================================================

    if not identified:

        print()
        print(
            "=" * 90
        )

        print(
            "NO COMPOUNDS IDENTIFIED"
        )

        print(
            "=" * 90
        )

        print()
        print(
            "The highest numerical candidate "
            "did not meet the required "
            "similarity threshold."
        )

        return None

    # ========================================================
    # PUBCHEM
    # ========================================================

    for compound in identified:

        compound["pubchem"] = (
            get_pubchem_information(
                compound["name"]
            )
        )

    # ========================================================
    # FINAL RESULTS
    # ========================================================

    print()
    print(
        "=" * 90
    )

    print(
        "FINAL MIXTURE RESULTS"
    )

    print(
        "=" * 90
    )

    summary_rows = []

    for position, compound in enumerate(
        identified,
        start=1
    ):

        pubchem = (
            compound.get(
                "pubchem"
            )
            or {}
        )

        print()
        print(
            f"{position}. "
            f"{compound['name']}"
        )

        print(
            f"   Combined score: "
            f"{compound['similarity']:.2f}%"
        )

        print(
            f"   IR score: "
            f"{compound['ir_score']:.2f}%"
        )

        print(
            f"   FFT score: "
            f"{compound['fft_score']:.2f}%"
        )

        print(
            f"   Diagnostic score: "
            f"{compound['diagnostic'] * 100:.2f}%"
        )

        print(
            f"   Nitrile score: "
            f"{compound['nitrile'] * 100:.2f}%"
        )

        print(
            f"   Correlation: "
            f"{compound['correlation'] * 100:.2f}%"
        )

        print(
            f"   Cosine: "
            f"{compound['cosine'] * 100:.2f}%"
        )

        print(
            f"   Spectral coefficient: "
            f"{compound['coefficient']:.6f}"
        )

        print(
            f"   Residual reduction: "
            f"{compound['residual_reduction'] * 100:.2f}%"
        )

        if pubchem:

            print(
                f"   PubChem CID: "
                f"{pubchem.get('cid', 'N/A')}"
            )

            print(
                f"   Formula: "
                f"{pubchem.get('formula', 'N/A')}"
            )

            print(
                f"   SMILES: "
                f"{pubchem.get('smiles', 'N/A')}"
            )

        summary_rows.append({

            "Rank":
                position,

            "Compound":
                compound["name"],

            "Combined Similarity (%)":
                round(
                    compound["similarity"],
                    2
                ),

            "IR Score (%)":
                round(
                    compound["ir_score"],
                    2
                ),

            "FFT Score (%)":
                round(
                    compound["fft_score"],
                    2
                ),

            "Diagnostic Score (%)":
                round(
                    compound["diagnostic"]
                    * 100,
                    2
                ),

            "Nitrile Score (%)":
                round(
                    compound["nitrile"]
                    * 100,
                    2
                ),

            "Correlation (%)":
                round(
                    compound["correlation"]
                    * 100,
                    2
                ),

            "Cosine Similarity (%)":
                round(
                    compound["cosine"]
                    * 100,
                    2
                ),

            "Derivative Similarity (%)":
                round(
                    compound["derivative"]
                    * 100,
                    2
                ),

            "MSE":
                compound["distance"],

            "RMSE":
                compound["rmse"],

            "Spectral Coefficient":
                compound["coefficient"],

            "Residual Reduction (%)":
                round(
                    compound[
                        "residual_reduction"
                    ]
                    * 100,
                    2
                ),

            "PubChem CID":
                pubchem.get(
                    "cid",
                    "N/A"
                ),

            "Molecular Formula":
                pubchem.get(
                    "formula",
                    "N/A"
                ),

            "SMILES":
                pubchem.get(
                    "smiles",
                    "N/A"
                ),
        })

    # ========================================================
    # SAVE RESULTS
    # ========================================================

    summary_df = pd.DataFrame(
        summary_rows
    )

    results_path = (
        mixture_ir_file.parent
        /
        f"{sample_name}_mixture_results.csv"
    )

    try:

        summary_df.to_csv(
            results_path,
            index=False
        )

        print()
        print(
            f"✅ Results saved: "
            f"{results_path}"
        )

    except Exception as error:

        print(
            f"⚠️ Could not save results: "
            f"{error}"
        )

    print()
    print(
        summary_df.to_string(
            index=False
        )
    )

    # ========================================================
    # STRUCTURES
    # ========================================================

    print()
    print(
        "Opening identified chemical structures..."
    )

    try:

        show_structure(
            identified
        )

    except Exception as error:

        print(
            f"⚠️ Structure window error: "
            f"{error}"
        )

    # ========================================================
    # FINAL MESSAGE
    # ========================================================

    print()
    print(
        "=" * 90
    )

    print(
        "IMPORTANT"
    )

    print(
        "=" * 90
    )

    print(
        "Similarity is a mathematical "
        "spectral score, not a probability."
    )

    print(
        "The spectral coefficient is not "
        "automatically concentration."
    )

    print(
        "Quantitative mixture composition "
        "requires calibrated reference mixtures."
    )

    print(
        "Unknown input:"
    )

    print(
        "  ✓ IR CSV"
    )

    print(
        "  ✓ FFT CSV"
    )

    print(
        "Database input:"
    )

    print(
        "  ✓ Database IR CSV"
    )

    print(
        "  ✓ Database FFT CSV"
    )

    print(
        "JCAMP-DX files are NOT used "
        "by the identification program."
    )

    print(
        "=" * 90
    )

    return identified


# ============================================================
# PUBLIC ENTRY POINT
# ============================================================

def Multiple():

    try:

        return identify_mixture_compounds()

    except KeyboardInterrupt:

        print(
            "\nProgram stopped by user."
        )

        return None

    except Exception as error:

        print()
        print(
            "=" * 90
        )

        print(
            "UNEXPECTED PROGRAM ERROR"
        )

        print(
            "=" * 90
        )

        print(
            f"{type(error).__name__}: "
            f"{error}"
        )

        import traceback

        traceback.print_exc()

        return None


# ============================================================
# START PROGRAM
# ============================================================

if __name__ == "__main__":

    Multiple()