
try:

    from databasis import (
        save_spectra_to_db,
        DB_PATH
    )

except ImportError as error:

    print()
    print("=" * 80)
    print("❌ DATABASE MODULE ERROR")
    print("=" * 80)
    print(
        "Could not import databasis.py"
    )
    print(
        f"Error: {error}"
    )
    print("=" * 80)

    save_spectra_to_db = None

    DB_PATH = None


# ============================================================
# STANDARD LIBRARY
# ============================================================

import re
import shutil
import platform

from pathlib import Path


# ============================================================
# NUMPY / PANDAS
# ============================================================

import numpy as np
import pandas as pd


# ============================================================
# MATPLOTLIB
# ============================================================

import matplotlib.pyplot as plt


# ============================================================
# TKINTER
# ============================================================

try:

    from tkinter import (
        Tk,
        filedialog,
        messagebox
    )

    TKINTER_AVAILABLE = True

except ImportError:

    TKINTER_AVAILABLE = False


# ============================================================
# OPTIONAL JCAMP PACKAGE
# ============================================================

try:

    import jcamp

    JCAMP_AVAILABLE = True

except ImportError:

    jcamp = None

    JCAMP_AVAILABLE = False


# ============================================================
# APPLICATION DIRECTORIES
# ============================================================

APPLICATION_FOLDER = (
    Path.home()
    /
    "IR_Spectra"
)

SPECTRA_FOLDER = (
    APPLICATION_FOLDER
    /
    "spectra"
)

CSV_FOLDER = (
    APPLICATION_FOLDER
    /
    "csv"
)

FFT_FOLDER = (
    APPLICATION_FOLDER
    /
    "fft"
)

PLOT_FOLDER = (
    APPLICATION_FOLDER
    /
    "plots"
)


# ============================================================
# CREATE APPLICATION DIRECTORIES
# ============================================================

def create_application_folders():

    folders = [

        APPLICATION_FOLDER,

        SPECTRA_FOLDER,

        CSV_FOLDER,

        FFT_FOLDER,

        PLOT_FOLDER

    ]

    for folder in folders:

        try:

            folder.mkdir(
                parents=True,
                exist_ok=True
            )

        except Exception as error:

            raise RuntimeError(
                f"Could not create folder "
                f"{folder}: {error}"
            )


# ============================================================
# SAFE FILE NAME
# ============================================================

def create_safe_filename(
    name
):

    name = str(
        name
    ).strip()

    if not name:

        name = "spectrum"

    # Replace invalid characters.

    name = re.sub(
        r'[<>:"/\\|?*]',
        "_",
        name
    )

    # Replace whitespace.

    name = re.sub(
        r"\s+",
        "_",
        name
    )

    # Remove repeated underscores.

    name = re.sub(
        r"_+",
        "_",
        name
    )

    # Remove leading/trailing dots.

    name = name.strip(
        ". "
    )

    if not name:

        name = "spectrum"

    return name


# ============================================================
# SELECT FILE
# ============================================================

def select_jdx_file():

    if not TKINTER_AVAILABLE:

        print()
        print("=" * 80)
        print("❌ TKINTER IS NOT AVAILABLE")
        print("=" * 80)

        print(
            "Python cannot open the file-selection "
            "window because Tkinter is unavailable."
        )

        return None

    root = None

    try:

        root = Tk()

        root.withdraw()

        try:

            root.attributes(
                "-topmost",
                True
            )

        except Exception:

            pass

        selected_file = (
            filedialog.askopenfilename(
                title=(
                    "Select JCAMP-DX IR Spectrum"
                ),
                filetypes=[

                    (
                        "JCAMP-DX files",
                        "*.jdx"
                    ),

                    (
                        "JCAMP files",
                        "*.dx"
                    ),

                    (
                        "JDX/DX files",
                        "*.jdx *.dx"
                    ),

                    (
                        "All files",
                        "*.*"
                    )

                ]
            )
        )

        return (
            Path(
                selected_file
            ).resolve()
            if selected_file
            else None
        )

    except Exception as error:

        print()
        print(
            f"❌ Could not open "
            f"file-selection window: "
            f"{error}"
        )

        return None

    finally:

        if root is not None:

            try:

                root.destroy()

            except Exception:

                pass


# ============================================================
# READ JCAMP HEADER VALUE
# ============================================================

def get_jcamp_header_value(
    lines,
    key
):

    key_upper = (
        key.upper()
    )

    prefix = (
        f"##{key_upper}="
    )

    for line in lines:

        stripped = (
            line.strip()
        )

        if stripped.upper().startswith(
            prefix
        ):

            return (
                stripped.split(
                    "=",
                    1
                )[1]
                .strip()
            )

    return None


# ============================================================
# FLOAT HEADER
# ============================================================

def get_float_header(
    lines,
    key,
    default=None
):

    value = get_jcamp_header_value(
        lines,
        key
    )

    if value is None:

        return default

    try:

        return float(
            value
        )

    except (
        ValueError,
        TypeError
    ):

        return default


# ============================================================
# INTEGER HEADER
# ============================================================

def get_int_header(
    lines,
    key,
    default=None
):

    value = get_jcamp_header_value(
        lines,
        key
    )

    if value is None:

        return default

    try:

        return int(
            float(value)
        )

    except (
        ValueError,
        TypeError
    ):

        return default


# ============================================================
# FIND XYDATA
# ============================================================

def find_xydata_index(
    lines
):

    for index, line in enumerate(
        lines
    ):

        if line.strip().upper().startswith(
            "##XYDATA="
        ):

            return index

    return None


# ============================================================
# DETERMINE Y UNITS
# ============================================================

def determine_y_units(
    lines
):

    value = get_jcamp_header_value(
        lines,
        "YUNITS"
    )

    if value is None:

        return ""

    return str(
        value
    ).strip()


# ============================================================
# DETERMINE SPECTRUM TYPE
# ============================================================

def determine_spectrum_type(
    y_units
):

    units = (
        str(
            y_units
        )
        .strip()
        .upper()
    )

    if "ABSORBANCE" in units:

        return "ABSORBANCE"

    if "TRANSMITTANCE" in units:

        return "TRANSMITTANCE"

    if "TRANSMISSION" in units:

        return "TRANSMITTANCE"

    return "UNKNOWN"


# ============================================================
# DECODE X++(Y..Y)
# ============================================================

def decode_xydata(
    lines
):

    xydata_index = find_xydata_index(
        lines
    )

    if xydata_index is None:

        raise ValueError(
            "JCAMP-DX file does not contain "
            "##XYDATA=."
        )

    xydata_header = (
        lines[
            xydata_index
        ]
        .strip()
    )

    print()
    print(
        f"XYDATA format: "
        f"{xydata_header}"
    )

    # --------------------------------------------------------
    # JCAMP scaling factors
    # --------------------------------------------------------

    x_factor = get_float_header(
        lines,
        "XFACTOR",
        1.0
    )

    y_factor = get_float_header(
        lines,
        "YFACTOR",
        1.0
    )

    delta_x = get_float_header(
        lines,
        "DELTAX",
        None
    )

    first_x = get_float_header(
        lines,
        "FIRSTX",
        None
    )

    last_x = get_float_header(
        lines,
        "LASTX",
        None
    )

    n_points = get_int_header(
        lines,
        "NPOINTS",
        None
    )

    # --------------------------------------------------------
    # Calculate DELTAX if missing
    # --------------------------------------------------------

    if delta_x is None:

        if (
            first_x is not None
            and
            last_x is not None
            and
            n_points is not None
            and
            n_points > 1
        ):

            delta_x = (
                last_x
                -
                first_x
            ) / (
                n_points - 1
            )

        else:

            delta_x = 1.0

    print()
    print(
        "JCAMP INFORMATION"
    )

    print(
        f"FIRSTX  : {first_x}"
    )

    print(
        f"LASTX   : {last_x}"
    )

    print(
        f"DELTAX  : {delta_x}"
    )

    print(
        f"NPOINTS : {n_points}"
    )

    print(
        f"XFACTOR : {x_factor}"
    )

    print(
        f"YFACTOR : {y_factor}"
    )

    # --------------------------------------------------------
    # Arrays
    # --------------------------------------------------------

    x_values = []
    y_values = []

    # --------------------------------------------------------
    # Read data
    # --------------------------------------------------------

    for line in lines[
        xydata_index + 1:
    ]:

        stripped = (
            line.strip()
        )

        if not stripped:

            continue

        upper = (
            stripped.upper()
        )

        # End of JCAMP data.

        if upper.startswith(
            "##END="
        ):

            break

        # Ignore metadata.

        if stripped.startswith(
            "##"
        ):

            continue

        parts = (
            stripped.split()
        )

        if len(parts) < 2:

            continue

        # ----------------------------------------------------
        # Starting X
        # ----------------------------------------------------

        try:

            start_x = float(
                parts[0]
            )

        except ValueError:

            continue

        # ----------------------------------------------------
        # Y values
        # ----------------------------------------------------

        row_y_values = []

        for value in parts[1:]:

            try:

                raw_y = float(
                    value
                )

                row_y_values.append(
                    raw_y
                )

            except ValueError:

                continue

        # ----------------------------------------------------
        # Generate X values
        # ----------------------------------------------------

        for index, raw_y in enumerate(
            row_y_values
        ):

            x = (
                start_x
                +
                index * delta_x
            )

            x = (
                x
                *
                x_factor
            )

            y = (
                raw_y
                *
                y_factor
            )

            x_values.append(
                x
            )

            y_values.append(
                y
            )

    if not x_values:

        raise ValueError(
            "No spectral data could be decoded "
            "from ##XYDATA=."
        )

    x_array = np.asarray(
        x_values,
        dtype=float
    )

    y_array = np.asarray(
        y_values,
        dtype=float
    )

    # --------------------------------------------------------
    # Respect NPOINTS if the decoder produced extras
    # --------------------------------------------------------

    if (
        n_points is not None
        and
        n_points > 0
        and
        len(x_array) > n_points
    ):

        x_array = (
            x_array[
                :n_points
            ]
        )

        y_array = (
            y_array[
                :n_points
            ]
        )

    print()
    print(
        f"Decoded points: "
        f"{len(x_array)}"
    )

    if n_points is not None:

        if len(x_array) == n_points:

            print(
                "✅ NPOINTS verified."
            )

        else:

            print(
                f"⚠️ Expected "
                f"{n_points} points, "
                f"decoded "
                f"{len(x_array)}."
            )

    return (
        x_array,
        y_array
    )


# ============================================================
# DECODE USING JCAMP PACKAGE
# ============================================================

def decode_using_jcamp(
    file_path
):

    if not JCAMP_AVAILABLE:

        raise ImportError(
            "The jcamp package is not installed."
        )

    data = jcamp.readfile(
        str(file_path)
    )

    x = np.asarray(
        data.get(
            "x",
            []
        ),
        dtype=float
    )

    y = np.asarray(
        data.get(
            "y",
            []
        ),
        dtype=float
    )

    if (
        len(x) == 0
        or
        len(y) == 0
    ):

        raise ValueError(
            "The jcamp package returned "
            "no spectral data."
        )

    if len(x) != len(y):

        raise ValueError(
            "JCamp X and Y arrays have "
            "different lengths."
        )

    return (
        x,
        y
    )


# ============================================================
# CONVERT TO TRANSMITTANCE
# ============================================================

def convert_to_transmittance(
    y_values,
    y_units
):

    values = np.asarray(
        y_values,
        dtype=float
    )

    units = (
        str(
            y_units
        )
        .strip()
        .upper()
    )

    # --------------------------------------------------------
    # Absorbance
    #
    # T = 10^(-A)
    # --------------------------------------------------------

    if "ABSORBANCE" in units:

        transmittance = (
            10.0
            **
            (
                -values
            )
        )

    # --------------------------------------------------------
    # Percentage transmittance
    # --------------------------------------------------------

    elif (
        "TRANSMITTANCE" in units
        and
        "%"
        in units
    ):

        transmittance = (
            values
            /
            100.0
        )

    # --------------------------------------------------------
    # Fractional transmittance
    # --------------------------------------------------------

    elif (
        "TRANSMITTANCE"
        in units
    ):

        transmittance = (
            values.copy()
        )

        # If values appear to be percentages despite
        # the header not containing %, convert them.

        finite_values = (
            transmittance[
                np.isfinite(
                    transmittance
                )
            ]
        )

        if (
            len(finite_values) > 0
            and
            np.nanmax(
                finite_values
            ) > 1.0
        ):

            transmittance /= 100.0

    # --------------------------------------------------------
    # Unknown
    #
    # Assume values are already transmittance.
    # --------------------------------------------------------

    else:

        print()
        print(
            "⚠️ YUNITS is unknown."
        )

        print(
            "Assuming Y values are "
            "transmittance."
        )

        transmittance = (
            values.copy()
        )

        finite_values = (
            transmittance[
                np.isfinite(
                    transmittance
                )
            ]
        )

        if (
            len(finite_values) > 0
            and
            np.nanmax(
                finite_values
            ) > 1.0
        ):

            print(
                "⚠️ Values exceed 1.0."
            )

            print(
                "Assuming they are "
                "percentage transmittance."
            )

            transmittance /= 100.0

    return np.asarray(
        transmittance,
        dtype=float
    )


# ============================================================
# CLEAN SPECTRUM
# ============================================================

def clean_spectrum(
    x,
    y
):

    x = np.asarray(
        x,
        dtype=float
    )

    y = np.asarray(
        y,
        dtype=float
    )

    if len(x) != len(y):

        raise ValueError(
            "X and Y arrays have "
            "different lengths."
        )

    # --------------------------------------------------------
    # Remove invalid values
    # --------------------------------------------------------

    valid = (
        np.isfinite(x)
        &
        np.isfinite(y)
    )

    x = x[
        valid
    ]

    y = y[
        valid
    ]

    if len(x) < 2:

        raise ValueError(
            "Fewer than two valid "
            "spectral points remain."
        )

    # --------------------------------------------------------
    # Create DataFrame
    # --------------------------------------------------------

    df = pd.DataFrame({

        "Wavenumber (cm⁻¹)":
            x,

        "Transmittance":
            y

    })

    # --------------------------------------------------------
    # Remove duplicates
    # --------------------------------------------------------

    df = df.drop_duplicates(
        subset=[
            "Wavenumber (cm⁻¹)"
        ]
    )

    # --------------------------------------------------------
    # Sort by wavenumber
    # --------------------------------------------------------

    df = df.sort_values(
        by="Wavenumber (cm⁻¹)"
    )

    df.reset_index(
        drop=True,
        inplace=True
    )

    # --------------------------------------------------------
    # Clean transmittance
    # --------------------------------------------------------

    df[
        "Transmittance"
    ] = pd.to_numeric(
        df[
            "Transmittance"
        ],
        errors="coerce"
    )

    df = df.dropna()

    # --------------------------------------------------------
    # Clip to physical transmittance range
    # --------------------------------------------------------

    df[
        "Transmittance"
    ] = np.clip(
        df[
            "Transmittance"
        ].to_numpy(
            dtype=float
        ),
        1e-6,
        1.0
    )

    # --------------------------------------------------------
    # DO NOT round wavenumbers here.
    #
    # Rounding before FFT can create duplicate X values and
    # can alter the spacing of the spectrum.
    # --------------------------------------------------------

    return df


# ============================================================
# RESAMPLE TO UNIFORM GRID
# ============================================================
#
# This is important for FFT.
#
# FFT assumes uniformly spaced samples.
#
# An IR spectrum may contain slightly irregular X positions.
#
# Therefore we interpolate the spectrum onto a uniform
# wavenumber grid before calculating the FFT.
#
# ============================================================

def resample_uniform_grid(
    df
):

    x = (
        df[
            "Wavenumber (cm⁻¹)"
        ]
        .to_numpy(
            dtype=float
        )
    )

    y = (
        df[
            "Transmittance"
        ]
        .to_numpy(
            dtype=float
        )
    )

    if len(x) < 2:

        raise ValueError(
            "Not enough points for "
            "uniform resampling."
        )

    # --------------------------------------------------------
    # Determine number of points
    # --------------------------------------------------------

    target_length = len(x)

    # --------------------------------------------------------
    # Uniform grid
    # --------------------------------------------------------

    uniform_x = np.linspace(
        x.min(),
        x.max(),
        target_length
    )

    # --------------------------------------------------------
    # Interpolate
    # --------------------------------------------------------

    uniform_y = np.interp(
        uniform_x,
        x,
        y
    )

    return (
        uniform_x,
        uniform_y
    )


# ============================================================
# CREATE FFT
# ============================================================

def create_fft(
    intensity
):

    intensity = np.asarray(
        intensity,
        dtype=float
    )

    if len(intensity) < 2:

        raise ValueError(
            "At least two points are "
            "required for FFT."
        )

    # --------------------------------------------------------
    # Replace invalid values
    # --------------------------------------------------------

    intensity = np.nan_to_num(
        intensity,
        nan=0.0,
        posinf=0.0,
        neginf=0.0
    )

    n = len(
        intensity
    )

    # ========================================================
    # REMOVE DC COMPONENT
    # ========================================================
    #
    # The mean of the signal produces a very large zero-
    # frequency component.
    #
    # For spectral-shape comparison we do not want the overall
    # baseline to dominate the FFT.
    #
    # ========================================================

    signal = (
        intensity
        -
        np.mean(
            intensity
        )
    )

    # ========================================================
    # ONE-SIDED REAL FFT
    # ========================================================

    fft_values = np.fft.rfft(
        signal
    )

    fft_magnitude = (
        np.abs(
            fft_values
        )
        /
        n
    )

    # ========================================================
    # AMPLITUDE CORRECTION
    # ========================================================
    #
    # For a real signal, positive-frequency components are
    # doubled except:
    #
    #   - DC component
    #   - Nyquist component for even N
    #
    # rfft() already returns only the positive-frequency side.
    #
    # ========================================================

    if n % 2 == 0:

        # Even N:
        #
        # index 0 = DC
        # last index = Nyquist

        if len(fft_magnitude) > 2:

            fft_magnitude[
                1:-1
            ] *= 2.0

    else:

        # Odd N:
        #
        # No Nyquist point.

        if len(fft_magnitude) > 1:

            fft_magnitude[
                1:
            ] *= 2.0

    # ========================================================
    # FORCE DC VALUE TO ZERO
    # ========================================================
    #
    # The signal was mean-centered, so DC should already be
    # essentially zero.
    #
    # Explicitly setting it to zero makes the output clearer
    # and prevents tiny numerical values from affecting
    # comparisons.
    #
    # ========================================================

    fft_magnitude[0] = 0.0

    # ========================================================
    # FREQUENCY INDEX
    # ========================================================

    frequency_index = np.arange(
        len(
            fft_magnitude
        ),
        dtype=int
    )

    # ========================================================
    # DATAFRAME
    # ========================================================

    fft_df = pd.DataFrame({

        "Frequency Index":
            frequency_index,

        "FFT Magnitude":
            fft_magnitude

    })

    return fft_df


# ============================================================
# PROCESS JCAMP-DX
# ============================================================

def process_jdx_file():

    print()
    print("=" * 80)
    print(
        "        JCAMP-DX IR + FOURIER TRANSFORM PROCESSOR"
    )
    print("=" * 80)

    print()
    print(
        "Operating system:"
    )

    print(
        platform.system()
    )

    # ========================================================
    # CREATE DIRECTORIES
    # ========================================================

    try:

        create_application_folders()

    except Exception as error:

        print()
        print(
            f"❌ Could not create application folders:"
        )

        print(
            error
        )

        return None

    # ========================================================
    # SELECT FILE
    # ========================================================

    print()
    print(
        "Select a JCAMP-DX IR spectrum."
    )

    file_path = select_jdx_file()

    if file_path is None:

        print()
        print(
            "⚠️ No JDX file selected."
        )

        return None

    # ========================================================
    # VERIFY FILE
    # ========================================================

    if not file_path.exists():

        print()
        print(
            "❌ Selected file does not exist:"
        )

        print(
            file_path
        )

        return None

    if not file_path.is_file():

        print()
        print(
            "❌ Selected path is not a file:"
        )

        print(
            file_path
        )

        return None

    # ========================================================
    # DISPLAY FILE
    # ========================================================

    print()
    print("=" * 80)
    print(
        "FILE SELECTED"
    )
    print("=" * 80)

    print(
        f"File: {file_path.name}"
    )

    print(
        f"Location: {file_path.parent}"
    )

    # ========================================================
    # COPY JDX TO APPLICATION FOLDER
    # ========================================================

    saved_jdx_path = (
        SPECTRA_FOLDER
        /
        file_path.name
    )

    try:

        if (
            file_path.resolve()
            !=
            saved_jdx_path.resolve()
        ):

            shutil.copy2(
                file_path,
                saved_jdx_path
            )

            print()
            print(
                "✅ JDX file copied to:"
            )

            print(
                saved_jdx_path
            )

        else:

            saved_jdx_path = file_path

    except Exception as error:

        print()
        print(
            f"⚠️ Could not copy JDX file:"
        )

        print(
            error
        )

        print(
            "Continuing with original file."
        )

        saved_jdx_path = file_path

    # ========================================================
    # READ FILE
    # ========================================================

    try:

        with open(
            saved_jdx_path,
            "r",
            encoding="utf-8",
            errors="replace"
        ) as file:

            lines = file.readlines()

    except Exception as error:

        print()
        print(
            "❌ Failed to read JCAMP-DX file:"
        )

        print(
            error
        )

        return None

    if not lines:

        print()
        print(
            "❌ The JDX file is empty."
        )

        return None

    # ========================================================
    # COMPOUND NAME
    # ========================================================

    print()

    compound_name = input(
        "Enter compound name "
        "(press Enter to use filename): "
    ).strip()

    if not compound_name:

        compound_name = (
            file_path.stem
        )

    safe_compound_name = (
        create_safe_filename(
            compound_name
        )
    )

    # ========================================================
    # Y UNITS
    # ========================================================

    y_units = determine_y_units(
        lines
    )

    spectrum_type = (
        determine_spectrum_type(
            y_units
        )
    )

    print()
    print(
        f"Y units: {y_units}"
    )

    print(
        f"Detected spectrum type: "
        f"{spectrum_type}"
    )

    # ========================================================
    # DECODE SPECTRUM
    # ========================================================

    decoded = None

    # --------------------------------------------------------
    # First attempt:
    #
    # Native JCAMP decoder.
    #
    # --------------------------------------------------------

    try:

        decoded = decode_xydata(
            lines
        )

        print()
        print(
            "✅ Spectrum decoded using "
            "native JCAMP-DX decoder."
        )

    except Exception as native_error:

        print()
        print(
            "⚠️ Native JCAMP-DX decoder failed:"
        )

        print(
            native_error
        )

        # ----------------------------------------------------
        # Second attempt:
        #
        # jcamp package.
        # ----------------------------------------------------

        if JCAMP_AVAILABLE:

            try:

                decoded = (
                    decode_using_jcamp(
                        saved_jdx_path
                    )
                )

                print()
                print(
                    "✅ Spectrum decoded using "
                    "jcamp package."
                )

            except Exception as jcamp_error:

                print()
                print(
                    "❌ jcamp decoder also failed:"
                )

                print(
                    jcamp_error
                )

        if decoded is None:

            print()
            print(
                "❌ Could not decode the "
                "JCAMP-DX spectrum."
            )

            return None

    # ========================================================
    # CONVERT DATA
    # ========================================================

    wavenumbers = (
        decoded[0]
    )

    raw_values = (
        decoded[1]
    )

    # ========================================================
    # CONVERT TO TRANSMITTANCE
    # ========================================================

    try:

        transmittance = (
            convert_to_transmittance(
                raw_values,
                y_units
            )
        )

    except Exception as error:

        print()
        print(
            "❌ Failed to convert spectrum "
            "to transmittance:"
        )

        print(
            error
        )

        return None

    # ========================================================
    # CLEAN
    # ========================================================

    try:

        df = clean_spectrum(
            wavenumbers,
            transmittance
        )

    except Exception as error:

        print()
        print(
            "❌ Failed to clean spectrum:"
        )

        print(
            error
        )

        return None

    # ========================================================
    # VERIFY DATA
    # ========================================================

    if len(df) < 10:

        print()
        print(
            "❌ Fewer than 10 usable spectral points "
            "remain after cleaning."
        )

        return None

    # ========================================================
    # DISPLAY SUMMARY
    # ========================================================

    print()
    print("=" * 80)
    print(
        "SPECTRUM SUMMARY"
    )
    print("=" * 80)

    print(
        f"Compound: "
        f"{compound_name}"
    )

    print(
        f"Data points: "
        f"{len(df)}"
    )

    print(
        f"Minimum wavenumber: "
        f"{df['Wavenumber (cm⁻¹)'].min():.3f}"
    )

    print(
        f"Maximum wavenumber: "
        f"{df['Wavenumber (cm⁻¹)'].max():.3f}"
    )

    print(
        f"Minimum transmittance: "
        f"{df['Transmittance'].min():.6f}"
    )

    print(
        f"Maximum transmittance: "
        f"{df['Transmittance'].max():.6f}"
    )

    # ========================================================
    # SAVE IR CSV
    # ========================================================

    csv_path = (
        CSV_FOLDER
        /
        f"{safe_compound_name}_ir.csv"
    )

    try:

        df.to_csv(
            csv_path,
            index=False
        )

        print()
        print(
            "✅ IR CSV saved:"
        )

        print(
            csv_path
        )

    except Exception as error:

        print()
        print(
            "❌ Failed to save IR CSV:"
        )

        print(
            error
        )

        return None

    # ========================================================
    # UNIFORM RESAMPLING FOR FFT
    # ========================================================

    try:

        (
            uniform_wavenumber,
            uniform_transmittance
        ) = resample_uniform_grid(
            df
        )

    except Exception as error:

        print()
        print(
            "❌ Failed to create "
            "uniform FFT grid:"
        )

        print(
            error
        )

        return None

    # ========================================================
    # CREATE FFT
    # ========================================================

    try:

        fft_df = create_fft(
            uniform_transmittance
        )

    except Exception as error:

        print()
        print(
            "❌ Failed to calculate FFT:"
        )

        print(
            error
        )

        return None

    # ========================================================
    # SAVE FFT CSV
    # ========================================================

    fft_csv_path = (
        FFT_FOLDER
        /
        f"{safe_compound_name}_fft.csv"
    )

    try:

        fft_df.to_csv(
            fft_csv_path,
            index=False
        )

        print()
        print(
            "✅ FFT CSV saved:"
        )

        print(
            fft_csv_path
        )

    except Exception as error:

        print()
        print(
            "❌ Failed to save FFT CSV:"
        )

        print(
            error
        )

        return None

    # ========================================================
    # DATABASE
    # ========================================================

    # Prefer databasecode.DB_PATH when available.
    #
    # This prevents Fourier.py from silently creating a
    # different database from the one used elsewhere.
    #
    if DB_PATH is not None:

        db_path = Path(
            DB_PATH
        ).resolve()

    else:

        db_path = (
            APPLICATION_FOLDER
            /
            "compound_spectra.db"
        )

    if save_spectra_to_db is not None:

        try:

            save_spectra_to_db(
                db_path,
                compound_name,
                df,
                fft_df
            )

            print()
            print(
                "✅ Spectrum saved to database:"
            )

            print(
                db_path
            )

        except Exception as error:

            print()
            print(
                "⚠️ Database save failed:"
            )

            print(
                error
            )

    else:

        print()
        print(
            "⚠️ Database save function is unavailable."
        )

    # ========================================================
    # PLOT
    # ========================================================

    plot_path = None

    try:

        fig, axes = plt.subplots(
            2,
            1,
            figsize=(
                12,
                9
            )
        )

        # ====================================================
        # IR PLOT
        # ====================================================

        ax1 = axes[0]

        ax1.plot(
            df[
                "Wavenumber (cm⁻¹)"
            ],
            df[
                "Transmittance"
            ],
            linewidth=1.2
        )

        # Conventional IR display:
        # high wavenumber on left.

        ax1.invert_xaxis()

        ax1.set_xlabel(
            "Wavenumber (cm⁻¹)"
        )

        ax1.set_ylabel(
            "Transmittance"
        )

        ax1.set_title(
            f"IR Spectrum - "
            f"{compound_name}"
        )

        ax1.grid(
            True,
            alpha=0.3
        )

        # ====================================================
        # FFT PLOT
        # ====================================================

        ax2 = axes[1]

        frequency_index = (
            fft_df[
                "Frequency Index"
            ]
            .to_numpy(
                dtype=float
            )
        )

        fft_magnitude = (
            fft_df[
                "FFT Magnitude"
            ]
            .to_numpy(
                dtype=float
            )
        )

        # Don't plot the zero-frequency point.
        #
        # It was deliberately removed from the FFT comparison.

        if len(
            frequency_index
        ) > 1:

            ax2.vlines(
                frequency_index[1:],
                0,
                fft_magnitude[1:],
                linewidth=0.8
            )

            ax2.scatter(
                frequency_index[1:],
                fft_magnitude[1:],
                s=10
            )

        ax2.set_xlabel(
            "Frequency Index"
        )

        ax2.set_ylabel(
            "FFT Magnitude"
        )

        ax2.set_title(
            f"FFT - "
            f"{compound_name}"
        )

        ax2.grid(
            True,
            alpha=0.3
        )

        # ----------------------------------------------------
        # Automatically show a useful FFT range.
        # ----------------------------------------------------

        maximum_display_frequency = min(
            150,
            max(
                1,
                len(
                    frequency_index
                ) - 1
            )
        )

        ax2.set_xlim(
            0,
            maximum_display_frequency
        )

        plt.tight_layout()

        # ====================================================
        # SAVE PLOT
        # ====================================================

        plot_path = (
            PLOT_FOLDER
            /
            f"{safe_compound_name}_IR_FFT.png"
        )

        plt.savefig(
            plot_path,
            dpi=300,
            bbox_inches="tight"
        )

        print()
        print(
            "✅ Plot saved:"
        )

        print(
            plot_path
        )

        # ====================================================
        # SHOW PLOT
        # ====================================================

        plt.show()

        plt.close(
            fig
        )

    except Exception as error:

        print()
        print(
            "⚠️ Plotting failed:"
        )

        print(
            error
        )

        try:

            plt.close(
                "all"
            )

        except Exception:

            pass

    # ========================================================
    # FINAL SUMMARY
    # ========================================================

    print()
    print("=" * 80)
    print(
        "PROCESSING COMPLETE"
    )
    print("=" * 80)

    print(
        f"Compound:"
    )

    print(
        compound_name
    )

    print()

    print(
        "JDX file:"
    )

    print(
        saved_jdx_path
    )

    print()

    print(
        "IR CSV:"
    )

    print(
        csv_path
    )

    print()

    print(
        "FFT CSV:"
    )

    print(
        fft_csv_path
    )

    print()

    if plot_path is not None:

        print(
            "Plot:"
        )

        print(
            plot_path
        )

    print()

    print(
        "Database:"
    )

    print(
        db_path
    )

    print("=" * 80)

    # ========================================================
    # RETURN
    # ========================================================

    return {

        "compound_name":
            compound_name,

        "jdx_file":
            str(
                saved_jdx_path
            ),

        "csv_file":
            str(
                csv_path
            ),

        "fft_file":
            str(
                fft_csv_path
            ),

        "plot_file":
            (
                str(
                    plot_path
                )
                if plot_path is not None
                else None
            ),

        "database_file":
            str(
                db_path
            ),

        "wavenumbers":
            df[
                "Wavenumber (cm⁻¹)"
            ].to_numpy(
                dtype=float
            ),

        "transmittance":
            df[
                "Transmittance"
            ].to_numpy(
                dtype=float
            ),

        "uniform_wavenumbers":
            uniform_wavenumber,

        "uniform_transmittance":
            uniform_transmittance,

        "dataframe":
            df,

        "fft_dataframe":
            fft_df

    }


# ============================================================
# PUBLIC ENTRY POINT
# ============================================================

def Fourier():

    """
    Public entry point for MainMenu.py.

    MainMenu.py can call:

        Fourier.Fourier()
    """

    try:

        return process_jdx_file()

    except KeyboardInterrupt:

        print()
        print()
        print(
            "Program stopped by user."
        )

        return None

    except Exception as error:

        print()
        print("=" * 80)
        print(
            "UNEXPECTED FOURIER ERROR"
        )
        print("=" * 80)

        print(
            f"{type(error).__name__}: "
            f"{error}"
        )

        print("=" * 80)

        return None


# ============================================================
# RUN DIRECTLY
# ============================================================

if __name__ == "__main__":

    Fourier()