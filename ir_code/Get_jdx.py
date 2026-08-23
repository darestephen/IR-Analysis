import re
import platform
import time
import requests
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from pathlib import Path
from urllib.parse import quote, urljoin

from bs4 import BeautifulSoup

from databasis import save_spectra_to_db


# ============================================================
# APPLICATION SETTINGS
# ============================================================

NIST_BASE_URL = "https://webbook.nist.gov"

HEADERS = {
    "User-Agent": (
        "Mozilla/5.0 (X11; Linux x86_64) "
        "AppleWebKit/537.36 "
        "(KHTML, like Gecko) "
        "Chrome/151.0 Safari/537.36"
    ),
    "Accept": (
        "text/html,application/xhtml+xml,"
        "application/xml;q=0.9,*/*;q=0.8"
    ),
    "Accept-Language": "en-US,en;q=0.9",
    "Connection": "keep-alive",
}


# ============================================================
# APPLICATION FOLDERS
# ============================================================

APPLICATION_FOLDER = Path.home() / "project"

SPECTRA_FOLDER = APPLICATION_FOLDER / "spectra"
CSV_FOLDER = APPLICATION_FOLDER / "csv"
FFT_FOLDER = APPLICATION_FOLDER / "fft"
PLOTS_FOLDER = APPLICATION_FOLDER / "plots"

DATABASE_PATH = APPLICATION_FOLDER / "compound_spectra.db"


# ============================================================
# CREATE DIRECTORIES
# ============================================================

def create_application_folders():

    folders = [
        APPLICATION_FOLDER,
        SPECTRA_FOLDER,
        CSV_FOLDER,
        FFT_FOLDER,
        PLOTS_FOLDER,
    ]

    for folder in folders:

        try:
            folder.mkdir(
                parents=True,
                exist_ok=True
            )

        except OSError as error:

            print()
            print("❌ Could not create folder:")
            print(folder)
            print(f"Error: {error}")

            raise


create_application_folders()


# ============================================================
# DISPLAY APPLICATION INFORMATION
# ============================================================

print()
print("=" * 75)
print("IR SPECTRUM APPLICATION")
print("=" * 75)

print(
    f"Operating system: {platform.system()}"
)

print(
    f"Application folder: {APPLICATION_FOLDER}"
)

print(
    f"Database: {DATABASE_PATH}"
)

print("=" * 75)


# ============================================================
# SAFE FILE NAME
# ============================================================

def make_safe_filename(compound_name):

    safe_name = str(
        compound_name
    ).strip()

    safe_name = re.sub(
        r'[<>:"/\\|?*]',
        "_",
        safe_name
    )

    safe_name = re.sub(
        r"\s+",
        "_",
        safe_name
    )

    safe_name = safe_name.strip(
        " ._"
    )

    if not safe_name:
        safe_name = "unknown_compound"

    return safe_name


def normalize_name(name):

    """
    Normalize a compound name for comparison.
    """

    name = str(
        name
    ).strip().lower()

    name = re.sub(
        r"\s+",
        " ",
        name
    )

    return name



def is_isotopologue(name):

    """
    Return True when a NIST result appears to be
    an isotope-labelled compound.

    Examples:

        [2H6]ethanol
        ethanol-d1
        ethanol-d5
        C2H3D3O
        ethanol-1,1-d2
    """

    name = str(
        name
    ).strip().lower()

    isotope_patterns = [

        r"\[[0-9]+h[0-9]*\]",

        r"\bd[0-9]+\b",

        r"-d[0-9]+",

        r"\bd[0-9]+$",

        r"\d+d\d+",

        r"\bd\d+[a-z]*",

    ]

    for pattern in isotope_patterns:

        if re.search(
            pattern,
            name,
            re.IGNORECASE
        ):

            return True

    # Explicit deuterium in molecular formula

    if "d" in name:

        # Avoid treating normal words containing d as isotope names
        if re.search(
            r"[a-z]\d*[d]\d*",
            name,
            re.IGNORECASE
        ):

            return True

    return False


def get_page(
    url,
    timeout=30,
    retries=3
):

    """
    Safely download a NIST page.

    Includes retry handling for HTTP 429.
    """

    for attempt in range(
        1,
        retries + 1
    ):

        try:

            response = requests.get(
                url,
                headers=HEADERS,
                timeout=timeout
            )

            if response.status_code == 429:

                wait_time = (
                    3 * attempt
                )

                print()
                print(
                    f"⚠️ NIST rate limit (429). "
                    f"Waiting {wait_time} seconds..."
                )

                time.sleep(
                    wait_time
                )

                continue

            response.raise_for_status()

            return response

        except requests.exceptions.RequestException as error:

            if attempt >= retries:

                print()
                print(
                    "❌ NIST request failed."
                )

                print(
                    f"URL: {url}"
                )

                print(
                    f"Error: {error}"
                )

                return None

            print()
            print(
                f"⚠️ Request failed. "
                f"Retry {attempt}/{retries}..."
            )

            time.sleep(
                2 * attempt
            )

    return None



def find_nist_id(compound_name):

    """
    Search NIST and select the parent compound.

    IMPORTANT:

    This function does NOT simply select the first C-number.

    It attempts to avoid isotopologues such as:

        [2H6]ethanol
        Ethanol-d1
        Ethanol-d5
        C2H3D3O

    For ethanol, it should select:

        C64175
    """

    compound_name = str(
        compound_name
    ).strip()

    if not compound_name:

        print(
            "❌ Empty compound name."
        )

        return None

    print()
    print("-" * 75)

    print(
        f"Searching NIST for: {compound_name}"
    )

    encoded_name = quote(
        compound_name,
        safe=""
    )

    search_url = (
        f"{NIST_BASE_URL}/cgi/cbook.cgi"
        f"?Name={encoded_name}"
        f"&Units=SI"
    )

    response = get_page(
        search_url
    )

    if response is None:

        return None

    soup = BeautifulSoup(
        response.text,
        "html.parser"
    )

    requested_name = normalize_name(
        compound_name
    )

    candidates = []

    # ========================================================
    # PARSE NIST LINKS
    # ========================================================

    for link in soup.find_all(
        "a",
        href=True
    ):

        href = link.get(
            "href",
            ""
        )

        match = re.search(
            r"[?&]ID=(C\d+)",
            href,
            re.IGNORECASE
        )

        if not match:
            continue

        nist_id = match.group(
            1
        ).upper()

        link_text = link.get_text(
            " ",
            strip=True
        )

        if not link_text:

            link_text = nist_id

        candidates.append(
            {
                "id": nist_id,
                "name": link_text,
                "href": urljoin(
                    NIST_BASE_URL,
                    href
                )
            }
        )


    unique_candidates = []

    seen = set()

    for candidate in candidates:

        key = (
            candidate["id"],
            candidate["name"]
        )

        if key not in seen:

            seen.add(key)

            unique_candidates.append(
                candidate
            )

    candidates = unique_candidates

    # ========================================================
    # DISPLAY CANDIDATES
    # ========================================================

    if candidates:

        print()
        print(
            f"NIST returned "
            f"{len(candidates)} candidate compound link(s)."
        )

        for candidate in candidates:

            print(
                f"   {candidate['id']} "
                f"-> {candidate['name']}"
            )

    else:

        print()
        print(
            "⚠️ No compound links were found."
        )

        return None

    # ========================================================
    # SCORE CANDIDATES
    # ========================================================

    scored = []

    for candidate in candidates:

        candidate_name = normalize_name(
            candidate["name"]
        )

        score = 0

        # ----------------------------------------------------
        # Exact name
        # ----------------------------------------------------

        if candidate_name == requested_name:

            score += 1000

        # ----------------------------------------------------
        # Exact requested name appears as a complete phrase
        # ----------------------------------------------------

        if re.search(
            r"\b" + re.escape(requested_name) + r"\b",
            candidate_name,
            re.IGNORECASE
        ):

            score += 500

        # ----------------------------------------------------
        # Avoid isotopologues
        # ----------------------------------------------------

        if is_isotopologue(
            candidate["name"]
        ):

            score -= 1000

        # ----------------------------------------------------
        # Prefer actual compound-name matches
        # ----------------------------------------------------

        if requested_name in candidate_name:

            score += 100

        # ----------------------------------------------------
        # Avoid generic data pages
        # ----------------------------------------------------

        generic_terms = [

            "gas phase thermochemistry data",
            "condensed phase thermochemistry data",
            "phase change data",
            "reaction thermochemistry data",
            "henry's law data",
            "gas phase ion energetics data",
            "ion clustering data",
            "ir spectrum",
            "mass spectrum",
            "gas chromatography",
            "switch to calorie-based units",

        ]

        if candidate_name in generic_terms:

            score -= 500

        scored.append(
            (
                score,
                candidate
            )
        )

    # ========================================================
    # SORT
    # ========================================================

    scored.sort(
        key=lambda item: item[0],
        reverse=True
    )

    best_score, best = scored[0]

    # ========================================================
    # IMPORTANT ETHANOL SAFETY CHECK
    # ========================================================

    # NIST ethanol parent compound is C64175.
    # This prevents accidental selection of isotope pages.

    if requested_name in {
        "ethanol",
        "ethyl alcohol",
        "ethylalcohol",
        "ethyl alcohol anhydrous",
    }:

        parent_candidate = None

        for candidate in candidates:

            candidate_name = normalize_name(
                candidate["name"]
            )

            if (
                candidate["id"] == "C64175"
                or candidate_name == "ethanol"
            ):

                parent_candidate = candidate

                break

        if parent_candidate is not None:

            best = parent_candidate

            print()
            print(
                "✅ Parent ethanol compound selected."
            )

    # ========================================================
    # OUTPUT
    # ========================================================

    print()
    print(
        "✅ Matching NIST compound found."
    )

    print(
        f"   Compound: {best['name']}"
    )

    print(
        f"   NIST ID: {best['id']}"
    )

    if is_isotopologue(
        best["name"]
    ):

        print()
        print(
            "⚠️ WARNING: selected result appears "
            "to be an isotopologue."
        )

    return best["id"]


# ============================================================
# FIND NIST IR PAGE
# ============================================================

def find_ir_pages(nist_id):

    """
    Find all NIST IR-SPEC pages for a compound.

    Example:

        C64175

    may have several IR spectra.
    """

    base_url = (
        f"{NIST_BASE_URL}/cgi/cbook.cgi"
        f"?ID={nist_id}"
        f"&Mask=2283"
    )

    response = get_page(
        base_url
    )

    if response is None:

        return []

    soup = BeautifulSoup(
        response.text,
        "html.parser"
    )

    pages = []

    seen = set()

    # ========================================================
    # SEARCH LINKS
    # ========================================================

    for link in soup.find_all(
        "a",
        href=True
    ):

        href = link.get(
            "href",
            ""
        )

        if "Type=IR-SPEC" not in href.upper():

            continue

        absolute_url = urljoin(
            NIST_BASE_URL,
            href
        )

        if absolute_url in seen:

            continue

        seen.add(
            absolute_url
        )

        pages.append(
            absolute_url
        )

    # ========================================================
    # DIRECT FALLBACK
    # ========================================================

    if not pages:

        direct_url = (
            f"{NIST_BASE_URL}/cgi/cbook.cgi"
            f"?ID={nist_id}"
            f"&Index=0"
            f"&Type=IR-SPEC"
        )

        response = get_page(
            direct_url
        )

        if response is not None:

            if "Infrared Spectrum" in response.text:

                pages.append(
                    direct_url
                )

    return pages


# ============================================================
# FIND JCAMP DOWNLOAD LINK
# ============================================================

def find_jcamp_link(
    ir_page_url,
    html
):

    """
    Find the actual JCAMP-DX download link
    on the NIST IR spectrum page.

    We do NOT construct:

        ?JCAMP=...

    ourselves.

    Instead we use NIST's actual download link.
    """

    soup = BeautifulSoup(
        html,
        "html.parser"
    )

    # ========================================================
    # SEARCH ALL LINKS
    # ========================================================

    possible_links = []

    for link in soup.find_all(
        "a",
        href=True
    ):

        href = link.get(
            "href",
            ""
        )

        text = link.get_text(
            " ",
            strip=True
        )

        combined = (
            f"{text} {href}"
        ).lower()

        score = 0

        if "jcamp" in combined:

            score += 100

        if "download" in combined:

            score += 50

        if "spectrum" in combined:

            score += 20

        if (
            "cgi/cbook.cgi" in href.lower()
            and (
                "jcamp" in href.lower()
                or "type=ir" in href.lower()
            )
        ):

            score += 20

        if score > 0:

            absolute_url = urljoin(
                ir_page_url,
                href
            )

            possible_links.append(
                (
                    score,
                    absolute_url,
                    text
                )
            )

    # ========================================================
    # SORT
    # ========================================================

    possible_links.sort(
        key=lambda item: item[0],
        reverse=True
    )

    # ========================================================
    # RETURN BEST
    # ========================================================

    if possible_links:

        best = possible_links[0]

        print()
        print(
            "✅ JCAMP-DX download link found."
        )

        print(
            f"   Link text: {best[2]}"
        )

        print(
            f"   URL: {best[1]}"
        )

        return best[1]

    # ========================================================
    # SEARCH RAW HTML
    # ========================================================

    jcamp_patterns = [

        r'href=["\']([^"\']*JCAMP[^"\']*)["\']',

        r'href=["\']([^"\']*jcamp[^"\']*)["\']',

    ]

    for pattern in jcamp_patterns:

        matches = re.findall(
            pattern,
            html,
            re.IGNORECASE
        )

        if matches:

            url = urljoin(
                ir_page_url,
                matches[0]
            )

            print()
            print(
                "✅ JCAMP-DX link found "
                "in page source."
            )

            print(
                url
            )

            return url

    return None


# ============================================================
# DOWNLOAD JCAMP
# ============================================================

def download_jcamp(
    compound_name,
    nist_id
):

    """
    Find the actual NIST IR page and download
    its JCAMP-DX spectrum.
    """

    print()
    print(
        "Downloading NIST IR spectrum..."
    )

    print(
        f"NIST ID: {nist_id}"
    )

    # ========================================================
    # FIND IR PAGES
    # ========================================================

    print()
    print(
        "Checking NIST IR spectrum page..."
    )

    ir_pages = find_ir_pages(
        nist_id
    )

    if not ir_pages:

        print()
        print(
            f"❌ No IR spectrum page found "
            f"for {nist_id}."
        )

        return None

    print()
    print(
        f"Found {len(ir_pages)} possible "
        f"IR spectrum page(s)."
    )

    for index, page in enumerate(
        ir_pages
    ):

        print(
            f"   IR page {index}: {page}"
        )

    # ========================================================
    # CHECK EACH IR PAGE
    # ========================================================

    for index, ir_page in enumerate(
        ir_pages
    ):

        print()
        print(
            f"Checking IR spectrum "
            f"{index + 1}/{len(ir_pages)}"
        )

        print(
            ir_page
        )

        response = get_page(
            ir_page
        )

        if response is None:

            continue

        # ====================================================
        # FIND JCAMP LINK
        # ====================================================

        jcamp_url = find_jcamp_link(
            ir_page,
            response.text
        )

        if jcamp_url is None:

            print()
            print(
                "⚠️ Could not find JCAMP link "
                "on this IR page."
            )

            continue

        # ====================================================
        # DOWNLOAD JCAMP
        # ====================================================

        print()
        print(
            "Downloading JCAMP-DX:"
        )

        print(
            jcamp_url
        )

        jcamp_response = get_page(
            jcamp_url
        )

        if jcamp_response is None:

            continue

        jcamp_data = (
            jcamp_response.text
        )

        # ====================================================
        # VERIFY JCAMP
        # ====================================================

        upper_data = (
            jcamp_data.upper()
        )

        if (
            "##JCAMP-DX=" not in upper_data
            and "##TITLE=" not in upper_data
            and "##XYDATA=" not in upper_data
        ):

            print()
            print(
                "⚠️ Response is not a valid JCAMP-DX file."
            )

            print()
            print(
                "Response preview:"
            )

            print(
                jcamp_data[:500]
            )

            continue

        # ====================================================
        # SAVE
        # ====================================================

        safe_name = make_safe_filename(
            compound_name
        )

        filepath = (
            SPECTRA_FOLDER
            / f"{safe_name}.jdx"
        )

        try:

            with open(
                filepath,
                "w",
                encoding="utf-8"
            ) as file:

                file.write(
                    jcamp_data
                )

        except OSError as error:

            print()
            print(
                "❌ Could not save JCAMP-DX."
            )

            print(
                f"Error: {error}"
            )

            return None

        print()
        print(
            "✅ JCAMP-DX download successful."
        )

        print(
            f"Saved to: {filepath}"
        )

        return filepath

    # ========================================================
    # FAILED
    # ========================================================

    print()
    print(
        "❌ NIST did not provide a usable "
        "JCAMP-DX file."
    )

    print(
        f"Compound: {compound_name}"
    )

    print(
        f"NIST ID: {nist_id}"
    )

    return None


# ============================================================
# JCAMP HEADER
# ============================================================

def get_jcamp_header(data):

    header = {}

    for raw_line in data.splitlines():

        line = raw_line.strip()

        if not line.startswith("##"):

            continue

        if "=" not in line:

            continue

        key, value = line[
            2:
        ].split(
            "=",
            1
        )

        key = key.strip().upper()

        value = value.strip()

        header[key] = value

    return header


# ============================================================
# DECODE JCAMP XYDATA
# ============================================================

def decode_xydata(
    data,
    header
):

    """
    Decode:

        ##XYDATA=(X++(Y..Y))

    Supports normal NIST JCAMP files.
    """

    # ========================================================
    # HEADER VALUES
    # ========================================================

    first_x = float(
        header.get(
            "FIRSTX",
            0
        )
    )

    last_x = float(
        header.get(
            "LASTX",
            0
        )
    )

    delta_x = float(
        header.get(
            "DELTAX",
            0
        )
    )

    n_points = int(
        float(
            header.get(
                "NPOINTS",
                0
            )
        )
    )

    x_factor = float(
        header.get(
            "XFACTOR",
            1
        )
    )

    y_factor = float(
        header.get(
            "YFACTOR",
            1
        )
    )

    print()
    print("=" * 75)
    print("JCAMP INFORMATION")
    print("=" * 75)

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

    # ========================================================
    # FIND XYDATA
    # ========================================================

    lines = data.splitlines()

    xy_start = None

    for index, line in enumerate(
        lines
    ):

        if line.strip().upper().startswith(
            "##XYDATA="
        ):

            xy_start = index + 1

            break

    if xy_start is None:

        raise ValueError(
            "##XYDATA= section was not found."
        )

    # ========================================================
    # DECODE
    # ========================================================

    x_values = []
    y_values = []

    current_x = None

    for line in lines[xy_start:]:

        line = line.strip()

        if not line:

            continue

        if line.upper().startswith(
            "##END="
        ):

            break

        if line.startswith("##"):

            break

        parts = line.split()

        if len(parts) < 2:

            continue

        try:

            starting_x = float(
                parts[0]
            )

        except ValueError:

            continue

        current_x = starting_x

        for position, value in enumerate(
            parts[1:]
        ):

            try:

                y_value = float(
                    value
                )

            except ValueError:

                continue

            x_value = (
                current_x
                + position * delta_x
            )

            x_values.append(
                x_value
            )

            y_values.append(
                y_value
            )

    # ========================================================
    # NUMPY
    # ========================================================

    x_values = np.asarray(
        x_values,
        dtype=float
    )

    y_values = np.asarray(
        y_values,
        dtype=float
    )

    # ========================================================
    # APPLY FACTORS
    # ========================================================

    x_values *= x_factor

    y_values *= y_factor

    # ========================================================
    # LIMIT NPOINTS
    # ========================================================

    if (
        n_points > 0
        and len(x_values) > n_points
    ):

        x_values = x_values[
            :n_points
        ]

        y_values = y_values[
            :n_points
        ]

    # ========================================================
    # VALIDATE
    # ========================================================

    print()
    print(
        f"Decoded points: {len(x_values)}"
    )

    if n_points > 0:

        if len(x_values) == n_points:

            print(
                "✅ NPOINTS verified."
            )

        else:

            print()
            print(
                "⚠️ NPOINTS mismatch."
            )

            print(
                f"NIST reports: {n_points}"
            )

            print(
                f"Decoded:       {len(x_values)}"
            )

    if len(x_values) < 10:

        raise ValueError(
            "Fewer than 10 spectral points decoded."
        )

    print()
    print(
        f"X range: "
        f"{x_values.min():.2f} - "
        f"{x_values.max():.2f} cm⁻¹"
    )

    print(
        f"Y range: "
        f"{y_values.min():.6f} - "
        f"{y_values.max():.6f}"
    )

    # ========================================================
    # LAST X CHECK
    # ========================================================

    if len(x_values) > 0:

        decoded_last = x_values[-1]

        if last_x != 0:

            difference = abs(
                decoded_last - last_x
            )

            if difference > max(
                abs(delta_x) * 2,
                1.0
            ):

                print()
                print(
                    "⚠️ LASTX differs from decoded "
                    "last X."
                )

    return (
        x_values,
        y_values
    )


# ============================================================
# DETERMINE SPECTRUM TYPE
# ============================================================

def determine_spectrum_type(header):

    yunits = header.get(
        "YUNITS",
        ""
    ).upper()

    print()
    print(
        f"YUNITS: "
        f"{yunits or 'Not specified'}"
    )

    if (
        "ABSORBANCE" in yunits
        or yunits.strip() in {
            "A",
            "AU"
        }
    ):

        return "ABSORBANCE"

    if (
        "TRANSMITTANCE" in yunits
        or "TRANSMISSION" in yunits
        or "%" in yunits
    ):

        return "TRANSMITTANCE"

    return "UNKNOWN"


# ============================================================
# STANDARDIZE TRANSMITTANCE
# ============================================================

def standardize_transmittance(
    y,
    spectrum_type
):

    y = np.asarray(
        y,
        dtype=float
    )

    if len(y) == 0:

        raise ValueError(
            "No Y values."
        )

    # ========================================================
    # ABSORBANCE
    # ========================================================

    if spectrum_type == "ABSORBANCE":

        print()
        print(
            "Converting Absorbance to Transmittance."
        )

        transmittance = (
            10.0 ** (-y)
        )

    # ========================================================
    # TRANSMITTANCE
    # ========================================================

    else:

        y_max = np.nanmax(
            y
        )

        y_min = np.nanmin(
            y
        )

        print()
        print(
            "Standardizing Transmittance."
        )

        print(
            f"Original Y range: "
            f"{y_min:.6f} - {y_max:.6f}"
        )

        if (
            y_max > 1.0
            and y_max <= 100.0
        ):

            print(
                "Interpreting Y as percent transmittance."
            )

            transmittance = (
                y / 100.0
            )

        else:

            print(
                "Interpreting Y as fractional transmittance."
            )

            transmittance = y.copy()

    # ========================================================
    # CLEAN
    # ========================================================

    transmittance = np.nan_to_num(
        transmittance,
        nan=0.0,
        posinf=1.0,
        neginf=0.0
    )

    transmittance = np.clip(
        transmittance,
        1e-6,
        1.0
    )

    return transmittance


# ============================================================
# CLEAN SPECTRUM
# ============================================================

def clean_spectrum(
    x,
    transmittance
):

    df = pd.DataFrame({

        "Wavenumber (cm⁻¹)": x,

        "Transmittance": transmittance

    })

    df = df.replace(
        [np.inf, -np.inf],
        np.nan
    )

    df = df.dropna()

    df = df.sort_values(
        "Wavenumber (cm⁻¹)"
    )

    df = df.drop_duplicates(
        subset=[
            "Wavenumber (cm⁻¹)"
        ],
        keep="first"
    )

    df = df.reset_index(
        drop=True
    )

    if len(df) < 10:

        raise ValueError(
            "Fewer than 10 valid points remain."
        )

    return df


# ============================================================
# SAVE IR CSV
# ============================================================

def save_ir_csv(
    compound_name,
    df
):

    safe_name = make_safe_filename(
        compound_name
    )

    csv_path = (
        CSV_FOLDER
        / f"{safe_name}_ir.csv"
    )

    try:

        df.to_csv(
            csv_path,
            index=False
        )

    except Exception as error:

        print()
        print(
            "❌ Failed to save IR CSV."
        )

        print(
            f"Error: {error}"
        )

        return None

    print()
    print(
        "✅ IR CSV saved:"
    )

    print(
        csv_path
    )

    return csv_path


# ============================================================
# SAVE DATABASE
# ============================================================

def save_to_database(
    compound_name,
    df
):

    try:

        save_spectra_to_db(
            DATABASE_PATH,
            compound_name,
            df,
            None
        )

        print()
        print(
            "✅ IR spectrum saved to database."
        )

        print(
            DATABASE_PATH
        )

        return True

    except Exception as error:

        print()
        print(
            "⚠️ Database save failed."
        )

        print(
            f"Error: {error}"
        )

        return False


# ============================================================
# PLOT SPECTRUM
# ============================================================

def plot_spectrum(
    compound_name,
    df
):

    safe_name = make_safe_filename(
        compound_name
    )

    plot_path = (
        PLOTS_FOLDER
        / f"{safe_name}.png"
    )

    x = df[
        "Wavenumber (cm⁻¹)"
    ].to_numpy()

    y = df[
        "Transmittance"
    ].to_numpy()

    try:

        plt.figure(
            figsize=(12, 6)
        )

        plt.plot(
            x,
            y,
            linewidth=1
        )

        plt.gca().invert_xaxis()

        plt.xlabel(
            "Wavenumber (cm⁻¹)"
        )

        plt.ylabel(
            "Transmittance"
        )

        plt.title(
            f"IR Spectrum - {compound_name}"
        )

        plt.grid(
            True,
            alpha=0.3
        )

        plt.tight_layout()

        plt.savefig(
            plot_path,
            dpi=300,
            bbox_inches="tight"
        )

        print()
        print(
            "✅ IR plot saved:"
        )

        print(
            plot_path
        )

        plt.show()

        plt.close()

        return plot_path

    except Exception as error:

        print()
        print(
            "⚠️ Could not create IR plot."
        )

        print(
            f"Error: {error}"
        )

        try:
            plt.close()
        except Exception:
            pass

        return None


# ============================================================
# LOAD JCAMP
# ============================================================

def load_jcamp(
    filepath
):

    try:

        with open(
            filepath,
            "r",
            encoding="utf-8",
            errors="replace"
        ) as file:

            data = file.read()

    except OSError as error:

        raise ValueError(
            f"Could not read JCAMP: {error}"
        )

    header = get_jcamp_header(
        data
    )

    title = header.get(
        "TITLE",
        ""
    )

    if title:

        print()
        print(
            f"JCAMP title: {title}"
        )

    spectrum_type = determine_spectrum_type(
        header
    )

    x, raw_y = decode_xydata(
        data,
        header
    )

    transmittance = standardize_transmittance(
        raw_y,
        spectrum_type
    )

    df = clean_spectrum(
        x,
        transmittance
    )

    return (
        df,
        header,
        spectrum_type
    )


# ============================================================
# PROCESS ONE COMPOUND
# ============================================================

def process_compound(
    compound_name
):

    print()
    print()
    print("#" * 75)

    print(
        f"PROCESSING: {compound_name}"
    )

    print(
        "#" * 75
    )

    # ========================================================
    # FIND NIST ID
    # ========================================================

    nist_id = find_nist_id(
        compound_name
    )

    if nist_id is None:

        return None

    # ========================================================
    # DOWNLOAD JCAMP
    # ========================================================

    jdx_path = download_jcamp(
        compound_name,
        nist_id
    )

    if jdx_path is None:

        return None

    # ========================================================
    # LOAD
    # ========================================================

    try:

        (
            df,
            header,
            spectrum_type
        ) = load_jcamp(
            jdx_path
        )

    except Exception as error:

        print()
        print(
            "❌ Failed to decode JCAMP."
        )

        print(
            f"Error: {error}"
        )

        return None

    # ========================================================
    # SUMMARY
    # ========================================================

    print()
    print("=" * 75)
    print("SPECTRUM SUMMARY")
    print("=" * 75)

    print(
        f"Compound: {compound_name}"
    )

    print(
        f"NIST ID: {nist_id}"
    )

    print(
        f"Original spectrum type: {spectrum_type}"
    )

    print(
        "Standardized Y units: Transmittance"
    )

    print(
        f"Data points: {len(df)}"
    )

    print(
        f"Wavenumber range: "
        f"{df['Wavenumber (cm⁻¹)'].min():.1f} - "
        f"{df['Wavenumber (cm⁻¹)'].max():.1f} cm⁻¹"
    )

    print(
        f"Transmittance range: "
        f"{df['Transmittance'].min():.6f} - "
        f"{df['Transmittance'].max():.6f}"
    )

    # ========================================================
    # SAVE CSV
    # ========================================================

    csv_path = save_ir_csv(
        compound_name,
        df
    )

    # ========================================================
    # DATABASE
    # ========================================================

    database_saved = save_to_database(
        compound_name,
        df
    )

    # ========================================================
    # PLOT
    # ========================================================

    plot_path = plot_spectrum(
        compound_name,
        df
    )

    # ========================================================
    # RETURN
    # ========================================================

    return {

        "compound": compound_name,

        "nist_id": nist_id,

        "jcamp_file": str(
            jdx_path
        ),

        "csv_file": (
            str(csv_path)
            if csv_path
            else None
        ),

        "plot_file": (
            str(plot_path)
            if plot_path
            else None
        ),

        "database_file": str(
            DATABASE_PATH
        ),

        "database_saved": database_saved,

        "spectrum_type": spectrum_type,

        "x": df[
            "Wavenumber (cm⁻¹)"
        ].to_numpy(),

        "y": df[
            "Transmittance"
        ].to_numpy(),

        "dataframe": df,

        "header": header
    }


# ============================================================
# PROCESS MULTIPLE COMPOUNDS
# ============================================================

def process_multiple_compounds():

    print()
    print("=" * 75)
    print(
        "NIST IR SPECTRUM DOWNLOADER"
    )
    print("=" * 75)

    print()
    print(
        "Enter compound names separated by commas."
    )

    compound_input = input(
        "\nCompound names: "
    ).strip()

    if not compound_input:

        print()
        print(
            "❌ No compound names entered."
        )

        return {}

    compounds = []

    for compound in compound_input.split(","):

        compound = compound.strip()

        if (
            compound
            and compound not in compounds
        ):

            compounds.append(
                compound
            )

    if not compounds:

        print()
        print(
            "❌ No valid compound names."
        )

        return {}

    results = {}

    for compound in compounds:

        try:

            result = process_compound(
                compound
            )

            results[compound] = result

        except KeyboardInterrupt:

            print()
            print(
                "⚠️ Processing interrupted."
            )

            results[compound] = None

            break

        except Exception as error:

            print()
            print(
                f"❌ Unexpected error processing "
                f"{compound}:"
            )

            print(
                error
            )

            results[compound] = None

    # ========================================================
    # SUMMARY
    # ========================================================

    print()
    print()
    print("=" * 75)
    print(
        "DOWNLOAD SUMMARY"
    )
    print("=" * 75)

    successful = 0
    failed = 0

    for compound, result in results.items():

        if result is None:

            failed += 1

            print()
            print(
                f"❌ {compound} --> FAILED"
            )

            continue

        successful += 1

        print()
        print(
            f"✅ {compound} --> SUCCESS"
        )

        print(
            f"    NIST ID: "
            f"{result['nist_id']}"
        )

        print(
            f"    Spectrum type: "
            f"{result['spectrum_type']}"
        )

        print(
            f"    Points: "
            f"{len(result['x'])}"
        )

        print(
            f"    JCAMP: "
            f"{result['jcamp_file']}"
        )

        print(
            f"    IR CSV: "
            f"{result['csv_file']}"
        )

        print(
            f"    Plot: "
            f"{result['plot_file']}"
        )

        print(
            f"    Database saved: "
            f"{result['database_saved']}"
        )

    print()
    print("-" * 75)

    print(
        f"Successful: {successful}"
    )

    print(
        f"Failed:     {failed}"
    )

    print(
        f"Database:   {DATABASE_PATH}"
    )

    print("=" * 75)

    return results


# ============================================================
# PROGRAM START
# ============================================================

if __name__ == "__main__":

    process_multiple_compounds()