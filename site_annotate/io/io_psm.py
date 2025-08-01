# io/io_psm.py

import pandas as pd
import janitor
import logging
from typing import Union, Sequence
import re

from pathlib import Path

from ..constants import RENAME, RENAME_SHORT
from ..constants import VALID_MODI_COLS, POSSIBLE_SITE_ID_COLS
from ..utils import data_generator

logger = logging.getLogger(__name__)


TMT_PATTERN = re.compile(r"(?:^|_)(1[2-3][0-9])([nc]?)\b", flags=re.IGNORECASE)


def convert_tmt_label(shorthand):
    """Convert shorthand TMT labels to standardized format."""
    # Remove 'TMT' or 'TMT_' prefix if present
    normalized_input = re.sub(r"^TMT[_]?", "", shorthand)

    # Explicit special case for 126
    if normalized_input == "126":
        return "TMT_126"

    # Handle N and C suffixes
    match = re.match(r"(\d+)[_]?([NC])$", normalized_input)
    if match:
        return f"TMT_{match.group(1)}_{match.group(2)}"

    # Default to "_N" if no suffix is present
    return f"TMT_{normalized_input}_N"


def extract_tmt_label(col: str) -> str | None:
    """
    Extract valid TMT label from column string, e.g., '127n' => 'TMT_127_N'
    Ignores other numbers.
    """
    match = TMT_PATTERN.search(col)
    if not match:
        return None
    num, suffix = match.groups()
    # label = f"{num}{suffix}"
    shorthand = f"{num}{suffix.upper()}" if suffix else num
    return convert_tmt_label(shorthand)


def validate_psm_file(df):
    exampledata = data_generator.generate_test_data(1)
    exampledata = janitor.clean_names(exampledata)
    # TODO check correctly
    setdiff = set(exampledata.columns) - set(df.columns)

    if "protein" in setdiff:
        raise ValueError(
            "`protein` not in psms file, this is used to assign isoform specific sites"
        )

    if len(set(df.columns) & set(VALID_MODI_COLS)) == 0:
        raise ValueError(
            f"""no modi column present
        should have a column of form {VALID_MODI_COLS[0]}
        """
        )

    if "protein_start" not in df.columns:
        raise ValueError(f"expected `protein_start` in input file")

    if "spectrum_file" not in df.columns:
        raise ValueError(f"expected `spectrum_file` in input file")

    return True
    # if 'peptide' not in


def get_rename_dict(sample_cols, protected=None):
    """
    Choose renaming strategy based on number of samples.
    """
    from ..constants import RENAME, RENAME_SHORT

    rename_source = RENAME_SHORT if len(sample_cols) < 10 else RENAME
    return update_rename(sample_cols, rename_source, protected)


def update_rename(cols, rename_mapping, protected=None):
    new_mapping = {}
    used_labels = set()

    for col in cols:
        if protected and col in protected:
            continue

        label = extract_tmt_label(col)
        if label and label not in used_labels:
            new_mapping[col] = label
            used_labels.add(label)
            continue  # skip fallback if confident TMT match

        # fallback: use provided mapping
        for key, val in rename_mapping.items():
            if key in col:
                renamed = f"{col}_{val}"
                if renamed not in new_mapping.values():
                    new_mapping[col] = renamed
                    break

    return new_mapping


# def best_match(key: str, candidates: list[str]) -> str | None:
#     """Try to find best matching candidate for key."""
#     norm_key = normalize_colname(key)
#     exact = [c for c in candidates if normalize_colname(c) == norm_key]
#     if exact:
#         return exact[0]
#     partial = [c for c in candidates if norm_key in normalize_colname(c)]
#     return partial[0] if partial else None


# def update_rename(cols: Sequence[str], rename_mapping: dict, protected=None) -> dict:
#     """
#     Build a renaming dict from existing column names and provided rename map.
#     Returns a new mapping dict: original_col -> renamed_col
#     """
#     result = {}
#     norm_cols = {normalize_colname(c): c for c in cols}
#     for key, new_label in rename_mapping.items():
#         if protected and key in protected:
#             continue

#         match = best_match(key, cols)
#         if match:
#             new_name = f"{match}_{new_label}"
#             result[match] = new_name
#         else:
#             logger.debug(f"No match for rename key '{key}' in provided columns.")

#     return result


def read_psm_file(psm_file: Union[str, Path]) -> pd.DataFrame:
    """Read and prepare a PSM file."""
    df = pd.read_csv(psm_file, sep="\t")
    df = janitor.clean_names(df)
    df = prepare_psm_file(df)
    validate_psm_file(df)
    return df


def prepare_psm_file(df: pd.DataFrame) -> pd.DataFrame:
    """
    Standardize and explode mapped_proteins column for downstream use.
    """
    # Fallback column assignments
    if "intensity" not in df.columns:
        for alt in ("area", "peakarea", "peak_area"):
            if alt in df.columns:
                df["intensity"] = df[alt]
                break

    if "mapped_proteins" not in df.columns and "alternative_proteins" in df.columns:
        df["mapped_proteins"] = df["alternative_proteins"]

    required_cols = ["peptide", "intensity", "protein", "mapped_proteins"]
    missing = [col for col in required_cols if col not in df.columns]
    if missing:
        raise ValueError(f"Invalid PSM file, missing columns: {missing}")

    # Rename sample columns using TMT mapping
    orig_cols = set(df.columns)
    rename_map = get_rename_dict(df.columns)

    rename_subset = {k: v for k, v in rename_map.items() if k in orig_cols}
    df = df.rename(columns=rename_subset)

    renamed = set(df.columns) - orig_cols
    if renamed:
        logger.info(f"Renamed columns: {rename_subset}")

    # Merge 'protein' and 'mapped_proteins', split and clean
    df["mapped_proteins"] = (
        df["protein"] + ", " + df["mapped_proteins"].fillna("").str.replace("@@", ", ")
    )
    df["mapped_proteins2"] = (
        df["mapped_proteins"]
        .apply(lambda x: x.split(", ") if isinstance(x, str) else [])
        .apply(lambda x: list(filter(None, set(x))))
    )

    df = df.explode("mapped_proteins2").reset_index(drop=True)
    df["protein"] = df["mapped_proteins2"]
    return df
