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

    # these get calculated later
    # if "protein_start" not in df.columns: 
    #     raise ValueError(f"expected `protein_start` in input file")

    # if "spectrum_file" not in df.columns:
    #     raise ValueError(f"expected `spectrum_file` in input file")

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


ALT_PEPTIDE_COLS = [
    "sequence",
    "stripped_sequence",
    # do not use these, these will not match the fasta sequence:
    #"sequence_modi",
    #"modified_peptide",
    #"modified_sequence",
]

ALT_INTENSITY_COLS = [
    "intensity",
    "precursor_area",
    "precursor_quantity",
    "sequence_area",
    "reporter_intensity",
]

ALT_PROTEIN_COLS = [
    "protein",
    "ensps_all",
    "protein_ids",
    "protein_group",
]


ALT_SPECTRUM_COLS = [
    "spectrum",
    "precursor_lib_index",
    "precursor.id",
]

ALT_SCORE_COLS = [
    "hyperscore",
    "xscore",
    "mass_evidence",
]

ALT_DELTA_MASS_COLS = [
    "delta_mass",
    "delta",
    "delta_ppm",
    "ms1_apex_mz_delta",
]


UNIMOD_TO_MOD_COLUMN = {
    "UniMod:21": "sty_79_9663",
    "UniMod:35": "m_15_9949",
    "UniMod:1": "k_42_0106",
    "UniMod:121": "k_114_0429",
}


def _assign_first_available(df: pd.DataFrame, target: str, candidates: list[str]) -> None:
    if target in df.columns:
        return
    for candidate in candidates:
        if candidate in df.columns:
            df[target] = df[candidate]
            return


def _coerce_intensity(df: pd.DataFrame) -> None:
    if "intensity" in df.columns:
        df["intensity"] = pd.to_numeric(df["intensity"], errors="coerce")


def _normalise_protein_columns(df: pd.DataFrame) -> None:
    if "protein" not in df.columns and "ensps_all" in df.columns:
        df["protein"] = df["ensps_all"].str.split(";").str[0]

    if "protein" not in df.columns and "protein_ids" in df.columns:
        df["protein"] = df["protein_ids"].str.split(";").str[0]

    if "protein" in df.columns:
        df["protein"] = df["protein"].astype(str).str.split(";").str[0]

    if "mapped_proteins" not in df.columns:
        if "ensps_all" in df.columns:
            df["mapped_proteins"] = df["ensps_all"].str.replace(";", ", ", regex=False)
        elif "protein_ids" in df.columns:
            df["mapped_proteins"] = df["protein_ids"].str.replace(";", ", ", regex=False)
        elif "protein" in df.columns:
            df["mapped_proteins"] = df["protein"]


def _sanitize_site_probs(value):
    if pd.isna(value):
        return pd.NA
    s = str(value).strip()
    if not s or s.lower() == "nan":
        return pd.NA
    s = re.sub(r"\d+$", "", s)
    s = re.sub(r"\(UniMod:\d+\)", "", s)
    s = s.replace("{", "(").replace("}", ")")
    return s


def _populate_mod_columns(df: pd.DataFrame) -> None:
    if "site_occupancy_probabilities" not in df.columns:
        return

    sanitized = df["site_occupancy_probabilities"].apply(_sanitize_site_probs)
    if sanitized.isna().all():
        return

    modified_series = df.get("modified_sequence")
    if modified_series is None:
        modified_series = df.get("modified_peptide")
    if modified_series is None:
        return

    unimod_lists = (
        modified_series.fillna("")
        .astype(str)
        .str.findall(r"UniMod:\d+")
    )

    observed_codes = {
        code for codes in unimod_lists if codes for code in codes
    }
    unknown_codes = observed_codes - set(UNIMOD_TO_MOD_COLUMN.keys())
    if unknown_codes:
        logger.debug("Unmapped UniMod codes encountered: %s", ", ".join(sorted(unknown_codes)))

    for unimod, column in UNIMOD_TO_MOD_COLUMN.items():
        mask = unimod_lists.apply(lambda codes: unimod in codes if codes else False)
        if not mask.any():
            continue

        if column not in df.columns:
            df[column] = pd.NA

        assign_mask = mask & sanitized.notna()
        if assign_mask.any():
            df.loc[assign_mask & df[column].isna(), column] = sanitized[assign_mask]

        best_col = f"{column}_best_localization"
        if best_col not in df.columns:
            df[best_col] = pd.NA

        def _best_prob(val):
            if pd.isna(val):
                return pd.NA
            matches = re.findall(r"\(([0-9]*\.?[0-9]+)\)", str(val))
            if not matches:
                return pd.NA
            best = max(float(m) for m in matches)
            return round(best, 3)

        df.loc[assign_mask, best_col] = df.loc[assign_mask, best_col].where(
            df.loc[assign_mask, best_col].notna(),
            df.loc[assign_mask, column].apply(_best_prob),
        )


def prepare_psm_file(df: pd.DataFrame) -> pd.DataFrame:
    """
    Standardize and explode mapped_proteins column for downstream use.
    """
    # Fallback column assignments for heterogeneous PSM exporters
    _assign_first_available(df, "peptide", ALT_PEPTIDE_COLS)
    _assign_first_available(df, "intensity", ALT_INTENSITY_COLS)
    _assign_first_available(df, "protein", ALT_PROTEIN_COLS)
    _assign_first_available(df, "spectrum", ALT_SPECTRUM_COLS)
    _assign_first_available(df, "score_msms", ALT_SCORE_COLS)
    _assign_first_available(df, "delta_mass", ALT_DELTA_MASS_COLS)

    if "score_msms" in df.columns:
        df["score_msms"] = pd.to_numeric(df["score_msms"], errors="coerce")
    if "delta_mass" in df.columns:
        df["delta_mass"] = pd.to_numeric(df["delta_mass"], errors="coerce")

    if "spectrum" in df.columns:
        df["spectrum"] = df["spectrum"].astype(str)

    if "mapped_proteins" not in df.columns and "alternative_proteins" in df.columns:
        df["mapped_proteins"] = df["alternative_proteins"]

    _normalise_protein_columns(df)
    _coerce_intensity(df)
    _populate_mod_columns(df)

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
