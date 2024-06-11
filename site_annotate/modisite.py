# modisite.py
import re
import pandas as pd
from collections import defaultdict

from .constants import VALID_MODI_COLS
from . import log

logger = log.get_logger(__file__)


def extract_positions(sequence):
    # Remove parentheses and numbers to calculate positions
    cleaned_sequence = re.sub(r"\(\d\.\d+\)", "", sequence)

    pattern = re.compile(r"(\w)\((\d\.?\d*)\)")
    matches = pattern.finditer(sequence)

    result = {}
    position = 1
    for match in matches:
        # import ipdb; ipdb.set_trace()
        aa, prob = match.groups()
        # Find the correct position in the cleaned sequence
        position = cleaned_sequence.find(aa, position - 1) + 1
        result[position] = {"AA": aa, "prob": float(prob)}
        # Increment position to search for the next occurrence
        position += 1

    return result


def reset_inner_index(res_series_of_dfs: pd.Series) -> pd.DataFrame:
    """
    Concatenates a series or dictionary of DataFrames into a single DataFrame, while adding an 'original_index' column to track the original index or key of each DataFrame in the input series or dictionary.

    This function resets the index of each DataFrame, drops the existing index, and appends the original index or key as a new column called 'original_index'. If the input series or dictionary is empty, the function returns an empty DataFrame.

    Parameters:
    res_series_of_dfs (pd.Series or dict): A Series or dictionary where each element is a DataFrame. The keys or indices of this series or dictionary are used to populate the 'original_index' column in the resulting DataFrame.

    Returns:
    pd.DataFrame: A concatenated DataFrame with the original indices or keys stored in the 'original_index' column. If the input is empty, returns an empty DataFrame.

    Examples:
    >>> df1 = pd.DataFrame({'A': [1, 2]})
    >>> df2 = pd.DataFrame({'A': [3, 4]})
    >>> series_of_dfs = pd.Series({0: df1, 1: df2})
    >>> reset_inner_index(series_of_dfs)
       A  original_index
    0  1               0
    1  2               0
    2  3               1
    3  4               1

    # res_series_of_dfs can be a simple dict of dataframes as well



    """
    df_list = []
    for idx, df in res_series_of_dfs.items():
        df_reset = df.reset_index(drop=True)
        df_reset["original_index"] = idx
        df_list.append(df_reset)
    if len(df_list) == 0:
        return pd.DataFrame()
    concatenated_df = pd.concat(df_list, ignore_index=True)
    return concatenated_df


def position_dict_to_df(position_dict):
    """
    Converts a dictionary where the keys represent positions and the values represent attributes associated with those positions, into a pandas DataFrame. The keys become a column in the DataFrame named 'position_relative'. The values, which should themselves be dictionary-like, are expanded into additional columns.

    Parameters:
    position_dict (dict): A dictionary where the keys are positions (typically integers or strings that can be interpreted as such) and the values are dictionaries containing data associated with these positions.

    Returns:
    pd.DataFrame: A DataFrame where each row corresponds to one key-value pair from the input dictionary. The 'position_relative' column represents the keys from the original dictionary.

    Examples:
    >>> position_dict = {1: {'AA': 'S', 'prob': 0.95}, 2: {'AA': 'T', 'prob': 0.85}}
    >>> position_dict_to_df(position_dict)
       position_relative AA   prob
    0                  1  S  0.95
    1                  2  T  0.85
    """
    # Convert the dictionary to a DataFrame
    df = pd.DataFrame(position_dict).T.reset_index(names=["position_relative"])
    return df


def create_15mer(sequence, position):
    sequence_padded = "_" * 7 + sequence + "_" * 7
    position_padded = position + 7
    position_padded = int(position_padded)  # attempt to have this guaranteed beforehand
    res = sequence_padded[position_padded - 7 : position_padded + 7 + 1]
    reslist = list(res)
    reslist[6] = reslist[6].lower()
    return "".join(reslist)


def quant_isobaric_site(psms_positions):
    tmtsum = psms_positions.filter(like="TMT").sum(axis=1)
    psms_positions["tmt_sum"] = tmtsum

    ratios = psms_positions.filter(like="TMT").div(psms_positions.tmt_sum, axis=0)
    _rename = {col: f"{col}_ratio" for col in ratios.columns}
    ratios = ratios.rename(columns=_rename)
    #

    intensity_dstr = ratios.mul(psms_positions["intensity"], axis=0)
    _rename = {col: col.strip("ratio") + "intensity" for col in intensity_dstr.columns}
    intensity_dstr = intensity_dstr.rename(columns=_rename)

    fullres = pd.concat([psms_positions, ratios, intensity_dstr], axis=1)

    return fullres


def extract_and_transform_data(df, col):
    """
    Extract data and transform into DataFrame with original indices.
    Extract non-null data from the specified column, apply position extraction and transformation to a DataFrame format,
    and reset index while preserving the original index in a new column.

    Parameters:
    df (pd.DataFrame): The DataFrame from which to extract data.
    col (str): The column name in df from which to extract non-null values.

    Returns:
    pd.DataFrame: A DataFrame with the extracted and transformed data, including the original index.

    """
    _data = df[~df[col].isna()][col]
    res = _data.apply(extract_positions)  # res is a pd.Series
    res_series_of_dfs = res.apply(position_dict_to_df)
    res_df = reset_inner_index(res_series_of_dfs)
    return res_df


def enhance_dataframe(res_df, df, seqinfo) -> pd.DataFrame:
    """
    Merge and enhance DataFrame with absolute positions and PSP data if available.

    Parameters:
    res_df (pd.DataFrame): The DataFrame to be enhanced.
    df (pd.DataFrame): The original DataFrame to merge with res_df.
    seqinfo (dict): A dictionary containing sequence information, potentially including PhosphoSitePlus data.

    Returns:
    pd.DataFrame: The enhanced DataFrame with additional positional data.
    """

    df_positions = pd.merge(res_df, df, left_on="original_index", right_index=True)
    # df_positions["position_absolut"] = (
    # df_positions["position_relative"] + df_positions["protein_start"]
    # )
    orig_sequence = seqinfo["sequence"]
    # import ipdb; ipdb.set_trace()
    df_positions["position_start"] = df_positions.apply(
        lambda x: orig_sequence.find(x["peptide"]) + 1, axis=1
    )
    df_positions["protein_start"] = df_positions.apply(
        lambda x: orig_sequence.find(x["peptide"]) + 1, axis=1
    )
    df_positions["position_absolut"] = (
        df_positions["position_relative"] + df_positions["position_start"] - 1
    )

    # pos = psp_sequence.find(x["peptide"])

    if "psp" in seqinfo:
        orig_sequence = seqinfo["sequence"]
        psp_sequence = seqinfo["psp"]["sequence"]

        df_positions["position_start_psp"] = df_positions.apply(
            lambda x: psp_sequence.find(x["peptide"]) + 1, axis=1
        )
        df_positions["protein_start_psp"] = df_positions.apply(
            lambda x: psp_sequence.find(x["peptide"]) + 1, axis=1
        )
        df_positions["position_absolut_psp"] = (
            df_positions["position_start_psp"] + df_positions["position_relative"] - 1
        )
        # minus 1 because e.g.:
        # position_relative == 1
        # position_start = 500
        # sequence.find(peptide) == 499
        # but 1 index converts the result to 500
        # therefore we need to subtract 1 after adding the (1 indexed) relative position
        # add two 1 indexed sets together requires subtraction of 1 to remain 1 indexed
    return df_positions


def compute_additional_attributes(df_positions, sequence) -> pd.DataFrame:
    """
    Compute and append additional attributes such as fifteenmer sequences based on absolute positions and protein length.

    Parameters:
    df_positions (pd.DataFrame): The DataFrame to compute additional attributes for.
    sequence (str): The sequence string used for attribute calculations.

    Returns:
    pd.DataFrame: The DataFrame with additional attributes.
    """

    df_positions["fifteenmer"] = df_positions.apply(
        lambda x: create_15mer(sequence, x["position_absolut"]), axis=1
    )
    df_positions["protein_length"] = len(sequence)
    return df_positions


def process_probability_and_filter(df_positions, col, cutoff=0.5, take_best=True):
    """Calculate max probability and filter data."""
    """
    Calculate the highest probability of localization per spectrum and peptide, and filter the positions
    based on these probabilities.

    Parameters:
    df_positions (pd.DataFrame): DataFrame containing probability data to process.
    col (str): The column name to use for calculating best localization probabilities.

    Returns:
    pd.DataFrame: The DataFrame with positions filtered based on calculated probabilities.
    """

    best_probability_col = col + "_best_localization"
    maxprob = df_positions.groupby(["spectrum", "peptide"])[best_probability_col].max()
    maxprob.name = "highest_prob"
    maxprob = maxprob.reset_index()
    df_positions = df_positions.merge(maxprob, on=["spectrum", "peptide"])
    if take_best:
        _res = df_positions[
            (df_positions.prob > cutoff)
            | (df_positions.prob >= df_positions["highest_prob"])
        ]
    else:
        _res = df_positions[(df_positions.prob > cutoff)]

    return _res


def reorder_nice(df):
    order = [
        "protein",
        "uniprot_id",
        "fifteenmer",
        "peptide",
        "modified_peptide",
        "protein_start",
        "protein_end",
        "protein_start_psp",
        "position_relative",
        "position_absolut",
        "position_absolut_psp",
        "AA",
        "prob",
        "original_index",
        "spectrum",
        "spectrum_file",
        "extended_peptide",
        "prev_aa",
        "next_aa",
        "peptide_length",
        "charge",
        "retention",
        "observed_mass",
        "calibrated_observed_mass",
        "observed_m_z",
        "calibrated_observed_m_z",
        "calculated_peptide_mass",
        "calculated_m_z",
        "delta_mass",
        "spectralsim",
        "rtscore",
        "expectation",
        "hyperscore",
        "nextscore",
        "peptideprophet_probability",
        "number_of_enzymatic_termini",
        "number_of_missed_cleavages",
        "intensity",
        "assigned_modifications",
        "observed_modifications",
        "m_15_9949",
        "m_15_9949_best_localization",
        "sty_79_9663",
        "sty_79_9663_best_localization",
        "purity",
        "protein_id",
        "entry_name",
        "gene",
        "protein_description",
        "mapped_genes",
        "mapped_proteins",
        "is_unique",
        "TMT_126",
        "TMT_127_N",
        "TMT_127_C",
        "TMT_128_N",
        "TMT_128_C",
        "TMT_129_N",
        "TMT_129_C",
        "TMT_130_N",
        "TMT_130_C",
        "TMT_131_N",
        "TMT_131_C",
        "TMT_132_N",
        "TMT_132_C",
        "TMT_133_N",
        "TMT_133_C",
        "TMT_134_N",
        "protein_length",
        "tmt_sum",
        "TMT_126_ratio",
        "TMT_127_N_ratio",
        "TMT_127_C_ratio",
        "TMT_128_N_ratio",
        "TMT_128_C_ratio",
        "TMT_129_N_ratio",
        "TMT_129_C_ratio",
        "TMT_130_N_ratio",
        "TMT_130_C_ratio",
        "TMT_131_N_ratio",
        "TMT_131_C_ratio",
        "TMT_132_N_ratio",
        "TMT_132_C_ratio",
        "TMT_133_N_ratio",
        "TMT_133_C_ratio",
        "TMT_134_N_ratio",
        "TMT_126_intensity",
        "TMT_127_N_intensity",
        "TMT_127_C_intensity",
        "TMT_128_N_intensity",
        "TMT_128_C_intensity",
        "TMT_129_N_intensity",
        "TMT_129_C_intensity",
        "TMT_130_N_intensity",
        "TMT_130_C_intensity",
        "TMT_131_N_intensity",
        "TMT_131_C_intensity",
        "TMT_132_N_intensity",
        "TMT_132_C_intensity",
        "TMT_133_N_intensity",
        "TMT_133_C_intensity",
        "TMT_134_N_intensity",
        "highest_prob",
    ]
    _order = [x for x in order if x in df.columns]
    _extra = set(df.columns) - set(_order)
    return df[_order]


def main(df: pd.DataFrame, seqinfo: dict, isobaric=True) -> dict:
    """
    df is a DataFrame with columns:
        - protein
        - protein_start
        - peptide
        - spectrum
        - intensity
        - columns of form "sty_79_9663", "k_42_0106", ... that are probabilities of modification at a given position
          (e.g. S(0.0416)ES(0.8863)AENHS(0.0405)Y(0.0316)AK )
        - (optional) TMT_126, TMT_127, ... columns that are intensities of isobaric tags

    seqinfo should/may have keys of:
        - id
        - description
        - sequence
        - geneid
        - taxon
        - symbol
        - ENSP (only protein specific designator supported now)
    """

    sequence = seqinfo["sequence"]
    RESULTS = dict()
    # import ipdb; ipdb.set_trace()

    for col in VALID_MODI_COLS:
        if col not in df.columns or len(df[~df[col].isna()]) == 0:
            continue

        # import ipdb; ipdb.set_trace()
        res_df = extract_and_transform_data(df, col)
        df_positions = enhance_dataframe(res_df, df, seqinfo)
        df_positions = compute_additional_attributes(df_positions, sequence)

        if isobaric:
            df_positions = quant_isobaric_site(df_positions)

        df_positions_filtered = process_probability_and_filter(df_positions, col)
        df_positions_filtered = reorder_nice(df_positions_filtered)
        RESULTS[col] = df_positions_filtered

    return RESULTS
