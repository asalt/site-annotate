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


def psp_assigner(x, psp_sequence):
    # this is a bit of a hack to ensure we don't match the first position of the protein
    # if the original "position_start" is a second/repeated peptide occurance in the protein
    start = int(max(x["protein_start"] - 4, 0))
    pos = psp_sequence[int(max(x["protein_start"] - 4, 0))].find(x["peptide"]) + start
    return pos


def main(df: pd.DataFrame, seqinfo: dict, isobaric=True):
    """
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
    # breakpoint()

    # VALID_MODI_COLS = [ # here as example, this is defined in constants.py
    #     "sty_79_9663",
    #     "k_42_0106",
    #     "k_43_0058",
    # ]
    for col in VALID_MODI_COLS:
        if col not in df.columns:
            continue
        if len(df[~df[col].isna()]) == 0:
            continue

        # at the moment, test_runner.test_flow tests the following routine
        _data = df[~df[col].isna()][col]
        res = _data.apply(extract_positions)  # res is a pd.Series
        res_series_of_dfs = res.apply(position_dict_to_df)
        res_df = reset_inner_index(res_series_of_dfs)
        #

        df_positions = pd.merge(res_df, df, left_on="original_index", right_index=True)
        df_positions["position_absolut"] = (
            df_positions["position_relative"] + df_positions["protein_start"] - 1
        )
        # here try phosphosite plus position map check
        if "psp" in seqinfo:  # this does
            psp_sequence = seqinfo["psp"]["sequence"]
            psp_position_start = df_positions.apply(
                lambda x: psp_assigner(x, psp_sequence),
                axis=1,
            )
            df_positions["psp_position_start"] = psp_position_start
            df_positions["position_absolut_psp"] = (
                df_positions["position_relative"]
                + df_positions["psp_position_start"]
                - 1
            )

        df_positions["fifteenmer"] = df_positions.apply(
            lambda x: create_15mer(sequence, x["position_absolut"]), axis=1
        )

        df_positions["protein_length"] = len(sequence)

        if isobaric:
            df_positions = quant_isobaric_site(df_positions)

        best_probability_col = col + "_best_localization"

        maxprob = df_positions.groupby(["spectrum", "peptide"])[
            best_probability_col
        ].max()
        maxprob.name = "highest_prob"
        maxprob = maxprob.reset_index()
        df_positions = df_positions.merge(maxprob, on=["spectrum", "peptide"])

        df_positions_filtered = df_positions[
            (df_positions.prob > 0.5)
            | (df_positions.prob >= df_positions["highest_prob"])
        ]

        RESULTS[col] = df_positions_filtered

    return RESULTS
