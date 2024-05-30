# modisite.py
import re
import pandas as pd
from collections import defaultdict

from . import log

logger = log.get_logger(__file__)


def quant_protein(**kwargs):
    pass


def modisite_quant():
    pass


def extract_positions(sequence):
    # Remove parentheses and numbers to calculate positions
    cleaned_sequence = re.sub(r"\(\d\.\d+\)", "", sequence)

    pattern = re.compile(r"(\w)\((\d\.\d+)\)")
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


def position_dict_to_df(position_dict):
    # Convert the dictionary to a DataFrame
    df = pd.DataFrame(position_dict).T.reset_index(names=["position"])
    return df


def reset_inner_index(res_df):
    df_list = []
    for idx, df in res_df.items():
        df_reset = df.reset_index(drop=True)
        df_reset["original_index"] = idx
        df_list.append(df_reset)
    concatenated_df = pd.concat(df_list, ignore_index=True)
    return concatenated_df


def create_15mer(sequence, position):
    sequence_padded = "_" * 7 + sequence + "_" * 7
    position_padded = position + 7
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


def main(df: pd.DataFrame, seqinfo: dict, isobaric=True):

    sequence = seqinfo["sequence"]
    VALID_COLS = ["sty_79_9663"]
    RESULTS = dict()
    # breakpoint()

    for col in VALID_COLS:
        if len(df[~df[col].isna()]) == 0:
            return

        res = df[~df[col].isna()][col].apply(extract_positions)
        res_df = res.apply(position_dict_to_df)
        res_df = reset_inner_index(res_df)

        df_positions = pd.merge(res_df, df, left_on="original_index", right_index=True)
        df_positions["position_absolut"] = (
            df_positions["position"] + df_positions["protein_start"] - 1
        )

        df_positions["fifteenmer"] = df_positions.apply(
            lambda x: create_15mer(sequence, x["position_absolut"]), axis=1
        )

        if isobaric:
            df_positions = quant_isobaric_site(df_positions)

        best_probability_col = col + "_best_localization"
        df_positions_filtered = df_positions[
            (df_positions.prob > 0.5)
            | (df_positions.prob == df_positions[best_probability_col])
        ]

        # psms_positions.groupby("15mer").apply(quant_protein)

        RESULTS[col] = df_positions_filtered

    return RESULTS
