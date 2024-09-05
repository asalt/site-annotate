# reduce.py
# reduce to site level
import re
from tokenize import group
import pandas as pd
from collections import defaultdict

import numpy as np

from tqdm import tqdm

from .constants import VALID_MODI_COLS
from . import log
from . import io

logger = log.get_logger(__file__)


tmt_pat = re.compile(r"TMT_(\d+)_intensity")


def make_nr(df_reduced) -> pd.DataFrame:
    #  not useable yet
    logger.info("makenr")

    mapping = io.get_isoform_hierarchy()

    df = pd.merge(
        df_reduced,
        mapping[["primary_select", "secondary_select", "protein_id"]],
        left_on="ENSP",
        right_on="protein_id",
        how="left",
    )

    if df_reduced.shape[0] != df.shape[0]:
        logger.error("merge failed, created more rows than started with")

    id_cols = [
        "fifteenmer",
        "ENSG",
        "AA",
        "geneid",
    ]  #'label']
    id_cols = [x for x in id_cols if x in df.columns]

    def condense_group(group, idxs, df, id_cols=id_cols):
        # idx = [0, 1]
        sel = df.loc[idxs]
        sel["protein"] = sel["protein"].fillna("")
        proteins = str.join(", ", map(str, sel["protein"].tolist()))
        aapos_list = str.join(", ", map(str, set(sel["position_absolut"].tolist())))

        try:
            primarysel = (
                str.join(
                    ",",
                    (
                        sel.apply(
                            lambda x: (
                                x["protein"] if x["primary_select"] == "Yes" else ""
                            ),
                            axis=1,
                        )
                    ),
                ).strip(",")
                or None
            )
        except Exception as e:
            pass

        secondarysel = (
            str.join(
                ",",
                (
                    sel.apply(
                        lambda x: (
                            x["protein"] if x["secondary_select"] == "Yes" else ""
                        ),
                        axis=1,
                    )
                ),
            ).strip(",")
            or None
        )
        #
        uniquesel = None
        bestsel = None
        if (
            bool(primarysel) == False
            and bool(secondarysel) == False
            and proteins.count(",") == 0
        ):
            uniquesel = proteins
            bestsel = None
        elif (
            bool(primarysel) == False
            and bool(secondarysel) == False
            and proteins.count(",") > 0
        ):
            bestsel = sel.sort_values(by="protein_length", ascending=False).iloc[0][
                "Protein"
            ]
            uniquesel = None

        finalsel = None
        reason = None
        if bool(primarysel):
            finalsel = primarysel
            reason = "primary"
        elif not bool(primarysel) and bool(secondarysel):
            finalsel = secondarysel
            reason = "secondary"
        elif uniquesel is not None:
            finalsel = uniquesel
            reason = "not primary or secondary, only one choice"
        elif bestsel is not None:
            finalsel = bestsel
            reason = "not primary or secondary, multiple choices, returning longest"

        #

        # str.join(", ", map(str, set(sel["Primary_select"].tolist())))
        # make unique before making string
        d = dict()
        d["Protein_Accessions"] = proteins
        d["AApos_List"] = aapos_list
        d["Primary_select"] = primarysel or ""
        d["Secondary_select"] = secondarysel or ""
        d["Unique_select"] = uniquesel or ""
        d["Best_select"] = bestsel or ""
        d["Protein_Accession_Selected"] = finalsel
        d["Reason_Selected"] = reason
        for k, v in zip(id_cols, group):
            d[k] = v

        # add quant values
        _cols = set(df.columns) - set(d.keys())
        _cols = [
            x
            for x in _cols
            if x not in ("AApos", "Protein", "Protein_id", "protein_length")
        ]
        # avoid some columns from the original table we do not want to add
        # in fact the only columns we want to add are the ones that start with TMT_ (or something else if it isn't a TMT experiment)

        for col in _cols:
            d[col] = sel.iloc[0][col]

        return d

    g = df.groupby(id_cols)
    groups = g.groups
    res = list()
    for group, idxs in groups.items():
        d = condense_group(group, idxs, df, id_cols=id_cols)
        res.append(d)


def _reduce_sites(df):

    common_cols = [
        "protein",
        "uniprot_id",
        "ENSP",
        #"mapped_proteins",
    ]

    # Group by 'fifteenmer' and other common columns
    groupby_cols = [
        "fifteenmer",
        "sitename",
    ] + [x for x in common_cols if x in df.columns]

    agg_dict2 = {
        "hyperscore": "max",
        "rtscore": "max",
        "delta_mass": "min",
        "highest_prob": "max",
        "spectrum": lambda x: "|".join(x),
    }
    _missing = [x for x in agg_dict2.keys() if x not in df.columns]
    for _m in _missing:
        logger.warning(f"column {_m} not found in df")
        if _m == "spectrum":
            df[_m] = ""
        else:
            df[_m] = np.nan

    if "intensity" not in df.columns:
        df["intensity"] = np.nan

    g = df.groupby(groupby_cols)
    # Summarize by calculating the sum of TMT intensities
    # agg1
    result1 = None
    if any([tmt_pat.match(col) for col in df.columns]):
        _agg_dict = {
            col: "sum" for col in df.columns if "TMT_" in col and "intensity" in col
        }
        # else:
        #     if 'intensity' not in df.columns:
        #         logger.error(f"no intensity or TMT intensity columns present")
        #         return
        #     _agg_dict = {col: "sum" for col in ('intensity',)}
        result1 = g.agg(_agg_dict)

    # agg 2

    result2 = g.agg(agg_dict2).rename(
        columns={
            k: f"{k}_best"
            for k in ["hyperscore", "rtscore", "delta_mass", "highest_prob"]
        }
    )
    result2 = result2.rename(columns={"spectrum": "spectra"})
    result2_1 = g.agg({"spectrum": lambda x: len(set(x))})
    result2_1 = result2_1.rename(columns={"spectrum": "nspectra"})
    result2 = result2.join(result2_1)

    # agg 3
    result3 = g.agg({"intensity": ["sum", "max"]})
    result3.columns = ["_".join(col).strip() for col in result3.columns.values]

    # result = pd.concat([result1, result2, result3], axis=1)
    result = pd.concat(
        [
            result2,
            result3,
            result1,
        ],
        axis=1,
    )

    # g2 = df.groupby(groupby_cols + ["spectrum"])
    # result2 = g2.agg({"intensity": "sum"})
    # result2 = result2.reset_index()

    result = result.reset_index()
    addl_cls = [
        "ENST",
        "ENSG",
        "geneid",
        "taxon",
        "symbol",
        "protein_start",
        "protein_end",
        "protein_start_psp",
        "position_relative",
        "position_absolut",
        "position_absolut_psp",
        # "sitename",
        "sitename2",
        "all_possible_positions",
        "all_possible_positions_psp",
        "AA",
        # "ENSP", # already added
    ]
    addl_cols = [x for x in addl_cls if x in df.columns]

    _meta = df.drop_duplicates(subset=groupby_cols)[groupby_cols + addl_cols]

    result_merged = _meta.merge(result, on=groupby_cols, how="outer")
    # result_merged2 = result_merged.merge(result2, on=groupby_cols, how='outer')

    return result_merged


def reduce_sites(data: dict, n_workers=None, **kwargs):

    newdata = dict()
    for key, df in tqdm(data.items()):
        logger.info(f"Reducing {key}")
        newdata[key] = _reduce_sites(df)

    return newdata
