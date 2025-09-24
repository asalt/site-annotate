# reduce.py
# reduce to site level
import re
from tokenize import group
import pandas as pd
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed

import numpy as np

from tqdm import tqdm


from .constants import VALID_MODI_COLS
from . import log

logger = log.get_logger(__file__)


# tmt_pat = re.compile(r".*TMT_(\d+)_intensity")


# Ensure uniqueness by appending .1, .2, etc., to duplicates
def make_unique(series):
    seen = {}
    unique = []
    for value in series:
        if value not in seen:
            seen[value] = 0
            unique.append(value)
        else:
            seen[value] += 1
            unique.append(f"{value}.{seen[value]}")
    return unique


def make_site_id(df):
    if "sitename2" in df:
        site_id = (
            df[["sitename2", "ENSP", "fifteenmer"]].fillna("").agg("_".join, axis=1)
        )
    elif "sitename" in df:  #
        site_id = (
            df[["sitename", "ENSP", "fifteenmer"]].fillna("").agg("_".join, axis=1)
        )
    site_id = make_unique(site_id)
    df["site_id"] = site_id
    return df


def group_map(groups, df, **kwargs):
    res = list()
    for group, idxs in tqdm(groups.items()):
        d = condense_group(group, idxs, df)
        res.append(d)
    return res


def get_id_cols(df):
    id_cols = [
        "fifteenmer",
        "ENSG",
        "AA",
        "geneid",
    ]  #'label']
    id_cols = [x for x in id_cols if x in df.columns]
    return id_cols


# def condense_group(group, idxs, df, id_cols=id_cols):
def condense_group(group, idxs, df):
    # idx = [0, 1]

    id_cols = get_id_cols(df)

    sel = df.loc[idxs]
    sel["protein_id"] = sel["protein_id"].fillna("")
    proteins = str.join(", ", map(str, sel["protein_id"].tolist()))
    aapos_list = str.join(", ", map(str, set(sel["position_absolut"].tolist())))

    try:
        primarysel = (
            str.join(
                ",",
                (
                    sel.apply(
                        lambda x: (
                            x["protein_id"] if x["primary_select"] == "Yes" else ""
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
                        x["protein_id"] if x["secondary_select"] == "Yes" else ""
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
        bestsel = sel.sort_values(by="ccds_length", ascending=False).iloc[0][
            "protein_id"
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

    # add rest of values
    _cols = set(df.columns) - set(d.keys())
    _cols = [
        x for x in _cols if x not in ("AApos", "Protein", "Protein_id", "ccds_length")
    ]
    # avoid some columns from the original table we do not want to add
    # in fact the only columns we want to add are the ones that start with TMT_ (or something else if it isn't a TMT experiment)

    for col in _cols:
        d[col] = sel.iloc[0][col]

    return d


def make_nr(df_reduced, cores=1) -> pd.DataFrame:

    logger.info("makenr")

    for col in ("ENSG", "ENSP"):
        if col not in df_reduced.columns:
            logger.warning(f"{col} not present, cannot make nr")
            return df_reduced

    from .io import io

    mapping = io.get_isoform_hierarchy()

    df = pd.merge(
        df_reduced,
        mapping[
            [
                "primary_select",
                "secondary_select",
                "protein_id",
                "gene_id",
                "ccds_length",
            ]
        ],
        left_on=["ENSP", "ENSG"],
        right_on=["protein_id", "gene_id"],
        how="left",
    )

    if df_reduced.shape[0] != df.shape[0]:
        logger.error("merge failed, created more rows than started with")

    id_cols = get_id_cols(df)
    g = df.groupby(id_cols)
    groups = g.groups

    batches = [df]
    if cores > 1:
        # batch and send to jobs
        group_keys = list(groups)
        batch_size = len(group_keys) // cores
        idxs = [x for x in range(0, len(group_keys), batch_size)]
        idxs[-1] = idxs[-1] + (len(group_keys) - idxs[-1])
        group_batches = [group_keys[a:b] for a, b in zip(idxs[0:-1], idxs[1:])]
        # indices = [[*g.indices[grp] for grp in group_batch] for group_batch in group_batches]
        indices = [
            [i for grp in group_batch for i in g.indices[grp]]
            for group_batch in group_batches
        ]
        batches = [df.loc[ixs] for ixs in indices]

    res = list()

    # for group, idxs in groups.items():
    #     d = condense_group(group, idxs, df)
    #     res.append(d)

    # res_df = pd.DataFrame.from_dict(res)
    # res_df = make_site_id(res_df)
    # res_df.index = res_df['site_id']

    # return res_df

    with ProcessPoolExecutor(max_workers=cores) as executor:
        futures = {
            executor.submit(group_map, subdf.groupby(id_cols).groups, subdf)
            for subdf in batches
        }

        for future in tqdm(
            as_completed(futures),
            total=len(futures),
            mininterval=0.4,
            smoothing=0.1,
        ):
            result = future.result()
            res.extend(result)
            # try:
            #     result = future.result()
            #     res.extend(result)
            # except Exception as e:
            #     # print(f"An error occurred: {e}")
            #     pass

    res = list(filter(None, res))

    res_df = pd.DataFrame.from_dict(res)
    res_df = make_site_id(res_df)
    res_df.index = res_df["site_id"]

    return res_df


def _reduce_sites(df):

    common_cols = [
        "protein",
        "protein_id",
        "uniprot_id",
        "ENSP",
        # "mapped_proteins",
    ]

    # Group by 'fifteenmer' and other common columns
    groupby_cols = [
        "fifteenmer",
        "sitename",
    ] + [x for x in common_cols if x in df.columns and not df[x].isnull().all()]

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
    # if any([tmt_pat.match(col) for col in df.columns]):
    if any("TMT" in x for x in df.columns):
        _agg_dict = {
            col: "sum"
            for col in df.columns
            if "TMT_" in col and "intensity" in col and "ratio" not in col
        }
        _agg_dict.update(
            {col: "mean" for col in df.columns if "TMT_" in col and "ratio" in col}
        )

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
    import pdb

    # agg 3
    result3 = g.agg({"intensity": ["sum", "max"]})
    result3.columns = ["_".join(col).strip() for col in result3.columns.values]

    # result = pd.concat([result1, result2, result3], axis=1)
    result = pd.concat(
        filter(
            lambda x: x is not None,
            [
                result2,
                result3,
                result1,
            ],
        ),
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
