# reduce.py
# reduce to site level
import re
from tokenize import group
import pandas as pd
from collections import defaultdict
from tqdm import tqdm

from .constants import VALID_MODI_COLS
from . import log

logger = log.get_logger(__file__)


def _reduce_sites(df):

    common_cols = [
        "protein",
        "uniprot_id",
        "ENSP",
    ]

    # Group by 'fifteenmer' and other common columns
    groupby_cols = ["fifteenmer"] + [x for x in common_cols if x in df.columns]

    g = df.groupby(groupby_cols)
    # Summarize by calculating the sum of TMT intensities
    result1 = g.agg(
        {col: "sum" for col in df.columns if "TMT_" in col and "intensity" in col}
    )
    result2 = g.agg(
        {
            "hyperscore": "max",
            "rtscore": "max",
            "delta_mass": "min",
            "highest_prob": "max",
        }
    ).rename(
        columns={
            k: f"{k}_best"
            for k in ["hyperscore", "rtscore", "delta_mass", "highest_prob"]
        }
    )

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
        "protein_start",
        "protein_end",
        "protein_start_psp",
        "position_relative",
        "position_absolut",
        "position_absolut_psp",
        "AA",
    ]
    addl_cols = [x for x in addl_cls if x in df.columns]

    _meta = df.drop_duplicates(subset=groupby_cols)[groupby_cols + addl_cols]

    result_merged = _meta.merge(result, on=groupby_cols, how="outer")
    # result_merged2 = result_merged.merge(result2, on=groupby_cols, how='outer')

    return result_merged


def reduce_sites(data: dict, n_workers=None):

    # if n_workers is None:
    #     n_workers = len(Client().ncores()) - 1  # Using Dask Client to manage workers
    newdata = dict()
    for key, df in tqdm(data.items()):
        logger.info(f"Reducing {key}")
        newdata[key] = _reduce_sites(df)

    # Group by 'fifteenmer' and other common columns

    return newdata
