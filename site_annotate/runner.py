# runner.py
import pandas as pd
import logging
from typing import Iterable, Tuple
from concurrent.futures import ThreadPoolExecutor, as_completed

import pyfaidx
from tqdm import tqdm

from . import io
from . import modisite
from .log import get_logger
from .constants import VALID_MODI_COLS

logger = get_logger(__file__)


def process_frame(key_frame, fa, fa_psp_ref=None):
    # logger.debug(f"processing {key_frame}")
    key, frame = key_frame
    # if 'Cont' not in key:
    #     print()
    # subfa = fa[(fa["id"] == key) & (~fa["id"].str.startswith(DECOY_FLAG))]
    # import ipdb; ipdb.set_trace()
    try:
        subfa = fa[key]
    except KeyError:
        logger.info(f"skipping key {key}, not found in fasta")
        return
    seqinfo = io.extract_info_from_header(subfa.name)
    # seqinfo["sequence"] = subfa.seq
    seqinfo["sequence"] = str(subfa)

    if fa_psp_ref is not None and "uniprot_id" in frame.columns:
        try:
            subfa_psp = fa_psp_ref[
                frame.uniprot_id.iloc[0]
            ]  # shoudl check to ensure only 1
        except KeyError:
            # logger.info(f"skipping key {key}, not found in psp fasta")
            pass
        seqinfo["psp"] = dict(
            name=subfa_psp.name,
            # sequence=subfa_psp.seq,
            sequence=str(subfa_psp),
            name_raw=subfa_psp.raw.split("\n")[0],
        )

    if len(subfa) == 0:
        logger.info(f"skipping key {key}, db mismatch")
        return
    # seqinfo = subfa.iloc[0].to_dict()
    res = modisite.main(frame, seqinfo)
    return res


def run_pipeline(
    # g: Iterable[Tuple[str, pd.DataFrame]],
    df: pd.DataFrame,
    fa: pyfaidx.Fasta,
    fa_psp_ref: pyfaidx.Fasta = None,
    cores=1,
) -> list:

    # import ipdb; ipdb.set_trace()
    modis_for_processing = set(df.columns) & set(VALID_MODI_COLS)

    g = df.groupby("protein")
    fullres = list()

    if cores == 1:
        for item in tqdm(g):
            res = process_frame(item, fa, fa_psp_ref)
            fullres.append(res)

    if cores > 1:
        with ThreadPoolExecutor() as executor:
            futures = {
                executor.submit(process_frame, item, fa, fa_psp_ref): item for item in g
            }

            for future in tqdm(
                as_completed(futures),
                total=len(futures),
                mininterval=0.4,
                smoothing=0.1,
            ):
                try:
                    result = future.result()
                    fullres.append(result)
                except Exception as e:

                    # print(f"An error occurred: {e}")
                    pass

    fullres = list(filter(None, fullres))
    return fullres
