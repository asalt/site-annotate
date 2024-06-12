# mapper.py

import os
from pathlib import Path
import re
import pandas as pd
from mygene import MyGeneInfo
import json
import sqlitedict

from site_annotate import io_external


from . import log

logger = log.get_logger(__file__)


def set_data_dir():
    # Path of the current script
    current_script_path = Path(__file__).resolve()
    # Path to the top-level data directory
    res = current_script_path.parent.parent / "data"
    logger.debug(f"data dir set to {res}")
    return res


data_dir = set_data_dir()
sqlitedict_filename = data_dir / "mappings.sqlite"


def get_db():
    return sqlitedict.SqliteDict(
        str(sqlitedict_filename),
        tablename="ensembl_uniprot_mapping",
        autocommit=False,
        encode=json.dumps,
        decode=json.loads,
    )


def update_db(db, items):
    for item in items:
        db[item["query"]] = item
    db.commit()


def find_ENSP(protein):
    ENSP_pattern = re.compile(r"(?<=ENSP\|)(.*?)(?=\|)")
    search_result = ENSP_pattern.search(protein)
    return search_result.group() if search_result else None


def fetch_uniprot_info_from_db(db, ENSP_values):
    return {ENSP: db[ENSP] for ENSP in ENSP_values if ENSP in db}


def fetch_uniprot_info_online(missing_ENSP):
    mg = MyGeneInfo()
    return mg.querymany(
        missing_ENSP,
        scopes="ensembl.protein",
        fields="uniprot.Swiss-Prot,uniprot.TrEMBL,name,symbol,other_names,entrezgene,taxid",
    )


def resolve_multi_uniprot(uniprot_list) -> str:
    from .io_external import read_psite_fasta

    fa = read_psite_fasta()  # this is cached/indexed so doesn't hurt to keep "reading"
    # could rename to "get_psite_fasta"

    _final_record = None
    _maxlen = 0
    for _id in uniprot_list:
        try:
            record = fa[_id]
        except KeyError:
            continue
        if len(record) > _maxlen:
            _final_record = _id
            _maxlen = len(record)
    return _final_record


def map_proteins_to_uniprot(df, final_items):
    lookup_dict = {}
    proteins = df["protein"].dropna().tolist()
    ENSPs = df["protein"].dropna().map(find_ENSP).tolist()

    # for protein, ENSP in df['protein'].dropna().map(find_ENSP).items():
    for protein, ENSP in zip(proteins, ENSPs):
        uniprot_info = final_items.get(ENSP, {}).get("uniprot")
        if uniprot_info:
            swissprot = uniprot_info.get("Swiss-Prot")
            if isinstance(swissprot, list):  # right now picks the longest
                picked = resolve_multi_uniprot(swissprot)
                swissprot = picked
            if swissprot:
                lookup_dict[protein] = swissprot
            else:
                logger.debug(f"No Swiss-Prot info for {protein}")
        else:
            logger.debug(f"No UniProt info for {protein}")
    df["uniprot_id"] = df["protein"].map(lookup_dict)
    return df


def add_uniprot(df):
    """
    df has a column `protein` that is a fasta header string
    with embedded ENSP values.
    e.g.: gene|ENSP|ENSP00000440235|more
    This tries to map ENSP values to uniprot IDs
    in the future could add support for other identifiers
    """
    if "uniprot_id" in df.columns or "acc_id" in df.columns:
        return (
            df.rename(columns={"acc_id": "uniprot_id"})
            if "acc_id" in df.columns
            else df
        )

    db = get_db()
    proteins_ENSP = {
        protein: find_ENSP(protein) for protein in df["protein"].dropna().unique()
    }
    final_items = fetch_uniprot_info_from_db(db, proteins_ENSP.values())

    missing_ENSP = set(proteins_ENSP.values()) - set(final_items.keys())
    # missing_ENSP = list(filter(missing_ENSP, None))
    missing_ENSP = [x for x in missing_ENSP if x is not None]

    if missing_ENSP:
        # import ipdb; ipdb.set_trace()
        online_info = fetch_uniprot_info_online(missing_ENSP)
        update_db(db, online_info)
        final_items.update({item["query"]: item for item in online_info})

    # xs = [ x for x in final_items.values() if isinstance(x.get('uniprot', {}).get('Swiss-Prot'), list) ]

    return map_proteins_to_uniprot(df, final_items)


def add_annotations(data):
    outdata = {}
    for k, df in data.items():
        if "uniprot_id" not in df.columns:
            outdata[k] = df
            logger.warning(f"{k} no uniprot id present, skipping")
            continue
        # data[k] = add_uniprot(df)
        psp_info = io_external.get_psiteplus_file(k)
        if psp_info is None:
            logger.warning(f"no info for {k}")
            continue

        df["_upper"] = df["fifteenmer"].str.upper()
        psp_info["_upper"] = psp_info["site_+_7_aa"].str.upper()

        dfm = df.merge(
            psp_info,
            left_on=["uniprot_id", "_upper"],
            right_on=["acc_id", "_upper"],
            how="left",
        )
        dfm = dfm[[x for x in [*psp_info.columns, *df.columns] if x in dfm.columns]]
        outdata[k] = dfm

    return outdata


def _xxadd_uniprot(df):
    """
    tries to extract ENSP value from `protein` column and then query mygene.info for uniprot info
    """
    if "uniprot_id" in df.columns:
        return df
    if "acc_id" in df.columns:
        df = df.rename(columns={"acc_id": "uniprot_id"})
        return df

    db = sqlitedict.SqliteDict(
        str(sqlitedict_filename),
        tablename="ensembl_uniprot_mapping",
        autocommit=False,  # too slow
        encode=json.dumps,
        decode=json.loads,
    )

    proteins = df.protein.unique()

    # gid_pattern = re.compile(r"(?<=geneid\|)(.*?)(?=\|)")
    ENSP_pattern = re.compile(r"(?<=ENSP\|)(.*?)(?=\|)")

    final_items = dict()  # maps ENSP to uniprot info
    proteins_ENSP = dict()  # maps long_protein string to ENSP value
    for ix, protein in enumerate(proteins):
        if pd.isna(protein):
            continue
        # if ix > 100:
        #     break
        searchres = ENSP_pattern.search(protein)
        if searchres is not None:
            ENSP_value = searchres.group()
        else:
            ENSP_value = None
        proteins_ENSP[protein] = ENSP_value
        # ====
        try:  # faster than checking if key is in db, maybe
            final_items[ENSP_value] = db[ENSP_value]
        except KeyError:
            pass

    # proteins = map(lambda x: x.split("|")[0], proteins)

    to_query = set(proteins_ENSP.values()) - set(final_items.keys())

    if len(to_query) != 0:
        mg = MyGeneInfo()

        mg_res = mg.querymany(
            proteins_ENSP,
            scopes="ensembl.protein",
            # fields="uniprot.Swiss-Prot,uniprot.TrEMBL,name,symbol,other_names,entrezgene,taxid,generif",
            fields="uniprot.Swiss-Prot,uniprot.TrEMBL,name,symbol,other_names,entrezgene,taxid",
        )

        for res in mg_res:
            # if "notfound" in res:
            #     continue
            db[res["query"]] = res
        db.commit()

        _ret_dict = {item["query"]: item for item in mg_res}
        final_items = final_items.update(_ret_dict)

    lookup_dict = dict()
    for key, ensp_value in proteins_ENSP.items():
        json_info = final_items.get(ensp_value)
        if json_info is None:
            continue
        uniprot_info = json_info.get("uniprot")
        if uniprot_info is None:
            logger.debug(f"no uniprot info for {key}")
            continue
        swissprot_value = uniprot_info.get("Swiss-Prot")
        if swissprot_value is None:
            logger.debug(f"no swissprot info for {key}")
            continue
        lookup_dict[key] = swissprot_value

    df["uniprot_id"] = df["protein"].map(lookup_dict)
    return df


def extract_keyvals_pipedsep(df, col="protein"):
    from .io import extract_info_from_header

    if col not in df.columns:
        return
    res = df[col].apply(extract_info_from_header)
    res_df = pd.DataFrame.from_records(res)
    df = df.join(res_df, how="left")
    return df
