import os
from pathlib import Path
import re
import pandas as pd
from mygene import MyGeneInfo
import json
import sqlitedict


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


def add_uniprot(df):
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
