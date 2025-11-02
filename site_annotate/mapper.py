# mapper.py

import os
from collections import defaultdict
from pathlib import Path
import re
import pandas as pd
import json
import sqlitedict
from mygene import MyGeneInfo

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
    """Extract the token after an 'ENSP|' or 'EN|' label if present; otherwise
    try to find an Ensembl protein-like ID (ENS..P<digits>) anywhere in the string.

    Note: Validation of the extracted token is done by the caller.
    """
    s = str(protein)
    m = re.search(r"(?<=ENSP\|)([^|]+)", s)
    if m:
        return m.group(1)
    m = re.search(r"(?<=EN\|)([^|]+)", s)
    if m:
        return m.group(1)
    m = re.search(r"(ENS[A-Z]*P\d{3,})", s)
    if m:
        return m.group(1)
    return None


from functools import cache
from tqdm import tqdm


MEM = dict()
    

def fetch_uniprot_info_from_db(db, ENSP_values):
    results = {}
    # import ipdb; ipdb.set_trace()
    for ENSP in tqdm(ENSP_values):
        if ENSP in MEM:
            results[ENSP] = MEM[ENSP]
            continue
        if ENSP not in db:
            continue

        record = db[ENSP]
        if isinstance(record, dict) and ENSP in record and "uniprot" not in record:
            # legacy layout where entry is nested under the ENSP key
            record = record[ENSP]

        results[ENSP] = record
        MEM[ENSP] = record

    return results

mg = MyGeneInfo()
def fetch_uniprot_info_online(missing_ENSP):

    if not missing_ENSP:
        return []

    try:
        return mg.querymany(
            missing_ENSP,
            scopes="ensembl.protein",
            fields="uniprot.Swiss-Prot,uniprot.TrEMBL,name,symbol,other_names,entrezgene,taxid",
        )
    except Exception as exc:  # pragma: no cover - defensive network fallback
        logger.warning("Failed to fetch UniProt info online: %s", exc)
        return []


def fetch_uniprot_by_ensg(ensg_ids):
    if not ensg_ids:
        return []
    try:
        return mg.querymany(
            list(ensg_ids),
            scopes="ensembl.gene",
            fields="uniprot.Swiss-Prot,uniprot.TrEMBL,name,symbol,entrezgene,taxid",
        )
    except Exception as exc:
        logger.warning("Failed to fetch UniProt by ENSG: %s", exc)
        return []


def fetch_uniprot_by_enst(enst_ids):
    if not enst_ids:
        return []
    try:
        return mg.querymany(
            list(enst_ids),
            scopes="ensembl.transcript",
            fields="uniprot.Swiss-Prot,uniprot.TrEMBL,name,symbol,entrezgene,taxid",
        )
    except Exception as exc:
        logger.warning("Failed to fetch UniProt by ENST: %s", exc)
        return []


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
    proteins = set(df["protein"].dropna().tolist())
    # ENSPs = df["protein"].dropna().map(find_ENSP).tolist()

    # for protein, ENSP in df['protein'].dropna().map(find_ENSP).items():
    for protein in proteins: 
        uniprot_info = (final_items.get(protein, {}) or dict()).get("uniprot")
        if uniprot_info:
            swissprot = uniprot_info.get("Swiss-Prot")
            if isinstance(swissprot, list):  # pick longest sequence
                swissprot = resolve_multi_uniprot(swissprot)
            if not swissprot:
                trembl = uniprot_info.get("TrEMBL")
                if isinstance(trembl, list):
                    trembl = resolve_multi_uniprot(trembl)
                if trembl:
                    lookup_dict[protein] = trembl
                    continue
            if swissprot:
                lookup_dict[protein] = swissprot
            else:
                logger.debug(f"No Swiss-Prot/TrEMBL info for {protein}")
        else:
            logger.debug(f"No ENSP mapping/UniProt info for {protein}")

    final_items.update(lookup_dict)

    # s = set([zz.strip() for yy in  [yy.split(',') for yy in xx[:5]] for zz in yy])
    # if "mapped_proteins" in df:
    #     mapped_proteins = df.mapped_proteins.apply(
    #         lambda s: ', '.join(
    #             (ens if (ens := find_ENSP(tok.strip())) else tok.strip().split('|',1)[0])
    #             for tok in str(s).split(',')
    #             if tok.strip()
    #         )
    #   )
    #     df['mapped_proteins'] = mapped_proteins


        # this one is supposed to be faster, untested:

        # s = df['mapped_proteins'].fillna('').astype(str).str.split(',')
        # long = s.explode().str.strip()
        # long = long[long.ne('')]

        # ens = long.map(find_ENSP)
        # mapped = ens.fillna(long.str.split('|', n=1).str[0])

        # mapped_proteins = mapped.groupby(level=0).apply(', '.join)

    if "uniprot_id" in df.columns:
        df["uniprot_id"] = df["uniprot_id"].fillna(df["protein"].map(lookup_dict))
    else:
        df["uniprot_id"] = df["protein"].map(lookup_dict)
    return df


def _is_valid_uniprot_id(val: str) -> bool:
    if pd.isna(val):
        return False
    s = str(val).strip()
    if not s:
        return False
    if any(ch in s for ch in (" ", "|", ",", ";", ":")):
        return False
    import re as _re
    # UniProt accession patterns (primary and isoform agnostic). Approximate but strict.
    # 6-char: 1 letter, 1 digit, 3 alnum, 1 digit; 10-char: mix of letters/digits.
    if _re.fullmatch(r"[OPQ][0-9][A-Z0-9]{3}[0-9]", s):
        return True
    if _re.fullmatch(r"[A-NR-Z][0-9][A-Z][A-Z0-9][A-Z0-9][0-9]", s):
        return True
    if _re.fullmatch(r"[A-Z0-9]{10}", s):
        return True
    return False


def add_uniprot(df, keycol="protein"):
    """
    df has a column `protein` that is a fasta header string
    with embedded ENSP values.
    e.g.: gene|ENSP|ENSP00000440235|more
    This tries to map ENSP values to uniprot IDs
    in the future could add support for other identifiers
    """
    # If acc_id exists but uniprot_id does not, rename; otherwise keep existing column values
    if "acc_id" in df.columns and "uniprot_id" not in df.columns:
        df = df.rename(columns={"acc_id": "uniprot_id"})

    df = df.copy()

    if "protein_ids" in df.columns:
        pre_pi = df["uniprot_id"].notna().sum() if "uniprot_id" in df.columns else 0
        if "uniprot_id" not in df.columns:
            df["uniprot_id"] = pd.Series(pd.NA, index=df.index, dtype="object")

        def _extract_primary_uniprot(val):
            if pd.isna(val):
                return pd.NA
            def pick_accession(token: str) -> str | None:
                token = token.strip()
                if not token:
                    return None
                # Recognize UniProt forms only; ignore Ensembl-like tokens completely
                if '|' in token:
                    toks = [t for t in token.split('|') if t]
                    # Look for sp|ACC| or tr|ACC| patterns
                    for i in range(len(toks) - 1):
                        left, mid = toks[i], toks[i + 1]
                        if left.lower() in ("sp", "tr") and _is_valid_uniprot_id(mid):
                            return mid
                    # Fallback: any subtoken that is a valid accession
                    for t in toks:
                        if _is_valid_uniprot_id(t.strip()):
                            return t.strip()
                    return None
                if _is_valid_uniprot_id(token):
                    return token
                return None

            if isinstance(val, list):
                for item in val:
                    if not item:
                        continue
                    acc = pick_accession(str(item))
                    if acc:
                        return acc
                return pd.NA
            text = str(val)
            if not text:
                return pd.NA
            for part in text.split(';'):
                part = part.strip()
                if not part:
                    continue
                acc = pick_accession(part)
                if acc:
                    return acc
            return pd.NA

        protein_ids_uniprot = df["protein_ids"].apply(_extract_primary_uniprot)
        df["uniprot_id"] = df["uniprot_id"].fillna(protein_ids_uniprot)

        post_pi = df["uniprot_id"].notna().sum()
        assigned_from_protein_ids = post_pi - pre_pi
        if assigned_from_protein_ids > 0:
            logger.info("UniProt IDs: %d assigned from Protein.Ids", assigned_from_protein_ids)

        if df["uniprot_id"].notna().all():
            return df

    # Build a mapping from the original 'protein' header to an ENSP token
    def _looks_like_ensembl_protein_id(s: str) -> bool:
        return bool(isinstance(s, str) and re.fullmatch(r"ENS[A-Z]*P\d{3,}", s))

    if keycol == "protein":
        proteins_ENSP = {}
        invalid_count = 0
        for protein in df["protein"].dropna().astype(str).unique():
            ens = find_ENSP(protein)
            if ens and not _looks_like_ensembl_protein_id(ens):
                ens = None
                invalid_count += 1
            proteins_ENSP[protein] = ens
        if invalid_count:
            logger.info("Headers with invalid ENSP tokens skipped: %d", invalid_count)
    else:
        proteins_ENSP = {}
        if keycol in df.columns and "protein" in df.columns:
            for header, token in zip(
                df["protein"].astype(str), df[keycol].astype(str)
            ):
                if pd.isna(header) or pd.isna(token):
                    continue
                extracted = find_ENSP(token)
                if extracted is None or not _looks_like_ensembl_protein_id(extracted):
                    continue
                proteins_ENSP[header] = extracted

    ensps = [val for val in proteins_ENSP.values() if val is not None]
    final_items = {}

    # import ipdb; ipdb.set_trace()
    if ensps:
        with get_db() as db:
            final_items = fetch_uniprot_info_from_db(db, ensps)
            db_hits = set(final_items.keys())

            missing_ENSP = sorted(set(ensps) - set(final_items.keys()))
            if missing_ENSP:
                online_info = fetch_uniprot_info_online(missing_ENSP)
                if online_info:
                    update_db(db, online_info)
                    final_items.update(
                        {
                            item["query"]: item
                            for item in online_info
                            if isinstance(item, dict) and "query" in item
                        }
                    )
                    online_hits = {item.get("query") for item in online_info if isinstance(item, dict)}
                else:
                    online_hits = set()
            else:
                online_hits = set()

            if db_hits:
                logger.info("UniProt cache hits: %d ENSP tokens", len(db_hits))
            if online_hits:
                logger.info("UniProt online fetch: %d ENSP tokens", len(online_hits))

    # xs = [ x for x in final_items.values() if isinstance(x.get('uniprot', {}).get('Swiss-Prot'), list) ]

    # Map back to original protein headers for downstream lookup
    final_items_original_value = {
        header: final_items.get(ensp)
        for header, ensp in proteins_ENSP.items()
        if ensp is not None
    }



    pre_count = df["uniprot_id"].notna().sum() if "uniprot_id" in df.columns else 0
    df = map_proteins_to_uniprot(df, final_items_original_value)
    # Fallback: try gene-level mapping for any proteins still missing UniProt
    if "uniprot_id" not in df.columns:
        df["uniprot_id"] = pd.NA
    missing_mask = df["uniprot_id"].isna()
    if missing_mask.any():
        # Attempt by NCBI geneid (preferred) or gene symbol (with species)
        cols = [c for c in ("protein", "geneid", "symbol", "taxon") if c in df.columns]
        candidates = df.loc[missing_mask, cols].copy()
        # Coerce numeric geneid if needed
        geneids = set()
        if "geneid" in candidates.columns:
            for val in candidates["geneid"].dropna().astype(str):
                s = str(val).strip()
                if s.isdigit():
                    geneids.add(int(s))
        if geneids:
            gi_info = fetch_uniprot_by_geneids(geneids)
        else:
            gi_info = []

        symbols = set()
        if "symbol" in candidates.columns:
            symbols = set(
                str(s).strip() for s in candidates["symbol"].dropna().astype(str)
                if str(s).strip() and str(s).strip().upper() != "NAN"
            )
        # Try symbols (human) regardless; MyGene dedupes
        sym_info = fetch_uniprot_by_symbols(symbols, species=9606) if symbols else []

        # Also try Ensembl gene/transcript fallbacks if present
        ensg_set = set()
        if "ENSG" in df.columns:
            ensg_set = set(
                str(x).strip() for x in df.loc[missing_mask, "ENSG"].dropna().astype(str)
                if str(x).strip()
            )
        enst_set = set()
        if "ENST" in df.columns:
            enst_set = set(
                str(x).strip() for x in df.loc[missing_mask, "ENST"].dropna().astype(str)
                if str(x).strip()
            )
        ensg_info = fetch_uniprot_by_ensg(ensg_set) if ensg_set else []
        enst_info = fetch_uniprot_by_enst(enst_set) if enst_set else []

        def _best_uniprot(rec: dict) -> str | None:
            unip = (rec or {}).get("uniprot") or {}
            sp = unip.get("Swiss-Prot")
            if isinstance(sp, list):
                sp = resolve_multi_uniprot(sp)
            if sp:
                return sp
            tr = unip.get("TrEMBL")
            if isinstance(tr, list):
                tr = resolve_multi_uniprot(tr)
            return tr or None

        gid_map = {}
        for item in gi_info or []:
            if not isinstance(item, dict):
                continue
            gid = item.get("query")
            acc = _best_uniprot(item)
            if gid and acc:
                gid_map[str(gid)] = acc

        sym_map = {}
        for item in sym_info or []:
            if not isinstance(item, dict):
                continue
            sym = item.get("query")
            acc = _best_uniprot(item)
            if sym and acc:
                sym_map[str(sym)] = acc

        # Build secondary maps
        ensg_map = {}
        for item in ensg_info or []:
            if not isinstance(item, dict):
                continue
            q = item.get("query")
            acc = _best_uniprot(item)
            if q and acc:
                ensg_map[str(q)] = acc

        enst_map = {}
        for item in enst_info or []:
            if not isinstance(item, dict):
                continue
            q = item.get("query")
            acc = _best_uniprot(item)
            if q and acc:
                enst_map[str(q)] = acc

        # Ensure dtype supports string assignments cleanly
        try:
            df["uniprot_id"] = df["uniprot_id"].astype("object")
        except Exception:
            pass
        filled = 0
        for idx, row in candidates.iterrows():
            if pd.notna(df.at[idx, "uniprot_id"]):
                continue
            gid = row.get("geneid")
            sym = row.get("symbol")
            acc = None
            if pd.notna(gid) and str(gid).isdigit():
                acc = gid_map.get(str(int(gid)))
            if not acc and pd.notna(sym):
                acc = sym_map.get(str(sym))
            if not acc and "ENSG" in df.columns:
                ensg = row.get("ENSG")
                if pd.notna(ensg):
                    acc = ensg_map.get(str(ensg))
            if not acc and "ENST" in df.columns:
                enst = row.get("ENST")
                if pd.notna(enst):
                    acc = enst_map.get(str(enst))
            if acc and _is_valid_uniprot_id(acc):
                df.at[idx, "uniprot_id"] = acc
                filled += 1
        if filled:
            logger.info("UniProt fallback by gene/symbol: %d rows", filled)
    # Sanitize obviously invalid UniProt IDs (accidental FASTA headers, etc.)
    invalid = df["uniprot_id"].notna() & ~df["uniprot_id"].apply(_is_valid_uniprot_id)
    if invalid.any():
        n_invalid = int(invalid.sum())
        logger.warning("Discarding %d invalid UniProt accessions (suspect headers)", n_invalid)
        df.loc[invalid, "uniprot_id"] = pd.NA
    post_count = df["uniprot_id"].notna().sum()
    assigned_from_mapping = post_count - pre_count
    if assigned_from_mapping > 0:
        logger.info("UniProt IDs: %d assigned via ENSP→UniProt mapping", assigned_from_mapping)
    # df['protein']= df.protein.map(proteins_ENSP)
    # if "mapped_proteins2" in df:
    #     df['mapped_proteins2'] = df.mapped_proteins2.map(proteins_ENSP)

    return df


def add_annotations(data):
    outdata = {}
    for k, df in data.items():
        has_uniprot = ("uniprot_id" in df.columns) and df["uniprot_id"].notna().any()
        # data[k] = add_uniprot(df)
        psp_info = io_external.get_psiteplus_file(k)
        if psp_info is None:
            logger.warning(f"no info for {k}")
            continue

        df["_upper"] = df["fifteenmer"].str.upper()
        df["uniprot_id"] = df["uniprot_id"].fillna("")
        psp_info["_upper"] = psp_info["site_+_7_aa"].str.upper()
        psp_info["acc_id"] = psp_info["acc_id"].fillna("")
        psp_info["psp_position"] = psp_info["mod_rsd"].str.extract(r"(\d+)")
        psp_info["psp_position"] = psp_info["psp_position"].astype(int)

        if has_uniprot:
            merge_keys_left = ["uniprot_id", "_upper"]
            merge_keys_right = ["acc_id", "_upper"]
            if "position_absolut_psp" in df.columns:
                merge_keys_left.append("position_absolut_psp")
                merge_keys_right.append("psp_position")
            else:
                logger.warning("No position_absolut_psp column; falling back to UniProt + 15-mer match for PSP mapping")
            dfm = df.merge(
                psp_info,
                left_on=merge_keys_left,
                right_on=merge_keys_right,
                how="left",
            )
        else:
            logger.warning(f"{k} has no UniProt IDs; trying symbol + 15-mer PSP mapping")
            if "symbol" not in df.columns:
                outdata[k] = df
                continue
            merge_keys_left = ["symbol", "_upper"]
            merge_keys_right = ["gene", "_upper"] if "gene" in psp_info.columns else ["gene", "_upper"]
            if "position_absolut_psp" in df.columns:
                merge_keys_left.append("position_absolut_psp")
                merge_keys_right.append("psp_position")
            elif "position_absolut" in df.columns:
                merge_keys_left.append("position_absolut")
                merge_keys_right.append("psp_position")
        dfm = df.merge(
            psp_info,
            left_on=merge_keys_left,
            right_on=merge_keys_right,
            how="left",
        )
        # If PSP supplied an accession and our uniprot_id is empty, backfill from PSP acc_id
        if "acc_id" in dfm.columns:
            if "uniprot_id" not in dfm.columns:
                dfm["uniprot_id"] = dfm["acc_id"]
            else:
                dfm["uniprot_id"] = dfm["uniprot_id"].fillna(dfm["acc_id"])
        dfm = dfm[[x for x in [*psp_info.columns, *df.columns] if x in dfm.columns]]
        outdata[k] = dfm
        if not len(df) == len(dfm):
            print(f"length of df is ", len(df))
            print(f"length of dfm is ", len(dfm))
            # This does not crash if not true, uses the overlap
            # this is by design
            # TODO: consider if this is a good idea
            # import ipdb; ipdb.set_trace()

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
    from .io.io import extract_info_from_header

    if col not in df.columns:
        return
    res = df[col].apply(extract_info_from_header)
    res_df = pd.DataFrame.from_records(res)
    df = df.join(res_df, how="left")
    return df


def build_fasta_token_index(fasta_df: pd.DataFrame) -> dict:
    token_cols = [col for col in ("ENSP", "ENST", "ENSG", "geneid", "symbol") if col in fasta_df.columns]
    index = defaultdict(set)

    proteins = fasta_df["protein"].astype(str).tolist()

    for key, header in zip(proteins, proteins):
        for token in header.split("|"):
            token = token.strip()
            if token:
                index[token].add(key)

    for col in token_cols:
        values = fasta_df[col].tolist()
        for key, val in zip(proteins, values):
            if pd.isna(val):
                continue
            for token in str(val).replace("|", ";").split(";"):
                token = token.strip()
                if token:
                    index[token].add(key)

    return {k: v for k, v in index.items()}


def build_fasta_uniprot_index(fasta_df: pd.DataFrame) -> dict:
    if "uniprot_id" not in fasta_df.columns:
        return {}

    index = defaultdict(set)
    for key, uni in zip(fasta_df["protein"].astype(str), fasta_df["uniprot_id"]):
        if pd.isna(uni):
            continue
        index[str(uni)].add(key)

    return {k: v for k, v in index.items()}


def build_uniprot_to_fasta_map(df: pd.DataFrame) -> dict:
    if df is None or df.empty:
        return {}

    index = defaultdict(set)
    for uni, key in zip(df["uniprot_id"], df["__fasta_key"]):
        if pd.isna(uni) or pd.isna(key):
            continue
        index[str(uni)].add(str(key))
    return {k: v for k, v in index.items()}


def build_fasta_synonym_headers(fasta_df: pd.DataFrame) -> list[str]:
    token_cols = [col for col in ("protein", "ENSP", "ENST", "ENSG", "geneid", "symbol") if col in fasta_df.columns]
    uniprot_present = "uniprot_id" in fasta_df.columns

    headers = []
    for _, row in fasta_df.iterrows():
        tokens = []
        for col in token_cols:
            val = row.get(col)
            if pd.isna(val):
                continue
            for token in re.split(r"[;|]", str(val)):
                token = token.strip()
                if token:
                    tokens.append(token)
        if uniprot_present and pd.notna(row.get("uniprot_id")):
            tokens.append(str(row["uniprot_id"]))

        tokens = [str(row["protein"])] + tokens + [str(row["protein"])]
        dedup = []
        seen = set()
        for token in tokens:
            token = token.strip()
            if token and token not in seen:
                dedup.append(token)
                seen.add(token)
        headers.append("|".join(dedup))

    return headers
def fetch_uniprot_by_geneids(geneids, species=9606):
    if not geneids:
        return []
    try:
        return mg.querymany(
            list(geneids),
            scopes="entrezgene",
            fields="uniprot.Swiss-Prot,uniprot.TrEMBL,name,symbol,entrezgene,taxid",
            species=species,
        )
    except Exception as exc:
        logger.warning("Failed to fetch UniProt by geneid: %s", exc)
        return []


def fetch_uniprot_by_symbols(symbols, species=9606):
    if not symbols:
        return []
    try:
        return mg.querymany(
            list(symbols),
            scopes="symbol",
            fields="uniprot.Swiss-Prot,uniprot.TrEMBL,name,symbol,entrezgene,taxid",
            species=species,
        )
    except Exception as exc:
        logger.warning("Failed to fetch UniProt by symbol: %s", exc)
        return []
