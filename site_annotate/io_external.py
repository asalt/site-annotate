import os
import re
from pathlib import Path
import requests
import pandas as pd

import pyfaidx
import janitor

from . import log
from .constants import PHOSPHOSITEPLUS_ANNOTATIONS

logger = log.get_logger(__file__)


def set_data_dir():
    # Path of the current script
    current_script_path = Path(__file__).resolve()
    # Path to the top-level data directory
    res = current_script_path.parent.parent / "data" / "phosphositeplus"
    logger.debug(f"data dir set to {res}")
    return res


data_dir = set_data_dir()


# print(data_dir)

FILE_CHOCIES = [
    "Phosphorylation_site_dataset",
    "Regulatory_sites",
]


def check_files_exist(file_names, directory=data_dir):
    """Check if all required files exist in the specified directory."""
    missing_files = []
    for file_name in file_names:
        if not os.path.exists(os.path.join(directory, file_name)):
            missing_files.append(file_name)
    if missing_files:
        print(f"Missing files: {missing_files}")
        return False
    return True


def load_data(file_name, directory="data"):
    """Load data from a specified file within the given directory."""
    file_path = os.path.join(directory, file_name)
    if os.path.exists(file_path):
        # Example of loading a CSV file, adjust the loading method as needed
        return pd.read_csv(file_path)
    else:
        print(f"File not found: {file_name}")
        return None


def get_psiteplus_file(file_or_abbv, **kwargs):

    skiprows = 3
    if "skiprows" in kwargs:
        skiprows = kwargs.pop("skiprows")

    if file_or_abbv in PHOSPHOSITEPLUS_ANNOTATIONS.keys():
        file = PHOSPHOSITEPLUS_ANNOTATIONS[file_or_abbv]
    elif file_or_abbv in PHOSPHOSITEPLUS_ANNOTATIONS.values():
        file = file_or_abbv
    else:
        logger.warning(f"file {file_or_abbv} not found in set constants, may fail")

    fullfile = data_dir / file

    df = pd.read_table(fullfile, skiprows=skiprows, **kwargs)
    df = janitor.clean_names(df)
    return df

    # import polars as pl

    # pl.Config.set_tbl_cols(10000)
    # psp = pl.read_csv(fullfile, separator="\t", skip_rows=3, infer_schema_length=10000)
    # colnames_series = pd.Series(index=psp.columns)
    # newnames = janitor.clean_names(colnames_series, axis="index")
    # mapping = {k: v for k, v in zip(colnames_series.index, newnames.index)}
    # psp = psp.rename(mapping)

    # return psp


def read_psite_fasta(file=None, skiprows=3) -> pyfaidx.Fasta:
    """
    custom func to handle phosphositeplus fasta files whereby the first few rows need to be skipped
    """

    if file is None:
        file = data_dir / "Phosphosite_seq.fasta"

    import tempfile

    falines = open(file, "rb").readlines()
    with tempfile.NamedTemporaryFile(mode="w") as tmpfile:
        tmpfile.writelines([x.decode() for x in falines[3:]])  # decode here not before
        fa = pyfastx.Fasta(tmpfile.name, key_func=lambda x: x.split("|")[-1])
    return fa


# # Example Usage
# required_files = ['phosphosite_data.csv', 'another_data_file.csv']
# if check_files_exist(required_files):
#     for file_name in required_files:
#         data = load_data(file_name)
#         # Process data as needed
# else:
#     print("Please ensure all required files are in the data directory.")


def maybe_add_uniprot(df):
    if "uniprot_id" in df.columns:
        return df
    if "acc_id" in df.columns:
        df = df.rename(columns={"acc_id": "uniprot_id"})
        return df

    proteins = df.protein.unique()

    gid_pattern = re.compile(r"(?<=geneid\|)(.*?)(?=\|)")
    ENSP_pattern = re.compile(r"(?<=ENSP\|)(.*?)(?=\|)")

    proteins_ENSP = list()
    for protein in proteins:
        searchres = ENSP_pattern.search(protein)
        if searchres is not None:
            res = searchres.group(0)

        proteins_ENSP.append(res)

    # proteins = map(lambda x: x.split("|")[1], proteins)
    from mygene import MyGeneInfo

    mg = MyGeneInfo()
