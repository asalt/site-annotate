# io.py
import os
import re
import glob
import pathlib
from functools import partial
from pathlib import Path
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed
import janitor
from Bio import SeqIO
from pyfaidx import Fasta

from ..io_external import read_psite_fasta
from ..utils import data_generator
from .. import mapper

from .. import log
from .. import reduce
from ..constants import VALID_MODI_COLS, POSSIBLE_SITE_ID_COLS, RENAME, RENAME_SHORT

logger = log.get_logger(__file__)


def set_data_dir():
    # Path of the current script
    current_script_path = Path(__file__).resolve()
    # Path to the top-level data directory
    res = current_script_path.parent.parent.parent / "data"
    logger.debug(f"data dir set to {res}")
    return res


data_dir = set_data_dir()


def conf_to_dataframe(conf_file, set_index=False, dtype: dict = None):
    # Initialize config parser
    import configparser

    config = configparser.ConfigParser()
    config.read(conf_file)

    # Prepare dictionary for storing data
    data = {}

    # Loop through each section to create rows
    for section in config.sections():
        section_data = config[section]
        row = {key: section_data[key] for key in section_data}
        row["name"] = section  # Include section name as 'name' column
        data[section] = row

    # Convert dictionary to DataFrame
    df = pd.DataFrame.from_dict(data, orient="index")
    if dtype:
        for col, typ in dtype.items():
            if col in df.columns:
                df[col] = df[col].astype(typ)

    # Optional: Set 'name' as index
    if set_index:
        df.set_index("name", inplace=True)

    return df


def _xxget_reader(file: str):
    if file.endswith(".tsv"):
        return pd.read_table
    elif file.endswith(".csv"):
        return pd.read_csv
    elif file.endswith(".xlsx"):
        return pd.read_excel
    elif file.endswith(".conf"):
        return conf_to_dataframe
    else:
        logger.error(f"do not know how to read file: {file}")
        return None


def get_reader(file: str, **kwargs):
    def wrapped_reader(reader):
        def cleaner(*args, **kwargs):
            df = reader(*args, **kwargs)
            if isinstance(df, pd.DataFrame):
                obj_cols = df.select_dtypes(include="object").columns.tolist()
                # this will break other things -
                # can add optional arguments to control name cleaning
                # df = df.clean_names(strip_underscores=True)  # Standardize column names
                # df = df.clean_names(
                #     strip_underscores=True, column_names=obj_cols, axis=None
                # )  # uses janitor cleaning names routine

                # df[col] = df[col].apply(lambda x: preprocess_string(x) if isinstance(x, str) else x)
            return df  # If it's not a DataFrame, return as is

        return cleaner

    if file.endswith(".tsv"):
        return wrapped_reader(pd.read_table)
    elif file.endswith(".csv"):
        return wrapped_reader(pd.read_csv)
    elif file.endswith(".xlsx"):
        return wrapped_reader(pd.read_excel)
    elif file.endswith(".conf"):
        return wrapped_reader(conf_to_dataframe)
    else:
        logger.error(f"do not know how to read file")
        raise ValueError("do not know how to read file")
        # return None


def convert_tmt_label(shorthand):
    """Convert shorthand TMT labels to standardized format."""
    # Remove 'TMT' or 'TMT_' prefix if present
    normalized_input = re.sub(r"^TMT[_]?", "", shorthand)

    # Explicit special case for 126
    if normalized_input == "126":
        return "TMT_126"

    # Handle N and C suffixes
    match = re.match(r"(\d+)[_]?([NC])$", normalized_input)
    if match:
        return f"TMT_{match.group(1)}_{match.group(2)}"

    # Default to "_N" if no suffix is present
    return f"TMT_{normalized_input}_N"


def find_expr_file(rec_run_search: str, data_dir):  # TODO fix this
    """
    data_dir should be absolute path by this point
    """
    search_pattern = os.path.join(data_dir, f"{rec_run_search}*reduced*mapped*tsv")
    results = glob.glob(search_pattern, recursive=True)

    if not results:  # try one more time
        logger.info(f"trying again with not reduced")
        search_pattern = os.path.join(data_dir, f"{rec_run_search}*tsv")
        results = glob.glob(search_pattern)

    if len(results) == 0:
        logger.warning(f"no files found for {rec_run_search} in {data_dir}, skipping")

    for result in results:
        logger.debug(f"found {result} for {rec_run_search}")
        # raise FileNotFoundError(f"no files found for {rec_run_search} in {data_dir}")
    # if len(results) > 1:
    #     #if any('mapped'  in x for x in results):
    #     #    results = [x for x in results if 'mapped' in x]
    #     results =
    #     logger.warning(
    #         f"Ambiguous, found multiple files for {rec_run_search}, {str.join(', ', results)}"
    #     )
    #     raise FileNotFoundError(
    #         f"Ambiguous, found multiple files for {rec_run_search}, {str.join(', ', results)}"
    #     )
    # if len(results) == 1:
    return results


def find_expr_files(rec_run_searches, data_dir):
    return {rrs: find_expr_file(rrs, data_dir) for rrs in rec_run_searches}


def validate_metadata(df: pd.DataFrame, default_runno=1, default_searchno=7):
    if "recno" not in df.columns:
        logger.error(f"recno not present in metadata file")
        raise ValueError(f"recno not present in metadata file")
    if "runno" not in df.columns:  # assume 1
        logger.info(f"runno not found, assuming {default_runno}")
        df["runno"] = default_runno
    if "searchno" not in df.columns:  # assume 1
        logger.info(f"searchno not found, assuming {default_searchno}")
        df["searchno"] = default_searchno

    for x in ("recno", "runno", "searchno"):
        df[x] = df[x].astype(str)

    df["rec_run_search"] = df.apply(
        lambda x: f"{x.recno}_{x.runno}_{x.searchno}", axis=1
    )

    if "label" in df.columns:
        renamer = get_rename_dict(
            df.label.unique(), #protected=df.label.unique().tolist()
        )

        if not renamer:
            raise ValueError("Error trying to match columns")


        df["label"] = df["label"].map(renamer)
        df["label"] = df["label"].apply(convert_tmt_label)
        # can add a check here to assert unique by label and rec_run_search

    if "name" not in df.columns:
        logger.info("name not found, using rec_run_search")
        df["name"] = df["rec_run_search"]
    df.index = df["name"]

    return df


_ISOFORM_DF = None


def get_isoform_hierarchy() -> pd.DataFrame:
    global _ISOFORM_DF
    if _ISOFORM_DF is not None:
        return _ISOFORM_DF
    target1 = data_dir / "GENCODE.M32.basic.CHR.protein.selection.mapping.txt"
    target2 = data_dir / "GENCODE.V42.basic.CHR.isoform.selection.mapping.txt"
    logger.info(f"Reading {target1}")
    df1 = pd.read_table(target1, sep="\t", low_memory=False)
    logger.info(f"Reading {target2}")
    df2 = pd.read_table(target2, sep="\t", low_memory=False)
    df = pd.concat([df1, df2])
    df = janitor.clean_names(df)
    _ISOFORM_DF = df
    return df


def get_rename_dict(sample_cols, protected=None):
    """
    chooses RENAME or RENAME_SHORT
    """
    if len(sample_cols) < 10:  # then tmt 6plex (or similar) no C isotopes
        return update_rename(sample_cols, RENAME_SHORT, protected)
    return update_rename(sample_cols, RENAME, protected)


def update_rename(cols, rename_mapping: dict = None, protected=None) -> dict:
    if rename_mapping is None:
        rename_mapping = RENAME
    new_vals = dict()
    for key in rename_mapping:
        if protected is not None and key in protected:
            continue

        matches = [x for x in cols if key in x]
        if len(matches) == 0:
            continue
        if len(matches) > 1:
            logger.warning(f"too many columns match a single key {matches} - {key}")
            continue
        matchval = matches[0]
        # new_vals[matchval] = matchval + "_" + rename_mapping[key]
        new_vals[matchval] = rename_mapping[key]
    rename_mapping.update(new_vals)
    return rename_mapping

    #


# def prepare_psm_file(df: pd.DataFrame) -> pd.DataFrame:
#     """Check if a DataFrame is a valid PSM file.
#     this makes use of the proteins and mapped_proteins columns in msfragger output
#     """
#     # Check if a DataFrame is a valid PSM file
#     required_cols = ["peptide", "intensity", "protein", "mapped_proteins"]
#
#     if "intensity" not in df.columns:
#         for _x in ("area", "peakarea", "peak_area"):
#             if _x in df.columns:
#                 df["intensity"] = df[_x]
#     if "mapped_proteins" not in df.columns:
#         if "alternative_proteins" in df.columns:
#             df["mapped_proteins"] = df["alternative_proteins"]
#
#     for required_col in required_cols:
#         if required_col not in df.columns:
#             raise ValueError(f"Invalid PSM file, missing {required_col} column")
#
#     # if not any(df.columns.str.startswith("TMT")):
#     #renamer = update_rename(df.columns, RENAME)
#     renamer = get_rename_dict(df.columns)
#     orig_cols = set(df.columns)
#     renamer_subset = {k:v for k,v in renamer.items() if k in orig_cols} # only reason to do this is for logging info
#
#     df = df.rename(columns=renamer_subset)
#     new_cols = set(df.columns) - orig_cols
#     if len(new_cols) > 0:
#         logger.info(f"renamed {renamer_subset}")
#
#     df["mapped_proteins"] = (
#         df["protein"] + ", " + df["mapped_proteins"].fillna("").str.replace("@@", ", ")
#     )  # mapped_proteins column does not contain the value in the protein column so we add it here
#     df["mapped_proteins2"] = df["mapped_proteins"].apply(
#         lambda x: x.split(", ") if isinstance(x, str) else []
#     )
#     df["mapped_proteins2"] = df["mapped_proteins2"].apply(
#         lambda x: list(filter(None, set(x)))
#     )
#
#     # # Step 2: Explode the 'mapped_proteins' column
#     # df_exploded = df[[*required_cols, "mapped_proteins2"]].explode("mapped_proteins2")
#     df_exploded = df.explode("mapped_proteins2")
#     df_exploded = df_exploded.reset_index(drop=True)
#     df_exploded["protein"] = df_exploded["mapped_proteins2"]
#
#     # return df
#     return df_exploded


# def read_psm_file(psm_file: str | pathlib.Path) -> pd.DataFrame:
#     """Read a PSM file into a DataFrame."""
#     # Read a PSM file into a DataFrame
#     df = pd.read_csv(psm_file, sep="\t")
#     df = janitor.clean_names(df)
#     df = prepare_psm_file(df)
#     validate_psm_file(df)
#     return df


def read_sequence_file_to_dataframe(file_path, file_format="fasta"):
    records = SeqIO.parse(file_path, file_format)
    data = []
    for record in records:
        data.append(
            {
                "id": record.id,
                "description": record.description,
                "sequence": str(record.seq),
            }
        )
    df = pd.DataFrame(data)
    return df


VALID_NAMES = ["ENSP", "ENST", "ENSG", "geneid", "taxon", "symbol"]


def extract_info_from_header(header: str):
    """
    extracts key-value pairs of identifier|value and returns as separate columns
    """
    # example = 'ENSP|ENSMUSP00000022222|ENST|ENSMUST00000022222|ENSG|ENSMUSG00000021709|geneid|59079|taxon|10090|symbol|Erbin|Erbin'
    pattern = r"(\w+)\|([^|]+)"
    # Find all matches
    matches = re.findall(pattern, header)

    # Convert matches to a dictionary
    result = {key: value for key, value in matches}
    filtered_result = {key: result[key] for key in result if key in VALID_NAMES}
    return filtered_result
    # return result


def read_fasta(file_path):
    df = read_sequence_file_to_dataframe(file_path, "fasta")
    extracted_info = df["id"].apply(extract_info_from_header)
    # Convert the Series of dictionaries to a DataFrame
    extracted_df = pd.DataFrame(extracted_info.tolist())
    # Concatenate the new DataFrame with the original one
    df = pd.concat([df, extracted_df], axis=1)

    return df


def validate_psm_file(df):
    exampledata = data_generator.generate_test_data(1)
    exampledata = janitor.clean_names(exampledata)
    # TODO check correctly
    setdiff = set(exampledata.columns) - set(df.columns)

    if "protein" in setdiff:
        raise ValueError(
            "`protein` not in psms file, this is used to assign isoform specific sites"
        )

    if len(set(df.columns) & set(VALID_MODI_COLS)) == 0:
        raise ValueError(
            f"""no modi column present
        should have a column of form {VALID_MODI_COLS[0]}
        """
        )

    if "protein_start" not in df.columns:
        raise ValueError(f"expected `protein_start` in input file")

    if "spectrum_file" not in df.columns:
        raise ValueError(f"expected `spectrum_file` in input file")

    return True
    # if 'peptide' not in


def validate_expr_files(rec_run_searches: dict, meta_df: pd.DataFrame):
    """
    type is Dict[str] -> list
    """
    results = dict()
    meta_df = meta_df.copy()
    for rrs, expr_files in rec_run_searches.items():
        rec, run, search = rrs.split("_")

        for expr_file in expr_files:
            _meta = meta_df[
                (
                    (meta_df.recno == rec)
                    & (meta_df.runno == run)
                    & (meta_df.searchno == search)
                )
            ]
            meta_df.loc[_meta.index, "expr_col"] = None
            meta_df.loc[_meta.index, "expr_file"] = expr_file
            df = get_reader(expr_file)(expr_file, nrows=5)

            # rename_dict = RENAME
            # this will fail if it is no longer called "sample"
            sample_columns = [x for x in df.columns if x.startswith("sample")]
            rename_dict = get_rename_dict(sample_columns)
            # if (
            #     sample_columns and len(sample_columns) < 10
            # ):  # then tmt 6plex (or similar) no C isotopes
            #     rename_dict = RENAME_SHORT
            # df = df.rename(columns=rename_dict)

            # if "intensity_sum" not in df.columns:
            #     raise ValueError(f"`intensity_sum` not found in {expr_file}")
            if "intensity_sum" not in df.columns:
                logger.warning(
                    f"`intensity_sum` not found in {expr_file}. Downstream analyis might fail"
                )
                # raise ValueError(f"`intensity_sum` not found in {expr_file}")

            if "label" in _meta.columns:
                for label in _meta["label"].tolist():
                    label_mapping = [
                        x for x in df.columns if label in x and "ratio" not in x
                    ]  #
                    # ratio is a col we arent interested in
                    if len(label_mapping) == 0:
                        logger.warning(f"could not find sample with label {label}")
                        continue
                    if len(label_mapping) > 1:
                        logger.warning(f"too many results for {label}")
                        continue
                    if len(label_mapping) == 1:
                        expr_col = label_mapping[0]

                    ix = meta_df[
                        (meta_df["rec_run_search"] == rrs) & (meta_df["label"] == label)
                    ].index
                    if (len(ix)) != 1:
                        raise ValueError()

                        1 + 1
                    ix = ix[0]
                    meta_df.loc[ix, "expr_col"] = expr_col
                    meta_df.loc[ix, "expr_file"] = expr_file
            else:
                expr_col = "intensity_sum"
                meta_df.loc[_meta.index, "expr_col"] = "intensity_sum"

            basename = os.path.split(str(expr_file))[-1]
            key = os.path.splitext(basename)[0]
            key = re.sub("site_.*", "site", key)
            key = key.lstrip(rrs)
            results[key] = meta_df.copy()  # copy again
    return results


def merge_metadata(
    metadata: pd.DataFrame,
    filepath="siteinfo_combined",
    taxon=None,
    make_nr=True,
    **kwargs,
):
    """
    merge and write gct
    """
    # Load metadata file

    # Find unique expression file paths
    unique_expr_files = metadata["expr_file"].dropna().unique()

    # Dictionary to hold processed DataFrames for each file
    processed_data = {}
    # combined_rdesc = None

    # Iterate over each unique expression file
    for expr_file in unique_expr_files:  # this is messy, begin to refactor
        # Load the expression data
        expr_data = pd.read_csv(
            expr_file, sep="\t"
        )  # Assuming TSV format for expr files which is true

        # this is tmt-integrator specific
        if (
            "Index" in expr_data.columns and "site_id" not in expr_data.columns
        ):  # Index is of form ProteinID_AApos from tmt-integrator
            expr_data["site_id"] = expr_data["Index"]
        if (
            "SequenceWindow" in expr_data.columns
            and "fifteenmer" not in expr_data.columns
        ):
            expr_data["fifteenmer"] = expr_data["SequenceWindow"]
        if "sitename" not in expr_data.columns:
            if "Gene" in expr_data.columns and "ProteinID" in expr_data.columns:
                expr_data["Gene"] = expr_data["Gene"].fillna("")
                expr_data["ProteinID"] = expr_data["ProteinID"].fillna("")
                expr_data["sitename"] = expr_data.apply(
                    lambda x: x["Gene"] + "_" + x["Index"].lstrip(x["ProteinID"] + "_"),
                    axis=1,
                )
            else:
                expr_data["sitename"] = expr_data["site_id"]

        #

        # default philosopher samplenames are of form experiment_sample_xx
        sample_columns = [
            x
            for x in expr_data.columns
            if (("sample_" in x or x.endswith("intensity")) and "ratio" not in x)
        ]
        # rename_dict = get_rename_dict(sample_columns, protected=metadata.expr_col.dropna().tolist())
        # expr_data = expr_data.rename(
        #    columns=rename_dict
        # )  # a universal renamer for various possible names

        if "site_id" not in expr_data.columns:
            expr_data = reduce.make_site_id(expr_data)

        expr_data = expr_data.set_index("site_id")

        # Filter metadata rows corresponding to the current expression file
        file_metadata = metadata[metadata["expr_file"] == expr_file]

        assert file_metadata["expr_col"].unique().shape[0] == len(file_metadata)
        mapper = file_metadata[["name", "expr_col"]].set_index("expr_col").to_dict()
        name_mapper = mapper["name"]
        expr_data = expr_data.rename(columns=name_mapper)

        for col in ["ms_lit", "lt_lit"]:
            if col in expr_data.columns:
                expr_data[[col]] = expr_data[[col]].fillna(0)
                expr_data[[col]] = expr_data[[col]].astype(int)

        # to_exclude = {
        #     *expr_data.columns.difference(name_mapper.values()),
        #     "spectra",
        #     "nspectra",
        #     *(col for col in expr_data.columns if "_upper" in col)
        # }

        # to_exclude = (set(expr_data.columns) & set(RENAME.values())) - set(
        #     name_mapper.keys()
        # )  # this is a set of columns that are in the expr_data but not in the metadata

        # to_exclude |= {"spectra", "nspectra"} # we also don't want this when merging multiple
        # to_exclude |= {x for x in expr_data.columns if '_upper' in x}

        # Store the processed DataFrame with unique renames in the dictionary
        # to_exclude = (
        #     set(expr_data.columns)
        #     | set(col for col in expr_data.columns if "_upper" in col)
        #     | {"spectra", "nspectra"}
        #     - set(name_mapper.values())
        #     - set(POSSIBLE_SITE_ID_COLS)
        # )
        # expr_data = expr_data.drop(columns=to_exclude)

        emat = expr_data[name_mapper.values()]
        # rdesc_cols = list(set(expr_data.columns) & set(POSSIBLE_SITE_ID_COLS))
        # rdesc_cols = list(set(expr_data.columns) - set(emat.columns))
        rdesc_cols = [x for x in POSSIBLE_SITE_ID_COLS if x in expr_data]
        rdesc = expr_data[rdesc_cols]
        # if combined_rdesc is None:
        #     combined_rdesc = rdesc
        # combined_rdesc = combined_rdesc.join(rdesc)

        processed_data[os.path.basename(expr_file)] = {"emat": emat, "rdesc": rdesc}
        # processed_data[os.path.basename(expr_file)] = expr_data
    # Combine all processed DataFrames (aligning by index) to create a single combined DataFrame

    combined_emat = pd.concat(
        [data["emat"] for data in processed_data.values()], axis=1
    )

    # Initialize an empty DataFrame
    combined_rdesc = None

    # Iterate through rdesc DataFrames and combine
    for rdesc in [data["rdesc"] for data in processed_data.values()]:
        if combined_rdesc is None:
            combined_rdesc = rdesc  # Start with the first DataFrame
        else:
            combined_rdesc = combined_rdesc.combine_first(
                rdesc
            )  # Merge, prioritizing non-NaN values

    emat_only = set(combined_emat.index) - set(combined_rdesc.index)
    rdesc_only = set(combined_rdesc.index) - set(combined_emat.index)
    if len(emat_only) > 0:
        raise ValueError()
    if len(rdesc_only) > 0:
        raise ValueError()

    combined_rdesc = combined_rdesc.loc[combined_emat.index]

    for col in ("taxon",):
        if col not in combined_rdesc:
            continue
        combined_rdesc[col] = combined_rdesc[col].astype(str)

    # combined_df = pd.concat(processed_data['emat'].values(), axis=1)
    # combined_rdesc = pd.concat(processed_data['rdesc'].values(), axis=1)
    # emat = combined_df[metadata["name"]]
    # with the correct index "site_id" set above this will work

    # get rid of this, resolve earlier
    # if "Index" in combined_df.columns:
    #     combined_df = combined_df.set_index("Index")
    # elif "index" in combined_df.columns:
    #     combined_df = combined_df.set_index("index")
    # elif "site_id" in combined_df.columns:
    #     combined_df = combined_df.set_index("site_id")
    # else:
    #     pass

    # combined_df.to_csv("combined_df.tsv", sep="\t", index=False)

    # rdesc = everything not in emat  but in orig mat

    # if "Gene" in combined_df.columns:  # from tmt-integrator output
    #     combined_df["Gene"] = combined_df["Gene"].fillna("")
    #     if "sitename" not in combined_df.columns:
    #         sitename = combined_df.apply(
    #             lambda x: x["Gene"] + "_" + str(x.name).split("_")[-1], axis=1
    #         )
    #         combined_df["sitename"] = sitename

    # # rdesc_cols = list(set(combined_df.columns) - set(emat.columns))
    # rdesc_cols = list(set(combined_df.columns) - set(emat.columns))
    # rdesc = None
    # if rdesc_cols:
    #     rdesc = combined_df[list(rdesc_cols)]

    cdesc = metadata.loc[combined_emat.columns]

    # subselect before writing gct
    if taxon:
        if "taxon" not in combined_rdesc:
            raise ValueError("No taxon info available")
        combined_rdesc = combined_rdesc[combined_rdesc.taxon == taxon]

    combined_emat = combined_emat.loc[combined_rdesc.index]

    if make_nr:
        subsel = reduce.make_nr(combined_rdesc, **kwargs)
        combined_rdesc_subsel = combined_rdesc.loc[subsel.index]
        combined_emat_subsel = combined_emat.loc[subsel.index]
        write_gct(
            combined_emat_subsel,
            cdesc=cdesc,
            rdesc=combined_rdesc,
            filename=str(filepath) + "_nr",
        )

        combined_df_nr = combined_rdesc_subsel.join(combined_emat_subsel)
        write_excel(
            combined_df_nr,
            filename=str(filepath) + "_nr",
            shape=combined_emat_subsel.shape,
            column_metadata=cdesc,
        )

    write_gct(combined_emat, cdesc=cdesc, rdesc=combined_rdesc, filename=filepath)

    combined_df = combined_rdesc.join(combined_emat)
    write_excel(
        combined_df,
        filename=filepath,
        shape=combined_emat.shape,
        column_metadata=cdesc,
    )

    return combined_df


def write_gct(emat, cdesc, rdesc=None, filename="siteinfo_combined"):
    """
    Write expression data to a .gct file format using separate dataframes for emat and cdesc.

    Parameters:
    - emat (pd.DataFrame): Expression matrix where rows are genes and columns are samples.
    - cdesc (pd.DataFrame): Column descriptions with the same columns as emat.
    - filename (str): Path to save the .gct file.
    - rdesc (pd.DataFrame, optional): Row descriptions (metadata) for each gene. Should align with rows in emat.
    """
    # Prepare metadata and dimension information
    num_rows, num_cols = emat.shape

    # Ensure cdesc has the correct columns and order
    if not all(cdesc.index == emat.columns):
        raise ValueError(
            "Columns of cdesc must match the columns of emat in the same order."
        )

    if rdesc is None:
        rdesc = pd.DataFrame(index=emat.index)
        rdesc["rdesc"] = rdesc.index

    # Prepare the header lines
    header_lines = [  # import ipdb; ipdb.set_trace()
        "#1.3",
        f"{num_rows}\t{num_cols}\t{rdesc.shape[1]}\t{cdesc.shape[1]}",
    ]

    # Open file and write header and column descriptions
    outname = f"{filename}_{num_rows}x{num_cols}.gct"
    print(f"Writing to {outname}")

    with open(outname, "w") as f:
        # Write header lines
        for line in header_lines:
            f.write(line + "\n")

        # Write column names
        col_header = ["id", *list(rdesc.columns), *list(emat.columns)]
        f.write("\t".join(col_header) + "\n")

        # Write column descriptions (cdesc)
        # import ipdb

        # ipdb.set_trace()
        for desc_col in cdesc.columns:
            desc_str = (
                f"{desc_col}\t"
                + "na\t" * len(rdesc.columns)
                + "\t".join(map(str, cdesc[desc_col]))
            )
            f.write(desc_str + "\n")

        # Write expression data with optional row descriptions
        for gene_name, row_data in emat.iterrows():
            rdesc_str = "\t".join(map(str, rdesc.loc[gene_name]))
            row_str = f"{gene_name}\t{rdesc_str}\t" + "\t".join(map(str, row_data))

            f.write(row_str + "\n")


def read_gct(filepath: str):
    """
    Read a GCT v1.3 file written by write_gct and return (emat, cdesc, rdesc).

    Returns
    - emat: DataFrame indexed by row id with sample columns
    - cdesc: DataFrame indexed by sample with column metadata columns
    - rdesc: DataFrame indexed by row id with row metadata columns
    """
    import math

    with open(filepath, "r") as f:
        header1 = f.readline().strip()
        if not header1.startswith("#1.3"):
            raise ValueError("Unsupported GCT version; expected #1.3")

        dims = f.readline().strip().split("\t")
        if len(dims) < 4:
            raise ValueError("Malformed GCT header dimensions line")
        n_rows, n_cols, n_rdesc, n_cdesc = map(int, dims[:4])

        col_header = f.readline().rstrip("\n").split("\t")
        if len(col_header) < (1 + n_rdesc + n_cols):
            raise ValueError("Malformed GCT column header")

        id_col = col_header[0]
        rdesc_cols = col_header[1 : 1 + n_rdesc]
        sample_cols = col_header[1 + n_rdesc : 1 + n_rdesc + n_cols]

        # Read cdesc rows
        cdesc_rows = []
        for _ in range(n_cdesc):
            line = f.readline()
            if not line:
                raise ValueError("Unexpected EOF while reading cdesc")
            parts = line.rstrip("\n").split("\t")
            if len(parts) < (1 + n_rdesc + n_cols):
                raise ValueError("Malformed cdesc row")
            desc_name = parts[0]
            values = parts[1 + n_rdesc : 1 + n_rdesc + n_cols]
            cdesc_rows.append((desc_name, values))

        if cdesc_rows:
            # Build DataFrame with index=sample, columns=cdesc fields
            cdesc_dict = {name: values for name, values in cdesc_rows}
            cdesc = pd.DataFrame(cdesc_dict, index=sample_cols)
        else:
            cdesc = pd.DataFrame(index=sample_cols)

        # Read expression + rdesc rows
        ids, rdesc_data, emat_data = [], [], []
        for _ in range(n_rows):
            line = f.readline()
            if not line:
                raise ValueError("Unexpected EOF while reading data rows")
            parts = line.rstrip("\n").split("\t")
            if len(parts) < (1 + n_rdesc + n_cols):
                raise ValueError("Malformed data row")
            rid = parts[0]
            ids.append(rid)
            rdesc_vals = parts[1 : 1 + n_rdesc]
            rdesc_data.append(rdesc_vals)
            expr_vals = parts[1 + n_rdesc : 1 + n_rdesc + n_cols]
            # coerce to float when possible
            def _to_num(x):
                try:
                    v = float(x)
                    if math.isfinite(v):
                        return v
                except Exception:
                    pass
                return pd.NA

            emat_data.append([_to_num(x) for x in expr_vals])

        emat = pd.DataFrame(emat_data, index=ids, columns=sample_cols)
        rdesc = pd.DataFrame(rdesc_data, index=ids, columns=rdesc_cols if rdesc_cols else [])

    return emat, cdesc, rdesc


def write_excel(
    df: pd.DataFrame,
    filename="siteinfo_combined",
    shape=None,
    sheet_name="data",
    column_metadata: pd.DataFrame | None = None,
):
    """Write the merged data to an Excel workbook, including column metadata."""
    filename = str(filename)
    if shape is None:
        shape = df.shape
    num_rows, num_cols = shape
    outname = f"{filename}_{num_rows}x{num_cols}.xlsx"

    excel_df = df.reset_index()
    if "index" in excel_df.columns:
        excel_df = excel_df.rename(columns={"index": "site_id"})

    sheets = [(sheet_name, excel_df, (1, 1))]

    if column_metadata is not None:
        column_metadata_sheet = column_metadata.copy()
        idx_name = column_metadata_sheet.index.name or "sample"
        column_metadata_sheet = column_metadata_sheet.reset_index()
        if "index" in column_metadata_sheet.columns:
            column_metadata_sheet = column_metadata_sheet.rename(
                columns={"index": idx_name}
            )
        sheets.append(("column_metadata", column_metadata_sheet, (1, 0)))

    try:
        import xlsxwriter  # noqa: F401

        engine = "xlsxwriter"
    except ImportError:
        engine = None

    if engine:
        with pd.ExcelWriter(outname, engine=engine) as writer:
            for sheet, sheet_df, freeze in sheets:
                sheet_df.to_excel(writer, index=False, sheet_name=sheet)
                worksheet = writer.sheets[sheet]
                if freeze is not None:
                    worksheet.freeze_panes(*freeze)
                worksheet.autofilter(0, 0, sheet_df.shape[0], sheet_df.shape[1] - 1)

                for col_idx, col_name in enumerate(sheet_df.columns):
                    series = sheet_df[col_name].fillna("")
                    lengths = [len(str(value)) for value in series]
                    max_len = max(lengths + [len(str(col_name))]) if lengths else len(
                        str(col_name)
                    )
                    max_len = min(max_len, 50)
                    worksheet.set_column(col_idx, col_idx, max_len + 2)
    else:
        logger.warning(
            "xlsxwriter not available; writing Excel output without custom formatting"
        )
        with pd.ExcelWriter(outname) as writer:
            for sheet, sheet_df, _ in sheets:
                sheet_df.to_excel(writer, index=False, sheet_name=sheet)


def load_and_validate_files(psm_path, fasta_path, uniprot_check):

    # try load phosphositeplus fasta
    # needs to be downloaded from https://www.phosphosite.org/staticDownloads manually (free non commercial)
    # fa_psp_ref = load_psite_fasta()

    # logger.info(f"Loading {psm_path}")
    # df = io.read_psm_file(psm_path)

    # fasta_data = Fasta(fasta_path)
    from . import io_psm

    with ThreadPoolExecutor(
        max_workers=3
    ) as executor:  # there's some significant postprocessing these funcs do that makes this worth it, I think
        # need to time test this with larger files
        # Submit the function to the executor
        psp_future = executor.submit(read_psite_fasta)
        df_future = executor.submit(io_psm.read_psm_file, psm_path)
        fasta_future = executor.submit(
            partial(lambda x: Fasta(x, key_function=lambda x: x.split(" ")[0])),
            fasta_path,
        )

        fasta_data = fasta_future.result()
        #

        fa_psp_ref = psp_future.result()

        logger.info(f"Loading {psm_path}")
        df = df_future.result()

        # logger.info(f"Loading {fasta_path}")

    df = mapper.extract_keyvals_pipedsep(df)
    df = mapper.add_uniprot(df)

    fasta_keys = [str(k) for k in fasta_data.keys()]
    fasta_df = pd.DataFrame({"protein": fasta_keys})
    fasta_df = mapper.extract_keyvals_pipedsep(fasta_df)
    # import ipdb; ipdb.set_trace()
    if "ENSP" in fasta_df:
        fasta_df = mapper.add_uniprot(fasta_df, keycol="ENSP")
    else:
        fasta_df = mapper.add_uniprot(fasta_df, keycol="protein")

    fasta_token_index = mapper.build_fasta_token_index(fasta_df)
    fasta_uniprot_index = mapper.build_fasta_uniprot_index(fasta_df)
    enhanced_headers = mapper.build_fasta_synonym_headers(fasta_df)
    header_lookup = dict(zip(enhanced_headers, fasta_keys))

    protein_id_vals = df.protein.astype(str).unique()
    protein_id_mappings, lowmatches = fast_token_match(protein_id_vals, enhanced_headers)

    df["__fasta_key"] = df.protein.map(protein_id_mappings).map(header_lookup)

    unmatched_mask = df["__fasta_key"].isna()
    if unmatched_mask.any():
        logger.info(
            "Attempting UniProt-based FASTA mapping for %d proteins",
            unmatched_mask.sum(),
        )
        snapshots = df.loc[~unmatched_mask & df["uniprot_id"].notna(), ["uniprot_id", "__fasta_key"]]
        matched_uniprot = mapper.build_uniprot_to_fasta_map(snapshots)

        rows_to_drop = []
        rows_to_add = []

        for idx in df.index[unmatched_mask]:
            row = df.loc[idx]
            candidates = set()

            for col in ("ENSP", "ENST", "ENSG", "geneid", "symbol"):
                if col in row and pd.notna(row[col]):
                    for token in str(row[col]).split(";"):
                        token = token.strip()
                        if token:
                            candidates.update(fasta_token_index.get(token, set()))

            uni = row.get("uniprot_id")
            if pd.notna(uni):
                candidates.update(matched_uniprot.get(str(uni), set()))
                candidates.update(fasta_uniprot_index.get(str(uni), set()))

            if len(candidates) == 1:
                df.at[idx, "__fasta_key"] = next(iter(candidates))
            elif len(candidates) > 1:
                rows_to_drop.append(idx)
                for key in sorted(candidates):
                    dup = row.copy()
                    dup["__fasta_key"] = key
                    rows_to_add.append(dup)

        if rows_to_drop:
            df = df.drop(index=rows_to_drop)
        if rows_to_add:
            df = pd.concat([df, pd.DataFrame(rows_to_add)], ignore_index=True)

    else:
        logger.info("Protein identifiers already match FASTA keys")

    still_unmatched = df["__fasta_key"].isna()
    if still_unmatched.any():
        logger.warning(
            "Dropping %d PSM rows that could not be aligned to provided FASTA",
            still_unmatched.sum(),
        )
        df = df.loc[~still_unmatched].reset_index(drop=True)

    df["protein"] = df.pop("__fasta_key")

    lookup_cols = [col for col in ("ENSP", "ENST", "ENSG", "geneid", "symbol", "uniprot_id") if col in fasta_df.columns]
    if lookup_cols:
        lookup = fasta_df.set_index("protein")
        pre_fastafill = df["uniprot_id"].notna().sum() if "uniprot_id" in df.columns else 0
        for col in lookup_cols:
            if col not in df.columns:
                df[col] = pd.NA
            else:
                df[col] = df[col].replace("", pd.NA)
            df[col] = df[col].fillna(df["protein"].map(lookup[col]))
        post_fastafill = df["uniprot_id"].notna().sum() if "uniprot_id" in df.columns else 0
        delta = post_fastafill - pre_fastafill
        if delta > 0:
            logger.info("UniProt IDs: %d filled from FASTA annotations", delta)
        # Sanity check: drop header-like or invalid UniProt IDs
        if "uniprot_id" in df.columns:
            from .. import mapper as _mapper
            invalid = df["uniprot_id"].notna() & ~df["uniprot_id"].apply(_mapper._is_valid_uniprot_id)
            if invalid.any():
                nbad = int(invalid.sum())
                logger.warning("Discarding %d invalid UniProt accessions from FASTA fill", nbad)
                df.loc[invalid, "uniprot_id"] = pd.NA

    if fa_psp_ref:
        logger.info("FASTA data loaded successfully.")
    else:
        logger.info("Failed to load FASTA data.")

    if uniprot_check:
        logger.info("adding uniprot info")
        df = mapper.add_uniprot(df)

    # df = df[ df.mapped_proteins.str.contains("ENSP00000359491") ]
    # "ENSP00000246785") ]
    return df, fasta_data, fa_psp_ref  #


def fast_token_match(
    protein_ids, fasta_headers, fullname_col="fullname", min_score=1, tokenizer=None
):
    """
    Universal fast token-based matching between protein_ids and FASTA headers using a sparse matrix.
    Returns a DataFrame: protein_id, best_fasta_fullname, overlap_score.
    - protein_ids: list or array of protein strings (e.g. 'ENSMUSP00000133117|Fga')
    - fasta_mapper: DataFrame with at least a 'fullname' column (the full FASTA header)
    - min_score: minimum required overlap to consider a match (default 1)
    """
    import numpy as np
    from scipy.sparse import csr_matrix

    def _tokenizer(x):
        return set(x.replace("-", "|").replace("_", "|").lower().split("|"))

    if tokenizer is None:
        tokenizer = _tokenizer

    prot_tokens = [tokenizer(p) for p in protein_ids]
    fasta_tokens = [tokenizer(f) for f in fasta_headers]

    # Build vocabulary
    all_tokens = set()
    for toks in prot_tokens + fasta_tokens:
        all_tokens.update(toks)
    all_tokens = sorted(all_tokens)
    token_idx = {tok: i for i, tok in enumerate(all_tokens)}

    # Encode proteins
    row, col = [], []
    for i, toks in enumerate(prot_tokens):
        for t in toks:
            if t in token_idx:
                row.append(i)
                col.append(token_idx[t])
    data = np.ones(len(row), dtype=np.uint8)
    X_prot = csr_matrix((data, (row, col)), shape=(len(protein_ids), len(all_tokens)))

    # Encode FASTA headers
    row, col = [], []
    for i, toks in enumerate(fasta_tokens):
        for t in toks:
            if t in token_idx:
                row.append(i)
                col.append(token_idx[t])
    data = np.ones(len(row), dtype=np.uint8)
    X_fasta = csr_matrix((data, (row, col)), shape=(len(fasta_tokens), len(all_tokens)))

    # Matrix multiplication
    similarity = X_prot @ X_fasta.T  # shape: (n_proteins, n_fasta)

    # Find best match per protein_id
    best_idx = similarity.argmax(axis=1).A.flatten()  # to 1D numpy
    best_score = similarity.max(axis=1).toarray().flatten()

    # matched_fasta = np.array(fasta_headers)[best_idx]
    matched_fasta = np.array(list(fasta_headers))[best_idx]

    result = pd.DataFrame(
        {
            "protein_id": protein_ids,
            "best_fasta_fullname": matched_fasta,  # .flaten(),
            "overlap_score": best_score,  # .flatten()
        }
    )

    # Optional: filter out weak matches
    # if min_score > 1:
    lowmatches = result[result["overlap_score"] < min_score].copy()
    if len(lowmatches) > 0:
        logger.warning(
            "%d low-overlap matches between protein column and FASTA header",
            len(lowmatches),
        )
        result.loc[lowmatches.index, "best_fasta_fullname"] = pd.NA

    fasta_header_set = set(str(h) for h in fasta_headers)
    not_in_fasta = ~result["best_fasta_fullname"].isin(fasta_header_set)
    if not_in_fasta.any():
        missing_df = result[not_in_fasta].copy()
        logger.warning(
            "%d matches not present in FASTA keys", missing_df.shape[0]
        )
        result.loc[not_in_fasta, "best_fasta_fullname"] = pd.NA
        lowmatches = pd.concat([lowmatches, missing_df], ignore_index=True)

    res_dict = result.set_index("protein_id")["best_fasta_fullname"].to_dict()

    return res_dict, lowmatches


from .io_psm import *
