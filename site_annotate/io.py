# io.py
import os
import re
import glob
import pathlib
from pathlib import Path
import pandas as pd
import janitor
from Bio import SeqIO

from .utils import data_generator

from . import log
from .constants import VALID_MODI_COLS, POSSIBLE_SITE_ID_COLS

logger = log.get_logger(__file__)

RENAME = {
    "sample_01": "TMT_126",
    "sample_02": "TMT_127_N",
    "sample_03": "TMT_127_C",
    "sample_04": "TMT_128_N",
    "sample_05": "TMT_128_C",
    "sample_06": "TMT_129_N",
    "sample_07": "TMT_129_C",
    "sample_08": "TMT_130_N",
    "sample_09": "TMT_130_C",
    "sample_10": "TMT_131_N",
    "sample_11": "TMT_131_C",
    "sample_12": "TMT_132_N",
    "sample_13": "TMT_132_C",
    "sample_14": "TMT_133_N",
    "sample_15": "TMT_133_C",
    "sample_16": "TMT_134_N",
    "sample_17": "TMT_134_C",
    "sample_18": "TMT_135_N",
    "ion_126_128": "TMT_126",
    "ion_127_125": "TMT_127_N",
    "ion_127_131": "TMT_127_C",
    "ion_128_128": "TMT_128_N",
    "ion_128_134": "TMT_128_C",
    "ion_129_131": "TMT_129_N",
    "ion_129_138": "TMT_129_C",
    "ion_130_135": "TMT_123_N",
    "ion_130_141": "TMT_130_C",
    "ion_131_138": "TMT_131_N",
}
RENAME.update({x.replace("_", "-"): y for x, y in RENAME.items()})
RENAME.update({x: y + "_intensity" for x, y in RENAME.items()})


def set_data_dir():
    # Path of the current script
    current_script_path = Path(__file__).resolve()
    # Path to the top-level data directory
    res = current_script_path.parent.parent / "data"
    logger.debug(f"data dir set to {res}")
    return res


data_dir = set_data_dir()


def conf_to_dataframe(conf_file, set_index=False):
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


def get_reader(file: str):
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
        return None


def convert_tmt_label(shorthand):
    # Normalize the input by removing 'TMT' or 'TMT_' if present
    normalized_input = re.sub(r"TMT[_]?", "", shorthand)

    # Conditional cases equivalent to R's case_when
    if normalized_input == "126":
        return "TMT_126"
    elif normalized_input == "131":
        return "TMT_131_N"
    elif normalized_input == "134":
        return "TMT_134_N"
    elif re.search(r"N$", normalized_input):
        return re.sub(r"(\d+)_?(N)", r"TMT_\1_N", normalized_input)
    elif re.search(r"C$", normalized_input):
        return re.sub(r"(\d+)_?(C)", r"TMT_\1_C", normalized_input)
    else:
        return f"TMT_{normalized_input}"


def find_expr_file(rec_run_search: str, data_dir):
    """
    data_dir should be absolute path by this point
    """
    search_pattern = os.path.join(data_dir, f"{rec_run_search}*reduced*tsv")
    results = glob.glob(search_pattern)

    if not results:  # try one more time
        logger.info(f"trying again with not reduced")
        search_pattern = os.path.join(data_dir, f"{rec_run_search}*tsv")
        results = glob.glob(search_pattern)

    if len(results) == 0:
        logger.warning(f"no files found for {rec_run_search} in {data_dir}")
        raise FileNotFoundError(f"no files found for {rec_run_search} in {data_dir}")
    if len(results) > 1:
        logger.warning(
            f"Ambiguous, found multiple files for {rec_run_search}, {str.join(', ', results)}"
        )
        raise FileNotFoundError(
            f"Ambiguous, found multiple files for {rec_run_search}, {str.join(', ', results)}"
        )
    if len(results) == 1:
        return results[0]


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
        df["label"] = df["label"].apply(convert_tmt_label)
        # can add a check here to assert unique by label and rec_run_search

    if "name" not in df.columns:
        logger.info("name not found, using rec_run_search")
        df["name"] = df["rec_run_search"]
    df.index = df["name"]

    return df


def get_isoform_hierarchy() -> pd.DataFrame:
    target1 = data_dir / "GENCODE.M32.basic.CHR.protein.selection.mapping.txt"
    target2 = data_dir / "GENCODE.V42.basic.CHR.isoform.selection.mapping.txt"
    logger.info(f"Reading {target1}")
    df1 = pd.read_table(target1, sep="\t", low_memory=False)
    logger.info(f"Reading {target2}")
    df2 = pd.read_table(target2, sep="\t", low_memory=False)
    df = pd.concat([df1, df2])
    df = janitor.clean_names(df)
    return df


def prepare_psm_file(df: pd.DataFrame) -> pd.DataFrame:
    """Check if a DataFrame is a valid PSM file.
    this makes use of the proteins and mapped_proteins columns in msfragger output
    """
    # Check if a DataFrame is a valid PSM file
    required_cols = ["peptide", "intensity", "protein", "mapped_proteins"]

    if "intensity" not in df.columns:
        for _x in ("area", "peakarea", "peak_area"):
            if _x in df.columns:
                df["intensity"] = df[_x]
    if "mapped_proteins" not in df.columns:
        if "alternative_proteins" in df.columns:
            df["mapped_proteins"] = df["alternative_proteins"]

    for required_col in required_cols:
        if required_col not in df.columns:
            raise ValueError(f"Invalid PSM file, missing {required_col} column")

    # if not any(df.columns.str.startswith("TMT")):
    df = df.rename(columns=RENAME)
    df["mapped_proteins"] = (
        df["protein"] + ", " + df["mapped_proteins"].fillna("").str.replace("@@", ", ")
    )  # mapped_proteins column does not contain the value in the protein column so we add it here
    df["mapped_proteins2"] = df["mapped_proteins"].apply(
        lambda x: x.split(", ") if isinstance(x, str) else []
    )
    df["mapped_proteins2"] = df["mapped_proteins2"].apply(
        lambda x: list(filter(None, set(x)))
    )

    # # Step 2: Explode the 'mapped_proteins' column
    # df_exploded = df[[*required_cols, "mapped_proteins2"]].explode("mapped_proteins2")
    df_exploded = df.explode("mapped_proteins2")
    df_exploded = df_exploded.reset_index(drop=True)
    df_exploded["protein"] = df_exploded["mapped_proteins2"]

    # return df
    return df_exploded


def read_psm_file(psm_file: str | pathlib.Path) -> pd.DataFrame:
    """Read a PSM file into a DataFrame."""
    # Read a PSM file into a DataFrame
    df = pd.read_csv(psm_file, sep="\t")
    df = janitor.clean_names(df)
    # df = df[df.protein.str.contains("Tcf12")]
    df = prepare_psm_file(df)
    validate_psm_file(df)
    return df


def rename_columns(df: pd.DataFrame, mapping: dict) -> pd.DataFrame:
    """Rename columns in a DataFrame."""
    # Rename columns in a DataFrame
    df = df.rename(columns=mapping)
    return df


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
    """ """
    for rrs, expr_file in rec_run_searches.items():
        rec, run, search = rrs.split("_")

        _meta = meta_df[
            (
                (meta_df.recno == rec)
                & (meta_df.runno == run)
                & (meta_df.searchno == search)
            )
        ]
        meta_df.loc[_meta.index, "expr_col"] = None
        meta_df.loc[_meta.index, "expr_file"] = expr_file
        df = pd.read_table(expr_file, nrows=5)
        df = df.rename(columns=RENAME)

        # if "intensity_sum" not in df.columns:
        #     raise ValueError(f"`intensity_sum` not found in {expr_file}")
        if "intensity_sum" not in df.columns:
            logger.warning(
                f"`intensity_sum` not found in {expr_file}. Downstream analyis might fail"
            )
            # raise ValueError(f"`intensity_sum` not found in {expr_file}")

        if "label" in _meta.columns:
            for label in _meta["label"].tolist():
                label_mapping = [x for x in df.columns if x.startswith(label)]
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
                    import ipdb

                    ipdb.set_trace()
                    1 + 1
                ix = ix[0]
                meta_df.loc[ix, "expr_col"] = expr_col
                meta_df.loc[ix, "expr_file"] = expr_file
        else:
            expr_col = "intensity_sum"
            meta_df.loc[_meta.index, "expr_col"] = "intensity_sum"
    return meta_df


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


def merge_metadata(metadata: pd.DataFrame, filepath="siteinfo_combined"):
    # Load metadata file

    # Find unique expression file paths
    unique_expr_files = metadata["expr_file"].unique()

    # Dictionary to hold processed DataFrames for each file
    processed_data = {}
    # combined_rdesc = None

    # Iterate over each unique expression file
    for expr_file in unique_expr_files:  # this is messy, begin to refactor
        # Load the expression data
        expr_data = pd.read_csv(
            expr_file, sep="\t"
        )  # Assuming TSV format for expr files which is true
        expr_data = expr_data.rename(
            columns=RENAME
        )  # a universal renamer for various possible names
        if "site_id" not in expr_data.columns:
            if "sitename2" in expr_data:
                site_id = (
                    expr_data[["sitename2", "ENSP", "fifteenmer"]]
                    .fillna("")
                    .agg("_".join, axis=1)
                )
            else:
                site_id = (
                    expr_data[["sitename", "ENSP", "fifteenmer"]]
                    .fillna("")
                    .agg("_".join, axis=1)
                )
            site_id = make_unique(site_id)
            expr_data["site_id"] = site_id
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
        # import ipdb; ipdb.set_trace()
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

    # import ipdb; ipdb.set_trace()

    # import ipdb; ipdb.set_trace()
    write_gct(combined_emat, cdesc=metadata, rdesc=combined_rdesc, filename=filepath)

    combined_df = combined_rdesc.join(combined_emat)

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
    header_lines = [
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
