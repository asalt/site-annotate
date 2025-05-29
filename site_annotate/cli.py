import sys
import os
import re
import glob
from functools import partial
from pprint import pprint
import subprocess
import functools
import logging
import pathlib
from pathlib import Path
import click
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm
import pandas as pd
from itertools import product

from rapidfuzz import process, fuzz

import pyfaidx
from pyfaidx import Fasta

import janitor

from . import log
from . import io
from . import io_external
from .io import get_reader, validate_metadata, load_and_validate_files
from . import modisite
from .utils import data_generator
#from .constants import VALID_MODI_COLS, get_all_columns
from .runner import run_pipeline
from . import mapper
from . import reduce
from . import tasks

logger = log.get_logger(__file__)

TEMPLATE_PATH = pathlib.Path(__file__).parent.parent / "scripts"  # 2 levels up
BASE_CONFIG = pathlib.Path(__file__).parent.parent / "config" / "base.toml"


@click.group(chain=True)
def main():
    pass


@main.command()
@click.option("-n", "--name", default="site-annotate.toml")
def get_config(name):
    import shutil

    new_file = pathlib.Path(".").absolute() / name
    if new_file.exists():
        logger.error(f"File {new_file} already exists")
        raise FileExistsError(f"File {new_file} already exists")
    print(f"Writing {new_file} ")
    shutil.copy(BASE_CONFIG, new_file)


def get_templates(TEMPLATE_PATH):
    """
    premade Rmd templates
    """

    if not TEMPLATE_PATH.exists():
        logger.error(f"Template path not found: {TEMPLATE_PATH}")
        raise
    REPORT_TEMPLATES = TEMPLATE_PATH.glob("*Rmd")
    REPORT_TEMPLATES = {x.stem: x for x in REPORT_TEMPLATES}
    # logger.info(REPORT_TEMPLATES)
    return REPORT_TEMPLATES


REPORT_TEMPLATES = get_templates(TEMPLATE_PATH)


def prepare_params(
    template, config, data_dir, output_dir, metadata, gct, root_dir, save_env
):
    params_dict = {}

    # Resolve paths
    if config:
        params_dict["config"] = str(pathlib.Path(config).absolute())
    if gct:
        if not os.path.exists(gct):
            raise FileNotFoundError(f"{gct} does not exist")
        params_dict["gct"] = str(pathlib.Path(gct).absolute())
    if data_dir:
        params_dict["data_dir"] = str(pathlib.Path(data_dir).absolute())
    if output_dir is None:
        output_dir = data_dir or "."
    output_dir = pathlib.Path(output_dir).absolute()
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    params_dict["output_dir"] = str(output_dir)
    params_dict["root_dir"] = str(root_dir)
    params_dict["save_env"] = str(save_env).lower()

    # Handle metadata
    if metadata and not gct:
        metadata_path = pathlib.Path(metadata).absolute()
        validated_meta = validate_meta(metadata_path, data_dir)
        validated_meta_path = metadata_path.parent / (
            metadata_path.stem + "_validated.tsv"
        )
        validated_meta.to_csv(validated_meta_path, sep="\t", index=False)
        params_dict["metadata"] = str(validated_meta_path)
        logger.info(f"Wrote validated metadata: {validated_meta_path}")

    return params_dict


@main.command()
@click.option(
    "-e",
    "--extended",
    is_flag=True,
    default=False,
    show_default=True,
    help="print extended info about templates",
)
def show_templates(extended):
    for k, v in REPORT_TEMPLATES.items():
        print(k, v)
    # end


def common_options(f):
    f = click.option(
        "-c",
        "--config",
        type=click.Path(exists=True, dir_okay=False),
        help=".toml file with additional parameters for report",
    )(f)
    f = click.option(
        "-m",
        "--metadata",
        type=click.Path(exists=True, dir_okay=False),
        help="tsv file with rec run search info used to associate with data",
    )(f)
    f = click.option(
        "-o",
        "--output-dir",
        type=click.Path(exists=False, file_okay=False, dir_okay=True),
        required=False,
        default=None,
        show_default=True,
        help="Output directory, defaults to data_dir if not explicitly set",
    )(f)
    f = click.option(
        "-d",
        "--data-dir",
        type=click.Path(exists=True, file_okay=False, dir_okay=True),
        required=True,
        default=pathlib.Path(".").absolute(),
        show_default=True,
    )(f)
    f = click.option(
        "-g",
        "--gct",
        type=click.Path(exists=True, file_okay=True, dir_okay=True),
        required=False,
        default=None,
        show_default=True,
    )(f)
    f = click.option(
        "-r",
        "--root-dir",
        type=click.Path(exists=True, file_okay=False, dir_okay=True),
        required=True,
        default=pathlib.Path(".").absolute(),
        show_default=True,
        help="Root directory, default is location of command invocation",
    )(f)
    f = click.option(
        "--save-env",
        type=click.BOOL,
        required=False,
        default=False,
        show_default=True,
        help="save r environ to output directory",
    )(f)
    return f


# @click.option(
#     "-d",
#     "--data-dir",
#     type=click.Path(exists=True, file_okay=False, dir_okay=True),
#     required=True,
#     default=pathlib.Path(".").absolute(),
#     show_default=True,
# )
# @click.option(
#     "-o",
#     "--
#     type=click.Path(exists=True, file_okay=False, dir_okay=True),
#     required=False,
#     default=None,
#     show_default=True,
#     help="Output directory, defaults to data_dir if not explictely set",
# )
# @click.option("-m", "--metadata", type=click.Path(exists=True, dir_okay=False))


def conf_to_dataframe(conf_file, set_index=False):
    # this should be easy to test

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


def validate_meta(
    metadata_file: pathlib.Path, data_dir: pathlib.Path, ensure_data=True
):
    """
    top level function to read metadata file and validate it
    """

    reader = io.get_reader(str(metadata_file))
    meta_df = reader(metadata_file, dtype={"label": "str"})
    meta_df = io.validate_metadata(meta_df)  # or .pipe

    rec_run_searches = meta_df.rec_run_search.unique()
    expr_data = io.find_expr_files(rec_run_searches, data_dir)
    expr_data = {k: v for (k, v) in expr_data.items() if v != None}
    if len(expr_data) == 0 and not ensure_data:  # no data
        return

    meta_df_final = io.validate_expr_files(expr_data, meta_df)
    logger.info(f"successfully validated {metadata_file}")

    return meta_df_final


@main.command()
@click.argument("metadata", type=click.Path(exists=True, dir_okay=False))
@click.argument("data_dir", type=click.Path(exists=True, file_okay=False))
def check_meta(metadata, data_dir):
    """
    used to check if metadata tabular file matches with tmt-integrator output file
    used internally in `merge_meta`
    looks for rec_run_search string in metadtaa file (e.g. 12345_1_1) and matching
    files (e.g. 12345_1_1_abundance_single-site_MD.tsv)
    a searchno and runno are not strictly necessary in the metadata, and will take on default values(1, 7) if not provided.
    but critically, the `rec_runno_searchno` pattern string needs to match
    TODO: allow partial matches and let runno and searchno take on default values if not provided
    """
    meta_validated = validate_meta(metadata, data_dir)
    print(meta_validated.to_string())
    print("Success")


@main.command()
@click.option(
    "--result-dir", type=click.Path(exists=False, file_okay=False, dir_okay=True)
)
@click.argument("metadata", type=click.Path(exists=True, dir_okay=False))
@click.argument(
    "data_dir",
    type=click.Path(exists=True, file_okay=False),
    default=click.Path("."),
    # show_default=True,
)
def merge_meta(result_dir, metadata, data_dir):
    r"""
    [DATA_DIR] is the directory where the data files are located

    default is current directory

    Merge expression matrix from fragpipe output:
    starting with this file:
        `tmt-report/abundance_single-site_MD.tsv`
        and a config file with recno 12345
    add the recno prefix and run:
        `tmt-report/12345_abundance_single-site_MD.tsv`
    and run:
        `site-annotate merge-meta ./path/to/meta.tsv ./siteinfo-dir`
    which will produce a gct file with metadata

    Note only works with a single "plex". (not true anymore)
    recno information is only available from the prefix (12345) appended to the front of the emat

    """
    data_dir = pathlib.Path(data_dir).absolute()

    # if result_dir is None:
    result_dir = result_dir or data_dir.parent
    result_dir = pathlib.Path(result_dir).absolute()

    metadata = pathlib.Path(metadata).absolute()
    # can add try/except here to provide better error msg
    meta_validated = validate_meta(metadata, data_dir)

    meta_validated_fname = metadata.parent / (metadata.stem + "_validated.tsv")
    meta_validated.to_csv(meta_validated_fname, sep="\t", index=False)

    outname = result_dir / (metadata.stem + "_siteinfo_combined")
    # ===
    logger.info(f"writing {meta_validated_fname}")
    res = io.merge_metadata(meta_validated, outname)  # writes gct file to output


def match_columns(sample_names: list, df: pd.DataFrame, threshold=97):
    """
    Matches `sample` column values from annotation_df to df.columns
    using fuzzy string matching with preprocessing.

    :param annotation_df: DataFrame containing a "sample" column.
    :param df: DataFrame whose columns need to be matched.
    :param threshold: Similarity threshold for a confident match (default=97%).
    :return: A dictionary mapping `sample` values to best matches in `df.columns`.
    """
    # Extract and preprocess column names and samples to match

    def preprocess_string(s):
        """
        Preprocesses a string by removing spaces, ...,  and normalizing case.
        """
        return re.sub(r"[\s]+", "", s).lower()

    samples = [preprocess_string(sample) for sample in sample_names]
    columns = [preprocess_string(col) for col in df.columns]

    # Dictionary to store matches
    matches = {}

    # Perform fuzzy matching
    for original_sample, processed_sample in zip(sample_names, samples):
        # Get best match using token_sort_ratio for better handling of word order
        match, score, _ = process.extractOne(
            processed_sample, columns, scorer=fuzz.token_sort_ratio
        )

        if (
            score / 100 >= threshold / 100
        ):  # Normalize score to 0-1 range for comparison
            # Find the original column name corresponding to the match
            original_match = df.columns[columns.index(match)]
            matches[original_sample] = original_match
        else:
            logger.warning(
                f"No high-confidence match for '{original_sample}' (Best: '{match}' with {score}%)"
            )
            matches[original_sample] = None  # No confident match

    return matches


@main.command()
@click.option(
    "-e",
    "--experiment-annotation",
    type=click.Path(exists=True, file_okay=True),
    required=True,
    help="Tabular file with columns: plex, channel, sample.",
)
@click.argument(
    "target",
    type=click.Path(exists=True, file_okay=True),
    default=".",
)
def add_experiment_annotation(experiment_annotation, target):
    """
    (experimental)
    For FragPipe TMT-integrator output.
    Merges FragPipe TMT-integrator output with an experiment annotation file.

    [experiment-annotation] is a tabular file generated from FragPipe.

    [TARGET] target file or directory

    :experiment-annotation:
    expected columns are : `plex`, `channel` and `sample`
    other columns such as `condition` and `replicate`, which come from FragPipe, are not used
    information about the "recno" should either be extractable from the `plex` column
    or a separate `recno` column can be manually added

    :target:
    example tmt-report/
    -> reads in abundance_xx_yy.tsv
        xx is a "level" such as gene, signle_site, protein.
            we are primarily interested in gene and signle_site
        we compare all of this to the config file, experiment_annotation)
        and attempt to merge.



    """
    # Validate inputs
    try:
        reader = get_reader(str(experiment_annotation))
        annotation_df = reader(experiment_annotation)
    except Exception as e:
        raise click.ClickException(f"Error reading annotation file: {e}")

    required_columns = {"plex", "channel", "sample"}
    if not required_columns.issubset(annotation_df.columns):
        raise click.ClickException(
            f"Annotation file missing required columns: {required_columns - set(annotation_df.columns)}"
        )

    click.echo(f"Processing target: {target}")
    click.echo(f"Loaded annotation file with {len(annotation_df)} rows.")

    # Placeholder for additional processing
    click.echo("Processing...")

    def search_files(base_dir, middle_part="gene"):
        """
        Searches for all files in the base_dir that match the middle_part pattern.

        :param base_dir: Directory to search in.
        :param middle_part: Middle part of the filename to match.
        :return: List of matching file paths.
        """
        # Define patterns
        prefixes = ["abundance", "ratio"]
        suffixes = ["GN.tsv", "MD.tsv", "None.tsv"]

        # Generate all possible filenames using product
        filenames = [
            f"{prefix}_{middle_part}_{suffix}"
            for prefix, suffix in product(prefixes, suffixes)
        ]

        # Check for file existence
        matching_files = [
            str(Path(base_dir) / filename)
            for filename in filenames
            if (Path(base_dir) / filename).exists()
        ]

        return matching_files

    files = search_files(target)

    for file in files:
        print(file)
        df = get_reader(file)(file)
        cid_cols = set(annotation_df["sample"]) & set(df.columns)
        # these should match
        # now do the merge and save the results
        sample_mapping = match_columns(annotation_df["sample"], df)
        # normalize names
        df = df.rename(columns={v: k for k, v in sample_mapping.items()})

        from .io import write_gct

        write_gct(
            df,  # this is not done yet
            cdesc=annotation_df,
        )


@main.command()
@common_options
@click.option(
    "-t",
    "--template",
    type=click.Choice(REPORT_TEMPLATES.keys()),
    required=False,
    default="analyze_modi",
    show_default=True,
)
@click.option("--interactive", default=False, show_default=True, is_flag=True)
@click.option("--save-env", default=False, show_default=True)
def report(
    template,
    config,
    data_dir,
    output_dir,
    metadata,
    gct,
    root_dir,
    interactive,
    save_env,
    **kwargs,
):

    params_dict = prepare_params(
        template, config, data_dir, output_dir, metadata, gct, root_dir, save_env
    )
    tasks.run_r_code_with_params(params_dict, interactive=interactive)
    return

    print(template)

    template_file = REPORT_TEMPLATES[template]

    params_dict = dict()

    if config is not None:
        config = pathlib.Path(config).absolute()
        params_dict["config"] = str(config)

    if gct is not None:
        if not os.path.exists(gct):
            logger.error(f"{gct} does not exist")
            raise FileNotFoundError(f"{gct} does not exist")
        params_dict["gct"] = pathlib.Path(gct).absolute()

    if data_dir is not None:
        data_dir = pathlib.Path(data_dir).absolute()
        params_dict["data_dir"] = str(data_dir)

    if metadata is not None and gct is None:
        metadata = pathlib.Path(metadata).absolute()
        # now we check if metadata is aligned with available experimental data
        meta_validated = validate_meta(metadata, data_dir)
        meta_validated_fname = metadata.parent / (metadata.stem + "_validated.tsv")
        meta_validated.to_csv(meta_validated_fname, sep="\t", index=False)

        # ===
        logger.info(f"wrote {meta_validated_fname}")
        params_dict["metadata"] = str(meta_validated_fname)
    if output_dir is None:
        output_dir = data_dir
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    output_dir = str(pathlib.Path(output_dir).absolute())
    #

    # Create a dictionary with all parameters
    for k, v in kwargs.items():
        params_dict[k] = v
    params_dict.update(
        {
            # "config": str(config),
            # "data_dir": str(data_dir),
            # "metadata": str(meta_validated_fname),
            "output_dir": str(output_dir),
            "root_dir": str(root_dir),
            "save_env": str.lower(str(save_env)),
        }
    )

    if (
        template == "merge_modis" or template == "heatmap"
    ):  # eventually will move this out
        params = (
            "list(" + ", ".join(f"'{k}' = '{v}'" for k, v in params_dict.items()) + ")"
        )
        cmd = [
            f"Rscript",
            "-e",
            f"""library(rmarkdown)
            rmarkdown::render("{template_file}",
            output_dir="{output_dir}",
            params={params},
            )""",
        ]
        logger.info(cmd)
        subprocess.run(cmd)
        return

    import rpy2.robjects as robjects

    for k, v in params_dict.items():
        robjects.r.assign(k, str(v))

    runfile = pathlib.Path(__file__).parent.parent / "R/run.R"
    robjects.r.assign("here_dir", str(runfile.parent))
    robjects.r("setwd(here_dir)")
    robjects.r.source(str(runfile))

    os.chdir(str(runfile.parent))

    robjects.r(
        """
        data_dir <- ifelse(exists("data_dir"), data_dir, '.')
        output_dir <- ifelse(exists("output_dir"), output_dir, '.')
        config_file <- ifelse(exists("config"), config, NA) # have to make it NA not null
        gct_file <- ifelse(exists("gct"), gct, NA) # have to make it NA not null
        save_env <- ifelse(exists("save_env"), save_env, FALSE)

        print(paste0('output dir is: ', output_dir))
        run(
            data_dir = data_dir,
            output_dir = output_dir,
            config_file = config_file,
            gct_file = gct_file,
            save_env = save_env,
               )
               """
    )
    return


# we are not using this
# can probably remove
@main.command()
@click.option(
    "-o",
    "--output-dir",
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
    required=False,
    default=None,
    show_default=True,
    help="Output directory, defaults to data_dir if not explicitly set",
)
@click.option(
    "-d",
    "--data-dir",
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
    required=True,
    default=pathlib.Path(".").absolute(),
    show_default=True,
)
@click.option("--modi-abbrev", default="p", show_default=True)
def reduce_sites(output_dir, data_dir, modi_abbrev, **kwargs):
    raise ValueError("depreciated")

    template_file = REPORT_TEMPLATES["site_annotate_post_not_reduced"]

    data_dir = pathlib.Path(data_dir).absolute()
    if output_dir is None:
        output_dir = data_dir

    params_dict = {
        "data_dir": str(data_dir),
        "modi_abbrev": modi_abbrev,
    }
    params = "list(" + ", ".join(f"'{k}' = '{v}'" for k, v in params_dict.items()) + ")"
    cmd = [
        f"Rscript",
        "-e",
        f"""library(rmarkdown)
        rmarkdown::render("{template_file}",
        output_dir="{output_dir}",
        params={params},
        )""",
    ]

    logger.info(cmd)
    subprocess.run(cmd)


@main.command()
@click.option("--cores", default=1, show_default=True)
@click.option(
    "-p",
    "--psms",
    type=click.Path(exists=True, dir_okay=False),
    multiple=True,
)
@click.option(
    "-o",
    "--output-dir",
    type=click.Path(exists=True, file_okay=False),
    default=None,
    show_default=True,
)
@click.option(
    "--prefix",
    type=str,
    default=None,
    show_default=True,
    help="prefix to append to beginning of files",
)
@click.option(
    "--uniprot-check/--no-uniprot-check", default=False, show_default=True, is_flag=True
)
@click.option(
    "-f", "--fasta", type=click.Path(exists=True, dir_okay=False), help="fasta file"
)
def run(cores, psms, output_dir, prefix, uniprot_check, fasta):
    """
    this is the main code for compiling site-level expression matrix from from PSMs table from FragPipe.
    for each modi of interest there should be a column of form `aas_mass_shift`.
    So for phospho a column of form `sty_79_9663`

    critical columns include:
        - spectrum
        - sty_79_9663
        - sty_79_9663_best_localization
    ...
    """
    if not psms:
        logger.warning("Warning: No PSM files provided.")
        return
    if len(psms) > 1:
        logger.error("more than 1 psms not supported yet")
        raise NotImplementedError("more than 1 psms not supported yet")

    df, fasta_data, fa_psp_ref = load_and_validate_files(psms[0], fasta, uniprot_check)
    fullres = run_pipeline(df, fasta_data, fa_psp_ref, cores)
    # type fullres is list of dicts
    # fullres[0].keys()
    # dict_keys(['sty_79_9663'])
    # fullres[0]['sty_79_9663'].columns
    # Index(['position_relative', 'AA', 'prob', 'original_index', 'spectrum', 'spectrum_file', 'peptide', 'modified_peptide',
    #     'extended_peptide', 'prev_aa', 'next_aa', 'peptide_length', 'charge', 'retention', 'observed_mass', 'calibrated_observed_mass',
    #     'observed_m_z', 'calibrated_observed_m_z', 'calculated_peptide_mass', 'calculated_m_z', 'delta_mass', 'spectralsim', 'rtscore',
    #     'expectation', 'hyperscore', 'nextscore', 'peptideprophet_probability', 'number_of_enzymatic_termini',
    #     'number_of_missed_cleavages', 'protein_start', 'protein_end', 'intensity', 'assigned_modifications', 'observed_modifications',
    #     'm_15_9949', 'm_15_9949_best_localization', 'sty_79_9663', 'sty_79_9663_best_localization', 'purity', 'is_unique', 'protein',
    #     'protein_id', 'entry_name', 'gene', 'protein_description', 'mapped_genes', 'mapped_proteins', 'TMT_126', 'TMT_127_N',
    #     'TMT_127_C', 'TMT_128_N', 'TMT_128_C', 'TMT_129_N', 'TMT_129_C', 'TMT_130_N', 'TMT_130_C', 'TMT_131_N', 'TMT_131_C',
    #     'TMT_132_N', 'TMT_132_C', 'TMT_133_N', 'TMT_133_C', 'TMT_134_N', 'position_absolut', 'fifteenmer', 'protein_length', 'tmt_sum',
    #     'TMT_126_ratio', 'TMT_127_N_ratio', 'TMT_127_C_ratio', 'TMT_128_N_ratio', 'TMT_128_C_ratio', 'TMT_129_N_ratio',
    #     'TMT_129_C_ratio', 'TMT_130_N_ratio', 'TMT_130_C_ratio', 'TMT_131_N_ratio', 'TMT_131_C_ratio', 'TMT_132_N_ratio',
    #     'TMT_132_C_ratio', 'TMT_133_N_ratio', 'TMT_133_C_ratio', 'TMT_134_N_ratio', 'TMT_126_intensity', 'TMT_127_N_intensity',
    #     'TMT_127_C_intensity', 'TMT_128_N_intensity', 'TMT_128_C_intensity', 'TMT_129_N_intensity', 'TMT_129_C_intensity',
    #     'TMT_130_N_intensity', 'TMT_130_C_intensity', 'TMT_131_N_intensity', 'TMT_131_C_intensity', 'TMT_132_N_intensity',
    #     'TMT_132_C_intensity', 'TMT_133_N_intensity', 'TMT_133_C_intensity', 'TMT_134_N_intensity', 'highest_prob'],
    #     dtype='object')
    finalres = process_results(fullres, decoy_label="rev_")
    # finalres is a dict with keys modi ( sty_79_9663, k_42_0106, ... )
    # and values the concatenated dataframe of all sites, along with tmt quant if applicable
    # this "not reduced" file can be huge ( 15g + )
    # save_results(finalres, psms[0])

    site_reduced = reduce.reduce_sites(
        finalres,
        decoy_label="rev_",
    )

    save_results(site_reduced, psms[0], name="site_annotation_reduced", prefix=prefix)

    site_reduced_mapped = mapper.add_annotations(site_reduced)
    save_results(site_reduced_mapped, psms[0], name="site_annotation_reduced_mapped", prefix=prefix)

    return



def process_results(fullres, decoy_label="rev_"):
    """
    fullres is a list of dicts
    """
    finalres = {}
    all_frames = defaultdict(list)
    for result_dict in fullres:
        for modi_id, df in result_dict.items():
            all_frames[modi_id].append(df)
        # for col in VALID_MODI_COLS:
        #     frames = [items.get(col) for items in fullres if items.get(col) is not None]
        # if frames:
        #     finalres[col] = pd.concat(frames, ignore_index=True)
        # else:
        #     pass  # not all need to be present
        # logger.warning(f"No results returned for {col}")
    for k, v in all_frames.items():
        _df = pd.concat(v)  # filter decoys here
        _df = _df[~_df.protein.str.startswith(decoy_label)]
        finalres[k] = _df
    return finalres


def save_results(finalres, input_file, name="site_annotation_notreduced", prefix=None):
    if prefix is None:
        prefix = ""
    infile = pathlib.Path(input_file)
    for key, val in finalres.items():
        # outfile = infile.parent / f"{key}_site_annotation_notreduced.tsv"
        outfile = infile.parent / f"{prefix}{key}_{name}.tsv"
        logger.info(f"Writing {outfile}")
        val.to_csv(outfile, sep="\t", index=False)
