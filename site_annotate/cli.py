import sys
import os
import re
import glob
import subprocess
import functools
import logging
import pathlib
import click
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed

import pyfaidx
from pyfaidx import Fasta

import janitor

from . import log
from . import io
from . import io_external
from . import modisite
from .utils import data_generator
from .constants import VALID_MODI_COLS, get_all_columns
from .runner import run_pipeline
from . import mapper
from . import reduce

logger = log.get_logger(__file__)


@click.group(chain=True)
def main():
    pass


TEMPLATE_PATH = pathlib.Path(__file__).parent.parent / "scripts"  # 2 levels up


def get_templates(TEMPLATE_PATH):

    if not TEMPLATE_PATH.exists():
        logger.error(f"Template path not found: {TEMPLATE_PATH}")
        raise
    REPORT_TEMPLATES = TEMPLATE_PATH.glob("*Rmd")
    REPORT_TEMPLATES = {x.stem: x for x in REPORT_TEMPLATES}
    # logger.info(REPORT_TEMPLATES)
    return REPORT_TEMPLATES


REPORT_TEMPLATES = get_templates(TEMPLATE_PATH)


@main.command()
@click.option("-e", "--extended", is_flag=True, default=False, show_default=True, help="print extended info about templates")
def show_templates(extended):
    for k, v in REPORT_TEMPLATES.items():
        print(k, v)


def common_options(f):
    f = click.option("-m", "--metadata", type=click.Path(exists=True, dir_okay=False))(
        f
    )
    f = click.option(
        "-o",
        "--output-dir",
        type=click.Path(exists=True, file_okay=False, dir_okay=True),
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

def get_reader(file:str ):
    if file.endswith('.tsv'):
        return pd.read_table
    elif file.endswith('.csv'):
        return pd.read_csv
    elif file.endswith('.xlsx'):
        return pd.read_excel
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


def validate_metadata(df: pd.DataFrame):
    if "recno" not in df.columns:
        logger.error(f"recno not present in metadata file")
        raise ValueError(f"recno not present in metadata file")
    if "runno"  not in df.columns: # assume 1
        logger.info("runno not found, assuming 1")
        df['runno'] = 1
    if "searchno"  not in df.columns: # assume 1
        logger.info("searchno not found, assuming 7")
        df['searchno'] = 7

    for x in ("recno", "runno", "searchno"):
        df[x] = df[x].astype(str)

    df['rec_run_search'] = df.apply(lambda x: f"{x.recno}_{x.runno}_{x.searchno}", axis=1)

    if "label" in df.columns:
        df["label"] = df["label"].apply(convert_tmt_label)
        # can add a check here to assert unique by label and rec_run_search

    return(df)

def find_expr_file(rec_run_search:str , data_dir):
    search_pattern = os.path.join(data_dir, f"{rec_run_search}*reduced*tsv")
    results = glob.glob(search_pattern)
    if len(results) == 0:
        logger.warning(f"no files found for {rec_run_search} in {data_dir}")
    if len(results) > 1:
        logger.warning(f"Ambiguous, found multiple files for {rec_run_search}, {str.join(', ', results)}")
    if len(results) == 1:
        return results[0]

def find_expr_files(rec_run_searches, data_dir):
    return { rrs: find_expr_file(rrs, data_dir) for rrs in rec_run_searches }

def validate_expr_files(rec_run_searches: dict, meta_df: pd.DataFrame):
    for (rrs, expr_file) in rec_run_searches.items():
        rec, run, search = rrs.split("_")

        _meta = meta_df[
                (
                    (meta_df.recno == rec) &
                    (meta_df.runno == run) & 
                    (meta_df.searchno == search)
                )]
        meta_df.loc[_meta.index, "expr_col"] = None
        meta_df.loc[_meta.index, "expr_file"] = expr_file
        _df = pd.read_table(expr_file, nrows=5)

        if "intensity_sum" not in _df.columns:
            raise ValueError(f"`intensity_sum` not found in {expr_file}")

        if "label" in _meta.columns:
            for label in _meta['label'].tolist():
                label_mapping = [ x for x in _df.columns if x.startswith(label) ]
                if len(label_mapping) == 0:
                    logger.warning(f"could not find sample with label {label}")
                    next
                if len(label_mapping) > 1:
                    logger.warning(f"too many results for {label}")
                    next
                if len(label_mapping) == 1:
                    expr_col = label_mapping[0]

                ix = meta_df[(meta_df['rec_run_search'] == rrs) & (meta_df['label'] == label)].index[0]
                meta_df.loc[ix, "expr_col"] = expr_col
                meta_df.loc[ix, "expr_file"] = expr_file
        else:
            expr_col = "intensity_sum"
            meta_df.loc[_meta.index, "expr_col"] = "intensity_sum"
        return meta_df


def validate_meta(metadata_file: pathlib.Path, data_dir: pathlib.Path):

    reader = get_reader(str(metadata_file))
    meta_df = reader(metadata_file) 
    meta_df = validate_metadata(meta_df) # or .pipe

    rec_run_searches = meta_df.rec_run_search.unique()
    expr_data = find_expr_files(rec_run_searches, data_dir)
    expr_data = { k:v for (k,v) in expr_data.items() if v != None }
    if len(expr_data) == 0: # no data
        return

    meta_df_final = validate_expr_files(expr_data, meta_df)
    return meta_df_final


@main.command()
@common_options
@click.option(
    "-t",
    "--template",
    type=click.Choice(REPORT_TEMPLATES.keys()),
    required=False,
    default="modi",
    show_default=True,
)
def report(template, data_dir, output_dir, metadata, **kwargs):
    print(template)
    template_file = REPORT_TEMPLATES[template]

    data_dir = pathlib.Path(data_dir).absolute()
    metadata = pathlib.Path(metadata).absolute()
    if output_dir is None:
        output_dir = data_dir

    # now we check if metadata is aligned with available experimental data
    # == move all of this to a separate function
    meta_validated = validate_meta(metadata, data_dir)
    meta_validated_fname = metadata.parent / (metadata.stem + "_validated.tsv")
    meta_validated.to_csv(meta_validated_fname, sep='\t', index=False)
    logger.info(f"wrote {meta_validated_fname}")
    #

    # Create a dictionary with all parameters
    params_dict = {
        "data_dir": str(data_dir),
        # "output_dir": str(output_dir),
        "metadata": str(meta_validated_fname),
    }
    for k, v in kwargs.items():
        params_dict[k] = v
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
    "--uniprot-check/--no-uniprot-check", default=False, show_default=True, is_flag=True
)
@click.option(
    "-f", "--fasta", type=click.Path(exists=True, dir_okay=False), help="fasta file"
)
def run(cores, psms, output_dir, uniprot_check, fasta):
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

    save_results(site_reduced, psms[0], name="site_annotation_reduced")

    site_reduced_mapped = mapper.add_annotations(site_reduced)
    save_results(site_reduced_mapped, psms[0], name="site_annotation_reduced_mapped")

    return


def load_psite_fasta():
    try:
        # Call the function that reads the FASTA file
        fa_psp_ref = io_external.read_psite_fasta()
        return fa_psp_ref
    except FileNotFoundError:
        logger.warning("PhosphositePlus fasta not found")
        return None


def load_and_validate_files(psm_path, fasta_path, uniprot_check):

    # try load phosphositeplus fasta
    # needs to be downloaded from https://www.phosphosite.org/staticDownloads manually (free non commercial)
    # fa_psp_ref = load_psite_fasta()

    # logger.info(f"Loading {psm_path}")
    # df = io.read_psm_file(psm_path)

    # fasta_data = Fasta(fasta_path)

    with ThreadPoolExecutor(
        max_workers=3
    ) as executor:  # there's some significant postprocessing these funcs do that makes this worth it, I think
        # need to time test this with larger files
        # Submit the function to the executor
        psp_future = executor.submit(load_psite_fasta)
        df_future = executor.submit(io.read_psm_file, psm_path)
        fasta_future = executor.submit(Fasta, fasta_path)
        #

        fa_psp_ref = psp_future.result()

        logger.info(f"Loading {psm_path}")
        df = df_future.result()

        logger.info(f"Loading {fasta_path}")
        fasta_data = fasta_future.result()

    if fa_psp_ref:
        logger.info("FASTA data loaded successfully.")
    else:
        logger.info("Failed to load FASTA data.")

    df = mapper.extract_keyvals_pipedsep(df)
    if uniprot_check:
        df = mapper.add_uniprot(df)

    return df, fasta_data, fa_psp_ref  #


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


def save_results(finalres, input_file, name="site_annotation_notreduced"):
    infile = pathlib.Path(input_file)
    for key, val in finalres.items():
        # outfile = infile.parent / f"{key}_site_annotation_notreduced.tsv"
        outfile = infile.parent / f"{key}_{name}.tsv"
        logger.info(f"Writing {outfile}")
        val.to_csv(outfile, sep="\t", index=False)
