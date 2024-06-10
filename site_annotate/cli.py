import sys
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
from .constants import VALID_MODI_COLS
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
def show_templates():
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


@main.command()
@common_options
@click.option(
    "-t",
    "--template",
    type=click.Choice(REPORT_TEMPLATES.keys()),
    required=False,
    default="site_annotate_post_not_reduced",
    show_default=True,
)
def report(template, data_dir, output_dir, metadata, **kwargs):
    print(template)
    template_file = REPORT_TEMPLATES[template]

    data_dir = pathlib.Path(data_dir).absolute()
    if output_dir is None:
        output_dir = data_dir

    # Create a dictionary with all parameters
    params_dict = {
        "data_dir": str(data_dir),
        # "output_dir": str(output_dir),
        "metadata": str(metadata),
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
    finalres = process_results(fullres)
    # finalres is a dict with keys modi ( sty_79_9663, k_42_0106, ... )
    # and values the concatenated dataframe of all sites, along with tmt quant if applicable
    save_results(finalres, psms[0])

    site_reduced = reduce.reduce_sites(finalres)
    save_results(site_reduced, psms[0], name="site_annotation_reduced")


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

    if uniprot_check:
        df = mapper.add_uniprot(df)

    return df, fasta_data, fa_psp_ref  #


def process_results(fullres):
    finalres = {}
    for col in VALID_MODI_COLS:
        frames = [items.get(col) for items in fullres if items.get(col) is not None]
        if frames:
            finalres[col] = pd.concat(frames, ignore_index=True)
        else:
            pass  # not all need to be present
            # logger.warning(f"No results returned for {col}")
    return finalres


def save_results(finalres, input_file, name="site_annotation_notreduced"):
    infile = pathlib.Path(input_file)
    for key, val in finalres.items():
        # outfile = infile.parent / f"{key}_site_annotation_notreduced.tsv"
        outfile = infile.parent / f"{key}_{name}.tsv"
        logger.info(f"Writing {outfile}")
        val.to_csv(outfile, sep="\t", index=False)
