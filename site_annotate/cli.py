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

import janitor

from . import log
from . import io
from . import io_external
from . import modisite
from .utils import data_generator
from .constants import VALID_MODI_COLS
from .runner import run_pipeline
from . import mapper

logger = log.get_logger(__file__)


@click.group(chain=True)
def main():
    pass


TEMPLATE_PATH = pathlib.Path(__file__).parent.parent / "scripts"  # 2 levels up

if not TEMPLATE_PATH.exists():
    logger.error(f"Template path not found: {TEMPLATE_PATH}")
    raise
REPORT_TEMPLATES = TEMPLATE_PATH.glob("*Rmd")
REPORT_TEMPLATES = {x.stem: x for x in REPORT_TEMPLATES}
# logger.info(REPORT_TEMPLATES)


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


# import ipdb
# ipdb.set_trace()


@main.command()
@click.option("--cores", default=1, show_default=True)
@click.option(
    "-p", "--psms", type=click.Path(exists=True, dir_okay=False), multiple=True
)
@click.option(
    "-o",
    "--output-dir",
    type=click.Path(exists=True, file_okay=False),
    default=None,
    show_default=True,
)
@click.option(
    "-f", "--fasta", type=click.Path(exists=True, dir_okay=False), help="fasta file"
)
def run(cores, psms, output_dir, fasta, **kwargs):

    DECOY_FLAG = "rev_"

    if len(psms) == 0:
        logging.warning("Warning: No PSM files provided")
        return

    # TODO expand for all psms
    logger.info(f"loading {psms[0]}")
    df = io.read_psm_file(psms[0])
    try:
        io.validate_psm_file(df)
    except Exception as e:
        raise e

    df = mapper.add_uniprot(df)  #

    logger.info(f"loading {fasta}")
    # fa = io.read_fasta(fasta)
    import pyfastx

    fa = pyfastx.Fasta(fasta)
    fa_psp_ref = io_external.read_psite_fasta()

    # g2 = df2.groupby("protein")
    # g3 = df3.group_by("protein")

    # breakpoint()

    # def process_frame(key_frame, fa, modisite):
    #     key, frame = key_frame
    #     # if 'Cont' not in key:
    #     #     print()
    #     #     import ipdb; ipdb.set_trace()
    #     # subfa = fa[(fa["id"] == key) & (~fa["id"].str.startswith(DECOY_FLAG))]
    #     try:
    #         subfa = fa[key]
    #     except KeyError:
    #         logger.info(f"skipping key {key}, not found in fasta")
    #         return
    #     seqinfo = io.extract_info_from_header(subfa.name)
    #     seqinfo["sequence"] = subfa.seq
    #     if len(subfa) == 0:
    #         logger.info(f"skipping key {key}, db mismatch")
    #         return
    #     # seqinfo = subfa.iloc[0].to_dict()
    #     res = modisite.main(frame, seqinfo)
    #     return res

    # fullres = list()

    # if cores == 1:
    #     for item in tqdm(g):
    #         res = process_frame(item, fa, modisite)
    #         fullres.append(res)

    # if cores > 1:
    #     with ThreadPoolExecutor() as executor:
    #         futures = {
    #             executor.submit(process_frame, item, fa, modisite): item for item in g
    #         }

    #         for future in tqdm(
    #             as_completed(futures),
    #             total=len(futures),
    #             mininterval=0.4,
    #             smoothing=0.1,
    #         ):
    #             try:
    #                 result = future.result()
    #                 fullres.append(result)
    #             except Exception as e:

    #                 # print(f"An error occurred: {e}")
    #                 pass

    # g = df.groupby("protein")
    # fullres = run_pipeline(g, fa, cores=cores)

    fullres = run_pipeline(df, fa, fa_psp_ref=fa_psp_ref, cores=cores)

    finalres = dict()
    for col in VALID_MODI_COLS:  # of form ["sty_79_9663", "k_42_0106", ...]
        frames = list()
        for items in fullres:
            vals = items.get(col)
            if vals is None:
                continue
            # vals = [x for x in vals if x is not None]
            if len(vals) == 0:
                logger.warning(f"no results returned for {col}")
                continue
            frames.append(vals)
        if len(frames) != 0:
            finalres[col] = pd.concat(frames, ignore_index=True)
    # res1 = [x.get("sty_79_9663") for x in fullres if x is not None]
    # res1 = [x for x in res1 if x is not None]
    # if len(res1) == 0:
    #     logger.error("no results returned")
    #     sys.exit(1)

    # fullres_df = pd.concat(res1, ignore_index=True)

    infile = pathlib.Path(psms[0])
    for key, val in finalres.items():
        outfile = infile.parent / f"{key}_site_annotation_notreduced.tsv"
        logger.info(f"writing {outfile}")
        val.to_csv(outfile, sep="\t", index=False)

    # with ProcessPoolExecutor() as executor:
    #     futures = {executor.submit(process_row, row): row for index, row in df.iterrows()}
    #     for future in as_completed(futures):
    #         peptide_to_proteins = future.result()
    #         for peptide, proteins in peptide_to_proteins.items():
    #             all_peptide_to_proteins[peptide].extend(proteins)


# # not using, for later
# def process_frame(key_frame, fa, modisite):
#     key, frame = key_frame
#     subfa = fa[fa["id"] == key]
#     assert len(subfa) == 1
#     seqinfo = subfa.iloc[0].to_dict()
#     res = modisite.main(frame, seqinfo)
#     return res


# def main_function(g, fa, modisite):
#     fullres = []

#     with ThreadPoolExecutor() as executor:
#         futures = {
#             executor.submit(process_frame, item, fa, modisite): item for item in g
#         }

#         for future in tqdm(as_completed(futures), total=len(futures)):
#             try:
#                 result = future.result()
#                 fullres.append(result)
#             except Exception as e:
#                 print(f"An error occurred: {e}")

#     return fullres
