import sys
import subprocess
import logging
import pathlib
import click
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed

from . import log
from . import io
from . import modisite

logger = log.get_logger(__file__)


@click.group(chain=True)
def main():
    pass


TEMPLATE_PATH = pathlib.Path(__file__).parent.parent / "scripts"  # 2 levels up
if not TEMPLATE_PATH.exists():
    logger.error(f"Template path not found: {TEMPLATE_PATH}")
    sys.exit(1)
REPORT_TEMPLATES = TEMPLATE_PATH.glob("*Rmd")
REPORT_TEMPLATES = {x.stem: x for x in REPORT_TEMPLATES}
# logger.info(REPORT_TEMPLATES)


@main.command()
@click.option(
    "-t",
    "--template",
    type=click.Choice(REPORT_TEMPLATES.keys()),
    required=False,
    default="site_annotate_post_not_reduced",
    show_default=True,
)
@click.option(
    "-d",
    "--data-dir",
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
    required=True,
    default=pathlib.Path(".").absolute(),
    show_default=True,
)
@click.option(
    "-o",
    "--output-dir",
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
    required=True,
    default=pathlib.Path(".").absolute(),
    show_default=True,
)
@click.option("-m", "--metadata", type=click.Path(exists=True, dir_okay=False))
def report(template, data_dir, output_dir, metadata):
    print(template)
    template_file = REPORT_TEMPLATES[template]

    # Create a dictionary with all parameters
    params_dict = {
        "data_dir": str(data_dir),
        # "output_dir": str(output_dir),
        "metadata": str(metadata),
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
    "-p", "--psms", type=click.Path(exists=True, dir_okay=False), multiple=True
)
@click.option("-o", "--out", default=None)
@click.option(
    "-f", "--fasta", type=click.Path(exists=True, dir_okay=False), help="fasta file"
)
def run(cores, psms, out, fasta, **kwargs):

    DECOY_FLAG = "rev_"

    if len(psms) == 0:
        logging.warning("Warning: No PSM files provided")
        return

    # TODO expand for all psms
    logger.info(f"loading {psms[0]}")
    df = io.read_psm_file(psms[0])

    logger.info(f"loading {fasta}")
    fa = io.read_fasta(fasta)

    g = df.groupby("protein")

    # breakpoint()

    def process_frame(key_frame, fa, modisite):
        key, frame = key_frame
        # if 'Cont' not in key:
        #     print()
        #     import ipdb; ipdb.set_trace()
        subfa = fa[(fa["id"] == key) & (~fa["id"].str.startswith(DECOY_FLAG))]
        if len(subfa) == 0:
            logger.info(f"skipping key {key}, db mismatch")
            return
        seqinfo = subfa.iloc[0].to_dict()
        res = modisite.main(frame, seqinfo)
        return res

    fullres = list()

    if cores == 1:
        for item in tqdm(g):
            res = process_frame(item, fa, modisite)
            fullres.append(res)

    if cores > 1:
        with ThreadPoolExecutor() as executor:
            futures = {
                executor.submit(process_frame, item, fa, modisite): item for item in g
            }

            for future in tqdm(as_completed(futures), total=len(futures)):
                try:
                    result = future.result()
                    fullres.append(result)
                except Exception as e:

                    # print(f"An error occurred: {e}")
                    pass

    res1 = [x["sty_79_9663"] for x in fullres if x is not None]
    if len(res1) == 0:
        logger.error("no results returned")
        sys.exit(1)
    fullres_df = pd.concat(res1, ignore_index=True)

    infile = pathlib.Path(psms[0])
    outfile = infile.parent / f"site_annotation_notcondensed.tsv"

    fullres_df.to_csv(outfile, sep="\t", index=False)

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
