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


@click.command()
@click.option("--cores", default=1, show_default=True)
@click.option(
    "-p", "--psms", type=click.Path(exists=True, dir_okay=False), multiple=True
)
@click.option("-o", "--out", default=None)
@click.option(
    "-f", "--fasta", type=click.Path(exists=True, dir_okay=False), help="fasta file"
)
def main(cores, psms, out, fasta, **kwargs):

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
    fullres = list()
    # for key, frame in tqdm(g):
    #     subfa = fa[fa["id"] == key]
    #     assert len(subfa) == 1
    #     seqinfo = subfa.iloc[0].to_dict()
    #     res = modisite.main(frame, seqinfo)
    #     fullres.append(res)

    def process_frame(key_frame, fa, modisite):
        key, frame = key_frame
        subfa = fa[fa["id"] == key]
        assert len(subfa) == 1
        seqinfo = subfa.iloc[0].to_dict()
        res = modisite.main(frame, seqinfo)
        return res

    fullres = []
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
