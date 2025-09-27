from pathlib import Path
import pandas as pd
from site_annotate import modisite
from site_annotate.modisite import reset_inner_index
from pytest import approx

import pyfaidx
from pyfaidx import Fasta

import site_annotate.io
import site_annotate.runner as runner


def test_isoform_case_two():
    """
    test case with two isoforms that map to the same gene but are non-overlapping
    """

    fa = Path(__file__).absolute().parent.parent / "testdata/test_isoforms_tcf12.fa"
    f = (
        Path(__file__).absolute().parent.parent / "testdata/test_isoforms_tcf12.tsv"
    )  # the names have already been cleaned with janitor
    df = pd.read_csv(f, sep="\t")
    df = site_annotate.io.prepare_psm_file(df)
    all(df[["spectrum", "peptide", "modified_peptide", "protein"]].value_counts() == 1)

    fa_fx = Fasta(fa)

    res = runner.run_pipeline(df, fa_fx, cores=1)
