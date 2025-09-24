from pathlib import Path
import pandas as pd

from pytest import approx
from pyfaidx import Fasta


import site_annotate.io

# from site_annotate.modisite import reset_inner_index
import site_annotate.runner as runner
from site_annotate import modisite
import site_annotate.reduce


def test_only_common_columns_no_grouping():
    # what does this function really test?
    df = pd.DataFrame(
        {
            "protein": ["p1", "p2", "p3"],
            "uniprot_id": ["u1", "u2", "u3"],
            "ENSP": ["e1", "e2", "e3"],
            "mapped_proteins": ["m1", "m2", "m3"],
            "fifteenmer": ["seq1", "seq2", "seq3"],
            "sitename": [1, 2, 3],
        }
    )
    reduced = site_annotate.reduce._reduce_sites(df)
    assert not reduced.empty, "Should not be empty with valid grouping columns."
    assert all(
        col in reduced.columns
        for col in ["protein", "uniprot_id", "ENSP"]  # , "mapped_proteins"]
    ), "All common columns should be in the result."


def test_aggregation_with_tmt_intensity():
    """ """

    df = pd.DataFrame(
        {
            "fifteenmer": ["seq1", "seq1", "seq2"],
            "sitename": [1, 1, 2],
            "TMT_1_intensity": [100, 200, 300],
            "TMT_2_intensity": [400, 500, 600],
            "hyperscore": [10, 20, 30],
            "rtscore": [5, 15, 25],
            "delta_mass": [0.1, 0.2, 0.3],
            "highest_prob": [0.9, 0.95, 0.99],
            "spectrum": ["s1", "s2", "s3"],
            "intensity": [1000, 2000, 3000],
        }
    )
    reduced = site_annotate.reduce._reduce_sites(df)
    #

    assert reduced.shape[0] == 2, "There should be two groups."
    assert all(
        reduced["intensity_sum"] == pd.Series([3000, 3000], index=reduced.index)
    ), "TMT intensities should be summed correctly."

    r0 = reduced[reduced.fifteenmer == "seq1"]
    r1 = reduced[reduced.fifteenmer == "seq2"]
    assert len(r0) == len(r1) == 1, "Each fifteenmer should have one row."

    assert r0["intensity_max"].iloc[0] == 2000, "Intensity should be max value."
    assert r1["intensity_max"].iloc[0] == 3000, "Intensity should be max value."

    assert r0["nspectra"][0] == 2, "Nspectra should count unique spectra."
    assert r0["hyperscore_best"][0] == 20, "Hyperscore should have the max value."
    assert r0["rtscore_best"][0] == 15, "Rtscore should have the max value."
    assert r0["delta_mass_best"][0] == 0.1, "Delta mass should have the min value."
    assert r0["highest_prob_best"][0] == 0.95, "Highest prob should have the max value."

    assert r0["TMT_1_intensity"][0] == 300, "TMT1 intensity should be summed."
    assert r0["TMT_2_intensity"][0] == 900, "TMT1 intensity should be summed."

    assert r1["nspectra"].iloc[0] == 1, "Nspectra should count unique spectra."
    assert r1["hyperscore_best"].iloc[0] == 30, "Hyperscore should have the max value."
    assert r1["rtscore_best"].iloc[0] == 25, "Rtscore should have the max value."
    assert r1["delta_mass_best"].iloc[0] == 0.3, "Delta mass should have the min value."
    assert (
        r1["highest_prob_best"].iloc[0] == 0.99
    ), "Highest prob should have the max value."


def test_reduce_1():

    f = (
        Path(__file__).absolute().parent.parent
        / "testdata/test_isoforms_tcf12_aftermodi.tsv"
    )  # the names have already been cleaned with janitor
    df = pd.read_csv(f, sep="\t")

    dfres = site_annotate.reduce._reduce_sites(df)

    query = dfres[dfres.fifteenmer == "AGGQAPSsPSYENSL"]
    assert len(query) == 3
    assert query.spectra.nunique() == 2

    query_distinct = query[query.ENSP == "ENSMUSP00000034755"]
    assert len(query_distinct) == 1

    assert not all(query.intensity_max == query_distinct.intensity_max.iloc[0])
    assert query.intensity_max.max() == query_distinct.intensity_max.iloc[0]
