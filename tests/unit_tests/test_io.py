# tests/test_io.py
import os
import pandas as pd
from io import StringIO
from pathlib import Path
import pytest
from click.testing import CliRunner
from site_annotate.io import io
from pytest_mock import MockerFixture


def test_read_psm_file():
    psm_file = Path(__file__).parent.parent / "testdata/psm.tsv"
    print(os.path.abspath(psm_file))
    df = io.read_psm_file(psm_file)

    # assert df.shape == (3, 3)
    # assert df.columns.tolist() == ["peptide", "protein", "score"]
    # assert df["score"].sum() == 0.6
    # assert df["peptide"].tolist() == ["PEPTIDE1", "PEPTIDE2", "PEPTIDE3"]
    # assert df["protein"].tolist() == ["PROTEIN1", "PROTEIN2", "PROTEIN3"]


# Sample FASTA content
fasta_content = """>ENSP|ENSP00000483084 description 1
ATGCGTA
>ENSP|ENSP00000492052 description 2
GCTAGCTA
"""


# Mock FASTA file as StringIO
@pytest.fixture
def mock_fasta_file():
    return StringIO(fasta_content)


def test_get_isoform_hierarchy():
    df = io.get_isoform_hierarchy()
    assert df is not None
    assert "gene_id" in df.columns
    assert "protein_id" in df.columns
    assert "primary_select" in df.columns
    assert "secondary_select" in df.columns


def test_read_sequence_file_to_dataframe(mock_fasta_file, mocker):
    mocker.patch("builtins.open", return_value=mock_fasta_file)
    df = io.read_sequence_file_to_dataframe("dummy_path.fasta", "fasta")

    assert not df.empty
    assert len(df) == 2
    assert list(df.columns) == ["id", "description", "sequence"]
    assert df.iloc[0]["id"] == "ENSP|ENSP00000483084"
    assert df.iloc[0]["sequence"] == "ATGCGTA"


def test_extract_info_from_header():
    header = "ENSP|ENSMUSP00000022222|ENST|ENSMUST00000022222|ENSG|ENSMUSG00000021709|geneid|59079|taxon|10090|symbol|Erbin|Erbin"
    result = io.extract_info_from_header(header)

    assert isinstance(result, dict)
    assert result == {
        "ENSP": "ENSMUSP00000022222",
        "ENST": "ENSMUST00000022222",
        "ENSG": "ENSMUSG00000021709",
        "geneid": "59079",
        "taxon": "10090",
        "symbol": "Erbin",
    }


def test_read_fasta(mock_fasta_file, mocker: MockerFixture):
    mocker.patch("builtins.open", return_value=mock_fasta_file)

    df = io.read_fasta("dummy_path.fasta")

    assert not df.empty
    # assert len(df) == 2
    # assert df.iloc[0]["ENSP"] == "ENSP00000483084"
    # assert df.iloc[0]["sequence"] == "ATGCGTA"
    assert "ENSP" in df.columns


def test_read_fasta_testdata():
    fasta_file = Path(__file__).parent.parent / "testdata/test.fa"
    if not fasta_file.exists():
        pytest.skip("Test data not found")
    df = io.read_fasta(fasta_file)
    assert "id" in df.columns
    assert "sequence" in df.columns
    assert "ENSP" in df.columns
    assert "geneid" in df.columns


def test_explode_mapped_proteins_testdata():
    """ """

    df = pd.DataFrame(
        {
            "spectrum": ["s1", "s2", "s3"],
            "peptide": ["p1", "p2", "p3"],
            "modified_peptide": ["p1", "p2", "p3"],
            "protein": ["P1", "P2", "P3"],
            "mapped_proteins": ["P1, P2", "P2", "P3"],
            "intensity": [1, 2, 3],
        }
    )

    df_new = io.prepare_psm_file(df)
    assert all(df_new.protein == df_new.mapped_proteins2)

    query = df[["spectrum", "peptide", "modified_peptide", "protein"]].value_counts()

    assert all(query == 1)


def test_explode_mapped_proteins_realdata():

    f = (
        Path(__file__).absolute().parent.parent / "testdata/test_isoforms_tcf12.tsv"
    )  # the names have already been cleaned with janitor

    df = pd.read_csv(f, sep="\t")
    df = io.prepare_psm_file(df)

    query = df[["spectrum", "peptide", "modified_peptide", "protein"]].value_counts()
    assert all(query == 1)
    # assert all(
    #     df[["spectrum", "peptide", "modified_peptide", "protein"]].value_counts() == 1
    # )


def test_prepare_psm_file_diann_like():

    df = pd.DataFrame(
        {
            "sequence": ["AAAK", "BBBK"],
            "sequence_modi": ["AAAK", "BBBK"],
            "precursor_area": [1000.0, 2000.0],
            "protein_ids": ["P12345;P67890", "P54321"],
            "ensps_all": [
                "ENSP|ENSP000001;ENSP|ENSP000002",
                "ENSP|ENSP000003",
            ],
            "precursor_lib_index": ["lib_a", "lib_b"],
            "modified_sequence": [
                "AAAS(UniMod:21)",
                "M(UniMod:35)BBK",
            ],
            "site_occupancy_probabilities": [
                "A{0.000000}A{0.000000}A{0.000000}S(UniMod:21){1.000000}2",
                "M(UniMod:35){0.750000}BBK",
            ],
        }
    )

    df_prepared = io.prepare_psm_file(df.copy())

    assert set(["peptide", "intensity", "protein", "mapped_proteins"]).issubset(
        df_prepared.columns
    )
    assert df_prepared["peptide"].iloc[0] == "AAAK"
    assert (
        df_prepared.loc[df_prepared["peptide"] == "BBBK", "intensity"].unique()[0]
        == 2000.0
    )
    assert set(df_prepared.loc[df_prepared["peptide"] == "AAAK", "protein"]) == {
        "ENSP|ENSP000001",
        "ENSP|ENSP000002",
    }
    mapped = (
        df_prepared.loc[df_prepared["peptide"] == "AAAK", "mapped_proteins2"]
        .explode()
        .dropna()
        .unique()
    )
    assert set(mapped) == {"ENSP|ENSP000001", "ENSP|ENSP000002"}

    sty_vals = df_prepared.loc[df_prepared["peptide"] == "AAAK", "sty_79_9663"]
    assert sty_vals.notna().any()
    assert sty_vals.dropna().iloc[0] == "A(0.000000)A(0.000000)A(0.000000)S(1.000000)"
    sty_best = df_prepared.loc[df_prepared["peptide"] == "AAAK", "sty_79_9663_best_localization"]
    assert sty_best.dropna().iloc[0] == 1.0

    ox_vals = df_prepared.loc[df_prepared["peptide"] == "BBBK", "m_15_9949"]
    assert ox_vals.notna().any()
    assert ox_vals.dropna().iloc[0] == "M(0.750000)BBK"
    ox_best = df_prepared.loc[df_prepared["peptide"] == "BBBK", "m_15_9949_best_localization"]
    assert ox_best.dropna().iloc[0] == 0.75

    spectrum_vals = df_prepared.loc[df_prepared["peptide"] == "AAAK", "spectrum"].unique()
    assert set(spectrum_vals) == {"lib_a"}
