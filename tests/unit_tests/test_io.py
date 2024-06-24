# tests/test_io.py
import os
import pandas as pd
from io import StringIO
from pathlib import Path
import pytest
from click.testing import CliRunner
import site_annotate.io
from pytest_mock import MockerFixture


def test_read_psm_file():
    psm_file = Path(__file__).parent.parent / "testdata/psm.tsv"
    print(os.path.abspath(psm_file))
    df = site_annotate.io.read_psm_file(psm_file)

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
    df = site_annotate.io.get_isoform_hierarchy()
    assert df is not None
    assert "gene_id" in df.columns
    assert "protein_id" in df.columns
    assert "primary_select" in df.columns
    assert "secondary_select" in df.columns


def test_read_sequence_file_to_dataframe(mock_fasta_file, mocker):
    mocker.patch("builtins.open", return_value=mock_fasta_file)
    df = site_annotate.io.read_sequence_file_to_dataframe("dummy_path.fasta", "fasta")

    assert not df.empty
    assert len(df) == 2
    assert list(df.columns) == ["id", "description", "sequence"]
    assert df.iloc[0]["id"] == "ENSP|ENSP00000483084"
    assert df.iloc[0]["sequence"] == "ATGCGTA"


def test_extract_info_from_header():
    header = "ENSP|ENSMUSP00000022222|ENST|ENSMUST00000022222|ENSG|ENSMUSG00000021709|geneid|59079|taxon|10090|symbol|Erbin|Erbin"
    result = site_annotate.io.extract_info_from_header(header)

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

    df = site_annotate.io.read_fasta("dummy_path.fasta")

    assert not df.empty
    # assert len(df) == 2
    # assert df.iloc[0]["ENSP"] == "ENSP00000483084"
    # assert df.iloc[0]["sequence"] == "ATGCGTA"
    assert "ENSP" in df.columns


def test_read_fasta_testdata():
    fasta_file = Path(__file__).parent.parent / "testdata/test.fa"
    if not fasta_file.exists():
        pytest.skip("Test data not found")
    df = site_annotate.io.read_fasta(fasta_file)
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

    df_new = site_annotate.io.prepare_psm_file(df)
    assert all(df_new.protein == df_new.mapped_proteins)

    query = df[["spectrum", "peptide", "modified_peptide", "protein"]].value_counts()

    assert all(query == 1)


def test_explode_mapped_proteins_realdata():
    f = "../testdata/test_isoforms_tcf12.tsv"  # the names have already been cleaned with janitor
    df = pd.read_csv(f, sep="\t")
    df = site_annotate.io.prepare_psm_file(df)

    query = df[["spectrum", "peptide", "modified_peptide", "protein"]].value_counts()
    assert all(query == 1)
    # assert all(
    #     df[["spectrum", "peptide", "modified_peptide", "protein"]].value_counts() == 1
    # )
