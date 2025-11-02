# mapper.py

import pandas as pd
import pytest
from unittest.mock import MagicMock, patch
from site_annotate.mapper import add_uniprot, find_ENSP


def test_find_ENSP():
    assert find_ENSP("gene|ENSP|XYZ123|more") == "XYZ123"
    assert find_ENSP("no_ensp_here") is None


@patch("site_annotate.mapper.get_db")
def test_add_uniprot(mock_get_db):
    df = pd.DataFrame(
        {"protein": ["gene|ENSP|ENSP00000440235|more", "another|ENSP|ABC789|data"]}
    )
    mock_db = mock_get_db.return_value.__enter__.return_value
    mock_db.__contains__.side_effect = lambda x: x in ["ENSP00000440235"]
    mock_db.__getitem__.side_effect = lambda x: {
        "ENSP00000440235": {"uniprot": {"Swiss-Prot": "Q8N1F7"}}
    }

    result_df = add_uniprot(df)
    assert "uniprot_id" in result_df.columns
    assert result_df.at[0, "uniprot_id"] == "Q8N1F7"
    assert pd.isna(result_df.at[1, "uniprot_id"])


def test_add_uniprot_from_protein_ids():
    df = pd.DataFrame(
        {
            "protein": [
                "gene|ENSP|ENSP00000440235|more",
                "gene|ENSP|ENSP00000440236|more",
            ],
            "protein_ids": ["Q8N1F7;Q9XXXX", pd.NA],
        }
    )

    result_df = add_uniprot(df)
    assert result_df.at[0, "uniprot_id"] == "Q8N1F7"
    assert pd.isna(result_df.at[1, "uniprot_id"])


@patch("site_annotate.mapper.get_db")
@patch("site_annotate.mapper.fetch_uniprot_info_online")
@patch("site_annotate.mapper.fetch_uniprot_by_geneids")
def test_add_uniprot_fallback_geneid(mock_by_gid, mock_online, mock_get_db):
    # Simulate no ENSP mapping in DB or online, but geneid resolves
    df = pd.DataFrame(
        {
            "protein": [
                "ENSP|ENSP00000483455|ENST|ENST00000610452|ENSG|ENSG00000204291|geneid|1306|taxon|9606|symbol|COL15A1|"
            ],
            "geneid": [1306],
            "symbol": ["COL15A1"],
            "taxon": [9606],
        }
    )

    mock_db = mock_get_db.return_value.__enter__.return_value
    mock_db.__contains__.side_effect = lambda x: False
    mock_online.return_value = []  # ENSP path returns nothing
    mock_by_gid.return_value = [
        {"query": 1306, "uniprot": {"Swiss-Prot": "Q8N1F7"}, "taxid": 9606}
    ]

    result_df = add_uniprot(df)
    assert result_df.loc[0, "uniprot_id"] == "Q8N1F7"


@patch("site_annotate.mapper.get_db")
@patch("site_annotate.mapper.fetch_uniprot_info_online")
@patch("site_annotate.mapper.fetch_uniprot_by_geneids")
@patch("site_annotate.mapper.fetch_uniprot_by_symbols")
@patch("site_annotate.mapper.fetch_uniprot_by_ensg")
@patch("site_annotate.mapper.fetch_uniprot_by_enst")
def test_add_uniprot_fallback_symbol(
    mock_by_enst, mock_by_ensg, mock_by_sym, mock_by_gid, mock_online, mock_get_db
):
    # Simulate no ENSP and no geneid mapping but symbol resolves
    df = pd.DataFrame(
        {
            "protein": [
                "ENSP|ENSP00000484086|ENST|ENST00000611833|ENSG|ENSG00000136098|geneid|NA|taxon|9606|symbol|NEK3|"
            ],
            "symbol": ["NEK3"],
            "taxon": [9606],
        }
    )

    mock_db = mock_get_db.return_value.__enter__.return_value
    mock_db.__contains__.side_effect = lambda x: False
    mock_online.return_value = []
    mock_by_gid.return_value = []
    mock_by_sym.return_value = [
        {"query": "NEK3", "uniprot": {"Swiss-Prot": "Q8N1F7"}, "taxid": 9606}
    ]
    mock_by_ensg.return_value = []
    mock_by_enst.return_value = []

    result_df = add_uniprot(df)
    assert result_df.loc[0, "uniprot_id"] == "Q8N1F7"
