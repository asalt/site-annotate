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
