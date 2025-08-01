from pathlib import Path

import pandas as pd
from site_annotate import modisite
from site_annotate.io import io
from site_annotate.modisite import reset_inner_index
from pytest import approx

from pyfaidx import Fasta

import site_annotate.runner as runner  # runner calls modisite
from site_annotate import modisite

# def test_explode_mapped_proteins_realdata():
#     from pyfaidx import Fasta
#     f = "../testdata/test_isoforms_tcf12.tsv"  # the names have already been cleaned with janitor
#     df = pd.read_csv(f, sep="\t")
#     df = site_annotate.io.prepare_psm_file(df)


def test_explode_mapped_proteins_realdata():

    f = (
        Path(__file__).absolute().parent.parent / "testdata/test_isoforms_tcf12.tsv"
    )  # the names have already been cleaned with janitor
    fa = Path(__file__).absolute().parent.parent / "testdata/test_isoforms_tcf12.fa"

    fa_fx = Fasta(fa)

    df = pd.read_csv(f, sep="\t")
    df = io.prepare_psm_file(df)

    res = runner.run_pipeline(df, fa_fx, cores=1)  # this uses modisite funcs
    resl = [x["sty_79_9663"] for x in res]
    dfl = pd.concat(resl)

    query = dfl[
        ["spectrum", "peptide", "modified_peptide", "protein", "position_absolut"]
    ].value_counts()

    # "AGGQAPSSPSYENSLHSLQSR"
    assert (query == 1).all()


def test_extract_positions_empty():
    sequence = "XXXXX"
    extracted_info = modisite.extract_positions(sequence)
    assert isinstance(extracted_info, dict)
    assert len(extracted_info) == 0


def test_extract_positions():
    sequence = "DGS(0.0032)GGAS(0.0033)GT(0.0037)LQPS(0.1980)S(0.1980)GGGS(0.1980)S(0.1980)NS(0.1980)RER"
    extracted_info = modisite.extract_positions(sequence)

    expected_result = {
        3: {"AA": "S", "prob": 0.0032},
        7: {"AA": "S", "prob": 0.0033},
        9: {"AA": "T", "prob": 0.0037},
        13: {"AA": "S", "prob": 0.198},
        14: {"AA": "S", "prob": 0.198},
        18: {"AA": "S", "prob": 0.198},
        19: {"AA": "S", "prob": 0.198},
        21: {"AA": "S", "prob": 0.198},
    }

    assert extracted_info == expected_result


def test_extract_positions_singleval():
    sequence = "abcd(1)e"
    extracted_info = modisite.extract_positions(sequence)
    assert extracted_info == {4: {"AA": "d", "prob": 1}}


def test_extract_positions_misc():
    sequence = "abcd(2.22.22)e"
    extracted_info = modisite.extract_positions(sequence)
    assert len(extracted_info) == 0


def test_create_15mer():
    sequence = "SGGRRKSASATSSSS"
    position = 1
    expected_result = "_______sGGRRKSA"
    result = modisite.create_15mer(sequence, position)
    assert result == expected_result


def test_create_15mer_position_is_float():
    sequence = "SGGRRKSASATSSSS"
    position = 1.0
    expected_result = "_______sGGRRKSA"
    result = modisite.create_15mer(sequence, position)
    assert result == expected_result


def test_flow():
    list_data = [
        "S(0.7959)T(0.1760)EDLS(0.0282)PQR",
        "S(0.8889)T(0.0618)EDLS(0.0493)PQR",
        "S(0.0416)ES(0.8863)AENHS(0.0405)Y(0.0316)AK",
        "S(0.0412)ES(0.9210)AENHS(0.0195)Y(0.0183)AK",
        "QS(0.1981)S(0.1981)VS(0.1981)S(0.1981)T(0.1981)AS(0.0057)VNLGDPT(0.0037)R",
        "AY(0.0958)S(0.8444)IDGPNT(0.0206)S(0.0206)RPQS(0.0187)AR",
        "AY(0.0508)S(0.8956)IDGPNT(0.0138)S(0.0138)RPQS(0.0132)ARPS(0.0129)INEIPER",
        "AY(0.0252)S(0.9106)IDGPNT(0.0167)S(0.0167)RPQS(0.0158)ARPS(0.0149)INEIPER",
        "AS(0.0084)S(0.0084)S(0.0085)ADVGIS(0.4675)KS(0.4675)T(0.0312)EDLS(0.0085)PQR",
        "AS(0.0071)S(0.0074)S(0.0074)ADVGIS(0.3237)KS(0.3237)T(0.3237)EDLS(0.0072)PQR",
        "AS(0.0085)S(0.0087)S(0.0087)ADVGIS(0.4691)KS(0.4691)T(0.0272)EDLS(0.0086)PQR",
        "S(0.0203)HS(0.6993)IT(0.2562)NMET(0.0243)GGLK",
        "S(0.0208)HS(0.8978)IT(0.0507)NMET(0.0307)GGLK",
        "S(0.4561)ES(0.4561)AENHS(0.0537)Y(0.0340)AK",
        "S(0.4731)ES(0.4731)AENHS(0.0328)Y(0.0210)AK",
    ]
    data = pd.Series(list_data)
    res = data.apply(modisite.extract_positions)
    assert len(res) == len(list_data) == len(data)

    assert res.apply(lambda x: isinstance(x, dict)).all()

    res_series_of_dfs = res.apply(modisite.position_dict_to_df)
    assert isinstance(res_series_of_dfs, pd.Series)

    assert res_series_of_dfs.apply(lambda x: isinstance(x, pd.DataFrame)).all()

    res_df = modisite.reset_inner_index(res_series_of_dfs)

    assert isinstance(res_df, pd.DataFrame)

    assert "original_index" in res_df.columns


def test_reset_inner_index():
    data = {
        1: pd.DataFrame({"a": [1, 2, 3]}),
        2: pd.DataFrame({"b": [4, 5, 6]}),
    }
    result = modisite.reset_inner_index(data)
    assert "original_index" in result.columns
    assert result["original_index"].nunique() == 2


#
def test_reset_inner_index_basic():
    # Creating a sample series of dataframes
    dfs = {
        19: pd.DataFrame(
            {
                "position_relative": [1, 2, 6],
                "AA": ["S", "T", "S"],
                "prob": [0.7959, 0.176, 0.0282],
            }
        ),
        20: pd.DataFrame(
            {"position_relative": [1, 2], "AA": ["S", "T"], "prob": [0.8889, 0.0618]}
        ),
    }

    # Expected DataFrame
    expected_df = pd.DataFrame(
        {
            "position_relative": [1, 2, 6, 1, 2],
            "AA": ["S", "T", "S", "S", "T"],
            "prob": [0.7959, 0.176, 0.0282, 0.8889, 0.0618],
            "original_index": [19, 19, 19, 20, 20],
        }
    )

    # Test the function
    result_df = reset_inner_index(dfs)
    pd.testing.assert_frame_equal(result_df, expected_df)


def test_reset_inner_index_empty():
    # Test with empty dict
    result_df = reset_inner_index(pd.Series())
    expected_df = pd.DataFrame()
    pd.testing.assert_frame_equal(result_df, expected_df)


def test_reset_inner_index_single():
    # Test with a single DataFrame
    dfs = {
        21: pd.DataFrame(
            {
                "position_relative": [1, 2, 6],
                "AA": ["S", "T", "S"],
                "prob": [0.7959, 0.176, 0.0282],
            }
        )
    }
    expected_df = pd.DataFrame(
        {
            "position_relative": [1, 2, 6],
            "AA": ["S", "T", "S"],
            "prob": [0.7959, 0.176, 0.0282],
            "original_index": [21, 21, 21],
        }
    )

    result_df = reset_inner_index(dfs)
    pd.testing.assert_frame_equal(result_df, expected_df)


def test_quant_isobaric_site():
    data = {
        "TMT_129_N": [18721.447],
        "TMT_129_C": [13258.230],
        "TMT_130_N": [16503.178],
        "TMT_130_C": [12118.164],
        "TMT_131_N": [17498.963],
        "TMT_131_C": [13208.370],
        "TMT_132_N": [11776.589],
        "TMT_133_N": [12446.593],
        "TMT_133_C": [11923.287],
        "TMT_134_N": [10962.315],
        "intensity": [10000],
    }
    psms_positions = pd.DataFrame(data)

    result_df = modisite.quant_isobaric_site(psms_positions)

    # Check that ratios and intensity distribution columns are present
    for col in [
        "TMT_129_N_ratio",
        "TMT_129_C_ratio",
        "TMT_130_N_ratio",
        "TMT_130_C_ratio",
        "TMT_131_N_ratio",
        "TMT_131_C_ratio",
        "TMT_132_N_ratio",
        "TMT_133_N_ratio",
        "TMT_133_C_ratio",
        "TMT_134_N_ratio",
    ]:
        assert col in result_df.columns

    for col in [
        "TMT_129_N_intensity",
        "TMT_129_C_intensity",
        "TMT_130_N_intensity",
        "TMT_130_C_intensity",
        "TMT_131_N_intensity",
        "TMT_131_C_intensity",
        "TMT_132_N_intensity",
        "TMT_133_N_intensity",
        "TMT_133_C_intensity",
        "TMT_134_N_intensity",
    ]:
        assert col in result_df.columns

    # Check that the sum of TMT columns is correct
    expected_sum = sum(data[col][0] for col in data if col.startswith("TMT_"))
    result_df["tmt_sum"].iloc[0] == approx(expected_sum, rel=1e-9)


def xx_test_main_basic():
    from unittest.mock import patch

    df = pd.DataFrame(
        {
            "sty_79_9663": [0.1, None, 0.3],
            "k_42_0106": [None, 0.2, 0.3],
            # other columns as needed...
        }
    )
    seqinfo = {
        "sequence": "MKLVT",
        "psp": {"sequence": "MKLVT"},
        "id": "protein1",
        "description": "example protein",
        "geneid": "gene123",
        "taxon": "9606",
        "symbol": "PR1",
        "ENSP": "ENSP00000369497",
    }

    # Mock the external functions called within `main`
    with patch("modisite.extract_positions") as mock_extract:
        with patch("modisite.position_dict_to_df") as mock_pos2df:
            with patch("modisite.reset_inner_index") as mock_reset:
                with patch("modisite.quant_isobaric_site") as mock_quant:
                    # Setup the mock returns
                    mock_extract.return_value = pd.Series([{"pos": 1}, {"pos": 2}])
                    mock_pos2df.return_value = pd.DataFrame(
                        {"position_relative": [1, 2], "original_index": [0, 0]}
                    )
                    mock_reset.return_value = pd.DataFrame(
                        {
                            "position_relative": [1, 2],
                            "original_index": [0, 0],
                            "prob": [0.7, 0.8],
                        }
                    )
                    mock_quant.return_value = pd.DataFrame(
                        {
                            "position_relative": [1, 2],
                            "original_index": [0, 0],
                            "prob": [0.7, 0.8],
                        }
                    )

                    # Run test
                    results = modisite.main(df, seqinfo)
                    # Add  assertions here to verify correct results
