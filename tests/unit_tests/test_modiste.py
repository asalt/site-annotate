import pandas as pd
from site_annotate import modisite
from pytest import approx

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
    assert extracted_info == {
            4: { "AA" : "d", "prob" : 1 }
            }

def test_extract_positions_misc():
    sequence = "abcd(2.22.22)e"
    extracted_info = modisite.extract_positions(sequence)
    assert len(extracted_info) == 0


def test_create_15mer():
    sequence = "SGGRRKSASATSSSS"
    position = 1
    expected_result = "______sGGRRKSAS"
    result = modisite.create_15mer(sequence, position)
    assert result == expected_result

def test_create_15mer_position_is_float():
    sequence = "SGGRRKSASATSSSS"
    position = 1.0 
    expected_result = "______sGGRRKSAS"
    result = modisite.create_15mer(sequence, position)
    assert result == expected_result


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
