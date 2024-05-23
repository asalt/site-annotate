from site_annotate import modisite


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
