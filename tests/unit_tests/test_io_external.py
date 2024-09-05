# tests/test_io_ext.py
import os
import tempfile
from pathlib import Path
import pytest
from io import StringIO
from pathlib import Path
import pytest
from click.testing import CliRunner
from pytest_mock import MockerFixture

import site_annotate.io_external
from site_annotate import io_external


def test_module_load():
    dir(site_annotate.io_external)


def test_get_data_dir():
    res = site_annotate.io_external.set_data_dir()
    res = Path(res)

    assert (
        res.resolve()
        == (Path(__file__).parent.parent.parent / "data" / "phosphositeplus").resolve()
    )


# def test_check_files_exists():
#     res = site_annotate.io_external.check_files_exist(["Phosphorylation_site_dataset"])
#     assert res == True


def test_load_file():
    file = "Phosphorylation_site_dataset"

    res = site_annotate.io_external.check_files_exist([file])
    if res == False:
        pytest.skip(f"file : {file} not found, skipping")
    expected_names = [
        "gene",
        "protein",
        "acc_id",
        "hu_chr_loc",
        "mod_rsd",
        "site_grp_id",
        "organism",
        "mw_kd",
        "domain",
        "site_+_7_aa",
        "lt_lit",
        "ms_lit",
        "ms_cst",
        "cst_cat#",
        "ambiguous_site",
    ]

    df = site_annotate.io_external.get_psiteplus_file(file, nrows=10)

    for name in expected_names:
        assert name in df.columns


TEST_FASTA = """>GN:Hist1h4m|H4H4|mouse|P62806
SGRGKGGKGLGKGGAKRHRKVLRDNIQGITKPAIRRLARRGGVKRISGLIYEETRGVLKV
FLENVIRDAVTYTEHAKRKTVTAMDVVYALKRQGRTLYGFGG
>GN:H4c2|H4|mouse|P62806
SGRGKGGKGLGKGGAKRHRKVLRDNIQGITKPAIRRLARRGGVKRISGLIYEETRGVLKV
FLENVIRDAVTYTEHAKRKTVTAMDVVYALKRQGRTLYGFGG
"""
TEST_FASTA_HEAD = """;051724
;Data extracted from PhosphoSitePlus(R), created by Cell Signaling Technology Inc. PhosphoSitePlus is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License. Attribution must be given in written, oral and digital presentations to PhosphoSitePlus, www.phosphosite.org. Written documents should additionally cite Hornbeck PV, Kornhauser JM, Tkachev S, Zhang B, Skrzypek E, Murray B, Latham V, Sullivan M (2012) PhosphoSitePlus: a comprehensive resource for investigating the structure and function of experimentally determined post-translational modifications in man and mouse. Nucleic Acids Res. 40, D261?70.; www.phosphosite.org.

>GN:Cbln1|CBLN1|mouse|Q9R171
MLGVVELLLLGTAWLAGPARGQNETEPIVLEGKCLVVCDSNPTSDPTGTALGISVRSGSA
KVAFSAIRSTNHEPSEMSNRTMIIYFDQVLVNIGNNFDSERSTFIAPRKGIYSFNFHVVK
VYNRQTIQVSLMLNGWPVISAFAGDQDVTREAASNGVLIQMEKGDRAYLKLERGNLMGGW
KYSTFSGFLVFPL
>GN:Cox7a2|COX7A2|mouse|P48771
MLRNLLALRQIAQRTISTTSRRHFENKVPEKQKLFQEDNGMPVHLKGGASDALLYRATMA
"""


def test_read_psite_fasta():
    # Create a temporary file to simulate the fasta file
    with tempfile.NamedTemporaryFile(mode="w+", delete=False) as tmpfile:
        tmpfile.write(TEST_FASTA_HEAD + TEST_FASTA)
        tmpfile.flush()  # Ensure all data is written

        # Run the function with the temporary file
        res = io_external.read_psite_fasta(tmpfile.name)

    # You can now perform asserts
    assert len(res) > 0, "No entries found in the FASTA data."

    # Optionally, check for specific content or properties
    assert "P62806" in res.keys(), "Expected key not found in FASTA entries."
