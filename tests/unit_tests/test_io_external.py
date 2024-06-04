# tests/test_io_ext.py
import os
from pathlib import Path
import pytest
from io import StringIO
from pathlib import Path
import pytest
from click.testing import CliRunner
import site_annotate.io_external
from pytest_mock import MockerFixture


def test_module_load():
    dir(site_annotate.io_external)


def test_get_data_dir():
    res = site_annotate.io_external.set_data_dir()
    res = Path(res)

    assert res.resolve() == (Path(__file__).parent.parent.parent / "data").resolve()


def test_check_files_exists():
    res = site_annotate.io_external.check_files_exist(["Phosphorylation_site_dataset"])
    assert res == True


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
