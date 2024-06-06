# tests/test_run.py
import traceback
import pytest
from pathlib import Path
from click.testing import CliRunner
from site_annotate.cli import main, run


def get_test_data():

    data_dir = Path(__file__).parent.parent / "testdata"
    assert data_dir.exists()
    psms_test_file = data_dir / "psm.tsv"
    assert psms_test_file.exists()

    fasta_test_file = data_dir / "test.fa"
    assert fasta_test_file.exists()

    return dict(psms=psms_test_file, fasta=fasta_test_file)


def test_get_test_data():
    get_test_data()  # will rase assertionerror if fail


# psms_file


def test_run():

    test_files = get_test_data()
    psms = test_files.get("psms")
    fa = test_files.get("fasta")

    runner = CliRunner()
    result = runner.invoke(run, ["--cores", "2", "--psms", psms, "--fasta", fa,
    "--no-uniprot-check"
    ])

    # assert result.exit_code == 0
    # import pdb; pdb.set_trace()
    try:
        assert result.exit_code == 0
        assert result.exception is None
    except AssertionError as e:

        tb = traceback.format_exception(*result.exc_info)
        print(''.join(tb))
        raise e

