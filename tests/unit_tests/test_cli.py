# tests/test_cli.py
import pytest
from click.testing import CliRunner
from site_annotate.cli import main, run


def test_greet():
    runner = CliRunner()
    result = runner.invoke(
        run,
        [
            "--cores",
            "2",
            # "--psms",
            # "file1.psm",
            # "--psms",
            # "file2.psm",
            # "--out",
            # "output.file",
            # "--fasta",
            # "file.fasta",
        ],
    )
    assert result.exit_code == 0


def test_missing_psms():
    runner = CliRunner()
    result = runner.invoke(
        run, ["--cores", "2", "--output-dir", ".", "--fasta", "file.fasta"]
    )
    assert result.exit_code != 0
    assert "Error" in result.output


# Add more tests as needed
