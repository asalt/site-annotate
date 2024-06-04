# tests/test_io_ext.py
import os
import pytest
from io import StringIO
from pathlib import Path
import pytest
from click.testing import CliRunner
import site_annotate.io_external
from pytest_mock import MockerFixture

def test_module_load():
    dir(site_annotate.io_external)

