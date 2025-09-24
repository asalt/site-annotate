import importlib
import pytest

MODULES = [
    "site_annotate",
    "site_annotate.cli",
    "site_annotate.io",
    "site_annotate.io.io",
    "site_annotate.io.io_psm",
    "site_annotate.mapper",
    "site_annotate.mapper_backend",
    "site_annotate.modisite",
    "site_annotate.runner",
    "site_annotate.tasks",
    "site_annotate.utils",
]


@pytest.mark.parametrize("module_name", MODULES)
def test_import_module(module_name):
    try:
        importlib.import_module(module_name)
    except ImportError as exc:  # pragma: no cover - skip optional deps
        if module_name == "site_annotate.tasks" and "rpy2" in str(exc):
            pytest.skip("rpy2 is optional for tasks module")
        raise


def test_io_reexports():
    from site_annotate.io import (
        prepare_psm_file,
        read_psm_file,
        validate_psm_file,
    )
    from site_annotate.io import io as io_module
    from site_annotate.io.io_psm import (
        prepare_psm_file as impl_prepare,
        read_psm_file as impl_read,
        validate_psm_file as impl_validate,
    )

    assert prepare_psm_file is impl_prepare
    assert read_psm_file is impl_read
    assert validate_psm_file is impl_validate
    assert hasattr(io_module, "read_fasta")
