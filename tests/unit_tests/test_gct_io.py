import pandas as pd
from pathlib import Path

from site_annotate.io.io import write_gct, read_gct


def test_write_read_gct(tmp_path: Path):
    emat = pd.DataFrame({
        'S1': [1.0, 2.0],
        'S2': [3.0, 4.0],
    }, index=['rowA', 'rowB'])
    cdesc = pd.DataFrame({
        'group': ['ctrl', 'treat']
    }, index=['S1', 'S2'])
    rdesc = pd.DataFrame({
        'symbol': ['A', 'B']
    }, index=['rowA', 'rowB'])

    out = tmp_path / 'mini'
    write_gct(emat, cdesc=cdesc, rdesc=rdesc, filename=str(out))

    # discover written gct path
    gct_files = list(tmp_path.glob('mini_*x*.gct'))
    assert gct_files, 'No GCT written'

    emat2, cdesc2, rdesc2 = read_gct(str(gct_files[0]))
    assert list(emat2.columns) == ['S1', 'S2']
    assert list(emat2.index) == ['rowA', 'rowB']
    assert float(emat2.loc['rowA', 'S1']) == 1.0
    assert list(cdesc2.index) == ['S1', 'S2']
    assert list(cdesc2.columns) == ['group']
    assert list(rdesc2.index) == ['rowA', 'rowB']
    assert list(rdesc2.columns) == ['symbol']

