import pytest

from ensembl_tui import _maf as elt_maf


def test_read(DATA_DIR):
    path = DATA_DIR / "sample.maf"
    blocks = list(elt_maf.parse(path))
    assert len(blocks) == 4
    block_ids = {b for b, *_ in blocks}
    assert block_ids == {20060000040557, 20060000042317, 20060000132559, 20060000102888}


def test_process_id_line():
    got = elt_maf.process_id_line("# id: 20060000042317 \n")
    assert got == 20060000042317


@pytest.mark.parametrize(
    "line",
    (
        "s pan_paniscus.11 2 7 + 13 ACTCTCCAGATGA",
        "s pan_paniscus.11 4 7 - 13 ACTCTCCAGATGA",
    ),
)
def test_process_maf_line_plus(line):
    n, s = elt_maf.process_maf_line(line)
    assert s == "ACTCTCCAGATGA"
    # maf is zero based
    assert n.start == 2
    assert n.stop == 2 + 7
