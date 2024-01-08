import pytest

from ensembl_lite import maf


def test_read(DATA_DIR):
    path = DATA_DIR / "sample.maf"
    blocks = list(maf.parse(path))
    assert len(blocks) == 4


@pytest.mark.parametrize(
    "line",
    (
        "s pan_paniscus.11 2 7 + 13 ACTCTCCAGATGA",
        "s pan_paniscus.11 4 7 - 13 ACTCTCCAGATGA",
    ),
)
def test_process_maf_line_plus(line):
    n, s = maf.process_maf_line(line)
    assert s == "ACTCTCCAGATGA"
    # maf is zero based
    assert n.start == 2
    assert n.end == 2 + 7
