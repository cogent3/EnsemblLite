from ensembl_cli import maf


def test_read(DATA_DIR):
    path = DATA_DIR / "sample.maf"
    blocks = list(maf.parse(path))
    assert len(blocks) == 4
