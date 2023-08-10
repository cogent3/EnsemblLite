from ensembl_cli.aligndb import AlignDb
from ensembl_cli.install import _load_one_align


def test_db_align(DATA_DIR, tmp_path):
    records = _load_one_align(DATA_DIR / "tiny.maf")
    outpath = tmp_path / "blah.sqlitedb"
    db = AlignDb(source=outpath)
    db.add_records(records)
    orig = len(db)
    db.close()
    got = AlignDb(source=outpath)
    assert len(got) == orig
