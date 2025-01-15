import numpy
import pytest

from ensembl_tui import _align as eti_align
from ensembl_tui import _maf as eti_maf
from ensembl_tui import _storage_mixin as eti_mixin


@pytest.fixture
def db_align(DATA_DIR, tmp_path):
    # apps are callable
    records = eti_maf.load_align_records()(  # pylint: disable=not-callable
        DATA_DIR / "tiny.maf",
    )
    outpath = tmp_path / "blah.sqlitedb"
    db = eti_align.AlignDb(source=outpath)
    db.add_records(records)
    db.close()
    return eti_align.AlignDb(source=outpath)


def test_db_align(db_align):
    orig = len(db_align)
    source = db_align.source
    db_align.close()
    got = eti_align.AlignDb(source=source)
    assert len(got) == orig


def test_db_align_add_records(db_align):
    gap_spans = numpy.array([[2, 5], [7, 1]], dtype=int)
    orig = {
        "source": "blah",
        "block_id": 42,
        "species": "human",
        "seqid": "1",
        "start": 22,
        "stop": 42,
        "strand": "-",
    }
    records = [eti_align.AlignRecord(gap_spans=gap_spans, **orig.copy())]
    db_align.add_records(records)
    got = next(
        iter(
            db_align._execute_sql(
                f"SELECT * from {db_align.table_name} where block_id=?",
                (42,),
            ),
        ),
    )
    got = {k: got[k] for k in orig if k != "id"}
    # gap_spans not stored in the db, but in a GapStore
    assert got == orig


@pytest.mark.parametrize("func", [str, repr])
def test_db_align_repr(db_align, func):
    func(db_align)


@pytest.mark.parametrize(
    "data",
    [numpy.array([], dtype=numpy.int32), numpy.array([0, 3], dtype=numpy.uint8)],
)
def test_array_blob_roundtrip(data):
    blob = eti_mixin.array_to_blob(data)
    assert isinstance(blob, bytes)
    inflated = eti_mixin.blob_to_array(blob)
    assert numpy.array_equal(inflated, data)
    assert inflated.dtype is data.dtype


@pytest.mark.parametrize(
    "data",
    [
        numpy.array([0, 3], dtype=numpy.uint8),
        eti_mixin.array_to_blob(numpy.array([0, 3], dtype=numpy.uint8)),
    ],
)
def test_blob_array(data):
    # handles array or bytes as input
    inflated = eti_mixin.blob_to_array(data)
    assert numpy.array_equal(inflated, numpy.array([0, 3], dtype=numpy.uint8))
