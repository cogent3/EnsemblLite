import pickle  # nosec B403

import numpy
import pytest
from cogent3 import load_table

from ensembl_tui import _align as elt_align
from ensembl_tui import _homology as elt_homology
from ensembl_tui import _maf as elt_maf
from ensembl_tui import _storage_mixin as elt_mixin


@pytest.fixture(scope="function")
def db_align(DATA_DIR, tmp_path):
    # apps are callable
    records = elt_maf.load_align_records()(  # pylint: disable=not-callable
        DATA_DIR / "tiny.maf",
    )
    outpath = tmp_path / "blah.sqlitedb"
    db = elt_align.AlignDb(source=outpath)
    db.add_records(records)
    db.close()
    return elt_align.AlignDb(source=outpath)


def test_db_align(db_align):
    orig = len(db_align)
    source = db_align.source
    db_align.close()
    got = elt_align.AlignDb(source=source)
    assert len(got) == orig


def test_db_align_add_records(db_align):
    gap_spans = numpy.array([[2, 5], [7, 1]], dtype=int)
    orig = dict(
        source="blah",
        block_id=42,
        species="human",
        seqid="1",
        start=22,
        stop=42,
        strand="-",
    )
    records = [elt_align.AlignRecord(gap_spans=gap_spans, **orig.copy())]
    db_align.add_records(records)
    got = list(
        db_align._execute_sql(
            f"SELECT * from {db_align.table_name} where block_id=?",
            (42,),
        ),
    )[0]
    got = {k: got[k] for k in orig if k != "id"}
    # gap_spans not stored in the db, but in a GapStore
    assert got == orig


@pytest.mark.parametrize("func", (str, repr))
def test_db_align_repr(db_align, func):
    func(db_align)


@pytest.fixture
def hom_dir(DATA_DIR, tmp_path):
    path = DATA_DIR / "small_protein_homologies.tsv.gz"
    table = load_table(path)
    outpath = tmp_path / "small_1.tsv.gz"
    table[:1].write(outpath)
    outpath = tmp_path / "small_2.tsv.gz"
    table[1:2].write(outpath)
    return tmp_path


def test_extract_homology_data(hom_dir):
    loader = elt_homology.load_homologies(
        {"gorilla_gorilla", "nomascus_leucogenys", "notamacropus_eugenii"},
    )
    records = []
    for result in loader.as_completed(hom_dir.glob("*.tsv.gz"), show_progress=False):
        records.extend(result.obj)
    assert len(records) == 2


def test_homology_db(hom_dir):
    loader = elt_homology.load_homologies(
        {"gorilla_gorilla", "nomascus_leucogenys", "notamacropus_eugenii"},
    )

    outpath = hom_dir / "species.sqlitedb"
    db = elt_homology.HomologyDb(source=outpath)
    grouped = [
        r.obj
        for r in loader.as_completed(hom_dir.glob("*.tsv.gz"), show_progress=False)
    ]
    got = elt_homology.merge_grouped(*grouped)
    num_genes = 0
    for rel_type, data in got.items():
        db.add_records(records=data, relationship_type=rel_type)
        for record in data:
            num_genes += len(record.gene_ids)
    assert len(db) == num_genes
    db.close()
    got = elt_homology.HomologyDb(source=outpath)
    assert len(got) == num_genes


def test_pickling_db(db_align):
    # should not fail
    pkl = pickle.dumps(db_align)  # nosec B301
    upkl = pickle.loads(pkl)  # nosec B301
    assert db_align.source == upkl.source


@pytest.mark.parametrize(
    "data",
    (numpy.array([], dtype=numpy.int32), numpy.array([0, 3], dtype=numpy.uint8)),
)
def test_array_blob_roundtrip(data):
    blob = elt_mixin.array_to_blob(data)
    assert isinstance(blob, bytes)
    inflated = elt_mixin.blob_to_array(blob)
    assert numpy.array_equal(inflated, data)
    assert inflated.dtype is data.dtype


@pytest.mark.parametrize(
    "data",
    (
        numpy.array([0, 3], dtype=numpy.uint8),
        elt_mixin.array_to_blob(numpy.array([0, 3], dtype=numpy.uint8)),
    ),
)
def test_blob_array(data):
    # handles array or bytes as input
    inflated = elt_mixin.blob_to_array(data)
    assert numpy.array_equal(inflated, numpy.array([0, 3], dtype=numpy.uint8))
