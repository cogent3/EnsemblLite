import shutil

import duckdb
import numpy
import pytest

import ensembl_tui._config as eti_config
import ensembl_tui._ingest_annotation as eti_db_ingest
import ensembl_tui._mysql_core_attr as eti_tables

TRANSLATION_SCHEMA = (
    "transcript_id INTEGER",
    "start_exon_id INTEGER",
    "seq_start INTEGER",
    "end_exon_id INTEGER",
    "seq_end INTEGER",
    "stable_id TEXT",
)

TRANSLATION_COLS = [c.split()[0] for c in TRANSLATION_SCHEMA]

EXON_VIEW_SCHEMA = (
    "transcript_id INTEGER",
    "exon_id INTEGER",
    "seqid TEXT",
    "start INTEGER",
    "stop INTEGER",
    "strand INTEGER",
    "rank INTEGER",
    "gene_id INTEGER",
    "transcript_stable_id TEXT",
    "cds_stable_id TEXT",
)

EXON_VIEW_COLS = [c.split()[0] for c in EXON_VIEW_SCHEMA]

SEQ_REGION_SCHEMA = (
    "seq_region_id INTEGER",
    "name TEXT",
)

SEQ_REGION_COLS = [c.split()[0] for c in SEQ_REGION_SCHEMA]

EXON_SCHEMA = (
    "exon_id INTEGER",
    "transcript_id INTEGER",
    "seq_region_id INTEGER",
    "seq_region_start INTEGER",
    "seq_region_end INTEGER",
    "seq_region_strand INTEGER",
)

EXON_COLS = [c.split()[0] for c in EXON_SCHEMA]

EXON_TRANSCRIPT_SCHEMA = (
    "transcript_id INTEGER",
    "exon_id INTEGER",
    "rank INTEGER",
)

EXON_TRANSCRIPT_COLS = [c.split()[0] for c in EXON_TRANSCRIPT_SCHEMA]

TRANSCRIPT_SCHEMA = (
    "transcript_id INTEGER",
    "gene_id INTEGER",
    "stable_id TEXT",
)

TRANSCRIPT_COLS = [c.split()[0] for c in TRANSCRIPT_SCHEMA]

GENE_SCHEMA = (
    "gene_id INTEGER",
    "stable_id TEXT",
)

GENE_COLS = [c.split()[0] for c in GENE_SCHEMA]


@pytest.fixture
def empty_db():
    return duckdb.connect(":memory:")


@pytest.fixture
def db_with_start_columns(empty_db):
    sql = """CREATE TABLE demo_table (
    name VARCHAR(2),
    seq_region_start INTEGER,
    seq_start INTEGER,
    )"""
    empty_db.sql(sql)
    return empty_db


@pytest.fixture
def db_without_start_columns(empty_db):
    sql = """CREATE TABLE demo_table (
    name VARCHAR(2),
    start_seq_region INTEGER NULL,
    start_seq INTEGER NULL,
    )"""
    empty_db.sql(sql)
    return empty_db


def test_get_start_column(db_with_start_columns):
    starts = eti_db_ingest.get_start_column(db_with_start_columns, "demo_table")
    assert starts == {"seq_region_start", "seq_start"}


def test_get_start_column_none(db_without_start_columns):
    starts = eti_db_ingest.get_start_column(db_without_start_columns, "demo_table")
    assert starts == set()


@pytest.fixture
def tsv_with_start_cols(tmp_path):
    tsv_path = tmp_path / "demo.tsv"
    tsv_path.write_text("ab\t1\t2\ncd\t3\t4\ned\t\\N\t5")
    return tsv_path


def test_import_mysql_table(db_with_start_columns, tsv_with_start_cols):
    eti_db_ingest.import_mysql_table(
        con=db_with_start_columns,
        mysql_dump_path=tsv_with_start_cols,
        table_name="demo_table",
    )
    rows = db_with_start_columns.sql(
        "SELECT seq_region_start,seq_start FROM demo_table",
    ).fetchall()
    assert rows == [(0, 1), (2, 3), (None, 4)]


@pytest.fixture
def downloaded_cfg(tmp_downloaded):
    return tmp_downloaded / "downloaded.cfg"


def test_write_parquet(downloaded_cfg):
    cfg = eti_config.read_config(downloaded_cfg)
    template_path = cfg.staging_template_path
    genome = "saccharomyces_cerevisiae"
    dump_path = cfg.staging_genomes / genome / "mysql" / "gene.txt.gz"
    assert dump_path.exists()
    dest_dir = downloaded_cfg.parent / "installed" / genome
    pqt_path = eti_db_ingest.write_parquet(
        db_templates=template_path,
        dump_path=dump_path,
        table_name="gene",
        dest_dir=dest_dir,
    )
    assert pqt_path.exists()
    assert pqt_path.suffix == ".parquet"
    db = duckdb.connect(":memory:")
    db.sql(f"CREATE TABLE gene AS SELECT * FROM read_parquet('{pqt_path}')")
    (got,) = db.sql(
        "SELECT COUNT(*) FROM gene WHERE biotype = 'protein_coding'",
    ).fetchone()
    assert got > 6_000


@pytest.mark.slow
def test_install_features(yeast_db):
    # this is a check on expectations rather than execution
    source = yeast_db.source
    pqts = list(source.glob("*.parquet"))
    # tables are seq_region, repeat, repeat_consensus, gene_attr, transcript_attr
    assert len(pqts) == 5


# fail to import if directories or files are missing
@pytest.fixture(params=[1, 2, 3])
def invalid_downloaded_cfg(downloaded_cfg, request):
    cfg = eti_config.read_config(downloaded_cfg)
    genome = "saccharomyces_cerevisiae"
    genome_root = cfg.staging_genomes / genome
    if request.param == 1:
        (genome_root / "mysql" / "gene.txt.gz").unlink()
    elif request.param == 2:
        shutil.rmtree(cfg.staging_template_path)
    else:
        shutil.rmtree(genome_root / "mysql")

    return downloaded_cfg


def test_install_features_invalid(invalid_downloaded_cfg):
    cfg = eti_config.read_config(invalid_downloaded_cfg)
    with pytest.raises(FileNotFoundError):  # noqa: PT012
        app = eti_db_ingest.mysql_dump_to_parquet(config=cfg)
        # an exception during call in finding mysql file
        app.main("saccharomyces_cerevisiae")


@pytest.mark.parametrize(
    "data",
    [
        numpy.array([1, 2, 3], dtype=numpy.int32),
        numpy.array([[1, 2, 3], [4, 5, 6]], dtype=numpy.int32),
        numpy.array([[1, 2, 3], [4, 5, 6]], dtype=numpy.uint8),
    ],
)
def test_array_blob_interconvert(data):
    blob = eti_tables.array_to_blob(data)
    assert isinstance(blob, bytes)
    back = eti_tables.blob_to_array(blob)
    assert numpy.array_equal(data, back)
    assert data.dtype == back.dtype


@pytest.fixture
def empty_db_ev_transcript():
    conn = duckdb.connect()
    conn.sql(
        f"""
        CREATE TABLE exon_view ({', '.join(EXON_VIEW_SCHEMA)});
        CREATE TABLE translation ({', '.join(TRANSLATION_SCHEMA)});
        """,
    )
    return conn


@pytest.fixture
def one_exon(empty_db_ev_transcript):
    conn = empty_db_ev_transcript
    conn.sql(
        """
        INSERT INTO exon_view (exon_id, transcript_id, rank, strand) VALUES (1, 1664, 1, 1);
        INSERT INTO translation (transcript_id, start_exon_id, end_exon_id, seq_start, seq_end) VALUES (1664, 1, 1, 1, 1);
        """,
    )
    return conn


def test_get_limiting_exons_one_exon(one_exon):
    all_exons = eti_tables.get_all_limit_exons(one_exon)
    lex = eti_tables.get_limit_exons(all_exons[1664])
    assert lex.start_rank == lex.stop_rank == 1
    assert lex.rel_start == 0
    assert lex.rel_stop == 1


def test_get_limiting_exons_no_transcript():
    with pytest.raises(ValueError):  # noqa: PT011
        eti_tables.get_limit_exons([])


@pytest.fixture
def two_exon(empty_db_ev_transcript):
    conn = empty_db_ev_transcript
    conn.sql(
        """
        INSERT INTO exon_view (exon_id, transcript_id, rank, strand) VALUES (858807, 269944, 4, 1);
        INSERT INTO exon_view (exon_id, transcript_id, rank, strand) VALUES (858815, 269944, 12, 1);
        INSERT INTO translation (transcript_id, start_exon_id, end_exon_id, seq_start, seq_end) VALUES (269944, 858807, 858815, 29, 54);
        """,
    )
    return conn


def test_get_limiting_exons_two_exons(two_exon):
    all_exons = eti_tables.get_all_limit_exons(two_exon)
    lex = eti_tables.get_limit_exons(all_exons[269944])
    assert lex.start_rank == 4
    assert lex.stop_rank == 12
    assert lex.rel_start == 28
    assert lex.rel_stop == 54


# these IDs are from human release 113
# start_exon != end_exon AND start exon.rank == 1:
#  transcript_id = 269936
# start_exon != end_exon AND start exon.rank > 1:
#  transcript_id = 269944
# start_exon != end_exon AND start exon.rank > 1 and strand == -1:
#  transcript_id = 271626


@pytest.fixture
def empty_ev_tr():
    conn = duckdb.connect()
    exon_view_sql = f"CREATE TABLE exon_view ({', '.join(EXON_VIEW_SCHEMA)})"
    conn.execute(exon_view_sql)

    translation_sql = f"CREATE TABLE translation ({', '.join(TRANSLATION_SCHEMA)})"
    conn.execute(translation_sql)
    return conn


@pytest.fixture
def four_exons(empty_ev_tr):
    conn = empty_ev_tr
    exon_view_data = [
        {
            "transcript_id": 11,
            "exon_id": 1,
            "seqid": "2",
            "start": 100,
            "stop": 200,
            "strand": 1,
            "rank": 1,
            "gene_id": 42,
            "transcript_stable_id": "a1",
            "cds_stable_id": "a1",
        },
        {
            "transcript_id": 11,
            "exon_id": 2,
            "seqid": "2",
            "start": 300,
            "stop": 400,
            "strand": 1,
            "rank": 2,
            "gene_id": 42,
            "transcript_stable_id": "a1",
            "cds_stable_id": "a1",
        },
        {
            "transcript_id": 11,
            "exon_id": 3,
            "seqid": "2",
            "start": 500,
            "stop": 600,
            "strand": 1,
            "rank": 3,
            "gene_id": 42,
            "transcript_stable_id": "a1",
            "cds_stable_id": "a1",
        },
        {
            "transcript_id": 11,
            "exon_id": 4,
            "seqid": "2",
            "start": 700,
            "stop": 800,
            "strand": 1,
            "rank": 4,
            "gene_id": 42,
            "transcript_stable_id": "a1",
            "cds_stable_id": "a1",
        },
    ]
    value_placeholder = "?, " * len(EXON_VIEW_COLS)
    conn.executemany(
        f"INSERT INTO exon_view VALUES ({value_placeholder})",
        parameters=[[r[c] for c in EXON_VIEW_COLS] for r in exon_view_data],
    )
    return conn


@pytest.fixture
def same_tr_cds(four_exons):
    conn = four_exons
    # rel start and end are interpreted as no change, so transcript and cds spans same

    translation_data = [
        {
            "transcript_id": 11,
            "start_exon_id": 1,
            "seq_start": 1,
            "end_exon_id": 4,
            "seq_end": 100,
            "stable_id": "a1",
        },
    ]
    conn.executemany(
        "INSERT INTO translation VALUES (?, ?, ?, ?, ?, ?)",
        parameters=[[r[c] for c in TRANSLATION_COLS] for r in translation_data],
    )
    return conn


def test_get_transcript_record(same_tr_cds):
    tr = next(iter(eti_tables.get_transcript_attr_records(same_tr_cds)))
    assert (tr.transcript_spans == tr.cds_spans).all()


@pytest.fixture
def diff_tr_cds(four_exons):
    # limiting exons are different, but spans are the same
    conn = four_exons
    # rel start and end are interpreted as no change, so shared transcript
    # and cds spans same
    translation_data = [
        {
            "transcript_id": 11,
            "start_exon_id": 2,
            "seq_start": 1,
            "end_exon_id": 3,
            "seq_end": 100,
        },
    ]
    value_placeholder = "?, " * len(TRANSLATION_COLS)
    conn.executemany(
        f"INSERT INTO translation VALUES ({value_placeholder})",
        parameters=[[r.get(c) for c in TRANSLATION_COLS] for r in translation_data],
    )
    return conn


def test_get_transcript_record_limit_exons(diff_tr_cds):
    tr = next(iter(eti_tables.get_transcript_attr_records(diff_tr_cds)))
    assert (tr.transcript_spans[1:-1] == tr.cds_spans).all()


@pytest.fixture
def tr_cds_rel_pos(four_exons):
    # limiting exons are different, but spans are the same
    conn = four_exons
    # rel start and end are non-zero, so start exon begins +1,
    # last exon ends -1
    translation_data = [
        {
            "transcript_id": 11,
            "start_exon_id": 1,
            "seq_start": 2,
            "end_exon_id": 4,
            "seq_end": 2,
        },
    ]
    value_placeholder = "?, " * len(TRANSLATION_COLS)
    conn.executemany(
        f"INSERT INTO translation VALUES ({value_placeholder})",
        parameters=[[r.get(c) for c in TRANSLATION_COLS] for r in translation_data],
    )
    return conn


def test_rel_start_ends_1(tr_cds_rel_pos):
    tr = next(iter(eti_tables.get_transcript_attr_records(tr_cds_rel_pos)))
    assert tr.start == 100
    assert tr.stop == 800
    assert (tr.transcript_spans[1:-1] == tr.cds_spans[1:-1]).all()
    assert not (tr.transcript_spans == tr.cds_spans).all()
    assert (tr.cds_spans[0] == (101, 200)).all()
    assert (tr.cds_spans[-1] == (700, 702)).all()


@pytest.fixture
def two_exons_minus_strand(empty_ev_tr):
    conn = empty_ev_tr
    exon_view_data = [
        {
            "transcript_id": 11,
            "exon_id": 1,
            "seqid": "3",
            "start": 100,
            "stop": 200,
            "strand": -1,
            "rank": 2,
            "gene_id": 42,
            "transcript_stable_id": "a1",
            "cds_stable_id": "a1",
        },
        {
            "transcript_id": 11,
            "exon_id": 2,
            "seqid": "3",
            "start": 300,
            "stop": 400,
            "strand": -1,
            "rank": 1,
            "gene_id": 42,
            "transcript_stable_id": "a1",
            "cds_stable_id": "a1",
        },
    ]
    conn.executemany(
        "INSERT INTO exon_view VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
        [[r[c] for c in EXON_VIEW_COLS] for r in exon_view_data],
    )
    return conn


@pytest.fixture
def tr_cds_rel_pos_minus(two_exons_minus_strand):
    # limiting exons are different, but spans are the same
    conn = two_exons_minus_strand
    # rel start and end are non-zero, so start exon begins +1,
    # last exon ends -1
    translation_data = [
        {
            "transcript_id": 11,
            "start_exon_id": 2,
            "seq_start": 5,
            "end_exon_id": 1,
            "seq_end": 10,
        },
    ]
    value_placeholder = "?, " * len(TRANSLATION_COLS)
    conn.executemany(
        f"INSERT INTO translation VALUES ({value_placeholder})",
        parameters=[[r.get(c) for c in TRANSLATION_COLS] for r in translation_data],
    )
    return conn


def test_rel_start_ends_2(tr_cds_rel_pos_minus):
    tr = next(iter(eti_tables.get_transcript_attr_records(tr_cds_rel_pos_minus)))
    assert tr.start == 100
    assert tr.stop == 400
    assert numpy.array_equal(tr.transcript_spans, [(100, 200), (300, 400)])
    assert numpy.array_equal(tr.cds_spans, [(200 - 10, 200), (300, 400 - 4)])


def test_no_cds_spans(four_exons):
    tr = next(iter(eti_tables.get_transcript_attr_records(four_exons)))
    assert tr.cds_spans is None
    record = tr.to_record(eti_tables.TRANSCRIPT_ATTR_COLS)
    assert record[-3] is None
    assert len(record) == 10


def single_tables_db():
    conn = duckdb.connect()
    sql = f"CREATE TABLE seq_region ({', '.join(SEQ_REGION_SCHEMA)})"
    conn.execute(sql)
    sql = f"CREATE TABLE exon ({', '.join(EXON_SCHEMA)})"
    conn.execute(sql)
    sql = f"CREATE TABLE exon_transcript ({', '.join(EXON_TRANSCRIPT_SCHEMA)})"
    conn.execute(sql)
    sql = f"CREATE TABLE transcript ({', '.join(TRANSCRIPT_SCHEMA)})"
    conn.execute(sql)
    sql = f"CREATE TABLE translation ({', '.join(TRANSLATION_SCHEMA)})"
    conn.execute(sql)
    return conn


@pytest.fixture
def mixed_data():
    conn = single_tables_db()
    # seq_region: seq_region_id, name
    sr_data = [{"seq_region_id": i, "name": str(i)} for i in range(1, 3)]
    conn.executemany(
        "INSERT INTO seq_region VALUES (?, ?)",
        parameters=[[r[c] for c in SEQ_REGION_COLS] for r in sr_data],
    )
    exon_data = [
        {
            "exon_id": 1,
            "transcript_id": 11,
            "seq_region_id": 1,
            "seq_region_start": 100,
            "seq_region_end": 200,
            "seq_region_strand": 1,
        },
        {
            "exon_id": 2,
            "transcript_id": 11,
            "seq_region_id": 1,
            "seq_region_start": 300,
            "seq_region_end": 400,
            "seq_region_strand": 1,
        },
        {
            "exon_id": 3,
            "transcript_id": 11,
            "seq_region_id": 1,
            "seq_region_start": 500,
            "seq_region_end": 600,
            "seq_region_strand": 1,
        },
        {
            "exon_id": 4,
            "transcript_id": 11,
            "seq_region_id": 1,
            "seq_region_start": 700,
            "seq_region_end": 800,
            "seq_region_strand": 1,
        },
        {
            "exon_id": 5,
            "transcript_id": 12,
            "seq_region_id": 2,
            "seq_region_start": 900,
            "seq_region_end": 1000,
            "seq_region_strand": -1,
        },
        {
            "exon_id": 6,
            "transcript_id": 12,
            "seq_region_id": 2,
            "seq_region_start": 1100,
            "seq_region_end": 1200,
            "seq_region_strand": -1,
        },
        {
            "exon_id": 7,
            "transcript_id": 12,
            "seq_region_id": 2,
            "seq_region_start": 1300,
            "seq_region_end": 1400,
            "seq_region_strand": -1,
        },
    ]
    conn.executemany(
        "INSERT INTO exon VALUES (?, ?, ?, ?, ?, ?)",
        parameters=[[r[c] for c in EXON_COLS] for r in exon_data],
    )
    et_data = [
        {"transcript_id": 11, "exon_id": 1, "rank": 1},
        {"transcript_id": 11, "exon_id": 2, "rank": 2},
        {"transcript_id": 11, "exon_id": 3, "rank": 3},
        {"transcript_id": 11, "exon_id": 4, "rank": 4},
        {"transcript_id": 12, "exon_id": 5, "rank": 3},
        {"transcript_id": 12, "exon_id": 6, "rank": 2},
        {"transcript_id": 12, "exon_id": 7, "rank": 1},
    ]
    conn.executemany(
        "INSERT INTO exon_transcript VALUES (?, ?, ?)",
        parameters=[[r[c] for c in EXON_TRANSCRIPT_COLS] for r in et_data],
    )
    tr_data = [
        {"transcript_id": 11, "gene_id": 42, "stable_id": "tr-01"},
        {"transcript_id": 12, "gene_id": 43, "stable_id": "tr-02"},
    ]
    conn.executemany(
        "INSERT INTO transcript VALUES (?, ?, ?)",
        parameters=[[r[c] for c in TRANSCRIPT_COLS] for r in tr_data],
    )

    tl_data = [
        {
            "transcript_id": 12,
            "start_exon_id": 7,
            "seq_start": 1,
            "end_exon_id": 5,
            "seq_end": 100,
            "stable_id": "pr-01",
        },
    ]

    conn.executemany(
        "INSERT INTO translation VALUES (?, ?, ?, ?, ?, ?)",
        parameters=[[r[c] for c in TRANSLATION_COLS] for r in tl_data],
    )
    return conn


def test_make_transcript_attr(mixed_data):
    conn = eti_tables.make_transcript_attr(mixed_data)
    sql = "SELECT * FROM transcript_attr"
    records = conn.sql(sql).fetchall()
    assert len(records) == 2
    (result,) = conn.sql(
        "SELECT cds_spans FROM transcript_attr WHERE transcript_id = 12",
    ).fetchone()
    spans = eti_tables.blob_to_array(result)
    expect = numpy.array([[900, 1000], [1100, 1200], [1300, 1400]], dtype=numpy.int32)
    assert numpy.allclose(spans, expect)
