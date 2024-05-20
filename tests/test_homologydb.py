import pytest

from cogent3 import load_table

from ensembl_lite._homologydb import (
    _HOMOLOGYDB_NAME,
    HomologyDb,
    HomologyRecord,
    grouped_related,
    homolog_group,
    merge_grouped,
    species_genes,
)
from ensembl_lite._install import LoadHomologies


def _make_expected_o2o(table):
    """return dict with keys stable ID's values list[tuple[str,str]]"""
    data = table.to_list(["gene_id_1", "gene_id_2"])

    result = {}
    for g1, g2 in data:
        value = {g1, g2}
        result[g1] = result.get(g1, set()) | value
        result[g2] = result.get(g2, set()) | value

    return result


@pytest.fixture
def o2o_db(DATA_DIR, tmp_dir):
    raw = DATA_DIR / "one2one_homologies.tsv"

    species = {
        "gorilla_gorilla",
        "macaca_mulatta",
        "microcebus_murinus",
        "homo_sapiens",
        "pongo_abelii",
        "pan_troglodytes",
        "macaca_fascicularis",
        "chlorocebus_sabaeus",
        "pan_paniscus",
    }
    loader = LoadHomologies(species)

    table = load_table(raw).get_columns(loader.src_cols)

    table = table.with_new_header(loader.src_cols, loader.dest_col)
    table = table.get_columns(["relationship", "gene_id_1", "gene_id_2"])
    data = loader([raw])
    # convert HomologyRecords into homology_group instances
    hom_groups = grouped_related(data)
    homdb = HomologyDb(tmp_dir / _HOMOLOGYDB_NAME)
    for rel_type, data in hom_groups.items():
        homdb.add_records(records=data, relationship_type=rel_type)
    return homdb, table


@pytest.mark.parametrize(
    "gene_id",
    (
        "ENSGGOG00000026757",
        "ENSGGOG00000025053",
        "ENSGGOG00000022688",
        "ENSGGOG00000026221",
        "ENSGGOG00000024015",
    ),
)
def test_hdb(o2o_db, gene_id):
    homdb, table = o2o_db
    expect = _make_expected_o2o(table)

    got = homdb.get_related_to(gene_id=gene_id, relationship_type="ortholog_one2one")
    assert set(got) == expect[gene_id]


@pytest.fixture
def orth_records():
    common = dict(relationship="ortholog_one2one")
    return [
        HomologyRecord(gene_id_1="1", gene_id_2="2", **common),  # grp 1
        HomologyRecord(gene_id_1="2", gene_id_2="3", **common),  # grp 1
        HomologyRecord(gene_id_1="4", gene_id_2="5", **common),
    ]


@pytest.fixture
def hom_records(orth_records):
    return orth_records + [
        HomologyRecord(
            gene_id_1="6",
            gene_id_2="7",
            relationship="ortholog_one2many",
        ),
    ]


def test_set_hom_set_value():
    r = HomologyRecord(
        gene_id_1="6",
        gene_id_2="7",
    )
    assert r.relationship is None
    r["relationship"] = "ortholog_one2many"
    assert r.relationship == "ortholog_one2many"


def test_hdb_get_related_groups(o2o_db):
    homdb, _ = o2o_db
    got = homdb.get_related_groups(relationship_type="ortholog_one2one")
    assert len(got) == 5


@pytest.fixture
def hom_hdb(hom_records):
    col_order = (
        "relationship",
        "gene_id_1",
        "gene_id_2",
    )
    records = [
        HomologyRecord(**dict(zip(col_order, [r[c] for c in col_order])))
        for r in hom_records
    ]
    groups = grouped_related(records)
    hdb = HomologyDb(source=":memory:")
    for rel_type, data in groups.items():
        hdb.add_records(records=data, relationship_type=rel_type)
    return hdb


def test_group_related(hom_records):
    orths = [r for r in hom_records if r["relationship"] == "ortholog_one2one"]
    related = grouped_related(orths)
    # the lambda is essential!
    got = sorted(
        related["ortholog_one2one"],
        key=lambda x: len(x),  # pylint: disable=unnecessary-lambda
        reverse=True,
    )
    expect = [
        homolog_group(relationship="ortholog_one2one", gene_ids=set("123")),
        homolog_group(
            relationship="ortholog_one2one",
            gene_ids=set("45"),
        ),
    ]
    assert got == expect


def test_homology_db(hom_hdb):
    got = sorted(
        hom_hdb.get_related_groups("ortholog_one2one"),
        key=lambda x: len(x),
        reverse=True,
    )

    expect = [
        homolog_group(
            relationship="ortholog_one2one",
            gene_ids={"1", "2", "3"},
        ),
        homolog_group(relationship="ortholog_one2one", gene_ids={"4", "5"}),
    ]
    assert got == expect


def test_species_genes_eq():
    a = species_genes(species="a", gene_ids=["1"])
    b = species_genes(species="a", gene_ids=["1"])
    assert a == b
    c = species_genes(species="a", gene_ids=["2"])
    assert a != c
    d = species_genes(species="b", gene_ids=["1"])
    assert a != d
    e = species_genes(species="b", gene_ids=["2"])
    assert a != e


def test_species_genes_hash():
    a = species_genes(species="a", gene_ids=["1"])
    assert hash(a) == hash("a")


@pytest.mark.parametrize("table_name", tuple(HomologyDb.table_names))
def test_indexing(o2o_db, table_name):
    db, _ = o2o_db
    # we just get the first column for the table
    col = HomologyDb._index_columns[table_name][0]
    expect = ("index", col, table_name)
    db.make_indexes()
    sql_template = (
        f"SELECT * FROM sqlite_master WHERE type = 'index' AND "  # nosec B608
        f"tbl_name = {table_name!r} and name = {col!r}"  # nosec B608
    )

    result = db._execute_sql(sql_template).fetchone()
    got = tuple(result)[:3]
    assert got == expect


@pytest.mark.parametrize("sp,geneids", (("abc", ()), ("abc", ["a", "b"])))
def test_species_genes_pickle_roundtrip(sp, geneids):
    import pickle  # nosec B403

    orig = species_genes(species=sp, gene_ids=geneids)
    got = pickle.loads(pickle.dumps(orig))  # nosec B301
    assert got == orig


def test_homolog_group_pickle_roundtrip():
    import pickle  # nosec B403

    orig = homolog_group(relationship="one2one", gene_ids={"1", "2", "3"})
    got = pickle.loads(pickle.dumps(orig))  # nosec B301
    assert got == orig


def test_homolog_group_union():
    a = homolog_group(relationship="one2one", gene_ids={"1", "2", "3"})
    b = homolog_group(relationship="one2one", gene_ids={"3", "4"})
    c = a | b
    assert c.gene_ids == {"1", "2", "3", "4"}


def test_homolog_group_union_invalid():
    a = homolog_group(relationship="one2one", gene_ids={"1", "2", "3"})
    b = homolog_group(relationship="one2many", gene_ids={"3", "4"})
    with pytest.raises(ValueError):
        _ = a | b


def test_merge_grouped():
    a1 = homolog_group(relationship="one2one", gene_ids={"1", "2", "3"})
    a2 = homolog_group(relationship="one2many", gene_ids={"3", "5", "6"})
    c = homolog_group(relationship="one2one", gene_ids={"3", "2", "4"})
    got = merge_grouped({"one2one": (a1,), "one2many": (a2,)}, {"one2one": (c,)})
    expect = {"one2one": (a1 | c,), "one2many": (a2,)}
    assert got == expect


def test_homdb_add_invalid_record():
    hom_db = HomologyDb()
    records = (
        homolog_group(relationship="one2one", gene_ids={"1", "2", "3"}),
        homolog_group(relationship="one2many", gene_ids={"3", "5", "6"}),
    )

    with pytest.raises(ValueError):
        hom_db.add_records(records=records, relationship_type="one2one")


@pytest.mark.parametrize(
    "gene_id,rel_type",
    (
        ("blah", "ortholog_one2one"),
        ("ENSMMUG00000065353", "ortholog_one2many"),
    ),
)
def test_homdb_get_related_to_non(o2o_db, gene_id, rel_type):
    db, _ = o2o_db
    assert not db.get_related_to(gene_id=gene_id, relationship_type=rel_type)
