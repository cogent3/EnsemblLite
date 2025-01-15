import pytest
from cogent3 import load_table

from ensembl_tui import _homology as eti_homology
from ensembl_tui import _ingest_homology as homol_ingest


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
    src_cols = (
        "homology_type",
        "species",
        "gene_stable_id",
        "homology_species",
        "homology_gene_stable_id",
    )
    dest_col = (
        "relationship",
        "species_1",
        "gene_id_1",
        "species_2",
        "gene_id_2",
    )
    table = load_table(raw).get_columns(src_cols)
    table = table.with_new_header(src_cols, dest_col)
    table = table.get_columns(["relationship", "gene_id_1", "gene_id_2"])
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

    loader = homol_ingest.load_homologies(species)
    agg = homol_ingest.make_homology_aggregator_db()
    hom_groups = loader(raw)  # pylint: disable=not-callable
    for rel_type, data in hom_groups.items():
        agg.add_records(records=data, relationship_type=rel_type)
    agg.commit()
    homol_ingest.write_homology_views(agg=agg, outdir=tmp_dir)
    homdb = eti_homology.HomologyDb(source=tmp_dir)
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
    assert got.gene_ids.keys() == expect[gene_id]


@pytest.fixture
def orth_records():
    return [
        ("ortholog_one2one", {"1": "sp1", "2": "sp2"}),  # grp 1
        ("ortholog_one2one", {"2": "sp2", "3": "sp3"}),  # grp 1
        ("ortholog_one2one", {"4": "sp1", "5": "sp3"}),
    ]


@pytest.fixture
def hom_records(orth_records):
    return orth_records + [("ortholog_one2many", {"6": "sp2", "7": "sp3"})]


def test_hdb_get_related_groups(o2o_db):
    homdb, _ = o2o_db
    got = homdb.get_related_groups(relationship_type="ortholog_one2one")
    assert len(got) == 5


@pytest.fixture
def hom_hdb(hom_records, tmp_dir):
    groups = homol_ingest.grouped_related(hom_records)
    agg = homol_ingest.make_homology_aggregator_db()
    for rel_type, data in groups.items():
        agg.add_records(records=data, relationship_type=rel_type)
    agg.commit()
    homol_ingest.write_homology_views(agg=agg, outdir=tmp_dir)
    return eti_homology.HomologyDb(source=tmp_dir)


def test_group_related(hom_records):
    orths = [r for r in hom_records if r[0] == "ortholog_one2one"]
    related = homol_ingest.grouped_related(orths)
    # the lambda is essential!
    got = sorted(
        related["ortholog_one2one"],
        key=lambda x: len(x),  # pylint: disable=unnecessary-lambda
        reverse=True,
    )
    expect = [
        eti_homology.homolog_group(
            relationship="ortholog_one2one",
            gene_ids={
                "1": "sp1",
                "2": "sp2",
                "3": "sp3",
            },
        ),
        eti_homology.homolog_group(
            relationship="ortholog_one2one",
            gene_ids={"4": "sp1", "5": "sp3"},
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
        eti_homology.homolog_group(
            relationship="ortholog_one2one",
            gene_ids={
                "1": "sp1",
                "2": "sp2",
                "3": "sp3",
            },
        ),
        eti_homology.homolog_group(
            relationship="ortholog_one2one",
            gene_ids={
                "4": "sp1",
                "5": "sp3",
            },
        ),
    ]
    assert got == expect


def test_homolog_group_pickle_roundtrip():
    import pickle  # nosec B403

    orig = eti_homology.homolog_group(
        relationship="one2one",
        gene_ids={
            "1": "sp1",
            "2": "sp2",
            "3": "sp3",
        },
    )
    got = pickle.loads(pickle.dumps(orig))  # nosec B301
    assert got == orig


def test_homolog_group_union():
    a = eti_homology.homolog_group(
        relationship="one2one",
        gene_ids={
            "1": "sp1",
            "2": "sp2",
            "3": "sp3",
        },
    )
    b = eti_homology.homolog_group(
        relationship="one2one",
        gene_ids={
            "3": "sp3",
            "4": "sp1",
        },
    )
    c = a | b
    assert c.gene_ids == {"1": "sp1", "2": "sp2", "3": "sp3", "4": "sp1"}


def test_homolog_group_union_invalid():
    a = eti_homology.homolog_group(
        relationship="one2one",
        gene_ids={"1", "2", "3"},
    )
    b = eti_homology.homolog_group(
        relationship="one2many",
        gene_ids={"3", "4"},
    )
    with pytest.raises(ValueError):
        _ = a | b


def test_homdb_add_invalid_record():
    agg = homol_ingest.make_homology_aggregator_db()
    records = (
        eti_homology.homolog_group(
            relationship="one2one",
            gene_ids={
                "1": "sp1",
                "2": "sp2",
                "3": "sp3",
            },
        ),
        eti_homology.homolog_group(
            relationship="one2many",
            gene_ids={
                "3": "sp3",
                "5": "sp4",
                "6": "sp4",
            },
        ),
    )

    with pytest.raises(ValueError):
        agg.add_records(records=records, relationship_type="one2one")


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


def test_load_homologies(DATA_DIR):
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

    loader = homol_ingest.load_homologies(species)
    got = loader(DATA_DIR / "one2one_homologies.tsv")  # pylint: disable=not-callable
    assert len(got["ortholog_one2one"]) == 5


def test_homdb_get_related_to(o2o_db):
    homdb, _ = o2o_db
    got = homdb.get_related_to(
        gene_id="ENSG00000198786",
        relationship_type="ortholog_one2one",
    )
    assert len(got) >= 9


def test_homdb_get_related_groups(o2o_db):
    homdb, _ = o2o_db
    got = homdb.get_related_groups(relationship_type="ortholog_one2one")
    assert len(got) == 5


def test_homdb_num_records(o2o_db):
    homdb, _ = o2o_db
    got = homdb.count_distinct(homology_type=True)
    assert got.shape == (1, 2)
    assert got.columns["homology_type"][0] == "ortholog_one2one"
    assert got.columns["count"][0] == 41
    got = homdb.count_distinct(species=True)
    assert got.shape == (9, 2)


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
    loader = homol_ingest.load_homologies(
        {"gorilla_gorilla", "nomascus_leucogenys", "notamacropus_eugenii"},
    )
    records = []
    for result in loader.as_completed(hom_dir.glob("*.tsv.gz"), show_progress=False):
        records.extend(result.obj)
    assert len(records) == 2
