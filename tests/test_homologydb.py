import pytest

from cogent3 import load_table

from ensembl_lite._homologydb import (
    _HOMOLOGYDB_NAME,
    HomologyDb,
    HomologyRecordType,
    grouped_related,
    id_by_species_group,
)
from ensembl_lite.install import LoadHomologies


def _make_expected_o2o(table):
    """return dict with keys stable ID's values list[tuple[str,str]]"""
    data = table.to_list(["species_1", "gene_id_1", "species_2", "gene_id_2"])

    result = {}
    for sp1, g1, sp2, g2 in data:
        value = {(sp1, g1), (sp2, g2)}
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

    table = table.with_new_header(loader.src_cols, loader.dest_col[:-1])

    data = loader([raw])
    homdb = HomologyDb(tmp_dir / _HOMOLOGYDB_NAME)
    homdb.add_records(records=data, col_order=loader.dest_col)
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
    assert got == expect[gene_id]


@pytest.fixture
def orth_records():
    common = dict(relationship="ortholog_one2one")
    return [
        HomologyRecordType(
            species_1="a", gene_id_1="1", species_2="b", gene_id_2="2", **common
        ),  # grp 1
        HomologyRecordType(
            species_1="b", gene_id_1="2", species_2="c", gene_id_2="3", **common
        ),  # grp 1
        HomologyRecordType(
            species_1="a", gene_id_1="4", species_2="b", gene_id_2="5", **common
        ),
    ]


@pytest.fixture
def hom_records(orth_records):
    return orth_records + [
        HomologyRecordType(
            species_1="a",
            gene_id_1="6",
            species_2="e",
            gene_id_2="7",
            relationship="ortholog_one2many",
        ),
    ]


def _reduced_to_ids(groups):
    return [[i for _, i in group] for group in groups]


def test_hdb_get_related_groups(o2o_db):
    homdb, _ = o2o_db
    got = homdb.get_related_groups(relationship_type="ortholog_one2one")
    groups = _reduced_to_ids(got)
    assert len(groups) == 5


@pytest.fixture
def hom_hdb(hom_records):
    col_order = (
        "relationship",
        "species_1",
        "gene_id_1",
        "species_2",
        "gene_id_2",
    )
    records = [[r[c] for c in col_order] for r in hom_records]
    hdb = HomologyDb(source=":memory:")
    hdb.add_records(records=records, col_order=col_order)
    return hdb


def test_group_related(hom_records):
    orths = [r for r in hom_records if r["relationship"] == "ortholog_one2one"]
    got = grouped_related(orths)
    expect = {
        frozenset([("a", "1"), ("b", "2"), ("c", "3")]),
        frozenset([("a", "4"), ("b", "5")]),
    }
    assert got == expect


def test_homology_db(hom_hdb):
    got = hom_hdb.get_related_groups("ortholog_one2one")
    expect = {
        frozenset([("a", "1"), ("b", "2"), ("c", "3")]),
        frozenset([("a", "4"), ("b", "5")]),
    }
    assert got == expect


def test_grouped_by_species(hom_hdb):
    got_species, got_gene_map = id_by_species_group(
        hom_hdb.get_related_groups("ortholog_one2one")
    )
    sp_map = {sp.species: set(sp.gene_ids) for sp in got_species}
    expected_species = {"a": {"1", "4"}, "b": {"2", "5"}, "c": {"3"}}
    assert sp_map == expected_species
    expected_groups = {("1", "2", "3"), ("4", "5")}
    got_groups = {}
    for g, i in got_gene_map.items():
        got_groups[i] = got_groups.get(i, []) + [g]
    got_groups = {tuple(sorted(g)) for g in got_groups.values()}
    assert got_groups == expected_groups
