import pytest

from ensembl_tui._name import EnsemblDbName


def test_cmp_name():
    """should validly compare names by attributes"""
    n1 = EnsemblDbName("homo_sapiens_core_46_36h")
    n2 = EnsemblDbName("homo_sapiens_core_46_36h")
    assert n1 == n2


def test_cmp_name_false():
    """should validly compare names by attributes"""
    n1 = EnsemblDbName("homo_sapiens_core_46_36h")
    n2 = EnsemblDbName("homo_sapiens_core_46_37h")
    assert n1 != n2


def test_name_without_build():
    """should correctly handle a db name without a build"""
    n = EnsemblDbName("pongo_pygmaeus_core_49_1")
    assert n.prefix == "pongo_pygmaeus"
    assert n.type == "core"
    assert n.build == "1"


def test_lt():
    """should validly compare names by attributes"""
    n1 = EnsemblDbName("homo_sapiens_core_46_36h")
    n2 = EnsemblDbName("homo_sapiens_core_46_37h")
    assert n1 < n2


def test_species_with_three_words_name():
    """should correctly parse a db name that contains a three words species name"""
    n = EnsemblDbName("mustela_putorius_furo_core_70_1")
    assert n.prefix == "mustela_putorius_furo"
    assert n.type == "core"
    assert n.build == "1"
    assert n.species == "Mustela putorius furo"
    n = EnsemblDbName("canis_lupus_familiaris_core_102_31")
    assert n.species == "Canis lupus familiaris"


def test_ensemblgenomes_names():
    """correctly handle the ensemblgenomes naming system"""
    n = EnsemblDbName("aedes_aegypti_core_5_58_1e")
    assert n.prefix == "aedes_aegypti"
    assert n.type == "core"
    assert n.release == "5"
    assert n.general_release == "58"
    assert n.build == "1e"
    n = EnsemblDbName("ensembl_compara_metazoa_6_59")
    assert n.release == "6"
    assert n.general_release == "59"
    assert n.type == "compara"


def test_invalid_name():
    with pytest.raises(ValueError):
        EnsemblDbName("ab cd _fg")


def test_repr():
    db = "aedes_aegypti_core_5_58_1e"
    n = EnsemblDbName(db)
    assert repr(n).startswith("db(prefix=")


def test_str():
    db = "aedes_aegypti_core_5_58_1e"
    n = EnsemblDbName(db)
    assert str(n) == db
