from configparser import ConfigParser
from random import shuffle

import pytest

from ensembl_tui import _config as eti_config
from ensembl_tui import _util as eti_util


@pytest.fixture(scope="function")
def compara_cfg(tmp_config):
    # we just add compara sections
    parser = ConfigParser()
    parser.read(eti_util.get_resource_path(tmp_config))
    parser.add_section("compara")
    alns = ",".join(("17_sauropsids.epc", "10_primates.epo"))
    parser.set("compara", "align_names", value=alns)
    with open(tmp_config, "w") as out:
        parser.write(out)

    return tmp_config


def test_parse_config(compara_cfg):
    cfg = eti_config.read_config(compara_cfg)
    assert set(cfg.align_names) == {"17_sauropsids.epc", "10_primates.epo"}


def test_load_ensembl_md5sum(DATA_DIR):
    got = eti_util.load_ensembl_md5sum(DATA_DIR / "sample-MD5SUM")
    assert len(got) == 3
    assert got["b.emf.gz"] == "3d9af835d9ed19975bd8b2046619a3a1"


def test_load_ensembl_checksum(DATA_DIR):
    got = eti_util.load_ensembl_checksum(DATA_DIR / "sample-CHECKSUMS")
    assert len(got) == 4  # README line is ignored
    assert got["c.fa.gz"] == (7242, 327577)


@pytest.fixture(scope="function")
def gorilla_cfg(tmp_config):
    # we add gorilla genome
    parser = ConfigParser()
    parser.read(eti_util.get_resource_path(tmp_config))
    parser.add_section("Gorilla")
    parser.set("Gorilla", "db", value="core")
    with open(tmp_config, "w") as out:
        parser.write(out)

    return tmp_config


def test_parse_config_gorilla(gorilla_cfg):
    # Gorilla has two synonyms, we need only one
    cfg = eti_config.read_config(gorilla_cfg)
    num_gorilla = sum(1 for k in cfg.species_dbs if "gorilla" in k)
    assert num_gorilla == 1


@pytest.mark.parametrize(
    "name",
    (
        "Gallus_gallus.bGalGall.mat.broiler.GRCg7b.dna_rm.primary_assembly.MT.fa.gz",
        "Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna_rm.primary_assembly.Z.fa.gz",
        "Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna_rm.toplevel.fa.gz",
        "Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna_sm.nonchromosomal.fa.gz",
        "Homo_sapiens.GRCh38.dna_rm.alt.fa.gz",
        "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz",
        "Homo_sapiens.GRCh38.dna.toplevel.fa.gz",
    ),
)
def test_invalid_seq(name):
    from ensembl_tui._download import valid_seq_file

    assert not valid_seq_file(name)


@pytest.mark.parametrize(
    "name",
    (
        "Homo_sapiens.GRCh38.dna.chromosome.Y.fa.gz",
        "Homo_sapiens.GRCh38.dna.nonchromosomal.fa.gz",
        "Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.primary_assembly.W.fa.gz",
        "Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.nonchromosomal.fa.gz",
        "Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.primary_assembly.MT.fa.gz",
    ),
)
def test_valid_seq(name):
    from ensembl_tui._download import valid_seq_file

    assert valid_seq_file(name)


@pytest.fixture(scope="function")
def just_compara_cfg(tmp_config):
    # no genomes!
    parser = ConfigParser()
    parser.read(tmp_config)
    parser.remove_section("Saccharomyces cerevisiae")
    parser.add_section("compara")
    parser.set("compara", "align_names", value="10_primates.epo")
    parser.set("compara", "tree_names", value="10_primates_EPO_default.nh")
    with open(tmp_config, "w") as out:
        parser.write(out)

    return tmp_config


@pytest.mark.internet
def test_just_compara(just_compara_cfg):
    # get species names from the alignment ref tree
    cfg = eti_config.read_config(just_compara_cfg)
    # 10 primates i the alignments, so we should have 10 db's
    assert len(cfg.species_dbs) == 10


def test_write_read_installed_config(tmp_config):
    config = eti_config.read_config(tmp_config)
    cfg_path = eti_config.write_installed_cfg(config)
    icfg = eti_config.read_installed_cfg(cfg_path.parent)
    assert icfg.release == config.release
    assert icfg.install_path == config.install_path


def test_match_align_tree(tmp_config):
    trees = [
        "pub/release-110/compara/species_trees/16_pig_breeds_EPO-Extended_default.nh",
        "pub/release-110/compara/species_trees/21_murinae_EPO_default.nh",
        "pub/release-110/compara/species_trees/39_fish_EPO_default.nh",
        "pub/release-110/compara/species_trees/65_amniota_vertebrates_Mercator-Pecan_default.nh",
    ]

    aligns = [
        "pub/release-110/maf/ensembl-compara/multiple_alignments/16_pig_breeds.epo_extended",
        "pub/release-110/maf/ensembl-compara/multiple_alignments/21_murinae.epo",
        "pub/release-110/maf/ensembl-compara/multiple_alignments/39_fish.epo",
        "pub/release-110/maf/ensembl-compara/multiple_alignments/65_amniotes.pecan",
    ]

    expect = dict(zip(aligns, trees, strict=False))
    shuffle(aligns)
    result = eti_util.trees_for_aligns(aligns, trees)
    assert result == expect


def test_missing_match_align_tree(tmp_config):
    trees = [
        "pub/release-110/compara/species_trees/16_pig_breeds_EPO-Extended_default.nh",
        "pub/release-110/compara/species_trees/21_murinae_EPO_default.nh",
        "pub/release-110/compara/species_trees/65_amniota_vertebrates_Mercator-Pecan_default.nh",
    ]

    aligns = [
        "pub/release-110/maf/ensembl-compara/multiple_alignments/16_pig_breeds.epo_extended",
        "pub/release-110/maf/ensembl-compara/multiple_alignments/21_murinae.epo",
        "pub/release-110/maf/ensembl-compara/multiple_alignments/39_fish.epo",
        "pub/release-110/maf/ensembl-compara/multiple_alignments/65_amniotes.pecan",
    ]
    with pytest.raises(ValueError):
        eti_util.trees_for_aligns(aligns, trees)


def test_config_update_invalid_species(tmp_config):
    config = eti_config.read_config(tmp_config)
    with pytest.raises(ValueError):
        config.update_species({"Micro bat": ["core"]})


def test_config_update_species(tmp_config):
    config = eti_config.read_config(tmp_config)
    config.update_species({"Human": ["core"]})
    assert len(list(config.db_names)) == 2
    assert set(config.db_names) == {"homo_sapiens", "saccharomyces_cerevisiae"}


@pytest.mark.internet
def test_cfg_to_dict(just_compara_cfg):
    cfg = eti_config.read_config(just_compara_cfg)
    data = cfg.to_dict()
    cfg.write()
    path = cfg.staging_path / eti_config.DOWNLOADED_CONFIG_NAME
    assert path.exists()
    got_cfg = eti_config.read_config(path)
    assert got_cfg.to_dict() == data


def test_blosc_apps():
    o = "ACGG" * 1000
    z = eti_util.eti_compress_it(o)
    assert isinstance(z, bytes)
    assert len(z) < len(o)
    assert eti_util.eti_decompress_it(z) == o


def test_get_sig_calc_func_invalid():
    with pytest.raises(NotImplementedError):
        eti_util.get_sig_calc_func(2)


def test_is_signature():
    assert not eti_util.is_signature("blah")


def test_exec_command():
    got = eti_util.exec_command("ls")
    assert isinstance(got, str)


def test_exec_command_fail(capsys):
    with pytest.raises(SystemExit):
        eti_util.exec_command("qwertyuiop")

    _ = capsys.readouterr()


@pytest.mark.parametrize("biotype", ("gene", "exon"))
def test_sanitise_stableid(biotype):
    identifier = "ENSG00012"
    stableid = f"{biotype}:{identifier}"
    got = eti_util.sanitise_stableid(stableid)
    assert got == identifier


@pytest.mark.parametrize("text", ["'primate'", '"primate"'])
def test_stripquotes(text):
    assert eti_util.strip_quotes(text) == "primate"


def test_unique_values():
    indexer = eti_util.unique_value_indexer()
    index = indexer("a")
    assert index == 1
    index = indexer("b")
    assert index == 2
    assert indexer("a") == 1


def test_unique_values_iter():
    indexer = eti_util.unique_value_indexer()
    indexer("a")
    indexer("b")
    got = list(indexer)
    assert got == [(1, "a"), (2, "b")]


def test_unique_values_tuples():
    indexer = eti_util.unique_value_indexer()
    i1 = indexer((1, "a"))
    i2 = indexer((2, "a"))
    i3 = indexer((2, "b"))
    assert indexer((1, "a")) == i1 == 1
    assert indexer((2, "a")) == i2 == 2
    assert indexer((2, "b")) == i3 == 3


def test_unique_values_tuple_iter():
    indexer = eti_util.unique_value_indexer()
    indexer((1, "a"))
    indexer((2, "a"))
    indexer((2, "b"))
    got = list(indexer)
    assert got == [(1, (1, "a")), (2, (2, "a")), (3, (2, "b"))]


@pytest.mark.parametrize(
    ("data", "expected"),
    [
        (
            (("cat1", {"a", "b"}), ("cat1", {"c", "d"})),
            [1, 2],
        ),  # diff indices
        (
            (("cat1", {"a", "b"}), ("cat1", {"c", "b"})),
            [1, 1],
        ),  # same indices
        (
            (("cat1", {"a", "b"}), ("cat2", {"a", "b"})),
            [1, 2],
        ),  # diff indices
        (
            (
                ("cat1", {"a", "b"}),
                ("cat2", {"a", "b"}),
                ("cat1", {"c", "b"}),
            ),
            [1, 2, 1],
        ),  # mix indices
    ][-1:],
)
def test_category_indexer(data, expected):
    indexer = eti_util.category_indexer()
    got = [indexer(cat, grp) for cat, grp in data]
    assert got == expected


def test_category_indexer_empty_vals():
    indexer = eti_util.category_indexer()
    with pytest.raises(ValueError):
        indexer("cat1", set())


def test_category_indexer_iter():
    indexer = eti_util.category_indexer()
    data = [
        ("c1", {"a", "b"}),
        ("c2", {"a", "b"}),
        ("c1", {"c", "b"}),
    ]
    expect = {
        (1, "c1", "a"),
        (1, "c1", "b"),
        (1, "c1", "c"),
        (2, "c2", "a"),
        (2, "c2", "b"),
    }
    _ = [indexer(cat, grp) for cat, grp in data]
    got = set(indexer)
    assert got == expect
