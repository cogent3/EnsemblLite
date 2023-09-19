from configparser import ConfigParser
from random import shuffle

import pytest

from ensembl_cli._config import (
    read_config,
    read_installed_cfg,
    write_installed_cfg,
)
from ensembl_cli.util import (
    get_resource_path,
    load_ensembl_checksum,
    load_ensembl_md5sum,
    trees_for_aligns,
)


@pytest.fixture(scope="function")
def compara_cfg(tmp_config):
    # we just add compara sections
    parser = ConfigParser()
    parser.read(get_resource_path(tmp_config))
    parser.add_section("compara")
    alns = ",".join(("17_sauropsids.epc", "10_primates.epo"))
    parser.set("compara", "align_names", value=alns)
    with open(tmp_config, "wt") as out:
        parser.write(out)

    yield tmp_config


def test_parse_config(compara_cfg):
    from ensembl_cli._config import read_config

    cfg = read_config(compara_cfg)
    assert set(cfg.align_names) == {"17_sauropsids.epc", "10_primates.epo"}


def test_load_ensembl_md5sum(DATA_DIR):
    got = load_ensembl_md5sum(DATA_DIR / "sample-MD5SUM")
    assert len(got) == 3
    assert got["b.emf.gz"] == "3d9af835d9ed19975bd8b2046619a3a1"


def test_load_ensembl_checksum(DATA_DIR):
    got = load_ensembl_checksum(DATA_DIR / "sample-CHECKSUMS")
    assert len(got) == 4  # README line is ignored
    assert got["c.fa.gz"] == (7242, 327577)


@pytest.fixture(scope="function")
def gorilla_cfg(tmp_config):
    # we add gorilla genome
    parser = ConfigParser()
    parser.read(get_resource_path(tmp_config))
    parser.add_section("Gorilla")
    parser.set("Gorilla", "db", value="core")
    with open(tmp_config, "wt") as out:
        parser.write(out)

    yield tmp_config


def test_parse_config_gorilla(gorilla_cfg):
    from ensembl_cli._config import read_config

    # Gorilla has two synonyms, we need only one
    cfg = read_config(gorilla_cfg)
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
    from ensembl_cli.download import valid_seq_file

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
    from ensembl_cli.download import valid_seq_file

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
    with open(tmp_config, "wt") as out:
        parser.write(out)

    yield tmp_config


def test_just_compara(just_compara_cfg):
    # get species names from the alignment ref tree
    cfg = read_config(just_compara_cfg)
    # 10 primates i the alignments, so we should have 10 db's
    assert len(cfg.species_dbs) == 10


def test_write_read_installed_config(tmp_config):
    config = read_config(tmp_config)
    cfg_path = write_installed_cfg(config)
    icfg = read_installed_cfg(cfg_path.parent)
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

    expect = dict(zip(aligns, trees))
    shuffle(aligns)
    result = trees_for_aligns(aligns, trees)
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
        trees_for_aligns(aligns, trees)


def test_config_update_invalid_species(tmp_config):
    config = read_config(tmp_config)
    with pytest.raises(ValueError):
        config.update_species({"Micro bat": ["core"]})


def test_config_update_species(tmp_config):
    config = read_config(tmp_config)
    config.update_species({"Human": ["core"]})
    assert len(list(config.db_names)) == 2
    assert set(config.db_names) == {"homo_sapiens", "saccharomyces_cerevisiae"}
