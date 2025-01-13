import configparser
import pathlib

import pytest

from ensembl_tui import _align as eti_align
from ensembl_tui import _config as eti_config
from ensembl_tui import _download as eti_download
from ensembl_tui import _util as eti_util


def test_installed_genome():
    cfg = eti_config.InstalledConfig(release="110", install_path="abcd")
    assert cfg.installed_genome("human") == pathlib.Path("abcd/genomes/homo_sapiens")


def test_installed_aligns():
    cfg = eti_config.InstalledConfig(release="110", install_path="abcd")
    assert cfg.aligns_path == pathlib.Path("abcd/compara/aligns")


def test_installed_homologies():
    cfg = eti_config.InstalledConfig(release="110", install_path="abcd")
    assert cfg.homologies_path == pathlib.Path("abcd/compara/homologies")


def test_read_installed(tmp_config, tmp_path):
    config = eti_config.read_config(tmp_config)
    outpath = eti_config.write_installed_cfg(config)
    got = eti_config.read_installed_cfg(outpath)
    assert str(got.installed_genome("human")) == str(
        got.install_path / "genomes/homo_sapiens",
    )


def test_installed_config_hash():
    ic = eti_config.InstalledConfig(release="11", install_path="abcd")
    assert hash(ic) == id(ic)
    v = {ic}
    assert len(v) == 1


@pytest.fixture
def installed_aligns(tmp_path):
    align_dir = tmp_path / eti_config._COMPARA_NAME / eti_config._ALIGNS_NAME
    align_dir.mkdir(parents=True, exist_ok=True)
    # make two alignment paths with similar names
    (align_dir / f"10_primates.epo.{eti_align.ALIGN_STORE_SUFFIX}").open(mode="w")
    (align_dir / f"24_primates.epo_extended.{eti_align.ALIGN_STORE_SUFFIX}").open(
        mode="w",
    )
    # and their associated HDF5 seqs
    (align_dir / f"10_primates.epo.{eti_align.GAP_STORE_SUFFIX}").open(mode="w")
    (align_dir / f"24_primates.epo_extended.{eti_align.GAP_STORE_SUFFIX}").open(
        mode="w",
    )

    return eti_config.InstalledConfig(release="11", install_path=tmp_path)


@pytest.mark.parametrize("pattern", ("10*", "1*prim*", "10_p*", "10_primates.epo"))
def test_get_alignment_path(installed_aligns, pattern):
    got = installed_aligns.path_to_alignment(pattern, eti_align.ALIGN_STORE_SUFFIX)
    assert got.name == f"10_primates.epo.{eti_align.ALIGN_STORE_SUFFIX}"


@pytest.mark.parametrize("pattern", ("10pri*", "blah-blah", ""))
def test_get_alignment_path_invalid(installed_aligns, pattern):
    assert (
        installed_aligns.path_to_alignment(pattern, eti_align.ALIGN_STORE_SUFFIX)
        is None
    )


@pytest.mark.parametrize("pattern", ("*pri*", "*epo*"))
def test_get_alignment_path_multiple(installed_aligns, pattern):
    with pytest.raises(ValueError):
        installed_aligns.path_to_alignment(pattern, eti_align.ALIGN_STORE_SUFFIX)


@pytest.fixture
def empty_cfg(tmp_dir, ENSEMBL_RELEASE_VERSION):
    parser = configparser.ConfigParser()
    parser.read(eti_util.get_resource_path("sample.cfg"))
    parser.remove_section("Caenorhabditis elegans")
    parser.remove_section("Saccharomyces cerevisiae")
    parser.remove_section("compara")
    parser.set("local path", "staging_path", value=str(tmp_dir / "staging"))
    parser.set("local path", "install_path", value=str(tmp_dir / "install"))
    parser.set("release", "release", value=ENSEMBL_RELEASE_VERSION)
    return tmp_dir, parser


@pytest.fixture
def cfg_just_aligns(empty_cfg):
    tmp_dir, parser = empty_cfg
    parser.add_section("compara")
    parser.set("compara", "align_names", value="10_primates.epo")
    download_cfg = tmp_dir / "download.cfg"
    with open(download_cfg, "w") as out:  # noqa: PTH123
        parser.write(out)

    return download_cfg


COMMON_NAMES = (
    "Crab-eating macaque",
    "Human",
    "Bonobo",
    "Chimpanzee",
    "Mouse Lemur",
    "Gorilla",
    "Gibbon",
    "Vervet-AGM",
    "Sumatran orangutan",
    "Macaque",
)


@pytest.fixture
def cfg_just_genomes(empty_cfg):
    tmp_dir, parser = empty_cfg
    parser.add_section("compara")
    parser.set("compara", "align_names", value="10_primates.epo")
    download_cfg = tmp_dir / "download.cfg"
    for name in COMMON_NAMES:
        parser.add_section(name)
        parser.set(name, "db", value="core")

    with open(download_cfg, "w") as out:  # noqa: PTH123
        parser.write(out)

    return download_cfg
    # we make a config using common names


@pytest.mark.internet
def test_read_config_compara_genomes(cfg_just_aligns):
    from ensembl_tui._species import Species

    config = eti_config.read_config(cfg_just_aligns)
    assert not config.species_dbs
    sp = eti_download.get_species_for_alignments(
        host=config.host,
        remote_path=config.remote_path,
        release=config.release,
        align_names=config.align_names,
    )
    expected = {Species.get_species_name(n) for n in COMMON_NAMES}
    assert set(sp.keys()) == expected


@pytest.mark.internet
def test_read_config_genomes(cfg_just_genomes):
    from ensembl_tui._species import Species

    config = eti_config.read_config(cfg_just_genomes)
    expected = {Species.get_species_name(n) for n in COMMON_NAMES}
    assert set(config.species_dbs.keys()) == expected
