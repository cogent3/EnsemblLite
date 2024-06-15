import pathlib

import pytest

from ensembl_lite._aligndb import _GAP_STORE_SUFFIX
from ensembl_lite._config import (
    _ALIGNS_NAME,
    _COMPARA_NAME,
    InstalledConfig,
    read_config,
    read_installed_cfg,
    write_installed_cfg,
)


def test_installed_genome():
    cfg = InstalledConfig(release=110, install_path="abcd")
    assert cfg.installed_genome("human") == pathlib.Path("abcd/genomes/homo_sapiens")


def test_installed_aligns():
    cfg = InstalledConfig(release=110, install_path="abcd")
    assert cfg.aligns_path == pathlib.Path("abcd/compara/aligns")


def test_installed_homologies():
    cfg = InstalledConfig(release=110, install_path="abcd")
    assert cfg.homologies_path == pathlib.Path("abcd/compara/homologies")


def test_read_installed(tmp_config, tmp_path):
    config = read_config(tmp_config)
    outpath = write_installed_cfg(config)
    got = read_installed_cfg(outpath)
    assert str(got.installed_genome("human")) == str(
        got.install_path / "genomes/homo_sapiens"
    )


def test_installed_config_hash():
    ic = InstalledConfig(release="11", install_path="abcd")
    assert hash(ic) == id(ic)
    v = {ic}
    assert len(v) == 1


@pytest.fixture
def installed_aligns(tmp_path):
    align_dir = tmp_path / _COMPARA_NAME / _ALIGNS_NAME
    align_dir.mkdir(parents=True, exist_ok=True)
    # make two alignment paths with similar names
    (align_dir / "10_primates.epo.sqlitedb").open(mode="w")
    (align_dir / "24_primates.epo_extended.sqlitedb").open(mode="w")
    # and their associated HDF5 seqs
    (align_dir / f"10_primates.epo.{_GAP_STORE_SUFFIX}").open(mode="w")
    (align_dir / f"24_primates.epo_extended.{_GAP_STORE_SUFFIX}").open(mode="w")

    return InstalledConfig(release="11", install_path=tmp_path)


@pytest.mark.parametrize("pattern", ("10*", "1*prim*", "10_p*", "*s.epo.*"))
def test_get_alignment_path(installed_aligns, pattern):
    got = installed_aligns.path_to_alignment(pattern)
    assert got.name == "10_primates.epo.sqlitedb"


@pytest.mark.parametrize("pattern", ("10pri*", "blah-blah", ""))
def test_get_alignment_path_invalid(installed_aligns, pattern):
    assert installed_aligns.path_to_alignment(pattern) is None


@pytest.mark.parametrize("pattern", ("*pri*", "*epo*"))
def test_get_alignment_path_multiple(installed_aligns, pattern):
    with pytest.raises(ValueError):
        installed_aligns.path_to_alignment(pattern)
