import pathlib

import pytest

from ensembl_tui import _align as elt_align
from ensembl_tui import _config as elt_config


def test_installed_genome():
    cfg = elt_config.InstalledConfig(release="110", install_path="abcd")
    assert cfg.installed_genome("human") == pathlib.Path("abcd/genomes/homo_sapiens")


def test_installed_aligns():
    cfg = elt_config.InstalledConfig(release="110", install_path="abcd")
    assert cfg.aligns_path == pathlib.Path("abcd/compara/aligns")


def test_installed_homologies():
    cfg = elt_config.InstalledConfig(release="110", install_path="abcd")
    assert cfg.homologies_path == pathlib.Path("abcd/compara/homologies")


def test_read_installed(tmp_config, tmp_path):
    config = elt_config.read_config(tmp_config)
    outpath = elt_config.write_installed_cfg(config)
    got = elt_config.read_installed_cfg(outpath)
    assert str(got.installed_genome("human")) == str(
        got.install_path / "genomes/homo_sapiens",
    )


def test_installed_config_hash():
    ic = elt_config.InstalledConfig(release="11", install_path="abcd")
    assert hash(ic) == id(ic)
    v = {ic}
    assert len(v) == 1


@pytest.fixture
def installed_aligns(tmp_path):
    align_dir = tmp_path / elt_config._COMPARA_NAME / elt_config._ALIGNS_NAME
    align_dir.mkdir(parents=True, exist_ok=True)
    # make two alignment paths with similar names
    (align_dir / f"10_primates.epo.{elt_align.ALIGN_STORE_SUFFIX}").open(mode="w")
    (align_dir / f"24_primates.epo_extended.{elt_align.ALIGN_STORE_SUFFIX}").open(
        mode="w",
    )
    # and their associated HDF5 seqs
    (align_dir / f"10_primates.epo.{elt_align.GAP_STORE_SUFFIX}").open(mode="w")
    (align_dir / f"24_primates.epo_extended.{elt_align.GAP_STORE_SUFFIX}").open(
        mode="w",
    )

    return elt_config.InstalledConfig(release="11", install_path=tmp_path)


@pytest.mark.parametrize("pattern", ("10*", "1*prim*", "10_p*", "10_primates.epo"))
def test_get_alignment_path(installed_aligns, pattern):
    got = installed_aligns.path_to_alignment(pattern, elt_align.ALIGN_STORE_SUFFIX)
    assert got.name == f"10_primates.epo.{elt_align.ALIGN_STORE_SUFFIX}"


@pytest.mark.parametrize("pattern", ("10pri*", "blah-blah", ""))
def test_get_alignment_path_invalid(installed_aligns, pattern):
    assert (
        installed_aligns.path_to_alignment(pattern, elt_align.ALIGN_STORE_SUFFIX)
        is None
    )


@pytest.mark.parametrize("pattern", ("*pri*", "*epo*"))
def test_get_alignment_path_multiple(installed_aligns, pattern):
    with pytest.raises(ValueError):
        installed_aligns.path_to_alignment(pattern, elt_align.ALIGN_STORE_SUFFIX)
