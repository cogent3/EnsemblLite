import pathlib

from ensembl_lite._config import (
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
