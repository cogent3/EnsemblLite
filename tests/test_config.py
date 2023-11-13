import pathlib

from ensembl_lite._config import InstalledConfig


def test_installed_genome():
    cfg = InstalledConfig(release=110, install_path="abcd")
    assert cfg.installed_genome("human") == pathlib.Path("abcd/genome/homo_sapiens")


def test_installed_aligns():
    cfg = InstalledConfig(release=110, install_path="abcd")
    assert cfg.install_aligns == pathlib.Path("abcd/compara/aligns")


def test_installed_homologies():
    cfg = InstalledConfig(release=110, install_path="abcd")
    assert cfg.install_homologies == pathlib.Path("abcd/compara/homologies")
