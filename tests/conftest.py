import pathlib

from configparser import ConfigParser

import pytest

from ensembl_cli.util import get_resource_path


@pytest.fixture(scope="session")
def DATA_DIR():
    return pathlib.Path(__file__).parent / "data"


@pytest.fixture(scope="function")
def tmp_dir(tmp_path_factory):
    return tmp_path_factory.mktemp("cli")


@pytest.fixture(scope="function")
def tmp_config(tmp_dir):
    # create a simpler download config
    # we want a very small test set
    parser = ConfigParser()
    parser.read(get_resource_path("ensembldb_download.cfg"))
    parser.remove_section("Caenorhabditis elegans")
    parser.remove_section("compara")
    parser.set("local path", "staging_path", value=str(tmp_dir / "staging"))
    parser.set("local path", "install_path", value=str(tmp_dir / "install"))
    download_cfg = tmp_dir / "download.cfg"
    with open(download_cfg, "wt") as out:
        parser.write(out)

    yield download_cfg
