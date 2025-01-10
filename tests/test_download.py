import pytest

import ensembl_tui._config as eti_config
import ensembl_tui._download as eti_download
from ensembl_tui import _mysql_core_attr as eti_db_attr


@pytest.mark.internet
def test_get_db_names():
    cfg = eti_config.read_config(eti_download.DEFAULT_CFG)
    db_names = eti_download.get_core_db_dirnames(cfg)
    assert db_names == {
        "caenorhabditis_elegans": "pub/release-112/mysql/caenorhabditis_elegans_core_112_282",
        "saccharomyces_cerevisiae": "pub/release-112/mysql/saccharomyces_cerevisiae_core_112_4",
    }


def test_make_dumpfiles():
    table_names = set(eti_db_attr.make_mysqldump_names())
    assert "CHECKSUMS" in table_names
    assert {
        n.split(".")[0] for n in table_names - {"CHECKSUMS"}
    } == eti_db_attr.get_all_tables()
