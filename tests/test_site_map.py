import pytest

from ensembl_cli._site_map import get_site_map


@pytest.mark.parametrize("site", ("ftp.ensembl.org", "ftp.ensemblgenomes.ebi.ac.uk"))
def test_correct_site(site):
    smap = get_site_map(site)
    assert smap.site == site


def test_standard_smp():
    sm = get_site_map("ftp.ensembl.org")
    assert str(sm._genomes_path) == "fasta/dna"
