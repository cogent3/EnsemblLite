import pytest

from ensembl_tui._site_map import get_site_map


@pytest.mark.parametrize("site", ("ftp.ensembl.org",))
def test_correct_site(site):
    smap = get_site_map(site)
    assert smap.site == site


def test_standard_smp():
    sm = get_site_map("ftp.ensembl.org")
    assert sm.get_seqs_path("abcd") == "fasta/abcd/dna"
    assert sm.get_annotations_path("abcd") == "gff3/abcd"
