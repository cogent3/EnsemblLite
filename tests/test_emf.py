import pytest

from ensembl_tui import _emf as eti_emf
from ensembl_tui import _name as eti_name


def test_load(DATA_DIR):
    path = DATA_DIR / "sample.emf"
    got = list(eti_emf.parse_emf(path))[0]
    expect = {
        eti_name.EmfName(
            "human",
            "4",
            "450000",
            "560000",
            "1",
            "(chr_length=201709)",
        ): "-TCGC",
        eti_name.EmfName(
            "mouse",
            "17",
            "780000",
            "790000",
            "-1",
            "(chr_length=201709)",
        ): "AT--G",
        eti_name.EmfName(
            "rat",
            "12",
            "879999",
            "889998",
            "1",
            "(chr_length=201709)",
        ): "AAA--",
    }
    assert got == expect


def test_unsupported_format(tmp_path, DATA_DIR):
    data = (DATA_DIR / "sample.emf").read_text().splitlines(keepends=True)
    data[0] = data[0].replace("compara", "resequencing")
    outpath = tmp_path / "sample.emf"
    outpath.write_text("".join(data))
    with pytest.raises(NotImplementedError):
        list(eti_emf.parse_emf(outpath))
