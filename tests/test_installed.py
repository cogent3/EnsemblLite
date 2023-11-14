# this will be used to test integrated features
import pytest

from cogent3 import load_annotations, load_seq

from ensembl_lite._config import (
    InstalledConfig,
    read_installed_cfg,
    write_installed_cfg,
)
from ensembl_lite._genomedb import (
    _ANNOTDB_NAME,
    _SEQDB_NAME,
    CompressedGenomeSeqsDb,
    get_seq_by_ids,
)


@pytest.fixture
def one_genome(DATA_DIR, tmp_dir):
    cfg = InstalledConfig(release="110", install_path=tmp_dir)
    # we're only making a genomes directory
    celegans = cfg.installed_genome("Caenorhabditis elegans")
    celegans.mkdir(parents=True, exist_ok=True)

    seqs_path = celegans / _SEQDB_NAME
    seqdb = CompressedGenomeSeqsDb(source=seqs_path, species=seqs_path.parent.name)
    input_seq = DATA_DIR / "c_elegans_WS199_shortened.fasta"
    seq = load_seq(
        input_seq,
        moltype="dna",
        label_to_name=lambda x: x.split()[0],
    )
    name = seq.name
    seqdb.add_records(records=[(name, str(seq))])
    seqdb.close()

    annot_path = celegans / _ANNOTDB_NAME
    input_ann = DATA_DIR / "c_elegans_WS199_shortened.gff3"
    ann_db = load_annotations(path=input_ann, write_path=annot_path)
    ann_db.db.close()
    seq = load_seq(input_seq, input_ann, moltype="dna")
    write_installed_cfg(cfg)
    return tmp_dir, seq


def test_installed_genomes(one_genome):
    inst, seq = one_genome
    config = read_installed_cfg(inst)
    got = config.list_genomes()
    assert got == ["caenorhabditis_elegans"]
