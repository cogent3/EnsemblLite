# this will be used to test integrated features
import pytest
from cogent3 import load_seq

from ensembl_tui import _config as elt_config
from ensembl_tui import _genome as elt_genome


@pytest.fixture
def one_genome(DATA_DIR, tmp_dir):
    cfg = elt_config.InstalledConfig(release="110", install_path=tmp_dir)
    # we're only making a genomes directory
    celegans = cfg.installed_genome("Caenorhabditis elegans")
    celegans.mkdir(parents=True, exist_ok=True)

    seqs_path = celegans / elt_genome.SEQ_STORE_NAME
    seqdb = elt_genome.SeqsDataHdf5(
        source=seqs_path,
        species=seqs_path.parent.name,
        mode="w",
    )
    input_seq = DATA_DIR / "c_elegans_WS199_shortened.fasta"
    seq = load_seq(
        input_seq,
        moltype="dna",
        label_to_name=lambda x: x.split()[0],
    )
    name = seq.name
    seqdb.add_records(records=[(name, str(seq))])
    seqdb.close()

    annot_path = celegans / elt_genome.ANNOT_STORE_NAME
    input_ann = DATA_DIR / "c_elegans_WS199_shortened.gff3"
    elt_genome.make_annotation_db((input_ann, annot_path))
    seq = load_seq(input_seq, input_ann, moltype="dna")
    elt_config.write_installed_cfg(cfg)
    return tmp_dir, seq


@pytest.mark.parametrize("make_seq_name", (False, True))
def test_get_genes(one_genome, make_seq_name):
    inst, seq = one_genome
    config = elt_config.read_installed_cfg(inst)
    species = "caenorhabditis_elegans"
    name = "WBGene00000138"
    cds_name = "B0019.1"
    if make_seq_name:
        # silly hack to make sure function applied
        make_seq_name = lambda x: x.name * 2  # noqa: E731

    gene = list(
        elt_genome.get_seqs_for_ids(
            config=config,
            species=species,
            names=[name],
            make_seq_name=make_seq_name,
        ),
    )[0]
    expect = [ft.get_slice() for ft in seq.get_features(name=f"CDS:{cds_name}")][0]
    assert gene.name == (
        cds_name * 2 if make_seq_name else f"caenorhabditis_elegans-{name}"
    )
    assert str(gene) == str(expect)


def test_installed_genomes(one_genome):
    inst, _ = one_genome
    config = elt_config.read_installed_cfg(inst)
    got = config.list_genomes()
    assert got == ["caenorhabditis_elegans"]
