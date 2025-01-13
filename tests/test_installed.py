# this will be used to test integrated features


from ensembl_tui import _config as eti_config
from ensembl_tui import _genome as eti_genome


def test_load_genome(small_install_path):
    config = eti_config.read_installed_cfg(small_install_path)
    species = "caenorhabditis_elegans"
    genome = eti_genome.load_genome(config=config, species=species)
    assert genome.species == species
    # directly interrogate the gene view
    stable_id = "WBGene00004893"
    gene = list(
        genome.get_cds(
            biotype="protein_coding",
            stable_id=stable_id,
        ),
    )
    assert len(gene) == 1


def test_get_genes(small_install_path):
    config = eti_config.read_installed_cfg(small_install_path)
    species = "caenorhabditis_elegans"
    name = "WBGene00004893"

    gene = next(
        iter(
            eti_genome.get_seqs_for_ids(
                config=config,
                species=species,
                names=[name],
            ),
        ),
    )
    # we check we can make a aa seq which has the correct start and end
    aa = str(gene.trim_stop_codon().get_translation(incomplete_ok=True))
    # expected values from ensembl.org
    assert aa.startswith("MTNSSEFTDVLQS")
    assert aa.endswith("TIMNRINYKLQ")


def test_installed_genomes(small_install_path):
    config = eti_config.read_installed_cfg(small_install_path)
    got = config.list_genomes()
    assert got == ["caenorhabditis_elegans", "saccharomyces_cerevisiae"]
