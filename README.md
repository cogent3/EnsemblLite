# ensembl_cli

## Installation

Suggest creating a conda environment or a python virtual environment, using python3.11. Then install directly into that environment from the GitHub repo as

```
$ python -m pip install "ensembl_cli  @ git+https://github.com/cogent3/ensembl_cli.git@develop"
```

Then run for the first time using

```
$ ensembl_cli tui
```

The first start takes a while as, behind the scenes, cogent3 is transpiling various functions into C and compiling them. Eventually, you get a very neat terminal interface you can click around in. To exit, make sure the "root" is selected on the left panel then `^+r`.

## Usage

The setup is (for now) controlled using a config file, defined in `ini` format. To get a starting template use the `exportrc` subcommand.

<!-- [[[cog
import cog
from ensembl_cli import cli
from click.testing import CliRunner
runner = CliRunner()
result = runner.invoke(cli.main, ["exportrc", "--help"])
help = result.output.replace("Usage: main", "Usage: ensembl_cli")
cog.out(
    "```\n{}\n```".format(help)
)
]]] -->
```
Usage: ensembl_cli exportrc [OPTIONS]

  exports sample config and species table to the nominated path

  setting an environment variable ENSEMBLDBRC with this path will force its
  contents to override the default ensembl_cli settings

Options:
  -o, --outpath PATH  path to directory to export all rc contents
  --help              Show this message and exit.

```
<!-- [[[end]]] -->

### Download

The `download` downloads, for the species indicated in the config file:

- genomes sequences as fasta format 
- annotations as gff3
- gene homologies for individual genomes in tsv format

Alignments indicated in the config file will be downloaded in `.maf` format.

Downloads are written to a local directory, specified in the config file. Downloads are done in parallel (using threads).

### Install

The downloaded files are (or will be) transformed into local sqlite3 databases which are written to a location specified in the config file.

From the maf alignment files, the "ancestral" sequences are discarded and for every aligned sequence only the gap data is stored (i.e. gap position and length) along with the genomic coordinates. These alignments will be reconstructable by combining this information with the whole genome sequence. (This approach reduces storage requirements ~5-fold).

Installation is done in parallel on multiple CPUs (since the data need to be decompressed on the fly).