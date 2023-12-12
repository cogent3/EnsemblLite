[![CI](https://github.com/cogent3/EnsemblLite/actions/workflows/testing_develop.yml/badge.svg)](https://github.com/cogent3/EnsemblLite/actions/workflows/testing_develop.yml)
[![CodeQL](https://github.com/cogent3/EnsemblLite/actions/workflows/codeql.yml/badge.svg)](https://github.com/cogent3/EnsemblLite/actions/workflows/codeql.yml)
[![Coverage Status](https://coveralls.io/repos/github/cogent3/EnsemblLite/badge.svg?branch=develop)](https://coveralls.io/github/cogent3/EnsemblLite?branch=develop)

# EnsemblLite

> **Warning**
> EnsemblLite is not ready for use! We will remove this notice when we are ready to post to PyPi at which point it will be ready for trialling. In the meantime, you can check the project progress towards being usable via the [EnsemblLite roadmap](https://github.com/cogent3/EnsemblLite/issues/38).

## A screencast of an early prototype

<details>
  <summary> 🎬 Very early proof-of-concept demo and plan for a new style terminal user interface </summary>
    <video src="https://user-images.githubusercontent.com/3102996/273427137-d3835f8b-8c0a-4370-a6e1-f8805f5dd320.mp4" controls="controls" style="max-height:640px">
    </video>
    
**NOTE:** the command line name has changed since this early version. See text below for the new name.

</details>

## Developer installs

Fork the repo and clone your fork to your local machine. In the terminal, create either a python virtual environment or a new conda environment and activate it. In that virtual environment

```
$ pip install flit
```

Then do the flit version of a "developer install". (It is basically creating a symlink to the repos source directory.)

```
$ flit install -s --python `which python`
```

## Installation

Suggest creating a conda environment or a python virtual environment, using python3.11. Then install directly into that environment from the GitHub repo as

```
$ python -m pip install "ensembl_lite @ git+https://github.com/cogent3/EnsemblLite.git@develop"
```

Then run for the first time using

```
$ elt tui
```

The first start takes a while as, behind the scenes, cogent3 is transpiling various functions into C and compiling them. Eventually, you get a very neat terminal interface you can click around in. To exit, make sure the "root" is selected on the left panel then `^+r`.

## Usage

The setup is (for now) controlled using a config file, defined in `ini` format. To get a starting template use the `exportrc` subcommand.

<!-- [[[cog
import cog
from ensembl_lite import cli
from click.testing import CliRunner
runner = CliRunner()
result = runner.invoke(cli.main, ["exportrc", "--help"])
help = result.output.replace("Usage: main", "Usage: elt")
cog.out(
    "```\n{}\n```".format(help)
)
]]] -->
```
Usage: elt exportrc [OPTIONS]

  exports sample config and species table to the nominated path

  setting an environment variable ENSEMBLDBRC with this path will force its
  contents to override the default ensembl_lite settings

Options:
  -o, --outpath PATH  path to directory to export all rc contents
  --help              Show this message and exit.

```
<!-- [[[end]]] -->

<details>
    <summary> Click to see a sample config file I've been using for development </summary>
    
Using this config, it takes approximately 16' to download (over a ~200MB/s WiFi connection) and ~45' to install on my M2 Macbook Pro (note the install is incomplete). (Note this step uses up to  10 CPU cores.)

```
[remote path]
host=ftp.ensembl.org
path=pub
[local path]
staging_path=~/Desktop/Outbox/ensembl_download
install_path=~/Desktop/Outbox/ensembl_install
[release]
release=110
[Mouse Lemur]
db=core
[Macaque]
db=core
[Gibbon]
db=core
[Orangutan]
db=core
[Bonobo]
db=core
[Human]
db=core
[Chimp]
db=core
[Gorilla]
db=core
[compara]
align_names=10_primates.epo
```
</details>

### Download

Downloads the species indicated in the config file:

- genomes sequences as fasta format 
- annotations as gff3
- gene homologies for individual genomes in tsv format

Alignments indicated in the config file will be downloaded in `.maf` format.

Downloads are written to a local directory, specified in the config file. Downloads are done in parallel (using threads).

### Install

"Installation" presently involves transforming downloaded files into local sqlite3 databases which are saved to the location specified in the config file.

From the maf alignment files, the "ancestral" sequences are discarded and for every aligned sequence only the gap data is stored (i.e. gap position and length) along with the genomic coordinates. These alignments will be reconstructable by combining this information with the whole genome sequence. (This approach reduces storage requirements ~5-fold).

Installation is done in parallel on multiple CPUs (since the data need to be decompressed on the fly).