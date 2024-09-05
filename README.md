[![CI](https://github.com/cogent3/ensembl_tui/actions/workflows/testing_develop.yml/badge.svg)](https://github.com/cogent3/ensembl_tui/actions/workflows/testing_develop.yml)
[![CodeQL](https://github.com/cogent3/ensembl_tui/actions/workflows/codeql.yml/badge.svg)](https://github.com/cogent3/ensembl_tui/actions/workflows/codeql.yml)
[![Coverage Status](https://coveralls.io/repos/github/cogent3/ensembl_tui/badge.svg?branch=develop)](https://coveralls.io/github/cogent3/ensembl_tui?branch=develop)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)

# ensembl-tui

ensembl-tui provides the `eti` terminal application for obtaining a subset of the data provided by Ensembl which can then be queried locally. You can have multiple such subsets on your machine, each corresponding to a different selection of species and data types.

> **Warning**
> ensembl-tui is in a preliminary phase of development with a limited feature set and incomplete test coverage! We currently **only support accessing data from the main ensembl.org** site. Please validate results against the web version. If you discover errors, please post a [bug report](https://github.com/cogent3/ensembl_tui/issues).

## Installing the software

<details>
  <summary>General user installation instructions</summary>

  ```
  $ pip install ensembl-tui
  ```

</details>

<details>
  <summary>Developer installation instructions</summary>
  Fork the repo and clone your fork to your local machine. In the terminal, create either a python virtual environment or a new conda environment and activate it. In that virtual environment

  ```
  $ pip install flit
  ```

  Then do the flit version of a "developer install". (It is basically creating a symlink to the repos source directory.)

  ```
  $ flit install -s --python `which python`
  ```
</details>

## Resources required to subset Ensembl data

Ensembl hosts some very large data sets. You need to have a machine with sufficient disk space to store the data you want to download. At present we do not have support for predicting how much storage would be required for a given selection of species and data types. You will need to experiment.

Some commands can be run in parallel but have moderate memory requirements. If you have a machine with limited RAM, you may need to reduce the number of parallel processes. Again, run some experiments.

## Getting setup

<details>
  <summary>Specifying what data you want to download and where to put it</summary>

  We use a plain text file to indicate the Ensembl domain, release and types of genomic data to download. Start by using the `exportrc` subcommand.

  <!-- [[[cog
  import cog
  from ensembl_tui import cli
  from click.testing import CliRunner
  runner = CliRunner()
  result = runner.invoke(cli.main, ["exportrc", "--help"])
  help = result.output.replace("Usage: main", "Usage: eti")
  cog.out(
      "```\n{}\n```".format(help)
  )
  ]]] -->
  ```
  Usage: eti exportrc [OPTIONS]

    exports sample config and species table to the nominated path

  Options:
    -o, --outpath PATH  Path to directory to export all rc contents.
    --help              Show this message and exit.

  ```
  <!-- [[[end]]] -->

  ```shell
  $ eti exportrc -o ~/Desktop/Outbox/ensembl_download
  ```
  This command creates a `ensembl_download` download directory and writes two plain text files into it:

  1. `species.tsv`: contains the Latin names, common names etc... of the species accessible at ensembl.org website.
  2. `sample.cfg`: a sample configuration file that you can edit to specify the data you want to download.

  The latter file includes comments on how to edit it in order to specify the genomic resources that you want.
</details>

<details>
  <summary>Downloading the data</summary>
  Downloads the data indicated in the config file to a local directory.

  <!-- [[[cog
  import cog
  from ensembl_tui import cli
  from click.testing import CliRunner
  runner = CliRunner()
  result = runner.invoke(cli.main, ["download", "--help"])
  help = result.output.replace("Usage: main", "Usage: eti")
  cog.out(
      "```\n{}\n```".format(help)
  )
  ]]] -->
  ```
  Usage: eti download [OPTIONS]

    download data from Ensembl's ftp site

  Options:
    -c, --configpath PATH  Path to config file specifying databases, (only species
                           or compara at present).
    -d, --debug            Maximum verbosity, and reduces number of downloads,
                           etc...
    -v, --verbose
    --help                 Show this message and exit.

  ```
  <!-- [[[end]]] -->

  For a config file named `config.cfg`, the download command would be:

  ```shell
  $ cd to/directory/with/config.cfg
  $ eti download -c config.cfg
  ```

  > **Note**
  > Downloads can be interrupted and resumed. The software deletes partially downloaded files.

The download creates a new `.cfg` file inside the download directory. This file is used by the `install` command.

</details>

<details>
  <summary>Installing the data</summary>
  
  <!-- [[[cog
  import cog
  from ensembl_tui import cli
  from click.testing import CliRunner
  runner = CliRunner()
  result = runner.invoke(cli.main, ["install", "--help"])
  help = result.output.replace("Usage: main", "Usage: eti")
  cog.out(
      "```\n{}\n```".format(help)
  )
  ]]] -->
  ```
  Usage: eti install [OPTIONS]

    create the local representations of the data

  Options:
    -d, --download PATH       Path to local download directory containing a cfg
                              file.
    -np, --num_procs INTEGER  Number of procs to use.  [default: 1]
    -f, --force_overwrite     Overwrite existing data.
    -v, --verbose
    --help                    Show this message and exit.

  ```
  <!-- [[[end]]] -->

The following command uses 2 CPUs and has been safe on systems with only 16GB of RAM for 10 primate genomes, including homology data and whole genome:

```shell
$ cd to/directory/with/downloaded_data
$ eti install -d downloaded_data -np 2
```

</details>

<details>
  <summary>Checking what has been installed</summary>
  
  <!-- [[[cog
  import cog
  from ensembl_tui import cli
  from click.testing import CliRunner
  runner = CliRunner()
  result = runner.invoke(cli.main, ["installed", "--help"])
  help = result.output.replace("Usage: main", "Usage: eti")
  cog.out(
      "```\n{}\n```".format(help)
  )
  ]]] -->
  ```
  Usage: eti installed [OPTIONS]

    show what is installed

  Options:
    -i, --installed TEXT  Path to root directory of an installation.  [required]
    --help                Show this message and exit.

  ```
  <!-- [[[end]]] -->

</details>

## Interrogating the data

We provide a conventional command line interface for querying the data with subcommands.

<details>
  <summary>The full list of subcommands</summary>

  You can get help on individual subcommands by running `eti <subcommand>` in the terminal.

  <!-- [[[cog
  import cog
  from ensembl_tui import cli
  from click.testing import CliRunner
  runner = CliRunner()
  result = runner.invoke(cli.main)
  help = result.output.replace("Usage: main", "Usage: eti")
  cog.out(
      "```\n{}\n```".format(help)
  )
  ]]] -->
  ```
  Usage: eti [OPTIONS] COMMAND [ARGS]...

    Tools for obtaining and interrogating subsets of https://ensembl.org genomic
    data.

  Options:
    --version  Show the version and exit.
    --help     Show this message and exit.

  Commands:
    alignments       export multiple alignments in fasta format for named genes
    download         download data from Ensembl's ftp site
    dump-genes       export meta-data table for genes from one species to...
    exportrc         exports sample config and species table to the nominated...
    homologs         exports CDS sequence data in fasta format for homology...
    install          create the local representations of the data
    installed        show what is installed
    species-summary  genome summary data for a species
    tui              Open Textual TUI.

  ```
  <!-- [[[end]]] -->

</details>

We also provide an experiment terminal user interface (TUI) that allows you to explore the data in a more interactive way. This is invoked with the `tui` subcommand.
