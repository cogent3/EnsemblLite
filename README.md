# ensembl_cli

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

  exports the rc directory to the nominated path

  setting an environment variable ENSEMBLDBRC with this path will force its
  contents to override the default ensembl_cli settings

Options:
  -o, --outpath PATH  path to directory to export all rc contents
  --help              Show this message and exit.

```
<!-- [[[end]]] -->

