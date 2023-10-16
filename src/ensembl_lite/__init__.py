from warnings import filterwarnings

from .species import Species


filterwarnings("ignore", message=".*MPI")
filterwarnings("ignore", message="Can't drop database.*")


__all__ = [
    "name",
    "species",
    "util",
    "Species",
]

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2023-, The ensembl cli project"
__credits__ = ["Gavin Huttley"]
__license__ = "BSD"
__version__ = "2023.7.9a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "alpha"
