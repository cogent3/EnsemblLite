"""Ensembl cli tools"""

from warnings import filterwarnings

filterwarnings("ignore", message=".*MPI")
filterwarnings("ignore", message="Can't drop database.*")


__version__ = "0.1a8"
