"""TaxoBench 2: reference-aware benchmarking of taxonomic classifiers.

The data model (contract.py) is the interface between tier builders, adapters
and the scorer.
"""
__version__ = "2.0.0"

RANKS = ("domain", "phylum", "class", "order", "family", "genus", "species")
PREFIX = {"domain": "d__", "phylum": "p__", "class": "c__", "order": "o__",
          "family": "f__", "genus": "g__", "species": "s__"}
