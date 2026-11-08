"""Classifier adapters (see base.py). Third-party adapters register under the
`taxobench.adapters` entry-point group."""
from .base import REGISTRY, Adapter, Parsed, get_adapter  # noqa: F401
