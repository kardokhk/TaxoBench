"""Shared helpers for the mock builders in mocks/scripts/."""

import os
import json
from pathlib import Path


def project_root() -> Path:
    """The repository root; data/ and the other relative paths of
    mocks/configs/pipeline_config.json resolve against it."""
    env = os.environ.get("TAXOBENCH_ROOT")
    if env:
        return Path(env).resolve()
    return Path(__file__).resolve().parents[2]


def configs_dir() -> Path:
    return Path(__file__).resolve().parent.parent / "configs"


def src_dir() -> Path:
    return Path(__file__).resolve().parent / "lib"


def load_pipeline_config():
    """Load mocks/configs/pipeline_config.json. Relative path values are resolved
    against project_root(); TAXOBENCH_<KEY> environment variables (for example
    TAXOBENCH_REFERENCE_DATABASE_DIR) override a key."""
    with open(configs_dir() / "pipeline_config.json") as fh:
        cfg = json.load(fh)
    for k, v in list(cfg.items()):
        env = os.environ.get("TAXOBENCH_" + k.upper())
        if env:
            cfg[k] = env
        elif isinstance(v, str) and k.endswith("_dir") and not os.path.isabs(v):
            cfg[k] = str(project_root() / v)
    return cfg


def load_mocks_config():
    with open(configs_dir() / "mocks.json") as fh:
        return json.load(fh)


def load_regions_config():
    with open(configs_dir() / "regions.json") as fh:
        return json.load(fh)


def add_src_to_path():
    import sys
    s = str(src_dir())
    if s not in sys.path:
        sys.path.insert(0, s)
