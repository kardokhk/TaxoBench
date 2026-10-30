"""Shared helpers for the 0X_*.py orchestration scripts."""

import os
import json
from pathlib import Path


def project_root() -> Path:
    env = os.environ.get("TAXOBENCH_ROOT")
    if env:
        return Path(env).resolve()
    return Path(__file__).resolve().parent.parent


def configs_dir() -> Path:
    return project_root() / "configs"


def src_dir() -> Path:
    return project_root() / "src"


def load_pipeline_config():
    """Load configs/pipeline_config.json. Relative path values are resolved
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


def conda_bin() -> str:
    """The conda executable used for `conda run`: TAXOBENCH_CONDA_BIN, else
    "conda_bin" in pipeline_config.json, else `conda` from PATH."""
    import shutil
    env = os.environ.get("TAXOBENCH_CONDA_BIN")
    if env:
        return env
    try:
        with open(configs_dir() / "pipeline_config.json") as fh:
            configured = json.load(fh).get("conda_bin")
    except (OSError, ValueError):
        configured = None
    return configured or shutil.which("conda") or "conda"


def load_classifiers_config(path=None):
    """Load a classifier config. Default: configs/classifiers_default.json.
    `path` may be a file path or a bare name under configs/ (for example
    "classifiers_paper.json")."""
    if path is None:
        p = configs_dir() / "classifiers_default.json"
    else:
        p = Path(path)
        if not p.exists() and (configs_dir() / p).exists():
            p = configs_dir() / p
    with open(p) as fh:
        cfg = json.load(fh)
    # A method_conda_envs value may name a key of conda_envs in
    # pipeline_config.json (for example "qiime2"), so the environment paths
    # are set in one place.
    named = load_pipeline_config().get("conda_envs", {})

    def _resolve(v):
        if isinstance(v, dict):
            return {k: _resolve(x) for k, x in v.items()}
        return named.get(v, v) if isinstance(v, str) else v

    if "method_conda_envs" in cfg:
        cfg["method_conda_envs"] = _resolve(cfg["method_conda_envs"])
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


def ref_dict(classifier, ref_root, region, taxonova_settings=None):
    """Reference files for one classifier and region, or None when the region
    has no reference directory. Shared by 03_assign.py and
    05_realworld_classify.py."""
    from os.path import join, exists
    region_dir = join(ref_root, region)
    if not os.path.isdir(region_dir):
        # Taxonova with a model id (for example "v4") needs no reference
        # files: the package downloads and verifies the model itself.
        if classifier == "taxonova":
            add_src_to_path()
            import generate_commands
            model = generate_commands.resolve_taxonova_model(
                region, region_dir, taxonova_settings or {}, ref_root=ref_root)
            if model and not (model.endswith(".pt") or os.sep in model):
                return {"seqs": None, "taxa": None, "merged": None, "model": model}
        return None
    if region == "RRN":
        seqs = join(region_dir, "grond_clean.fasta")
        taxa = join(region_dir, "grond_clean.tsv")
        merged = seqs
    else:
        seqs = join(region_dir, "sequences.fasta")
        taxa = join(region_dir, "taxonomy.tsv")
        merged = join(region_dir, "merged.fasta")
        if not exists(merged):
            merged = None
    base = {"seqs": seqs, "taxa": taxa, "merged": merged}

    if classifier == "taxonova":
        add_src_to_path()
        import generate_commands
        base["model"] = generate_commands.resolve_taxonova_model(
            region, region_dir, taxonova_settings or {}, ref_root=ref_root)
    elif classifier == "naive-bayes":
        import glob
        models = sorted(glob.glob(join(region_dir, "naive-bayes", "*sklearn*qza")))
        base["model"] = models[0] if models else None
    elif classifier == "idtaxa":
        import glob
        m = (sorted(glob.glob(join(region_dir, "idtaxa", "*.rds"))) or
             sorted(glob.glob(join(region_dir, "idtaxa", "*.RData"))))
        base["model"] = m[0] if m else None
    elif classifier == "spingo":
        import glob
        m = sorted(glob.glob(join(region_dir, "spingo", "*.fasta")))
        base["database"] = m[0] if m else None
    elif classifier == "kraken2":
        kdir = join(region_dir, "kraken2")
        ok = all(exists(join(kdir, f)) for f in ("hash.k2d", "opts.k2d", "taxo.k2d"))
        base["database"] = kdir if ok else None
    elif classifier == "mapseq":
        mdir = join(region_dir, "mapseq")
        fa = join(mdir, "reference.fasta")
        tx = join(mdir, "reference.tax")
        base["database_fasta"] = fa if exists(fa) else None
        base["database_tax"]   = tx if exists(tx) else None
    return base
