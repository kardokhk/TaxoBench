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


# ---------------------------------------------------------------------------
# Fail-closed run guards (2026-10-31). A rerun used to leave the previous
# taxonomy_final.tsv in place when the classifier failed, exit 0, and let
# 04_evaluate.py score the old file as the new run.

#: Files a classifier command writes into its output directory. They are
#: removed before a run so that nothing from an earlier run can survive a
#: failed one.
RUN_OUTPUTS = (
    "taxonomy_final.tsv", "taxonomy_final.tsv.failed", "taxonomy_raw.tsv", "taxonomy_raw.mseq",
    "taxonomy_classification_out.qza", "taxonomy_search_out.qza",
    "kraken2_classification.txt", "kraken2_report.txt",
    "kraken2_unclassified.fasta", "kraken2_classified.fasta",
)
RUN_STATUS = "run.status.json"
_HEADER_IDS = ("id", "Feature ID", "Feature.ID", "taxonomy", "Taxon", "#OTU ID")


def clear_run_outputs(outdir):
    """Delete the outputs of any earlier run in `outdir` (and its status)."""
    outdir = Path(outdir)
    removed = []
    for name in RUN_OUTPUTS + (RUN_STATUS,):
        fp = outdir / name
        if fp.exists() or fp.is_symlink():
            fp.unlink()
            removed.append(name)
    return removed


def fasta_ids(fp):
    """Sequence IDs of a FASTA file, in order."""
    out = []
    with open(fp) as fh:
        for line in fh:
            if line.startswith(">"):
                out.append(line[1:].split()[0])
    return out


def table_ids(fp):
    """First-column IDs of a taxonomy TSV, header rows skipped, in order."""
    out = []
    with open(fp) as fh:
        for line in fh:
            if not line.strip():
                continue
            first = line.rstrip("\n").split("\t")[0]
            if first in _HEADER_IDS:
                continue
            out.append(first)
    return out


def coverage_problems(got_ids, want_ids):
    """Empty list when `got_ids` covers `want_ids` exactly, once each;
    otherwise human-readable problems (duplicates, missing, unexpected)."""
    got_ids = list(got_ids)
    got, want = set(got_ids), set(want_ids)
    probs = []
    dup = len(got_ids) - len(got)
    if dup:
        probs.append(f"{dup} duplicated ID(s)")
    if want - got:
        ex = sorted(want - got)[:3]
        probs.append(f"{len(want - got)} of {len(want)} query IDs missing (e.g. {ex})")
    if got - want:
        ex = sorted(got - want)[:3]
        probs.append(f"{len(got - want)} unexpected ID(s) (e.g. {ex})")
    return probs


def write_run_status(outdir, **fields):
    """Write run.status.json beside a run's outputs."""
    import time
    fields.setdefault("written", time.time())
    with open(Path(outdir) / RUN_STATUS, "w") as fh:
        json.dump(fields, fh, indent=1, sort_keys=True)


def read_run_status(outdir):
    fp = Path(outdir) / RUN_STATUS
    if not fp.exists():
        return None
    with open(fp) as fh:
        return json.load(fh)


def check_run(outdir, query_fasta, started, returncode):
    """Problems with a finished run (empty list when it is usable): non-zero
    exit, no taxonomy_final.tsv, a taxonomy_final.tsv older than the start of
    the run, or an output that does not cover the query IDs exactly."""
    outdir = Path(outdir)
    final = outdir / "taxonomy_final.tsv"
    probs = []
    if returncode != 0:
        probs.append(f"exit status {returncode}")
    if not final.exists():
        probs.append("no taxonomy_final.tsv written")
        return probs
    if final.stat().st_mtime < started - 1:
        probs.append("taxonomy_final.tsv predates this run (stale)")
    if query_fasta is not None and os.path.exists(query_fasta):
        probs += coverage_problems(table_ids(final), fasta_ids(query_fasta))
    return probs


def run_guarded(cmd, outdir, log_fp, query_fasta):
    """Clear old outputs, run `cmd` in a shell, check the result and write
    run.status.json. Returns (cmd, ok, problems)."""
    import subprocess
    import time
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    clear_run_outputs(outdir)
    Path(log_fp).parent.mkdir(parents=True, exist_ok=True)
    started = time.time()
    with open(log_fp, "w") as fh:
        ret = subprocess.run(cmd, shell=True, stdout=fh, stderr=subprocess.STDOUT)
    probs = check_run(outdir, query_fasta, started, ret.returncode)
    final = outdir / "taxonomy_final.tsv"
    if probs and final.exists():
        # Keep a failed run's output for diagnosis, but not under the name the
        # evaluator reads.
        final.rename(outdir / "taxonomy_final.tsv.failed")
    write_run_status(outdir, returncode=ret.returncode, started=started,
                     finished=time.time(), ok=not probs, problems=probs,
                     query_fasta=str(query_fasta) if query_fasta else None)
    return cmd, not probs, probs
