"""`taxobench run`: execute classifier runs and write contract tables.

A run spec (JSON) describes one run:

    {
      "run_id": "taxonova-26.11.2__gg2-2024.09-V4__default",   # optional: derived
      "adapter": "taxonova",               # parser, see taxobench.adapters
      "classifier": "taxonova", "version": "26.11.2",
      "ref_id": "gg2-2024.09/V4",
      "taxa": "<tier dir with taxa.tsv>" | "<id<TAB>lineage file>",
      "query_fasta": "queries.fasta",
      "outdir": "runs/<run_id>",
      "command": ["taxonova", "classify", "{query}", "--model", "{model}", "-o", "{outdir}/raw.tsv"],
      "raw": "raw.tsv"  (or {"output": ..., "report": ...} relative to outdir),
      "model": "path/to/model.pt",         # optional; sha256 and md5 recorded
      "parameters": {...}, "parse_params": {...},
      "runner": {"type": "local"} | {"type": "apptainer", "image": "x.sif", "binds": [...]}
               | {"type": "conda", "prefix": "/env", "exe": "micromamba"},
      "threads": 8
    }

With no "command", the raw output must already exist (ingesting precomputed or
archived outputs). The run is fail-closed: the output directory is cleared
first, a non-zero exit, a missing raw output or predictions that do not cover
the query IDs exactly once stop it, and nothing is written as a result unless
every check passed (results are moved into place last).
"""
from __future__ import annotations

import hashlib
import json
import os
import platform
import shlex
import shutil
import subprocess
import time
from pathlib import Path

import pandas as pd

from . import __version__
from .adapters import get_adapter
from .contract import ContractError, read_table, sha256_file, write_table
from .taxonomy import Reference, load_reference_lineages

RESULT_FILES = ("predictions.tsv", "predictions.tsv.json", "runs.tsv", "runs.tsv.json", "run.json")


def fasta_ids(fp):
    with open(fp) as fh:
        return [ln[1:].split()[0] for ln in fh if ln.startswith(">")]


def _md5(fp):
    h = hashlib.md5()
    with open(fp, "rb") as fh:
        for c in iter(lambda: fh.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def load_reference(spec) -> Reference:
    p = Path(spec["taxa"])
    if p.is_dir() or p.name.startswith("taxa"):
        taxa = read_table(p, "taxa", verify=spec.get("verify", True))
        return Reference.from_taxa(taxa, spec["ref_id"])
    return load_reference_lineages(p, spec["ref_id"])


def derive_run_id(spec) -> str:
    params = json.dumps(spec.get("parameters", {}), sort_keys=True)
    model = spec.get("model_sha256", "")
    h = hashlib.sha256(f"{spec['adapter']}|{spec['version']}|{spec['ref_id']}|{params}|{model}".encode()).hexdigest()[:10]
    return f"{spec.get('classifier', spec['adapter'])}-{spec['version']}__{spec['ref_id'].replace('/', '-')}__{h}"


def _wrap(argv, runner, outdir):
    t = (runner or {}).get("type", "local")
    if t == "local":
        return argv
    if t == "apptainer":
        binds = []
        for b in runner.get("binds", []) + [str(outdir)]:
            binds += ["--bind", b]
        return [runner.get("exe", "apptainer"), "exec", "--cleanenv", *binds, runner["image"], *argv]
    if t == "conda":
        return [runner.get("exe", os.environ.get("TAXOBENCH_CONDA_BIN", "conda")), "run", "-p",
                runner["prefix"], *argv]
    raise ValueError(f"unknown runner type {t!r}")


def run_one(spec: dict, log=print) -> Path:
    spec = dict(spec)
    if spec.get("model"):
        spec["model_sha256"] = sha256_file(spec["model"])
        spec["model_md5"] = _md5(spec["model"])
    spec.setdefault("classifier", spec["adapter"])
    run_id = spec.get("run_id") or derive_run_id(spec)
    outdir = Path(spec.get("outdir") or Path("runs") / run_id)
    outdir.mkdir(parents=True, exist_ok=True)
    for name in RESULT_FILES:
        (outdir / name).unlink(missing_ok=True)
    adapter = get_adapter(spec["adapter"])
    query = Path(spec["query_fasta"])
    qids = fasta_ids(query)
    if len(qids) != len(set(qids)):
        raise ContractError(f"{query}: duplicated query IDs")
    status = {"run_id": run_id, "taxobench_version": __version__, "host": platform.node(),
              "start": time.time(), "spec": spec, "adapter": adapter.describe()}
    if spec.get("command"):
        raw_names = spec["raw"] if isinstance(spec["raw"], dict) else {"raw": spec["raw"]}
        for v in raw_names.values():
            (outdir / v).unlink(missing_ok=True)
        fmt = {"query": str(query), "outdir": str(outdir), "threads": str(spec.get("threads", 1)),
               "model": str(spec.get("model", "")), **{k: str(v) for k, v in spec.get("parameters", {}).items()}}
        argv = [a.format(**fmt) for a in spec["command"]]
        argv = _wrap(argv, spec.get("runner"), outdir)
        status["argv"] = argv
        log("run: " + " ".join(shlex.quote(a) for a in argv))
        with open(outdir / "run.log", "w") as fh:
            ret = subprocess.run(argv, stdout=fh, stderr=subprocess.STDOUT,
                                 stdin=subprocess.DEVNULL)
        status["returncode"] = ret.returncode
        if ret.returncode != 0:
            status["ok"] = False
            (outdir / "run.failed.json").write_text(json.dumps(status, indent=1, default=str))
            raise RuntimeError(f"{run_id}: command exited {ret.returncode} (log {outdir/'run.log'})")
    raw = spec["raw"]
    raw_path = {k: outdir / v for k, v in raw.items()} if isinstance(raw, dict) else outdir / raw
    if isinstance(raw, str) and Path(raw).is_absolute():
        raw_path = Path(raw)
    for p in (raw_path.values() if isinstance(raw_path, dict) else [raw_path]):
        if not Path(p).exists():
            raise ContractError(f"{run_id}: raw output {p} missing")
        if spec.get("command") and Path(p).stat().st_mtime < status["start"] - 1:
            raise ContractError(f"{run_id}: raw output {p} predates the run (stale)")
    status["raw_sha256"] = {str(p): sha256_file(p) for p in
                            (raw_path.values() if isinstance(raw_path, dict) else [raw_path])}
    parsed = adapter.parse(raw_path, **spec.get("parse_params", {}))
    ref = load_reference(spec)
    preds = adapter.to_predictions(parsed, ref, run_id, query_ids=qids)
    runs = pd.DataFrame([{
        "run_id": run_id, "classifier": spec["classifier"], "version": spec["version"],
        "ref_id": spec["ref_id"], "adapter": adapter.name, "adapter_version": adapter.version,
        "confidence_kind": adapter.confidence_kind, "model_md5": spec.get("model_md5", ""),
        "model_sha256": spec.get("model_sha256", ""),
        "parameters": json.dumps(spec.get("parameters", {}), sort_keys=True),
        "thresholds": json.dumps(spec.get("thresholds", {}), sort_keys=True),
        "command": " ".join(shlex.quote(a) for a in status.get("argv", [])),
        "date": time.strftime("%Y-%m-%dT%H:%M:%S"), "code_commit": spec.get("code_commit", ""),
        "taxobench_version": __version__,
    }])
    tmp = outdir / ".tmp"
    if tmp.exists():
        shutil.rmtree(tmp)
    tmp.mkdir()
    write_table(preds, tmp / "predictions.tsv", "predictions")
    write_table(runs, tmp / "runs.tsv", "runs")
    counts = preds.groupby(["rank", "call"]).size().unstack(fill_value=0).to_dict("index")
    status.update(end=time.time(), ok=True, n_queries=len(qids), calls=counts)
    (tmp / "run.json").write_text(json.dumps(status, indent=1, default=str))
    for name in RESULT_FILES:
        (tmp / name).rename(outdir / name)
    tmp.rmdir()
    (outdir / "run.failed.json").unlink(missing_ok=True)
    log(f"{run_id}: {len(qids)} queries; oov at species "
        f"{counts.get('species', {}).get('oov', 0)}; written to {outdir}")
    return outdir


def expand_plan(plan: dict) -> list:
    """A plan is {"defaults": {...}, "runs": [spec, ...]}; each spec inherits
    the defaults. Run ids are derived where absent."""
    out = []
    for s in plan["runs"]:
        spec = dict(plan.get("defaults", {}), **s)
        spec.setdefault("classifier", spec["adapter"])
        if not spec.get("run_id"):
            spec["run_id"] = derive_run_id(spec)
        spec.setdefault("outdir", str(Path(plan.get("outroot", "runs")) / spec["run_id"]))
        out.append(spec)
    ids = [s["run_id"] for s in out]
    if len(ids) != len(set(ids)):
        raise ContractError("plan has duplicated run ids")
    return out


def collect(run_dirs, out_dir) -> Path:
    """Merge the per-run tables into one runs.tsv and predictions.tsv."""
    runs, preds = [], []
    for d in run_dirs:
        runs.append(read_table(d, "runs"))
        preds.append(read_table(d, "predictions"))
    out = Path(out_dir)
    write_table(pd.concat(runs, ignore_index=True), out / "runs.tsv", "runs")
    write_table(pd.concat(preds, ignore_index=True), out / "predictions.tsv", "predictions")
    return out
