"""Snakemake scheduling for `taxobench run --snakemake`.

Snakemake only schedules: each job is `taxobench run --spec <spec.json>`,
which does the fail-closed work (and the Apptainer or conda wrapping named in
the spec). On a cluster, running Snakemake inside one node allocation packs
many small runs into its cores.
"""
from __future__ import annotations

import json
import shutil
import subprocess
from importlib import resources
from pathlib import Path


def snakefile() -> Path:
    return Path(str(resources.files("taxobench").joinpath("Snakefile")))


def write_specs(specs, workdir) -> Path:
    wd = Path(workdir)
    (wd / "specs").mkdir(parents=True, exist_ok=True)
    for s in specs:
        (wd / "specs" / f"{s['run_id']}.json").write_text(json.dumps(s, indent=1))
    cfg = {"runs": {s["run_id"]: {"spec": str((wd / "specs" / f"{s['run_id']}.json").resolve()),
                                  "outdir": str(Path(s["outdir"]).resolve()),
                                  "threads": int(s.get("threads", 1))} for s in specs}}
    fp = wd / "config.json"
    fp.write_text(json.dumps(cfg, indent=1))
    return fp


def run_snakemake(specs, workdir, cores=1, extra=()) -> int:
    exe = shutil.which("snakemake")
    if exe is None:
        print("snakemake not found; install taxobench[workflow]")
        return 1
    cfg = write_specs(specs, workdir)
    cmd = [exe, "-s", str(snakefile()), "--configfile", str(cfg), "--cores", str(cores),
           "--directory", str(Path(workdir).resolve()), "--rerun-incomplete", "--keep-going", *extra]
    print(" ".join(cmd))
    return subprocess.run(cmd).returncode
