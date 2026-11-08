"""`taxobench import-mocks | run | ingest | score | report | compare | validate | adapters`.

Every subcommand exits non-zero on any failure.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from . import __version__


def _err(msg):
    print(f"taxobench: error: {msg}", file=sys.stderr)
    return 1


def cmd_import_mocks(a):
    from .mocks import build
    cells = [tuple(c.split(":", 1)) for c in a.cell] or None
    build(a.data_dir, a.ref_dir, a.out, cells=cells, tier=a.tier_name)
    return 0


def cmd_run(a):
    from .run import expand_plan, run_one
    if a.spec:
        specs = [json.loads(Path(a.spec).read_text())]
    else:
        specs = expand_plan(json.loads(Path(a.plan).read_text()))
    if a.dry_run:
        for s in specs:
            print(json.dumps(s))
        return 0
    if a.snakemake:
        from .workflow import run_snakemake
        return run_snakemake(specs, a.workdir, cores=a.cores, extra=a.snakemake_args)
    failed = []
    for s in specs:
        try:
            run_one(s)
        except Exception as e:          # noqa: BLE001  (report every failure, then exit 1)
            print(f"FAILED {s.get('run_id')}: {e}", file=sys.stderr)
            failed.append(s.get("run_id"))
    return 1 if failed else 0


def cmd_ingest(a):
    from .run import run_one
    spec = {"adapter": a.adapter, "classifier": a.classifier or a.adapter, "version": a.version,
            "ref_id": a.ref_id, "taxa": a.taxa, "query_fasta": a.query_fasta,
            "raw": str(Path(a.raw).resolve()), "outdir": a.out, "parse_params": dict(
                kv.split("=", 1) for kv in a.param), "verify": not a.no_verify}
    if a.run_id:
        spec["run_id"] = a.run_id
    run_one(spec)
    return 0


def cmd_score(a):
    from .score import Context, score, write_scores
    ctx = Context.load(a.tier, a.runs, verify=not a.no_verify)
    df = score(ctx, rules=tuple(a.group_rules), n_boot=a.n_boot, seed=a.seed, split=a.split)
    fp = write_scores(df, a.out)
    print(f"wrote {fp} ({len(df):,} rows)")
    return 0


def cmd_report(a):
    import pandas as pd
    from .report import to_markdown
    s = pd.read_csv(Path(a.scores), sep="\t")
    md = to_markdown(s, ranks=tuple(a.ranks))
    if a.out:
        Path(a.out).write_text(md)
    else:
        print(md)
    return 0


def cmd_compare(a):
    from .compare import compare, load_plan, write
    from .score import Context
    ctx = Context.load(a.tier, a.runs, verify=not a.no_verify)
    out, guard, meta = compare(ctx, a.a, a.b, plan=load_plan(a.plan), n_boot=a.n_boot,
                               seed=a.seed, split=a.split, rules=tuple(a.group_rules))
    d = write(out, guard, meta, a.out)
    print(f"wrote {d}/compare.md and compare.tsv")
    if len(guard) and not guard["passed"].all():
        print("guardrail failed", file=sys.stderr)
        return 3
    return 0


def cmd_validate(a):
    from .contract import TierTables, RunTables
    for d in a.tier or []:
        t = TierTables.load(d, verify=not a.no_verify)
        print(f"tier {d}: {len(t.queries)} queries, {len(t.truth)} truth rows, "
              f"{t.taxa['ref_id'].nunique()} references: ok")
    for d in a.runs or []:
        r = RunTables.load(d, verify=not a.no_verify)
        print(f"runs {d}: {len(r.runs)} runs, {len(r.predictions)} prediction rows: ok")
    return 0


def cmd_adapters(a):
    from .adapters.base import REGISTRY, get_adapter
    get_adapter("taxonova")
    for name in sorted(REGISTRY):
        ad = REGISTRY[name]()
        print(f"{name}\tv{ad.version}\tconfidence={ad.confidence_kind}\t{ad.confidence_note}")
    return 0


def build_parser():
    p = argparse.ArgumentParser(prog="taxobench", description=__doc__.splitlines()[0])
    p.add_argument("--version", action="version", version=f"taxobench {__version__}")
    sub = p.add_subparsers(dest="cmd", required=True)

    m = sub.add_parser("import-mocks", help="mock builder output -> tier tables and query FASTA files")
    m.add_argument("--data-dir", required=True, help="directory holding mock_community/<mock>/<region>/")
    m.add_argument("--ref-dir", required=True, help="reference root holding <region>/taxonomy.tsv "
                                                    "(RRN/grond_clean.tsv for the operon)")
    m.add_argument("--cell", action="append", default=[], help="mock:region (repeatable; default: all)")
    m.add_argument("--tier-name", default="mocks")
    m.add_argument("--out", required=True)
    m.set_defaults(func=cmd_import_mocks)

    r = sub.add_parser("run", help="run classifiers from a plan or one spec (fail closed)")
    g = r.add_mutually_exclusive_group(required=True)
    g.add_argument("--plan")
    g.add_argument("--spec")
    r.add_argument("--dry-run", action="store_true")
    r.add_argument("--snakemake", action="store_true", help="schedule with Snakemake")
    r.add_argument("--workdir", default="taxobench-work")
    r.add_argument("--cores", type=int, default=1)
    r.add_argument("--snakemake-args", nargs=argparse.REMAINDER, default=[])
    r.set_defaults(func=cmd_run)

    i = sub.add_parser("ingest", help="parse an existing raw output into contract tables")
    i.add_argument("--adapter", required=True)
    i.add_argument("--classifier")
    i.add_argument("--version", required=True)
    i.add_argument("--ref-id", required=True)
    i.add_argument("--taxa", required=True, help="tier dir with taxa.tsv, or id<TAB>lineage file")
    i.add_argument("--query-fasta", required=True)
    i.add_argument("--raw", required=True)
    i.add_argument("--run-id")
    i.add_argument("--param", action="append", default=[], help="parser parameter key=value")
    i.add_argument("--no-verify", action="store_true")
    i.add_argument("--out", required=True)
    i.set_defaults(func=cmd_ingest)

    s = sub.add_parser("score", help="score runs against tier tables")
    s.add_argument("--tier", nargs="+", required=True)
    s.add_argument("--runs", nargs="+", required=True)
    s.add_argument("--group-rules", nargs="+", default=["G1", "G2", "G3"])
    s.add_argument("--split", default="eval")
    s.add_argument("--n-boot", type=int, default=2000)
    s.add_argument("--seed", type=int, default=20261031)
    s.add_argument("--no-verify", action="store_true", help="skip sha256 sidecar checks")
    s.add_argument("--out", required=True)
    s.set_defaults(func=cmd_score)

    rp = sub.add_parser("report", help="markdown summary of a scores.tsv")
    rp.add_argument("scores")
    rp.add_argument("--ranks", nargs="+", default=["genus", "species"])
    rp.add_argument("--out")
    rp.set_defaults(func=cmd_report)

    c = sub.add_parser("compare", help="paired comparison of two runs")
    c.add_argument("--a", required=True, help="run_id of the current release")
    c.add_argument("--b", required=True, help="run_id of the candidate")
    c.add_argument("--tier", nargs="+", required=True)
    c.add_argument("--runs", nargs="+", required=True)
    c.add_argument("--plan", help="JSON with the confirmatory family and guardrails")
    c.add_argument("--group-rules", nargs="+", default=["G1", "G2"])
    c.add_argument("--split", default="eval")
    c.add_argument("--n-boot", type=int, default=10000)
    c.add_argument("--seed", type=int, default=20261031)
    c.add_argument("--no-verify", action="store_true")
    c.add_argument("--out", required=True)
    c.set_defaults(func=cmd_compare)

    v = sub.add_parser("validate", help="check tier or run tables against the contract")
    v.add_argument("--tier", nargs="*")
    v.add_argument("--runs", nargs="*")
    v.add_argument("--no-verify", action="store_true")
    v.set_defaults(func=cmd_validate)

    ad = sub.add_parser("adapters", help="list installed adapters")
    ad.set_defaults(func=cmd_adapters)
    return p


def main(argv=None):
    a = build_parser().parse_args(argv)
    try:
        return a.func(a) or 0
    except Exception as e:              # noqa: BLE001
        if "TAXOBENCH_DEBUG" in __import__("os").environ:
            raise
        return _err(f"{type(e).__name__}: {e}")


if __name__ == "__main__":
    sys.exit(main())
