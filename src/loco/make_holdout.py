#!/usr/bin/env python3
"""LOCO holdout generator.

Given a per-region reference DB (sequences.fasta + taxonomy.tsv + merged.fasta),
produce:

  retained_ref/
      sequences.fasta         # ref sequences with held-out taxa removed
      taxonomy.tsv            # matching taxonomy
      merged.fasta            # taxonova-style "id<TAB>tax" header fasta
  novel_queries.fasta         # held-out sequences (true taxon is novel)
  novel_truth.tsv             # ground-truth taxonomy of novel queries
  novel_meta.tsv              # which rank each novel query is novel at (the
                              # FIRST rank at which its true taxon is absent)
  control_queries.fasta       # sample of retained sequences (true taxon stays)
  control_truth.tsv           # ground-truth taxonomy of controls
  holdout_manifest.json       # bookkeeping: per-rank held-out taxa lists,
                              # seed, sizes, etc.

Stratified holdout design — for each region, ONE holdout that produces queries
novel at every rank simultaneously (see docs/loco/README.md).
"""
import argparse
import json
import random
from collections import Counter, defaultdict
from pathlib import Path


LEVELS = ["domain", "phylum", "class", "order", "family", "genus", "species"]
LEVEL_TO_INT = {n: i for i, n in enumerate(LEVELS)}


def parse_taxonomy(tax_str):
    """Return list[7] of stripped tax strings, padded with '' if short."""
    parts = [p.strip() for p in tax_str.split(";")]
    while len(parts) < 7:
        parts.append("")
    return parts[:7]


def lineage_at(parts, level_int):
    """Lineage prefix up to (and incl) `level_int`, joined with '; '."""
    return "; ".join(parts[: level_int + 1])


def load_reference(taxonomy_tsv, sequences_fasta):
    """Yield (id, tax_str, sequence) tuples in fasta order."""
    tax_map = {}
    with open(taxonomy_tsv) as fh:
        for i, line in enumerate(fh):
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if i == 0 and parts[0] in ("Feature ID", "id"):
                continue
            if len(parts) >= 2:
                tax_map[parts[0]] = parts[1]
    cur_id, cur_seq = None, []
    with open(sequences_fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur_id is not None and cur_id in tax_map:
                    yield cur_id, tax_map[cur_id], "".join(cur_seq)
                # support either ">id" or ">id\ttaxonomy"
                cur_id = line[1:].rstrip("\n").split("\t", 1)[0].split()[0]
                cur_seq = []
            else:
                cur_seq.append(line.rstrip("\n"))
        if cur_id is not None and cur_id in tax_map:
            yield cur_id, tax_map[cur_id], "".join(cur_seq)


def select_holdout_taxa(records, rank_targets, min_seqs, seed, fixed_ranks=None):
    """records: list[(id, tax_parts, seq)]. Returns dict[level_int] -> set of
    holdout lineage strings at that level. Stratified across the taxonomy.

    fixed_ranks: optional dict[rank_name] -> iterable of lineage strings. Those
    ranks are NOT redrawn; the given lineages are used verbatim. Used to hold
    the higher-rank strata constant across replicates so that only the species
    draw varies. Without it, redrawing all five strata moves the retained
    reference by -41 % to +20 % between seeds, which swamps the species-level
    sampling variability the replicates are meant to measure.
    """
    rng = random.Random(seed)
    fixed_ranks = fixed_ranks or {}
    # Build counts of sequences per lineage at each rank.
    counts_by_rank = defaultdict(Counter)
    for _, parts, _ in records:
        for level_name, lvl in LEVEL_TO_INT.items():
            counts_by_rank[lvl][lineage_at(parts, lvl)] += 1

    holdout_lineages = {}
    already_held = set()  # union of (lvl, lineage) so we don't double-hold
    # IMPORTANT: process from species (deepest) UP to class (shallowest) so a
    # genus held-out doesn't accidentally also be a held-out species lineage.
    for rank_name in ["species", "genus", "family", "order", "class"]:
        lvl = LEVEL_TO_INT[rank_name]
        target = rank_targets[rank_name]
        if rank_name in fixed_ranks:
            # Carried over from a previous replicate. Verify every lineage is
            # still present in this reference rather than silently holding out
            # nothing, which would inflate the retained reference.
            chosen = set(fixed_ranks[rank_name])
            present = set(counts_by_rank[lvl])
            absent = sorted(chosen - present)
            if absent:
                raise SystemExit(
                    "[loco/01] FATAL: %d fixed %s lineages are absent from this "
                    "reference, e.g. %s" % (len(absent), rank_name, absent[:3]))
            for lin in chosen:
                already_held.add((lvl, lin))
            holdout_lineages[lvl] = chosen
            continue
        # Eligible: has >= min_seqs sequences AND its parent rank has at least
        # 2 child taxa (else removing it strips the parent entirely).
        parent_lvl = lvl - 1
        children_per_parent = defaultdict(set)
        for lin, _ in counts_by_rank[lvl].items():
            parent_lin = "; ".join(lin.split("; ")[: parent_lvl + 1])
            children_per_parent[parent_lin].add(lin)
        eligible = []
        for lin, n in counts_by_rank[lvl].items():
            if n < min_seqs:
                continue
            # skip empty taxa labels (e.g. "g__")
            tail = lin.split("; ")[-1]
            if tail.endswith("__") or tail == "":
                continue
            parent_lin = "; ".join(lin.split("; ")[: parent_lvl + 1])
            if len(children_per_parent[parent_lin]) < 2:
                continue
            # avoid descendants of already-held parent ranks
            ancestor_lineages = ["; ".join(lin.split("; ")[: l + 1]) for l in range(lvl)]
            if any((al, l) in already_held for l, al in enumerate(ancestor_lineages)):
                continue
            eligible.append(lin)
        rng.shuffle(eligible)
        chosen = set(eligible[:target])
        for lin in chosen:
            already_held.add((lvl, lin))
        holdout_lineages[lvl] = chosen
    return holdout_lineages


def first_novel_rank(parts, holdout_lineages):
    """Return the deepest level at which `parts`'s lineage is in holdout."""
    for lvl in sorted(holdout_lineages.keys()):
        if lineage_at(parts, lvl) in holdout_lineages[lvl]:
            return lvl
    return None


def write_fasta(path, records, with_tax_header=False):
    with open(path, "w") as fh:
        for rec in records:
            if with_tax_header:
                _id, tax_str, seq = rec
                fh.write(f">{_id}\t{tax_str}\n{seq}\n")
            else:
                _id, _tax, seq = rec
                fh.write(f">{_id}\n{seq}\n")


def write_taxonomy_tsv(path, records, header=True):
    with open(path, "w") as fh:
        if header:
            fh.write("Feature ID\tTaxon\n")
        for _id, tax_str, _seq in records:
            fh.write(f"{_id}\t{tax_str}\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--region", required=True,
                    choices=["V4", "V3-V4", "V5-V7", "V1-V9", "RRN"])
    ap.add_argument("--ref-dir", required=True,
                    help="Source per-region ref dir under data/ref_dbs/<region>/")
    ap.add_argument("--out-dir", required=True,
                    help="Output dir, e.g. data/loco/<region>/")
    ap.add_argument("--seed", type=int, default=2026)
    ap.add_argument("--n-species",  type=int, default=200)
    ap.add_argument("--n-genera",   type=int, default=100)
    ap.add_argument("--n-families", type=int, default=50)
    ap.add_argument("--n-orders",   type=int, default=25)
    ap.add_argument("--n-classes",  type=int, default=10)
    ap.add_argument("--min-seqs-per-taxon", type=int, default=4)
    ap.add_argument("--fixed-higher-ranks", default=None, metavar="MANIFEST",
                    help="Path to a previous replicate's holdout_manifest.json. "
                         "Its genus, family, order and class hold-outs are reused "
                         "verbatim and only the species stratum is redrawn, so that "
                         "replicates differ in the species draw alone and the "
                         "retained reference stays comparable.")
    ap.add_argument("--n-controls", type=int, default=2000,
                    help="Random sample of retained sequences used as MC controls")
    args = ap.parse_args()

    ref_dir = Path(args.ref_dir)
    if args.region == "RRN":
        tax_tsv = ref_dir / "grond_clean.tsv"
        seq_fa  = ref_dir / "grond_clean.fasta"
    else:
        tax_tsv = ref_dir / "taxonomy.tsv"
        seq_fa  = ref_dir / "sequences.fasta"
    assert tax_tsv.exists() and seq_fa.exists(), f"missing {tax_tsv} or {seq_fa}"

    print(f"[loco] loading {tax_tsv} + {seq_fa}")
    records = []
    for _id, tax_str, seq in load_reference(tax_tsv, seq_fa):
        if not seq:
            continue
        parts = parse_taxonomy(tax_str)
        records.append((_id, parts, tax_str, seq))
    print(f"[loco] loaded {len(records):,} sequences")

    rank_targets = {
        "species":  args.n_species,
        "genus":    args.n_genera,
        "family":   args.n_families,
        "order":    args.n_orders,
        "class":    args.n_classes,
    }
    fixed_ranks = None
    fixed_src = None
    if args.fixed_higher_ranks:
        with open(args.fixed_higher_ranks) as fh:
            prev = json.load(fh)
        fixed_ranks = {r: prev["holdout_taxa"][r]
                       for r in ("genus", "family", "order", "class")}
        fixed_src = args.fixed_higher_ranks
        print("[loco/01] higher ranks fixed from %s: %s"
              % (fixed_src, {r: len(v) for r, v in fixed_ranks.items()}))

    holdout = select_holdout_taxa(
        [(r[0], r[1], r[3]) for r in records],
        rank_targets, args.min_seqs_per_taxon, args.seed,
        fixed_ranks=fixed_ranks,
    )
    print("[loco] holdout taxa per rank:")
    for lvl in sorted(holdout):
        print(f"        {LEVELS[lvl]:>8}: {len(holdout[lvl])} taxa")

    # Split sequences into novel / retained
    novel = []
    retained = []
    for _id, parts, tax_str, seq in records:
        novel_lvl = first_novel_rank(parts, holdout)
        if novel_lvl is not None:
            novel.append((_id, parts, tax_str, seq, novel_lvl))
        else:
            retained.append((_id, parts, tax_str, seq))
    print(f"[loco] novel={len(novel):,}  retained={len(retained):,}")

    # MC controls: random sample of retained sequences
    rng = random.Random(args.seed + 1)
    n_controls = min(args.n_controls, len(retained))
    controls = rng.sample(retained, n_controls)
    control_ids = {r[0] for r in controls}
    # Retained reference excludes the controls (so when we classify them, they
    # are not trivially matched to themselves).
    retained_ref = [r for r in retained if r[0] not in control_ids]
    print(f"[loco] retained_ref={len(retained_ref):,}  controls={len(controls):,}")

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    ref_out = out / "retained_ref"
    ref_out.mkdir(exist_ok=True)

    # Reference fasta + taxonomy
    write_fasta(ref_out / "sequences.fasta",
                [(r[0], r[2], r[3]) for r in retained_ref], with_tax_header=False)
    write_taxonomy_tsv(ref_out / "taxonomy.tsv",
                [(r[0], r[2], r[3]) for r in retained_ref], header=True)
    write_fasta(ref_out / "merged.fasta",
                [(r[0], r[2], r[3]) for r in retained_ref], with_tax_header=True)

    # Novel queries
    write_fasta(out / "novel_queries.fasta",
                [(r[0], r[2], r[3]) for r in novel], with_tax_header=False)
    write_taxonomy_tsv(out / "novel_truth.tsv",
                [(r[0], r[2], r[3]) for r in novel], header=True)
    with open(out / "novel_meta.tsv", "w") as fh:
        fh.write("Feature ID\tnovel_rank\tnovel_rank_int\ttrue_lineage\n")
        for _id, parts, tax_str, _seq, lvl in novel:
            fh.write(f"{_id}\t{LEVELS[lvl]}\t{lvl}\t{tax_str}\n")

    # Controls
    write_fasta(out / "control_queries.fasta",
                [(r[0], r[2], r[3]) for r in controls], with_tax_header=False)
    write_taxonomy_tsv(out / "control_truth.tsv",
                [(r[0], r[2], r[3]) for r in controls], header=True)

    # Manifest
    manifest = {
        "region": args.region,
        "seed": args.seed,
        "fixed_higher_ranks_from": fixed_src,
        "rank_targets": rank_targets,
        "min_seqs_per_taxon": args.min_seqs_per_taxon,
        "n_total_records": len(records),
        "n_retained_ref": len(retained_ref),
        "n_novel_queries": len(novel),
        "n_controls": len(controls),
        "holdout_taxa_count": {LEVELS[l]: len(taxa) for l, taxa in holdout.items()},
        "holdout_taxa": {LEVELS[l]: sorted(taxa) for l, taxa in holdout.items()},
        "novel_queries_per_rank": dict(Counter(LEVELS[r[4]] for r in novel)),
    }
    with open(out / "holdout_manifest.json", "w") as fh:
        json.dump(manifest, fh, indent=2)
    print(f"[loco] wrote manifest to {out/'holdout_manifest.json'}")
    print(f"[loco] novel queries by first-novel rank: {manifest['novel_queries_per_rank']}")


if __name__ == "__main__":
    main()
