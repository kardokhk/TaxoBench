"""Weighted (normalized) UniFrac on the GG2 2024.09 backbone phylogeny.

Designed to slot into `scripts/05_realworld_evaluate.py` and emit a new
`weighted_unifrac` column alongside jaccard / bray_curtis at species rank.

Implementation:
  - One-time `bp.parse_newick` of the GG2 phylogeny (~30 s, vs ~280 s
    for skbio.TreeNode.read). Cached in module state.
  - Per (study × classifier × param) build a `biom.Table` (tip × sample)
    where each species's per-sample mass is distributed uniformly across
    the GG2 backbone tips annotated with that species. Tips not present
    in the BP tree are dropped (libssu requires full tree-table
    representability).
  - Call `unifrac.ssu_inmem(table, tree, "weighted_normalized_fp32")`
    once per evaluation and extract per-pair distances from the
    returned DistanceMatrix.

Bracket predictions `s__[a|b|c]` are expanded to their members and the
ASV's mass is split equally across the union of bracket-member tips.

Species names not present in GG2 taxonomy contribute zero mass; this
matches BC's behaviour (those keys never overlap with truth either way).
"""
from __future__ import annotations
import json
import re
import sys
import time
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np

_GROUP_RE = re.compile(r"([dpcofgs]__)\[([^\]]+)\]")
_SPECIES_PREFIX = "s__"

def _reference_dir() -> Path:
    """reference_database_dir from configs/pipeline_config.json, honouring
    TAXOBENCH_REFERENCE_DATABASE_DIR."""
    import os
    env = os.environ.get("TAXOBENCH_REFERENCE_DATABASE_DIR")
    if env:
        return Path(env)
    root = Path(__file__).resolve().parents[1]
    try:
        cfg = json.loads((root / "configs" / "pipeline_config.json").read_text())
        ref = Path(cfg.get("reference_database_dir", "data/ref_dbs"))
    except (OSError, ValueError):
        ref = Path("data/ref_dbs")
    return ref if ref.is_absolute() else root / ref


_DEFAULT_GG2_DIR = _reference_dir() / "Greengenes2_backbone"
_DEFAULT_TAXONOMY = _DEFAULT_GG2_DIR / "2024.09_taxonomy.tsv"
_DEFAULT_TREE = _DEFAULT_GG2_DIR / "2024.09.phylogeny.id.nwk"
_CACHE_DIR = _DEFAULT_GG2_DIR / ".unifrac_cache"


# ----- species → tip mapping ------------------------------------------------ #


def _species_key(tax_string: str) -> str:
    if not tax_string:
        return ""
    for token in tax_string.split(";"):
        token = token.strip()
        if token.startswith(_SPECIES_PREFIX):
            return token
    return ""


def species_to_tips(taxonomy_tsv: Path = _DEFAULT_TAXONOMY) -> Dict[str, List[str]]:
    """{species_token: [backbone_id, ...]} from GG2 backbone taxonomy."""
    _CACHE_DIR.mkdir(exist_ok=True)
    cache_fp = _CACHE_DIR / "species_to_tips.json"
    if cache_fp.exists():
        try:
            with cache_fp.open() as f:
                return json.load(f)
        except Exception:
            pass

    out: Dict[str, List[str]] = defaultdict(list)
    with open(taxonomy_tsv) as f:
        header = f.readline()
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            try:
                fid, tax = line.split("\t", 1)
            except ValueError:
                continue
            sp = _species_key(tax)
            if not sp or sp == _SPECIES_PREFIX:
                continue
            out[sp].append(fid)
    out_d = dict(out)
    try:
        with cache_fp.open("w") as f:
            json.dump(out_d, f)
    except Exception:
        pass
    return out_d


# ----- tree (BP) cache ------------------------------------------------------ #

_TREE: Dict[str, object] = {"bp": None, "tip_set": None, "path": None}


def _load_bp_tree(tree_fp: Path = _DEFAULT_TREE):
    """Return (bp.BP tree, set_of_tip_names). Cached once per process."""
    if _TREE["bp"] is not None and _TREE["path"] == str(tree_fp):
        return _TREE["bp"], _TREE["tip_set"]
    import bp
    t0 = time.time()
    with open(tree_fp) as f:
        nwk = f.read()
    tree = bp.parse_newick(nwk)
    # Enumerate all tip names: BP.B[i] == 1 indicates an "open" (we want tips
    # which are positions with B[i]==1 and B[i+1]==0 in BP encoding). Use
    # tree.name(i) to fetch labels; iterate all positions and keep nonempty,
    # non-internal names that are leaves.
    tip_set = set()
    n = len(tree.B)
    for i in range(n):
        if tree.isleaf(i):
            name = tree.name(i)
            if name:
                tip_set.add(name)
    sys.stderr.write(
        f"[_unifrac_metric] BP tree loaded in {time.time()-t0:.1f}s; "
        f"tips={len(tip_set):,}\n"
    )
    _TREE["bp"] = tree
    _TREE["tip_set"] = tip_set
    _TREE["path"] = str(tree_fp)
    return tree, tip_set


# ----- bucket expansion + matrix build -------------------------------------- #


def _expand_species_bucket(bucket: str) -> List[str]:
    sp_token = _species_key(bucket)
    if not sp_token:
        return []
    m = _GROUP_RE.match(sp_token)
    if m:
        return [_SPECIES_PREFIX + x for x in m.group(2).split("|") if x]
    return [sp_token]


# ----- bracket scoring: best-of-N nearest member ---------------------------- #
#
# A bracket `s__[a|b|c]` whose truth-set intersection is empty (per
# _resolve_group_keys) should be scored as ONE wrong prediction at the
# bracket member's tip-set that's closest to truth. This mirrors Jaccard's
# bracket semantics (any matching member → TP) extended to UniFrac: any
# bracket-member-near-truth → low UniFrac penalty. The bracket is never
# more-penalised than its closest-to-truth member alone.
#
# "Closest" uses taxonomic-rank distance (cheap proxy for tree-distance;
# GG2's molecular tree is taxonomy-aligned by construction).


def _full_tax_split(tax_str: str) -> List[str]:
    """Return per-rank tokens, stripped of whitespace."""
    return [t.strip() for t in (tax_str or "").split(";")]


def _expand_bracket_full_tax(tax_str: str) -> List[str]:
    """Given a full-tax bucket like 'd__A; ...; s__[a|b|c]', return list of
    single-species full-tax strings: ['d__A; ...; s__a', ...].
    Non-bracket strings return [tax_str]."""
    parts = _full_tax_split(tax_str)
    if not parts:
        return [tax_str]
    last = parts[-1]
    m = _GROUP_RE.match(last)
    if not m:
        return [tax_str]
    prefix = m.group(1)
    members = [x for x in m.group(2).split("|") if x]
    out = []
    for memb in members:
        new_last = prefix + memb
        out.append("; ".join(parts[:-1] + [new_last]))
    return out


def _tax_rank_distance(tax_a: str, tax_b: str) -> int:
    """Number of divergent ranks between two full-tax strings (lower=closer).
    0 = identical species; 1 = same genus diff species; ... ; 7 = no shared rank."""
    a = _full_tax_split(tax_a)
    b = _full_tax_split(tax_b)
    n = min(len(a), len(b))
    if n == 0:
        return 99
    # Walk from species-rank (last) up; first match defines the divergence.
    for i in range(n - 1, -1, -1):
        if a[i] and a[i] == b[i]:
            return n - 1 - i
    return n


def _select_best_bracket_member(
    bucket_full_tax: str,
    truth_freqs: Dict[str, float],
) -> str:
    """Pick the bracket-member full-tax string with the smallest
    taxonomic-rank distance to any truth species. Ties broken by the
    truth species with the largest abundance (mass-weighted preference),
    then lexically. If no truth species, return the first member.
    """
    members = _expand_bracket_full_tax(bucket_full_tax)
    if len(members) <= 1 or not truth_freqs:
        return members[0] if members else bucket_full_tax
    best = None
    best_score = None  # (dist, -truth_mass, lex_index)
    for mi, m in enumerate(members):
        for t, tm in truth_freqs.items():
            d = _tax_rank_distance(m, t)
            sc = (d, -float(tm), mi)
            if best_score is None or sc < best_score:
                best_score = sc
                best = m
    return best or members[0]


def _resolve_unmatched_brackets(
    amp_freqs: Dict[str, Dict[str, float]],
    wgs_freqs: Dict[str, Dict[str, float]],
    pairing: Dict[str, str],
) -> Dict[str, Dict[str, float]]:
    """Substitute every unmatched-bracket bucket key with the best
    bracket member (nearest to paired-truth by tax-rank distance).
    Matched brackets are already resolved upstream by _resolve_group_keys."""
    out: Dict[str, Dict[str, float]] = {}
    for sid, fd in amp_freqs.items():
        truth_for_sid = wgs_freqs.get(pairing.get(sid, ""), {})
        new_fd: Dict[str, float] = {}
        for bucket, freq in fd.items():
            if "[" in bucket:
                chosen = _select_best_bracket_member(bucket, truth_for_sid)
                new_fd[chosen] = new_fd.get(chosen, 0.0) + freq
            else:
                new_fd[bucket] = new_fd.get(bucket, 0.0) + freq
        out[sid] = new_fd
    return out


def _build_biom(
    sample_freqs: Dict[str, Dict[str, float]],
    sp2tips: Dict[str, List[str]],
    tip_set: set,
):
    """Construct a biom.Table (tips × samples) for the given samples.

    Tips not present in the phylogeny's tip_set are dropped. Species
    without any in-tree tip contribute zero mass (consistent with BC).
    """
    import biom
    from scipy.sparse import lil_matrix

    sample_ids = list(sample_freqs.keys())
    used_tips: set = set()
    expanded: Dict[str, List[Tuple[List[str], float]]] = {}
    for sid in sample_ids:
        per_sample = []
        for bucket, freq in sample_freqs[sid].items():
            members = _expand_species_bucket(bucket)
            tips_union = []
            for sp in members:
                tips_union.extend(t for t in sp2tips.get(sp, []) if t in tip_set)
            if not tips_union:
                continue
            used_tips.update(tips_union)
            per_sample.append((tips_union, freq))
        expanded[sid] = per_sample

    if not used_tips:
        return None, sample_ids, []

    tip_ids = sorted(used_tips)
    tip_index = {t: i for i, t in enumerate(tip_ids)}
    mat = lil_matrix((len(tip_ids), len(sample_ids)), dtype=np.float64)
    for si, sid in enumerate(sample_ids):
        for tips_union, freq in expanded[sid]:
            share = freq / len(tips_union)
            for t in tips_union:
                mat[tip_index[t], si] += share
    table = biom.Table(mat.tocsr(), observation_ids=tip_ids, sample_ids=sample_ids)
    return table, sample_ids, tip_ids


# ----- public batched UniFrac ----------------------------------------------- #


def weighted_unifrac_pairs(
    amp_freqs: Dict[str, Dict[str, float]],
    wgs_freqs: Dict[str, Dict[str, float]],
    pairing: Dict[str, str],
    sp2tips: Dict[str, List[str]] | None = None,
    tree_fp: Path = _DEFAULT_TREE,
    bypass_tips: bool = True,
) -> Dict[Tuple[str, str], float]:
    """Batched weighted-normalized UniFrac for a study's amp/wgs pairs.

    `bypass_tips=True` enables libssu's tip-bypass approximation (~2× speed,
    near-zero accuracy loss; standard for downstream Qiime2 use).

    Returns {(amp_id, wgs_id): distance}. NaN if a pair has no in-tree mass.
    """
    if sp2tips is None:
        sp2tips = species_to_tips()

    import unifrac
    tree, tip_set = _load_bp_tree(tree_fp)

    amp_used = {sid: amp_freqs[sid] for sid in pairing if sid in amp_freqs}
    # Best-of-N bracket resolution against paired-truth (see module docstring).
    amp_resolved = _resolve_unmatched_brackets(amp_used, wgs_freqs, pairing)
    wgs_used = {pairing[sid]: wgs_freqs.get(pairing[sid], {})
                for sid in amp_used}

    all_freqs: Dict[str, Dict[str, float]] = {}
    for sid, fd in amp_resolved.items():
        all_freqs[f"AMP::{sid}"] = fd
    for wsid, fd in wgs_used.items():
        all_freqs[f"WGS::{wsid}"] = fd

    table, sample_ids, tip_ids = _build_biom(all_freqs, sp2tips, tip_set)
    if table is None:
        return {(a, pairing[a]): float("nan") for a in amp_used}

    dm = unifrac.ssu_inmem(
        table, tree, "weighted_normalized_fp32",
        variance_adjust=False, alpha=1.0,
        bypass_tips=bypass_tips, n_substeps=1,
    )
    dm_df = dm.to_data_frame()
    out: Dict[Tuple[str, str], float] = {}
    for amp_id, wgs_id in pairing.items():
        a = f"AMP::{amp_id}"
        w = f"WGS::{wgs_id}"
        if a in dm_df.index and w in dm_df.columns:
            out[(amp_id, wgs_id)] = float(dm_df.loc[a, w])
        else:
            out[(amp_id, wgs_id)] = float("nan")
    return out
