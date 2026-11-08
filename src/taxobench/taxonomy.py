"""Reference vocabularies, lineage parsing and name resolution.

One group parser for the whole package: a token is a group
when it is wrapped in one pair of brackets that closes at its last character,
and its members are split on ``|`` at bracket depth 0. So
``[[Ruminococcus] gnavus group|Blautia]`` has two members and
``[Ruminococcus] gnavus group`` is one ordinary (SILVA-style) name.

Name resolution never guesses. A name resolves when exactly one taxon of
that rank and name exists in the run's reference, after narrowing by the
resolved parent when several share the name. Anything else is ``oov``. An
adapter may supply a repair function proposing alternative spellings; a repair
is accepted only when it, too, leads to exactly one taxon.
"""
from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field
from typing import Callable, Iterable

import pandas as pd

from . import PREFIX, RANKS
from .contract import ContractError

RANK_INDEX = {r: i for i, r in enumerate(RANKS)}
PREFIX_TO_RANK = {v: k for k, v in PREFIX.items()}

#: Tokens that mean "no call" at a rank, after the rank prefix is removed
#: (compared case-insensitively). Tool-specific forms are normalised by the
#: adapters; these are the shared ones.
NOCALL_TOKENS = {"", "unassigned", "unclassified", "unknown", "ambiguous",
                 "unidentified", "na", "n/a", "none", "no blast hit"}


#: Non-specific words (open nomenclature): a name containing one of them, or a
#: word starting with uncultured / unidentified / metagenome / endosymbiont,
#: may label distinct lineages and is never merged across paths. Source: the
#: placeholder rule of the species-F1 campaign (PREPROCESSING_SPEC.md rule 1;
#: preproc/lib/preproc.py OPEN_WORDS, 2026-10-31).
NON_SPECIFIC_WORDS = {"sp", "sp.", "bacterium", "uncultured", "unidentified", "metagenome",
                      "environmental", "unknown", "unclassified", "cf.", "cf", "aff.", "aff",
                      "str.", "n.", "group", "clone", "endosymbiont", "symbiont"}
NON_SPECIFIC_PREFIXES = ("uncultured", "unidentified", "metagenome", "endosymbiont")


def is_non_specific(name: str) -> bool:
    words = [w.lower() for w in name.replace("_", " ").split()]
    return any(w in NON_SPECIFIC_WORDS or w.startswith(NON_SPECIFIC_PREFIXES) for w in words)


def is_nocall(name: str) -> bool:
    n = name.strip().lower()
    return n in NOCALL_TOKENS or n.startswith("unclassified_") or n.startswith("unclassified ")


def split_group(token: str) -> list[str]:
    """Members of a group token, or [token] when it is not a group."""
    t = token.strip()
    if len(t) < 2 or t[0] != "[" or t[-1] != "]":
        return [t]
    depth = 0
    for i, ch in enumerate(t):
        if ch == "[":
            depth += 1
        elif ch == "]":
            depth -= 1
            if depth == 0 and i != len(t) - 1:
                return [t]            # e.g. "[Ruminococcus] gnavus group"
    inner = t[1:-1]
    members, depth, cur = [], 0, []
    for ch in inner:
        if ch == "[":
            depth += 1
        elif ch == "]":
            depth -= 1
        if ch == "|" and depth == 0:
            members.append("".join(cur).strip())
            cur = []
        else:
            cur.append(ch)
    members.append("".join(cur).strip())
    return [m for m in members if m]


def parse_lineage(label: str, strict_prefix=True) -> dict:
    """``d__A; p__B; ...`` -> {rank: token}. Unknown prefixes raise.

    Tokens without any rank prefix are an error when strict_prefix is set;
    adapters for tools that write positional lineages convert them first.
    """
    out = {}
    if label is None:
        return out
    label = label.strip()
    if not label or is_nocall(label):
        return out
    for part in label.split(";"):
        part = part.strip()
        if not part:
            continue
        pre = part[:3]
        if pre in PREFIX_TO_RANK:
            rank = PREFIX_TO_RANK[pre]
            if rank in out:
                raise ContractError(f"rank {rank} twice in lineage {label!r}")
            out[rank] = part[3:].strip()
        elif len(part) >= 3 and part[1:3] == "__":
            raise ContractError(f"unknown rank prefix {pre!r} in {label!r}")
        elif strict_prefix:
            raise ContractError(f"token without rank prefix {part!r} in {label!r}")
    return out


@dataclass
class Resolved:
    """Resolution of one call at one rank."""
    call: str                      # named | group | nocall | oov
    taxon_ids: tuple = ()
    note: str = ""


@dataclass
class Reference:
    ref_id: str
    rank: dict = field(default_factory=dict)        # taxon_id -> rank
    name: dict = field(default_factory=dict)        # taxon_id -> name
    parent: dict = field(default_factory=dict)      # taxon_id -> parent id or ""
    is_group: dict = field(default_factory=dict)    # taxon_id -> bool
    index: dict = field(default_factory=lambda: defaultdict(list))   # (rank, name) -> [ids]
    _atoms: dict = field(default_factory=dict)
    cache: dict = field(default_factory=dict)      # adapter repair indexes
    _last_candidates: tuple = ()
    homonym_groups: set = field(default_factory=set)
    in_reference: dict = field(default_factory=dict)   # taxon_id -> bool (default True)
    unexpanded_groups: set = field(default_factory=set)

    @classmethod
    def from_taxa(cls, taxa: pd.DataFrame, ref_id: str) -> "Reference":
        sub = taxa[taxa["ref_id"] == ref_id]
        if sub.empty:
            raise ContractError(f"no taxa rows for ref_id {ref_id!r}")
        ref = cls(ref_id=ref_id)
        grp = sub["is_group"] if "is_group" in sub.columns else pd.Series("0", index=sub.index)
        inr = sub["in_reference"] if "in_reference" in sub.columns else pd.Series("1", index=sub.index)
        for tid, rk, nm, par, g, ir in zip(sub["taxon_id"], sub["rank"], sub["name"],
                                           sub["parent_taxon_id"], grp, inr):
            ref.in_reference[tid] = str(ir) not in ("0", "false", "False")
            ref.rank[tid] = rk
            ref.name[tid] = nm
            ref.parent[tid] = par
            ref.is_group[tid] = str(g) in ("1", "true", "True")
            ref.index[(rk, nm)].append(tid)
        for tid, par in ref.parent.items():
            if par and par not in ref.rank:
                raise ContractError(f"taxa {ref_id}: parent {par!r} of {tid!r} is not a taxon")
        return ref

    @classmethod
    def from_lineages(cls, lineages: Iterable[str], ref_id: str) -> "Reference":
        """Build a reference from full lineage strings (for example a
        Greengenes2 taxonomy.tsv). taxon_id = rank-prefixed path."""
        rows = {}
        for lab in lineages:
            lin = parse_lineage(lab)
            path, parent = [], ""
            for r in RANKS:
                nm = lin.get(r, "")
                if is_nocall(nm):
                    break
                path.append(PREFIX[r] + nm)
                tid = ";".join(path).replace("|", "/")   # "|" separates ids in the tables
                if tid not in rows:
                    rows[tid] = (ref_id, r, tid, nm, parent,
                                 "1" if len(split_group(nm)) > 1 else "0")
                parent = tid
        df = pd.DataFrame(list(rows.values()),
                          columns=["ref_id", "rank", "taxon_id", "name", "parent_taxon_id", "is_group"])
        return cls.from_taxa(df, ref_id)

    # ------------------------------------------------------------------
    def ancestor_at(self, tid: str, rank: str):
        while tid:
            if self.rank[tid] == rank:
                return tid
            tid = self.parent[tid]
        return None

    def path(self, tid: str) -> list:
        out = []
        while tid:
            out.append(tid)
            tid = self.parent[tid]
        return out[::-1]

    def group_members(self, tid: str) -> tuple:
        """Members a taxon stands for, each as a frozenset of candidate taxa.

        A non-group taxon is one member, itself (a placeholder taxon, see
        below, may stand for its homonyms). A group taxon (is_group) is split
        into its bracketed members; each member resolves to the taxon of that
        rank and name, narrowed by the group's parent when several exist.

        Homonymous paths: when several candidates remain, or the member is a
        placeholder taxon (in_reference = 0, as a tier may write group members), the
        member is resolved by name within the rank and stands for every
        candidate (counted as n_homonym_resolved), unless the guard applies:
        the name is non-specific (NON_SPECIFIC_WORDS) or its candidates sit
        under more than one lineage at the rank above, with group labels on
        the parents expanded (Pantoea under Erwiniaceae and under
        [Erwiniaceae|X] is one lineage; endosymbionts under six families are
        distinct). A guarded member, or a member with no candidate, leaves the
        group unexpanded, standing for itself (counted as n_group_unexpanded).
        """
        if tid in self._atoms:
            return self._atoms[tid]
        res = (self._by_name_if_placeholder(tid),)
        if self.is_group.get(tid):
            members = split_group(self.name[tid])
            rank = self.rank[tid]
            sets, homonym = [], False
            for m in members:
                every = [x for x in self.index.get((rank, m), []) if not self.is_group[x]]
                c = every
                if len(c) > 1:
                    c2 = [x for x in c if self.parent[x] == self.parent[tid]]
                    if len(c2) == 1:
                        c = c2
                if c and (len(c) > 1 or not self.in_reference.get(c[0], True)) and len(every) > 1:
                    if not self.one_lineage(every, m):
                        sets = None
                        break
                    c, homonym = every, True
                if not c:
                    sets = None
                    break
                sets.append(frozenset(c))
            if sets is None:
                self.unexpanded_groups.add(tid)
            else:
                uniq = []
                for x in sets:
                    if x not in uniq:
                        uniq.append(x)
                res = tuple(uniq)
                if homonym:
                    self.homonym_groups.add(tid)
        self._atoms[tid] = res
        return res

    def one_lineage(self, cands, name: str) -> bool:
        """False when merging same-name candidates would join distinct
        lineages: a non-specific name, or parents (group labels expanded to
        their member names) that fall into more than one connected set."""
        if len(cands) <= 1:
            return True
        if is_non_specific(name):
            return False
        sets = [frozenset(split_group(self.name[self.parent[x]])) if self.parent[x] else frozenset([""])
                for x in cands]
        comps: list = []
        for s in sets:
            hit = [c for c in comps if c & s]
            merged = set(s).union(*hit) if hit else set(s)
            comps = [c for c in comps if not (c & s)] + [merged]
        return len(comps) == 1

    def _by_name_if_placeholder(self, tid: str) -> frozenset:
        """A non-group placeholder taxon (in_reference = 0) stands for every
        taxon with its name at its rank, unless the guard of group_members()
        applies (then it stands for itself and counts as unexpanded)."""
        if self.is_group.get(tid) or self.in_reference.get(tid, True):
            return frozenset([tid])
        wide = [x for x in self.index.get((self.rank[tid], self.name[tid]), []) if not self.is_group[x]]
        if len(wide) > 1:
            if not self.one_lineage(wide, self.name[tid]):
                self.unexpanded_groups.add(tid)
                return frozenset([tid])
            self.homonym_groups.add(tid)
            return frozenset(wide)
        return frozenset([tid])

    def members(self, ids: Iterable[str]) -> tuple:
        """Member sets of a call or truth (one or more taxon_ids), deduplicated."""
        out = []
        for t in ids:
            if t not in self.rank:
                raise ContractError(f"{self.ref_id}: unknown taxon_id {t!r}")
            for m in self.group_members(t):
                if m not in out:
                    out.append(m)
        return tuple(out)

    def atoms(self, tid: str) -> frozenset:
        """All taxa a taxon stands for (union of its member sets)."""
        return frozenset().union(*self.group_members(tid))

    def expand(self, ids: Iterable[str]) -> frozenset:
        out = set()
        for t in ids:
            if t not in self.rank:
                raise ContractError(f"{self.ref_id}: unknown taxon_id {t!r}")
            out |= self.atoms(t)
        return frozenset(out)

    # ------------------------------------------------------------------
    def _lookup(self, rank, name, parents, repair):
        cands = self.index.get((rank, name), [])
        if not cands and repair is not None:
            alts = set()
            for alt in repair(rank, name, self):
                alts.update(self.index.get((rank, alt), []))
            cands = sorted(alts)
        if len(cands) > 1 and parents:
            narrowed = [c for c in cands if self.parent[c] in parents]
            if narrowed:
                cands = narrowed
        if len(cands) > 1:
            atomic = [c for c in cands if not self.is_group[c]]
            if len(atomic) == 1:
                cands = atomic
        if len(cands) == 1:
            return cands[0], ""
        self._last_candidates = tuple(cands)
        return None, ("ambiguous" if cands else "absent")

    def resolve_lineage(self, lineage: dict,
                        repair: Callable | None = None) -> dict:
        """{rank: token} -> {rank: Resolved} for all seven ranks.

        A rank without a token, or with a no-call token, is nocall. A lower
        rank named below a nocall rank is resolved on its own (no parent
        narrowing); adapters enforce hierarchy where the tool requires it.
        """
        out = {}
        ambiguous = {}
        parents: set = set()
        for r in RANKS:
            tok = lineage.get(r, "")
            if is_nocall(tok):
                out[r] = Resolved("nocall")
                parents = set()
                continue
            members = split_group(tok)
            ids, notes = [], []
            if len(members) > 1:
                # a group label that is itself a taxon of the reference (a tier may
                # write a classifier's training groups to `taxa`) resolves to that taxon; its
                # members are expanded by group_members()
                g = [t for t in self.index.get((r, tok.strip()), []) if self.is_group[t]]
                if len(g) > 1 and parents:
                    g = [t for t in g if self.parent[t] in parents] or g
                if len(g) == 1:
                    out[r] = Resolved("group" if len(self.members(g)) > 1 else "named", (g[0],))
                    parents = {g[0], self.parent[g[0]]}
                    continue
            for m in members:
                self._last_candidates = ()
                tid, why = self._lookup(r, m, parents, repair)
                if tid is None:
                    notes.append(f"{m}: {why}")
                    if why == "ambiguous" and len(members) == 1:
                        ambiguous[r] = self._last_candidates
                else:
                    ids.append(tid)
            if notes or not ids:
                out[r] = Resolved("oov", (), "; ".join(notes))
                parents = set()
                continue
            ids = tuple(sorted(set(ids)))
            out[r] = Resolved("group" if len(self.members(ids)) > 1 else "named", ids)
            parents = set(ids) | {self.parent[t] for t in ids if self.is_group[t]}
        # An ambiguous name is settled by the call's own deeper rank when that
        # rank resolved to one taxon whose ancestor is one of the candidates
        # (for example a suffix-less genus above a fully named species).
        for r in reversed(RANKS):
            if r not in ambiguous:
                continue
            below = [out[q] for q in RANKS[RANKS.index(r) + 1:]
                     if out[q].call == "named" and len(out[q].taxon_ids) == 1]
            if below:
                anc = self.ancestor_at(below[0].taxon_ids[0], r)
                if anc in ambiguous[r]:
                    out[r] = Resolved("named", (anc,), "settled by the deeper rank")
        return out


def load_reference_lineages(fp, ref_id: str) -> Reference:
    """Reference from a two-column ``id<TAB>lineage`` file (header optional)."""
    labs = []
    with open(fp) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2 or parts[0] in ("Feature ID", "id", "#OTU ID"):
                continue
            labs.append(parts[1])
    return Reference.from_lineages(labs, ref_id)
