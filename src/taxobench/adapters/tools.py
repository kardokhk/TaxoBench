"""Built-in adapters. Each parses the tool's own raw output (not the TaxoBench 1
`taxonomy_final.tsv`, except `taxonomy-tsv`, which reads that format for QIIME 2
tools and for precomputed outputs).

Changes from the TaxoBench 1 parsers, by design:
- MAPseq: a rank below the cutoff ends the lineage; no lower rank survives it.
- Kraken 2: the per-taxon value the TaxoBench 1 parser called confidence is the clade's
  share of all reads in the run (report percentage / 100). It is kept, but marked
  ``run_share``: not a per-read value and not a probability. It never enters
  curves or calibration.
- SPINGO and Kraken 2 names are repaired only when the repair is unique in the
  run's reference; otherwise the call is `oov`.
"""
from __future__ import annotations

import csv
import re
from pathlib import Path

from .. import RANKS
from ..contract import ContractError
from ..taxonomy import Reference, is_nocall
from .base import Adapter, Parsed, register, split_prefixed, to_float


def _read_tsv(fp, header=True):
    with open(fp, newline="") as fh:
        rows = list(csv.reader(fh, delimiter="\t", quoting=csv.QUOTE_NONE))
    rows = [r for r in rows if r and any(x.strip() for x in r)]
    if header:
        return rows[0], rows[1:]
    return None, rows


@register
class TaxonovaAdapter(Adapter):
    """`taxonova classification|classify` raw table: sequence_id, then
    <rank>, <rank>_confidence for the seven ranks. "Unclassified"/"Unknown"
    at a rank is no call there and below."""
    name = "taxonova"
    version = "1"
    tool = "taxonova"
    confidence_kind = "probability"

    def parse(self, raw, **params):
        hdr, rows = _read_tsv(raw)
        col = {h: i for i, h in enumerate(hdr)}
        need = ["sequence_id"] + list(RANKS) + [f"{r}_confidence" for r in RANKS]
        miss = [c for c in need if c not in col]
        if miss:
            raise ContractError(f"taxonova raw {raw}: missing columns {miss}")
        out = []
        for row in rows:
            lin, conf = {}, {}
            for r in RANKS:
                nm = row[col[r]].strip()
                if nm.lower() in ("unclassified", "unknown", ""):
                    break
                lin[r] = nm
                conf[r] = to_float(row[col[f"{r}_confidence"]])
            out.append(Parsed(row[col["sequence_id"]], lin, conf,
                              "; ".join(row[col[r]] for r in RANKS)))
        return out


@register
class TaxonomyTsvAdapter(Adapter):
    """QIIME 2 style `Feature ID, Taxon, Confidence|Consensus` table: BLAST+,
    VSEARCH and naive Bayes through q2-feature-classifier, and any precomputed
    output in that format (including the TaxoBench 1 `taxonomy_final.tsv`).

    One confidence value applies to every named rank; a `; `-separated value
    list is read per rank."""
    name = "taxonomy-tsv"
    version = "1"
    tool = "precomputed"
    confidence_kind = "score"

    def __init__(self, confidence_kind=None):
        if confidence_kind:
            self.confidence_kind = confidence_kind

    def parse(self, raw, **params):
        hdr, rows = _read_tsv(raw)
        if hdr[0] not in ("Feature ID", "Feature.ID", "id", "#OTU ID", "query_id"):
            raise ContractError(f"{raw}: unexpected header {hdr[:3]}")
        out = []
        for row in rows:
            qid, label = row[0], (row[1] if len(row) > 1 else "")
            lin, order = split_prefixed(label)
            lin = {r: t for r, t in lin.items() if not is_nocall(t)}
            conf = {}
            if len(row) > 2 and row[2].strip():
                vals = [to_float(x) for x in row[2].split(";")]
                if len(vals) == 1:
                    conf = {r: vals[0] for r in lin}
                else:
                    conf = {r: v for r, v in zip(order, vals) if r in lin}
            out.append(Parsed(qid, lin, conf, label))
        return out


register(type("QiimeBlastAdapter", (TaxonomyTsvAdapter,), {"name": "blast", "tool": "blast"}))
register(type("QiimeVsearchAdapter", (TaxonomyTsvAdapter,), {"name": "vsearch", "tool": "vsearch"}))
register(type("QiimeNaiveBayesAdapter", (TaxonomyTsvAdapter,),
              {"name": "naive-bayes", "tool": "naive-bayes", "confidence_kind": "probability"}))


@register
class IdtaxaAdapter(Adapter):
    """`run_idtaxa_classification.R` output: Feature.ID, Taxon (starting with
    Root; `unclassified_*` for abstentions), Confidence (per rank)."""
    name = "idtaxa"
    version = "1"
    tool = "idtaxa"
    confidence_kind = "score"

    def parse(self, raw, **params):
        hdr, rows = _read_tsv(raw)
        out = []
        for row in rows:
            qid, tax, conf_s = row[0], row[1], (row[2] if len(row) > 2 else "")
            toks = [t.strip() for t in tax.split(";")]
            confs = [to_float(c) for c in conf_s.split(";")] if conf_s else []
            if toks and toks[0].endswith("Root"):
                toks, confs = toks[1:], confs[1:]
            lin, conf = {}, {}
            for i, r in enumerate(RANKS):
                if i >= len(toks):
                    break
                t = re.sub(r"^([a-z]__)+", "", toks[i]).strip()
                if is_nocall(t):
                    break
                lin[r] = t
                if i < len(confs):
                    c = confs[i]
                    conf[r] = c / 100.0 if c > 1.0 else c
            out.append(Parsed(qid, lin, conf, tax))
        return out


def _underscore_forms(name: str) -> set:
    return {name.replace("_", " ")}


@register
class SpingoAdapter(Adapter):
    """SPINGO 1.3 raw lines: id, score, genus, genus_conf, species, species_conf.

    SPINGO writes species names with a space where GTDB-style references have
    an underscore (`Enterococcus J nangangensis` for `Enterococcus_J
    nangangensis`). A name absent from the reference is replaced by the
    reference name equal to it after `_` -> space, only when exactly one such
    name exists. Ranks above genus come from the reference path of the genus.
    """
    name = "spingo"
    version = "1"
    tool = "spingo"
    confidence_kind = "score"
    complete_from_reference = True

    def parse(self, raw, **params):
        _, rows = _read_tsv(raw, header=False)
        out = []
        for row in rows:
            if len(row) != 6:
                raise ContractError(f"SPINGO {raw}: line with {len(row)} columns, expected 6: {row[:2]}")
            qid, _, g, gc, s, sc = row
            lin, conf = {}, {}
            if not is_nocall(g):
                lin["genus"], conf["genus"] = g.strip(), to_float(gc)
                if not is_nocall(s):
                    lin["species"], conf["species"] = s.strip(), to_float(sc)
            out.append(Parsed(qid, lin, conf, "\t".join(row[2:])))
        return out

    def repair(self, rank, name, ref: Reference):
        key = ("_spingo_norm", rank)
        idx = ref.cache.get(key)
        if idx is None:
            idx = {}
            for (rk, nm) in ref.index:
                if rk == rank:
                    idx.setdefault(nm.replace("_", " "), set()).add(nm)
            ref.cache[key] = idx
        cands = idx.get(name.replace("_", " "), set())
        return sorted(cands) if len(cands) == 1 else ()


_GTDB_SUFFIX = [re.compile(r"_\d+$"), re.compile(r"(_[A-Z]+)?(_\d+)?$")]


def _stripped(name: str) -> set:
    return {p.sub("", name) for p in _GTDB_SUFFIX} | {name}


@register
class Kraken2Adapter(Adapter):
    """Kraken 2 `--output` (with `--use-names`) plus `--report`.

    The lineage of a read comes from the report tree (ranks D, P, C, O, F, G, S;
    sub-ranks skipped). Kraken 2 databases built from GTDB-style taxonomies can
    drop the numeric suffix of a genus (`Dialister` for `Dialister_36337`); a
    name absent from the reference is matched to reference names with such
    suffixes removed, and accepted only when one taxon remains (after narrowing
    by the resolved family).

    Confidence: the report percentage of the assigned clade divided by 100,
    kept for continuity with TaxoBench 1 and marked run_share: it is the clade's share
    of all reads in the run, not a per-read value and not a probability.
    """
    name = "kraken2"
    version = "1"
    tool = "kraken2"
    confidence_kind = "run_share"
    confidence_note = "clade share of all reads in the run (report %/100); not a probability"
    RANK_CODE = {"D": "domain", "P": "phylum", "C": "class", "O": "order",
                 "F": "family", "G": "genus", "S": "species"}

    def parse(self, raw, **params):
        if isinstance(raw, dict):
            out_fp, rep_fp = Path(raw["output"]), Path(raw["report"])
        else:
            raw = Path(raw)
            out_fp, rep_fp = raw / "kraken2_classification.txt", raw / "kraken2_report.txt"
        tree = self._report(rep_fp)
        out = []
        with open(out_fp) as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 3:
                    continue
                status, qid, tx = parts[0], parts[1], parts[2]
                m = re.match(r"^(.*)\(taxid (\d+)\)\s*$", tx)
                taxid = int(m.group(2)) if m else int(tx) if tx.strip().isdigit() else 0
                lin, conf = {}, {}
                if status == "C" and taxid != 0:
                    if taxid not in tree:
                        raise ContractError(f"kraken2: taxid {taxid} of {qid} not in report {rep_fp}")
                    t = taxid
                    while t:
                        code, nm, share, par = tree[t]
                        r = self.RANK_CODE.get(code)
                        if r and r not in lin:
                            lin[r] = nm
                            conf[r] = share
                        t = par
                out.append(Parsed(qid, lin, conf, tx))
        return out

    @staticmethod
    def _report(fp):
        tree, stack = {}, []
        with open(fp) as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 6:
                    continue
                pct, code, taxid, nm = float(parts[0]), parts[3].strip(), int(parts[4]), parts[5]
                depth = len(nm) - len(nm.lstrip(" "))
                while stack and stack[-1][0] >= depth:
                    stack.pop()
                par = stack[-1][1] if stack else 0
                tree[taxid] = (code, nm.strip(), pct / 100.0, par)
                stack.append((depth, taxid))
        return tree

    def repair(self, rank, name, ref: Reference):
        key = ("_k2_norm", rank)
        idx = ref.cache.get(key)
        if idx is None:
            idx = {}
            for (rk, nm) in ref.index:
                if rk != rank:
                    continue
                if rank == "species" and " " in nm:
                    g, ep = nm.split(" ", 1)
                    forms = {f + " " + ep for f in _stripped(g)}
                else:
                    forms = _stripped(nm)
                for f in forms:
                    idx.setdefault(f, set()).add(nm)
            ref.cache[key] = idx
        return sorted(idx.get(name, set()))


@register
class MapseqAdapter(Adapter):
    """MAPseq 2 `.mseq` output. Columns 14 onwards hold, per rank, the taxon,
    combined_cf and score_cf. The confidence column (default score_cf;
    `confidence_column="combined_cf"` selects the other) is compared with
    `cutoff`: the first rank below it, or without a taxon, ends the lineage
    (TaxoBench 1 kept lower ranks after a dropped one)."""
    name = "mapseq"
    version = "1"
    tool = "mapseq"
    confidence_kind = "score"
    FIRST = 13

    def parse(self, raw, cutoff=0.0, confidence_column="score_cf", **params):
        off = {"combined_cf": 1, "score_cf": 2}[confidence_column]
        out = []
        with open(raw) as fh:
            for line in fh:
                if not line.strip() or line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                qid = parts[0]
                lin, conf = {}, {}
                for i, r in enumerate(RANKS):
                    b = self.FIRST + 3 * i
                    if b + 2 >= len(parts):
                        break
                    nm = parts[b].strip()
                    c = to_float(parts[b + off])
                    if is_nocall(nm) or not (c >= float(cutoff)):
                        break
                    lin[r], conf[r] = nm, c
                out.append(Parsed(qid, lin, conf, "\t".join(parts[self.FIRST:])))
        return out


@register
class SpingoFinalAdapter(TaxonomyTsvAdapter):
    """SPINGO calls already converted to the TaxoBench 1 `taxonomy_final.tsv`
    format (empty upper ranks), with the SPINGO name repair and reference
    completion."""
    name = "spingo-final"
    tool = "spingo"
    confidence_kind = "score"
    complete_from_reference = True
    repair = SpingoAdapter.repair


@register
class Kraken2FinalAdapter(TaxonomyTsvAdapter):
    """Kraken 2 calls in the TaxoBench 1 `taxonomy_final.tsv` format, with the unique
    suffix repair. Its confidence column is the run share (see Kraken2Adapter)."""
    name = "kraken2-final"
    tool = "kraken2"
    confidence_kind = "run_share"
    confidence_note = Kraken2Adapter.confidence_note
    repair = Kraken2Adapter.repair


@register
class MapseqCombinedAdapter(MapseqAdapter):
    """Labelled variant: the cutoff applies to combined_cf instead of score_cf."""
    name = "mapseq-combined"

    def parse(self, raw, cutoff=0.0, confidence_column="combined_cf", **params):
        return super().parse(raw, cutoff=cutoff, confidence_column=confidence_column, **params)
