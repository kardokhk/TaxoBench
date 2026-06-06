#!/usr/bin/env python3
"""Step 01b — prepare the MBARC-26 mock community for the v_mocks panel.

Adds:
  - mbarc26  (Singer 2016 MBARC-26 mock; 26 species, 16S from RefSeq NR_)

ZymoBIOMICS D6300 was evaluated and dropped: chkMocks D6300 16S sub-region
sequences match gg2 sister-species labels (e.g., S. aureus → argenteus,
B. subtilis → inaquosorum, S. enterica → bongori, P. aeruginosa →
paraeruginosa) at higher identity than canonical-name gg2 references, so the
gg2-trained classifier panel systematically miscalls species. Not a pipeline
issue — a reference-DB nomenclature mismatch.

CAMI2 was evaluated and rejected — it is a shotgun-only challenge with no
amplicon-level truth and no curated 16S reference.

Outputs (taxobench convention, mirroring mockrobiota-16S):
  data/mock_community/input_data/<mock>_seqs.fasta
  data/mock_community/input_data/<mock>_taxa.txt   (Feature ID \t Taxon)

Then ready for 02_generate.py (V4/V3-V4/V5-V7) and
scripts/02c_generate_v1v9_direct.py (V1-V9).

Network reach to NCBI eutils is required.
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

# ---------------------------------------------------------------------------
TAXOBENCH = Path("/path/to/taxobench")
INPUT_DIR = TAXOBENCH / "data" / "mock_community" / "input_data"

# MBARC-26 strain set (Singer 2016, Sci. Data) — gg2-style taxonomies for the
# 23 bacteria + 3 archaea. Species labels are the canonical binomials Singer
# used; higher ranks follow Greengenes2 r2024.09.
MBARC26 = [
    # (genus_species, gg2-style 7-rank taxonomy, NCBI search hint for 16S)
    ("Clostridium perfringens",
     "d__Bacteria; p__Bacillota_A; c__Clostridia; o__Clostridiales; "
     "f__Clostridiaceae; g__Clostridium; s__Clostridium perfringens",
     "Clostridium perfringens"),
    ("Hungateiclostridium thermocellum",
     "d__Bacteria; p__Bacillota_A; c__Clostridia; o__Eubacteriales; "
     "f__Hungateiclostridiaceae; g__Hungateiclostridium; "
     "s__Hungateiclostridium thermocellum",
     "Clostridium thermocellum"),
    ("Coraliomargarita akajimensis",
     "d__Bacteria; p__Verrucomicrobiota; c__Opitutae; o__Puniceicoccales; "
     "f__Puniceicoccaceae; g__Coraliomargarita; s__Coraliomargarita akajimensis",
     "Coraliomargarita akajimensis"),
    ("Corynebacterium glutamicum",
     "d__Bacteria; p__Actinomycetota; c__Actinomycetes; o__Mycobacteriales; "
     "f__Corynebacteriaceae; g__Corynebacterium; s__Corynebacterium glutamicum",
     "Corynebacterium glutamicum"),
    ("Desulfosporosinus acidiphilus",
     "d__Bacteria; p__Bacillota_E; c__Desulfitobacteriia; o__Desulfitobacteriales; "
     "f__Desulfitobacteriaceae; g__Desulfosporosinus; "
     "s__Desulfosporosinus acidiphilus",
     "Desulfosporosinus acidiphilus"),
    ("Desulfosporosinus meridiei",
     "d__Bacteria; p__Bacillota_E; c__Desulfitobacteriia; o__Desulfitobacteriales; "
     "f__Desulfitobacteriaceae; g__Desulfosporosinus; "
     "s__Desulfosporosinus meridiei",
     "Desulfosporosinus meridiei"),
    ("Desulfotomaculum gibsoniae",
     "d__Bacteria; p__Bacillota_F; c__Desulfotomaculia; o__Desulfotomaculales; "
     "f__Desulfotomaculaceae; g__Desulfotomaculum; "
     "s__Desulfotomaculum gibsoniae",
     "Desulfotomaculum gibsoniae"),
    ("Echinicola vietnamensis",
     "d__Bacteria; p__Bacteroidota; c__Cytophagia; o__Cytophagales; "
     "f__Cyclobacteriaceae; g__Echinicola; s__Echinicola vietnamensis",
     "Echinicola vietnamensis"),
    ("Escherichia coli",
     "d__Bacteria; p__Pseudomonadota; c__Gammaproteobacteria; "
     "o__Enterobacterales; f__Enterobacteriaceae; g__Escherichia; "
     "s__Escherichia coli",
     "Escherichia coli"),
    ("Fervidobacterium pennivorans",
     "d__Bacteria; p__Thermotogota; c__Thermotogae; o__Thermotogales; "
     "f__Fervidobacteriaceae; g__Fervidobacterium; "
     "s__Fervidobacterium pennivorans",
     "Fervidobacterium pennivorans"),
    ("Frateuria aurantia",
     "d__Bacteria; p__Pseudomonadota; c__Gammaproteobacteria; o__Xanthomonadales; "
     "f__Rhodanobacteraceae; g__Frateuria; s__Frateuria aurantia",
     "Frateuria aurantia"),
    ("Halovivax ruber",
     "d__Archaea; p__Halobacteriota; c__Halobacteria; o__Natrialbales; "
     "f__Natrialbaceae; g__Halovivax; s__Halovivax ruber",
     "Halovivax ruber XH-70"),
    ("Hirschia baltica",
     "d__Bacteria; p__Pseudomonadota; c__Alphaproteobacteria; o__Hyphomonadales; "
     "f__Hyphomonadaceae; g__Hirschia; s__Hirschia baltica",
     "Hirschia baltica"),
    ("Meiothermus silvanus",
     "d__Bacteria; p__Deinococcota; c__Deinococci; o__Thermales; "
     "f__Thermaceae; g__Meiothermus; s__Meiothermus silvanus",
     "Meiothermus silvanus"),
    ("Natronobacterium gregoryi",
     "d__Archaea; p__Halobacteriota; c__Halobacteria; o__Natrialbales; "
     "f__Natrialbaceae; g__Natronobacterium; s__Natronobacterium gregoryi",
     "Natronobacterium gregoryi"),
    ("Natronococcus occultus",
     "d__Archaea; p__Halobacteriota; c__Halobacteria; o__Natrialbales; "
     "f__Natrialbaceae; g__Natronococcus; s__Natronococcus occultus",
     "Natronococcus occultus"),
    ("Nocardiopsis dassonvillei",
     "d__Bacteria; p__Actinomycetota; c__Actinomycetes; o__Streptosporangiales; "
     "f__Nocardiopsaceae; g__Nocardiopsis; s__Nocardiopsis dassonvillei",
     "Nocardiopsis dassonvillei"),
    ("Olsenella uli",
     "d__Bacteria; p__Actinomycetota; c__Coriobacteriia; o__Coriobacteriales; "
     "f__Atopobiaceae; g__Olsenella; s__Olsenella uli",
     "Olsenella uli"),
    ("Pseudomonas stutzeri",
     "d__Bacteria; p__Pseudomonadota; c__Gammaproteobacteria; o__Pseudomonadales; "
     "f__Pseudomonadaceae; g__Stutzerimonas; s__Stutzerimonas stutzeri",
     "Pseudomonas stutzeri"),
    ("Salmonella bongori",
     "d__Bacteria; p__Pseudomonadota; c__Gammaproteobacteria; o__Enterobacterales; "
     "f__Enterobacteriaceae; g__Salmonella; s__Salmonella bongori",
     "Salmonella bongori"),
    ("Salmonella enterica",
     "d__Bacteria; p__Pseudomonadota; c__Gammaproteobacteria; o__Enterobacterales; "
     "f__Enterobacteriaceae; g__Salmonella; s__Salmonella enterica",
     "Salmonella enterica subsp. arizonae"),
    ("Segniliparus rotundus",
     "d__Bacteria; p__Actinomycetota; c__Actinomycetes; o__Mycobacteriales; "
     "f__Segniliparaceae; g__Segniliparus; s__Segniliparus rotundus",
     "Segniliparus rotundus"),
    ("Sediminispirochaeta smaragdinae",
     "d__Bacteria; p__Spirochaetota; c__Spirochaetia; o__Spirochaetales; "
     "f__Sediminispirochaetaceae; g__Sediminispirochaeta; "
     "s__Sediminispirochaeta smaragdinae",
     "Spirochaeta smaragdinae"),
    ("Streptococcus pyogenes",
     "d__Bacteria; p__Bacillota; c__Bacilli; o__Lactobacillales; "
     "f__Streptococcaceae; g__Streptococcus; s__Streptococcus pyogenes",
     "Streptococcus pyogenes"),
    ("Terriglobus roseus",
     "d__Bacteria; p__Acidobacteriota; c__Terriglobia; o__Terriglobales; "
     "f__Acidobacteriaceae; g__Terriglobus; s__Terriglobus roseus",
     "Terriglobus roseus"),
    ("Thermobacillus composti",
     "d__Bacteria; p__Bacillota; c__Bacilli; o__Paenibacillales; "
     "f__Paenibacillaceae; g__Thermobacillus; s__Thermobacillus composti",
     "Thermobacillus composti"),
]


# ---------------------------------------------------------------------------
def parse_fasta(text: str) -> list[tuple[str, str]]:
    out: list[tuple[str, str]] = []
    cur, buf = None, []
    for line in text.splitlines():
        if line.startswith(">"):
            if cur is not None:
                out.append((cur, "".join(buf).upper()))
            cur = line[1:].strip()
            buf = []
        else:
            buf.append(line.strip())
    if cur is not None:
        out.append((cur, "".join(buf).upper()))
    return out


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        for h, s in records:
            f.write(f">{h}\n{s}\n")


def write_taxa(path: Path, records: list[tuple[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        f.write("Feature ID\tTaxon\n")
        for fid, tax in records:
            f.write(f"{fid}\t{tax}\n")


# ---------------------------------------------------------------------------
def fetch_one_16s(query_species: str, retries: int = 3) -> tuple[str, str] | None:
    """Search RefSeq NR_ records for a full-length 16S sequence for `query_species`.

    Returns (accession, sequence) of the first NR_ hit. Length 1300-1800 bp.
    """
    from Bio import Entrez  # lazy import — only when needed
    Entrez.email = "taxobench@projects.local"

    queries = [
        f'"{query_species}"[Organism] AND (16S ribosomal RNA[Title]) '
        f'AND 1300:1800[SLEN] AND refseq[filter]',
        f'"{query_species}"[Organism] AND 16S[Title] AND 1300:1800[SLEN]',
        f'"{query_species}"[All Fields] AND 16S ribosomal RNA[Title] '
        f'AND 1300:1800[SLEN]',
    ]
    for q in queries:
        for attempt in range(retries):
            try:
                h = Entrez.esearch(db="nucleotide", term=q, retmax=5)
                r = Entrez.read(h); h.close()
                ids = r.get("IdList", [])
                if not ids:
                    break
                # prefer NR_ accessions; fall back to first hit
                summary_h = Entrez.esummary(db="nucleotide", id=",".join(ids))
                summaries = Entrez.read(summary_h); summary_h.close()
                pick = None
                for s in summaries:
                    acc = s.get("AccessionVersion") or s.get("Caption", "")
                    if acc.startswith("NR_"):
                        pick = (acc, s["Id"]); break
                if pick is None:
                    pick = (summaries[0].get("AccessionVersion", ""),
                            summaries[0]["Id"])
                acc, uid = pick
                fetch_h = Entrez.efetch(db="nucleotide", id=uid,
                                        rettype="fasta", retmode="text")
                fa = fetch_h.read(); fetch_h.close()
                recs = parse_fasta(fa)
                if recs and len(recs[0][1]) >= 1200:
                    return acc, recs[0][1]
                break
            except Exception as e:
                print(f"      retry {attempt+1}/{retries} for {q!r}: {e}")
                time.sleep(2 ** attempt)
        # try next query
    return None


def prepare_mbarc26() -> None:
    print("\n=== MBARC-26 (Singer 2016) ===")
    print(f"  fetching 16S for {len(MBARC26)} species via NCBI eutils...")
    seq_recs: list[tuple[str, str]] = []
    tax_recs: list[tuple[str, str]] = []
    miss: list[str] = []
    for i, (label, tax, query) in enumerate(MBARC26, start=1):
        print(f"  [{i:2d}/26] {label} (query: {query!r})", flush=True)
        result = fetch_one_16s(query)
        if result is None:
            print(f"     ! NO HIT for {label}")
            miss.append(label)
            continue
        acc, seq = result
        fid = f"MBARC_{i:02d}_{acc}"
        seq_recs.append((fid, seq))
        tax_recs.append((fid, tax))
        print(f"     + {acc} ({len(seq)} bp)")
        time.sleep(0.4)  # be polite to NCBI
    print(f"\n  fetched {len(seq_recs)}/{len(MBARC26)} 16S sequences")
    if miss:
        print(f"  ! missing: {miss}")
        sys.exit(2)

    write_fasta(INPUT_DIR / "mbarc26_seqs.fasta", seq_recs)
    write_taxa(INPUT_DIR / "mbarc26_taxa.txt", tax_recs)
    print(f"  -> {INPUT_DIR / 'mbarc26_seqs.fasta'}")
    print(f"  -> {INPUT_DIR / 'mbarc26_taxa.txt'}")


def main():
    argparse.ArgumentParser(description=__doc__).parse_args()
    INPUT_DIR.mkdir(parents=True, exist_ok=True)
    prepare_mbarc26()
    print("\ndone.")


if __name__ == "__main__":
    main()
