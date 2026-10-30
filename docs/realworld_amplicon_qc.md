# Real-world amplicon preprocessing: paired-end audit (2026-05-04)

A dated audit of how four real-world studies (Mantri, Virtanen, Cao (coral) and Zuo)
were primer-trimmed, filtered and denoised with DADA2, and of an error found in the
paired-end handling of Zuo. The per-study preprocessing scripts belong to the archived
benchmark data and are not yet in this repository.

## Summary

| Study | Region | Read layout | Pipeline | Outcome |
|---|---|---|---|---|
| Mantri | V3-V4 | MiSeq 2x300 | cutadapt, DADA2 paired (`truncLen` 270/200), merge | correct; ASVs 402 to 427 bp |
| Virtanen | V3-V4 | MiSeq 2x300 | Trim Galore, DADA2 paired (270/260), merge | correct; ASVs 442 to 466 bp |
| Cao (coral) | V4 | iSeq 2x150 with primers | cutadapt (primers removed), DADA2 single-end (125) | kept; the pair cannot merge (below) |
| Zuo | V4 | MiSeq 2x150, primers already removed | cutadapt (quality only), DADA2 single-end (110) | error; rerun paired-end (below) |

## Reference values

| Region | Primers (length) | Insert length, primers excluded | Source |
|---|---|---|---|
| V4 | 515F (19 nt) / 806R (20 nt) | about 252 to 253 bp | Earth Microbiome Project; Illumina application note |
| V3-V4 | 341F (17 nt) / 785R (21 nt) | about 444 to 460 bp | Klindworth et al. 2013 |

DADA2 `mergePairs()` defaults: `minOverlap = 12`, `maxMismatch = 0`.

The V4 insert length excludes the primers. Once cutadapt removes the primers from a
150 bp read, the remaining 131 or 130 bp cover part of the 252 bp insert, not part of
252 + 39 bp. An earlier draft of this audit used the wrong convention, which gave an
optimistic overlap for the coral study. Both layouts below were checked on the reads.

## Per study

**Mantri.** Reads after cutadapt: R1 about 284 bp, R2 about 280 bp. With
`truncLen = c(270, 200)` the overlap is 270 + 200 - 410 = 60 bp. Observed ASVs: 402 to
427 bp.

**Virtanen.** Trim Galore (primers and quality), DADA2 paired with
`truncLen = c(270, 260)`, `mergePairs`. About 92% of filtered reads merged; 1,294 ASVs
of 442 to 466 bp.

**Cao (coral).** cutadapt removes the 515FY/806RB primers (19 and 20 nt), then DADA2
runs single-end with `truncLen = 125`. A paired-end merge is not possible with the
DADA2 defaults:

```
after cutadapt:  R1 about 131 bp, R2 about 130 bp
V4 insert:       about 252 bp
overlap:         131 + 130 - 252 = 9 bp, below the 12 bp minimum
```

A paired-end run with `truncLen = c(140, 140)` lost every read in `filterAndTrim`,
because most reads are shorter than 140 bp after primer removal; a lower `truncLen`
reduces the overlap further. Lowering `minOverlap` to 8 risks false merges and was not
adopted. The single-end ASVs (125 bp, the first half of V4) are what every classifier
in the panel receives.

**Zuo.** The original run used single-end DADA2 with `truncLen = 110`, on the
assumption that primer removal leaves about 111 bp per read. The submitted reads had
already been primer-trimmed, so each read is 150 bp of insert:

```
overlap:  150 + 150 - 252 = 48 bp, above the 12 bp minimum
```

The study was rerun with DADA2 paired-end, `truncLen = c(145, 140)` and the
`mergePairs` defaults, after the same quality-only cutadapt step:

| | Single-end (`truncLen` 110) | Paired-end |
|---|---|---|
| ASVs | 1,769 | 3,218 |
| ASV length | 110 bp (all) | 252 to 253 bp (2,633 of 3,218 at 253 bp) |
| Samples | 42 | 42 |

`rep_seqs_zuo.fasta` and the DADA2 tables were replaced by the paired-end output; the
sample keys of `asv_table_ids.csv` still match `paired_samples.csv` (all 42 samples
checked). The Mantri, Virtanen and coral outputs did not change.

## Why the two V4 studies differ

Both are V4 2x150. In the coral study the first 19 or 20 bases of each read are primer,
so the pair covers 261 bp of a 252 bp insert (9 bp overlap). In the Zuo study the reads
start in the insert, so the pair covers 300 bp (48 bp overlap).

## References

- DADA2 `mergePairs` defaults: <https://rdrr.io/bioc/dada2/man/mergePairs.html>;
  <https://github.com/benjjneb/dada2/issues/852>
- Earth Microbiome Project 16S protocol (515F/806R):
  <https://earthmicrobiome.ucsd.edu/protocols-and-standards/16s/>
- Illumina 16S metagenomic sequencing application note:
  <https://www.illumina.com/content/dam/illumina-marketing/documents/products/appnotes/appnote_miseq_16S.pdf>
- Klindworth A, et al. Evaluation of general 16S ribosomal RNA gene PCR primers for
  classical and next-generation sequencing-based diversity studies. Nucleic Acids Res.
  2013;41(1):e1. doi:10.1093/nar/gks808
