# Real-world amplicon trimming/denoising QC: paired-end audit (2026-05-04)

A dated audit of the real-world preprocessing. The batch scripts it names
(`slurm/*.sbatch`) belong to the benchmark data tree and are not part of this
repository.

This note audits how each of the four real-world studies (virtanen, mantri,
coral, zuo) handles primer trimming, quality filtering, and DADA2 denoising,
and documents the bug found in the paired-end handling.

---

## TL;DR

| Study    | Region | Read layout                   | Pipeline                                         | Verdict                                |
|----------|--------|-------------------------------|--------------------------------------------------|----------------------------------------|
| mantri   | V3-V4  | MiSeq 2×300                   | cutadapt → DADA2 paired (270/200) → merge        | yes correct (ASVs ~402–427 bp)          |
| virtanen | V3-V4  | MiSeq 2×300                   | trim_galore → DADA2 paired (270/260) → merge     | yes correct (ASVs 442–466 bp)           |
| coral    | V4     | iSeq **2×150 with primers**   | cutadapt (primers stripped) → DADA2 single-end (125) | yes kept (paired-end is mathematically infeasible, see §3) |
| **zuo**  | V4     | MiSeq 2×150 (**pre-trimmed**) | cutadapt-Q only → DADA2 **single-end** (110)     | no paired-end SHOULD merge, fixed     |

---

## 1. Reference numbers used below

| Region | Primers (length)            | Insert length          | Source |
|--------|-----------------------------|------------------------|--------|
| V4     | 515F (19 nt) / 806R (20 nt) | ~252–253 bp (already excludes primers) | Earth Microbiome Project; Illumina app note |
| V3-V4  | 341F (17 nt) / 785R (21 nt) | ~444–460 bp            | Klindworth et al. 2013 |

DADA2 `mergePairs()` defaults: `minOverlap = 12`, `maxMismatch = 0`.

> **Important convention for V4.** "V4 is ~252 bp" refers to the length of
> the amplified insert *between* the primer-binding sites; the primers
> themselves are not included. So once cutadapt removes the primers from a
> 150 bp read, the remaining 131 bp / 130 bp covers a portion of the *252
> bp* insert (not a portion of "252 + 39"). An earlier draft of this audit
> used the wrong convention, which led to over-optimistic merge math for
> coral. Both layouts below have been re-verified empirically.

---

## 2. Per-study findings

### 2.1 mantri yes correct (kept as reference)

Sampled lengths after cutadapt: R1 ≈ 284 bp, R2 ≈ 280 bp.
DADA2 `truncLen = c(270, 200)`. Merge math: 270 + 200 − ~410 = **60 bp** ≥ 12 OK.
Observed ASV peak 402–427 bp.

### 2.2 virtanen yes correct (kept as reference)

Pipeline at `data/real_world_input/virtanen/dada2_output/dada2_analysis.R`:
trim_galore (primer + quality) → DADA2 paired `truncLen = c(270, 260)` →
mergePairs. ~92 % of filtered reads merge. 1 294 ASVs at 442–466 bp.

### 2.3 coral: paired-end infeasible with default DADA2 yes keep single-end

`slurm/coral_amplicon_pipeline.sbatch` runs cutadapt to strip the
515FY/806RB primers (19+20 nt) and then DADA2 single-end with
`truncLen = 125`.

The original justification ("9 bp overlap < 12") is correct under the
right insert convention:

```
post-cutadapt:  R1 ≈ 131 bp, R2 ≈ 130 bp
V4 insert:      ~252 bp
overlap:        131 + 130 − 252  =  9 bp  <  12  no
```

I verified this empirically: a paired-end variant with `truncLen=c(140,140)`
dropped *every* read in `filterAndTrim` (because most post-cutadapt reads
are already shorter than 140 bp), and lower `truncLen` values would push
the available overlap below DADA2's minimum. The Earth Microbiome 515F/806R
+ 2×150 layout is well known not to merge cleanly without bespoke
parameters.

Possible future paths (not adopted now to keep the panel uniform with
defaults):
- run DADA2 with `mergePairs(..., minOverlap=8)`, risks false merges;
- skip cutadapt and use `filterAndTrim(..., trimLeft=c(19,20),
  truncLen=c(140,140))`, same overlap budget, no gain;
- ignore primers, use full 150 bp, overlap 150+150−290 = 10 bp, still
  below default.

The original single-end (`truncLen = 125`) is the right baseline. ASVs are
125 bp covering the front half of V4, that is what every classifier in
the panel sees, so the comparison is fair.

### 2.4 zuo: paired-end was disabled in error yes fixed

`slurm/zuo_amplicon_pipeline.sbatch` reasons:

> V4 252bp amplicon + 2x150 MiSeq + 39bp primers leaves ~111bp post-trim
> reads with no overlap.

…and runs DADA2 single-end with `truncLen = 110`.

**The math is wrong.** Zuo reads were already primer-trimmed by the
submitter (the script itself notes this and skips the primer step in
cutadapt). So the reads start at the V4 insert and are the full 150 bp:

```
R1 + R2  − insert  =  150 + 150 − 252  =  48 bp  ≥ 12  OK
```

**Fix applied** (`slurm/zuo_amplicon_pipeline_paired.sbatch`, run
2026-05-04 as job 9168, COMPLETED in 23 min):

- Same download + Q-only cutadapt as before.
- DADA2 paired-end with `truncLen = c(145, 140)`.
- `mergePairs` (defaults).

Result:

| Metric          | Before (single-end, truncLen=110) | After (paired-end) |
|-----------------|------------------------------------|--------------------|
| Total ASVs      | 1 769                              | **3 218**          |
| ASV length      | 110 bp (all)                       | **252–253 bp** (peak; 2633 of 3218 at 253) |
| Samples         | 42                                 | 42                 |

The new `rep_seqs_zuo.fasta` is in place at
`data/real_world_input/rep_seqs_zuo.fasta`. The original DADA2 outputs in
`zuo/dada2_output/` were overwritten by the new run; the asv_table_ids.csv
keys still match `metadata/paired_samples.csv` (verified all 42 samples).

---

## 3. Why merging is impossible for coral but trivial for zuo

Both studies are V4 2×150. The difference is *what the 150 bp covers*:

- **coral (raw + primers)**: read[0:19/20] is primer; only ~131 bp of each
  read is insert. R1+R2 covers 261/252 of the insert ⇒ overlap = 9 bp.
- **zuo (submitter pre-trimmed)**: read[0:150] is already insert.
  R1+R2 covers 300/252 of the insert ⇒ overlap = 48 bp.

This is also why `truncLen=140` failed for coral: most reads were
already shorter than 140 bp post-cutadapt and `filterAndTrim` dropped
them all.

---

## 4. Effect on downstream evaluation

For zuo: the new ASVs cover the **full V4 variable region** instead of
just the first ~110 bp, restoring ~140 bp of variable-region signal.
Species-level scores for zuo are expected to rise across all 8
classifiers; the relative ordering of classifiers is the comparable axis
for the manuscript.

For coral, virtanen, mantri: no change.

---

## 5. References (web)

- DADA2 `mergePairs` defaults `minOverlap=12`, `maxMismatch=0`,
  https://rdrr.io/bioc/dada2/man/mergePairs.html ;
  https://github.com/benjjneb/dada2/issues/852
- Earth Microbiome Project 16S protocol (515F/806R; insert ≈ 252 bp),
  https://earthmicrobiome.ucsd.edu/protocols-and-standards/16s/
- Illumina V4 amplicon length / 2×150 limitations,
  https://www.illumina.com/content/dam/illumina-marketing/documents/products/appnotes/appnote_miseq_16S.pdf
- V3-V4 amplicon length (Klindworth et al. 2013, Nucleic Acids Res.)
