#!/usr/bin/env Rscript
# Train an IDTAXA model (DECIPHER::LearnTaxa) for LOCO experiments.
# Usage: Rscript train_idtaxa.R <fasta> <taxonomy_tsv> <output_RData>

suppressMessages({
  library(DECIPHER)
  library(Biostrings)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) stop("Usage: train_idtaxa.R <fasta> <taxonomy.tsv> <out.RData>")
fa <- args[1]; tx <- args[2]; out <- args[3]

cat("[idtaxa] reading fasta:", fa, "\n")
seqs <- readDNAStringSet(fa)

cat("[idtaxa] reading taxonomy:", tx, "\n")
df <- read.delim(tx, sep="\t", header=TRUE, stringsAsFactors=FALSE,
                 check.names=FALSE, quote="")
# Some files are headerless or have "Feature ID"; handle both.
if (!"Feature ID" %in% colnames(df) && !"id" %in% colnames(df)) {
  colnames(df)[1:2] <- c("Feature ID", "Taxon")
}
id_col <- if ("Feature ID" %in% colnames(df)) "Feature ID" else "id"
tax_map <- setNames(df$Taxon, df[[id_col]])
seq_names <- names(seqs)
# Strip any trailing tax info from fasta headers (defensive)
seq_names <- sapply(seq_names, function(x) strsplit(x, "\t", fixed=TRUE)[[1]][1])
names(seqs) <- seq_names
keep <- seq_names %in% names(tax_map)
seqs <- seqs[keep]
labels <- tax_map[names(seqs)]
# Filter out sequences with no usable k-mers (too short or all-N) and missing labels
slens <- width(seqs)
keep <- slens >= 50 & !is.na(labels) & labels != ""
seqs <- seqs[keep]; labels <- labels[keep]
cat("[idtaxa] dropped", sum(!keep), "sequences with width<50 or NA label\n")
# Filter sequences that are nearly all N (DECIPHER's k-mer needs unambiguous bases)
acgt_frac <- letterFrequency(seqs, c("A","C","G","T"), as.prob=TRUE)
n_acgt <- rowSums(acgt_frac)
keep2 <- n_acgt >= 0.5
seqs <- seqs[keep2]; labels <- labels[keep2]
cat("[idtaxa] dropped", sum(!keep2), "sequences with <50% A/C/G/T content\n")
# IDTAXA expects "Root;Domain;Phylum;..." format
clean_labels <- paste0("Root;", gsub("; ", ";", labels))
cat("[idtaxa] training LearnTaxa on", length(seqs), "sequences...\n")
trainingSet <- LearnTaxa(seqs, clean_labels, verbose=TRUE)
cat("[idtaxa] saving to", out, "\n")
save(trainingSet, file=out)
cat("[idtaxa] done\n")
