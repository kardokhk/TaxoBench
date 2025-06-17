#!/usr/bin/env Rscript

# idTaxa classification script for taxonova benchmarking framework
# Usage: Rscript run_idtaxa_classification.R <input_fasta> <model_file> <output_tsv> <confidence_threshold>
# model_file can be either .rds or .RData format

library(DECIPHER)
library(Biostrings)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4) {
  stop("Usage: Rscript run_idtaxa_classification.R <input_fasta> <model_file> <output_tsv> <confidence_threshold>")
}

input_fasta <- args[1]
model_rds <- args[2]
output_tsv <- args[3]
confidence_threshold <- as.numeric(args[4])

cat("Running idTaxa classification...\n")
cat("Input FASTA:", input_fasta, "\n")
cat("Model file:", model_rds, "\n")
cat("Output file:", output_tsv, "\n")
cat("Confidence threshold:", confidence_threshold, "\n")

# Load the pre-trained idTaxa model
cat("Loading idTaxa model...\n")

# Check file extension and load accordingly
if (grepl("\\.rds$", model_rds, ignore.case = TRUE)) {
  trainingSet <- readRDS(model_rds)
} else if (grepl("\\.RData$", model_rds, ignore.case = TRUE)) {
  # Load RData file (assumes the training set object is named 'trainingSet')
  load(model_rds)
  if (!exists("trainingSet")) {
    stop("RData file must contain an object named 'trainingSet'")
  }
} else {
  stop("Model file must be either .rds or .RData format")
}

# Read input sequences
cat("Reading input sequences...\n")
seqs <- readDNAStringSet(input_fasta)

# Run idTaxa classification
cat("Running classification...\n")
ids <- IdTaxa(seqs, trainingSet, strand="top", threshold=confidence_threshold, processors=NULL)

# Process results and create output dataframe
cat("Processing results...\n")
results <- data.frame(
  "Feature ID" = names(seqs),
  "Taxon" = character(length(ids)),
  "Confidence" = character(length(ids)),
  stringsAsFactors = FALSE
)

for (i in seq_along(ids)) {
  if (length(ids[[i]]$taxon) > 0) {
    # Get the taxonomic classification as-is from idTaxa
    taxon_levels <- ids[[i]]$taxon
    confidence_levels <- ids[[i]]$confidence
    
    # For raw output, preserve the original idTaxa format
    # Join taxonomy levels with semicolons (idTaxa returns them without prefixes usually)
    taxonomy_str = paste(taxon_levels, collapse="; ")
    
    # Convert confidence from percentage to decimal and format
    confidence_decimal <- confidence_levels / 100
    confidence_str <- paste(sprintf("%.4f", confidence_decimal), collapse="; ")
    
    results[i, "Taxon"] <- taxonomy_str
    results[i, "Confidence"] <- confidence_str
    
  } else {
    # No classification - unclassified
    results[i, "Taxon"] <- "Unclassified"
    results[i, "Confidence"] <- "0.0000"
  }
}

# Write results to output file
cat("Writing results to", output_tsv, "\n")
write.table(results, file=output_tsv, sep="\t", quote=FALSE, row.names=FALSE)

cat("idTaxa classification completed successfully!\n")
