
args <- commandArgs(trailingOnly = TRUE)

# packages
suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(DECIPHER)
})  
set.seed(123)




# get trained classifier
trainingSet <- args[1]
trainingSet <- readRDS(trainingSet)

# get test seqs
test_seq  <- args[2]
test_seq <- readDNAStringSet(test_seq)

# perform classification
confidence_threshold <- args[4]
n_cpus <- args[5]
ids <- IdTaxa(test = test_seq, 
              trainingSet = trainingSet, 
              type = "extended", 
              strand = "both", 
              threshold = confidence_threshold,
              processors = n_cpus,
              fullLength = 0.99, # default is 0 (if using v_region classifiers we are good to go with 0.99)
              verbose = TRUE)
 

# Transform to a dataframe
n_seq <- length(ids)
df_rows <- list()

# Go through all the elements of the list
for(i in 1:n_seq){
  sequence_id <- names(ids[i])
  taxonomy <- ids[[i]]$taxon
  confidence <- ids[[i]]$confidence
  rank <- ids[[i]]$rank
  # mutate taxonomy to unclassified if matched
  if (taxonomy[2] == "unclassified_Root") {
    taxonomy <- rep("Unclassified", length(rank))
  }
  df_rows[[i]] = data.frame(sequence_id, taxonomy, confidence, taxo_level=rank)
}

# Bind rows of all seqs
df <- rbindlist(df_rows) %>% filter(taxo_level !="rootrank")

# To wide dataframe and then mutate
df <- df %>%
  pivot_wider(
    names_from = taxo_level, 
    values_from = c(taxonomy, confidence),
    names_glue = "{taxo_level}_{.value}"
  ) %>%
  rename_with(~ gsub("_taxonomy$", "", .x))

# Desired order
desired_order <- c("sequence_id", 
                   "domain", "domain_confidence", 
                   "phylum", "phylum_confidence", 
                   "class", "class_confidence",
                   "order", "order_confidence",
                   "family", "family_confidence",
                   "genus", "genus_confidence",
                   "species", "species_confidence")

# Reorder columns
df <- df %>% select(any_of(desired_order))

# finalize df
df <- df %>%
  mutate(across(where(is.character), ~ replace_na(.x, "Unclassified"))) %>%
  mutate(across(where(is.numeric), ~ replace_na(.x, 0))) %>%
  mutate("domain_confidence" = if_else(domain != "Unclassified", domain_confidence, 0))


# save to txt file
out_df <- args[3]
fwrite(x = df, file = out_df, sep = "\t")