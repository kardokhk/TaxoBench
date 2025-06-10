
args <- commandArgs(trailingOnly = TRUE)

# packages
suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(DECIPHER)
})  
set.seed(123)



# function to create taxid file
decipher_build_taxid <- function(groups, taxo_levels) {

  # Define the ranks and create a parent-child relationship
  ranks <- taxo_levels
  
  # get taxonomy
  taxonomy <- data.frame(taxo_groups = groups) |>
    mutate(taxo_groups = str_replace(taxo_groups, "Root;", "")) |>
    tidyr::separate(taxo_groups, sep = "; ", into = taxo_levels) |>
    dplyr::distinct()

  # Process the data and write to a new file
  taxid <- taxonomy %>%
    pivot_longer(
      cols = all_of(ranks),
      names_to = "rank",
      values_to = "name"
    )  %>%
    filter(name != "") %>%
    mutate(
      ID = row_number(),
      parent_name = if_else(rank == ranks[1], NA, lag(name)),
      parent_rank = if_else(rank == ranks[1], NA, lag(rank))
    ) %>%
    select(
      ID,
      parent_name,
      parent_rank,
      name,
      rank
    )|> 
    distinct(name, rank, .keep_all = TRUE)

  taxid_parents <- taxid |>
    select(parent_ID = ID, parent_name = name, parent_rank = rank)

  taxid <- left_join(taxid, taxid_parents, by = c("parent_name", "parent_rank")) |>
    mutate(rank_number = match(rank, taxo_levels),
           parent_ID = if_else(is.na(parent_ID), 0, parent_ID)) |>
    select(Index=ID, Name=name, Parent=parent_ID, Level=rank_number, Rank=rank) 

  root_taxid <- data.frame(Index=0, Name="Root", Parent = -1, Level=0, Rank= "rootrank")
  taxid <- bind_rows(root_taxid, taxid)
  
  return(taxid)
}


## import training sequences
ref_seqs_path <- args[1]
seqs <- readDNAStringSet(ref_seqs_path)
names(seqs) <- sapply(X = strsplit(names(seqs), "\t"), FUN = function(x) { x[2] })

# remove any gaps
seqs <- RemoveGaps(seqs)

# reorient seqs
seqs <- OrientNucleotides(seqs)

# taxonomic level
taxo_levels = c("domain", "phylum", "class", "order", "family", "genus", "species")

# Get taxonomy assignments
groups <- names(seqs)
groupCounts <- table(groups)
u_groups <- names(groupCounts)
taxid <- decipher_build_taxid(u_groups, taxo_levels)

# add "Root;" to taxonomy lineage
seqs_name <- paste0("Root; ", names(seqs))

# train the classifier
trainingSet <- LearnTaxa(train = seqs, taxonomy = seqs_name, rank = taxid)

# save trained classifer
trainset_path <- args[2]
saveRDS(trainingSet, trainset_path)

