
args <- commandArgs(trailingOnly = TRUE)

# packages
library(data.table)
library(tidyverse)
set.seed(123)


## Import GG2.2024 taxonomy
gg2 <- data.table::fread(input = args[1], sep = "\t", check.names = FALSE)

## Import results
# MMSEQS = query,target,fident,alnlen,qlen,tlen,qcov,tcov,evalue,bits
df <- data.table::fread(input = args[2], sep = "\t", check.names = FALSE)

# join taxonomy info
df <- df %>%
  rename("query"=V1, "target"=V2, "id"=V3, "alnlen"=V4, "qlen"=V5, "rlen"=V6, "qcov"=V7, "rcov"=V8, "evalue"=V9, "bits"=V10) %>%
  inner_join(gg2, by = c("target"="Feature ID"))



### Identify if there are any best hits for a given query that has hits where a query sequence has 100% identity over >100% coverage
exact_match_pass <- df %>%
  filter(id == 1.00, qcov == 1.00) %>%
  group_by(query) %>%
  mutate(n_taxa = n()) %>%
  filter(n_taxa == 1) %>%
  ungroup() %>%
  distinct(Taxon, .keep_all = TRUE) %>%
  mutate(confidence = "100% identity + 100% coverage + unambiguous taxonomy") %>%
  select(query, Taxon, confidence)

### Identify if there are any best hits for a given query that has a only hit where a query sequence has 100% identity over >99% coverage
best_hits_pass_1 <- df %>%
  filter(!query %in% exact_match_pass$query) %>%
  filter(id == 1.00, qcov >= 0.99) %>%
  group_by(query) %>%
  mutate(n_valid_hits = n()) %>%
  filter(n_valid_hits == 1) %>%
  mutate(confidence = "100% identity + >99% coverage + only best hit") %>%
  select(query, Taxon, confidence)

### Identify if there are any best hits for a given query that has a only hit where a query sequence has >99% identity over >99% coverage
best_hits_pass_2 <- df %>%
  filter(!query %in% c(exact_match_pass$query, best_hits_pass_1$query)) %>%
  filter(id >= 0.99, qcov >= 0.99) %>%
  group_by(query) %>%
  mutate(n_valid_hits = n()) %>%
  filter(n_valid_hits == 1) %>%
  mutate(confidence = ">99% identity + >99% coverage + only best hit") %>%
  select(query, Taxon, confidence)

### Identify if there are any hits for a given query that has a query sequence with >99% identity over >99% coverage and first top hit starts with ^G000
# ^G00 genomes represent well-characterized genomes with reliable taxonomic information derived from WOL.
best_hits_pass_3 <- df %>%
  filter(!query %in% c(exact_match_pass$query, best_hits_pass_1$query, best_hits_pass_2$query)) %>%
  filter(id >= 0.99, qcov >= 0.99) %>%
  group_by(query) %>%
  mutate(
    n_valid_hits = n(),
    is_G00 = grepl("^G00", target)
  ) %>%
  arrange(desc(bits), desc(is_G00)) %>%
  slice_head(n = 3) %>%
  filter(grepl("^G00", target)) %>%
  select(-is_G00) %>%
  slice_head(n = 1) %>%
  mutate(confidence = ">99% identity + >99% coverage + WOL genome in top hits") %>%
  select(query, Taxon, confidence)

## Majority rule attempt [identity >= 97% and qcov >= 90% and consensus assignment based on >=60% majority rule]
majority_rule_pass <- df %>%
  filter(!query %in% c(exact_match_pass$query, best_hits_pass_1$query, best_hits_pass_2$query, best_hits_pass_3$query)) %>%
  filter(id >= 0.97, qcov >= 0.90) %>%
  group_by(query) %>%
  mutate(
    min_id = min(id),
    max_id = max(id),
    min_qcov = min(qcov),
    max_qcov = max(qcov),
    min_alnlen = min(alnlen),
    max_alnlen = max(alnlen)
  ) %>%
  group_by(query, Taxon) %>%
  summarise(
    No_matched = n(),
    min_id = first(min_id),
    max_id = first(max_id),
    min_qcov = first(min_qcov),
    max_qcov = first(max_qcov),
    min_alnlen = first(min_alnlen),
    max_alnlen = first(max_alnlen),
    .groups = "drop_last"
  ) %>%
  group_by(query) %>%
  mutate(
    No_total = sum(No_matched),
    Freq = No_matched / No_total
  ) %>%
  filter(Freq == max(Freq)) %>%
  filter(Freq >= 0.60) %>%
  ungroup() %>%
  mutate(
    identity_range = paste0(round(min_id, 1), "-", round(max_id, 1)),
    qcov_range = paste0(round(min_qcov, 1), "-", round(max_qcov, 1)),
    alnlen_range = paste0(min_alnlen, "-", max_alnlen),
    confidence = "Majority rule > 60%"
  ) %>%
  select(query, Taxon, confidence)

### Identify if there are any best hits for a given query For each query where second-best hit is at least N bits lower [bits[1] - bits[2] >= N]
# bits score if log2s scaled.. so difference of 10 means 1024X significant
# taxon aware bits difference
bits_difference_score <- 10  # bit score difference threshold
best_hits_pass <- df %>%
  filter(!query %in% c(majority_rule_pass$query, exact_match_pass$query, best_hits_pass_1$query, best_hits_pass_2$query, best_hits_pass_3$query)) %>%
  filter(id >= 0.97, qcov >= 0.90) %>%
  group_by(query) %>%
  arrange(desc(bits)) %>%
  mutate(
    species_name = str_extract(Taxon, "s__[^; ]+"),
    is_unclassified = is.na(species_name) | species_name == "s__",
    species_name = ifelse(is_unclassified, NA, species_name),
    top_species = first(species_name),
    is_diff = species_name != top_species & !is.na(species_name)
  ) %>%
  # Keep the top hit and the first taxon-different hit (excluding unclassified)
  filter(row_number() == 1 | is_diff) %>%
  slice_head(n = 2) %>%
  summarise(
    top_hit = first(target),
    top_bits = first(bits),
    second_bits = ifelse(n() > 1, bits[2], NA),
    Taxon = first(Taxon),
    .groups = "drop"
  ) %>%
  mutate(
    bit_diff = top_bits - second_bits,
    confident = ifelse(is.na(second_bits) | bit_diff >= bits_difference_score, TRUE, FALSE)
  ) %>%
  filter(confident) %>%
  mutate(confidence = "Top hit ≥ 10-bit score above first differing taxon") %>%
  select(query, Taxon, confidence)


## not classified with any of the above rules
not_classified_mocks <- df %>% 
  filter(!query %in% c(exact_match_pass$query, best_hits_pass_1$query, best_hits_pass_2$query, best_hits_pass_3$query, 
                       majority_rule_pass$query, best_hits_pass$query))
cat("Following mock sequences could not be classified reliably:-\n")
cat(" ", not_classified_mocks %>% pull(query) %>% unique() %>% sort())

## combine results
final <- bind_rows(exact_match_pass, best_hits_pass_1, best_hits_pass_2, best_hits_pass_3, majority_rule_pass, best_hits_pass)


## write
data.table::fwrite(x = final, file = args[3], sep = "\t")
