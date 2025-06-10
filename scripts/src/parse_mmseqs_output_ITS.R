
args <- commandArgs(trailingOnly = TRUE)

# packages
library(data.table)
library(tidyverse)
set.seed(123)

## Import UNITE taxonomy
unite <- data.table::fread(input = args[1], sep = "\t", check.names = FALSE, header = FALSE)
colnames(unite) <- c("Feature ID", "Taxon")

## Import results
# VSEARCH = query+target+id+alnlen+qcov+tcov+mism+opens+qlo+qhi+tlo+thi+evalue+bits
# MMSEQS = query,target,fident,alnlen,qlen,tlen,qcov,tcov,evalue,bits
df <- data.table::fread(input = args[2], sep = "\t", check.names = FALSE)

# join taxonomy info
df <- df %>%
  rename("query"=V1, "target"=V2, "id"=V3, "alnlen"=V4, "qlen"=V5, "rlen"=V6, "qcov"=V7, "rcov"=V8, "evalue"=V9, "bits"=V10) %>%
  inner_join(unite, by = c("target"="Feature ID"))

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

### Identify if there are any best hits for a given query For each query where second-best hit is at least N bits lower [bits[1] - bits[2] >= N]
# bits score if log2s scaled.. so difference of 30 means very very significant
bits_difference_score <- 30  # bit score difference threshold
best_hits_pass <- df %>%
  filter(!query %in% exact_match_pass$query) %>%
  group_by(query) %>%
  arrange(desc(bits)) %>%
  mutate(rank = row_number()) %>%
  # Keep only top 2 hits to compare
  filter(rank <= 2) %>%
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
  mutate(confidence = "Top hit ≥ 30-bit score above next best match") %>%
  select(query, Taxon, confidence)

## combine results
final <- bind_rows(exact_match_pass, best_hits_pass)

## write
data.table::fwrite(x = final, file = args[3], sep = "\t")
