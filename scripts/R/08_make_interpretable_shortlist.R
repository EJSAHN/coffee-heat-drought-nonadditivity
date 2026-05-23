library(readr)
library(dplyr)
library(stringr)

outdir <- "results/paper_pack_v1/tables"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

read_panel <- function(group) {
  path <- file.path("results/next/tables", paste0(group, "_candidate_panel_top20_per_class.tsv"))
  readr::read_tsv(path, show_col_types = FALSE) %>%
    mutate(group = group)
}

read_enrich <- function(group) {
  path <- file.path("results/next/tables", paste0(group, "_module_enrichment.tsv"))
  readr::read_tsv(path, show_col_types = FALSE) %>%
    mutate(group = group)
}

read_class <- function(group) {
  path <- file.path("results/next/tables", paste0(group, "_nonadditivity_class_summary_pct.tsv"))
  readr::read_tsv(path, show_col_types = FALSE) %>%
    mutate(group = group)
}

groups <- c("arabica_Icatu", "canephora_CL153")

panels <- bind_rows(lapply(groups, read_panel))
enrich <- bind_rows(lapply(groups, read_enrich))
classes <- bind_rows(lapply(groups, read_class))

interpretable <- panels %>%
  mutate(
    product_clean = str_replace_all(product, "%2C", ","),
    is_unannotated = is.na(product_clean) |
      str_detect(str_to_lower(product_clean), "unannotated|uncharacterized|hypothetical|unknown")
  ) %>%
  filter(
    nonadd_class != "not_significant_or_small",
    !is.na(padj),
    padj < 0.05
  ) %>%
  mutate(
    abs_nonadd = abs(nonadditive_score),
    priority_score = -log10(padj + 1e-300) + abs_nonadd
  ) %>%
  arrange(group, target_stage, nonadd_class, desc(priority_score))

interpretable_non_other <- interpretable %>%
  filter(functional_module != "other_or_unannotated" | !is_unannotated)

top_by_module <- interpretable_non_other %>%
  group_by(group, target_stage, nonadd_class, functional_module) %>%
  slice_max(priority_score, n = 5, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(group, target_stage, nonadd_class, functional_module, desc(priority_score))

top_overall <- interpretable_non_other %>%
  group_by(group, target_stage) %>%
  slice_max(priority_score, n = 25, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(group, target_stage, desc(priority_score))

sig_modules <- enrich %>%
  mutate(sig_level = case_when(
    padj < 0.05 ~ "FDR<0.05",
    padj < 0.10 ~ "FDR<0.10",
    TRUE ~ "not_FDR_significant"
  )) %>%
  arrange(group, target_stage, padj)

readr::write_tsv(interpretable, file.path(outdir, "paper_all_candidate_panel_ranked.tsv"))
readr::write_tsv(interpretable_non_other, file.path(outdir, "paper_interpretable_candidate_panel.tsv"))
readr::write_tsv(top_by_module, file.path(outdir, "paper_top5_candidates_per_module.tsv"))
readr::write_tsv(top_overall, file.path(outdir, "paper_top25_candidates_per_stage.tsv"))
readr::write_tsv(sig_modules, file.path(outdir, "paper_module_enrichment_ranked.tsv"))
readr::write_tsv(classes, file.path(outdir, "paper_nonadditivity_class_numbers.tsv"))

message("DONE paper shortlist tables")