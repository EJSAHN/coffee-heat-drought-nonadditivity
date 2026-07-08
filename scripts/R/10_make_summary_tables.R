suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

base <- 'results/analysis_tables/tables'
dir.create(base, recursive = TRUE, showWarnings = FALSE)

class_tbl <- readr::read_tsv(file.path(base, 'analysis_class_summary_long.tsv'), show_col_types = FALSE)
stage_totals <- readr::read_tsv(file.path(base, 'analysis_stage_nonadditive_totals.tsv'), show_col_types = FALSE)

class_within <- class_tbl %>%
  dplyr::left_join(stage_totals, by = c('group', 'target_stage')) %>%
  dplyr::mutate(
    within_stage_pct = 100 * n / nonadditive_genes,
    target_stage = factor(target_stage, levels = c('T37', 'T42', 'REC14'))
  ) %>%
  dplyr::arrange(group, target_stage, nonadd_class)

readr::write_tsv(class_within, file.path(base, 'analysis_class_summary_within_stage_pct.tsv'))
message('DONE summary tables')
