suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

message('== QC table synthesis ==')
dir.create('results/next/qc', recursive = TRUE, showWarnings = FALSE)
dir.create('results/analysis_tables/qc', recursive = TRUE, showWarnings = FALSE)

meta <- readr::read_tsv('config/sample_metadata_curated.tsv', show_col_types = FALSE) %>%
  dplyr::filter(tolower(include) %in% c('yes','y','true','1'))
map <- readr::read_csv('results/tables/salmon_mapping_rates.csv', show_col_types = FALSE)
qc <- meta %>%
  dplyr::left_join(map, by = 'sample_id') %>%
  dplyr::mutate(percent_mapped = as.numeric(percent_mapped))

low <- qc %>%
  dplyr::filter(percent_mapped < 50) %>%
  dplyr::arrange(percent_mapped)

readr::write_tsv(qc, 'results/next/qc/sample_qc_with_metadata.tsv')
readr::write_tsv(low, 'results/next/qc/low_mapping_samples_lt50.tsv')
readr::write_tsv(qc, 'results/analysis_tables/qc/sample_qc_with_metadata.tsv')
readr::write_tsv(low, 'results/analysis_tables/qc/low_mapping_samples_lt50.tsv')

message('== DONE QC table synthesis ==')
