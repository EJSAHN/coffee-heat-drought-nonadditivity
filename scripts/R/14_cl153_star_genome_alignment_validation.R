suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(DESeq2)
  library(openxlsx)
})

# Rebind tidyverse verbs after Bioconductor packages are loaded.
select <- dplyr::select
filter <- dplyr::filter
mutate <- dplyr::mutate
arrange <- dplyr::arrange
count <- dplyr::count
summarise <- dplyr::summarise
group_by <- dplyr::group_by
ungroup <- dplyr::ungroup
left_join <- dplyr::left_join
bind_rows <- dplyr::bind_rows
recode <- dplyr::recode
pivot_wider <- tidyr::pivot_wider
pivot_longer <- tidyr::pivot_longer

message("== CL153 STAR genome-alignment validation summary")

out_dir <- "results/qc_validation/cl153_star_genome_alignment"
star_dir <- file.path(out_dir, "star_alignments")
manifest_path <- file.path(out_dir, "cl153_star_manifest.tsv")
if (!file.exists(manifest_path)) stop("Missing manifest: ", manifest_path)
manifest <- readr::read_tsv(manifest_path, show_col_types = FALSE) %>%
  mutate(
    sample_id = as.character(sample_id),
    water = factor(water, levels = c("WW", "SWD")),
    stage = factor(stage, levels = c("T25", "T37", "T42", "REC14"))
  )

parse_star_log <- function(sample_id) {
  p <- file.path(star_dir, sample_id, "Log.final.out")
  if (!file.exists(p)) stop("Missing STAR Log.final.out for ", sample_id)
  x <- readLines(p, warn = FALSE)
  get_val <- function(pattern) {
    line <- x[grepl(pattern, x, fixed = TRUE)]
    if (!length(line)) return(NA_character_)
    sub(".*\\|\\s*", "", line[1])
  }
  pct_num <- function(z) suppressWarnings(as.numeric(gsub("%", "", z)))
  tibble::tibble(
    sample_id = sample_id,
    input_reads = suppressWarnings(as.numeric(get_val("Number of input reads"))),
    avg_input_read_length = suppressWarnings(as.numeric(get_val("Average input read length"))),
    uniquely_mapped_reads_number = suppressWarnings(as.numeric(get_val("Uniquely mapped reads number"))),
    uniquely_mapped_reads_pct = pct_num(get_val("Uniquely mapped reads %")),
    multi_mapped_reads_number = suppressWarnings(as.numeric(get_val("Number of reads mapped to multiple loci"))),
    multi_mapped_reads_pct = pct_num(get_val("% of reads mapped to multiple loci")),
    too_many_loci_pct = pct_num(get_val("% of reads mapped to too many loci")),
    unmapped_too_short_pct = pct_num(get_val("% of reads unmapped: too short")),
    unmapped_other_pct = pct_num(get_val("% of reads unmapped: other")),
    total_mapped_pct_approx = uniquely_mapped_reads_pct + multi_mapped_reads_pct
  )
}

mapping <- bind_rows(lapply(manifest$sample_id, parse_star_log)) %>%
  left_join(manifest %>% select(sample_id, run_accession, water, stage, replicate), by = "sample_id") %>%
  select(sample_id, run_accession, water, stage, replicate, everything())

read_gene_counts <- function(sample_id) {
  p <- file.path(star_dir, sample_id, "ReadsPerGene.out.tab")
  if (!file.exists(p)) stop("Missing ReadsPerGene.out.tab for ", sample_id)
  d <- readr::read_tsv(p, col_names = c("gene_id", "unstranded", "stranded_forward", "stranded_reverse"), show_col_types = FALSE)
  d %>%
    filter(!grepl("^N_", gene_id)) %>%
    transmute(gene_id = as.character(gene_id), !!sample_id := as.integer(unstranded))
}

counts_list <- lapply(manifest$sample_id, read_gene_counts)
count_df <- Reduce(function(x, y) full_join(x, y, by = "gene_id"), counts_list) %>%
  mutate(across(-gene_id, ~replace_na(.x, 0L)))
count_mat <- as.data.frame(count_df)
rownames(count_mat) <- count_mat$gene_id
count_mat$gene_id <- NULL
count_mat <- as.matrix(count_mat[, manifest$sample_id, drop = FALSE])
storage.mode(count_mat) <- "integer"

coldata <- manifest %>%
  select(sample_id, water, stage, replicate) %>%
  as.data.frame()
rownames(coldata) <- coldata$sample_id
coldata$sample_id <- NULL

# Low-count filter matching the primary workflow in spirit.
keep <- rowSums(count_mat >= 5) >= max(2, ceiling(0.10 * ncol(count_mat)))
count_mat_f <- count_mat[keep, , drop = FALSE]
message("Features kept in STAR gene-count model: ", nrow(count_mat_f))

dds <- DESeqDataSetFromMatrix(countData = count_mat_f, colData = coldata, design = ~ water + stage + water:stage)
dds <- DESeq(dds)
norm_counts <- counts(dds, normalized = TRUE)
log_norm <- log2(norm_counts + 1)

# Helper mean per cell
cell_mean <- function(water, stage) {
  idx <- rownames(colData(dds))[colData(dds)$water == water & colData(dds)$stage == stage]
  rowMeans(log_norm[, idx, drop = FALSE])
}

ww_t25 <- cell_mean("WW", "T25")
swd_t25 <- cell_mean("SWD", "T25")

res_stage <- function(stage_name) {
  interaction_name <- paste0("waterSWD.stage", stage_name)
  rn <- resultsNames(dds)
  if (!(interaction_name %in% rn)) {
    # DESeq2 may encode ':' in names in unusual versions; search robustly.
    hit <- rn[grepl("water.*SWD", rn) & grepl(stage_name, rn)]
    if (length(hit) != 1) stop("Could not identify interaction coefficient for ", stage_name, ". resultsNames: ", paste(rn, collapse=", "))
    interaction_name <- hit
  }
  res <- as.data.frame(results(dds, name = interaction_name)) %>%
    tibble::rownames_to_column("gene_id")
  ww_s <- cell_mean("WW", stage_name)
  swd_s <- cell_mean("SWD", stage_name)
  nonadd <- (swd_s - ww_t25) - (ww_s - ww_t25) - (swd_t25 - ww_t25)
  combo <- swd_s - ww_t25
  heat <- ww_s - ww_t25
  drought <- swd_t25 - ww_t25

  res %>%
    mutate(
      target_stage = stage_name,
      nonadditivity_score = nonadd[gene_id],
      combo_effect = combo[gene_id],
      heat_effect = heat[gene_id],
      drought_effect = drought[gene_id],
      is_nonadditive = !is.na(padj) & padj < 0.05 & abs(nonadditivity_score) >= 1,
      class = case_when(
        is_nonadditive & abs(combo_effect) >= 1 & abs(heat_effect) < 0.5 & abs(drought_effect) < 0.5 ~ "Emergent Only",
        is_nonadditive & nonadditivity_score > 0 ~ "Synergistic Positive",
        is_nonadditive & nonadditivity_score < 0 ~ "Antagonistic or Buffered",
        TRUE ~ "Not significant or small"
      )
    )
}

all_res <- bind_rows(lapply(c("T37", "T42", "REC14"), res_stage))
counts_by_stage <- all_res %>%
  filter(is_nonadditive) %>%
  count(target_stage, class, name = "n") %>%
  tidyr::complete(target_stage = c("T37","T42","REC14"),
                  class = c("Synergistic Positive","Antagonistic or Buffered","Emergent Only"),
                  fill = list(n = 0)) %>%
  group_by(target_stage) %>%
  mutate(total_nonadditive = sum(n), pct = ifelse(total_nonadditive > 0, 100*n/total_nonadditive, 0)) %>%
  ungroup() %>%
  arrange(target_stage, class)

total_counts <- counts_by_stage %>%
  group_by(target_stage) %>%
  summarise(star_gene_count_total = unique(total_nonadditive), .groups = "drop")

# Pull primary Salmon all-stage counts if available.
# The low-mapping sensitivity table may be either long-format
# (model / target_stage / total_nonadditive) or wide-format
# (target_stage / all_24_libraries / drop_5_libraries_lt50pct).
primary_path <- "results/qc_validation/cl153_low_mapping_exclusion/nonadditive_gene_count_comparison.tsv"
if (file.exists(primary_path)) {
  primary <- readr::read_tsv(primary_path, show_col_types = FALSE)

  if (all(c("model", "target_stage", "total_nonadditive") %in% names(primary))) {
    primary_tot <- primary %>%
      filter(model == "all_24_libraries") %>%
      transmute(target_stage, salmon_tximport_total = total_nonadditive)
  } else if (all(c("target_stage", "all_24_libraries") %in% names(primary))) {
    primary_tot <- primary %>%
      transmute(target_stage, salmon_tximport_total = all_24_libraries)
  } else {
    warning("Could not recognize primary Salmon comparison table columns: ", paste(names(primary), collapse = ", "))
    primary_tot <- tibble::tibble(target_stage = character(), salmon_tximport_total = numeric())
  }

  compare <- total_counts %>%
    left_join(primary_tot, by = "target_stage") %>%
    mutate(
      star_minus_salmon = ifelse(is.na(salmon_tximport_total), NA_real_, star_gene_count_total - salmon_tximport_total),
      star_to_salmon_ratio = ifelse(is.na(salmon_tximport_total) | salmon_tximport_total == 0, NA_real_, star_gene_count_total / salmon_tximport_total)
    )
} else {
  compare <- total_counts
}

summary <- tibble::tibble(
  metric = c("n_libraries", "n_features_kept", "median_STAR_total_mapped_pct_approx", "min_STAR_total_mapped_pct_approx", "median_STAR_unique_pct", "min_STAR_unique_pct"),
  value = c(nrow(manifest), nrow(count_mat_f),
            median(mapping$total_mapped_pct_approx, na.rm = TRUE),
            min(mapping$total_mapped_pct_approx, na.rm = TRUE),
            median(mapping$uniquely_mapped_reads_pct, na.rm = TRUE),
            min(mapping$uniquely_mapped_reads_pct, na.rm = TRUE))
)

readr::write_tsv(mapping, file.path(out_dir, "star_mapping_rate_summary.tsv"))
readr::write_tsv(counts_by_stage, file.path(out_dir, "star_nonadditive_class_counts.tsv"))
readr::write_tsv(compare, file.path(out_dir, "star_vs_salmon_nonadditive_count_comparison.tsv"))
readr::write_tsv(summary, file.path(out_dir, "star_reference_summary.tsv"))
readr::write_tsv(all_res, file.path(out_dir, "star_all_interaction_results.tsv"))

wb <- openxlsx::createWorkbook()
add_sheet <- function(name, df) {
  openxlsx::addWorksheet(wb, name)
  openxlsx::writeData(wb, name, df)
  openxlsx::freezePane(wb, name, firstRow = TRUE)
}
add_sheet("README", tibble::tibble(
  item = c("Purpose", "Reference framework", "Counting mode", "Interpretation"),
  description = c(
    "CL153 genome-alignment validation using the same AUK_PRJEB4211_v1 assembly framework.",
    "STAR alignment to Ensembl Plants/Coffee Genome Hub-derived C. canephora AUK_PRJEB4211_v1 genome and GTF converted from GFF3.",
    "STAR GeneCounts unstranded counts; DESeq2 water + stage + water:stage model.",
    "This analysis evaluates whether transcript-level Salmon mapping-rate estimates materially alter CL153 non-additive gene counts."
  )
))
add_sheet("STAR_mapping_summary", mapping)
add_sheet("STAR_nonadd_counts", counts_by_stage)
add_sheet("STAR_vs_Salmon_counts", compare)
add_sheet("STAR_summary", summary)
openxlsx::saveWorkbook(wb, file.path(out_dir, "CL153_STAR_genome_alignment_validation.xlsx"), overwrite = TRUE)

message("== DONE CL153 STAR genome-alignment validation summary")
message("Outputs written to: ", out_dir)
