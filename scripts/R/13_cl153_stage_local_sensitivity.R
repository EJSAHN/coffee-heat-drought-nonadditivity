#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tximport)
  library(DESeq2)
})

# Robust tidyverse aliases: Bioconductor packages can mask dplyr verbs.
select <- dplyr::select
filter <- dplyr::filter
mutate <- dplyr::mutate
arrange <- dplyr::arrange
slice <- dplyr::slice
count <- dplyr::count
summarise <- dplyr::summarise
group_by <- dplyr::group_by
ungroup <- dplyr::ungroup
left_join <- dplyr::left_join
bind_rows <- dplyr::bind_rows
complete <- tidyr::complete
pivot_wider <- tidyr::pivot_wider

message("== CL153 stage-local interaction analysis")

root <- getwd()
out_root <- "results/qc_validation/cl153_stage_local"
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

water_levels <- c("WW", "SWD")
stage_levels_all <- c("T25", "T37", "T42", "REC14")
target_stages <- c("T37", "T42", "REC14")
nonadditive_classes <- c("synergistic_positive", "antagonistic_or_buffered", "emergent_only", "not_significant_or_small")

read_metadata <- function() {
  meta_path <- if (file.exists("config/sample_metadata_curated.tsv")) {
    "config/sample_metadata_curated.tsv"
  } else {
    "config/sample_metadata_auto.tsv"
  }
  if (!file.exists(meta_path)) stop("Missing metadata: run metadata step first")

  meta <- readr::read_tsv(meta_path, show_col_types = FALSE)
  required <- c("sample_id", "run_accession", "index_species", "genotype", "water", "stage", "replicate", "include")
  missing <- setdiff(required, names(meta))
  if (length(missing)) stop("Metadata missing columns: ", paste(missing, collapse = ", "))

  meta %>%
    mutate(include = tolower(as.character(include)),
           quant_file = file.path("results", "salmon", sample_id, "quant.sf")) %>%
    filter(include %in% c("yes", "y", "true", "1")) %>%
    filter(index_species == "canephora" | genotype == "CL153") %>%
    filter(file.exists(quant_file)) %>%
    mutate(water = as.character(water),
           stage = as.character(stage),
           replicate = as.character(replicate))
}

extract_json_number <- function(txt, key) {
  pat <- paste0('"', key, '"[[:space:]]*:[[:space:]]*([0-9.]+)')
  m <- regexpr(pat, txt, perl = TRUE)
  if (m[1] < 0) return(NA_real_)
  hit <- regmatches(txt, m)
  as.numeric(sub(pat, "\\1", hit, perl = TRUE))
}

read_mapping_rates <- function(samples) {
  tsv_path <- "results/analysis_tables/qc/sample_qc_with_metadata.tsv"
  csv_path <- "results/tables/salmon_mapping_rates.csv"

  if (file.exists(tsv_path)) {
    qc <- readr::read_tsv(tsv_path, show_col_types = FALSE)
    if (all(c("sample_id", "percent_mapped") %in% names(qc))) {
      return(qc %>% select(any_of(c("sample_id", "num_processed", "num_mapped", "percent_mapped"))))
    }
  }

  if (file.exists(csv_path)) {
    qc <- readr::read_csv(csv_path, show_col_types = FALSE)
    if (all(c("sample_id", "percent_mapped") %in% names(qc))) {
      return(qc %>% select(any_of(c("sample_id", "num_processed", "num_mapped", "percent_mapped"))))
    }
  }

  rows <- lapply(samples, function(s) {
    f <- file.path("results", "salmon", s, "aux_info", "meta_info.json")
    if (!file.exists(f)) {
      return(tibble(sample_id = s, num_processed = NA_real_, num_mapped = NA_real_, percent_mapped = NA_real_))
    }
    txt <- paste(readLines(f, warn = FALSE), collapse = "\n")
    tibble(sample_id = s,
           num_processed = extract_json_number(txt, "num_processed"),
           num_mapped = extract_json_number(txt, "num_mapped"),
           percent_mapped = extract_json_number(txt, "percent_mapped"))
  })
  bind_rows(rows)
}

safe_group_mean <- function(mat, samples) {
  samples <- intersect(samples, colnames(mat))
  if (length(samples) == 0) return(rep(NA_real_, nrow(mat)))
  rowMeans(mat[, samples, drop = FALSE], na.rm = TRUE)
}

find_interaction_coef <- function(coefs, stage) {
  hits <- coefs[grepl("water", coefs, ignore.case = TRUE) &
                  grepl("SWD", coefs, ignore.case = TRUE) &
                  grepl("stage", coefs, ignore.case = TRUE) &
                  grepl(stage, coefs, ignore.case = TRUE)]
  if (length(hits) >= 1) return(hits[1])
  NA_character_
}

classify_nonadd <- function(nonadd, combo, heat, drought, padj) {
  out <- rep("not_significant_or_small", length(nonadd))
  sig <- !is.na(padj) & padj < 0.05 & !is.na(nonadd) & abs(nonadd) >= 1
  out[sig & nonadd > 0] <- "synergistic_positive"
  out[sig & nonadd < 0] <- "antagonistic_or_buffered"
  emergent <- sig & !is.na(combo) & abs(combo) >= 1 & !is.na(heat) & abs(heat) < 0.5 & !is.na(drought) & abs(drought) < 0.5
  out[emergent] <- "emergent_only"
  out
}

pretty_class <- function(x) {
  dplyr::recode(x,
                synergistic_positive = "Synergistic Positive",
                antagonistic_or_buffered = "Antagonistic or Buffered",
                emergent_only = "Emergent Only",
                not_significant_or_small = "Not significant or small",
                .default = x)
}

run_stage_local_model <- function(meta_all, target_stage, out_dir) {
  message("\n== Stage-local CL153 model: T25 + ", target_stage)

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  meta <- meta_all %>%
    filter(stage %in% c("T25", target_stage)) %>%
    mutate(stage = factor(stage, levels = c("T25", target_stage)),
           water = factor(water, levels = water_levels),
           low_mapping = !is.na(percent_mapped) & percent_mapped < 50)

  balance <- meta %>%
    count(water, stage, name = "n_libraries") %>%
    complete(water = factor(water_levels, levels = water_levels),
             stage = factor(c("T25", target_stage), levels = c("T25", target_stage)),
             fill = list(n_libraries = 0)) %>%
    arrange(stage, water) %>%
    mutate(target_stage = target_stage)
  write_tsv(balance, file.path(out_dir, "sample_balance.tsv"))

  low_map <- meta %>% filter(low_mapping) %>%
    select(sample_id, run_accession, water, stage, replicate, percent_mapped)
  write_tsv(low_map, file.path(out_dir, "low_mapping_in_stage_model.tsv"))

  tx2gene_path <- "data/processed/annotations/tx2gene_canephora.tsv"
  if (!file.exists(tx2gene_path)) stop("Missing tx2gene: ", tx2gene_path)
  tx2gene <- readr::read_tsv(tx2gene_path, show_col_types = FALSE) %>% as.data.frame()
  tx2gene[[1]] <- sub("\\.[0-9]+$", "", tx2gene[[1]])
  tx2gene <- tx2gene[!duplicated(tx2gene[[1]]), , drop = FALSE]

  files <- meta$quant_file
  names(files) <- meta$sample_id

  txi <- tximport(files,
                  type = "salmon",
                  tx2gene = tx2gene,
                  ignoreTxVersion = TRUE,
                  dropInfReps = TRUE,
                  countsFromAbundance = "lengthScaledTPM")

  coldata <- meta %>%
    select(sample_id, run_accession, genotype, water, stage, replicate, index_species, percent_mapped, low_mapping) %>%
    mutate(condition = factor(paste(water, stage, sep = "_"))) %>%
    as.data.frame()
  rownames(coldata) <- coldata$sample_id

  keep <- rowSums(txi$counts >= 5) >= max(2, floor(ncol(txi$counts) * 0.1))
  txi$counts <- txi$counts[keep, , drop = FALSE]
  txi$abundance <- txi$abundance[keep, , drop = FALSE]
  txi$length <- txi$length[keep, , drop = FALSE]
  message("Features kept for stage-local ", target_stage, ": ", sum(keep))

  dds <- DESeqDataSetFromTximport(txi, colData = coldata, design = ~ water + stage + water:stage)
  dds <- DESeq(dds)
  saveRDS(dds, file.path(out_dir, "dds.rds"))

  coef_names <- resultsNames(dds)
  coef <- find_interaction_coef(coef_names, target_stage)
  if (is.na(coef)) stop("Could not find interaction coefficient for ", target_stage, ". Coefficients: ", paste(coef_names, collapse = ", "))

  res <- results(dds, name = coef)
  res_tbl <- as.data.frame(res) %>%
    tibble::rownames_to_column("gene_id") %>%
    mutate(target_stage = target_stage, coefficient = coef)
  write_tsv(res_tbl, file.path(out_dir, paste0("interaction_", target_stage, ".tsv")))

  norm_counts <- counts(dds, normalized = TRUE)
  log_norm <- log2(norm_counts + 1)

  # Manual non-additivity within this two-stage model.
  samples <- list()
  for (w in water_levels) {
    for (st in c("T25", target_stage)) {
      key <- paste(w, st, sep = "_")
      samples[[key]] <- coldata$sample_id[coldata$water == w & coldata$stage == st]
    }
  }
  means <- sapply(names(samples), function(k) safe_group_mean(log_norm, samples[[k]]))
  rownames(means) <- rownames(log_norm)

  baseline <- means[, "WW_T25"]
  drought25 <- means[, "SWD_T25"] - baseline
  heat <- means[, paste("WW", target_stage, sep = "_")] - baseline
  combo <- means[, paste("SWD", target_stage, sep = "_")] - baseline
  nonadd <- combo - heat - drought25
  rec_mem <- if (target_stage == "REC14") means[, paste("SWD", target_stage, sep = "_")] - means[, paste("WW", target_stage, sep = "_")] else rep(NA_real_, length(nonadd))

  nonadd_tbl <- tibble(gene_id = rownames(means),
                       target_stage = target_stage,
                       effect_WW_stage_vs_WW_T25 = heat,
                       effect_SWD_T25_vs_WW_T25 = drought25,
                       effect_SWD_stage_vs_WW_T25 = combo,
                       nonadditive_score = nonadd,
                       recovery_memory_score = rec_mem) %>%
    left_join(res_tbl %>% select(gene_id, log2FoldChange, lfcSE, stat, pvalue, padj, coefficient), by = "gene_id") %>%
    mutate(nonadd_class = classify_nonadd(nonadditive_score, effect_SWD_stage_vs_WW_T25, effect_WW_stage_vs_WW_T25, effect_SWD_T25_vs_WW_T25, padj)) %>%
    arrange(padj, desc(abs(nonadditive_score)))
  write_tsv(nonadd_tbl, file.path(out_dir, paste0("stage_local_nonadditivity_", target_stage, ".tsv")))

  class_summary <- nonadd_tbl %>%
    count(target_stage, nonadd_class, name = "n") %>%
    complete(target_stage = target_stage, nonadd_class = nonadditive_classes, fill = list(n = 0)) %>%
    mutate(nonadd_class = factor(nonadd_class, levels = nonadditive_classes)) %>%
    group_by(target_stage) %>%
    mutate(nonadditive_genes = sum(n[nonadd_class != "not_significant_or_small"]),
           within_nonadditive_pct = ifelse(nonadd_class == "not_significant_or_small", NA_real_, 100 * n / pmax(nonadditive_genes, 1))) %>%
    ungroup() %>%
    mutate(class_label = pretty_class(as.character(nonadd_class)),
           model = paste0("stage_local_", target_stage),
           n_libraries = nrow(meta),
           low_mapping_libraries_in_model = sum(meta$low_mapping, na.rm = TRUE),
           features_kept = sum(keep))
  write_tsv(class_summary, file.path(out_dir, "class_summary.tsv"))

  # Dispersion summary.
  disp_tbl <- tibble(
    gene_id = rownames(dds),
    baseMean = rowMeans(counts(dds, normalized = TRUE)),
    dispGeneEst = mcols(dds)$dispGeneEst,
    dispersion = dispersions(dds),
    dispFit = mcols(dds)$dispFit
  )
  disp_summary <- disp_tbl %>%
    summarise(target_stage = target_stage,
              model = paste0("stage_local_", target_stage),
              n_features = dplyr::n(),
              median_baseMean = median(baseMean, na.rm = TRUE),
              median_dispersion = median(dispersion, na.rm = TRUE),
              q25_dispersion = quantile(dispersion, 0.25, na.rm = TRUE),
              q75_dispersion = quantile(dispersion, 0.75, na.rm = TRUE))
  write_tsv(disp_summary, file.path(out_dir, "dispersion_summary.tsv"))

  list(target_stage = target_stage,
       balance = balance,
       low_mapping = low_map,
       nonadd = nonadd_tbl,
       class_summary = class_summary,
       dispersion_summary = disp_summary)
}

meta <- read_metadata()
mapping <- read_mapping_rates(meta$sample_id)
meta <- meta %>%
  left_join(mapping, by = "sample_id") %>%
  mutate(percent_mapped = as.numeric(percent_mapped),
         low_mapping = !is.na(percent_mapped) & percent_mapped < 50)

write_tsv(meta %>% select(sample_id, run_accession, water, stage, replicate, percent_mapped, low_mapping),
          file.path(out_root, "cl153_samples_for_stage_local_sensitivity.tsv"))

stage_results <- lapply(target_stages, function(st) {
  run_stage_local_model(meta, st, file.path(out_root, paste0("stage_local_", st)))
})

class_all <- bind_rows(lapply(stage_results, `[[`, "class_summary"))
write_tsv(class_all, file.path(out_root, "stage_local_class_summary_all.tsv"))

stage_counts <- class_all %>%
  filter(nonadd_class != "not_significant_or_small") %>%
  group_by(model, target_stage, n_libraries, low_mapping_libraries_in_model, features_kept) %>%
  summarise(stage_local_nonadditive_genes = sum(n), .groups = "drop")
write_tsv(stage_counts, file.path(out_root, "stage_local_nonadditive_gene_counts.tsv"))

balance_all <- bind_rows(lapply(stage_results, `[[`, "balance"))
write_tsv(balance_all, file.path(out_root, "stage_local_sample_balance.tsv"))

low_mapping_all <- bind_rows(lapply(stage_results, function(x) {
  x$low_mapping %>% mutate(target_stage_model = x$target_stage)
}))
write_tsv(low_mapping_all, file.path(out_root, "stage_local_low_mapping_libraries.tsv"))

disp_all <- bind_rows(lapply(stage_results, `[[`, "dispersion_summary"))
write_tsv(disp_all, file.path(out_root, "stage_local_dispersion_summary.tsv"))

# Compare against existing global all-24 and drop-5 outputs if present.
global_count_path <- "results/qc_validation/cl153_low_mapping_exclusion/nonadditive_gene_count_comparison.tsv"
if (file.exists(global_count_path)) {
  global_counts <- readr::read_tsv(global_count_path, show_col_types = FALSE)
  compare <- stage_counts %>%
    select(target_stage, stage_local_nonadditive_genes, n_libraries, low_mapping_libraries_in_model, features_kept) %>%
    left_join(global_counts, by = "target_stage") %>%
    mutate(stage_local_vs_all24_delta = stage_local_nonadditive_genes - all_24_libraries,
           stage_local_vs_all24_ratio = ifelse(all_24_libraries > 0, stage_local_nonadditive_genes / all_24_libraries, NA_real_),
           stage_local_vs_drop5_delta = stage_local_nonadditive_genes - drop_5_libraries_lt50pct,
           stage_local_vs_drop5_ratio = ifelse(drop_5_libraries_lt50pct > 0, stage_local_nonadditive_genes / drop_5_libraries_lt50pct, NA_real_))
  write_tsv(compare, file.path(out_root, "stage_local_vs_global_count_comparison.tsv"))
} else {
  compare <- stage_counts
}

# Concise T37 stage-local summary.
t37_row <- compare %>% filter(target_stage == "T37") %>% dplyr::slice(1)
t37_summary <- tibble(
  item = c(
    "Purpose",
    "T37 stage-local libraries",
    "T37 stage-local low-mapping libraries",
    "Primary all-stage/all-library T37 non-additive genes",
    "Drop-low-mapping all-stage T37 non-additive genes",
    "Stage-local T25+T37 T37 non-additive genes",
    "Interpretation"
  ),
  value = c(
    "Evaluates CL153 T37 non-additivity when the DESeq2 model excludes all T42/REC14 libraries, including the five low-mapping libraries.",
    as.character(t37_row$n_libraries),
    as.character(t37_row$low_mapping_libraries_in_model),
    as.character(t37_row$all_24_libraries),
    as.character(t37_row$drop_5_libraries_lt50pct),
    as.character(t37_row$stage_local_nonadditive_genes),
    "Reports whether the T25+T37 model recovers a large early CL153 interaction response when all T42/REC14 libraries are excluded."
  )
)
write_tsv(t37_summary, file.path(out_root, "stage_local_T37_summary.tsv"))

if (requireNamespace("openxlsx", quietly = TRUE)) {
  wb <- openxlsx::createWorkbook()
  add_sheet <- function(name, dat) {
    nm <- substr(name, 1, 31)
    openxlsx::addWorksheet(wb, nm)
    openxlsx::writeData(wb, nm, dat)
  }
  add_sheet("README", tibble(
    field = c("Purpose", "Key model", "Question addressed"),
    value = c("Stage-local CL153 sensitivity models using T25 + one target stage at a time.",
              "The T37 stage-local model uses only T25 and T37 samples, excluding all T42/REC14 libraries and therefore excluding all five low-mapping libraries.",
              "Whether later-stage low-mapping libraries influence the CL153 T37 interaction estimate.")
  ))
  add_sheet("Stage_local_counts", stage_counts)
  add_sheet("Global_comparison", compare)
  add_sheet("Sample_balance", balance_all)
  add_sheet("Low_mapping_in_models", low_mapping_all)
  add_sheet("Class_summary", class_all)
  add_sheet("Dispersion_summary", disp_all)
  add_sheet("T37_summary", t37_summary)
  openxlsx::saveWorkbook(wb, file.path(out_root, "CL153_stage_local_sensitivity.xlsx"), overwrite = TRUE)
}

message("== DONE CL153 stage-local sensitivity")
message("Outputs written to: ", out_root)
