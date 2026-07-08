#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tximport)
  library(DESeq2)
})

# Rebind tidyverse verbs after Bioconductor packages are loaded.
# DESeq2/AnnotationDbi can mask dplyr verbs such as count() and select().
# Keeping these aliases makes the script robust across R/Bioconductor versions.
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
complete <- tidyr::complete
pivot_wider <- tidyr::pivot_wider
pivot_longer <- tidyr::pivot_longer


message("== CL153 low-mapping library exclusion analysis")
root <- getwd()
out_root <- "results/qc_validation/cl153_low_mapping_exclusion"
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

stage_levels <- c("T25", "T37", "T42", "REC14")
water_levels <- c("WW", "SWD")

# -----------------------------
# Utility functions
# -----------------------------
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
    mutate(water = factor(water, levels = water_levels),
           stage = factor(stage, levels = stage_levels),
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
  csv_path <- "results/tables/salmon_mapping_rates.csv"
  tsv_path <- "results/analysis_tables/qc/sample_qc_with_metadata.tsv"
  if (file.exists(tsv_path)) {
    qc <- readr::read_tsv(tsv_path, show_col_types = FALSE)
    if (all(c("sample_id", "percent_mapped") %in% names(qc))) {
      return(qc %>% select(sample_id, num_processed, num_mapped, percent_mapped))
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

nonadditive_classes <- c("synergistic_positive", "antagonistic_or_buffered", "emergent_only", "not_significant_or_small")

run_model <- function(meta, label, out_dir) {
  message("\n== Running CL153 model: ", label, "  n=", nrow(meta))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  condition_counts <- meta %>%
    count(water, stage, name = "n_libraries") %>%
    complete(water = factor(water_levels, levels = water_levels),
             stage = factor(stage_levels, levels = stage_levels),
             fill = list(n_libraries = 0)) %>%
    arrange(water, stage) %>%
    mutate(model = label)
  write_tsv(condition_counts, file.path(out_dir, "sample_balance.tsv"))

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
    mutate(water = factor(water, levels = water_levels),
           stage = factor(stage, levels = stage_levels),
           condition = factor(paste(water, stage, sep = "_"))) %>%
    as.data.frame()
  rownames(coldata) <- coldata$sample_id

  keep <- rowSums(txi$counts >= 5) >= max(2, floor(ncol(txi$counts) * 0.1))
  txi$counts <- txi$counts[keep, , drop = FALSE]
  txi$abundance <- txi$abundance[keep, , drop = FALSE]
  txi$length <- txi$length[keep, , drop = FALSE]
  message("Genes/transcript-groups kept for ", label, ": ", sum(keep))

  dds <- DESeqDataSetFromTximport(txi, colData = coldata, design = ~ water + stage + water:stage)
  dds <- DESeq(dds)
  vsd <- vst(dds, blind = FALSE)

  write_tsv(tibble(coefficient = resultsNames(dds)), file.path(out_dir, "deseq2_coefficients.tsv"))
  norm_counts <- counts(dds, normalized = TRUE)
  write_tsv(as.data.frame(norm_counts) %>% tibble::rownames_to_column("gene_id"), file.path(out_dir, "normalized_counts.tsv"))
  saveRDS(dds, file.path(out_dir, "dds.rds"))
  saveRDS(vsd, file.path(out_dir, "vst.rds"))

  # Dispersion summary for evaluating model stability.
  disp_tbl <- tibble(
    gene_id = rownames(dds),
    baseMean = rowMeans(counts(dds, normalized = TRUE)),
    dispGeneEst = mcols(dds)$dispGeneEst,
    dispersion = dispersions(dds),
    dispFit = mcols(dds)$dispFit
  )
  write_tsv(disp_tbl, file.path(out_dir, "dispersion_estimates.tsv"))
  disp_summary <- disp_tbl %>%
    summarise(model = label,
              n_features = dplyr::n(),
              median_baseMean = median(baseMean, na.rm = TRUE),
              median_dispersion = median(dispersion, na.rm = TRUE),
              q25_dispersion = quantile(dispersion, 0.25, na.rm = TRUE),
              q75_dispersion = quantile(dispersion, 0.75, na.rm = TRUE),
              median_dispGeneEst = median(dispGeneEst, na.rm = TRUE))
  write_tsv(disp_summary, file.path(out_dir, "dispersion_summary.tsv"))

  # Interaction results.
  coef_names <- resultsNames(dds)
  interaction_tables <- list()
  for (st in c("T37", "T42", "REC14")) {
    coef <- find_interaction_coef(coef_names, st)
    if (is.na(coef)) {
      message("No interaction coefficient found for ", st, " in ", label)
      next
    }
    res <- results(dds, name = coef)
    tbl <- as.data.frame(res) %>%
      tibble::rownames_to_column("gene_id") %>%
      mutate(model = label, target_stage = st, coefficient = coef)
    write_tsv(tbl, file.path(out_dir, paste0("interaction_", st, ".tsv")))
    interaction_tables[[st]] <- tbl
  }

  # Manual stress arithmetic.
  log_norm <- log2(norm_counts + 1)
  sample_sets <- list()
  for (w in water_levels) {
    for (st in stage_levels) {
      key <- paste(w, st, sep = "_")
      sample_sets[[key]] <- coldata$sample_id[coldata$water == w & coldata$stage == st]
    }
  }
  means <- sapply(names(sample_sets), function(k) safe_group_mean(log_norm, sample_sets[[k]]))
  rownames(means) <- rownames(log_norm)

  nonadd_all <- list()
  baseline <- means[, "WW_T25"]
  drought25 <- means[, "SWD_T25"] - baseline
  for (st in c("T37", "T42", "REC14")) {
    ww_key <- paste("WW", st, sep = "_")
    swd_key <- paste("SWD", st, sep = "_")
    if (!(ww_key %in% colnames(means)) || !(swd_key %in% colnames(means))) next
    heat_or_stage <- means[, ww_key] - baseline
    combo <- means[, swd_key] - baseline
    nonadd <- combo - heat_or_stage - drought25
    rec_mem <- if (st == "REC14") means[, swd_key] - means[, ww_key] else rep(NA_real_, length(nonadd))
    tbl <- tibble(gene_id = rownames(means),
                  model = label,
                  target_stage = st,
                  effect_WW_stage_vs_WW_T25 = heat_or_stage,
                  effect_SWD_T25_vs_WW_T25 = drought25,
                  effect_SWD_stage_vs_WW_T25 = combo,
                  nonadditive_score = nonadd,
                  recovery_memory_score = rec_mem)
    if (!is.null(interaction_tables[[st]])) {
      tbl <- tbl %>% left_join(interaction_tables[[st]] %>% select(gene_id, log2FoldChange, lfcSE, stat, pvalue, padj, coefficient), by = "gene_id")
    } else {
      tbl <- tbl %>% mutate(log2FoldChange = NA_real_, lfcSE = NA_real_, stat = NA_real_, pvalue = NA_real_, padj = NA_real_, coefficient = NA_character_)
    }
    tbl <- tbl %>%
      mutate(nonadd_class = classify_nonadd(nonadditive_score, effect_SWD_stage_vs_WW_T25, effect_WW_stage_vs_WW_T25, effect_SWD_T25_vs_WW_T25, padj)) %>%
      arrange(padj, desc(abs(nonadditive_score)))
    nonadd_all[[st]] <- tbl
    write_tsv(tbl, file.path(out_dir, paste0("nonadditivity_", st, ".tsv")))
  }

  all_tbl <- bind_rows(nonadd_all)
  write_tsv(all_tbl, file.path(out_dir, "nonadditivity_all_stages.tsv"))

  # PCA coordinates.
  mat <- assay(vsd)
  vars <- matrixStats::rowVars(mat)
  keep_pca <- order(vars, decreasing = TRUE)[seq_len(min(500, length(vars)))]
  pcs <- prcomp(t(mat[keep_pca, , drop = FALSE]), scale. = FALSE)
  pca_df <- as.data.frame(pcs$x[, 1:2]) %>%
    tibble::rownames_to_column("sample_id") %>%
    left_join(coldata %>% tibble::rownames_to_column("sample_id0") %>% select(-sample_id0), by = "sample_id") %>%
    mutate(model = label,
           PC1_percent = 100 * summary(pcs)$importance[2, 1],
           PC2_percent = 100 * summary(pcs)$importance[2, 2])
  write_tsv(pca_df, file.path(out_dir, "pca_coordinates.tsv"))

  class_summary <- all_tbl %>%
    mutate(nonadd_class = factor(nonadd_class, levels = nonadditive_classes)) %>%
    count(model, target_stage, nonadd_class, name = "n") %>%
    complete(model, target_stage = c("T37", "T42", "REC14"), nonadd_class = nonadditive_classes, fill = list(n = 0)) %>%
    group_by(model, target_stage) %>%
    mutate(stage_total_all_classes = sum(n),
           nonadditive_genes = sum(n[nonadd_class != "not_significant_or_small"]),
           within_nonadditive_pct = ifelse(nonadd_class == "not_significant_or_small", NA_real_, 100 * n / pmax(nonadditive_genes, 1))) %>%
    ungroup() %>%
    mutate(class_label = pretty_class(as.character(nonadd_class)))
  write_tsv(class_summary, file.path(out_dir, "nonadditivity_class_summary.tsv"))

  # Return compact object.
  list(label = label,
       out_dir = out_dir,
       dds = dds,
       vsd = vsd,
       sample_balance = condition_counts,
       dispersion_summary = disp_summary,
       nonadd = all_tbl,
       class_summary = class_summary,
       pca = pca_df)
}

# -----------------------------
# Run sensitivity models
# -----------------------------
meta <- read_metadata()
if (nrow(meta) != 24) {
  message("WARNING: expected 24 CL153 libraries with quant.sf, found ", nrow(meta), ". Continuing with available libraries.")
}

mapping <- read_mapping_rates(meta$sample_id)
meta <- meta %>%
  left_join(mapping, by = "sample_id") %>%
  mutate(percent_mapped = as.numeric(percent_mapped),
         low_mapping = !is.na(percent_mapped) & percent_mapped < 50)

low_map_tbl <- meta %>%
  filter(low_mapping) %>%
  select(sample_id, run_accession, genotype, water, stage, replicate, num_processed, num_mapped, percent_mapped) %>%
  arrange(percent_mapped)
write_tsv(low_map_tbl, file.path(out_root, "low_mapping_libraries_lt50.tsv"))

sample_qc <- meta %>%
  select(sample_id, run_accession, genotype, water, stage, replicate, num_processed, num_mapped, percent_mapped, low_mapping, quant_file) %>%
  arrange(stage, water, replicate)
write_tsv(sample_qc, file.path(out_root, "cl153_sample_qc_for_sensitivity.tsv"))

model_all <- run_model(meta, "all_24_libraries", file.path(out_root, "model_all_24_libraries"))
meta_drop <- meta %>% filter(!low_mapping)
model_drop <- run_model(meta_drop, "drop_5_libraries_lt50pct", file.path(out_root, "model_drop_5_libraries_lt50pct"))

# -----------------------------
# Compare model outputs
# -----------------------------
class_all <- bind_rows(model_all$class_summary, model_drop$class_summary)
write_tsv(class_all, file.path(out_root, "class_summary_all_vs_drop.tsv"))

nonadd_counts <- class_all %>%
  filter(nonadd_class != "not_significant_or_small") %>%
  group_by(model, target_stage) %>%
  summarise(nonadditive_genes = sum(n), .groups = "drop") %>%
  complete(model = c("all_24_libraries", "drop_5_libraries_lt50pct"),
           target_stage = c("T37", "T42", "REC14"),
           fill = list(nonadditive_genes = 0))

count_comparison <- nonadd_counts %>%
  pivot_wider(names_from = model, values_from = nonadditive_genes) %>%
  mutate(delta_drop_minus_all = drop_5_libraries_lt50pct - all_24_libraries,
         ratio_drop_over_all = ifelse(all_24_libraries > 0, drop_5_libraries_lt50pct / all_24_libraries, NA_real_))
write_tsv(count_comparison, file.path(out_root, "nonadditive_gene_count_comparison.tsv"))

# Overlap/Jaccard for non-additive gene sets by stage.
overlap_tbl <- lapply(c("T37", "T42", "REC14"), function(st) {
  a <- model_all$nonadd %>% filter(target_stage == st, nonadd_class != "not_significant_or_small") %>% pull(gene_id) %>% unique()
  b <- model_drop$nonadd %>% filter(target_stage == st, nonadd_class != "not_significant_or_small") %>% pull(gene_id) %>% unique()
  tibble(target_stage = st,
         all_24_n = length(a),
         drop_5_n = length(b),
         overlap_n = length(intersect(a, b)),
         union_n = length(union(a, b)),
         jaccard = ifelse(length(union(a, b)) > 0, length(intersect(a, b)) / length(union(a, b)), NA_real_))
}) %>% bind_rows()
write_tsv(overlap_tbl, file.path(out_root, "nonadditive_gene_set_overlap.tsv"))

balance_all <- bind_rows(model_all$sample_balance, model_drop$sample_balance)
write_tsv(balance_all, file.path(out_root, "sample_balance_all_vs_drop.tsv"))

disp_all <- bind_rows(model_all$dispersion_summary, model_drop$dispersion_summary)
write_tsv(disp_all, file.path(out_root, "dispersion_summary_all_vs_drop.tsv"))

pca_all <- bind_rows(model_all$pca, model_drop$pca)
write_tsv(pca_all, file.path(out_root, "pca_coordinates_all_vs_drop.tsv"))

# Concise QC summary table.
qc_summary <- tibble(
  item = c("Low-mapping libraries flagged", "Primary CL153 libraries", "Sensitivity CL153 libraries after exclusion", "Primary T37 non-additive genes", "Drop-low-mapping T37 non-additive genes", "Primary T42 non-additive genes", "Drop-low-mapping T42 non-additive genes", "Primary REC14 non-additive genes", "Drop-low-mapping REC14 non-additive genes"),
  value = c(nrow(low_map_tbl), nrow(meta), nrow(meta_drop),
            count_comparison$all_24_libraries[count_comparison$target_stage == "T37"],
            count_comparison$drop_5_libraries_lt50pct[count_comparison$target_stage == "T37"],
            count_comparison$all_24_libraries[count_comparison$target_stage == "T42"],
            count_comparison$drop_5_libraries_lt50pct[count_comparison$target_stage == "T42"],
            count_comparison$all_24_libraries[count_comparison$target_stage == "REC14"],
            count_comparison$drop_5_libraries_lt50pct[count_comparison$target_stage == "REC14"])
)
write_tsv(qc_summary, file.path(out_root, "cl153_low_mapping_summary.tsv"))

# -----------------------------
# Workbook
# -----------------------------
if (requireNamespace("openxlsx", quietly = TRUE)) {
  wb <- openxlsx::createWorkbook()
  add_sheet <- function(name, dat) {
    openxlsx::addWorksheet(wb, substr(name, 1, 31))
    openxlsx::writeData(wb, substr(name, 1, 31), dat)
  }
  add_sheet("README", tibble(
    field = c("Purpose", "Primary model", "Sensitivity model", "Question addressed"),
    value = c("CL153 low-mapping library exclusion sensitivity analysis",
              "All 24 CL153 libraries retained",
              "Five CL153 libraries with Salmon mapping rate <50% excluded",
              "Whether low-mapping T42/REC14 libraries altered dispersion estimates and non-additive gene counts")
  ))
  add_sheet("Low_mapping_libraries", low_map_tbl)
  add_sheet("Sample_balance", balance_all)
  add_sheet("Mapping_QC", sample_qc)
  add_sheet("Nonadd_count_comparison", count_comparison)
  add_sheet("Class_summary", class_all)
  add_sheet("Gene_set_overlap", overlap_tbl)
  add_sheet("Dispersion_summary", disp_all)
  add_sheet("PCA_coordinates", pca_all)
  add_sheet("QC_summary", qc_summary)
  openxlsx::saveWorkbook(wb, file.path(out_root, "CL153_low_mapping_sensitivity.xlsx"), overwrite = TRUE)
}


message("== DONE CL153 low-mapping sensitivity")
message("Outputs written to: ", out_root)
