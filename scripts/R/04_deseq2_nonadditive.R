suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tximport)
  library(DESeq2)
})

message("== Coffee DESeq2 non-additivity analysis")
root <- getwd()
meta_path <- if (file.exists("config/sample_metadata_curated.tsv")) "config/sample_metadata_curated.tsv" else "config/sample_metadata_auto.tsv"
if (!file.exists(meta_path)) stop("Missing metadata: run metadata step first")
meta <- read_tsv(meta_path, show_col_types = FALSE)

required_cols <- c("sample_id","run_accession","index_species","genotype","water","stage","fastq_1","salmon_index","include")
missing_cols <- setdiff(required_cols, names(meta))
if (length(missing_cols)) stop("Metadata missing columns: ", paste(missing_cols, collapse=", "))

meta <- meta %>%
  mutate(include = tolower(include),
         quant_file = file.path("results", "salmon", sample_id, "quant.sf")) %>%
  filter(include %in% c("yes","y","true","1")) %>%
  filter(!if_any(c(index_species, genotype, water, stage), ~ .x == "UNPARSED")) %>%
  filter(file.exists(quant_file))

if (nrow(meta) == 0) stop("No usable samples with quant.sf found. Run Salmon quant first and curate metadata.")

message("Samples with quant.sf: ", nrow(meta))
print(meta %>% dplyr::count(index_species, genotype, water, stage), n = 100)

stage_levels <- c("T25", "T37", "T42", "REC14")
water_levels <- c("WW", "SWD")

tx2gene_for <- function(index_species) {
  if (index_species == "arabica") return("data/processed/annotations/tx2gene_arabica.tsv")
  if (index_species == "canephora") return("data/processed/annotations/tx2gene_canephora.tsv")
  stop("Unknown index_species: ", index_species)
}

find_interaction_coef <- function(coefs, stage) {
  # DESeq2 coefficient names vary by version/order. Find a coefficient containing water, SWD, stage, and target stage.
  hits <- coefs[grepl("water", coefs, ignore.case=TRUE) &
                  grepl("SWD", coefs, ignore.case=TRUE) &
                  grepl("stage", coefs, ignore.case=TRUE) &
                  grepl(stage, coefs, ignore.case=TRUE)]
  if (length(hits) >= 1) return(hits[1])
  return(NA_character_)
}

safe_group_mean <- function(mat, samples) {
  samples <- intersect(samples, colnames(mat))
  if (length(samples) == 0) return(rep(NA_real_, nrow(mat)))
  rowMeans(mat[, samples, drop=FALSE], na.rm=TRUE)
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

for (grp in unique(meta$index_species)) {
  m <- meta %>% filter(index_species == grp)
  # Keep genotype in the group name. Usually arabica=Icatu and canephora=CL153.
  gname <- paste(unique(m$genotype), collapse="_")
  analysis_group <- paste(grp, gname, sep="_") %>% str_replace_all("[^A-Za-z0-9_.-]+", "_")
  message("\n== Analysis group: ", analysis_group)

  condition_table <- m %>% dplyr::count(water, stage) %>% arrange(water, stage)
  print(condition_table, n = 100)

  needed <- expand.grid(water=water_levels, stage=stage_levels, stringsAsFactors = FALSE)
  available <- condition_table %>% mutate(key = paste(water, stage, sep="_")) %>% pull(key)
  missing_conditions <- setdiff(paste(needed$water, needed$stage, sep="_"), available)
  if (length(missing_conditions) > 0) {
    message("WARNING: missing conditions in ", analysis_group, ": ", paste(missing_conditions, collapse=", "))
    message("Continuing with available conditions; some non-additivity terms may be skipped.")
  }

  tx2gene_path <- tx2gene_for(grp)
  if (!file.exists(tx2gene_path)) stop("Missing tx2gene: ", tx2gene_path, ". Run index step first.")
  tx2gene <- readr::read_tsv(tx2gene_path, show_col_types = FALSE)
  tx2gene <- as.data.frame(tx2gene)
  tx2gene[[1]] <- sub("\\.[0-9]+$", "", tx2gene[[1]])
  tx2gene <- tx2gene[!duplicated(tx2gene[[1]]), , drop = FALSE]

  files <- m$quant_file
  names(files) <- m$sample_id

  txi <- tximport(
    files,
    type = "salmon",
    tx2gene = tx2gene,
    ignoreTxVersion = TRUE,
    dropInfReps = TRUE,
    countsFromAbundance = "lengthScaledTPM"
  )
  coldata <- m %>%
    select(sample_id, run_accession, genotype, water, stage, replicate, index_species) %>%
    mutate(water = factor(water, levels=water_levels),
           stage = factor(stage, levels=stage_levels),
           condition = factor(paste(water, stage, sep="_"))) %>%
    as.data.frame()
  rownames(coldata) <- coldata$sample_id

  # Remove genes with almost no counts.
  keep <- rowSums(txi$counts >= 5) >= max(2, floor(ncol(txi$counts) * 0.1))
  txi$counts <- txi$counts[keep, , drop=FALSE]
  txi$abundance <- txi$abundance[keep, , drop=FALSE]
  txi$length <- txi$length[keep, , drop=FALSE]
  message("Genes/transcript-groups kept: ", sum(keep))

  dds <- DESeqDataSetFromTximport(txi, colData=coldata, design= ~ water + stage + water:stage)
  dds <- DESeq(dds)

  outdir <- file.path("results", "deseq2", analysis_group)
  dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
  dir.create("results/tables", recursive=TRUE, showWarnings=FALSE)

  write_tsv(tibble(coefficient = resultsNames(dds)), file.path(outdir, "deseq2_coefficients.tsv"))

  norm_counts <- counts(dds, normalized=TRUE)
  write_tsv(as.data.frame(norm_counts) %>% tibble::rownames_to_column("gene_id"), file.path(outdir, "normalized_counts.tsv"))
  # Variance-stabilized matrix retained for downstream inspection.
  vsd <- vst(dds, blind=FALSE)

  # Interaction coefficients from DESeq2.
  coef_names <- resultsNames(dds)
  interaction_tables <- list()
  for (st in c("T37", "T42", "REC14")) {
    coef <- find_interaction_coef(coef_names, st)
    if (is.na(coef)) {
      message("No interaction coefficient found for ", st, " in ", analysis_group)
      next
    }
    res <- results(dds, name=coef)
    tbl <- as.data.frame(res) %>%
      tibble::rownames_to_column("gene_id") %>%
      mutate(target_stage = st, coefficient = coef)
    write_tsv(tbl, file.path(outdir, paste0("interaction_", st, ".tsv")))
    interaction_tables[[st]] <- tbl
  }

  # Manual stress arithmetic on normalized log counts.
  log_norm <- log2(norm_counts + 1)
  sample_sets <- list()
  for (w in water_levels) {
    for (st in stage_levels) {
      key <- paste(w, st, sep="_")
      sample_sets[[key]] <- coldata$sample_id[coldata$water == w & coldata$stage == st]
    }
  }
  means <- sapply(names(sample_sets), function(k) safe_group_mean(log_norm, sample_sets[[k]]))
  rownames(means) <- rownames(log_norm)

  nonadd_all <- list()
  baseline <- means[, "WW_T25"]
  drought25 <- means[, "SWD_T25"] - baseline
  for (st in c("T37", "T42", "REC14")) {
    ww_key <- paste("WW", st, sep="_")
    swd_key <- paste("SWD", st, sep="_")
    if (!(ww_key %in% colnames(means)) || !(swd_key %in% colnames(means))) next

    heat_or_stage <- means[, ww_key] - baseline
    combo <- means[, swd_key] - baseline
    nonadd <- combo - heat_or_stage - drought25
    rec_mem <- if (st == "REC14") means[, swd_key] - means[, ww_key] else rep(NA_real_, length(nonadd))

    tbl <- tibble(
      gene_id = rownames(means),
      target_stage = st,
      effect_WW_stage_vs_WW_T25 = heat_or_stage,
      effect_SWD_T25_vs_WW_T25 = drought25,
      effect_SWD_stage_vs_WW_T25 = combo,
      nonadditive_score = nonadd,
      recovery_memory_score = rec_mem
    )
    if (!is.null(interaction_tables[[st]])) {
      tbl <- tbl %>% left_join(interaction_tables[[st]] %>% select(gene_id, log2FoldChange, lfcSE, stat, pvalue, padj, coefficient), by="gene_id")
    } else {
      tbl <- tbl %>% mutate(log2FoldChange=NA_real_, lfcSE=NA_real_, stat=NA_real_, pvalue=NA_real_, padj=NA_real_, coefficient=NA_character_)
    }
    tbl <- tbl %>% mutate(nonadd_class = classify_nonadd(nonadditive_score, effect_SWD_stage_vs_WW_T25, effect_WW_stage_vs_WW_T25, effect_SWD_T25_vs_WW_T25, padj)) %>%
      arrange(padj, desc(abs(nonadditive_score)))
    nonadd_all[[st]] <- tbl
    write_tsv(tbl, file.path(outdir, paste0("nonadditivity_", st, ".tsv")))
  }

  if (length(nonadd_all)) {
    all_tbl <- bind_rows(nonadd_all)
    write_tsv(all_tbl, file.path("results/tables", paste0(analysis_group, "_nonadditivity_scores.tsv")))

    top <- all_tbl %>% filter(!is.na(padj), padj < 0.05) %>% arrange(padj, desc(abs(nonadditive_score))) %>% head(50)
    write_tsv(top, file.path("results/tables", paste0(analysis_group, "_top50_interaction_candidates.tsv")))

    summary_tbl <- all_tbl %>% dplyr::count(target_stage, nonadd_class) %>% arrange(target_stage, desc(n))
    write_tsv(summary_tbl, file.path("results/tables", paste0(analysis_group, "_nonadditivity_class_summary.tsv")))
  }

  saveRDS(dds, file.path(outdir, "dds.rds"))
  saveRDS(vsd, file.path(outdir, "vst.rds"))
  message("== Finished ", analysis_group)
}

message("== All done")
