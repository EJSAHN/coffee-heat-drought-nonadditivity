library(readr)
library(dplyr)
library(stringr)

base <- "results/paper_pack_v1/tables"
qcbase <- "results/paper_pack_v1/qc"
outdir <- "results/paper_pack_v1"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

classes <- read_tsv(file.path(base, "paper_nonadditivity_class_numbers.tsv"), show_col_types = FALSE)
enrich  <- read_tsv(file.path(base, "paper_module_enrichment_ranked.tsv"), show_col_types = FALSE)
cand    <- read_tsv(file.path(base, "paper_interpretable_candidate_panel.tsv"), show_col_types = FALSE)
topmod  <- read_tsv(file.path(base, "paper_top5_candidates_per_module.tsv"), show_col_types = FALSE)
qc      <- read_tsv(file.path(qcbase, "sample_qc_with_metadata.tsv"), show_col_types = FALSE)

clean_product <- function(x) {
  x %>%
    str_replace_all("%2C", ",") %>%
    str_replace_all("%3B", ";") %>%
    str_replace_all("Source:Projected from ", "projected from ")
}

bad_product_regex <- "uncharacterized|hypothetical|unknown|predicted:.*uncharacterized|rrna|ribosomal rna|ncrna|snoRNA|cylicin"

hc <- cand %>%
  mutate(
    product_clean = clean_product(product_clean),
    product_lower = str_to_lower(product_clean),
    high_conf_product = !str_detect(product_lower, bad_product_regex),
    abs_nonadd = abs(nonadditive_score)
  ) %>%
  filter(
    high_conf_product,
    functional_module != "other_or_unannotated",
    !is.na(padj),
    padj < 0.05,
    abs_nonadd >= 1
  ) %>%
  arrange(group, target_stage, functional_module, desc(priority_score))

hc_top <- hc %>%
  group_by(group, target_stage, functional_module) %>%
  slice_max(priority_score, n = 3, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(group, target_stage, functional_module, desc(priority_score))

stage_summary <- classes %>%
  filter(nonadd_class != "not_significant_or_small") %>%
  group_by(group, target_stage) %>%
  summarise(nonadditive_genes = sum(n), .groups = "drop") %>%
  arrange(group, target_stage)

class_wide <- classes %>%
  filter(nonadd_class != "not_significant_or_small") %>%
  select(group, target_stage, nonadd_class, n, pct) %>%
  arrange(group, target_stage, nonadd_class)

enrich_sig <- enrich %>%
  filter(padj < 0.10) %>%
  arrange(padj)

qc_low <- qc %>%
  filter(percent_mapped < 50) %>%
  arrange(percent_mapped)

write_tsv(hc, file.path(base, "paper_high_confidence_candidates.tsv"))
write_tsv(hc_top, file.path(base, "paper_high_confidence_top3_per_module.tsv"))
write_tsv(stage_summary, file.path(base, "paper_stage_nonadditive_totals.tsv"))
write_tsv(class_wide, file.path(base, "paper_class_summary_long.tsv"))
write_tsv(enrich_sig, file.path(base, "paper_significant_or_suggestive_modules.tsv"))
write_tsv(qc_low, file.path(qcbase, "paper_low_mapping_samples_lt50.tsv"))

md <- c(
"# Coffee multistress paper digest v1",
"",
"## Main message",
"Coffee heat-drought responses are non-additive and genotype-specific. Icatu shows a strong T37 interaction burst, whereas CL153 shows delayed T42/REC14 non-additivity and a recovery-stage photosynthesis signal.",
"",
"## Non-additive gene totals",
capture.output(print(stage_summary)),
"",
"## Non-additive class counts",
capture.output(print(class_wide)),
"",
"## Significant or suggestive functional modules",
capture.output(print(enrich_sig)),
"",
"## High-confidence candidate examples",
capture.output(print(hc_top %>% select(group, target_stage, functional_module, nonadd_class, gene_id, product_clean, nonadditive_score, padj, priority_score), n = 100)),
"",
"## QC caveat: low-mapping samples",
capture.output(print(qc_low %>% select(sample_id, genotype, water, stage, replicate, num_processed, num_mapped, percent_mapped), n = 100)),
"",
"## Recommended manuscript framing",
"- Primary claim: genotype-specific non-additive stress arithmetic.",
"- Strongest numeric pattern: Icatu T37 burst vs CL153 T42/REC14 delayed response.",
"- Strongest functional enrichment: CL153 REC14 photosynthesis.",
"- Candidate-level Icatu story: sugar transport, thylakoid assembly, DnaJ/chaperone, ATG7/ATG8, VPS60.2, E3/ERAD ubiquitin pathway.",
"- QC note: CL153 low-mapping samples require cautious interpretation and PCA-supported retention."
)

writeLines(md, file.path(outdir, "paper_digest_v1.md"))

message("DONE paper digest v1")