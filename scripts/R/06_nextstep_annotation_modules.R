Sys.setenv(VROOM_CONNECTION_SIZE = "52428800")
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
})

message('== Next-step annotation and module summaries ==')

dir.create('results/next/annotations', recursive = TRUE, showWarnings = FALSE)
dir.create('results/next/tables', recursive = TRUE, showWarnings = FALSE)

strip_ver <- function(x) sub('\\.[0-9]+$', '', as.character(x))
clean_id <- function(x) {
  x <- as.character(x)
  x <- sub('^transcript:', '', x)
  x <- sub('^gene:', '', x)
  x <- sub('^rna-', '', x)
  x <- sub('^cds-', '', x)
  x
}

parse_attrs_one <- function(attr) {
  parts <- strsplit(attr, ';', fixed = TRUE)[[1]]
  out <- list()
  for (p in parts) {
    if (!nzchar(p)) next
    kv <- strsplit(p, '=', fixed = TRUE)[[1]]
    if (length(kv) >= 2) out[[kv[1]]] <- paste(kv[-1], collapse='=')
  }
  out
}

parse_gff_annotation <- function(gff_path, label) {
  if (!file.exists(gff_path)) {
    message('Missing GFF for ', label, ': ', gff_path)
    return(tibble())
  }
  message('Reading GFF: ', gff_path)
  gff <- readr::read_tsv(
    gff_path,
    comment = '#',
    col_names = c('seqid','source','type','start','end','score','strand','phase','attributes'),
    col_types = cols(.default = col_character()),
    progress = FALSE
  )
  gff <- gff %>% filter(type %in% c('mRNA','transcript','gene'))
  if (nrow(gff) == 0) return(tibble())

  rows <- lapply(seq_len(nrow(gff)), function(i) {
    a <- parse_attrs_one(gff$attributes[i])
    id <- a[['ID']] %||% ''
    name <- a[['Name']] %||% ''
    gene <- a[['gene']] %||% a[['gene_id']] %||% a[['Parent']] %||% ''
    product <- a[['product']] %||% a[['description']] %||% a[['Note']] %||% ''
    dbxref <- a[['Dbxref']] %||% ''

    # Pull GenBank accession if present, e.g. Dbxref=GenBank:XM_027204879.1,GeneID:...
    genbank <- ''
    m <- stringr::str_match(dbxref, 'GenBank:([^,;]+)')
    if (!is.na(m[1,2])) genbank <- m[1,2]

    tibble(
      feature_type = gff$type[i],
      feature_id = clean_id(id),
      feature_id_nover = strip_ver(clean_id(id)),
      name = clean_id(name),
      name_nover = strip_ver(clean_id(name)),
      genbank = clean_id(genbank),
      genbank_nover = strip_ver(clean_id(genbank)),
      gene_symbol = gene,
      product = product,
      dbxref = dbxref
    )
  })
  ann <- bind_rows(rows)

  keys <- bind_rows(
    ann %>% transmute(gene_id_key = feature_id_nover, feature_id, gene_symbol, product, dbxref, source = paste0(label, '_gff_feature')),
    ann %>% transmute(gene_id_key = name_nover,    feature_id, gene_symbol, product, dbxref, source = paste0(label, '_gff_name')),
    ann %>% transmute(gene_id_key = genbank_nover, feature_id, gene_symbol, product, dbxref, source = paste0(label, '_gff_genbank'))
  ) %>%
    filter(!is.na(gene_id_key), gene_id_key != '') %>%
    group_by(gene_id_key) %>%
    summarise(
      feature_id = first(feature_id[feature_id != ''], default = first(feature_id)),
      gene_symbol = first(gene_symbol[gene_symbol != ''], default = ''),
      product = first(product[product != ''], default = ''),
      dbxref = first(dbxref[dbxref != ''], default = ''),
      source = paste(unique(source), collapse=';'),
      .groups = 'drop'
    )
  keys
}

parse_fasta_desc <- function(fasta_path, label) {
  if (!file.exists(fasta_path)) {
    message('Missing FASTA for ', label, ': ', fasta_path)
    return(tibble())
  }
  message('Reading FASTA headers: ', fasta_path)
  con <- if (grepl('\\.gz$', fasta_path)) gzfile(fasta_path, 'rt') else file(fasta_path, 'rt')
  on.exit(close(con), add = TRUE)
  headers <- character()
  while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
    if (startsWith(line, '>')) headers <- c(headers, substring(line, 2))
  }
  if (!length(headers)) return(tibble())
  id <- sub('\\s.*$', '', headers)
  desc <- sub('^[^ ]+\\s*', '', headers)
  tibble(
    gene_id_key = strip_ver(clean_id(id)),
    feature_id = clean_id(id),
    gene_symbol = '',
    product = desc,
    dbxref = '',
    source = paste0(label, '_fasta')
  ) %>% filter(gene_id_key != '') %>% distinct(gene_id_key, .keep_all = TRUE)
}

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

arabica_gff <- 'data/ref/arabica_Cara_1_0_NCBI/GCF_003713225.1_Cara_1.0_genomic.gff.gz'
arabica_rna <- 'data/ref/arabica_Cara_1_0_NCBI/GCF_003713225.1_Cara_1.0_rna.fna.gz'
cane_gff <- 'data/ref/canephora_AUK_PRJEB4211_v1_EnsemblPlants/Coffea_canephora.AUK_PRJEB4211_v1.62.gff3.gz'
cane_pep <- 'data/ref/canephora_AUK_PRJEB4211_v1_EnsemblPlants/Coffea_canephora.AUK_PRJEB4211_v1.pep.all.fa.gz'

ann_arabica <- bind_rows(
  parse_gff_annotation(arabica_gff, 'arabica'),
  parse_fasta_desc(arabica_rna, 'arabica_rna')
) %>% group_by(gene_id_key) %>% summarise(
  feature_id = first(feature_id[feature_id != ''], default = first(feature_id)),
  gene_symbol = first(gene_symbol[gene_symbol != ''], default = ''),
  product = first(product[product != ''], default = ''),
  dbxref = first(dbxref[dbxref != ''], default = ''),
  source = paste(unique(source), collapse=';'),
  .groups = 'drop'
)
ann_cane <- bind_rows(
  parse_gff_annotation(cane_gff, 'canephora'),
  parse_fasta_desc(cane_pep, 'canephora_pep')
) %>% group_by(gene_id_key) %>% summarise(
  feature_id = first(feature_id[feature_id != ''], default = first(feature_id)),
  gene_symbol = first(gene_symbol[gene_symbol != ''], default = ''),
  product = first(product[product != ''], default = ''),
  dbxref = first(dbxref[dbxref != ''], default = ''),
  source = paste(unique(source), collapse=';'),
  .groups = 'drop'
)

write_tsv(ann_arabica, 'results/next/annotations/annotation_arabica.tsv')
write_tsv(ann_cane, 'results/next/annotations/annotation_canephora.tsv')
message('Arabica annotation rows: ', nrow(ann_arabica))
message('Canephora annotation rows: ', nrow(ann_cane))

assign_module <- function(product) {
  p <- tolower(ifelse(is.na(product), '', product))
  case_when(
    str_detect(p, 'photosystem|chlorophyll|rubisco|light-harvesting|lhc|thylakoid|ferredoxin|plastocyanin') ~ 'photosynthesis',
    str_detect(p, 'heat shock|hsp|chaperon|dna[jk]|dnak|clp|groel') ~ 'proteostasis_HSP_chaperone',
    str_detect(p, 'dehydrin|late embryogenesis|lea|osmoprotect|galactinol|raffinose') ~ 'dehydration_LEA_osmoprotection',
    str_detect(p, 'aquaporin|nodulin|pip|tip') ~ 'water_transport_aquaporin',
    str_detect(p, 'abscisic|aba|snrk|pp2c|bzip|areb|abi') ~ 'ABA_signaling',
    str_detect(p, 'peroxidase|superoxide|glutathione|thioredoxin|ascorbate|catalase|oxidoreductase|ros') ~ 'ROS_redox',
    str_detect(p, 'lipid|fatty acid|lipoxygenase|phospholipase|desaturase|wax|cutin') ~ 'lipid_membrane',
    str_detect(p, 'protease|aspartic|peptidase|vacuol|ubiquitin|autophagy|proteasome|saposin') ~ 'protease_vacuole_ubiquitin',
    str_detect(p, 'wrky|nac|myb|bhlh|erf|ethylene|jasmon|salicy|auxin|gibberellin|transcription factor') ~ 'TF_hormone_crosstalk',
    TRUE ~ 'other_or_unannotated'
  )
}

process_group <- function(group_name, ann) {
  score_path <- paste0('results/tables/', group_name, '_nonadditivity_scores.tsv')
  top_path <- paste0('results/tables/', group_name, '_top50_interaction_candidates.tsv')
  class_path <- paste0('results/tables/', group_name, '_nonadditivity_class_summary.tsv')
  if (!file.exists(score_path)) {
    message('Missing scores: ', score_path)
    return(NULL)
  }
  scores <- read_tsv(score_path, show_col_types = FALSE) %>%
    mutate(gene_id_key = strip_ver(gene_id)) %>%
    left_join(ann, by = 'gene_id_key') %>%
    mutate(product = if_else(is.na(product) | product == '', 'unannotated', product),
           functional_module = assign_module(product))
  write_tsv(scores, paste0('results/next/tables/', group_name, '_nonadditivity_scores_annotated.tsv'))

  if (file.exists(top_path)) {
    top <- read_tsv(top_path, show_col_types = FALSE) %>%
      mutate(gene_id_key = strip_ver(gene_id)) %>%
      left_join(ann, by = 'gene_id_key') %>%
      mutate(product = if_else(is.na(product) | product == '', 'unannotated', product),
             functional_module = assign_module(product))
    write_tsv(top, paste0('results/next/tables/', group_name, '_top50_interaction_candidates_annotated.tsv'))
  }

  class <- read_tsv(class_path, show_col_types = FALSE) %>%
    group_by(target_stage) %>%
    mutate(total = sum(n), pct = 100*n/total) %>%
    ungroup()
  write_tsv(class, paste0('results/next/tables/', group_name, '_nonadditivity_class_summary_pct.tsv'))

  sig <- scores %>% filter(nonadd_class != 'not_significant_or_small')
  module_summary <- sig %>%
    count(target_stage, nonadd_class, functional_module, sort = TRUE) %>%
    group_by(target_stage, nonadd_class) %>%
    mutate(pct = 100*n/sum(n)) %>%
    ungroup()
  write_tsv(module_summary, paste0('results/next/tables/', group_name, '_module_summary.tsv'))

  # simple module enrichment: each target stage significant non-additive vs background
  modules <- sort(unique(scores$functional_module))
  enrich_rows <- list()
  for (st in sort(unique(scores$target_stage))) {
    bg <- scores %>% filter(target_stage == st)
    sig_st <- bg %>% filter(nonadd_class != 'not_significant_or_small')
    N <- nrow(bg); Ksig <- nrow(sig_st)
    for (mod in modules) {
      M <- sum(bg$functional_module == mod)
      k <- sum(sig_st$functional_module == mod)
      p <- if (Ksig > 0 && M > 0) phyper(k-1, M, N-M, Ksig, lower.tail = FALSE) else NA_real_
      enrich_rows[[length(enrich_rows)+1]] <- tibble(
        target_stage = st, functional_module = mod,
        k_sig_module = k, K_sig = Ksig, M_bg_module = M, N_bg = N,
        odds_proxy = (k / max(Ksig,1)) / (M / max(N,1)),
        pvalue = p
      )
    }
  }
  enrich <- bind_rows(enrich_rows) %>%
    mutate(padj = p.adjust(pvalue, method = 'BH')) %>%
    arrange(padj, desc(odds_proxy))
  write_tsv(enrich, paste0('results/next/tables/', group_name, '_module_enrichment.tsv'))

  # top candidates per stage/class/module
  candidates <- scores %>%
    filter(nonadd_class != 'not_significant_or_small') %>%
    arrange(padj, desc(abs(nonadditive_score))) %>%
    group_by(target_stage, nonadd_class) %>%
    slice_head(n = 20) %>% ungroup() %>%
    select(gene_id, target_stage, nonadd_class, functional_module, product, gene_symbol, nonadditive_score, recovery_memory_score, log2FoldChange, padj)
  write_tsv(candidates, paste0('results/next/tables/', group_name, '_candidate_panel_top20_per_class.tsv'))
}

process_group('arabica_Icatu', ann_arabica)
process_group('canephora_CL153', ann_cane)

message('== DONE next-step annotation and module summaries ==')
