suppressPackageStartupMessages({
  library(openxlsx)
})

message('== Building analysis table workbook ==')
out_dir <- 'results/paper_pack_v1/tables'
qc_dir <- 'results/paper_pack_v1/qc'
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

files <- c(
  list.files(qc_dir, pattern = '\\.(tsv|csv)$', full.names = TRUE),
  list.files(out_dir, pattern = '\\.(tsv|csv)$', full.names = TRUE)
)
files <- unique(files[file.exists(files)])

if (length(files) == 0) {
  stop('No table files found for workbook generation.')
}

wb <- openxlsx::createWorkbook()
used <- character()
for (f in files) {
  df <- if (grepl('\\.csv$', f, ignore.case = TRUE)) {
    read.csv(f, check.names = FALSE)
  } else {
    read.delim(f, check.names = FALSE)
  }
  nm <- tools::file_path_sans_ext(basename(f))
  nm <- gsub('[^A-Za-z0-9_]', '_', nm)
  base <- substr(nm, 1, 28)
  sheet <- base
  i <- 1
  while (sheet %in% used) {
    suffix <- paste0('_', i)
    sheet <- paste0(substr(base, 1, 31 - nchar(suffix)), suffix)
    i <- i + 1
  }
  used <- c(used, sheet)
  openxlsx::addWorksheet(wb, sheet)
  openxlsx::writeData(wb, sheet, df)
}

out <- 'results/paper_pack_v1/coffee_analysis_tables.xlsx'
openxlsx::saveWorkbook(wb, out, overwrite = TRUE)
message('Workbook written: ', out)
message('== DONE analysis table workbook ==')
