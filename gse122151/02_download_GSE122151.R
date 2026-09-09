# ============================================================
# GSE122151 — STEP 2: DOWNLOAD OFFICIAL GEO DATA
# ============================================================

project_root <- getwd()

raw_dir <- file.path(
  project_root,
  "data/raw/GSE122151"
)

dir.create(
  raw_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# Download GEO supplementary files
# ------------------------------------------------------------

GEOquery::getGEOSuppFiles(
  GEO = "GSE122151",
  baseDir = raw_dir,
  makeDirectory = TRUE
)

# ------------------------------------------------------------
# Locate official FPKM table
# ------------------------------------------------------------

fpkm_candidates <- list.files(
  raw_dir,
  pattern = "GSE122151_Gene-fpkm-table\\.txt(\\.gz)?$",
  recursive = TRUE,
  full.names = TRUE
)

if (length(fpkm_candidates) == 0) {

  stop(
    "Could not find GSE122151_Gene-fpkm-table.txt.gz ",
    "after downloading GEO supplementary files."
  )
}

if (length(fpkm_candidates) > 1) {

  stop(
    "More than one GSE122151 FPKM table was found:\n",
    paste(fpkm_candidates, collapse = "\n")
  )
}

fpkm_file <- fpkm_candidates[1]

cat(
  "Official GEO FPKM file:\n",
  fpkm_file,
  "\n"
)

# ------------------------------------------------------------
# Save file manifest
# ------------------------------------------------------------

manifest <- data.frame(
  file = fpkm_file,
  file_size_bytes = file.info(fpkm_file)$size,
  modification_time = file.info(fpkm_file)$mtime,
  stringsAsFactors = FALSE
)

write.csv(
  manifest,
  file.path(
    project_root,
    "results/GSE122151/raw_file_manifest.csv"
  ),
  row.names = FALSE
)

cat("\nDOWNLOAD COMPLETE\n")
