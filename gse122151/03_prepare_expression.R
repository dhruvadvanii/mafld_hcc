# ============================================================
# GSE122151 — STEP 3: PREPARE EXPRESSION MATRIX
# ============================================================
#
# Input:
#   Official GEO GSE122151 Gene-FPKM table
#
# Processing:
#   1. Identify the six HepaRG samples
#   2. Extract gene-level FPKM
#   3. Remove invalid gene identifiers
#   4. Collapse duplicate gene symbols using median FPKM
#   5. Transform to log2(FPKM + 1)
#   6. Save expression matrices
#   7. Save exact sample metadata
#
# ============================================================

project_root <- getwd()

# ------------------------------------------------------------
# Directories
# ------------------------------------------------------------

processed_dir <- file.path(
  project_root,
  "data/processed/GSE122151"
)

results_dir <- file.path(
  project_root,
  "results/GSE122151"
)

dir.create(
  processed_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  results_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# Locate official GEO FPKM file
# ------------------------------------------------------------

fpkm_files <- list.files(
  file.path(
    project_root,
    "data/raw/GSE122151"
  ),
  pattern = "GSE122151_Gene-fpkm-table\\.txt(\\.gz)?$",
  recursive = TRUE,
  full.names = TRUE
)

if (length(fpkm_files) == 0) {

  stop(
    "Could not find the official GSE122151 Gene-FPKM table."
  )
}

if (length(fpkm_files) > 1) {

  stop(
    "More than one GSE122151 Gene-FPKM table was found:\n",
    paste(
      fpkm_files,
      collapse = "\n"
    )
  )
}

fpkm_file <- fpkm_files[1]

cat(
  "\nUsing GSE122151 FPKM file:\n",
  fpkm_file,
  "\n\n"
)

# ------------------------------------------------------------
# Read FPKM table
# ------------------------------------------------------------

fpkm <- data.table::fread(
  fpkm_file,
  data.table = FALSE,
  check.names = FALSE
)

cat(
  "Raw table dimensions:\n",
  nrow(fpkm),
  "rows x",
  ncol(fpkm),
  "columns\n\n"
)

cat(
  "Available columns:\n",
  paste(
    colnames(fpkm),
    collapse = ", "
  ),
  "\n\n"
)

# ------------------------------------------------------------
# Identify gene identifier column
# ------------------------------------------------------------
#
# The actual GSE122151 file contains:
# gene
# gene_name
# gene_type
# chr
# start
# end
# strand
# length
#
# We use gene_name as the primary expression identifier
# because the analysis is performed at the gene-symbol level.
#
# ------------------------------------------------------------

if ("gene_name" %in% colnames(fpkm)) {

  gene_col <- "gene_name"

} else if ("gene" %in% colnames(fpkm)) {

  gene_col <- "gene"

} else {

  stop(
    "Could not identify a gene identifier column."
  )
}

cat(
  "Gene identifier column:",
  gene_col,
  "\n\n"
)

# ------------------------------------------------------------
# Identify EXACT six GSE122151 sample columns
# ------------------------------------------------------------
#
# These are the actual column names in the GEO processed file.
#
# ------------------------------------------------------------

sample_cols <- c(
  "HepaRG1_control_S7",
  "HepaRG2_control_S10",
  "HepaRG3_control_S11",
  "HepaRG5_FA_S12",
  "HepaRG6_FA_S13",
  "HepaRG7_FA_S14"
)

# ------------------------------------------------------------
# Check that all six columns exist
# ------------------------------------------------------------

missing_samples <- setdiff(
  sample_cols,
  colnames(fpkm)
)

if (length(missing_samples) > 0) {

  stop(
    "The following expected GSE122151 sample columns are missing:\n",
    paste(
      missing_samples,
      collapse = "\n"
    ),
    "\n\nAvailable columns are:\n",
    paste(
      colnames(fpkm),
      collapse = ", "
    )
  )
}

# ------------------------------------------------------------
# Check that there are exactly six samples
# ------------------------------------------------------------

if (length(sample_cols) != 6) {

  stop(
    "Internal error: expected exactly six sample columns."
  )
}

cat(
  "Six GSE122151 sample columns successfully identified:\n",
  paste(
    sample_cols,
    collapse = "\n"
  ),
  "\n\n"
)

# ------------------------------------------------------------
# Extract gene identifiers
# ------------------------------------------------------------

genes <- toupper(
  trimws(
    as.character(
      fpkm[[gene_col]]
    )
  )
)

# ------------------------------------------------------------
# Remove missing / empty gene identifiers
# ------------------------------------------------------------

valid_gene <- (
  !is.na(genes) &
  genes != "" &
  genes != "NA"
)

if (sum(!valid_gene) > 0) {

  cat(
    "Removing",
    sum(!valid_gene),
    "rows with missing/empty gene identifiers.\n\n"
  )
}

genes <- genes[valid_gene]

# ------------------------------------------------------------
# Extract expression matrix
# ------------------------------------------------------------

expr_raw <- as.matrix(
  fpkm[
    valid_gene,
    sample_cols,
    drop = FALSE
  ]
)

# ------------------------------------------------------------
# Convert expression matrix to numeric
# ------------------------------------------------------------

suppressWarnings(
  storage.mode(expr_raw) <- "numeric"
)

# ------------------------------------------------------------
# Check for conversion problems
# ------------------------------------------------------------

if (anyNA(expr_raw)) {

  na_count <- sum(
    is.na(expr_raw)
  )

  stop(
    "Expression matrix contains ",
    na_count,
    " NA values after numeric conversion."
  )
}

if (any(!is.finite(expr_raw))) {

  stop(
    "Expression matrix contains non-finite values."
  )
}

if (any(expr_raw < 0)) {

  stop(
    "Negative FPKM values detected. ",
    "This is inconsistent with the expected FPKM input."
  )
}

# ------------------------------------------------------------
# Assign gene identifiers
# ------------------------------------------------------------

rownames(expr_raw) <- genes

colnames(expr_raw) <- sample_cols

# ------------------------------------------------------------
# Count duplicate gene identifiers
# ------------------------------------------------------------

n_duplicate_rows <- sum(
  duplicated(
    rownames(expr_raw)
  )
)

n_unique_genes_before_collapse <- length(
  unique(
    rownames(expr_raw)
  )
)

cat(
  "Rows after removing invalid genes:",
  nrow(expr_raw),
  "\n"
)

cat(
  "Unique gene identifiers:",
  n_unique_genes_before_collapse,
  "\n"
)

cat(
  "Duplicate gene rows:",
  n_duplicate_rows,
  "\n\n"
)

# ------------------------------------------------------------
# Collapse duplicate gene symbols
# ------------------------------------------------------------
#
# If multiple rows correspond to the same gene symbol,
# take the median FPKM across those rows.
#
# ------------------------------------------------------------

gene_index <- split(
  seq_len(nrow(expr_raw)),
  rownames(expr_raw)
)

expr_collapsed <- do.call(
  rbind,
  lapply(
    gene_index,
    function(idx) {

      if (length(idx) == 1) {

        as.numeric(
          expr_raw[idx, ]
        )

      } else {

        apply(
          expr_raw[idx, , drop = FALSE],
          2,
          median,
          na.rm = TRUE
        )
      }
    }
  )
)

expr_collapsed <- as.matrix(
  expr_collapsed
)

rownames(expr_collapsed) <- names(
  gene_index
)

colnames(expr_collapsed) <- sample_cols

storage.mode(expr_collapsed) <- "numeric"

# ------------------------------------------------------------
# Validate collapsed matrix
# ------------------------------------------------------------

if (ncol(expr_collapsed) != 6) {

  stop(
    "After duplicate-gene collapsing, the expression matrix ",
    "does not contain exactly six samples."
  )
}

if (any(!is.finite(expr_collapsed))) {

  stop(
    "Non-finite values detected after duplicate-gene collapsing."
  )
}

if (any(expr_collapsed < 0)) {

  stop(
    "Negative values detected after duplicate-gene collapsing."
  )
}

# ------------------------------------------------------------
# log2(FPKM + 1)
# ------------------------------------------------------------

expr_log2 <- log2(
  expr_collapsed + 1
)

# ------------------------------------------------------------
# Final expression matrix validation
# ------------------------------------------------------------

if (any(!is.finite(expr_log2))) {

  stop(
    "Non-finite values detected in log2(FPKM + 1) matrix."
  )
}

if (nrow(expr_log2) < 1000) {

  stop(
    "Unexpectedly low number of genes after preprocessing: ",
    nrow(expr_log2)
  )
}

# ------------------------------------------------------------
# Save FPKM matrix as RDS
# ------------------------------------------------------------

saveRDS(
  expr_collapsed,
  file = file.path(
    processed_dir,
    "GSE122151_FPKM_gene_symbol.rds"
  )
)

# ------------------------------------------------------------
# Save log2(FPKM + 1) matrix as RDS
# ------------------------------------------------------------

saveRDS(
  expr_log2,
  file = file.path(
    processed_dir,
    "GSE122151_log2_FPKM_plus1.rds"
  )
)

# ------------------------------------------------------------
# Save FPKM matrix as CSV
# ------------------------------------------------------------

fpkm_export <- data.frame(
  gene = rownames(expr_collapsed),
  expr_collapsed,
  check.names = FALSE
)

write.csv(
  fpkm_export,
  file = file.path(
    processed_dir,
    "GSE122151_FPKM_gene_symbol.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)

# ------------------------------------------------------------
# Save log2 matrix as CSV
# ------------------------------------------------------------

log2_export <- data.frame(
  gene = rownames(expr_log2),
  expr_log2,
  check.names = FALSE
)

write.csv(
  log2_export,
  file = file.path(
    processed_dir,
    "GSE122151_log2_FPKM_plus1.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)

# ------------------------------------------------------------
# Exact GSM/sample mapping
# ------------------------------------------------------------
#
# Based on the GSE122151 sample annotations:
#
# GSM3456189 = HepaRG1_control_S7
# GSM3456190 = HepaRG2_control_S10
# GSM3456191 = HepaRG3_control_S11
# GSM3456192 = HepaRG5_FA_S12
# GSM3456193 = HepaRG6_FA_S13
# GSM3456194 = HepaRG7_FA_S14
#
# ------------------------------------------------------------

sample_metadata <- data.frame(

  GSM = c(
    "GSM3456189",
    "GSM3456190",
    "GSM3456191",
    "GSM3456192",
    "GSM3456193",
    "GSM3456194"
  ),

  sample_name = c(
    "HepaRG1_control_S7",
    "HepaRG2_control_S10",
    "HepaRG3_control_S11",
    "HepaRG5_FA_S12",
    "HepaRG6_FA_S13",
    "HepaRG7_FA_S14"
  ),

  condition = c(
    "Control",
    "Control",
    "Control",
    "Lipid_Loaded",
    "Lipid_Loaded",
    "Lipid_Loaded"
  ),

  replicate = c(
    1,
    2,
    3,
    1,
    2,
    3
  ),

  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# Make sample names the row names
# ------------------------------------------------------------

rownames(sample_metadata) <-
  sample_metadata$sample_name

# ------------------------------------------------------------
# Validate metadata against expression matrix
# ------------------------------------------------------------

if (
  !identical(
    sample_metadata$sample_name,
    colnames(expr_log2)
  )
) {

  stop(
    "Sample metadata order does not exactly match ",
    "expression matrix column order."
  )
}

# ------------------------------------------------------------
# Save metadata
# ------------------------------------------------------------

write.csv(
  sample_metadata,
  file = file.path(
    processed_dir,
    "sample_metadata.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)

saveRDS(
  sample_metadata,
  file = file.path(
    processed_dir,
    "sample_metadata.rds"
  )
)

# ------------------------------------------------------------
# Create preprocessing summary
# ------------------------------------------------------------

qc_summary <- data.frame(

  metric = c(
    "Raw rows",
    "Raw columns",
    "Rows after invalid gene removal",
    "Unique genes after collapsing",
    "Duplicate gene rows collapsed",
    "Number of samples",
    "Control samples",
    "Lipid-loaded samples"
  ),

  value = c(
    nrow(fpkm),
    ncol(fpkm),
    nrow(expr_raw),
    nrow(expr_collapsed),
    n_duplicate_rows,
    ncol(expr_log2),
    sum(
      sample_metadata$condition == "Control"
    ),
    sum(
      sample_metadata$condition == "Lipid_Loaded"
    )
  ),

  stringsAsFactors = FALSE
)

write.csv(
  qc_summary,
  file = file.path(
    results_dir,
    "expression_preprocessing_summary.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)

# ------------------------------------------------------------
# Final console report
# ------------------------------------------------------------

cat(
  "\n============================================================\n"
)

cat(
  "GSE122151 EXPRESSION PREPROCESSING COMPLETE\n"
)

cat(
  "============================================================\n\n"
)

cat(
  "Final genes:",
  nrow(expr_log2),
  "\n"
)

cat(
  "Final samples:",
  ncol(expr_log2),
  "\n"
)

cat(
  "Control:",
  sum(
    sample_metadata$condition == "Control"
  ),
  "\n"
)

cat(
  "Lipid-loaded:",
  sum(
    sample_metadata$condition == "Lipid_Loaded"
  ),
  "\n\n"
)

cat(
  "Expression transformation:\n"
)

cat(
  "log2(FPKM + 1)\n\n"
)

cat(
  "Saved files:\n"
)

cat(
  "- GSE122151_FPKM_gene_symbol.rds\n"
)

cat(
  "- GSE122151_FPKM_gene_symbol.csv\n"
)

cat(
  "- GSE122151_log2_FPKM_plus1.rds\n"
)

cat(
  "- GSE122151_log2_FPKM_plus1.csv\n"
)

cat(
  "- sample_metadata.csv\n"
)

cat(
  "- sample_metadata.rds\n"
)

cat(
  "- expression_preprocessing_summary.csv\n\n"
)

cat(
  "============================================================\n"
)
