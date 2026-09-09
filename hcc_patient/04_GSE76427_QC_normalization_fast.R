# ============================================================
# MAFLD-HCC MANUSCRIPT PIPELINE
# STEP 4 — GSE76427 FAST QC + NORMALIZATION
# ============================================================
#
# PURPOSE:
#   Convert validated GSE76427 Illumina signal data into an
#   analysis-ready gene-level expression matrix.
#
# PIPELINE:
#
#   AVG_Signal
#       ↓
#   log2 transformation
#       ↓
#   Detection P-value filtering
#       ↓
#   Quantile normalization
#       ↓
#   Probe → gene mapping
#       ↓
#   Median collapse of multiple probes per gene
#       ↓
#   Final gene-level expression matrix
#
# This version intentionally omits:
#   - PCA
#   - density plots
#   - boxplots
#   - expensive visualization steps
#
# These can be performed separately later if required.
#
# ============================================================


# ------------------------------------------------------------
# 1. PROJECT DIRECTORY
# ------------------------------------------------------------

project_dir <- "C:/Users/Dhruv Advani/OneDrive/Desktop/AIIMS/Research/MAFLD_HCC_Pipeline"

setwd(project_dir)

cat("\nProject directory:\n")
cat(normalizePath(getwd()), "\n\n")


# ------------------------------------------------------------
# 2. DIRECTORIES
# ------------------------------------------------------------

processed_dir <- file.path(
  "data",
  "processed",
  "GSE76427"
)

results_dir <- "results"

if (!dir.exists(processed_dir)) {
  dir.create(
    processed_dir,
    recursive = TRUE
  )
}

if (!dir.exists(results_dir)) {
  dir.create(
    results_dir,
    recursive = TRUE
  )
}


# ------------------------------------------------------------
# 3. REQUIRED PACKAGES
# ------------------------------------------------------------

required_packages <- c(
  "limma"
)

missing_packages <- required_packages[
  !sapply(
    required_packages,
    requireNamespace,
    quietly = TRUE
  )
]

if (length(missing_packages) > 0) {

  stop(
    "Missing required package(s): ",
    paste(
      missing_packages,
      collapse = ", "
    ),
    "\nInstall before continuing."
  )
}

library(limma)


# ------------------------------------------------------------
# 4. LOAD STEP 3 OBJECT
# ------------------------------------------------------------

input_rds <- file.path(
  processed_dir,
  "GSE76427_raw_expression.rds"
)

if (!file.exists(input_rds)) {

  stop(
    "Step 3 output not found:\n",
    normalizePath(
      input_rds,
      mustWork = FALSE
    )
  )
}

cat("Loading Step 3 expression object...\n")

gse_raw <- readRDS(
  input_rds
)


# ------------------------------------------------------------
# 5. EXTRACT MATRICES
# ------------------------------------------------------------

signal_matrix <- gse_raw$signal_matrix

detection_matrix <- gse_raw$detection_matrix

probe_mapping <- gse_raw$probe_mapping


cat("\n============================================================\n")
cat("INPUT\n")
cat("============================================================\n\n")

cat(
  "Samples: ",
  ncol(signal_matrix),
  "\n",
  sep = ""
)

cat(
  "Probes: ",
  nrow(signal_matrix),
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------
# 6. VERIFY MATRICES
# ------------------------------------------------------------

if (!identical(
  dim(signal_matrix),
  dim(detection_matrix)
)) {

  stop(
    "Signal and detection matrices have different dimensions."
  )
}

if (!identical(
  rownames(signal_matrix),
  rownames(detection_matrix)
)) {

  stop(
    "Signal and detection matrices have different probe order."
  )
}

if (!identical(
  colnames(signal_matrix),
  colnames(detection_matrix)
)) {

  stop(
    "Signal and detection matrices have different sample order."
  )
}


# ------------------------------------------------------------
# 7. LOG2 TRANSFORMATION
# ------------------------------------------------------------
#
# The original AVG_Signal values are positive intensity
# measurements and strongly right-skewed.
#
# Therefore:
#
#   log2(AVG_Signal)
#
# No pseudocount is required because Step 3 showed that
# there are no zero values.
# ------------------------------------------------------------

cat("Applying log2 transformation...\n")

if (any(
  signal_matrix <= 0,
  na.rm = TRUE
)) {

  stop(
    "Non-positive signal values detected. ",
    "Do not proceed with this version."
  )
}

log2_signal <- log2(
  signal_matrix
)

cat("Log2 transformation completed.\n\n")


# ------------------------------------------------------------
# 8. DETECTION P-VALUE FILTER
# ------------------------------------------------------------
#
# Retain probes detected at:
#
#   P < 0.05
#
# in at least 20% of samples.
#
# For 167 samples:
#
#   ceiling(167 × 0.20) = 34 samples
# ------------------------------------------------------------

detection_threshold <- 0.05

minimum_detected_samples <- ceiling(
  ncol(detection_matrix) * 0.20
)


detected_samples_per_probe <- rowSums(
  detection_matrix <
    detection_threshold,
  na.rm = TRUE
)


probe_keep <- (
  detected_samples_per_probe >=
    minimum_detected_samples
)


cat("============================================================\n")
cat("DETECTION FILTER\n")
cat("============================================================\n\n")

cat(
  "Detection P threshold: < ",
  detection_threshold,
  "\n",
  sep = ""
)

cat(
  "Minimum detected samples: ",
  minimum_detected_samples,
  "\n",
  sep = ""
)

cat(
  "Probes before filtering: ",
  nrow(signal_matrix),
  "\n",
  sep = ""
)

cat(
  "Probes retained: ",
  sum(probe_keep),
  "\n",
  sep = ""
)

cat(
  "Probes removed: ",
  sum(!probe_keep),
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------
# 9. FILTER PROBE MATRIX
# ------------------------------------------------------------

log2_filtered <- log2_signal[
  probe_keep,
  ,
  drop = FALSE
]

detection_filtered <- detection_matrix[
  probe_keep,
  ,
  drop = FALSE
]


probe_mapping_filtered <- probe_mapping[
  match(
    rownames(log2_filtered),
    probe_mapping$probe_id
  ),
  ,
  drop = FALSE
]


# ------------------------------------------------------------
# 10. QUANTILE NORMALIZATION
# ------------------------------------------------------------
#
# limma::normalizeBetweenArrays()
#
# is applied to the log2-transformed, detection-filtered
# single-channel microarray data.
# ------------------------------------------------------------

cat("============================================================\n")
cat("QUANTILE NORMALIZATION\n")
cat("============================================================\n\n")

cat(
  "Running limma::normalizeBetweenArrays()...\n"
)

normalized_probe_matrix <- normalizeBetweenArrays(
  log2_filtered,
  method = "quantile"
)

cat(
  "Quantile normalization completed.\n\n"
)


# ------------------------------------------------------------
# 11. SAVE NORMALIZED PROBE MATRIX
# ------------------------------------------------------------

saveRDS(
  normalized_probe_matrix,
  file = file.path(
    processed_dir,
    "GSE76427_normalized_probe_matrix.rds"
  )
)

cat(
  "Normalized probe matrix saved.\n\n"
)


# ------------------------------------------------------------
# 12. MAP PROBES TO GENE SYMBOLS
# ------------------------------------------------------------

cat("============================================================\n")
cat("PROBE → GENE MAPPING\n")
cat("============================================================\n\n")


probe_symbols <- trimws(
  probe_mapping_filtered$symbol
)


probe_symbols[
  probe_symbols == ""
] <- NA


valid_mapping <- !is.na(
  probe_symbols
)


cat(
  "Retained probes: ",
  nrow(normalized_probe_matrix),
  "\n",
  sep = ""
)

cat(
  "Probes with gene symbols: ",
  sum(valid_mapping),
  "\n",
  sep = ""
)

cat(
  "Probes without gene symbols: ",
  sum(!valid_mapping),
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------
# 13. KEEP ONLY MAPPED PROBES
# ------------------------------------------------------------

mapped_expression <- normalized_probe_matrix[
  valid_mapping,
  ,
  drop = FALSE
]

mapped_symbols <- probe_symbols[
  valid_mapping
]


# ------------------------------------------------------------
# 14. COLLAPSE MULTIPLE PROBES PER GENE
# ------------------------------------------------------------
#
# Multiple probes representing the same gene are collapsed
# using the median normalized expression.
#
# This is done efficiently using split() rather than a
# gene-by-gene loop.
# ------------------------------------------------------------

cat(
  "Collapsing multiple probes per gene using median...\n"
)


gene_groups <- split(
  seq_along(mapped_symbols),
  mapped_symbols
)


gene_expression <- vapply(
  gene_groups,
  function(idx) {

    apply(
      mapped_expression[
        idx,
        ,
        drop = FALSE
      ],
      2,
      median,
      na.rm = TRUE
    )

  },
  numeric(
    ncol(mapped_expression)
  )
)


gene_expression <- t(
  gene_expression
)


# ------------------------------------------------------------
# 15. FIX DIMENSION NAMES
# ------------------------------------------------------------

rownames(gene_expression) <- names(
  gene_groups
)

colnames(gene_expression) <- colnames(
  mapped_expression
)


# ------------------------------------------------------------
# 16. FINAL MATRIX VALIDATION
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("FINAL GENE-LEVEL MATRIX\n")
cat("============================================================\n\n")

cat(
  "Genes: ",
  nrow(gene_expression),
  "\n",
  sep = ""
)

cat(
  "Samples: ",
  ncol(gene_expression),
  "\n\n",
  sep = ""
)


if (anyDuplicated(
  rownames(gene_expression)
) > 0) {

  stop(
    "Duplicate gene symbols remain."
  )
}


if (any(
  !is.finite(gene_expression)
)) {

  stop(
    "Non-finite values detected in gene expression matrix."
  )
}


if (!identical(
  colnames(gene_expression),
  colnames(signal_matrix)
)) {

  stop(
    "Sample names changed during processing."
  )
}


# ------------------------------------------------------------
# 17. GENE-LEVEL SUMMARY
# ------------------------------------------------------------

probe_count_per_gene <- table(
  mapped_symbols
)


gene_mapping_summary <- data.frame(

  metric = c(
    "initial_probes",
    "probes_after_detection_filter",
    "mapped_probes",
    "unmapped_probes",
    "unique_genes",
    "genes_with_multiple_probes"
  ),

  value = c(

    nrow(signal_matrix),

    nrow(normalized_probe_matrix),

    sum(valid_mapping),

    sum(!valid_mapping),

    nrow(gene_expression),

    sum(
      probe_count_per_gene > 1
    )
  ),

  stringsAsFactors = FALSE
)


write.csv(
  gene_mapping_summary,
  file = file.path(
    results_dir,
    "04_GSE76427_gene_mapping_summary.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 18. SAMPLE-LEVEL QC SUMMARY
# ------------------------------------------------------------
#
# We record summary statistics without generating expensive
# plots.
# ------------------------------------------------------------

sample_qc <- data.frame(

  sample =
    colnames(signal_matrix),

  raw_median =
    apply(
      signal_matrix,
      2,
      median,
      na.rm = TRUE
    ),

  log2_median =
    apply(
      log2_signal,
      2,
      median,
      na.rm = TRUE
    ),

  normalized_median =
    apply(
      normalized_probe_matrix,
      2,
      median,
      na.rm = TRUE
    ),

  percent_detected =
    100 *
    colMeans(
      detection_matrix <
        detection_threshold,
      na.rm = TRUE
    ),

  stringsAsFactors = FALSE
)


write.csv(
  sample_qc,
  file = file.path(
    results_dir,
    "04_GSE76427_sample_QC.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 19. FILTERING SUMMARY
# ------------------------------------------------------------

filtering_summary <- data.frame(

  metric = c(
    "samples",
    "initial_probes",
    "detection_threshold",
    "minimum_detection_fraction",
    "minimum_detected_samples",
    "probes_retained",
    "probes_removed",
    "percent_retained"
  ),

  value = c(

    ncol(signal_matrix),

    nrow(signal_matrix),

    detection_threshold,

    0.20,

    minimum_detected_samples,

    sum(probe_keep),

    sum(!probe_keep),

    100 *
      sum(probe_keep) /
      nrow(signal_matrix)
  ),

  stringsAsFactors = FALSE
)


write.csv(
  filtering_summary,
  file = file.path(
    results_dir,
    "04_GSE76427_filtering_summary.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 20. SAVE FINAL GENE-LEVEL OBJECT
# ------------------------------------------------------------

gse76427_gene_expression <- list(

  created_at =
    Sys.time(),

  source =
    "GSE76427_non-normalized.txt",

  platform =
    "GPL10558 — Illumina HumanHT-12 v4",

  n_samples =
    ncol(gene_expression),

  n_genes =
    nrow(gene_expression),

  transformation =
    "log2",

  normalization =
    "quantile",

  detection_p_threshold =
    detection_threshold,

  minimum_detection_fraction =
    0.20,

  minimum_detected_samples =
    minimum_detected_samples,

  probe_collapsing =
    "median",

  sample_names =
    colnames(gene_expression),

  gene_expression =
    gene_expression
)


final_rds <- file.path(
  processed_dir,
  "GSE76427_gene_expression.rds"
)


saveRDS(
  gse76427_gene_expression,
  file = final_rds
)


# ------------------------------------------------------------
# 21. WRITE PROCESSING REPORT
# ------------------------------------------------------------

report_lines <- c(

  "============================================================",
  "GSE76427 FAST PROCESSING REPORT",
  "============================================================",
  "",

  paste(
    "Date:",
    Sys.time()
  ),

  "",

  "INPUT",
  "------------------------------------------------------------",

  paste(
    "Samples:",
    ncol(signal_matrix)
  ),

  paste(
    "Initial probes:",
    nrow(signal_matrix)
  ),

  "",

  "TRANSFORMATION",
  "------------------------------------------------------------",

  "log2 transformation: YES",

  "",

  "DETECTION FILTER",
  "------------------------------------------------------------",

  paste(
    "Detection P threshold:",
    detection_threshold
  ),

  paste(
    "Minimum detection fraction:",
    "20%"
  ),

  paste(
    "Minimum detected samples:",
    minimum_detected_samples
  ),

  paste(
    "Probes retained:",
    sum(probe_keep)
  ),

  paste(
    "Probes removed:",
    sum(!probe_keep)
  ),

  "",

  "NORMALIZATION",
  "------------------------------------------------------------",

  "Method: limma quantile normalization",

  "",

  "GENE MAPPING",
  "------------------------------------------------------------",

  paste(
    "Mapped probes:",
    sum(valid_mapping)
  ),

  paste(
    "Unique genes:",
    nrow(gene_expression)
  ),

  paste(
    "Genes with multiple probes:",
    sum(
      probe_count_per_gene > 1
    )
  ),

  "Multiple probes collapsed using median.",

  "",

  "FINAL MATRIX",
  "------------------------------------------------------------",

  paste(
    "Genes:",
    nrow(gene_expression)
  ),

  paste(
    "Samples:",
    ncol(gene_expression)
  ),

  "",

  "QC",
  "------------------------------------------------------------",

  "Sample-level QC summary: YES",

  "PCA: NOT PERFORMED",

  "Density plots: NOT PERFORMED",

  "Boxplots: NOT PERFORMED",

  "Automatic sample removal: NO",

  "",

  "FINAL OUTPUT",
  "------------------------------------------------------------",

  final_rds
)


writeLines(
  report_lines,
  con = file.path(
    results_dir,
    "04_GSE76427_QC_report.txt"
  )
)


# ------------------------------------------------------------
# 22. SAVE PROCESSING PARAMETERS
# ------------------------------------------------------------

processing_parameters <- list(

  created_at =
    Sys.time(),

  input =
    normalizePath(
      input_rds
    ),

  samples =
    ncol(signal_matrix),

  initial_probes =
    nrow(signal_matrix),

  transformation =
    "log2",

  detection_p_threshold =
    detection_threshold,

  minimum_detection_fraction =
    0.20,

  minimum_detected_samples =
    minimum_detected_samples,

  probes_retained =
    sum(probe_keep),

  normalization =
    "quantile",

  probe_collapsing =
    "median",

  final_genes =
    nrow(gene_expression),

  final_samples =
    ncol(gene_expression)
)


saveRDS(
  processing_parameters,
  file = file.path(
    results_dir,
    "04_GSE76427_processing_parameters.rds"
  )
)


# ------------------------------------------------------------
# 23. CLEAN MEMORY
# ------------------------------------------------------------

rm(
  gse_raw,
  signal_matrix,
  detection_matrix,
  log2_signal,
  log2_filtered,
  detection_filtered,
  normalized_probe_matrix,
  mapped_expression
)

gc()


# ------------------------------------------------------------
# 24. FINAL STATUS
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("STEP 4 COMPLETED SUCCESSFULLY\n")
cat("============================================================\n\n")

cat(
  "Samples: ",
  ncol(gene_expression),
  "\n",
  sep = ""
)

cat(
  "Initial probes: ",
  nrow(signal_matrix),
  "\n",
  sep = ""
)

cat(
  "Probes retained after detection filtering: ",
  sum(probe_keep),
  "\n",
  sep = ""
)

cat(
  "Final genes: ",
  nrow(gene_expression),
  "\n\n",
  sep = ""
)

cat(
  "Transformation: log2\n"
)

cat(
  "Normalization: quantile\n"
)

cat(
  "Detection filter: P < 0.05 in >=20% samples\n"
)

cat(
  "Probe collapsing: median\n\n"
)

cat(
  "FINAL OUTPUT:\n",
  final_rds,
  "\n\n",
  sep = ""
)

cat(
  "QC summary:\n",
  file.path(
    results_dir,
    "04_GSE76427_sample_QC.csv"
  ),
  "\n"
)

cat(
  "\nNo samples were automatically removed.\n"
)

cat(
  "\nNext step: review the numerical QC summary and proceed\n",
  "to phenotype/sample matching before pathway scoring.\n"
)
