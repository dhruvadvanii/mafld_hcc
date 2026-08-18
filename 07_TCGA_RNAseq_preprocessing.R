# ============================================================
# MAFLD-HCC MANUSCRIPT PIPELINE
# STEP 7 — TCGA-LIHC RNA-seq PREPROCESSING
# ============================================================
#
# PURPOSE:
#   Convert TCGA-LIHC STAR-aligned RNA-seq counts into an
#   analysis-ready expression matrix for downstream GSVA.
#
# DISCOVERY COHORT:
#   ALL 371 TCGA-LIHC PRIMARY TUMOR SAMPLES
#
# PIPELINE:
#
#   STAR counts
#       ↓
#   371 primary tumors
#       ↓
#   Remove genes with zero counts in all samples
#       ↓
#   CPM filtering:
#       CPM >= 1 in >=20% of samples
#       ↓
#   edgeR TMM normalization
#       ↓
#   limma-voom transformation
#       ↓
#   Analysis-ready expression matrix
#
# IMPORTANT:
#   - No n=163 filtering
#   - No MAF required
#   - No GSVA
#   - No survival-based exclusions
#   - No automatic sample removal
#
# RATIONALE FOR FILTER:
#
#   This is a pathway/GSVA study rather than a conventional
#   differential-expression analysis. Therefore, we avoid the
#   overly stringent default filterByExpr() behavior when no
#   biological group/design has been specified.
#
#   A gene is retained if:
#
#       CPM >= 1
#
#   in at least:
#
#       20% of the 371 TCGA tumors
#
#   Minimum number of samples:
#
#       ceiling(371 * 0.20) = 75
#
# ============================================================


# ------------------------------------------------------------
# 1. SET PROJECT DIRECTORY
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
  "TCGA-LIHC"
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
  "edgeR",
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
    "The following required packages are missing:\n",
    paste(
      missing_packages,
      collapse = ", "
    ),
    "\n\nInstall them before running Step 7."
  )
}


library(edgeR)
library(limma)


# ------------------------------------------------------------
# 4. LOAD STEP 6 TCGA OBJECT
# ------------------------------------------------------------

input_rds <- file.path(
  processed_dir,
  "TCGA_LIHC_raw_expression.rds"
)


if (!file.exists(input_rds)) {

  stop(
    "Step 6 TCGA expression object not found:\n",
    normalizePath(
      input_rds,
      mustWork = FALSE
    ),
    "\n\nRun Step 6 first."
  )
}


cat(
  "Loading TCGA primary-tumor count matrix...\n"
)


tcga_object <- readRDS(
  input_rds
)


counts <- tcga_object$expression_counts


# ------------------------------------------------------------
# 5. INITIAL VALIDATION
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("INITIAL TCGA COUNT MATRIX\n")
cat("============================================================\n\n")


cat(
  "Genes: ",
  nrow(counts),
  "\n",
  sep = ""
)


cat(
  "Samples: ",
  ncol(counts),
  "\n\n",
  sep = ""
)


# Expected discovery cohort

if (
  ncol(counts) != 371
) {

  stop(
    "Expected 371 primary tumor samples, found ",
    ncol(counts),
    "."
  )
}


if (
  is.null(rownames(counts))
) {

  stop(
    "Gene identifiers are missing."
  )
}


if (
  is.null(colnames(counts))
) {

  stop(
    "Sample identifiers are missing."
  )
}


# ------------------------------------------------------------
# 6. COUNT MATRIX VALIDATION
# ------------------------------------------------------------

cat(
  "Checking count matrix values...\n"
)


if (anyNA(counts)) {

  stop(
    "NA values detected in count matrix."
  )
}


if (any(!is.finite(counts))) {

  stop(
    "Non-finite values detected in count matrix."
  )
}


if (any(counts < 0)) {

  stop(
    "Negative count values detected."
  )
}


integer_check <- all(
  counts == floor(counts)
)


cat(
  "All values integer-valued: ",
  integer_check,
  "\n",
  sep = ""
)


if (!integer_check) {

  warning(
    "Some count values are not integer-valued. ",
    "Inspect the input before proceeding."
  )
}


# ------------------------------------------------------------
# 7. REMOVE GENES WITH ZERO COUNTS IN ALL SAMPLES
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("ZERO-COUNT FILTER\n")
cat("============================================================\n\n")


zero_count_genes <- rowSums(
  counts
) == 0


cat(
  "Genes before zero-count filter: ",
  nrow(counts),
  "\n",
  sep = ""
)


cat(
  "Genes with zero counts in all samples: ",
  sum(zero_count_genes),
  "\n",
  sep = ""
)


counts_nonzero <- counts[
  !zero_count_genes,
  ,
  drop = FALSE
]


cat(
  "Genes remaining: ",
  nrow(counts_nonzero),
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------
# 8. CREATE edgeR DGEList
# ------------------------------------------------------------

cat(
  "Creating edgeR DGEList...\n"
)


dge <- DGEList(
  counts = counts_nonzero
)


# ------------------------------------------------------------
# 9. CALCULATE CPM
# ------------------------------------------------------------
#
# CPM = counts per million mapped reads.
#
# We use CPM rather than raw counts because library sizes
# differ between TCGA samples.
# ------------------------------------------------------------

cat(
  "Calculating counts per million (CPM)...\n"
)


cpm_matrix <- cpm(
  dge,
  normalized.lib.sizes = FALSE
)


# ------------------------------------------------------------
# 10. LOW-EXPRESSION FILTER
# ------------------------------------------------------------
#
# Primary criterion:
#
#   CPM >= 1
#
# in at least 20% of the 371 samples.
#
# 20% of 371 = 74.2
#
# Therefore at least 75 samples.
# ------------------------------------------------------------

cpm_threshold <- 1

minimum_fraction <- 0.20

minimum_samples <- ceiling(
  ncol(dge) *
    minimum_fraction
)


cat("\n")
cat("============================================================\n")
cat("LOW-EXPRESSION FILTER\n")
cat("============================================================\n\n")


cat(
  "CPM threshold: >= ",
  cpm_threshold,
  "\n",
  sep = ""
)


cat(
  "Minimum sample fraction: ",
  minimum_fraction * 100,
  "%\n",
  sep = ""
)


cat(
  "Minimum number of samples: ",
  minimum_samples,
  "\n\n",
  sep = ""
)


# Number of samples in which each gene reaches CPM >= 1

samples_above_cpm <- rowSums(
  cpm_matrix >= cpm_threshold
)


keep_genes <- (
  samples_above_cpm >=
    minimum_samples
)


cat(
  "Genes before expression filter: ",
  nrow(dge),
  "\n",
  sep = ""
)


cat(
  "Genes retained: ",
  sum(keep_genes),
  "\n",
  sep = ""
)


cat(
  "Genes removed: ",
  sum(!keep_genes),
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------
# 11. APPLY EXPRESSION FILTER
# ------------------------------------------------------------

dge_filtered <- dge[
  keep_genes,
  ,
  keep.lib.sizes = FALSE
]


# ------------------------------------------------------------
# 12. RE-CALCULATE CPM AFTER FILTERING
# ------------------------------------------------------------

cpm_filtered <- cpm(
  dge_filtered,
  normalized.lib.sizes = FALSE
)


# ------------------------------------------------------------
# 13. TMM NORMALIZATION
# ------------------------------------------------------------
#
# TMM = Trimmed Mean of M-values.
#
# This adjusts for differences in RNA-seq library size and
# composition between samples.
# ------------------------------------------------------------

cat("============================================================\n")
cat("TMM NORMALIZATION\n")
cat("============================================================\n\n")


cat(
  "Running edgeR::calcNormFactors()...\n"
)


dge_filtered <- calcNormFactors(
  dge_filtered,
  method = "TMM"
)


cat(
  "TMM normalization completed.\n\n"
)


# ------------------------------------------------------------
# 14. LIBRARY SIZE SUMMARY
# ------------------------------------------------------------

library_summary <- data.frame(

  sample =
    colnames(dge_filtered),

  raw_library_size =
    dge_filtered$samples$lib.size,

  TMM_normalization_factor =
    dge_filtered$samples$norm.factors,

  effective_library_size =
    dge_filtered$samples$lib.size *
    dge_filtered$samples$norm.factors,

  stringsAsFactors = FALSE
)


library_summary$relative_effective_library_size <-
  library_summary$effective_library_size /
  median(
    library_summary$effective_library_size
  )


write.csv(
  library_summary,
  file = file.path(
    results_dir,
    "07_TCGA_library_size_summary.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 15. LIMMA-VOOM TRANSFORMATION
# ------------------------------------------------------------
#
# We use an intercept-only design because this step is not
# performing differential expression.
#
# The purpose of voom here is to generate a normalized,
# log-expression matrix suitable for downstream continuous
# expression/pathway analyses.
# ------------------------------------------------------------

cat("============================================================\n")
cat("LIMMA-VOOM TRANSFORMATION\n")
cat("============================================================\n\n")


cat(
  "Running limma::voom()...\n"
)


voom_design <- matrix(
  1,
  nrow = ncol(dge_filtered),
  ncol = 1
)


colnames(
  voom_design
) <- "Intercept"


rownames(
  voom_design
) <- colnames(
  dge_filtered
)


voom_object <- voom(
  dge_filtered,
  design = voom_design,
  plot = FALSE
)


cat(
  "Voom transformation completed.\n\n"
)


voom_expression <- voom_object$E


# ------------------------------------------------------------
# 16. VERIFY VOOM MATRIX
# ------------------------------------------------------------

cat(
  "============================================================\n"
)

cat(
  "FINAL VOOM MATRIX\n"
)

cat(
  "============================================================\n\n"
)


cat(
  "Genes: ",
  nrow(voom_expression),
  "\n",
  sep = ""
)


cat(
  "Samples: ",
  ncol(voom_expression),
  "\n\n",
  sep = ""
)


if (anyNA(voom_expression)) {

  stop(
    "NA values detected in voom matrix."
  )
}


if (any(!is.finite(voom_expression))) {

  stop(
    "Non-finite values detected in voom matrix."
  )
}


if (!identical(
  colnames(voom_expression),
  colnames(counts)
)) {

  stop(
    "Sample order changed during preprocessing."
  )
}


if (!identical(
  rownames(voom_expression),
  rownames(dge_filtered$counts)
)) {

  stop(
    "Gene order changed unexpectedly."
  )
}


# ------------------------------------------------------------
# 17. EXPRESSION SUMMARY
# ------------------------------------------------------------

expression_summary <- data.frame(

  metric = c(

    "initial_genes",

    "zero_count_genes_removed",

    "genes_after_zero_filter",

    "CPM_threshold",

    "minimum_fraction",

    "minimum_samples",

    "genes_after_CPM_filter",

    "genes_removed_by_CPM_filter",

    "samples",

    "NA_values_in_voom",

    "non_finite_values_in_voom"
  ),

  value = c(

    nrow(counts),

    sum(zero_count_genes),

    nrow(counts_nonzero),

    cpm_threshold,

    minimum_fraction,

    minimum_samples,

    nrow(dge_filtered),

    sum(!keep_genes),

    ncol(voom_expression),

    sum(
      is.na(
        voom_expression
      )
    ),

    sum(
      !is.finite(
        voom_expression
      )
    )
  ),

  stringsAsFactors = FALSE
)


write.csv(
  expression_summary,
  file = file.path(
    results_dir,
    "07_TCGA_expression_summary.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 18. FILTERING SUMMARY
# ------------------------------------------------------------

filtering_summary <- data.frame(

  stage = c(

    "Initial STAR count matrix",

    "After zero-count filtering",

    "After CPM filtering"
  ),

  genes = c(

    nrow(counts),

    nrow(counts_nonzero),

    nrow(dge_filtered)
  ),

  samples = c(

    ncol(counts),

    ncol(counts_nonzero),

    ncol(dge_filtered)
  ),

  stringsAsFactors = FALSE
)


write.csv(
  filtering_summary,
  file = file.path(
    results_dir,
    "07_TCGA_filtering_summary.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 19. SAVE FILTERED RAW COUNTS
# ------------------------------------------------------------
#
# This is still count-scale data.
# ------------------------------------------------------------

filtered_counts_object <- list(

  created_at =
    Sys.time(),

  source =
    "TCGA-LIHC STAR counts",

  cohort =
    "371 primary tumor samples",

  filtering_method =
    "CPM >= 1 in >=20% of samples",

  CPM_threshold =
    cpm_threshold,

  minimum_fraction =
    minimum_fraction,

  minimum_samples =
    minimum_samples,

  genes_after_filtering =
    nrow(dge_filtered),

  samples =
    ncol(dge_filtered),

  counts =
    dge_filtered$counts
)


saveRDS(
  filtered_counts_object,
  file = file.path(
    processed_dir,
    "TCGA_LIHC_filtered_counts.rds"
  )
)


# ------------------------------------------------------------
# 20. SAVE VOOM EXPRESSION
# ------------------------------------------------------------

voom_expression_object <- list(

  created_at =
    Sys.time(),

  source =
    "TCGA-LIHC STAR counts",

  cohort =
    "371 primary tumor samples",

  filtering =
    "CPM >= 1 in >=20% of samples",

  normalization =
    "edgeR TMM",

  transformation =
    "limma voom",

  n_genes =
    nrow(voom_expression),

  n_samples =
    ncol(voom_expression),

  sample_names =
    colnames(voom_expression),

  gene_ids =
    rownames(voom_expression),

  expression =
    voom_expression
)


saveRDS(
  voom_expression_object,
  file = file.path(
    processed_dir,
    "TCGA_LIHC_voom_expression.rds"
  )
)


# ------------------------------------------------------------
# 21. SAVE COMPLETE VOOM OBJECT
# ------------------------------------------------------------

saveRDS(
  voom_object,
  file = file.path(
    processed_dir,
    "TCGA_LIHC_voom_object.rds"
  )
)


# ------------------------------------------------------------
# 22. SAVE PROCESSING PARAMETERS
# ------------------------------------------------------------

processing_parameters <- list(

  created_at =
    Sys.time(),

  initial_genes =
    nrow(counts),

  zero_count_genes_removed =
    sum(zero_count_genes),

  genes_after_zero_filter =
    nrow(counts_nonzero),

  CPM_threshold =
    cpm_threshold,

  minimum_fraction =
    minimum_fraction,

  minimum_samples =
    minimum_samples,

  genes_after_CPM_filter =
    nrow(dge_filtered),

  samples =
    ncol(voom_expression),

  filtering =
    "CPM >= 1 in >=20% of samples",

  normalization =
    "edgeR TMM",

  transformation =
    "limma voom",

  n163_filtering =
    FALSE,

  maf_required =
    FALSE,

  GSVA_performed =
    FALSE
)


saveRDS(
  processing_parameters,
  file = file.path(
    results_dir,
    "07_TCGA_processing_parameters.rds"
  )
)


# ------------------------------------------------------------
# 23. WRITE PROCESSING REPORT
# ------------------------------------------------------------

report_lines <- c(

  "============================================================",

  "TCGA-LIHC RNA-seq PREPROCESSING REPORT",

  "============================================================",

  "",

  paste(
    "Date:",
    Sys.time()
  ),

  "",

  "DISCOVERY COHORT",

  "------------------------------------------------------------",

  "All TCGA-LIHC primary tumors",

  paste(
    "Samples:",
    ncol(voom_expression)
  ),

  "",

  "GENE FILTERING",

  "------------------------------------------------------------",

  paste(
    "Initial genes:",
    nrow(counts)
  ),

  paste(
    "Zero-count genes removed:",
    sum(zero_count_genes)
  ),

  paste(
    "Genes after zero-count filtering:",
    nrow(counts_nonzero)
  ),

  paste(
    "CPM threshold:",
    cpm_threshold
  ),

  paste(
    "Minimum fraction:",
    minimum_fraction
  ),

  paste(
    "Minimum samples:",
    minimum_samples
  ),

  paste(
    "Genes retained after CPM filter:",
    nrow(dge_filtered)
  ),

  paste(
    "Genes removed by CPM filter:",
    sum(!keep_genes)
  ),

  "",

  "NORMALIZATION",

  "------------------------------------------------------------",

  "edgeR TMM",

  "",

  "TRANSFORMATION",

  "------------------------------------------------------------",

  "limma-voom",

  "Intercept-only design used.",

  "",

  "QUALITY CONTROL",

  "------------------------------------------------------------",

  paste(
    "NA values:",
    sum(
      is.na(
        voom_expression
      )
    )
  ),

  paste(
    "Non-finite values:",
    sum(
      !is.finite(
        voom_expression
      )
    )
  ),

  "",

  "IMPORTANT",

  "------------------------------------------------------------",

  "All 371 primary tumors retained.",

  "No n=163 filtering performed.",

  "No survival-based exclusions performed.",

  "No automatic sample removal performed.",

  "MAF not required.",

  "GSVA not performed.",

  "",

  "NEXT STEP",

  "------------------------------------------------------------",

  "Audit gene identifiers and pathway gene-set coverage."
)


writeLines(
  report_lines,
  con = file.path(
    results_dir,
    "07_TCGA_preprocessing_report.txt"
  )
)


# ------------------------------------------------------------
# 24. CLEAN MEMORY
# ------------------------------------------------------------

rm(
  tcga_object,
  counts,
  counts_nonzero,
  cpm_matrix,
  cpm_filtered,
  dge,
  dge_filtered
)

gc()


# ------------------------------------------------------------
# 25. FINAL STATUS
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("STEP 7 COMPLETED SUCCESSFULLY\n")
cat("============================================================\n\n")


cat(
  "TCGA primary tumor samples: ",
  ncol(voom_expression),
  "\n",
  sep = ""
)


cat(
  "Initial genes: ",
  processing_parameters$initial_genes,
  "\n",
  sep = ""
)


cat(
  "Genes after zero-count removal: ",
  processing_parameters$genes_after_zero_filter,
  "\n",
  sep = ""
)


cat(
  "CPM threshold: >= ",
  processing_parameters$CPM_threshold,
  "\n",
  sep = ""
)


cat(
  "Minimum samples: ",
  processing_parameters$minimum_samples,
  "\n",
  sep = ""
)


cat(
  "Genes after CPM filtering: ",
  processing_parameters$genes_after_CPM_filter,
  "\n\n",
  sep = ""
)


cat(
  "Normalization: edgeR TMM\n"
)


cat(
  "Transformation: limma-voom\n"
)


cat(
  "n=163 filtering: NOT PERFORMED\n"
)


cat(
  "MAF: NOT REQUIRED\n"
)


cat(
  "GSVA: NOT PERFORMED\n\n"
)


cat(
  "Final TCGA expression object:\n",
  file.path(
    processed_dir,
    "TCGA_LIHC_voom_expression.rds"
  ),
  "\n\n",
  sep = ""
)


cat(
  "Next step: gene identifier and pathway gene-set audit.\n"
)
