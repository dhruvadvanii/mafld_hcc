# ============================================================
# GSE122151 — STEP 4: SAMPLE QUALITY CONTROL
# ============================================================
#
# Input:
#   data/processed/GSE122151/GSE122151_log2_FPKM_plus1.rds
#   data/processed/GSE122151/sample_metadata.rds
#
# QC analyses:
#   1. Sample-level expression summary
#   2. Expression distribution
#   3. Spearman sample correlation
#   4. PCA
#   5. Hierarchical clustering
#
# Every figure is generated from a corresponding saved
# CSV file in results/GSE122151/.
#
# ============================================================

# ------------------------------------------------------------
# STEP 0 — Clean R environment
# ------------------------------------------------------------

rm(list = ls())
gc()

options(
  stringsAsFactors = FALSE,
  scipen = 999,
  warn = 1
)

# ------------------------------------------------------------
# STEP 1 — Project directories
# ------------------------------------------------------------

project_root <- getwd()

processed_dir <- file.path(
  project_root,
  "data/processed/GSE122151"
)

results_dir <- file.path(
  project_root,
  "results/GSE122151"
)

figures_dir <- file.path(
  project_root,
  "figures/GSE122151"
)

dir.create(
  results_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  figures_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# STEP 2 — Check required input files
# ------------------------------------------------------------

required_files <- c(
  file.path(
    processed_dir,
    "GSE122151_log2_FPKM_plus1.rds"
  ),
  file.path(
    processed_dir,
    "sample_metadata.rds"
  )
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {

  stop(
    "The following required input files are missing:\n",
    paste(
      missing_files,
      collapse = "\n"
    ),
    "\n\nRun Step 3 before running Step 4."
  )
}

# ------------------------------------------------------------
# STEP 3 — Load expression matrix
# ------------------------------------------------------------

expr_log2 <- readRDS(
  file.path(
    processed_dir,
    "GSE122151_log2_FPKM_plus1.rds"
  )
)

# ------------------------------------------------------------
# STEP 4 — Load sample metadata
# ------------------------------------------------------------

sample_metadata <- readRDS(
  file.path(
    processed_dir,
    "sample_metadata.rds"
  )
)

# ------------------------------------------------------------
# STEP 5 — Validate expression matrix
# ------------------------------------------------------------

expr_log2 <- as.matrix(
  expr_log2
)

storage.mode(expr_log2) <- "numeric"

if (is.null(rownames(expr_log2))) {

  stop(
    "Expression matrix does not have gene identifiers as row names."
  )
}

if (is.null(colnames(expr_log2))) {

  stop(
    "Expression matrix does not have sample identifiers as column names."
  )
}

if (anyDuplicated(rownames(expr_log2)) > 0) {

  stop(
    "Expression matrix contains duplicated gene identifiers."
  )
}

if (anyDuplicated(colnames(expr_log2)) > 0) {

  stop(
    "Expression matrix contains duplicated sample identifiers."
  )
}

if (any(!is.finite(expr_log2))) {

  stop(
    "Expression matrix contains NA, NaN, Inf, or -Inf values."
  )
}

# ------------------------------------------------------------
# STEP 6 — Validate metadata structure
# ------------------------------------------------------------

required_metadata_columns <- c(
  "GSM",
  "sample_name",
  "condition",
  "replicate"
)

missing_metadata_columns <- setdiff(
  required_metadata_columns,
  colnames(sample_metadata)
)

if (length(missing_metadata_columns) > 0) {

  stop(
    "Sample metadata is missing required columns:\n",
    paste(
      missing_metadata_columns,
      collapse = ", "
    )
  )
}

# ------------------------------------------------------------
# STEP 7 — Validate metadata uniqueness
# ------------------------------------------------------------

if (
  anyDuplicated(
    sample_metadata$GSM
  ) > 0
) {

  stop(
    "Duplicate GSM identifiers found in sample metadata."
  )
}

if (
  anyDuplicated(
    sample_metadata$sample_name
  ) > 0
) {

  stop(
    "Duplicate sample_name identifiers found in sample metadata."
  )
}

# ------------------------------------------------------------
# STEP 8 — Check expected six samples
# ------------------------------------------------------------

expected_samples <- c(
  "HepaRG1_control_S7",
  "HepaRG2_control_S10",
  "HepaRG3_control_S11",
  "HepaRG5_FA_S12",
  "HepaRG6_FA_S13",
  "HepaRG7_FA_S14"
)

if (
  !identical(
    colnames(expr_log2),
    expected_samples
  )
) {

  stop(
    "The expression matrix sample order does not match ",
    "the expected GSE122151 sample order.\n\n",
    "Expected:\n",
    paste(
      expected_samples,
      collapse = "\n"
    ),
    "\n\nFound:\n",
    paste(
      colnames(expr_log2),
      collapse = "\n"
    )
  )
}

# ------------------------------------------------------------
# STEP 9 — Reorder metadata to expression matrix
# ------------------------------------------------------------
#
# IMPORTANT:
# We use sample_name for matching.
#
# We do NOT use rownames_to_column().
# This avoids the previous duplicated-GSM error.
#
# ------------------------------------------------------------

metadata_match <- match(
  colnames(expr_log2),
  sample_metadata$sample_name
)

if (anyNA(metadata_match)) {

  stop(
    "One or more expression-matrix samples could not be ",
    "matched to sample metadata."
  )
}

sample_metadata <- sample_metadata[
  metadata_match,
  ,
  drop = FALSE
]

# ------------------------------------------------------------
# STEP 10 — Final metadata/order validation
# ------------------------------------------------------------

if (
  !identical(
    sample_metadata$sample_name,
    colnames(expr_log2)
  )
) {

  stop(
    "Metadata sample order does not match expression matrix."
  )
}

if (
  !identical(
    sample_metadata$condition,
    c(
      "Control",
      "Control",
      "Control",
      "Lipid_Loaded",
      "Lipid_Loaded",
      "Lipid_Loaded"
    )
  )
) {

  stop(
    "Unexpected sample condition assignment."
  )
}

# ------------------------------------------------------------
# STEP 11 — Save validated metadata
# ------------------------------------------------------------

write.csv(
  sample_metadata,
  file = file.path(
    results_dir,
    "validated_sample_metadata.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)

# ------------------------------------------------------------
# STEP 12 — Sample-level expression QC
# ------------------------------------------------------------

sample_qc <- data.frame(

  GSM = sample_metadata$GSM,

  sample_name = sample_metadata$sample_name,

  condition = sample_metadata$condition,

  replicate = sample_metadata$replicate,

  mean_expression =
    colMeans(expr_log2),

  median_expression =
    apply(
      expr_log2,
      2,
      median
    ),

  sd_expression =
    apply(
      expr_log2,
      2,
      sd
    ),

  min_expression =
    apply(
      expr_log2,
      2,
      min
    ),

  max_expression =
    apply(
      expr_log2,
      2,
      max
    ),

  detected_genes =
    colSums(
      expr_log2 > 0
    ),

  genes_above_1 =
    colSums(
      expr_log2 > 1
    ),

  genes_above_2 =
    colSums(
      expr_log2 > 2
    ),

  stringsAsFactors = FALSE
)

write.csv(
  sample_qc,
  file = file.path(
    results_dir,
    "sample_QC.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)

# ------------------------------------------------------------
# STEP 13 — Expression distribution data
# ------------------------------------------------------------

expression_distribution_data <- data.frame(
  gene = rep(
    rownames(expr_log2),
    times = ncol(expr_log2)
  ),
  sample_name = rep(
    colnames(expr_log2),
    each = nrow(expr_log2)
  ),
  expression = as.vector(
    expr_log2
  ),
  stringsAsFactors = FALSE
)

expression_distribution_data <- merge(
  expression_distribution_data,
  sample_metadata[
    ,
    c(
      "GSM",
      "sample_name",
      "condition",
      "replicate"
    )
  ],
  by = "sample_name",
  all.x = TRUE,
  sort = FALSE
)

# Reorder columns

expression_distribution_data <- expression_distribution_data[
  ,
  c(
    "gene",
    "GSM",
    "sample_name",
    "condition",
    "replicate",
    "expression"
  )
]

write.csv(
  expression_distribution_data,
  file = file.path(
    results_dir,
    "FIG_expression_distribution_data.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)

# ------------------------------------------------------------
# STEP 14 — Expression distribution figure
# ------------------------------------------------------------

p_expression <- ggplot2::ggplot(
  expression_distribution_data,
  ggplot2::aes(
    x = sample_name,
    y = expression
  )
) +

  ggplot2::geom_boxplot(
    outlier.size = 0.1
  ) +

  ggplot2::theme_bw() +

  ggplot2::theme(
    axis.text.x =
      ggplot2::element_text(
        angle = 45,
        hjust = 1
      )
  ) +

  ggplot2::labs(
    title =
      "GSE122151 expression distribution",

    x =
      NULL,

    y =
      "log2(FPKM + 1)"
  )

ggplot2::ggsave(
  filename = file.path(
    figures_dir,
    "Figure_S1_expression_distribution.png"
  ),
  plot = p_expression,
  width = 9,
  height = 6,
  dpi = 300
)

# ------------------------------------------------------------
# STEP 15 — Spearman correlation matrix
# ------------------------------------------------------------

sample_cor <- cor(
  expr_log2,
  method = "spearman",
  use = "pairwise.complete.obs"
)

if (
  any(!is.finite(sample_cor))
) {

  stop(
    "Non-finite values detected in sample correlation matrix."
  )
}

write.csv(
  sample_cor,
  file = file.path(
    results_dir,
    "FIG_sample_correlation_matrix.csv"
  ),
  row.names = TRUE,
  quote = FALSE
)

# ------------------------------------------------------------
# STEP 16 — Correlation figure
# ------------------------------------------------------------

png(
  filename = file.path(
    figures_dir,
    "Figure_S2_sample_correlation.png"
  ),
  width = 2200,
  height = 2000,
  res = 300
)

pheatmap::pheatmap(
  sample_cor,

  clustering_distance_rows =
    "correlation",

  clustering_distance_cols =
    "correlation",

  main =
    "Spearman correlation between samples"
)

dev.off()

# ------------------------------------------------------------
# STEP 17 — PCA
# ------------------------------------------------------------

pca <- stats::prcomp(
  t(expr_log2),
  center = TRUE,
  scale. = FALSE
)

# ------------------------------------------------------------
# Calculate PCA variance
# ------------------------------------------------------------

pca_variance <- (
  pca$sdev ^ 2
) / sum(
  pca$sdev ^ 2
)

# ------------------------------------------------------------
# PCA data
# ------------------------------------------------------------

pca_data <- data.frame(

  sample_name =
    rownames(pca$x),

  PC1 =
    pca$x[, "PC1"],

  PC2 =
    pca$x[, "PC2"],

  PC3 =
    if (
      "PC3" %in% colnames(pca$x)
    ) {
      pca$x[, "PC3"]
    } else {
      NA_real_
    },

  stringsAsFactors = FALSE
)

pca_data <- merge(
  pca_data,
  sample_metadata[
    ,
    c(
      "GSM",
      "sample_name",
      "condition",
      "replicate"
    )
  ],
  by = "sample_name",
  all.x = TRUE,
  sort = FALSE
)

# ------------------------------------------------------------
# Reorder columns
# ------------------------------------------------------------

pca_data <- pca_data[
  ,
  c(
    "GSM",
    "sample_name",
    "condition",
    "replicate",
    "PC1",
    "PC2",
    "PC3"
  )
]

write.csv(
  pca_data,
  file = file.path(
    results_dir,
    "FIG_PCA_data.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)

# ------------------------------------------------------------
# PCA variance table
# ------------------------------------------------------------

pca_variance_data <- data.frame(

  principal_component =
    paste0(
      "PC",
      seq_along(pca_variance)
    ),

  variance_fraction =
    pca_variance,

  variance_percent =
    100 * pca_variance,

  cumulative_variance_percent =
    100 * cumsum(
      pca_variance
    ),

  stringsAsFactors = FALSE
)

write.csv(
  pca_variance_data,
  file = file.path(
    results_dir,
    "PCA_variance_explained.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)

# ------------------------------------------------------------
# PCA figure
# ------------------------------------------------------------

p_pca <- ggplot2::ggplot(
  pca_data,
  ggplot2::aes(
    x = PC1,
    y = PC2,
    shape = condition
  )
) +

  ggplot2::geom_point(
    size = 4
  ) +

  ggplot2::geom_text(
    ggplot2::aes(
      label = sample_name
    ),
    vjust = -0.8,
    size = 3
  ) +

  ggplot2::theme_bw() +

  ggplot2::labs(

    title =
      "Principal component analysis",

    x =
      paste0(
        "PC1 (",
        round(
          pca_variance[1] * 100,
          1
        ),
        "%)"
      ),

    y =
      paste0(
        "PC2 (",
        round(
          pca_variance[2] * 100,
          1
        ),
        "%)"
      ),

    shape =
      "Condition"
  )

ggplot2::ggsave(
  filename = file.path(
    figures_dir,
    "Figure_S3_PCA.png"
  ),
  plot = p_pca,
  width = 8,
  height = 6,
  dpi = 300
)

# ------------------------------------------------------------
# STEP 18 — Hierarchical clustering
# ------------------------------------------------------------

sample_distance <- dist(
  t(expr_log2),
  method = "euclidean"
)

hc <- hclust(
  sample_distance,
  method = "complete"
)

# ------------------------------------------------------------
# Save distance matrix
# ------------------------------------------------------------

sample_distance_matrix <- as.matrix(
  sample_distance
)

write.csv(
  sample_distance_matrix,
  file = file.path(
    results_dir,
    "FIG_sample_euclidean_distance_matrix.csv"
  ),
  row.names = TRUE,
  quote = FALSE
)

# ------------------------------------------------------------
# Save clustering order
# ------------------------------------------------------------

hc_order_data <- data.frame(

  cluster_order =
    seq_along(hc$order),

  sample_name =
    hc$labels[hc$order],

  GSM =
    sample_metadata$GSM[
      match(
        hc$labels[hc$order],
        sample_metadata$sample_name
      )
    ],

  condition =
    sample_metadata$condition[
      match(
        hc$labels[hc$order],
        sample_metadata$sample_name
      )
    ],

  stringsAsFactors = FALSE
)

write.csv(
  hc_order_data,
  file = file.path(
    results_dir,
    "FIG_hierarchical_clustering_order.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)

# ------------------------------------------------------------
# Hierarchical clustering figure
# ------------------------------------------------------------

png(
  filename = file.path(
    figures_dir,
    "Figure_S4_hierarchical_clustering.png"
  ),
  width = 2200,
  height = 1800,
  res = 300
)

plot(
  hc,
  main =
    "Hierarchical clustering of GSE122151 samples",
  xlab = "",
  sub = "",
  cex = 0.9
)

dev.off()

# ------------------------------------------------------------
# STEP 19 — Combined QC summary
# ------------------------------------------------------------

qc_summary <- data.frame(

  metric = c(
    "Number of genes",
    "Number of samples",
    "Control samples",
    "Lipid-loaded samples",
    "Minimum expression value",
    "Maximum expression value",
    "Mean sample expression",
    "Median sample expression"
  ),

  value = c(

    nrow(expr_log2),

    ncol(expr_log2),

    sum(
      sample_metadata$condition ==
        "Control"
    ),

    sum(
      sample_metadata$condition ==
        "Lipid_Loaded"
    ),

    min(
      expr_log2
    ),

    max(
      expr_log2
    ),

    mean(
      expr_log2
    ),

    median(
      expr_log2
    )
  ),

  stringsAsFactors = FALSE
)

write.csv(
  qc_summary,
  file = file.path(
    results_dir,
    "QC_summary.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)

# ------------------------------------------------------------
# STEP 20 — Save complete QC object
# ------------------------------------------------------------

qc_object <- list(

  expression_matrix =
    expr_log2,

  sample_metadata =
    sample_metadata,

  sample_qc =
    sample_qc,

  correlation_matrix =
    sample_cor,

  pca =
    pca,

  pca_data =
    pca_data,

  pca_variance =
    pca_variance_data,

  hierarchical_clustering =
    hc,

  distance_matrix =
    sample_distance_matrix

)

saveRDS(
  qc_object,
  file = file.path(
    processed_dir,
    "GSE122151_QC_object.rds"
  )
)

# ------------------------------------------------------------
# STEP 21 — Final validation
# ------------------------------------------------------------

if (
  nrow(expr_log2) < 1000
) {

  stop(
    "QC validation failed: fewer than 1,000 genes."
  )
}

if (
  ncol(expr_log2) != 6
) {

  stop(
    "QC validation failed: expected six samples."
  )
}

if (
  sum(
    sample_metadata$condition ==
      "Control"
  ) != 3
) {

  stop(
    "QC validation failed: expected three control samples."
  )
}

if (
  sum(
    sample_metadata$condition ==
      "Lipid_Loaded"
  ) != 3
) {

  stop(
    "QC validation failed: expected three lipid-loaded samples."
  )
}

# ------------------------------------------------------------
# STEP 22 — Completion message
# ------------------------------------------------------------

cat(
  "\n============================================================\n"
)

cat(
  "GSE122151 STEP 4 — SAMPLE QC COMPLETE\n"
)

cat(
  "============================================================\n\n"
)

cat(
  "Genes:",
  nrow(expr_log2),
  "\n"
)

cat(
  "Samples:",
  ncol(expr_log2),
  "\n"
)

cat(
  "Control:",
  sum(
    sample_metadata$condition ==
      "Control"
  ),
  "\n"
)

cat(
  "Lipid-loaded:",
  sum(
    sample_metadata$condition ==
      "Lipid_Loaded"
  ),
  "\n\n"
)

cat(
  "QC outputs saved to:\n",
  results_dir,
  "\n\n"
)

cat(
  "QC figures saved to:\n",
  figures_dir,
  "\n\n"
)

cat(
  "No duplicated GSM column construction was used.\n"
)

cat(
  "============================================================\n"
)
