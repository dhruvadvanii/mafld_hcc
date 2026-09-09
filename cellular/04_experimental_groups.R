# ============================================================
# GSE122151 EXPERIMENT
# STEP 5 — DEFINE EXPERIMENTAL SAMPLES AND GROUPS
# ============================================================

cat("\n============================================================\n")
cat("GSE122151 — STEP 5: SAMPLE AND GROUP DEFINITION\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. DEFINE THE SIX EXPERIMENTAL SAMPLES
# ------------------------------------------------------------
# These are the six HepaRG samples present in the supplied
# GSE122151 FPKM table.
#
# Control:
#   HepaRG1_control_S7
#   HepaRG2_control_S10
#   HepaRG3_control_S11
#
# Lipid-loaded / FA:
#   HepaRG5_FA_S12
#   HepaRG6_FA_S13
#   HepaRG7_FA_S14

control_samples <- c(
  "HepaRG1_control_S7",
  "HepaRG2_control_S10",
  "HepaRG3_control_S11"
)

lipid_loaded_samples <- c(
  "HepaRG5_FA_S12",
  "HepaRG6_FA_S13",
  "HepaRG7_FA_S14"
)

gse122151_samples <- c(
  control_samples,
  lipid_loaded_samples
)


# ------------------------------------------------------------
# 2. VERIFY THAT ALL SIX SAMPLES EXIST
# ------------------------------------------------------------

missing_samples <- setdiff(
  gse122151_samples,
  colnames(gse122151_raw)
)

if (length(missing_samples) > 0) {

  stop(
    "The following expected GSE122151 samples were not found:\n",
    paste(missing_samples, collapse = "\n")
  )

}


# ------------------------------------------------------------
# 3. VERIFY THAT THE SIX SAMPLES ARE UNIQUE
# ------------------------------------------------------------

if (length(unique(gse122151_samples)) != 6) {

  stop(
    "The GSE122151 sample list does not contain exactly six unique samples."
  )

}


# ------------------------------------------------------------
# 4. CREATE SAMPLE METADATA
# ------------------------------------------------------------

gse122151_sample_metadata <- data.frame(
  sample = gse122151_samples,
  group = c(
    rep("Control", length(control_samples)),
    rep("Lipid_Loaded", length(lipid_loaded_samples))
  ),
  stringsAsFactors = FALSE
)

rownames(gse122151_sample_metadata) <-
  gse122151_sample_metadata$sample


# ------------------------------------------------------------
# 5. CONVERT GROUP TO FACTOR
# ------------------------------------------------------------
# Control is deliberately specified as the reference group.

gse122151_sample_metadata$group <- factor(
  gse122151_sample_metadata$group,
  levels = c(
    "Control",
    "Lipid_Loaded"
  )
)


# ------------------------------------------------------------
# 6. VERIFY GROUP COUNTS
# ------------------------------------------------------------

group_counts <- table(
  gse122151_sample_metadata$group
)

cat("============================================================\n")
cat("GROUP COUNTS\n")
cat("============================================================\n\n")

print(group_counts)

if (
  !identical(
    as.integer(group_counts),
    c(3L, 3L)
  )
) {

  stop(
    "Expected exactly 3 Control and 3 Lipid_Loaded samples."
  )

}


# ------------------------------------------------------------
# 7. DISPLAY SAMPLE METADATA
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("SAMPLE METADATA\n")
cat("============================================================\n\n")

print(
  gse122151_sample_metadata
)


# ------------------------------------------------------------
# 8. EXTRACT ONLY THE SIX EXPERIMENTAL SAMPLES
# ------------------------------------------------------------
# Gene annotation columns are deliberately retained in the
# raw object but excluded from the expression matrix.

gene_annotation_columns <- c(
  "gene_name",
  "gene_type",
  "chr",
  "start",
  "end",
  "strand",
  "length"
)

missing_annotation_columns <- setdiff(
  gene_annotation_columns,
  colnames(gse122151_raw)
)

if (length(missing_annotation_columns) > 0) {

  stop(
    "Expected gene annotation columns are missing:\n",
    paste(
      missing_annotation_columns,
      collapse = "\n"
    )
  )

}


gse122151_expression_raw <- gse122151_raw[
  ,
  gse122151_samples,
  drop = FALSE
]


# ------------------------------------------------------------
# 9. EXTRACT GENE IDENTIFIERS
# ------------------------------------------------------------

gse122151_gene_ids <- as.character(
  gse122151_raw$gene_name
)


# ------------------------------------------------------------
# 10. VERIFY EXPRESSION MATRIX DIMENSIONS
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("EXTRACTED EXPRESSION MATRIX\n")
cat("============================================================\n\n")

cat(
  "Number of genes/rows:",
  nrow(gse122151_expression_raw),
  "\n"
)

cat(
  "Number of samples:",
  ncol(gse122151_expression_raw),
  "\n"
)

cat("\nSample order:\n")
print(colnames(gse122151_expression_raw))


if (ncol(gse122151_expression_raw) != 6) {

  stop(
    "The extracted GSE122151 expression matrix does not contain exactly six samples."
  )

}


# ------------------------------------------------------------
# 11. VERIFY SAMPLE ORDER AGAINST METADATA
# ------------------------------------------------------------

if (
  !identical(
    colnames(gse122151_expression_raw),
    rownames(gse122151_sample_metadata)
  )
) {

  stop(
    "Expression matrix sample order does not match sample metadata."
  )

}


# ------------------------------------------------------------
# 12. SAVE SAMPLE METADATA
# ------------------------------------------------------------

write.csv(
  gse122151_sample_metadata,
  file = file.path(
    results_dir,
    "05_GSE122151_sample_metadata.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 13. SAVE RAW SIX-SAMPLE EXPRESSION MATRIX
# ------------------------------------------------------------

saveRDS(
  gse122151_expression_raw,
  file = file.path(
    results_dir,
    "05_GSE122151_expression_raw_six_samples.rds"
  )
)


# ------------------------------------------------------------
# 14. SAVE GENE IDENTIFIERS
# ------------------------------------------------------------

write.csv(
  data.frame(
    gene_name = gse122151_gene_ids,
    stringsAsFactors = FALSE
  ),
  file = file.path(
    results_dir,
    "05_GSE122151_gene_identifiers.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 15. FINAL MESSAGE
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP 5 COMPLETE\n")
cat("============================================================\n\n")

cat("Six experimental samples successfully defined:\n")
cat("  Control      :", length(control_samples), "samples\n")
cat("  Lipid_Loaded :", length(lipid_loaded_samples), "samples\n\n")

cat("The expression matrix has been restricted to these six\n")
cat("samples only.\n\n")

cat("No log transformation was performed.\n")
cat("No gene filtering was performed.\n")
cat("No duplicate genes were collapsed.\n")
cat("No GSVA was performed.\n")
cat("No statistical comparison was performed.\n")

cat("\n============================================================\n")
