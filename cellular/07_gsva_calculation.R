# ============================================================
# STEP 9 — GSVA FOR ALL 17 GENE SETS
# ============================================================

library(GSVA)

# ------------------------------------------------------------
# 1. Prepare expression matrix
# ------------------------------------------------------------

expr <- as.matrix(gse122151_expression_log2)

storage.mode(expr) <- "numeric"

rownames(expr) <- toupper(trimws(rownames(expr)))

# Remove duplicated gene symbols
expr <- expr[!duplicated(rownames(expr)), , drop = FALSE]


# ------------------------------------------------------------
# 2. Prepare all 17 gene sets
# ------------------------------------------------------------

gene_sets <- gene_set_list

gene_sets <- lapply(gene_sets, function(x) {
  
  x <- toupper(trimws(x))
  x <- x[!is.na(x) & x != ""]
  
  unique(x)
})


# Match gene sets to genes present in the expression matrix
gene_sets <- lapply(gene_sets, function(x) {
  intersect(x, rownames(expr))
})


# Keep pathway names
names(gene_sets) <- names(gene_set_list)


# ------------------------------------------------------------
# 3. Remove only completely empty pathways
# ------------------------------------------------------------

gene_sets <- gene_sets[lengths(gene_sets) > 0]


# ------------------------------------------------------------
# 4. Run GSVA
# ------------------------------------------------------------

gsva_param <- gsvaParam(
  exprData = expr,
  geneSets = gene_sets,
  kcdf = "Gaussian",
  minSize = 1,
  maxSize = Inf
)

gsva_result <- gsva(
  gsva_param,
  verbose = TRUE
)


# ------------------------------------------------------------
# 5. Extract GSVA scores
# ------------------------------------------------------------

# In the installed GSVA version, the result is already a matrix
gse122151_gsva_matrix <- as.matrix(gsva_result)


# ------------------------------------------------------------
# 6. Convert to data frame
# ------------------------------------------------------------

gse122151_gsva_scores <- as.data.frame(
  gse122151_gsva_matrix,
  check.names = FALSE
)

gse122151_gsva_scores$pathway <- rownames(
  gse122151_gsva_scores
)

# Put pathway first
gse122151_gsva_scores <- gse122151_gsva_scores[
  c(
    "pathway",
    setdiff(
      names(gse122151_gsva_scores),
      "pathway"
    )
  )
]


# ------------------------------------------------------------
# 7. Display results
# ------------------------------------------------------------

print(gse122151_gsva_scores)

cat(
  "\nGSVA completed successfully.\n",
  "Pathways:", nrow(gse122151_gsva_scores), "\n",
  "Samples:", ncol(gse122151_gsva_scores) - 1, "\n"
)


# ------------------------------------------------------------
# 8. Save GSVA results
# ------------------------------------------------------------

results_dir <- file.path(
  project_root,
  "results",
  "GSE122151"
)

if (!dir.exists(results_dir)) {
  dir.create(
    results_dir,
    recursive = TRUE
  )
}


write.csv(
  gse122151_gsva_scores,
  file.path(
    results_dir,
    "09_GSE122151_GSVA_all_17_pathways.csv"
  ),
  row.names = FALSE
)


saveRDS(
  gse122151_gsva_scores,
  file.path(
    results_dir,
    "09_GSE122151_GSVA_all_17_pathways.rds"
  )
)


cat(
  "\nSaved:\n",
  file.path(
    results_dir,
    "09_GSE122151_GSVA_all_17_pathways.csv"
  ),
  "\n"
)
