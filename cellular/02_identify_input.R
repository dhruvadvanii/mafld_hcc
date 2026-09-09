# ============================================================
# GSE122151 EXPERIMENT
# STEP 2 — IDENTIFY INPUT OBJECTS
# ============================================================

cat("\n============================================================\n")
cat("GSE122151 — STEP 2: INPUT OBJECT IDENTIFICATION\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. HELPER FUNCTION
# ------------------------------------------------------------
# Safely describe an object without modifying it.

describe_object <- function(object_name) {

  x <- get(object_name, envir = .GlobalEnv)

  cat("\n------------------------------------------------------------\n")
  cat("OBJECT:", object_name, "\n")
  cat("CLASS :", paste(class(x), collapse = ", "), "\n")

  if (is.matrix(x) || is.data.frame(x)) {
    cat(
      "DIMENSIONS:",
      paste(dim(x), collapse = " x "),
      "\n"
    )
  }

  if (is.list(x)) {
    cat("LIST LENGTH:", length(x), "\n")

    if (!is.null(names(x))) {
      cat(
        "FIRST NAMES:",
        paste(head(names(x), 10), collapse = ", "),
        "\n"
      )
    }
  }

  if (is.factor(x)) {
    cat(
      "LEVELS:",
      paste(levels(x), collapse = ", "),
      "\n"
    )
  }

  cat("------------------------------------------------------------\n")
}


# ------------------------------------------------------------
# 2. INSPECT POTENTIAL GENE-SET OBJECTS
# ------------------------------------------------------------
# We deliberately exclude all objects whose names clearly
# belong to TCGA, GSE76427, GSE122151 legacy analyses,
# or other specific cohorts.

gene_set_objects_to_inspect <- c(
  "gene_set_list",
  "gene_set_objects",
  "final_gene_sets",
  "pathway_gene_sets_raw",
  "pathway_genes",
  "pathway_set",
  "ddr_pathways",
  "ddr_pathway",
  "hr_gmt",
  "nhej_gmt",
  "mmr_pathway_name",
  "dsb_gmt"
)

gene_set_objects_to_inspect <- gene_set_objects_to_inspect[
  gene_set_objects_to_inspect %in% ls(envir = .GlobalEnv)
]

cat("\n============================================================\n")
cat("POTENTIAL GENE-SET OBJECTS\n")
cat("============================================================\n")

if (length(gene_set_objects_to_inspect) == 0) {

  cat("No candidate gene-set objects found.\n")

} else {

  for (obj in gene_set_objects_to_inspect) {
    describe_object(obj)
  }

}


# ------------------------------------------------------------
# 3. INSPECT POTENTIAL EXPRESSION OBJECTS
# ------------------------------------------------------------
# Do NOT inspect or use TCGA/GSE76427 objects here.

expression_objects_to_inspect <- c(
  "current_expression",
  "gse_expression",
  "gse_expr",
  "gse_gsva_expr",
  "gse_symbol_expression",
  "gene_expression",
  "collapsed_expression",
  "final_matrix",
  "candidate_matrix",
  "potential_expression",
  "expression_objects",
  "expression_structure"
)

expression_objects_to_inspect <- expression_objects_to_inspect[
  expression_objects_to_inspect %in% ls(envir = .GlobalEnv)
]

cat("\n============================================================\n")
cat("POTENTIAL EXPRESSION OBJECTS\n")
cat("============================================================\n")

if (length(expression_objects_to_inspect) == 0) {

  cat("No candidate expression objects found.\n")

} else {

  for (obj in expression_objects_to_inspect) {
    describe_object(obj)
  }

}


# ------------------------------------------------------------
# 4. INSPECT POTENTIAL SAMPLE / GROUP OBJECTS
# ------------------------------------------------------------

metadata_objects_to_inspect <- c(
  "biological_group",
  "group_assignment",
  "group_defining_data",
  "sample_manifest",
  "sample_manifest_tumor",
  "sample_names",
  "sample_columns",
  "sample_summary",
  "sample_qc"
)

metadata_objects_to_inspect <- metadata_objects_to_inspect[
  metadata_objects_to_inspect %in% ls(envir = .GlobalEnv)
]

cat("\n============================================================\n")
cat("POTENTIAL SAMPLE / GROUP OBJECTS\n")
cat("============================================================\n")

if (length(metadata_objects_to_inspect) == 0) {

  cat("No candidate metadata objects found.\n")

} else {

  for (obj in metadata_objects_to_inspect) {
    describe_object(obj)
  }

}


# ------------------------------------------------------------
# 5. INSPECT GENE-SET CONTENT
# ------------------------------------------------------------
# For list-like gene-set objects, show pathway names and
# gene counts. Nothing is modified.

cat("\n============================================================\n")
cat("GENE-SET CONTENT AUDIT\n")
cat("============================================================\n")

gene_set_summary <- NULL

for (obj in gene_set_objects_to_inspect) {

  x <- get(obj, envir = .GlobalEnv)

  if (is.list(x) && length(x) > 0) {

    is_gene_set_list <- all(
      vapply(
        x,
        function(z) is.character(z) || is.factor(z),
        logical(1)
      )
    )

    if (is_gene_set_list) {

      gene_counts <- vapply(
        x,
        function(z) length(unique(as.character(z))),
        integer(1)
      )

      gene_set_summary <- rbind(
        gene_set_summary,
        data.frame(
          object = obj,
          n_gene_sets = length(x),
          min_genes = min(gene_counts),
          median_genes = median(gene_counts),
          max_genes = max(gene_counts),
          stringsAsFactors = FALSE
        )
      )

      cat("\nOBJECT:", obj, "\n")
      cat("Number of gene sets:", length(x), "\n")
      cat("First pathway names:\n")
      print(head(names(x), 20))
      cat("Gene counts:\n")
      print(summary(gene_counts))
    }
  }
}


# ------------------------------------------------------------
# 6. SAVE THE STRUCTURAL AUDIT
# ------------------------------------------------------------

if (!is.null(gene_set_summary)) {

  write.csv(
    gene_set_summary,
    file = file.path(
      results_dir,
      "02_gene_set_object_audit.csv"
    ),
    row.names = FALSE
  )

}


# ------------------------------------------------------------
# 7. SAVE CANDIDATE OBJECT NAMES
# ------------------------------------------------------------

candidate_object_audit <- list(
  gene_sets = gene_set_objects_to_inspect,
  expression = expression_objects_to_inspect,
  metadata = metadata_objects_to_inspect
)

saveRDS(
  candidate_object_audit,
  file = file.path(
    results_dir,
    "02_candidate_input_objects.rds"
  )
)


# ------------------------------------------------------------
# 8. FINAL MESSAGE
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP 2 COMPLETE\n")
cat("============================================================\n")
cat("\nNo input object has been modified.\n")
cat("No TCGA object has been used.\n")
cat("No GSE76427 object has been used.\n")
cat("No legacy GSE122151 object has been used.\n")
cat("\nThe next step will select the correct inputs based on\n")
cat("the structural audit printed above.\n")
cat("============================================================\n")
