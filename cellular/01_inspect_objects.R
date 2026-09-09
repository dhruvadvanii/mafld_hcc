###############################################################################
# MAFLD-HCC Manuscript Pipeline
# GSE122151 Experimental Analysis
#
# STEP 1 — Project Setup and Object Audit
#
# Objective:
#   Establish a reproducible analysis environment and verify the existing
#   GSE122151 expression and gene-set objects before performing any analysis.
#
# Experimental design:
#   3 Control samples
#   3 Lipid_Loaded samples
#
# Primary analysis:
#   Prespecified DDR gene sets from downloaded GMT files
#   Control vs Lipid_Loaded
#
# No pathway selection is performed in this step.
###############################################################################


# =============================================================================
# 1. PROJECT DIRECTORY
# =============================================================================

# Set this to the root directory of the MAFLD-HCC project on the machine
# running the analysis.

project_root <- "C:/Users/Dhruv Advani/OneDrive/Desktop/AIIMS/Research/MAFLD_HCC_Pipeline"

if (!dir.exists(project_root)) {
  stop(
    "Project directory does not exist: ",
    project_root,
    "\nPlease update 'project_root' to the root of the MAFLD_HCC project."
  )
}

setwd(project_root)

cat("Project root:\n", getwd(), "\n\n")


# =============================================================================
# 2. REQUIRED DIRECTORIES
# =============================================================================

required_dirs <- c(
  "data",
  "data/reference",
  "data/reference/gene_sets",
  "scripts",
  "results",
  "figures"
)

missing_dirs <- required_dirs[!dir.exists(required_dirs)]

if (length(missing_dirs) > 0) {
  stop(
    "The following required project directories are missing:\n",
    paste(missing_dirs, collapse = "\n")
  )
}

cat("Required project directories verified.\n\n")


# =============================================================================
# 3. REQUIRED R PACKAGES
# =============================================================================

required_packages <- c(
  "GSVA",
  "limma",
  "edgeR",
  "matrixStats"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0) {
  stop(
    "The following required R packages are not installed:\n",
    paste(missing_packages, collapse = "\n"),
    "\nInstall them before continuing."
  )
}

suppressPackageStartupMessages({
  library(GSVA)
  library(limma)
  library(edgeR)
  library(matrixStats)
})

cat("Required R packages loaded successfully.\n\n")


# =============================================================================
# 4. CREATE STEP-SPECIFIC RESULT DIRECTORY
# =============================================================================

step_dir <- file.path(
  "results",
  "GSE122151"
)

if (!dir.exists(step_dir)) {
  dir.create(step_dir, recursive = TRUE)
}

cat("Step-specific result directory:\n", step_dir, "\n\n")


# =============================================================================
# 5. SAVE R SESSION INFORMATION
# =============================================================================

capture.output(
  sessionInfo(),
  file = file.path(
    step_dir,
    "01_sessionInfo.txt"
  )
)

cat("R session information saved.\n\n")


# =============================================================================
# 6. IDENTIFY EXISTING GSE122151 OBJECTS
# =============================================================================

# The existing analysis workspace contains objects from previous analyses.
# We do not assume that every similarly named object belongs to the current
# experiment. Instead, we explicitly identify the objects relevant to GSE122151.

existing_objects <- ls(envir = .GlobalEnv)

gse122151_object_names <- existing_objects[
  grepl("^gse12_", existing_objects)
]

if (length(gse122151_object_names) == 0) {
  stop(
    "No objects beginning with 'gse12_' were found in the current R session.\n",
    "Load the existing GSE122151 analysis workspace before continuing."
  )
}

cat("GSE122151-related objects detected:\n")
cat(paste(gse122151_object_names, collapse = "\n"))
cat("\n\n")


# =============================================================================
# 7. OBJECT INVENTORY
# =============================================================================

object_inventory <- data.frame(
  object = gse122151_object_names,
  class = vapply(
    gse122151_object_names,
    function(x) paste(class(get(x)), collapse = "; "),
    character(1)
  ),
  mode = vapply(
    gse122151_object_names,
    function(x) mode(get(x)),
    character(1)
  ),
  stringsAsFactors = FALSE
)

print(object_inventory)

write.csv(
  object_inventory,
  file = file.path(
    step_dir,
    "01_GSE122151_object_inventory.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# 8. VERIFY PRIMARY EXPRESSION OBJECT
# =============================================================================

if (!exists("gse12_expression", envir = .GlobalEnv)) {
  stop(
    "The expected GSE122151 expression object 'gse12_expression' was not found."
  )
}

gse12_expression <- get(
  "gse12_expression",
  envir = .GlobalEnv
)

if (!is.matrix(gse12_expression) &&
    !is.data.frame(gse12_expression)) {
  stop(
    "'gse12_expression' must be a matrix or data.frame."
  )
}

gse12_expression <- as.matrix(gse12_expression)

storage.mode(gse12_expression) <- "numeric"


# =============================================================================
# 9. EXPRESSION MATRIX STRUCTURAL QC
# =============================================================================

if (is.null(rownames(gse12_expression))) {
  stop(
    "'gse12_expression' does not have gene identifiers as row names."
  )
}

if (is.null(colnames(gse12_expression))) {
  stop(
    "'gse12_expression' does not have sample identifiers as column names."
  )
}

if (nrow(gse12_expression) == 0 ||
    ncol(gse12_expression) == 0) {
  stop(
    "'gse12_expression' is empty."
  )
}

if (anyDuplicated(rownames(gse12_expression)) > 0) {
  stop(
    "'gse12_expression' contains duplicated gene identifiers."
  )
}

if (anyDuplicated(colnames(gse12_expression)) > 0) {
  stop(
    "'gse12_expression' contains duplicated sample identifiers."
  )
}

if (any(!is.finite(gse12_expression))) {
  stop(
    "'gse12_expression' contains NA, NaN, Inf, or -Inf values."
  )
}

cat("GSE122151 expression matrix dimensions:\n")
cat(
  "Genes:   ", nrow(gse12_expression), "\n",
  "Samples: ", ncol(gse12_expression), "\n\n",
  sep = ""
)

cat("Sample identifiers:\n")
cat(paste(colnames(gse12_expression), collapse = "\n"))
cat("\n\n")


# =============================================================================
# 10. VERIFY THE SIX-SAMPLE EXPERIMENTAL DESIGN
# =============================================================================

if (ncol(gse12_expression) != 6) {
  stop(
    "The GSE122151 expression matrix does not contain exactly 6 samples.\n",
    "Expected: 3 Control + 3 Lipid_Loaded."
  )
}

if (!exists("gse12_targets", envir = .GlobalEnv)) {
  stop(
    "The expected GSE122151 sample metadata object 'gse12_targets' ",
    "was not found."
  )
}

gse12_targets <- get(
  "gse12_targets",
  envir = .GlobalEnv
)

if (!is.data.frame(gse12_targets)) {
  stop(
    "'gse12_targets' must be a data.frame."
  )
}

cat("GSE122151 sample metadata dimensions:\n")
print(dim(gse12_targets))
cat("\n")


# =============================================================================
# 11. RECORD THE EXISTING EXPERIMENTAL GROUP INFORMATION
# =============================================================================

cat("Existing GSE122151 sample metadata:\n")
print(gse12_targets)
cat("\n")


# =============================================================================
# 12. VERIFY THE PRESPECIFIED GENE-SET OBJECT
# =============================================================================

if (!exists("gse12_gene_sets", envir = .GlobalEnv)) {
  stop(
    "The expected GSE122151 gene-set object 'gse12_gene_sets' was not found."
  )
}

gse12_gene_sets <- get(
  "gse12_gene_sets",
  envir = .GlobalEnv
)

if (!is.list(gse12_gene_sets)) {
  stop(
    "'gse12_gene_sets' must be a list of gene sets."
  )
}

if (length(gse12_gene_sets) == 0) {
  stop(
    "'gse12_gene_sets' contains no gene sets."
  )
}

if (is.null(names(gse12_gene_sets))) {
  stop(
    "'gse12_gene_sets' must have pathway names."
  )
}

cat(
  "Number of prespecified GSE122151 gene sets detected: ",
  length(gse12_gene_sets),
  "\n\n",
  sep = ""
)

cat("Prespecified gene-set names:\n")
cat(paste(names(gse12_gene_sets), collapse = "\n"))
cat("\n\n")


# =============================================================================
# 13. VERIFY THAT GENE SETS ARE NON-EMPTY
# =============================================================================

gene_set_sizes <- vapply(
  gse12_gene_sets,
  length,
  integer(1)
)

if (any(gene_set_sizes == 0)) {
  empty_sets <- names(gene_set_sizes)[gene_set_sizes == 0]

  stop(
    "The following prespecified gene sets are empty:\n",
    paste(empty_sets, collapse = "\n")
  )
}

gene_set_inventory <- data.frame(
  pathway = names(gse12_gene_sets),
  genes_in_gmt_object = as.integer(gene_set_sizes),
  stringsAsFactors = FALSE
)

print(gene_set_inventory)

write.csv(
  gene_set_inventory,
  file = file.path(
    step_dir,
    "01_GSE122151_gene_set_inventory.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# 14. VERIFY EXPRESSION SCALE
# =============================================================================

expression_summary <- data.frame(
  minimum = min(gse12_expression),
  maximum = max(gse12_expression),
  median = median(gse12_expression),
  mean = mean(gse12_expression),
  stringsAsFactors = FALSE
)

print(expression_summary)

write.csv(
  expression_summary,
  file = file.path(
    step_dir,
    "01_GSE122151_expression_summary.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# 15. SAVE THE STEP 1 AUDIT OBJECT
# =============================================================================

step1_audit <- list(
  project_root = normalizePath(getwd(), winslash = "/", mustWork = TRUE),
  analysis = "GSE122151 experimental DDR-GSVA analysis",
  expression_object = "gse12_expression",
  expression_dimensions = dim(gse12_expression),
  sample_ids = colnames(gse12_expression),
  gene_count = nrow(gse12_expression),
  sample_count = ncol(gse12_expression),
  gene_set_object = "gse12_gene_sets",
  number_of_gene_sets = length(gse12_gene_sets),
  gene_set_names = names(gse12_gene_sets),
  gene_set_sizes = gene_set_sizes,
  metadata_object = "gse12_targets",
  timestamp = Sys.time()
)

saveRDS(
  step1_audit,
  file = file.path(
    step_dir,
    "01_GSE122151_setup_audit.rds"
  )
)


# =============================================================================
# 16. FINAL STEP 1 MESSAGE
# =============================================================================

cat(
  "\n============================================================\n",
  "STEP 1 COMPLETE\n",
  "============================================================\n",
  "GSE122151 expression object verified: YES\n",
  "Six-sample design verified: ",
  ncol(gse12_expression) == 6,
  "\n",
  "Prespecified gene-set object verified: YES\n",
  "Gene-set inventory saved: YES\n",
  "Session information saved: YES\n",
  "Setup audit saved: YES\n",
  "============================================================\n",
  sep = ""
)
