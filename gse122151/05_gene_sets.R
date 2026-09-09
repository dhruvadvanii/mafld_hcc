# ============================================================
# GSE122151 — STEP 5: GENE-SET PREPARATION
# ============================================================
#
# Purpose:
#   1. Retrieve Hallmark gene sets
#   2. Retrieve Reactome gene sets
#   3. Standardize gene symbols
#   4. Assess representation in GSE122151
#   5. Define usable gene sets for GSVA
#   6. Save all gene sets and coverage information
#
# Input:
#   data/processed/GSE122151/
#       GSE122151_log2_FPKM_plus1.rds
#
# Output:
#   data/processed/GSE122151/
#       Hallmark_gene_sets.rds
#       Reactome_gene_sets.rds
#       usable_Hallmark_gene_sets.rds
#       usable_Reactome_gene_sets.rds
#
#   results/GSE122151/
#       gene_set_coverage.csv
#       Hallmark_gene_set_summary.csv
#       Reactome_gene_set_summary.csv
#       gene_set_source_metadata.csv
#       GSE122151_gene_set_summary.xlsx
#
# ============================================================


# ------------------------------------------------------------
# STEP 0 — Clean environment
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
# STEP 2 — Check required packages
# ------------------------------------------------------------

required_packages <- c(
  "msigdbr",
  "dplyr",
  "tidyr",
  "tibble",
  "openxlsx"
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
      collapse = "\n"
    ),
    "\n\nRun Step 1 before running Step 5."
  )
}


# ------------------------------------------------------------
# STEP 3 — Load packages
# ------------------------------------------------------------

suppressPackageStartupMessages({

  library(msigdbr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(openxlsx)

})


# ------------------------------------------------------------
# STEP 4 — Check expression matrix
# ------------------------------------------------------------

expression_file <- file.path(
  processed_dir,
  "GSE122151_log2_FPKM_plus1.rds"
)

if (!file.exists(expression_file)) {

  stop(
    "Could not find:\n",
    expression_file,
    "\n\nRun Step 3 first."
  )
}

expr_log2 <- readRDS(
  expression_file
)

expr_log2 <- as.matrix(
  expr_log2
)

storage.mode(expr_log2) <- "numeric"


# ------------------------------------------------------------
# STEP 5 — Validate expression matrix
# ------------------------------------------------------------

if (is.null(rownames(expr_log2))) {

  stop(
    "Expression matrix does not contain gene identifiers ",
    "as row names."
  )
}

if (anyDuplicated(rownames(expr_log2)) > 0) {

  stop(
    "Expression matrix contains duplicated gene identifiers."
  )
}

if (any(!is.finite(expr_log2))) {

  stop(
    "Expression matrix contains non-finite values."
  )
}

expression_genes <- unique(
  toupper(
    trimws(
      rownames(expr_log2)
    )
  )
)

expression_genes <- expression_genes[
  !is.na(expression_genes) &
    expression_genes != ""
]


cat(
  "\nNumber of genes in GSE122151 expression matrix:",
  length(expression_genes),
  "\n\n"
)


# ============================================================
# STEP 6 — msigdbr helper functions
# ============================================================

# ------------------------------------------------------------
# The msigdbr API has changed across package versions.
#
# This helper attempts the current API first and then falls
# back to older parameter names where necessary.
# ------------------------------------------------------------

get_hallmark_msigdbr <- function() {

  # Current API
  result <- tryCatch(

    msigdbr::msigdbr(
      species = "Homo sapiens",
      collection = "H"
    ),

    error = function(e) {
      NULL
    }

  )

  # Older API fallback
  if (is.null(result)) {

    result <- tryCatch(

      msigdbr::msigdbr(
        species = "Homo sapiens",
        category = "H"
      ),

      error = function(e) {
        NULL
      }

    )
  }

  if (is.null(result)) {

    stop(
      "Could not retrieve the Hallmark collection using ",
      "the installed msigdbr version."
    )
  }

  result
}


# ------------------------------------------------------------
# Reactome retrieval
# ------------------------------------------------------------

get_reactome_msigdbr <- function() {

  # Current API
  result <- tryCatch(

    msigdbr::msigdbr(
      species = "Homo sapiens",
      collection = "C2",
      subcollection = "CP:REACTOME"
    ),

    error = function(e) {
      NULL
    }

  )

  # Older API fallback
  if (is.null(result)) {

    result <- tryCatch(

      msigdbr::msigdbr(
        species = "Homo sapiens",
        category = "C2",
        subcategory = "CP:REACTOME"
      ),

      error = function(e) {
        NULL
      }

    )
  }

  if (is.null(result)) {

    stop(
      "Could not retrieve the Reactome collection using ",
      "the installed msigdbr version."
    )
  }

  result
}


# ============================================================
# STEP 7 — Retrieve Hallmark
# ============================================================

cat(
  "Retrieving Hallmark gene sets...\n"
)

hallmark_raw <- get_hallmark_msigdbr()

cat(
  "Hallmark rows retrieved:",
  nrow(hallmark_raw),
  "\n"
)


# ============================================================
# STEP 8 — Retrieve Reactome
# ============================================================

cat(
  "Retrieving Reactome gene sets...\n"
)

reactome_raw <- get_reactome_msigdbr()

cat(
  "Reactome rows retrieved:",
  nrow(reactome_raw),
  "\n\n"
)


# ============================================================
# STEP 9 — Identify required msigdbr columns
# ============================================================

# ------------------------------------------------------------
# Gene-symbol column
# ------------------------------------------------------------

possible_gene_columns <- c(
  "gene_symbol",
  "gene_symbol.x"
)

hallmark_gene_column <- possible_gene_columns[
  possible_gene_columns %in%
    colnames(hallmark_raw)
][1]

reactome_gene_column <- possible_gene_columns[
  possible_gene_columns %in%
    colnames(reactome_raw)
][1]


if (is.na(hallmark_gene_column)) {

  stop(
    "Could not identify the gene-symbol column in the ",
    "Hallmark msigdbr output."
  )
}

if (is.na(reactome_gene_column)) {

  stop(
    "Could not identify the gene-symbol column in the ",
    "Reactome msigdbr output."
  )
}


# ------------------------------------------------------------
# Gene-set name column
# ------------------------------------------------------------

possible_set_columns <- c(
  "gs_name",
  "gs_name.x"
)

hallmark_set_column <- possible_set_columns[
  possible_set_columns %in%
    colnames(hallmark_raw)
][1]

reactome_set_column <- possible_set_columns[
  possible_set_columns %in%
    colnames(reactome_raw)
][1]


if (is.na(hallmark_set_column)) {

  stop(
    "Could not identify the gene-set name column in the ",
    "Hallmark msigdbr output."
  )
}

if (is.na(reactome_set_column)) {

  stop(
    "Could not identify the gene-set name column in the ",
    "Reactome msigdbr output."
  )
}


# ============================================================
# STEP 10 — Standardize Hallmark table
# ============================================================

hallmark_tbl <- hallmark_raw %>%

  dplyr::transmute(

    gs_name =
      .data[[hallmark_set_column]],

    gene_symbol =
      toupper(
        trimws(
          .data[[hallmark_gene_column]]
        )
      )

  ) %>%

  filter(
    !is.na(gs_name),
    gs_name != "",
    !is.na(gene_symbol),
    gene_symbol != ""
  ) %>%

  distinct()


# ============================================================
# STEP 11 — Standardize Reactome table
# ============================================================

reactome_tbl <- reactome_raw %>%

  dplyr::transmute(

    gs_name =
      .data[[reactome_set_column]],

    gene_symbol =
      toupper(
        trimws(
          .data[[reactome_gene_column]]
        )
      )

  ) %>%

  filter(
    !is.na(gs_name),
    gs_name != "",
    !is.na(gene_symbol),
    gene_symbol != ""
  ) %>%

  distinct()


# ============================================================
# STEP 12 — Convert to gene-set lists
# ============================================================

hallmark_sets <- split(
  hallmark_tbl$gene_symbol,
  hallmark_tbl$gs_name
)

reactome_sets <- split(
  reactome_tbl$gene_symbol,
  reactome_tbl$gs_name
)


# ------------------------------------------------------------
# Remove duplicate genes within individual gene sets
# ------------------------------------------------------------

hallmark_sets <- lapply(
  hallmark_sets,
  unique
)

reactome_sets <- lapply(
  reactome_sets,
  unique
)


# ------------------------------------------------------------
# Remove empty sets
# ------------------------------------------------------------

hallmark_sets <- hallmark_sets[
  lengths(hallmark_sets) > 0
]

reactome_sets <- reactome_sets[
  lengths(reactome_sets) > 0
]


# ============================================================
# STEP 13 — Basic gene-set validation
# ============================================================

if (length(hallmark_sets) == 0) {

  stop(
    "No Hallmark gene sets were successfully generated."
  )
}

if (length(reactome_sets) == 0) {

  stop(
    "No Reactome gene sets were successfully generated."
  )
}


cat(
  "\nHallmark gene sets:",
  length(hallmark_sets),
  "\n"
)

cat(
  "Reactome gene sets:",
  length(reactome_sets),
  "\n\n"
)


# ============================================================
# STEP 14 — Gene-set coverage
# ============================================================

calculate_coverage <- function(
    gene_sets,
    expression_genes,
    collection
) {

  result <- data.frame(

    collection =
      collection,

    pathway =
      names(gene_sets),

    total_genes =
      lengths(gene_sets),

    genes_present =
      sapply(
        gene_sets,
        function(x) {
          sum(
            x %in% expression_genes
          )
        }
      ),

    stringsAsFactors = FALSE
  )

  result$coverage_percent <-
    100 *
    result$genes_present /
    result$total_genes

  result
}


hallmark_coverage <- calculate_coverage(
  gene_sets = hallmark_sets,
  expression_genes = expression_genes,
  collection = "Hallmark"
)


reactome_coverage <- calculate_coverage(
  gene_sets = reactome_sets,
  expression_genes = expression_genes,
  collection = "Reactome"
)


gene_set_coverage <- dplyr::bind_rows(
  hallmark_coverage,
  reactome_coverage
)


# ============================================================
# STEP 15 — Define usable gene sets for GSVA
# ============================================================
#
# We require at least five genes from each gene set to be
# represented in the GSE122151 expression matrix.
#
# This does NOT alter the biological gene sets themselves.
# It only prevents poorly represented sets from being used
# for GSVA.
#
# ============================================================

minimum_genes_present <- 5

hallmark_usable_names <-
  hallmark_coverage$pathway[
    hallmark_coverage$genes_present >=
      minimum_genes_present
  ]


reactome_usable_names <-
  reactome_coverage$pathway[
    reactome_coverage$genes_present >=
      minimum_genes_present
  ]


usable_hallmark <- hallmark_sets[
  hallmark_usable_names
]

usable_reactome <- reactome_sets[
  reactome_usable_names
]


# ------------------------------------------------------------
# Validate usable sets
# ------------------------------------------------------------

if (length(usable_hallmark) == 0) {

  stop(
    "No Hallmark pathways have at least ",
    minimum_genes_present,
    " genes represented in GSE122151."
  )
}

if (length(usable_reactome) == 0) {

  stop(
    "No Reactome pathways have at least ",
    minimum_genes_present,
    " genes represented in GSE122151."
  )
}


# ============================================================
# STEP 16 — Add usability flag to coverage table
# ============================================================

gene_set_coverage <- gene_set_coverage %>%

  mutate(

    minimum_genes_for_GSVA =
      minimum_genes_present,

    usable_for_GSVA =
      genes_present >=
        minimum_genes_present

  )


# ============================================================
# STEP 17 — Gene-set summary tables
# ============================================================

hallmark_summary <- hallmark_coverage %>%

  mutate(

    usable_for_GSVA =
      genes_present >=
        minimum_genes_present

  ) %>%

  arrange(
    desc(genes_present),
    pathway
  )


reactome_summary <- reactome_coverage %>%

  mutate(

    usable_for_GSVA =
      genes_present >=
        minimum_genes_present

  ) %>%

  arrange(
    desc(genes_present),
    pathway
  )


# ============================================================
# STEP 18 — Save complete gene-set objects
# ============================================================

saveRDS(
  hallmark_sets,
  file = file.path(
    processed_dir,
    "Hallmark_gene_sets.rds"
  )
)

saveRDS(
  reactome_sets,
  file = file.path(
    processed_dir,
    "Reactome_gene_sets.rds"
  )
)

saveRDS(
  usable_hallmark,
  file = file.path(
    processed_dir,
    "usable_Hallmark_gene_sets.rds"
  )
)

saveRDS(
  usable_reactome,
  file = file.path(
    processed_dir,
    "usable_Reactome_gene_sets.rds"
  )
)


# ============================================================
# STEP 19 — Save gene sets as CSV
# ============================================================

hallmark_long <- hallmark_tbl %>%

  arrange(
    gs_name,
    gene_symbol
  )


reactome_long <- reactome_tbl %>%

  arrange(
    gs_name,
    gene_symbol
  )


write.csv(
  hallmark_long,
  file = file.path(
    results_dir,
    "Hallmark_gene_sets_long.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)


write.csv(
  reactome_long,
  file = file.path(
    results_dir,
    "Reactome_gene_sets_long.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)


# ============================================================
# STEP 20 — Save coverage tables
# ============================================================

write.csv(
  gene_set_coverage,
  file = file.path(
    results_dir,
    "gene_set_coverage.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)


write.csv(
  hallmark_summary,
  file = file.path(
    results_dir,
    "Hallmark_gene_set_summary.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)


write.csv(
  reactome_summary,
  file = file.path(
    results_dir,
    "Reactome_gene_set_summary.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)


# ============================================================
# STEP 21 — Save source/version information
# ============================================================

msigdbr_version <- as.character(
  packageVersion("msigdbr")
)

gene_set_source_metadata <- data.frame(

  source =
    c(
      "MSigDB Hallmark",
      "MSigDB Reactome"
    ),

  collection =
    c(
      "H",
      "C2"
    ),

  subcollection =
    c(
      NA,
      "CP:REACTOME"
    ),

  species =
    c(
      "Homo sapiens",
      "Homo sapiens"
    ),

  msigdbr_version =
    c(
      msigdbr_version,
      msigdbr_version
    ),

  total_gene_sets =
    c(
      length(hallmark_sets),
      length(reactome_sets)
    ),

  usable_gene_sets =
    c(
      length(usable_hallmark),
      length(usable_reactome)
    ),

  minimum_genes_present_for_GSVA =
    c(
      minimum_genes_present,
      minimum_genes_present
    ),

  stringsAsFactors = FALSE
)


write.csv(
  gene_set_source_metadata,
  file = file.path(
    results_dir,
    "gene_set_source_metadata.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)


# ============================================================
# STEP 22 — Master Excel export
# ============================================================

xlsx_file <- file.path(
  results_dir,
  "GSE122151_gene_set_summary.xlsx"
)

wb <- openxlsx::createWorkbook()


# ------------------------------------------------------------
# Hallmark summary
# ------------------------------------------------------------

openxlsx::addWorksheet(
  wb,
  "Hallmark_summary"
)

openxlsx::writeData(
  wb,
  "Hallmark_summary",
  hallmark_summary
)

openxlsx::freezePane(
  wb,
  "Hallmark_summary",
  firstRow = TRUE
)


# ------------------------------------------------------------
# Reactome summary
# ------------------------------------------------------------

openxlsx::addWorksheet(
  wb,
  "Reactome_summary"
)

openxlsx::writeData(
  wb,
  "Reactome_summary",
  reactome_summary
)

openxlsx::freezePane(
  wb,
  "Reactome_summary",
  firstRow = TRUE
)


# ------------------------------------------------------------
# Combined coverage
# ------------------------------------------------------------

openxlsx::addWorksheet(
  wb,
  "Gene_set_coverage"
)

openxlsx::writeData(
  wb,
  "Gene_set_coverage",
  gene_set_coverage
)

openxlsx::freezePane(
  wb,
  "Gene_set_coverage",
  firstRow = TRUE
)


# ------------------------------------------------------------
# Source metadata
# ------------------------------------------------------------

openxlsx::addWorksheet(
  wb,
  "Source_metadata"
)

openxlsx::writeData(
  wb,
  "Source_metadata",
  gene_set_source_metadata
)


# ------------------------------------------------------------
# Save workbook
# ------------------------------------------------------------

openxlsx::saveWorkbook(
  wb,
  xlsx_file,
  overwrite = TRUE
)


# ============================================================
# STEP 23 — Save an analysis object
# ============================================================

gene_set_object <- list(

  hallmark_sets =
    hallmark_sets,

  reactome_sets =
    reactome_sets,

  usable_hallmark =
    usable_hallmark,

  usable_reactome =
    usable_reactome,

  hallmark_coverage =
    hallmark_coverage,

  reactome_coverage =
    reactome_coverage,

  combined_coverage =
    gene_set_coverage,

  expression_genes =
    expression_genes,

  minimum_genes_present =
    minimum_genes_present,

  msigdbr_version =
    msigdbr_version

)


saveRDS(
  gene_set_object,
  file = file.path(
    processed_dir,
    "GSE122151_gene_set_object.rds"
  )
)


# ============================================================
# STEP 24 — Final validation
# ============================================================

# Confirm Hallmark sets contain genes
if (
  length(hallmark_sets) == 0
) {

  stop(
    "Final validation failed: no Hallmark gene sets."
  )
}


# Confirm Reactome sets contain genes
if (
  length(reactome_sets) == 0
) {

  stop(
    "Final validation failed: no Reactome gene sets."
  )
}


# Confirm usable Hallmark sets
if (
  length(usable_hallmark) == 0
) {

  stop(
    "Final validation failed: no usable Hallmark gene sets."
  )
}


# Confirm usable Reactome sets
if (
  length(usable_reactome) == 0
) {

  stop(
    "Final validation failed: no usable Reactome gene sets."
  )
}


# Confirm all usable sets meet threshold
hallmark_minimum_check <- min(
  sapply(
    usable_hallmark,
    function(x) {
      sum(
        x %in% expression_genes
      )
    }
  )
)


reactome_minimum_check <- min(
  sapply(
    usable_reactome,
    function(x) {
      sum(
        x %in% expression_genes
      )
    }
  )
)


if (
  hallmark_minimum_check <
    minimum_genes_present
) {

  stop(
    "Final validation failed: a Hallmark gene set has ",
    "fewer than the required number of represented genes."
  )
}


if (
  reactome_minimum_check <
    minimum_genes_present
) {

  stop(
    "Final validation failed: a Reactome gene set has ",
    "fewer than the required number of represented genes."
  )
}


# ============================================================
# STEP 25 — Console report
# ============================================================

cat(
  "\n============================================================\n"
)

cat(
  "GSE122151 STEP 5 — GENE-SET PREPARATION COMPLETE\n"
)

cat(
  "============================================================\n\n"
)

cat(
  "msigdbr version:",
  msigdbr_version,
  "\n\n"
)

cat(
  "Hallmark gene sets retrieved:",
  length(hallmark_sets),
  "\n"
)

cat(
  "Hallmark gene sets usable for GSVA:",
  length(usable_hallmark),
  "\n\n"
)

cat(
  "Reactome gene sets retrieved:",
  length(reactome_sets),
  "\n"
)

cat(
  "Reactome gene sets usable for GSVA:",
  length(usable_reactome),
  "\n\n"
)

cat(
  "Minimum represented genes required:",
  minimum_genes_present,
  "\n\n"
)

cat(
  "Saved to:\n",
  processed_dir,
  "\n\n"
)

cat(
  "Results saved to:\n",
  results_dir,
  "\n\n"
)

cat(
  "Master Excel:\n",
  xlsx_file,
  "\n\n"
)

cat(
  "============================================================\n"
)
