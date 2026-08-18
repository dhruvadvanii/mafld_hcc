# ============================================================
# MAFLD-HCC MANUSCRIPT PIPELINE
# STEP 6 — TCGA-LIHC SAMPLE + EXPRESSION VALIDATION
# ============================================================
#
# PURPOSE:
#   Validate and construct the TCGA-LIHC discovery expression
#   cohort before transformation, normalization and GSVA.
#
# IMPORTANT:
#   - No GSVA
#   - No pathway scoring
#   - No MAF required
#   - No differential expression
#   - No survival analysis
#   - No filtering to the historical n=163 cohort yet
#
# INPUT:
#   data/raw/TCGA-LIHC/TCGA-LIHC.star_counts.tsv
#   data/raw/TCGA-LIHC/TCGA-LIHC.clinical.tsv
#   data/raw/TCGA-LIHC/TCGA-LIHC.survival.tsv
#
# OUTPUT:
#   data/processed/TCGA-LIHC/
#       TCGA_LIHC_raw_expression.rds
#       TCGA_LIHC_sample_manifest.rds
#
#   results/
#       06_TCGA_expression_structure.csv
#       06_TCGA_sample_matching.csv
#       06_TCGA_sample_manifest.csv
#       06_TCGA_sample_summary.csv
#       06_TCGA_validation_report.txt
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
# 3. INPUT FILES
# ------------------------------------------------------------

expression_file <- file.path(
  "data",
  "raw",
  "TCGA-LIHC",
  "TCGA-LIHC.star_counts.tsv"
)

clinical_file <- file.path(
  "data",
  "raw",
  "TCGA-LIHC",
  "TCGA-LIHC.clinical.tsv"
)

survival_file <- file.path(
  "data",
  "raw",
  "TCGA-LIHC",
  "TCGA-LIHC.survival.tsv"
)


# ------------------------------------------------------------
# 4. CHECK INPUTS
# ------------------------------------------------------------

input_files <- c(
  expression = expression_file,
  clinical = clinical_file,
  survival = survival_file
)

input_status <- data.frame(

  file =
    names(input_files),

  path =
    unname(input_files),

  exists =
    file.exists(
      unname(input_files)
    ),

  size_MB =
    NA_real_,

  stringsAsFactors = FALSE
)


for (i in seq_len(
  nrow(input_status)
)) {

  if (
    input_status$exists[i]
  ) {

    input_status$size_MB[i] <-
      file.info(
        input_status$path[i]
      )$size / 1024^2

    cat(
      "[FOUND] ",
      input_status$file[i],
      " | ",
      round(
        input_status$size_MB[i],
        2
      ),
      " MB\n",
      sep = ""
    )

  } else {

    cat(
      "[MISSING] ",
      input_status$file[i],
      "\n",
      sep = ""
    )
  }
}


missing_inputs <- input_status$file[
  !input_status$exists
]


if (
  length(missing_inputs) > 0
) {

  stop(
    "Missing required TCGA input(s):\n",
    paste(
      missing_inputs,
      collapse = "\n"
    )
  )
}


cat("\n")


# ------------------------------------------------------------
# 5. REQUIRE data.table
# ------------------------------------------------------------

if (
  !requireNamespace(
    "data.table",
    quietly = TRUE
  )
) {

  stop(
    "Package 'data.table' is required.\n",
    "Install it with:\n",
    "install.packages('data.table')"
  )
}

library(data.table)


# ------------------------------------------------------------
# 6. INSPECT EXPRESSION HEADER
# ------------------------------------------------------------

cat("============================================================\n")
cat("TCGA-LIHC EXPRESSION STRUCTURE\n")
cat("============================================================\n\n")


expression_header <- names(
  fread(
    expression_file,
    nrows = 0,
    sep = "\t",
    check.names = FALSE
  )
)


cat(
  "Total expression columns: ",
  length(expression_header),
  "\n",
  sep = ""
)


cat(
  "First 10 column names:\n"
)

print(
  head(
    expression_header,
    10
  )
)

cat("\n")


# ------------------------------------------------------------
# 7. READ EXPRESSION MATRIX
# ------------------------------------------------------------

cat(
  "Reading TCGA STAR-count expression matrix...\n",
  "This may take some time.\n\n"
)


tcga_raw <- fread(
  expression_file,
  sep = "\t",
  header = TRUE,
  check.names = FALSE,
  data.table = TRUE,
  showProgress = TRUE
)


cat("\nExpression dimensions:\n")

cat(
  "Rows: ",
  nrow(tcga_raw),
  "\n",
  sep = ""
)

cat(
  "Columns: ",
  ncol(tcga_raw),
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------
# 8. IDENTIFY GENE ID COLUMN
# ------------------------------------------------------------

possible_gene_columns <- c(
  "gene_id",
  "Gene",
  "gene",
  "Ensembl",
  "Ensembl_ID",
  "ID"
)


gene_column_candidates <- intersect(
  possible_gene_columns,
  names(tcga_raw)
)


if (
  length(gene_column_candidates) == 0
) {

  stop(
    "Could not identify the gene identifier column.\n",
    "Available columns include:\n",
    paste(
      head(
        names(tcga_raw),
        20
      ),
      collapse = ", "
    )
  )
}


gene_column <- gene_column_candidates[1]


cat(
  "Gene identifier column: ",
  gene_column,
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------
# 9. EXTRACT SAMPLE COLUMNS
# ------------------------------------------------------------
#
# TCGA sample barcodes generally begin with TCGA-.
# ------------------------------------------------------------

sample_columns <- grep(
  "^TCGA-[A-Z0-9]+-[A-Z0-9]+",
  names(tcga_raw)
)


if (
  length(sample_columns) == 0
) {

  stop(
    "No TCGA sample columns were detected."
  )
}


sample_barcodes <- names(
  tcga_raw
)[sample_columns]


cat(
  "TCGA sample columns detected: ",
  length(sample_barcodes),
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------
# 10. EXTRACT GENE IDS
# ------------------------------------------------------------

gene_ids_raw <- tcga_raw[[
  gene_column
]]


gene_ids_raw <- as.character(
  gene_ids_raw
)


# ------------------------------------------------------------
# 11. REMOVE ENSEMBL VERSION SUFFIX
# ------------------------------------------------------------
#
# Example:
#
# ENSG00000141510.17
#
# becomes:
#
# ENSG00000141510
#
# This makes identifiers compatible with standard gene
# annotation resources.
# ------------------------------------------------------------

gene_ids <- sub(
  "\\..*$",
  "",
  gene_ids_raw
)


# ------------------------------------------------------------
# 12. CHECK GENE ID STRUCTURE
# ------------------------------------------------------------

ensembl_like <- grepl(
  "^ENSG[0-9]+$",
  gene_ids
)


cat(
  "Ensembl-like gene identifiers: ",
  sum(ensembl_like),
  "\n",
  sep = ""
)

cat(
  "Non-Ensembl-like identifiers: ",
  sum(!ensembl_like),
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------
# 13. REMOVE DUPLICATE GENE IDS
# ------------------------------------------------------------
#
# STAR-count files can contain duplicate/non-standard rows.
#
# We do NOT yet perform biological filtering.
#
# Duplicate Ensembl IDs will be aggregated by sum.
# ------------------------------------------------------------

expression_values <- as.matrix(
  tcga_raw[
    ,
    sample_columns,
    with = FALSE
  ]
)


storage.mode(
  expression_values
) <- "numeric"


rownames(
  expression_values
) <- gene_ids


# ------------------------------------------------------------
# 14. AGGREGATE DUPLICATE GENE IDS
# ------------------------------------------------------------

duplicate_gene_ids <- duplicated(
  gene_ids
)


cat(
  "Duplicate gene-ID rows: ",
  sum(duplicate_gene_ids),
  "\n",
  sep = ""
)


if (
  any(duplicate_gene_ids)
) {

  cat(
    "Aggregating duplicate gene IDs by sum...\n"
  )

  expression_values <- rowsum(
    expression_values,
    group = gene_ids,
    reorder = FALSE
  )

} else {

  cat(
    "No duplicate gene IDs detected.\n"
  )
}


cat(
  "Unique gene rows after aggregation: ",
  nrow(expression_values),
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------
# 15. SAMPLE BARCODE LEVEL
# ------------------------------------------------------------
#
# TCGA barcodes:
#
# TCGA-XX-YYYY-01A...
#
# The first three components identify the participant,
# while the sample type is encoded in the fourth component.
#
# We construct:
#
#   participant_id = TCGA-XX-YYYY
#   sample_type_code = first two digits of sample component
# ------------------------------------------------------------

participant_id <- sub(
  "^((TCGA-[A-Z0-9]+-[A-Z0-9]+)).*$",
  "\\1",
  sample_barcodes
)


sample_type_component <- sub(
  "^TCGA-[A-Z0-9]+-[A-Z0-9]+-([0-9]{2}).*$",
  "\\1",
  sample_barcodes
)


sample_type_code <- suppressWarnings(
  as.integer(
    sample_type_component
  )
)


# ------------------------------------------------------------
# 16. DEFINE PRIMARY TUMOR SAMPLES
# ------------------------------------------------------------
#
# TCGA sample type:
#
#   01 = Primary Solid Tumor
#
# This is the appropriate starting definition for the
# transcriptomic discovery cohort.
# ------------------------------------------------------------

is_primary_tumor <- (
  sample_type_code == 1
)


cat("============================================================\n")
cat("TCGA SAMPLE-TYPE CLASSIFICATION\n")
cat("============================================================\n\n")


sample_type_summary <- table(
  sample_type_code,
  useNA = "ifany"
)


print(
  sample_type_summary
)

cat("\n")


cat(
  "Primary tumor samples identified: ",
  sum(
    is_primary_tumor,
    na.rm = TRUE
  ),
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------
# 17. CREATE SAMPLE MANIFEST
# ------------------------------------------------------------

sample_manifest <- data.frame(

  expression_column =
    sample_columns,

  sample_barcode =
    sample_barcodes,

  participant_id =
    participant_id,

  sample_type_code =
    sample_type_code,

  is_primary_tumor =
    is_primary_tumor,

  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 18. CHECK SAMPLE BARCODE UNIQUENESS
# ------------------------------------------------------------

if (
  anyDuplicated(
    sample_manifest$sample_barcode
  ) > 0
) {

  stop(
    "Duplicate TCGA sample barcodes detected."
  )
}


# ------------------------------------------------------------
# 19. CREATE PRIMARY TUMOR EXPRESSION MATRIX
# ------------------------------------------------------------

tumor_sample_indices <- which(
  sample_manifest$is_primary_tumor
)


tumor_sample_barcodes <-
  sample_manifest$sample_barcode[
    tumor_sample_indices
  ]


tumor_expression <- expression_values[
  ,
  tumor_sample_indices,
  drop = FALSE
]


colnames(
  tumor_expression
) <- tumor_sample_barcodes


# ------------------------------------------------------------
# 20. VERIFY TUMOR MATRIX
# ------------------------------------------------------------

cat(
  "Primary tumor expression matrix:\n"
)

cat(
  "Genes: ",
  nrow(tumor_expression),
  "\n",
  sep = ""
)

cat(
  "Samples: ",
  ncol(tumor_expression),
  "\n\n",
  sep = ""
)


if (
  ncol(tumor_expression) == 0
) {

  stop(
    "No primary tumor samples were identified."
  )
}


# ------------------------------------------------------------
# 21. READ CLINICAL DATA
# ------------------------------------------------------------

cat(
  "Reading TCGA clinical data...\n"
)


clinical <- fread(
  clinical_file,
  sep = "\t",
  header = TRUE,
  check.names = FALSE,
  data.table = TRUE,
  showProgress = TRUE
)


cat(
  "Clinical dimensions: ",
  nrow(clinical),
  " x ",
  ncol(clinical),
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------
# 22. IDENTIFY CLINICAL SAMPLE/PATIENT IDENTIFIER
# ------------------------------------------------------------

clinical_names <- names(
  clinical
)


clinical_id_candidates <- c(
  "submitter_id",
  "bcr_patient_barcode",
  "case_submitter_id",
  "patient_id",
  "sample_id"
)


clinical_id_column <- intersect(
  clinical_id_candidates,
  clinical_names
)


if (
  length(clinical_id_column) == 0
) {

  stop(
    "Could not identify the patient/sample identifier ",
    "column in TCGA clinical data.\n\nAvailable columns:\n",
    paste(
      head(
        clinical_names,
        30
      ),
      collapse = ", "
    )
  )
}


clinical_id_column <- clinical_id_column[1]


cat(
  "Clinical identifier column: ",
  clinical_id_column,
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------
# 23. EXTRACT CLINICAL PATIENT IDS
# ------------------------------------------------------------

clinical_patient_id <- as.character(
  clinical[[
    clinical_id_column
  ]]
)


clinical_patient_id <- sub(
  "^((TCGA-[A-Z0-9]+-[A-Z0-9]+)).*$",
  "\\1",
  clinical_patient_id
)


# ------------------------------------------------------------
# 24. MATCH TCGA TUMOR PARTICIPANTS TO CLINICAL DATA
# ------------------------------------------------------------

tumor_participant_ids <- sub(
  "^((TCGA-[A-Z0-9]+-[A-Z0-9]+)).*$",
  "\\1",
  colnames(tumor_expression)
)


clinical_match <- match(
  tumor_participant_ids,
  clinical_patient_id
)


clinical_matched <- !is.na(
  clinical_match
)


cat(
  "Primary tumor samples: ",
  length(tumor_participant_ids),
  "\n",
  sep = ""
)

cat(
  "Tumor samples with clinical match: ",
  sum(clinical_matched),
  "\n",
  sep = ""
)

cat(
  "Tumor samples without clinical match: ",
  sum(!clinical_matched),
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------
# 25. READ SURVIVAL DATA
# ------------------------------------------------------------

cat(
  "Reading TCGA survival data...\n"
)


survival <- fread(
  survival_file,
  sep = "\t",
  header = TRUE,
  check.names = FALSE,
  data.table = TRUE,
  showProgress = TRUE
)


cat(
  "Survival dimensions: ",
  nrow(survival),
  " x ",
  ncol(survival),
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------
# 26. IDENTIFY SURVIVAL PATIENT ID
# ------------------------------------------------------------

survival_names <- names(
  survival
)


survival_id_candidates <- c(
  "sample",
  "sample_id",
  "submitter_id",
  "bcr_patient_barcode",
  "patient_id",
  "case_submitter_id"
)


survival_id_column <- intersect(
  survival_id_candidates,
  survival_names
)


if (
  length(survival_id_column) == 0
) {

  stop(
    "Could not identify patient/sample ID in survival file.\n",
    "Available columns:\n",
    paste(
      head(
        survival_names,
        30
      ),
      collapse = ", "
    )
  )
}


survival_id_column <-
  survival_id_column[1]


cat(
  "Survival identifier column: ",
  survival_id_column,
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------
# 27. MATCH SURVIVAL DATA
# ------------------------------------------------------------

survival_patient_id <- as.character(
  survival[[
    survival_id_column
  ]]
)


survival_patient_id <- sub(
  "^((TCGA-[A-Z0-9]+-[A-Z0-9]+)).*$",
  "\\1",
  survival_patient_id
)


survival_match <- match(
  tumor_participant_ids,
  survival_patient_id
)


survival_matched <- !is.na(
  survival_match
)


cat(
  "Tumor samples with survival match: ",
  sum(survival_matched),
  "\n",
  sep = ""
)

cat(
  "Tumor samples without survival match: ",
  sum(!survival_matched),
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------
# 28. COMBINED MATCHING STATUS
# ------------------------------------------------------------

sample_manifest_tumor <- sample_manifest[
  sample_manifest$is_primary_tumor,
  ,
  drop = FALSE
]


sample_manifest_tumor$clinical_match <-
  clinical_matched


sample_manifest_tumor$survival_match <-
  survival_matched


sample_manifest_tumor$clinical_row <-
  clinical_match


sample_manifest_tumor$survival_row <-
  survival_match


# ------------------------------------------------------------
# 29. SAVE SAMPLE MANIFEST
# ------------------------------------------------------------

write.csv(
  sample_manifest,
  file = file.path(
    results_dir,
    "06_TCGA_sample_manifest.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 30. MATCHING SUMMARY
# ------------------------------------------------------------

matching_summary <- data.frame(

  metric = c(
    "total_expression_sample_columns",
    "primary_tumor_samples",
    "clinical_matches",
    "clinical_unmatched",
    "survival_matches",
    "survival_unmatched",
    "clinical_and_survival_matches"
  ),

  value = c(

    length(sample_barcodes),

    ncol(tumor_expression),

    sum(clinical_matched),

    sum(!clinical_matched),

    sum(survival_matched),

    sum(!survival_matched),

    sum(
      clinical_matched &
      survival_matched
    )
  ),

  stringsAsFactors = FALSE
)


write.csv(
  matching_summary,
  file = file.path(
    results_dir,
    "06_TCGA_sample_matching.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 31. SAMPLE TYPE SUMMARY
# ------------------------------------------------------------

sample_summary <- data.frame(

  sample_type_code =
    names(
      sample_type_summary
    ),

  n =
    as.integer(
      sample_type_summary
    ),

  stringsAsFactors = FALSE
)


write.csv(
  sample_summary,
  file = file.path(
    results_dir,
    "06_TCGA_sample_summary.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 32. EXPRESSION STRUCTURE SUMMARY
# ------------------------------------------------------------

expression_structure <- data.frame(

  metric = c(
    "expression_rows_initial",
    "expression_columns_initial",
    "gene_identifier_column",
    "TCGA_sample_columns",
    "unique_gene_rows_after_aggregation",
    "primary_tumor_samples"
  ),

  value = c(

    nrow(tcga_raw),

    ncol(tcga_raw),

    gene_column,

    length(sample_barcodes),

    nrow(expression_values),

    ncol(tumor_expression)
  ),

  stringsAsFactors = FALSE
)


write.csv(
  expression_structure,
  file = file.path(
    results_dir,
    "06_TCGA_expression_structure.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 33. SAVE RAW PRIMARY-TUMOR EXPRESSION OBJECT
# ------------------------------------------------------------
#
# IMPORTANT:
#   This object is still COUNT DATA.
#
# No normalization or transformation has yet been performed.
# ------------------------------------------------------------

tcga_expression_object <- list(

  created_at =
    Sys.time(),

  source =
    "TCGA-LIHC STAR counts",

  gene_identifier =
    gene_column,

  n_genes =
    nrow(tumor_expression),

  n_samples =
    ncol(tumor_expression),

  sample_names =
    colnames(tumor_expression),

  expression_counts =
    tumor_expression
)


saveRDS(
  tcga_expression_object,
  file = file.path(
    processed_dir,
    "TCGA_LIHC_raw_expression.rds"
  )
)


# ------------------------------------------------------------
# 34. SAVE SAMPLE MANIFEST OBJECT
# ------------------------------------------------------------

tcga_manifest_object <- list(

  created_at =
    Sys.time(),

  total_expression_samples =
    length(sample_barcodes),

  primary_tumor_samples =
    ncol(tumor_expression),

  sample_manifest =
    sample_manifest,

  tumor_manifest =
    sample_manifest_tumor
)


saveRDS(
  tcga_manifest_object,
  file = file.path(
    processed_dir,
    "TCGA_LIHC_sample_manifest.rds"
  )
)


# ------------------------------------------------------------
# 35. WRITE VALIDATION REPORT
# ------------------------------------------------------------

report_lines <- c(

  "============================================================",
  "TCGA-LIHC SAMPLE + EXPRESSION VALIDATION REPORT",
  "============================================================",
  "",

  paste(
    "Date:",
    Sys.time()
  ),

  "",

  "EXPRESSION",
  "------------------------------------------------------------",

  paste(
    "Initial expression rows:",
    nrow(tcga_raw)
  ),

  paste(
    "Initial expression columns:",
    ncol(tcga_raw)
  ),

  paste(
    "Gene identifier column:",
    gene_column
  ),

  paste(
    "TCGA sample columns:",
    length(sample_barcodes)
  ),

  paste(
    "Unique gene rows:",
    nrow(expression_values)
  ),

  "",

  "SAMPLE TYPE",
  "------------------------------------------------------------",

  paste(
    "Primary tumor samples:",
    ncol(tumor_expression)
  ),

  "",

  "CLINICAL MATCHING",
  "------------------------------------------------------------",

  paste(
    "Clinical matches:",
    sum(clinical_matched)
  ),

  paste(
    "Clinical unmatched:",
    sum(!clinical_matched)
  ),

  "",

  "SURVIVAL MATCHING",
  "------------------------------------------------------------",

  paste(
    "Survival matches:",
    sum(survival_matched)
  ),

  paste(
    "Survival unmatched:",
    sum(!survival_matched)
  ),

  "",

  "CURRENT PROCESSING STATUS",
  "------------------------------------------------------------",

  "Expression transformation: NOT PERFORMED",

  "Normalization: NOT PERFORMED",

  "Gene filtering: NOT PERFORMED",

  "GSVA: NOT PERFORMED",

  "MAF: NOT REQUIRED",

  "Historical n=163 restriction: NOT APPLIED",

  "",

  "NEXT STEP",
  "------------------------------------------------------------",

  "Determine the appropriate RNA-seq transformation and",
  "gene-level filtering strategy before GSVA."
)


writeLines(
  report_lines,
  con = file.path(
    results_dir,
    "06_TCGA_validation_report.txt"
  )
)


# ------------------------------------------------------------
# 36. CLEAN MEMORY
# ------------------------------------------------------------

rm(
  tcga_raw,
  expression_values
)

gc()


# ------------------------------------------------------------
# 37. FINAL STATUS
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("STEP 6 COMPLETED SUCCESSFULLY\n")
cat("============================================================\n\n")

cat(
  "TCGA expression samples detected: ",
  length(sample_barcodes),
  "\n",
  sep = ""
)

cat(
  "Primary tumor samples: ",
  ncol(tumor_expression),
  "\n",
  sep = ""
)

cat(
  "Unique gene rows: ",
  nrow(tumor_expression),
  "\n",
  sep = ""
)

cat(
  "Clinical matches: ",
  sum(clinical_matched),
  "\n",
  sep = ""
)

cat(
  "Survival matches: ",
  sum(survival_matched),
  "\n\n",
  sep = ""
)

cat(
  "Expression transformation: NOT PERFORMED\n"
)

cat(
  "Normalization: NOT PERFORMED\n"
)

cat(
  "GSVA: NOT PERFORMED\n"
)

cat(
  "MAF: NOT REQUIRED\n\n"
)

cat(
  "Primary tumor expression object saved to:\n",
  file.path(
    processed_dir,
    "TCGA_LIHC_raw_expression.rds"
  ),
  "\n\n",
  sep = ""
)

cat(
  "Next step: TCGA RNA-seq preprocessing.\n"
)
