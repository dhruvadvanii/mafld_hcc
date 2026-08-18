# ============================================================
# MAFLD-HCC MANUSCRIPT PIPELINE
# STEP 5 — GSE76427 SAMPLE / PHENOTYPE MANIFEST
# ============================================================
#
# PURPOSE:
#   Define the validation cohort correctly before any pathway
#   scoring or correlation analysis.
#
# GSE76427:
#   167 total samples
#   115 primary HCC tumor samples (PT)
#    52 adjacent non-tumor samples (ANTT)
#
# VALIDATION COHORT FOR THIS MANUSCRIPT:
#   Primary HCC tumor samples only
#   Expected n = 115
#
# INPUT:
#   data/processed/GSE76427/GSE76427_gene_expression.rds
#
# OUTPUT:
#   data/processed/GSE76427/
#       GSE76427_tumor_expression.rds
#       GSE76427_sample_manifest.rds
#
#   results/
#       05_GSE76427_sample_manifest.csv
#       05_GSE76427_cohort_summary.csv
#       05_GSE76427_sample_manifest_report.txt
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
# 2. DEFINE DIRECTORIES
# ------------------------------------------------------------

processed_dir <- file.path(
  "data",
  "processed",
  "GSE76427"
)

results_dir <- "results"


# ------------------------------------------------------------
# 3. CHECK INPUT
# ------------------------------------------------------------

input_rds <- file.path(
  processed_dir,
  "GSE76427_gene_expression.rds"
)

if (!file.exists(input_rds)) {

  stop(
    "Final GSE76427 gene-expression matrix not found:\n",
    normalizePath(
      input_rds,
      mustWork = FALSE
    ),
    "\nRun Step 4 first."
  )
}


# ------------------------------------------------------------
# 4. LOAD FINAL GENE EXPRESSION MATRIX
# ------------------------------------------------------------

cat(
  "Loading final GSE76427 gene-expression matrix...\n"
)

gse_object <- readRDS(
  input_rds
)

gene_expression <- gse_object$gene_expression


cat("\nInput matrix:\n")

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


# ------------------------------------------------------------
# 5. EXTRACT SAMPLE NAMES
# ------------------------------------------------------------

sample_names <- colnames(
  gene_expression
)


if (is.null(sample_names)) {

  stop(
    "The gene-expression matrix has no sample names."
  )
}


# ------------------------------------------------------------
# 6. CLASSIFY SAMPLE TYPE
# ------------------------------------------------------------
#
# GSE76427 sample names use:
#
#   PT##   = Primary tumor
#   ANTT## = Adjacent non-tumor tissue
#
# This classification is based on the official GEO sample
# naming scheme for GSE76427.
# ------------------------------------------------------------

sample_type <- ifelse(

  grepl(
    "^PT[0-9]+$",
    sample_names,
    ignore.case = FALSE
  ),

  "Primary HCC tumor",

  ifelse(

    grepl(
      "^ANTT[0-9]+$",
      sample_names,
      ignore.case = FALSE
    ),

    "Adjacent non-tumor",

    NA_character_
  )
)


# ------------------------------------------------------------
# 7. EXTRACT PATIENT NUMBER
# ------------------------------------------------------------

patient_number <- rep(
  NA_integer_,
  length(sample_names)
)


pt_index <- grepl(
  "^PT[0-9]+$",
  sample_names
)

antt_index <- grepl(
  "^ANTT[0-9]+$",
  sample_names
)


patient_number[pt_index] <-
  as.integer(
    sub(
      "^PT",
      "",
      sample_names[pt_index]
    )
  )


patient_number[antt_index] <-
  as.integer(
    sub(
      "^ANTT",
      "",
      sample_names[antt_index]
    )
  )


# ------------------------------------------------------------
# 8. CREATE SAMPLE MANIFEST
# ------------------------------------------------------------

sample_manifest <- data.frame(

  expression_column =
    seq_along(sample_names),

  sample_name =
    sample_names,

  patient_number =
    patient_number,

  sample_type =
    sample_type,

  is_primary_tumor =
    sample_type ==
    "Primary HCC tumor",

  is_adjacent_non_tumor =
    sample_type ==
    "Adjacent non-tumor",

  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 9. CHECK FOR UNCLASSIFIED SAMPLES
# ------------------------------------------------------------

unclassified <- sample_manifest[
  is.na(
    sample_manifest$sample_type
  ),
  ,
  drop = FALSE
]


if (nrow(unclassified) > 0) {

  cat(
    "UNCLASSIFIED SAMPLES:\n"
  )

  print(
    unclassified
  )

  stop(
    "\nOne or more samples could not be classified."
  )
}


# ------------------------------------------------------------
# 10. COUNT COHORTS
# ------------------------------------------------------------

tumor_n <- sum(
  sample_manifest$is_primary_tumor
)

antt_n <- sum(
  sample_manifest$is_adjacent_non_tumor
)

total_n <- nrow(
  sample_manifest
)


cat("============================================================\n")
cat("GSE76427 COHORT COMPOSITION\n")
cat("============================================================\n\n")

cat(
  "Total samples: ",
  total_n,
  "\n",
  sep = ""
)

cat(
  "Primary HCC tumor: ",
  tumor_n,
  "\n",
  sep = ""
)

cat(
  "Adjacent non-tumor: ",
  antt_n,
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------
# 11. VERIFY EXPECTED COHORT SIZE
# ------------------------------------------------------------

if (total_n != 167) {

  stop(
    "Expected 167 total GSE76427 samples, found ",
    total_n,
    "."
  )
}


if (tumor_n != 115) {

  stop(
    "Expected 115 primary HCC tumor samples, found ",
    tumor_n,
    "."
  )
}


if (antt_n != 52) {

  stop(
    "Expected 52 adjacent non-tumor samples, found ",
    antt_n,
    "."
  )
}


if (
  tumor_n +
  antt_n !=
  total_n
) {

  stop(
    "Cohort counts do not sum to total sample count."
  )
}


# ------------------------------------------------------------
# 12. CHECK SAMPLE NAME UNIQUENESS
# ------------------------------------------------------------

if (
  anyDuplicated(
    sample_manifest$sample_name
  ) > 0
) {

  stop(
    "Duplicate sample names detected."
  )
}


# ------------------------------------------------------------
# 13. CHECK EXPRESSION/MANIFEST ORDER
# ------------------------------------------------------------

if (!identical(
  colnames(gene_expression),
  sample_manifest$sample_name
)) {

  stop(
    "Sample manifest order does not match expression matrix."
  )
}


cat(
  "Expression/manifest sample order: VERIFIED\n\n"
)


# ------------------------------------------------------------
# 14. CHECK PATIENT NUMBER DUPLICATION
# ------------------------------------------------------------
#
# PT and ANTT samples from the same patient should have the
# same patient number.
#
# This allows us to identify the expected paired structure.
# ------------------------------------------------------------

patient_sample_counts <- table(
  sample_manifest$patient_number
)


paired_patient_numbers <- names(
  patient_sample_counts[
    patient_sample_counts == 2
  ]
)


single_sample_patient_numbers <- names(
  patient_sample_counts[
    patient_sample_counts == 1
  ]
)


cat(
  "Patients represented by both PT and ANTT samples: ",
  length(paired_patient_numbers),
  "\n",
  sep = ""
)

cat(
  "Patient numbers represented by one sample only: ",
  length(single_sample_patient_numbers),
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------
# 15. CHECK THAT TWO-SAMPLE PATIENTS ARE ACTUALLY PT + ANTT
# ------------------------------------------------------------

paired_patient_check <- sample_manifest[
  sample_manifest$patient_number %in%
    as.integer(paired_patient_numbers),
  ,
  drop = FALSE
]


patient_pair_valid <- TRUE


if (nrow(paired_patient_check) > 0) {

  pair_types <- split(
    paired_patient_check$sample_type,
    paired_patient_check$patient_number
  )

  for (p in names(pair_types)) {

    types <- sort(
      unique(
        pair_types[[p]]
      )
    )

    if (!identical(
      types,
      sort(
        c(
          "Adjacent non-tumor",
          "Primary HCC tumor"
        )
      )
    )) {

      patient_pair_valid <- FALSE

      cat(
        "Unexpected patient pairing for patient ",
        p,
        ":\n",
        sep = ""
      )

      print(types)
    }
  }
}


if (!patient_pair_valid) {

  stop(
    "One or more paired patients do not have exactly ",
    "one PT and one ANTT sample."
  )
}


# ------------------------------------------------------------
# 16. SAVE FULL SAMPLE MANIFEST
# ------------------------------------------------------------

manifest_csv <- file.path(
  results_dir,
  "05_GSE76427_sample_manifest.csv"
)

write.csv(
  sample_manifest,
  file = manifest_csv,
  row.names = FALSE
)


# ------------------------------------------------------------
# 17. CREATE COHORT SUMMARY
# ------------------------------------------------------------

cohort_summary <- data.frame(

  cohort = c(
    "All GSE76427 samples",
    "Primary HCC tumor",
    "Adjacent non-tumor"
  ),

  n = c(
    total_n,
    tumor_n,
    antt_n
  ),

  percent = c(
    100,
    100 * tumor_n / total_n,
    100 * antt_n / total_n
  ),

  stringsAsFactors = FALSE
)


write.csv(
  cohort_summary,
  file = file.path(
    results_dir,
    "05_GSE76427_cohort_summary.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 18. CREATE TUMOR-ONLY EXPRESSION MATRIX
# ------------------------------------------------------------
#
# THIS is the expression matrix that will be used for the
# external validation analyses.
#
# n = 115 primary HCC tumors.
# ------------------------------------------------------------

tumor_samples <- sample_manifest[
  sample_manifest$is_primary_tumor,
  "sample_name"
]


tumor_expression <- gene_expression[
  ,
  tumor_samples,
  drop = FALSE
]


# ------------------------------------------------------------
# 19. VERIFY TUMOR MATRIX
# ------------------------------------------------------------

if (
  ncol(tumor_expression) !=
  115
) {

  stop(
    "Tumor-only matrix does not contain exactly 115 samples."
  )
}


if (!identical(
  colnames(tumor_expression),
  tumor_samples
)) {

  stop(
    "Tumor-only matrix sample order is incorrect."
  )
}


# ------------------------------------------------------------
# 20. SAVE TUMOR-ONLY OBJECT
# ------------------------------------------------------------

tumor_expression_object <- list(

  created_at =
    Sys.time(),

  source =
    "GSE76427",

  cohort =
    "Primary HCC tumor",

  n_samples =
    ncol(tumor_expression),

  n_genes =
    nrow(tumor_expression),

  sample_names =
    colnames(tumor_expression),

  gene_expression =
    tumor_expression
)


tumor_rds <- file.path(
  processed_dir,
  "GSE76427_tumor_expression.rds"
)


saveRDS(
  tumor_expression_object,
  file = tumor_rds
)


# ------------------------------------------------------------
# 21. SAVE SAMPLE MANIFEST OBJECT
# ------------------------------------------------------------

manifest_object <- list(

  created_at =
    Sys.time(),

  total_samples =
    total_n,

  primary_tumor_samples =
    tumor_n,

  adjacent_non_tumor_samples =
    antt_n,

  sample_manifest =
    sample_manifest
)


saveRDS(
  manifest_object,
  file = file.path(
    processed_dir,
    "GSE76427_sample_manifest.rds"
  )
)


# ------------------------------------------------------------
# 22. WRITE REPORT
# ------------------------------------------------------------

report_lines <- c(

  "============================================================",
  "GSE76427 SAMPLE / PHENOTYPE MANIFEST REPORT",
  "============================================================",
  "",

  paste(
    "Date:",
    Sys.time()
  ),

  "",

  "COHORT",
  "------------------------------------------------------------",

  paste(
    "Total samples:",
    total_n
  ),

  paste(
    "Primary HCC tumor:",
    tumor_n
  ),

  paste(
    "Adjacent non-tumor:",
    antt_n
  ),

  "",

  "SAMPLE CLASSIFICATION",
  "------------------------------------------------------------",

  "PT## → Primary HCC tumor",

  "ANTT## → Adjacent non-tumor",

  "Unclassified samples: 0",

  "",

  "VALIDATION COHORT",
  "------------------------------------------------------------",

  "Primary HCC tumor samples only",

  paste(
    "Validation cohort n =",
    tumor_n
  ),

  "",

  "PAIRING",
  "------------------------------------------------------------",

  paste(
    "Patients represented by both PT and ANTT:",
    length(paired_patient_numbers)
  ),

  paste(
    "Patients represented by one sample:",
    length(single_sample_patient_numbers)
  ),

  paste(
    "PT/ANTT pairing validation:",
    ifelse(
      patient_pair_valid,
      "PASSED",
      "FAILED"
    )
  ),

  "",

  "EXPRESSION",
  "------------------------------------------------------------",

  paste(
    "Genes:",
    nrow(tumor_expression)
  ),

  paste(
    "Tumor samples:",
    ncol(tumor_expression)
  ),

  "",

  "NEXT STEP",
  "------------------------------------------------------------",

  "Calculate GSVA pathway scores in the tumor-only",
  "validation cohort."
)


writeLines(
  report_lines,
  con = file.path(
    results_dir,
    "05_GSE76427_sample_manifest_report.txt"
  )
)


# ------------------------------------------------------------
# 23. FINAL STATUS
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("STEP 5 COMPLETED SUCCESSFULLY\n")
cat("============================================================\n\n")

cat(
  "Total GSE76427 samples: ",
  total_n,
  "\n",
  sep = ""
)

cat(
  "Primary HCC tumor samples: ",
  tumor_n,
  "\n",
  sep = ""
)

cat(
  "Adjacent non-tumor samples: ",
  antt_n,
  "\n",
  sep = ""
)

cat(
  "Unclassified samples: 0\n"
)

cat(
  "Expression/manifest order: VERIFIED\n"
)

cat(
  "PT/ANTT pairing: VERIFIED\n\n"
)

cat(
  "Validation cohort: PRIMARY HCC TUMORS ONLY\n"
)

cat(
  "Validation cohort n = ",
  tumor_n,
  "\n\n",
  sep = ""
)

cat(
  "Tumor-only expression matrix saved to:\n",
  tumor_rds,
  "\n\n",
  sep = ""
)

cat(
  "Next step: GSVA pathway scoring.\n"
)
