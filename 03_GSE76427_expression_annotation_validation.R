# ============================================================
# MAFLD-HCC MANUSCRIPT PIPELINE
# STEP 3 — GSE76427 EXPRESSION + ANNOTATION VALIDATION
# ============================================================
#
# PURPOSE:
#   Validate and separate the GSE76427 Illumina expression
#   signals and detection P-values, and validate the supplied
#   GPL10558 annotation files.
#
# IMPORTANT:
#   - Raw input files are NOT modified.
#   - No normalization is performed.
#   - No log transformation is performed.
#   - No probe filtering is performed.
#   - No probe-to-gene collapsing is performed.
#
# INPUT:
#   data/raw/GSE76427/GSE76427_non-normalized.txt
#   data/raw/GSE76427/GPL10558_HumanHT-12_V4_0_R1_15002873_B.txt
#   data/raw/GSE76427/GPL10558_HumanHT-12_V4_0_R2_15002873_B.txt
#
# OUTPUT:
#   results/03_GSE76427_expression_structure.csv
#   results/03_GSE76427_annotation_comparison.csv
#   results/03_GSE76427_probe_mapping_summary.csv
#   results/03_GSE76427_value_distribution.csv
#   results/03_GSE76427_detection_pvalue_summary.csv
#   results/03_GSE76427_validation_report.txt
#   results/03_GSE76427_raw_expression.rds
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
# 2. CREATE REQUIRED DIRECTORIES
# ------------------------------------------------------------

results_dir <- "results"

processed_dir <- file.path(
  "data",
  "processed",
  "GSE76427"
)

if (!dir.exists(results_dir)) {
  dir.create(
    results_dir,
    recursive = TRUE
  )
}

if (!dir.exists(processed_dir)) {
  dir.create(
    processed_dir,
    recursive = TRUE
  )
}


# ------------------------------------------------------------
# 3. DEFINE INPUT FILES
# ------------------------------------------------------------

gse_dir <- file.path(
  "data",
  "raw",
  "GSE76427"
)

expression_file <- file.path(
  gse_dir,
  "GSE76427_non-normalized.txt"
)

annotation_r1_file <- file.path(
  gse_dir,
  "GPL10558_HumanHT-12_V4_0_R1_15002873_B.txt"
)

annotation_r2_file <- file.path(
  gse_dir,
  "GPL10558_HumanHT-12_V4_0_R2_15002873_B.txt"
)


# ------------------------------------------------------------
# 4. CHECK INPUT FILES
# ------------------------------------------------------------

input_files <- c(
  expression = expression_file,
  annotation_R1 = annotation_r1_file,
  annotation_R2 = annotation_r2_file
)

input_status <- data.frame(
  file = names(input_files),
  path = unname(input_files),
  exists = file.exists(unname(input_files)),
  size_MB = NA_real_,
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(input_status))) {

  if (input_status$exists[i]) {

    input_status$size_MB[i] <-
      file.info(input_status$path[i])$size / 1024^2

    cat(
      "[FOUND] ",
      input_status$file[i],
      " | ",
      round(input_status$size_MB[i], 2),
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

if (length(missing_inputs) > 0) {

  stop(
    "Required GSE76427 input file(s) missing:\n",
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

if (!requireNamespace(
  "data.table",
  quietly = TRUE
)) {

  stop(
    "The 'data.table' package is required for Step 3.\n",
    "Install it once using:\n",
    "install.packages('data.table')"
  )
}

library(data.table)


# ------------------------------------------------------------
# 6. READ EXPRESSION FILE HEADER ONLY
# ------------------------------------------------------------

cat("============================================================\n")
cat("EXPRESSION FILE STRUCTURE\n")
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
  "Total columns detected: ",
  length(expression_header),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 7. IDENTIFY SIGNAL AND DETECTION P-VALUE COLUMNS
# ------------------------------------------------------------

signal_columns <- grep(
  "\\.AVG_Signal$",
  expression_header
)

detection_columns <- grep(
  "\\.Detection Pval$",
  expression_header
)

id_columns <- which(
  expression_header == "ID_REF"
)


cat(
  "ID_REF columns: ",
  length(id_columns),
  "\n",
  sep = ""
)

cat(
  "AVG_Signal columns: ",
  length(signal_columns),
  "\n",
  sep = ""
)

cat(
  "Detection Pval columns: ",
  length(detection_columns),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 8. STRUCTURAL VALIDATION
# ------------------------------------------------------------

if (length(id_columns) != 1) {

  stop(
    "Expected exactly one ID_REF column."
  )
}


if (length(signal_columns) !=
    length(detection_columns)) {

  stop(
    "Number of AVG_Signal columns does not equal number ",
    "of Detection Pval columns."
  )
}


expected_samples <- 167

if (length(signal_columns) !=
    expected_samples) {

  stop(
    "Expected ",
    expected_samples,
    " samples based on the known GSE76427 structure, ",
    "but detected ",
    length(signal_columns),
    "."
  )
}


expected_total_columns <-
  1 + 2 * expected_samples


if (length(expression_header) !=
    expected_total_columns) {

  stop(
    "Expected ",
    expected_total_columns,
    " total columns but detected ",
    length(expression_header),
    "."
  )
}


# ------------------------------------------------------------
# 9. EXTRACT SAMPLE NAMES
# ------------------------------------------------------------

signal_names <- expression_header[
  signal_columns
]

detection_names <- expression_header[
  detection_columns
]


signal_samples <- sub(
  "\\.AVG_Signal$",
  "",
  signal_names
)

detection_samples <- sub(
  "\\.Detection Pval$",
  "",
  detection_names
)


if (!identical(
  signal_samples,
  detection_samples
)) {

  stop(
    "Sample ordering mismatch between AVG_Signal and ",
    "Detection Pval columns."
  )
}


sample_names <- signal_samples


# ------------------------------------------------------------
# 10. SAVE EXPRESSION STRUCTURE SUMMARY
# ------------------------------------------------------------

expression_structure <- data.frame(

  total_columns =
    length(expression_header),

  id_ref_columns =
    length(id_columns),

  avg_signal_columns =
    length(signal_columns),

  detection_pval_columns =
    length(detection_columns),

  expected_samples =
    expected_samples,

  detected_samples =
    length(sample_names),

  sample_order_verified =
    TRUE,

  stringsAsFactors = FALSE
)


write.csv(
  expression_structure,
  file = file.path(
    results_dir,
    "03_GSE76427_expression_structure.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 11. READ EXPRESSION FILE
# ------------------------------------------------------------
#
# We now load the complete file.
#
# There are 335 columns:
#   1 ID_REF
#   167 AVG_Signal
#   167 Detection Pval
#
# This is large but manageable with data.table.
# ------------------------------------------------------------

cat("\n")
cat("Reading GSE76427 expression file...\n")
cat("This may take some time.\n\n")

gse_raw <- fread(
  expression_file,
  sep = "\t",
  header = TRUE,
  check.names = FALSE,
  data.table = TRUE,
  showProgress = TRUE
)


# ------------------------------------------------------------
# 12. VERIFY LOADED DIMENSIONS
# ------------------------------------------------------------

cat("\nLoaded dimensions:\n")
cat(
  "Rows: ",
  nrow(gse_raw),
  "\n",
  sep = ""
)

cat(
  "Columns: ",
  ncol(gse_raw),
  "\n",
  sep = ""
)


if (ncol(gse_raw) !=
    expected_total_columns) {

  stop(
    "Loaded expression matrix has unexpected number of columns."
  )
}


# ------------------------------------------------------------
# 13. VERIFY PROBE IDS
# ------------------------------------------------------------

probe_ids <- gse_raw[[
  "ID_REF"
]]


cat(
  "Number of probe rows: ",
  length(probe_ids),
  "\n",
  sep = ""
)

cat(
  "Unique probe IDs: ",
  length(unique(probe_ids)),
  "\n",
  sep = ""
)


if (anyDuplicated(probe_ids) > 0) {

  warning(
    "Duplicate probe IDs detected in expression file."
  )
}


# ------------------------------------------------------------
# 14. EXTRACT SIGNAL MATRIX
# ------------------------------------------------------------

signal_matrix <- as.matrix(
  gse_raw[
    ,
    ..signal_columns
  ]
)

rownames(signal_matrix) <- probe_ids

colnames(signal_matrix) <- sample_names


# ------------------------------------------------------------
# 15. EXTRACT DETECTION P-VALUE MATRIX
# ------------------------------------------------------------

detection_matrix <- as.matrix(
  gse_raw[
    ,
    ..detection_columns
  ]
)

rownames(detection_matrix) <- probe_ids

colnames(detection_matrix) <- sample_names


# ------------------------------------------------------------
# 16. FORCE NUMERIC VALUES
# ------------------------------------------------------------

signal_matrix <- apply(
  signal_matrix,
  2,
  as.numeric
)

rownames(signal_matrix) <- probe_ids

colnames(signal_matrix) <- sample_names


detection_matrix <- apply(
  detection_matrix,
  2,
  as.numeric
)

rownames(detection_matrix) <- probe_ids

colnames(detection_matrix) <- sample_names


# ------------------------------------------------------------
# 17. CHECK NUMERIC CONVERSION
# ------------------------------------------------------------

signal_na_count <- sum(
  is.na(signal_matrix)
)

detection_na_count <- sum(
  is.na(detection_matrix)
)


cat("\nNumeric validation:\n")

cat(
  "Signal NA values: ",
  signal_na_count,
  "\n",
  sep = ""
)

cat(
  "Detection P-value NA values: ",
  detection_na_count,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 18. SIGNAL VALUE SUMMARY
# ------------------------------------------------------------

signal_values <- as.numeric(
  signal_matrix
)

signal_values_finite <- signal_values[
  is.finite(signal_values)
]


signal_summary <- data.frame(

  metric = c(
    "min",
    "1_percentile",
    "5_percentile",
    "25_percentile",
    "median",
    "75_percentile",
    "95_percentile",
    "99_percentile",
    "max",
    "mean",
    "sd",
    "NA_count",
    "negative_count",
    "zero_count"
  ),

  value = c(

    min(
      signal_values_finite,
      na.rm = TRUE
    ),

    quantile(
      signal_values_finite,
      0.01,
      na.rm = TRUE
    ),

    quantile(
      signal_values_finite,
      0.05,
      na.rm = TRUE
    ),

    quantile(
      signal_values_finite,
      0.25,
      na.rm = TRUE
    ),

    median(
      signal_values_finite,
      na.rm = TRUE
    ),

    quantile(
      signal_values_finite,
      0.75,
      na.rm = TRUE
    ),

    quantile(
      signal_values_finite,
      0.95,
      na.rm = TRUE
    ),

    quantile(
      signal_values_finite,
      0.99,
      na.rm = TRUE
    ),

    max(
      signal_values_finite,
      na.rm = TRUE
    ),

    mean(
      signal_values_finite,
      na.rm = TRUE
    ),

    sd(
      signal_values_finite,
      na.rm = TRUE
    ),

    sum(
      is.na(signal_values)
    ),

    sum(
      signal_values < 0,
      na.rm = TRUE
    ),

    sum(
      signal_values == 0,
      na.rm = TRUE
    )
  ),

  stringsAsFactors = FALSE
)


write.csv(
  signal_summary,
  file = file.path(
    results_dir,
    "03_GSE76427_value_distribution.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 19. DETECTION P-VALUE SUMMARY
# ------------------------------------------------------------

detection_values <- as.numeric(
  detection_matrix
)

detection_values_finite <- detection_values[
  is.finite(detection_values)
]


detection_summary <- data.frame(

  metric = c(
    "min",
    "1_percentile",
    "5_percentile",
    "25_percentile",
    "median",
    "75_percentile",
    "95_percentile",
    "99_percentile",
    "max",
    "mean",
    "sd",
    "NA_count",
    "values_between_0_and_1",
    "values_below_0",
    "values_above_1"
  ),

  value = c(

    min(
      detection_values_finite,
      na.rm = TRUE
    ),

    quantile(
      detection_values_finite,
      0.01,
      na.rm = TRUE
    ),

    quantile(
      detection_values_finite,
      0.05,
      na.rm = TRUE
    ),

    quantile(
      detection_values_finite,
      0.25,
      na.rm = TRUE
    ),

    median(
      detection_values_finite,
      na.rm = TRUE
    ),

    quantile(
      detection_values_finite,
      0.75,
      na.rm = TRUE
    ),

    quantile(
      detection_values_finite,
      0.95,
      na.rm = TRUE
    ),

    quantile(
      detection_values_finite,
      0.99,
      na.rm = TRUE
    ),

    max(
      detection_values_finite,
      na.rm = TRUE
    ),

    mean(
      detection_values_finite,
      na.rm = TRUE
    ),

    sd(
      detection_values_finite,
      na.rm = TRUE
    ),

    sum(
      is.na(detection_values)
    ),

    sum(
      detection_values_finite >= 0 &
      detection_values_finite <= 1
    ),

    sum(
      detection_values_finite < 0
    ),

    sum(
      detection_values_finite > 1
    )
  ),

  stringsAsFactors = FALSE
)


write.csv(
  detection_summary,
  file = file.path(
    results_dir,
    "03_GSE76427_detection_pvalue_summary.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 20. READ GPL10558 R1 ANNOTATION
# ------------------------------------------------------------

cat("\n")
cat("Reading GPL10558 R1 annotation...\n")

annotation_r1 <- fread(
  annotation_r1_file,
  sep = "\t",
  header = TRUE,
  check.names = FALSE,
  data.table = TRUE,
  skip = "[Probes]",
  showProgress = TRUE
)


# ------------------------------------------------------------
# 21. READ GPL10558 R2 ANNOTATION
# ------------------------------------------------------------

cat("Reading GPL10558 R2 annotation...\n")

annotation_r2 <- fread(
  annotation_r2_file,
  sep = "\t",
  header = TRUE,
  check.names = FALSE,
  data.table = TRUE,
  skip = "[Probes]",
  showProgress = TRUE
)


# ------------------------------------------------------------
# 22. ANNOTATION STRUCTURE
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("GPL10558 ANNOTATION STRUCTURE\n")
cat("============================================================\n\n")

cat(
  "R1 rows: ",
  nrow(annotation_r1),
  "\n",
  sep = ""
)

cat(
  "R1 columns: ",
  ncol(annotation_r1),
  "\n",
  sep = ""
)

cat(
  "R2 rows: ",
  nrow(annotation_r2),
  "\n",
  sep = ""
)

cat(
  "R2 columns: ",
  ncol(annotation_r2),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 23. VERIFY REQUIRED ANNOTATION COLUMNS
# ------------------------------------------------------------

required_annotation_columns <- c(
  "Probe_Id",
  "Symbol",
  "Entrez_Gene_ID",
  "RefSeq_ID",
  "Accession"
)


missing_r1 <- setdiff(
  required_annotation_columns,
  names(annotation_r1)
)

missing_r2 <- setdiff(
  required_annotation_columns,
  names(annotation_r2)
)


if (length(missing_r1) > 0) {

  stop(
    "R1 annotation is missing required columns:\n",
    paste(
      missing_r1,
      collapse = ", "
    )
  )
}


if (length(missing_r2) > 0) {

  stop(
    "R2 annotation is missing required columns:\n",
    paste(
      missing_r2,
      collapse = ", "
    )
  )
}


# ------------------------------------------------------------
# 24. COMPARE R1 AND R2 ANNOTATIONS
# ------------------------------------------------------------

r1_probe_ids <- annotation_r1$Probe_Id

r2_probe_ids <- annotation_r2$Probe_Id


common_annotation_probes <- intersect(
  r1_probe_ids,
  r2_probe_ids
)


r1_only <- setdiff(
  r1_probe_ids,
  r2_probe_ids
)

r2_only <- setdiff(
  r2_probe_ids,
  r1_probe_ids
)


# Compare Symbol assignments for common probes
r1_symbols <- annotation_r1[
  match(
    common_annotation_probes,
    annotation_r1$Probe_Id
  ),
  Symbol
]

r2_symbols <- annotation_r2[
  match(
    common_annotation_probes,
    annotation_r2$Probe_Id
  ),
  Symbol
]


same_symbol <- (
  r1_symbols == r2_symbols
) |
  (
    is.na(r1_symbols) &
    is.na(r2_symbols)
)


annotation_comparison <- data.frame(

  r1_probes =
    length(r1_probe_ids),

  r2_probes =
    length(r2_probe_ids),

  common_probes =
    length(common_annotation_probes),

  r1_only_probes =
    length(r1_only),

  r2_only_probes =
    length(r2_only),

  common_probes_same_symbol =
    sum(
      same_symbol,
      na.rm = TRUE
    ),

  common_probes_different_symbol =
    sum(
      !same_symbol,
      na.rm = TRUE
    ),

  stringsAsFactors = FALSE
)


write.csv(
  annotation_comparison,
  file = file.path(
    results_dir,
    "03_GSE76427_annotation_comparison.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 25. SELECT R1 AS PRIMARY ANNOTATION SOURCE
# ------------------------------------------------------------
#
# We use R1 provisionally because it is the first supplied
# GPL10558 annotation file.
#
# This is NOT yet the final probe-collapsing decision.
#
# R1 vs R2 concordance will be reported before proceeding.
# ------------------------------------------------------------

annotation <- annotation_r1


# ------------------------------------------------------------
# 26. MATCH EXPRESSION PROBES TO ANNOTATION
# ------------------------------------------------------------

annotation_match <- match(
  probe_ids,
  annotation$Probe_Id
)


probe_mapping <- data.frame(

  probe_id =
    probe_ids,

  symbol =
    annotation$Symbol[
      annotation_match
    ],

  entrez_gene_id =
    annotation$Entrez_Gene_ID[
      annotation_match
    ],

  refseq_id =
    annotation$RefSeq_ID[
      annotation_match
    ],

  accession =
    annotation$Accession[
      annotation_match
    ],

  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 27. CLEAN EMPTY SYMBOLS
# ------------------------------------------------------------

probe_mapping$symbol <- trimws(
  probe_mapping$symbol
)

probe_mapping$symbol[
  probe_mapping$symbol == ""
] <- NA


# ------------------------------------------------------------
# 28. PROBE MAPPING SUMMARY
# ------------------------------------------------------------

total_expression_probes <-
  nrow(probe_mapping)


mapped_symbol_probes <- sum(
  !is.na(probe_mapping$symbol)
)


unique_symbols <- length(
  unique(
    probe_mapping$symbol[
      !is.na(probe_mapping$symbol)
    ]
  )
)


unmapped_probes <-
  total_expression_probes -
  mapped_symbol_probes


probe_mapping_summary <- data.frame(

  metric = c(
    "total_expression_probes",
    "probes_with_gene_symbol",
    "probes_without_gene_symbol",
    "unique_gene_symbols",
    "percent_probes_with_gene_symbol"
  ),

  value = c(

    total_expression_probes,

    mapped_symbol_probes,

    unmapped_probes,

    unique_symbols,

    100 *
      mapped_symbol_probes /
      total_expression_probes
  ),

  stringsAsFactors = FALSE
)


write.csv(
  probe_mapping_summary,
  file = file.path(
    results_dir,
    "03_GSE76427_probe_mapping_summary.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 29. CHECK DUPLICATE PROBE → GENE MAPPINGS
# ------------------------------------------------------------

symbol_counts <- table(
  probe_mapping$symbol,
  useNA = "no"
)


duplicate_symbol_count <- sum(
  symbol_counts > 1
)


maximum_probes_per_gene <- ifelse(
  length(symbol_counts) > 0,
  max(symbol_counts),
  NA
)


cat("\n")
cat("Probe → gene mapping:\n")

cat(
  "Total probes: ",
  total_expression_probes,
  "\n",
  sep = ""
)

cat(
  "Probes with gene symbol: ",
  mapped_symbol_probes,
  "\n",
  sep = ""
)

cat(
  "Unique gene symbols: ",
  unique_symbols,
  "\n",
  sep = ""
)

cat(
  "Genes represented by >1 probe: ",
  duplicate_symbol_count,
  "\n",
  sep = ""
)

cat(
  "Maximum probes representing one gene: ",
  maximum_probes_per_gene,
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------
# 30. DETECTION P-VALUE QC
# ------------------------------------------------------------
#
# We do NOT filter probes at this stage.
#
# We simply quantify detection rates.
# ------------------------------------------------------------

detection_threshold <- 0.05


detected_fraction_per_sample <- colMeans(
  detection_matrix < detection_threshold,
  na.rm = TRUE
)


detection_sample_summary <- data.frame(

  sample =
    names(detected_fraction_per_sample),

  fraction_probes_detected =
    as.numeric(
      detected_fraction_per_sample
    ),

  percent_probes_detected =
    100 *
    as.numeric(
      detected_fraction_per_sample
    ),

  stringsAsFactors = FALSE
)


write.csv(
  detection_sample_summary,
  file = file.path(
    results_dir,
    "03_GSE76427_detection_by_sample.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 31. SIGNAL SUMMARY BY SAMPLE
# ------------------------------------------------------------

sample_signal_summary <- data.frame(

  sample =
    colnames(signal_matrix),

  mean =
    colMeans(
      signal_matrix,
      na.rm = TRUE
    ),

  median =
    apply(
      signal_matrix,
      2,
      median,
      na.rm = TRUE
    ),

  sd =
    apply(
      signal_matrix,
      2,
      sd,
      na.rm = TRUE
    ),

  min =
    apply(
      signal_matrix,
      2,
      min,
      na.rm = TRUE
    ),

  max =
    apply(
      signal_matrix,
      2,
      max,
      na.rm = TRUE
    ),

  stringsAsFactors = FALSE
)


write.csv(
  sample_signal_summary,
  file = file.path(
    results_dir,
    "03_GSE76427_signal_by_sample.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 32. SAVE RAW-DERIVED EXPRESSION OBJECT
# ------------------------------------------------------------
#
# These are still raw/non-normalized values.
#
# Saving them prevents repeated reading of the 138 MB file
# in subsequent steps.
# ------------------------------------------------------------

gse76427_raw_expression <- list(

  expression_file =
    normalizePath(expression_file),

  annotation_file =
    normalizePath(annotation_r1_file),

  created_at =
    Sys.time(),

  platform =
    "GPL10558 — Illumina HumanHT-12 v4",

  n_samples =
    ncol(signal_matrix),

  n_probes =
    nrow(signal_matrix),

  sample_names =
    sample_names,

  probe_ids =
    probe_ids,

  signal_matrix =
    signal_matrix,

  detection_matrix =
    detection_matrix,

  probe_mapping =
    probe_mapping,

  detection_threshold =
    detection_threshold
)


saveRDS(
  gse76427_raw_expression,
  file = file.path(
    processed_dir,
    "GSE76427_raw_expression.rds"
  )
)


# ------------------------------------------------------------
# 33. WRITE VALIDATION REPORT
# ------------------------------------------------------------

report_file <- file.path(
  results_dir,
  "03_GSE76427_validation_report.txt"
)

report_lines <- c(

  "============================================================",
  "GSE76427 EXPRESSION + ANNOTATION VALIDATION REPORT",
  "============================================================",
  "",

  paste(
    "Analysis date:",
    Sys.time()
  ),

  "",

  paste(
    "Expression file:",
    normalizePath(expression_file)
  ),

  paste(
    "Expression file size (MB):",
    round(
      file.info(expression_file)$size / 1024^2,
      2
    )
  ),

  "",

  "EXPRESSION STRUCTURE",
  "------------------------------------------------------------",

  paste(
    "Total columns:",
    ncol(gse_raw)
  ),

  paste(
    "Expression signal columns:",
    length(signal_columns)
  ),

  paste(
    "Detection P-value columns:",
    length(detection_columns)
  ),

  paste(
    "Samples:",
    ncol(signal_matrix)
  ),

  paste(
    "Probes:",
    nrow(signal_matrix)
  ),

  paste(
    "Signal/detection sample ordering verified:",
    identical(
      signal_samples,
      detection_samples
    )
  ),

  "",

  "SIGNAL VALUES",
  "------------------------------------------------------------",

  paste(
    "Minimum:",
    min(
      signal_values_finite,
      na.rm = TRUE
    )
  ),

  paste(
    "Median:",
    median(
      signal_values_finite,
      na.rm = TRUE
    )
  ),

  paste(
    "Maximum:",
    max(
      signal_values_finite,
      na.rm = TRUE
    )
  ),

  paste(
    "NA values:",
    signal_na_count
  ),

  paste(
    "Negative values:",
    sum(
      signal_values < 0,
      na.rm = TRUE
    )
  ),

  "",

  "DETECTION P-VALUES",
  "------------------------------------------------------------",

  paste(
    "Minimum:",
    min(
      detection_values_finite,
      na.rm = TRUE
    )
  ),

  paste(
    "Median:",
    median(
      detection_values_finite,
      na.rm = TRUE
    )
  ),

  paste(
    "Maximum:",
    max(
      detection_values_finite,
      na.rm = TRUE
    )
  ),

  paste(
    "Values within [0,1]:",
    sum(
      detection_values_finite >= 0 &
      detection_values_finite <= 1
    )
  ),

  "",

  "ANNOTATION",
  "------------------------------------------------------------",

  paste(
    "R1 annotation probes:",
    nrow(annotation_r1)
  ),

  paste(
    "R2 annotation probes:",
    nrow(annotation_r2)
  ),

  paste(
    "Common R1/R2 probes:",
    length(common_annotation_probes)
  ),

  paste(
    "Expression probes with gene symbol:",
    mapped_symbol_probes
  ),

  paste(
    "Unique gene symbols:",
    unique_symbols
  ),

  paste(
    "Genes represented by >1 probe:",
    duplicate_symbol_count
  ),

  "",

  "PROCESSING STATUS",
  "------------------------------------------------------------",

  "Normalization: NOT PERFORMED",

  "Log transformation: NOT PERFORMED",

  "Probe filtering: NOT PERFORMED",

  "Probe collapsing: NOT PERFORMED",

  "Raw input files modified: NO",

  "",

  "NEXT STEP:",
  "Determine the appropriate normalization/transformation",
  "strategy and define probe-level QC before generating",
  "the final gene-level expression matrix."
)


writeLines(
  report_lines,
  con = report_file
)


# ------------------------------------------------------------
# 34. CLEAN MEMORY
# ------------------------------------------------------------

rm(
  gse_raw
)

gc()


# ------------------------------------------------------------
# 35. FINAL STATUS
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("STEP 3 COMPLETED SUCCESSFULLY\n")
cat("============================================================\n\n")

cat(
  "GSE76427 samples verified: ",
  ncol(signal_matrix),
  "\n",
  sep = ""
)

cat(
  "GSE76427 probes loaded: ",
  nrow(signal_matrix),
  "\n",
  sep = ""
)

cat(
  "AVG_Signal / Detection Pval pairing: VERIFIED\n"
)

cat(
  "GPL10558 annotation loaded: YES\n"
)

cat(
  "R1/R2 annotation comparison: COMPLETED\n"
)

cat(
  "Probe → gene mapping audit: COMPLETED\n"
)

cat(
  "Normalization: NOT PERFORMED\n"
)

cat(
  "Probe filtering: NOT PERFORMED\n\n"
)

cat("Important output:\n")

cat(
  " - data/processed/GSE76427/GSE76427_raw_expression.rds\n"
)

cat("\nResults:\n")

cat(
  " - results/03_GSE76427_expression_structure.csv\n"
)

cat(
  " - results/03_GSE76427_annotation_comparison.csv\n"
)

cat(
  " - results/03_GSE76427_probe_mapping_summary.csv\n"
)

cat(
  " - results/03_GSE76427_value_distribution.csv\n"
)

cat(
  " - results/03_GSE76427_detection_pvalue_summary.csv\n"
)

cat(
  " - results/03_GSE76427_detection_by_sample.csv\n"
)

cat(
  " - results/03_GSE76427_signal_by_sample.csv\n"
)

cat(
  " - results/03_GSE76427_validation_report.txt\n"
)

cat("\n")
cat(
  "NEXT STEP: determine normalization and QC strategy.\n"
)
