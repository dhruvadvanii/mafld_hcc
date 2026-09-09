# ============================================================
# MAFLD-HCC MANUSCRIPT PIPELINE
# STEP 2 — GSE76427 INPUT AUDIT
# ============================================================
#
# PURPOSE:
#   Inspect the downloaded GSE76427 files before any
#   expression processing, normalization, or probe mapping.
#
# IMPORTANT:
#   This step does NOT modify raw files.
#   This step does NOT normalize expression data.
#   This step does NOT perform probe-to-gene conversion.
#
# OUTPUTS:
#   results/02_GSE76427_file_audit.csv
#   results/02_GSE76427_tabular_file_summary.csv
#   results/02_GSE76427_file_previews.txt
#   results/02_GSE76427_input_audit.rds
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

gse_dir <- file.path(
  "data",
  "raw",
  "GSE76427"
)

results_dir <- "results"

if (!dir.exists(gse_dir)) {
  stop(
    "GSE76427 directory does not exist:\n",
    normalizePath(
      file.path("data", "raw", "GSE76427"),
      mustWork = FALSE
    )
  )
}

if (!dir.exists(results_dir)) {
  dir.create(
    results_dir,
    recursive = TRUE
  )
}


# ------------------------------------------------------------
# 3. LIST ALL FILES
# ------------------------------------------------------------

cat("============================================================\n")
cat("GSE76427 FILE INVENTORY\n")
cat("============================================================\n\n")

gse_files <- list.files(
  gse_dir,
  recursive = TRUE,
  full.names = TRUE,
  include.dirs = FALSE
)

if (length(gse_files) == 0) {
  stop(
    "No files were found in:\n",
    normalizePath(gse_dir)
  )
}

cat(
  "Number of files found: ",
  length(gse_files),
  "\n\n",
  sep = ""
)

for (i in seq_along(gse_files)) {

  f <- gse_files[i]

  size_mb <- file.info(f)$size / 1024^2

  cat(
    sprintf(
      "%02d. %-60s %10.2f MB\n",
      i,
      basename(f),
      size_mb
    )
  )
}

cat("\n")


# ------------------------------------------------------------
# 4. CREATE FILE AUDIT TABLE
# ------------------------------------------------------------

get_extension <- function(x) {

  x <- basename(x)

  if (grepl("\\.gz$", x, ignore.case = TRUE)) {
    x <- sub("\\.gz$", "", x, ignore.case = TRUE)
  }

  if (grepl("\\.tar$", x, ignore.case = TRUE)) {
    x <- sub("\\.tar$", "", x, ignore.case = TRUE)
  }

  ext <- tools::file_ext(x)

  if (ext == "") {
    return(NA_character_)
  }

  ext
}


get_compression <- function(x) {

  if (grepl("\\.gz$", x, ignore.case = TRUE)) {
    return("gzip")
  }

  if (grepl("\\.bz2$", x, ignore.case = TRUE)) {
    return("bzip2")
  }

  if (grepl("\\.zip$", x, ignore.case = TRUE)) {
    return("zip")
  }

  return(NA_character_)
}


file_audit <- data.frame(
  file_name = basename(gse_files),
  relative_path = gse_files,
  extension = sapply(gse_files, get_extension),
  compression = sapply(gse_files, get_compression),
  size_MB = file.info(gse_files)$size / 1024^2,
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 5. CLASSIFY FILES BY NAME
# ------------------------------------------------------------

file_audit$possible_expression <- grepl(
  paste(
    c(
      "expr",
      "expression",
      "matrix",
      "series_matrix",
      "normalized",
      "normalised",
      "data_table",
      "table"
    ),
    collapse = "|"
  ),
  file_audit$file_name,
  ignore.case = TRUE
)


file_audit$possible_annotation <- grepl(
  paste(
    c(
      "annot",
      "annotation",
      "platform",
      "probe",
      "GPL",
      "mapping",
      "gene"
    ),
    collapse = "|"
  ),
  file_audit$file_name,
  ignore.case = TRUE
)


file_audit$possible_metadata <- grepl(
  paste(
    c(
      "clinical",
      "pheno",
      "metadata",
      "sample",
      "characteristic",
      "survival"
    ),
    collapse = "|"
  ),
  file_audit$file_name,
  ignore.case = TRUE
)


file_audit$possible_geo_file <- grepl(
  paste(
    c(
      "GSE",
      "GPL",
      "GSM",
      "series_matrix",
      "soft"
    ),
    collapse = "|"
  ),
  file_audit$file_name,
  ignore.case = TRUE
)


# ------------------------------------------------------------
# 6. SAVE INITIAL FILE AUDIT
# ------------------------------------------------------------

write.csv(
  file_audit,
  file = file.path(
    results_dir,
    "02_GSE76427_file_audit.csv"
  ),
  row.names = FALSE
)

cat(
  "Initial file audit saved:\n",
  file.path(
    normalizePath(results_dir),
    "02_GSE76427_file_audit.csv"
  ),
  "\n\n"
)


# ------------------------------------------------------------
# 7. PREVIEW TEXT-BASED FILES
# ------------------------------------------------------------
#
# We inspect only the first few lines.
# Large files are NOT loaded into memory at this stage.
#
# This is particularly important because GSE files can be
# very large.
# ------------------------------------------------------------

cat("============================================================\n")
cat("FILE PREVIEWS\n")
cat("============================================================\n\n")


preview_file <- file.path(
  results_dir,
  "02_GSE76427_file_previews.txt"
)

preview_connection <- file(
  preview_file,
  open = "wt"
)

writeLines(
  "============================================================",
  preview_connection
)

writeLines(
  "GSE76427 FILE PREVIEWS",
  preview_connection
)

writeLines(
  "============================================================",
  preview_connection
)

writeLines(
  "",
  preview_connection
)


# ------------------------------------------------------------
# 8. FUNCTION TO READ FIRST LINES
# ------------------------------------------------------------

read_preview <- function(
    file,
    n = 15
) {

  con <- NULL

  result <- tryCatch({

    if (grepl(
      "\\.gz$",
      file,
      ignore.case = TRUE
    )) {

      con <- gzfile(
        file,
        open = "rt"
      )

    } else {

      con <- file(
        file,
        open = "rt"
      )
    }

    lines <- readLines(
      con,
      n = n,
      warn = FALSE
    )

    close(con)
    con <- NULL

    lines

  }, error = function(e) {

    if (!is.null(con)) {
      try(close(con), silent = TRUE)
    }

    paste0(
      "[ERROR READING FILE: ",
      conditionMessage(e),
      "]"
    )
  })

  result
}


# ------------------------------------------------------------
# 9. PREVIEW EACH FILE
# ------------------------------------------------------------

for (f in gse_files) {

  writeLines(
    paste0(
      "\n\n============================================================\n",
      "FILE: ",
      basename(f),
      "\n",
      "PATH: ",
      normalizePath(f),
      "\n",
      "SIZE: ",
      round(file.info(f)$size / 1024^2, 2),
      " MB\n",
      "============================================================"
    ),
    preview_connection
  )

  # Binary files are not meaningfully previewed
  binary_ext <- c(
    "xlsx",
    "xls",
    "rds",
    "rda",
    "RData",
    "zip",
    "tar"
  )

  ext <- tolower(get_extension(f))

  if (ext %in% tolower(binary_ext)) {

    writeLines(
      "[Binary/archive file — preview skipped]",
      preview_connection
    )

    next
  }

  lines <- read_preview(
    f,
    n = 15
  )

  writeLines(
    lines,
    preview_connection
  )
}


close(preview_connection)


cat(
  "File previews saved:\n",
  file.path(
    normalizePath(results_dir),
    "02_GSE76427_file_previews.txt"
  ),
  "\n\n"
)


# ------------------------------------------------------------
# 10. IDENTIFY TABULAR FILE CANDIDATES
# ------------------------------------------------------------

tabular_extensions <- c(
  "txt",
  "tsv",
  "csv",
  "gz"
)

tabular_candidates <- gse_files[
  sapply(
    gse_files,
    function(x) {

      ext <- tolower(get_extension(x))

      ext %in% tabular_extensions
    }
  )
]


cat("============================================================\n")
cat("TABULAR FILE CANDIDATES\n")
cat("============================================================\n\n")

if (length(tabular_candidates) == 0) {

  cat(
    "No obvious tabular files were detected.\n\n"
  )

} else {

  for (f in tabular_candidates) {

    cat(
      " - ",
      basename(f),
      "\n",
      sep = ""
    )
  }

  cat("\n")
}


# ------------------------------------------------------------
# 11. FUNCTION TO GUESS DELIMITER
# ------------------------------------------------------------

guess_delimiter <- function(lines) {

  lines <- lines[
    nzchar(trimws(lines))
  ]

  if (length(lines) == 0) {
    return(NA_character_)
  }

  candidates <- c(
    "\t",
    ",",
    ";"
  )

  counts <- sapply(
    candidates,
    function(delim) {

      mean(
        lengths(
          strsplit(
            lines,
            split = delim,
            fixed = TRUE
          )
        )
      )
    }
  )

  best <- which.max(counts)

  if (length(best) == 0 ||
      counts[best] <= 1) {

    return(NA_character_)
  }

  candidates[best]
}


# ------------------------------------------------------------
# 12. FUNCTION TO TEST TABULAR STRUCTURE
# ------------------------------------------------------------

inspect_tabular_file <- function(f) {

  result <- list(
    file_name = basename(f),
    path = normalizePath(f),
    delimiter = NA_character_,
    n_columns = NA_integer_,
    n_rows_sampled = NA_integer_,
    first_column = NA_character_,
    possible_matrix = FALSE,
    possible_annotation = FALSE,
    possible_metadata = FALSE,
    error = NA_character_
  )


  # Read a limited number of lines
  lines <- tryCatch(
    read_preview(
      f,
      n = 30
    ),
    error = function(e) {
      return(
        paste0(
          "[ERROR: ",
          conditionMessage(e),
          "]"
        )
      )
    }
  )


  if (length(lines) == 1 &&
      grepl("^\\[ERROR", lines[1])) {

    result$error <- lines[1]

    return(
      as.data.frame(
        result,
        stringsAsFactors = FALSE
      )
    )
  }


  # Remove empty lines
  lines <- lines[
    nzchar(trimws(lines))
  ]


  if (length(lines) == 0) {

    result$error <- "No non-empty lines found."

    return(
      as.data.frame(
        result,
        stringsAsFactors = FALSE
      )
    )
  }


  # GEO series matrix files often contain metadata lines
  # beginning with ! or # before the actual matrix.
  data_lines <- lines[
    !grepl(
      "^[!#]",
      trimws(lines)
    )
  ]


  if (length(data_lines) == 0) {

    result$error <-
      "Only metadata/comment lines detected."

    return(
      as.data.frame(
        result,
        stringsAsFactors = FALSE
      )
    )
  }


  delim <- guess_delimiter(
    data_lines
  )

  result$delimiter <- delim


  if (is.na(delim)) {

    result$error <-
      "Could not identify a tabular delimiter."

    return(
      as.data.frame(
        result,
        stringsAsFactors = FALSE
      )
    )
  }


  split_lines <- strsplit(
    data_lines,
    split = delim,
    fixed = TRUE
  )


  column_counts <- lengths(
    split_lines
  )


  result$n_columns <- max(
    column_counts,
    na.rm = TRUE
  )

  result$n_rows_sampled <- length(
    data_lines
  )


  first_row <- split_lines[[1]]

  if (length(first_row) > 0) {

    result$first_column <- first_row[1]
  }


  # ----------------------------------------------------------
  # Heuristic classification
  # ----------------------------------------------------------

  filename_lower <- tolower(
    basename(f)
  )


  # Expression matrix:
  # Usually many columns and many rows.
  #
  # We deliberately use broad heuristics here.
  # Final identification will happen in Step 3.
  #
  result$possible_matrix <-
    result$n_columns >= 10 &&
    (
      grepl(
        "expr|expression|matrix|series_matrix|normalized|normalised",
        filename_lower
      ) ||
      result$n_columns >= 50
    )


  result$possible_annotation <-
    grepl(
      "annot|platform|probe|mapping|GPL|gene",
      basename(f),
      ignore.case = TRUE
    )


  result$possible_metadata <-
    grepl(
      "pheno|clinical|metadata|sample|characteristic|survival",
      basename(f),
      ignore.case = TRUE
    )


  as.data.frame(
    result,
    stringsAsFactors = FALSE
  )
}


# ------------------------------------------------------------
# 13. INSPECT ALL TABULAR CANDIDATES
# ------------------------------------------------------------

if (length(tabular_candidates) > 0) {

  tabular_summary_list <- lapply(
    tabular_candidates,
    inspect_tabular_file
  )

  tabular_summary <- do.call(
    rbind,
    tabular_summary_list
  )

} else {

  tabular_summary <- data.frame(
    file_name = character(),
    path = character(),
    delimiter = character(),
    n_columns = integer(),
    n_rows_sampled = integer(),
    first_column = character(),
    possible_matrix = logical(),
    possible_annotation = logical(),
    possible_metadata = logical(),
    error = character(),
    stringsAsFactors = FALSE
  )
}


# ------------------------------------------------------------
# 14. SAVE TABULAR SUMMARY
# ------------------------------------------------------------

write.csv(
  tabular_summary,
  file = file.path(
    results_dir,
    "02_GSE76427_tabular_file_summary.csv"
  ),
  row.names = FALSE
)

cat(
  "Tabular file summary saved:\n",
  file.path(
    normalizePath(results_dir),
    "02_GSE76427_tabular_file_summary.csv"
  ),
  "\n\n"
)


# ------------------------------------------------------------
# 15. PRINT POTENTIAL EXPRESSION FILES
# ------------------------------------------------------------

cat("============================================================\n")
cat("POTENTIAL EXPRESSION FILES\n")
cat("============================================================\n\n")

potential_expression <- tabular_summary[
  tabular_summary$possible_matrix,
  ,
  drop = FALSE
]

if (nrow(potential_expression) == 0) {

  cat(
    "No obvious expression matrix candidate was automatically detected.\n\n"
  )

} else {

  for (i in seq_len(nrow(potential_expression))) {

    cat(
      " - ",
      potential_expression$file_name[i],
      " | columns: ",
      potential_expression$n_columns[i],
      " | delimiter: ",
      potential_expression$delimiter[i],
      "\n",
      sep = ""
    )
  }

  cat("\n")
}


# ------------------------------------------------------------
# 16. PRINT POTENTIAL ANNOTATION FILES
# ------------------------------------------------------------

cat("============================================================\n")
cat("POTENTIAL ANNOTATION FILES\n")
cat("============================================================\n\n")

potential_annotation <- tabular_summary[
  tabular_summary$possible_annotation,
  ,
  drop = FALSE
]

if (nrow(potential_annotation) == 0) {

  cat(
    "No obvious annotation file candidate was detected.\n\n"
  )

} else {

  for (i in seq_len(nrow(potential_annotation))) {

    cat(
      " - ",
      potential_annotation$file_name[i],
      "\n",
      sep = ""
    )
  }

  cat("\n")
}


# ------------------------------------------------------------
# 17. SAVE COMPLETE AUDIT OBJECT
# ------------------------------------------------------------

gse76427_audit <- list(

  project_directory =
    normalizePath(getwd()),

  gse_directory =
    normalizePath(gse_dir),

  audited_at =
    Sys.time(),

  files =
    gse_files,

  file_audit =
    file_audit,

  tabular_summary =
    tabular_summary,

  potential_expression =
    potential_expression,

  potential_annotation =
    potential_annotation
)


saveRDS(
  gse76427_audit,
  file = file.path(
    results_dir,
    "02_GSE76427_input_audit.rds"
  )
)


# ------------------------------------------------------------
# 18. FINAL STATUS
# ------------------------------------------------------------

cat("============================================================\n")
cat("STEP 2 COMPLETED SUCCESSFULLY\n")
cat("============================================================\n\n")

cat("GSE76427 files inventoried: YES\n")
cat("Raw files modified: NO\n")
cat("Expression matrix loaded: NO\n")
cat("Normalization performed: NO\n")
cat("Probe-to-gene mapping performed: NO\n")
cat("Annotation selected: NO\n\n")

cat("Outputs saved:\n")
cat(
  " - results/02_GSE76427_file_audit.csv\n"
)
cat(
  " - results/02_GSE76427_tabular_file_summary.csv\n"
)
cat(
  " - results/02_GSE76427_file_previews.txt\n"
)
cat(
  " - results/02_GSE76427_input_audit.rds\n\n"
)

cat(
  "NEXT STEP:\n",
  "Review the generated file audit and previews to identify\n",
  "the exact expression matrix and probe annotation source.\n",
  "Do NOT begin normalization until this identification is confirmed.\n"
)
