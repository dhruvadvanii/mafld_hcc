# ============================================================
# MAFLD-HCC MANUSCRIPT PIPELINE
# STEP 1 — PROJECT SETUP + INPUT VALIDATION
# ============================================================

# ------------------------------------------------------------
# 1. SET PROJECT DIRECTORY
# ------------------------------------------------------------

project_dir <- "C:/Users/Dhruv Advani/OneDrive/Desktop/AIIMS/Research/MAFLD_HCC_Pipeline"

setwd(project_dir)

cat("\nProject directory:\n")
cat(normalizePath(getwd()), "\n\n")


# ------------------------------------------------------------
# 2. VERIFY EXISTING PROJECT STRUCTURE
# ------------------------------------------------------------

required_dirs <- c(
  "data",
  "data/raw",
  "data/raw/TCGA-LIHC",
  "data/raw/GSE76427",
  "data/reference",
  "figures",
  "scripts",
  "results",
  "reformat code"
)

missing_dirs <- required_dirs[!dir.exists(required_dirs)]

if (length(missing_dirs) > 0) {

  cat("The following expected directories are missing:\n")
  cat(paste0(" - ", missing_dirs, "\n"), sep = "")

  stop(
    "Project structure validation failed. ",
    "Check the missing directories before continuing."
  )
}

cat("Project structure verified successfully.\n\n")


# ------------------------------------------------------------
# 3. DEFINE INPUT FILES
# ------------------------------------------------------------

expected_files <- c(

  # ---------------- TCGA-LIHC ----------------

  tcga_clinical =
    file.path(
      "data/raw/TCGA-LIHC",
      "TCGA-LIHC.clinical.tsv"
    ),

  tcga_expression =
    file.path(
      "data/raw/TCGA-LIHC",
      "TCGA-LIHC.star_counts.tsv"
    ),

  tcga_survival =
    file.path(
      "data/raw/TCGA-LIHC",
      "TCGA-LIHC.survival.tsv"
    ),

  # ---------------- GSE76427 ----------------

  gse76427_directory =
    file.path(
      "data/raw/GSE76427"
    ),

  # ---------------- Reference ----------------

  gene_sets =
    file.path(
      "data/reference",
      "gene_sets.RData"
    ),

  m3_public =
    file.path(
      "data/reference",
      "m3.public"
    )
)


# ------------------------------------------------------------
# 4. CHECK TCGA-LIHC INPUT FILES
# ------------------------------------------------------------

cat("============================================================\n")
cat("TCGA-LIHC INPUT FILES\n")
cat("============================================================\n\n")

tcga_files <- expected_files[
  c(
    "tcga_clinical",
    "tcga_expression",
    "tcga_survival"
  )
]

tcga_status <- data.frame(
  file = names(tcga_files),
  path = unname(tcga_files),
  exists = file.exists(unname(tcga_files)),
  size_MB = NA_real_,
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(tcga_status))) {

  if (tcga_status$exists[i]) {

    tcga_status$size_MB[i] <-
      file.info(tcga_status$path[i])$size / 1024^2

    cat(
      "[FOUND]   ",
      tcga_status$file[i],
      " | ",
      round(tcga_status$size_MB[i], 2),
      " MB\n",
      sep = ""
    )

  } else {

    cat(
      "[MISSING] ",
      tcga_status$file[i],
      " | ",
      tcga_status$path[i],
      "\n",
      sep = ""
    )
  }
}

cat("\n")


# ------------------------------------------------------------
# 5. VERIFY GSE76427 DIRECTORY
# ------------------------------------------------------------

cat("============================================================\n")
cat("GSE76427 INPUT DIRECTORY\n")
cat("============================================================\n\n")

gse_dir <- expected_files["gse76427_directory"]

if (!dir.exists(gse_dir)) {

  stop(
    "GSE76427 directory was not found at: ",
    gse_dir
  )
}

gse_files <- list.files(
  gse_dir,
  recursive = TRUE,
  full.names = TRUE,
  include.dirs = FALSE
)

if (length(gse_files) == 0) {

  stop(
    "GSE76427 directory exists but contains no files."
  )
}

cat(
  "Found ",
  length(gse_files),
  " file(s) in GSE76427 directory:\n\n",
  sep = ""
)

for (f in gse_files) {

  size_mb <- file.info(f)$size / 1024^2

  cat(
    " - ",
    basename(f),
    " | ",
    round(size_mb, 2),
    " MB\n",
    sep = ""
  )
}

cat("\n")


# ------------------------------------------------------------
# 6. CHECK REFERENCE FILES
# ------------------------------------------------------------
# These are reported but are NOT required for Step 1.
# We will determine their necessity in the relevant
# downstream analysis step.

# NOTE:
# m3.public is intentionally NOT required at this stage.

cat("============================================================\n")
cat("REFERENCE FILE STATUS\n")
cat("============================================================\n\n")

reference_files <- expected_files[
  c(
    "gene_sets",
    "m3_public"
  )
]

reference_status <- data.frame(
  file = names(reference_files),
  path = unname(reference_files),
  exists = file.exists(unname(reference_files)),
  size_MB = NA_real_,
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(reference_status))) {

  if (reference_status$exists[i]) {

    reference_status$size_MB[i] <-
      file.info(reference_status$path[i])$size / 1024^2

    cat(
      "[FOUND]   ",
      reference_status$file[i],
      " | ",
      round(reference_status$size_MB[i], 2),
      " MB\n",
      sep = ""
    )

  } else {

    cat(
      "[NOT YET PRESENT] ",
      reference_status$file[i],
      " | ",
      reference_status$path[i],
      "\n",
      sep = ""
    )
  }
}

cat("\n")


# ------------------------------------------------------------
# 7. VERIFY CORE TCGA INPUTS
# ------------------------------------------------------------

core_tcga <- c(
  "tcga_clinical",
  "tcga_expression",
  "tcga_survival"
)

missing_core <- core_tcga[
  !tcga_status$exists[
    match(core_tcga, tcga_status$file)
  ]
]

if (length(missing_core) > 0) {

  cat("============================================================\n")
  cat("CORE TCGA INPUTS MISSING\n")
  cat("============================================================\n\n")

  cat(
    paste0(
      " - ",
      missing_core,
      "\n"
    ),
    sep = ""
  )

  stop(
    "\nStep 1 stopped because one or more core TCGA files are missing."
  )
}


# ------------------------------------------------------------
# 8. SAVE TCGA INPUT MANIFEST
# ------------------------------------------------------------

write.csv(
  tcga_status,
  file = file.path(
    "results",
    "01_TCGA_input_manifest.csv"
  ),
  row.names = FALSE
)

cat(
  "TCGA input manifest saved to:\n",
  file.path(
    normalizePath("results"),
    "01_TCGA_input_manifest.csv"
  ),
  "\n\n"
)


# ------------------------------------------------------------
# 9. SAVE GSE76427 FILE MANIFEST
# ------------------------------------------------------------

gse_manifest <- data.frame(
  file = basename(gse_files),
  path = normalizePath(gse_files),
  size_MB = file.info(gse_files)$size / 1024^2,
  stringsAsFactors = FALSE
)

write.csv(
  gse_manifest,
  file = file.path(
    "results",
    "01_GSE76427_file_manifest.csv"
  ),
  row.names = FALSE
)

cat(
  "GSE76427 file manifest saved to:\n",
  file.path(
    normalizePath("results"),
    "01_GSE76427_file_manifest.csv"
  ),
  "\n\n"
)


# ------------------------------------------------------------
# 10. SAVE REFERENCE FILE MANIFEST
# ------------------------------------------------------------

write.csv(
  reference_status,
  file = file.path(
    "results",
    "01_reference_file_manifest.csv"
  ),
  row.names = FALSE
)

cat(
  "Reference file manifest saved to:\n",
  file.path(
    normalizePath("results"),
    "01_reference_file_manifest.csv"
  ),
  "\n\n"
)


# ------------------------------------------------------------
# 11. RECORD R SESSION INFORMATION
# ------------------------------------------------------------

capture.output(
  sessionInfo(),
  file = file.path(
    "results",
    "01_sessionInfo.txt"
  )
)

cat(
  "R session information saved.\n\n"
)


# ------------------------------------------------------------
# 12. SAVE PROJECT STATE
# ------------------------------------------------------------

project_state <- list(
  project_directory = normalizePath(getwd()),
  created_at = Sys.time(),
  expected_files = expected_files,
  tcga_status = tcga_status,
  reference_status = reference_status,
  gse_manifest = gse_manifest
)

saveRDS(
  project_state,
  file = file.path(
    "results",
    "01_project_state.rds"
  )
)

cat(
  "Project state saved.\n\n"
)


# ------------------------------------------------------------
# 13. FINAL STATUS
# ------------------------------------------------------------

cat("============================================================\n")
cat("STEP 1 COMPLETED SUCCESSFULLY\n")
cat("============================================================\n\n")

cat("Project structure: VERIFIED\n")
cat("TCGA-LIHC core files: VERIFIED\n")
cat("GSE76427 directory: VERIFIED\n")
cat("Reference files: STATUS RECORDED\n")
cat("TCGA MAF: NOT INCLUDED\n")
cat("Input manifests: SAVED\n")
cat("Session information: SAVED\n")
cat("Project state: SAVED\n\n")

cat(
  "Next step: inspect the exact GSE76427 files and determine\n",
  "which file is the expression matrix and how its probes\n",
  "map to genes.\n"
)
