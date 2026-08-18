# ============================================================
# STEP 9 — FINAL PATHWAY COLLECTION, GENE-SET QC, GSVA
# ============================================================
#
# MAFLD-HCC in silico analysis
#
# PURPOSE:
#   1. Load Step 8 harmonized matrices
#   2. Load common gene universe
#   3. Load the 17 curated GMT pathways
#   4. Perform effective gene-set QC
#   5. Apply predefined minimum effective gene-set size
#   6. Construct final pathway collection
#   7. Run GSVA separately in TCGA-LIHC and GSE76427
#   8. Perform pathway-score QC
#   9. Save pathway-score matrices and QC outputs
#
# IMPORTANT:
#   - GMT files are located recursively under:
#       data/reference/gene_sets/
#   - No arbitrary 50% coverage rule is used.
#   - Minimum effective gene-set size = 10 genes.
#   - Gene identifiers are SYMBOLS.
#   - The SAME final gene-set definitions are used in both
#     cohorts.
#   - The script will stop rather than silently substitute
#     an incorrect pathway.
#
# ============================================================


# ============================================================
# 9.0 — INITIAL SETUP
# ============================================================

options(stringsAsFactors = FALSE)

cat("\n============================================================\n")
cat("STEP 9 — FINAL PATHWAY COLLECTION + GSVA\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# Required packages
# ------------------------------------------------------------

required_packages <- c(
  "GSEABase",
  "GSVA"
)

for (pkg in required_packages) {

  if (!requireNamespace(pkg, quietly = TRUE)) {

    stop(
      paste0(
        "Required package '",
        pkg,
        "' is not installed.\n",
        "Install it before running Step 9."
      )
    )

  }
}


# ============================================================
# 9.1 — DIRECTORY STRUCTURE
# ============================================================

data_dir <- "data"

processed_dir <- file.path(
  data_dir,
  "processed"
)

# IMPORTANT:
# GMT files are stored here, NOT in data/gmt/.
gmt_root_dir <- file.path(
  data_dir,
  "reference",
  "gene_sets"
)

results_dir <- "results"


dir.create(
  data_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

dir.create(
  processed_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

dir.create(
  results_dir,
  showWarnings = FALSE,
  recursive = TRUE
)


# ------------------------------------------------------------
# Confirm GMT root exists
# ------------------------------------------------------------

if (!dir.exists(gmt_root_dir)) {

  stop(
    paste0(
      "\nGMT root directory not found:\n",
      normalizePath(
        gmt_root_dir,
        winslash = "/",
        mustWork = FALSE
      ),
      "\n\nExpected structure:\n",
      "data/reference/gene_sets/\n",
      "with nested folders such as lipid/, hr/, ddr/, etc."
    )
  )

}


cat(
  "GMT root directory:\n",
  normalizePath(
    gmt_root_dir,
    winslash = "/",
    mustWork = FALSE
  ),
  "\n\n"
)


# ============================================================
# 9.2 — PREDEFINED EFFECTIVE GENE-SET SIZE
# ============================================================

# A pathway must have at least 10 genes represented in:
#
#   1. TCGA-LIHC
#   2. GSE76427
#   3. common gene universe
#
# This is the predefined effective-size criterion.
#
# IMPORTANT:
# This is NOT a 50% coverage rule.

MIN_EFFECTIVE_GENE_SET_SIZE <- 10


cat(
  "Minimum effective gene-set size:",
  MIN_EFFECTIVE_GENE_SET_SIZE,
  "genes\n\n"
)


# ============================================================
# 9.3 — LOAD STEP 8 HARMONIZED MATRICES
# ============================================================

tcga_file <- file.path(
  processed_dir,
  "TCGA-LIHC",
  "TCGA_LIHC_voom_symbol_expression_harmonized.rds"
)


gse_file <- file.path(
  processed_dir,
  "GSE76427",
  "GSE76427_tumor_symbol_expression_harmonized.rds"
)


common_file <- file.path(
  processed_dir,
  "common_gene_universe.rds"
)


required_files <- c(
  tcga_file,
  gse_file,
  common_file
)


missing_files <- required_files[
  !file.exists(required_files)
]


if (length(missing_files) > 0) {

  stop(
    paste0(
      "\nThe following required Step 8 files are missing:\n\n",
      paste(
        missing_files,
        collapse = "\n"
      ),
      "\n\nRun Step 8 successfully before running Step 9."
    )
  )

}


# ------------------------------------------------------------
# Load objects
# ------------------------------------------------------------

tcga_expr <- readRDS(
  tcga_file
)


gse_expr <- readRDS(
  gse_file
)


common_genes <- readRDS(
  common_file
)


# ============================================================
# 9.4 — STANDARDIZE EXPRESSION MATRICES
# ============================================================

if (!is.matrix(tcga_expr)) {

  tcga_expr <- as.matrix(
    tcga_expr
  )

}


if (!is.matrix(gse_expr)) {

  gse_expr <- as.matrix(
    gse_expr
  )

}


storage.mode(tcga_expr) <- "numeric"

storage.mode(gse_expr) <- "numeric"


# ------------------------------------------------------------
# Standardize gene symbols
# ------------------------------------------------------------

rownames(tcga_expr) <- toupper(
  trimws(
    rownames(tcga_expr)
  )
)


rownames(gse_expr) <- toupper(
  trimws(
    rownames(gse_expr)
  )
)


# ------------------------------------------------------------
# Remove invalid gene names
# ------------------------------------------------------------

tcga_keep <- (
  !is.na(
    rownames(tcga_expr)
  )
  &
  rownames(tcga_expr) != ""
)


gse_keep <- (
  !is.na(
    rownames(gse_expr)
  )
  &
  rownames(gse_expr) != ""
)


tcga_expr <- tcga_expr[
  tcga_keep,
  ,
  drop = FALSE
]


gse_expr <- gse_expr[
  gse_keep,
  ,
  drop = FALSE
]


# ============================================================
# 9.5 — COLLAPSE DUPLICATE GENE SYMBOLS
# ============================================================

collapse_duplicate_genes <- function(
  expr_matrix
) {

  gene_names <- rownames(
    expr_matrix
  )


  duplicated_genes <- unique(
    gene_names[
      duplicated(gene_names)
    ]
  )


  if (
    length(duplicated_genes) == 0
  ) {

    return(
      expr_matrix
    )

  }


  keep_indices <- rep(
    TRUE,
    nrow(expr_matrix)
  )


  for (g in duplicated_genes) {

    idx <- which(
      gene_names == g
    )


    mean_expression <- rowMeans(
      expr_matrix[
        idx,
        ,
        drop = FALSE
      ],
      na.rm = TRUE
    )


    best_idx <- idx[
      which.max(
        mean_expression
      )
    ]


    keep_indices[idx] <- FALSE

    keep_indices[best_idx] <- TRUE

  }


  expr_matrix[
    keep_indices,
    ,
    drop = FALSE
  ]

}


tcga_expr <- collapse_duplicate_genes(
  tcga_expr
)


gse_expr <- collapse_duplicate_genes(
  gse_expr
)


# ============================================================
# 9.6 — STANDARDIZE COMMON GENE UNIVERSE
# ============================================================

if (
  is.data.frame(common_genes)
) {

  possible_gene_columns <- c(
    "gene",
    "Gene",
    "symbol",
    "Symbol",
    "gene_symbol",
    "GeneSymbol"
  )


  found_column <- intersect(
    possible_gene_columns,
    colnames(common_genes)
  )


  if (
    length(found_column) == 0
  ) {

    common_genes <- common_genes[[1]]

  } else {

    common_genes <- common_genes[[found_column[1]]]

  }

}


common_genes <- unique(
  toupper(
    trimws(
      as.character(
        common_genes
      )
    )
  )
)


common_genes <- common_genes[
  !is.na(
    common_genes
  )
  &
  common_genes != ""
]


# ============================================================
# 9.7 — COMMON GENE UNIVERSE AUDIT
# ============================================================

tcga_genes <- rownames(
  tcga_expr
)


gse_genes <- rownames(
  gse_expr
)


common_genes_observed <- Reduce(
  intersect,
  list(
    common_genes,
    tcga_genes,
    gse_genes
  )
)


common_genes_observed <- sort(
  unique(
    common_genes_observed
  )
)


cat(
  "\n------------------------------------------------------------\n"
)

cat(
  "STEP 9 INPUT AUDIT\n"
)

cat(
  "------------------------------------------------------------\n"
)


cat(
  "TCGA-LIHC:",
  nrow(tcga_expr),
  "genes x",
  ncol(tcga_expr),
  "samples\n"
)


cat(
  "GSE76427:",
  nrow(gse_expr),
  "genes x",
  ncol(gse_expr),
  "samples\n"
)


cat(
  "Step 8 common gene universe:",
  length(common_genes),
  "genes\n"
)


cat(
  "Common genes present in both cohorts:",
  length(common_genes_observed),
  "\n"
)


if (
  length(common_genes_observed) == 0
) {

  stop(
    "The common gene universe contains no genes present in both cohorts."
  )

}


# ============================================================
# 9.8 — DEFINE THE 17 CURATED PATHWAYS
# ============================================================

target_pathways <- c(

  # ----------------------------------------------------------
  # Lipid metabolism
  # ----------------------------------------------------------

  "HALLMARK_FATTY_ACID_METABOLISM",

  "CHOLESTEROL_BIOSYNTHESIS",

  "FATTY_ACYL_COA_BIOSYNTHESIS",

  "MITOCHONDRIAL_FATTY_ACID_BETA_OXIDATION",

  "PEROXISOMAL_LIPID_METABOLISM",

  "SREBP_CHOLESTEROL_REGULATION",

  "PPAR_ALPHA_LIPID_REGULATION",

  "TRIGLYCERIDE_METABOLISM",


  # ----------------------------------------------------------
  # DDR / repair
  # ----------------------------------------------------------

  "HOMOLOGOUS_RECOMBINATION",

  "MISMATCH_REPAIR",

  "BASE_EXCISION_REPAIR",

  "DOUBLE_STRAND_BREAK_REPAIR",

  "NON_HOMOLOGOUS_END_JOINING",

  "HALLMARK_DNA_REPAIR",


  # ----------------------------------------------------------
  # Checkpoint / proliferation
  # ----------------------------------------------------------

  "P53_PATHWAY",

  "P53_DEPENDENT_G1_DNA_DAMAGE_RESPONSE",

  "G2M_CHECKPOINT"
)


if (
  length(target_pathways) != 17
) {

  stop(
    paste0(
      "Expected 17 target pathways but found ",
      length(target_pathways),
      "."
    )
  )

}


# ============================================================
# 9.9 — LOAD STEP 8 GMT INVENTORY
# ============================================================

inventory_file <- file.path(
  results_dir,
  "08_gene_set_file_inventory.csv"
)


if (
  !file.exists(inventory_file)
) {

  stop(
    paste0(
      "Step 8 GMT inventory not found:\n",
      inventory_file,
      "\n\nRun Step 8 first."
    )
  )

}


gene_set_inventory_df <- read.csv(
  inventory_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


required_inventory_columns <- c(
  "gmt_file",
  "pathway_name"
)


missing_inventory_columns <- setdiff(
  required_inventory_columns,
  colnames(gene_set_inventory_df)
)


if (
  length(missing_inventory_columns) > 0
) {

  stop(
    paste0(
      "The Step 8 GMT inventory is missing:\n",
      paste(
        missing_inventory_columns,
        collapse = ", "
      )
    )
  )

}


# ============================================================
# 9.10 — NORMALIZE PATHWAY NAMES
# ============================================================

normalize_pathway_name <- function(
  x
) {

  x <- toupper(
    trimws(
      as.character(
        x
      )
    )
  )


  x <- gsub(
    "&",
    " AND ",
    x
  )


  x <- gsub(
    "[^A-Z0-9]+",
    "_",
    x
  )


  x <- gsub(
    "_+",
    "_",
    x
  )


  x <- gsub(
    "^_|_$",
    "",
    x
  )


  x

}


gene_set_inventory_df$pathway_key <-
  normalize_pathway_name(
    gene_set_inventory_df$pathway_name
  )


gene_set_inventory_df$file_key <-
  normalize_pathway_name(
    tools::file_path_sans_ext(
      basename(
        gene_set_inventory_df$gmt_file
      )
    )
  )


# ============================================================
# 9.11 — DISCOVER ALL GMT FILES RECURSIVELY
# ============================================================
#
# THIS IS THE IMPORTANT CORRECTION.
#
# The GMT files are under:
#
# data/reference/gene_sets/
#
# and may be nested several folders deep.
#
# Therefore we search recursively.
#
# ============================================================

all_gmt_files <- list.files(
  gmt_root_dir,
  pattern = "\\.gmt$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)


if (
  length(all_gmt_files) == 0
) {

  stop(
    paste0(
      "\nNo .gmt files were found under:\n",
      normalizePath(
        gmt_root_dir,
        winslash = "/",
        mustWork = FALSE
      ),
      "\n\nCheck your data/reference/gene_sets/ directory."
    )
  )

}


cat(
  "\n============================================================\n"
)

cat(
  "GMT FILES DISCOVERED\n"
)

cat(
  "============================================================\n\n"
)


cat(
  "Total GMT files found:",
  length(all_gmt_files),
  "\n\n"
)


for (
  i in seq_along(
    all_gmt_files
  )
) {

  cat(
    sprintf(
      "%03d. %s\n",
      i,
      normalizePath(
        all_gmt_files[i],
        winslash = "/",
        mustWork = FALSE
      )
    )
  )

}


# ============================================================
# 9.12 — CREATE SEARCHABLE GMT FILE INVENTORY
# ============================================================

gmt_discovery_df <- data.frame(

  gmt_path =
    normalizePath(
      all_gmt_files,
      winslash = "/",
      mustWork = FALSE
    ),

  gmt_file =
    basename(
      all_gmt_files
    ),

  file_key =
    normalize_pathway_name(
      tools::file_path_sans_ext(
        basename(
          all_gmt_files
        )
      )
    ),

  stringsAsFactors = FALSE

)


# ============================================================
# 9.13 — RESOLVE GMT FILE FROM STEP 8 INVENTORY
# ============================================================
#
# This function tries:
#
#   1. exact full path
#   2. exact basename
#   3. normalized basename
#   4. normalized pathway-name match
#
# Because the GMT files are recursively discovered first,
# the function no longer depends on where inside
# data/reference/gene_sets/ the file is located.
#
# ============================================================

resolve_gmt_path <- function(
  inventory_file_value,
  pathway_name,
  discovery_df
) {

  inventory_file_value <- trimws(
    as.character(
      inventory_file_value
    )
  )


  # ----------------------------------------------------------
  # 1. If Step 8 already contains a valid full path
  # ----------------------------------------------------------

  if (
    file.exists(
      inventory_file_value
    )
  ) {

    return(
      normalizePath(
        inventory_file_value,
        winslash = "/",
        mustWork = TRUE
      )
    )

  }


  # ----------------------------------------------------------
  # 2. Exact basename match
  # ----------------------------------------------------------

  basename_value <- basename(
    inventory_file_value
  )


  exact_basename_idx <- which(
    tolower(
      discovery_df$gmt_file
    )
    ==
    tolower(
      basename_value
    )
  )


  if (
    length(exact_basename_idx) == 1
  ) {

    return(
      discovery_df$gmt_path[
        exact_basename_idx
      ]
    )

  }


  if (
    length(exact_basename_idx) > 1
  ) {

    stop(
      paste0(
        "Multiple GMT files have the same basename: ",
        basename_value
      )
    )

  }


  # ----------------------------------------------------------
  # 3. Normalized filename match
  # ----------------------------------------------------------

  normalized_inventory_file <-
    normalize_pathway_name(
      tools::file_path_sans_ext(
        basename_value
      )
    )


  normalized_idx <- which(
    discovery_df$file_key ==
      normalized_inventory_file
  )


  if (
    length(normalized_idx) == 1
  ) {

    return(
      discovery_df$gmt_path[
        normalized_idx
      ]
    )

  }


  if (
    length(normalized_idx) > 1
  ) {

    stop(
      paste0(
        "Multiple GMT files matched normalized filename: ",
        normalized_inventory_file
      )
    )

  }


  # ----------------------------------------------------------
  # 4. Match using pathway name
  # ----------------------------------------------------------

  pathway_key <- normalize_pathway_name(
    pathway_name
  )


  pathway_idx <- which(
    discovery_df$file_key ==
      pathway_key
  )


  if (
    length(pathway_idx) == 1
  ) {

    return(
      discovery_df$gmt_path[
        pathway_idx
      ]
    )

  }


  # ----------------------------------------------------------
  # 5. No match
  # ----------------------------------------------------------

  return(
    NA_character_
  )

}


# ============================================================
# 9.14 — MATCH 17 CURATED PATHWAYS TO STEP 8 INVENTORY
# ============================================================

# The Step 8 inventory is still the authority for determining
# which biological pathway each GMT represents.
#
# This section uses the pathway matching logic already
# established in Step 8/previous Step 9.
#
# ============================================================


# ------------------------------------------------------------
# Explicit aliases
# ------------------------------------------------------------

pathway_aliases <- list(

  HALLMARK_FATTY_ACID_METABOLISM = c(
    "HALLMARK_FATTY_ACID_METABOLISM"
  ),

  CHOLESTEROL_BIOSYNTHESIS = c(
    "CHOLESTEROL_BIOSYNTHESIS",
    "CHOLESTEROL_BIOSYNTHESIS_PATHWAY",
    "REACTOME_CHOLESTEROL_BIOSYNTHESIS"
  ),

  FATTY_ACYL_COA_BIOSYNTHESIS = c(
    "FATTY_ACYL_COA_BIOSYNTHESIS",
    "FATTY_ACYL_COA_SYNTHESIS"
  ),

  MITOCHONDRIAL_FATTY_ACID_BETA_OXIDATION = c(
    "MITOCHONDRIAL_FATTY_ACID_BETA_OXIDATION",
    "MITOCHONDRIAL_FATTY_ACID_OXIDATION"
  ),

  PEROXISOMAL_LIPID_METABOLISM = c(
    "PEROXISOMAL_LIPID_METABOLISM"
  ),

  SREBP_CHOLESTEROL_REGULATION = c(
    "SREBP_CHOLESTEROL_REGULATION",
    "SREBP_REGULATION",
    "SREBP_MEDIATED_REGULATION_OF_CHOLESTEROL"
  ),

  PPAR_ALPHA_LIPID_REGULATION = c(
    "PPAR_ALPHA_LIPID_REGULATION",
    "PPAR_ALPHA_REGULATION",
    "PPARA_LIPID_REGULATION"
  ),

  TRIGLYCERIDE_METABOLISM = c(
    "TRIGLYCERIDE_METABOLISM",
    "TRIACYLGLYCEROL_METABOLISM"
  ),

  HOMOLOGOUS_RECOMBINATION = c(
    "HOMOLOGOUS_RECOMBINATION",
    "REACTOME_HOMOLOGOUS_RECOMBINATION"
  ),

  MISMATCH_REPAIR = c(
    "MISMATCH_REPAIR",
    "DNA_MISMATCH_REPAIR",
    "REACTOME_MISMATCH_REPAIR"
  ),

  BASE_EXCISION_REPAIR = c(
    "BASE_EXCISION_REPAIR",
    "DNA_BASE_EXCISION_REPAIR",
    "REACTOME_BASE_EXCISION_REPAIR"
  ),

  DOUBLE_STRAND_BREAK_REPAIR = c(
    "DOUBLE_STRAND_BREAK_REPAIR",
    "DOUBLE_STRAND_BREAK_REPAIR_PATHWAY",
    "REACTOME_DOUBLE_STRAND_BREAK_REPAIR"
  ),

  NON_HOMOLOGOUS_END_JOINING = c(
    "NON_HOMOLOGOUS_END_JOINING",
    "NON_HOMOLOGOUS_END_JOINING_NHEJ",
    "NHEJ",
    "REACTOME_NON_HOMOLOGOUS_END_JOINING",
    "REACTOME_NON_HOMOLOGOUS_END_JOINING_NHEJ"
  ),

  HALLMARK_DNA_REPAIR = c(
    "HALLMARK_DNA_REPAIR"
  ),

  P53_PATHWAY = c(
    "HALLMARK_P53_PATHWAY",
    "P53_PATHWAY",
    "P53_SIGNALING_PATHWAY",
    "P53_SIGNALING"
  ),

  P53_DEPENDENT_G1_DNA_DAMAGE_RESPONSE = c(
    "P53_DEPENDENT_G1_DNA_DAMAGE_RESPONSE",
    "REACTOME_P53_DEPENDENT_G1_DNA_DAMAGE_RESPONSE"
  ),

  G2M_CHECKPOINT = c(
    "G2M_CHECKPOINT",
    "G2_M_CHECKPOINT",
    "HALLMARK_G2M_CHECKPOINT"
  )
)


# ============================================================
# 9.15 — ROBUST PATHWAY MATCHING FUNCTION
# ============================================================

find_pathway_candidates <- function(
  target_name,
  inventory_df
) {

  target_key <- normalize_pathway_name(
    target_name
  )


  # ----------------------------------------------------------
  # Exact match
  # ----------------------------------------------------------

  exact_idx <- which(
    inventory_df$pathway_key ==
      target_key
    |
    inventory_df$file_key ==
      target_key
  )


  if (
    length(exact_idx) > 0
  ) {

    return(
      unique(
        exact_idx
      )
    )

  }


  # ----------------------------------------------------------
  # Alias match
  # ----------------------------------------------------------

  aliases <- pathway_aliases[[target_name]]


  if (
    !is.null(aliases)
  ) {

    alias_keys <- normalize_pathway_name(
      aliases
    )


    alias_idx <- which(
      inventory_df$pathway_key %in%
        alias_keys
      |
      inventory_df$file_key %in%
        alias_keys
    )


    if (
      length(alias_idx) > 0
    ) {

      return(
        unique(
          alias_idx
        )
      )

    }

  }


  # ----------------------------------------------------------
  # Special cholesterol matching
  # ----------------------------------------------------------

  if (
    target_name ==
      "CHOLESTEROL_BIOSYNTHESIS"
  ) {

    candidates <- which(

      grepl(
        "CHOLESTEROL",
        inventory_df$pathway_key,
        fixed = TRUE
      )

      &

      grepl(
        "BIOSYNTH",
        inventory_df$pathway_key
      )

    )


    if (
      length(candidates) == 0
    ) {

      candidates <- which(

        grepl(
          "CHOLESTEROL",
          inventory_df$file_key,
          fixed = TRUE
        )

        &

        grepl(
          "BIOSYNTH",
          inventory_df$file_key
        )

      )

    }


    return(
      unique(
        candidates
      )
    )

  }


  # ----------------------------------------------------------
  # Special p53 matching
  # ----------------------------------------------------------

  if (
    target_name ==
      "P53_PATHWAY"
  ) {

    candidates <- which(

      inventory_df$pathway_key ==
        "HALLMARK_P53_PATHWAY"

      |

      inventory_df$file_key ==
        "HALLMARK_P53_PATHWAY"

    )


    if (
      length(candidates) > 0
    ) {

      return(
        unique(
          candidates
        )
      )

    }


    candidates <- which(

      grepl(
        "^P53",
        inventory_df$pathway_key
      )

      &

      !grepl(
        "G1",
        inventory_df$pathway_key
      )

      &

      !grepl(
        "DNA_DAMAGE",
        inventory_df$pathway_key
      )

    )


    return(
      unique(
        candidates
      )
    )

  }


  # ----------------------------------------------------------
  # General controlled matching
  # ----------------------------------------------------------

  pattern_map <- list(

    HALLMARK_FATTY_ACID_METABOLISM =
      "FATTY_ACID_METABOLISM",

    FATTY_ACYL_COA_BIOSYNTHESIS =
      "ACYL.*COA.*(BIO|SYNTH)",

    MITOCHONDRIAL_FATTY_ACID_BETA_OXIDATION =
      "MITOCHONDRIAL.*FATTY.*OXIDATION",

    PEROXISOMAL_LIPID_METABOLISM =
      "PEROXISOM.*LIPID",

    SREBP_CHOLESTEROL_REGULATION =
      "SREBP",

    PPAR_ALPHA_LIPID_REGULATION =
      "PPAR.*ALPHA|PPARA",

    TRIGLYCERIDE_METABOLISM =
      "TRIGLYCERIDE|TRIACYLGLYCEROL",

    HOMOLOGOUS_RECOMBINATION =
      "HOMOLOGOUS.*RECOMBINATION",

    MISMATCH_REPAIR =
      "MISMATCH.*REPAIR",

    BASE_EXCISION_REPAIR =
      "BASE.*EXCISION.*REPAIR",

    DOUBLE_STRAND_BREAK_REPAIR =
      "DOUBLE.*STRAND.*BREAK.*REPAIR",

    NON_HOMOLOGOUS_END_JOINING =
      "NON.*HOMOLOGOUS.*END.*JOIN|^NHEJ$",

    HALLMARK_DNA_REPAIR =
      "HALLMARK.*DNA.*REPAIR",

    P53_DEPENDENT_G1_DNA_DAMAGE_RESPONSE =
      "P53.*G1.*DNA.*DAMAGE",

    G2M_CHECKPOINT =
      "G2.*M.*CHECKPOINT|G2M.*CHECKPOINT"
  )


  pattern <- pattern_map[[target_name]]


  if (
    is.null(pattern)
  ) {

    return(
      integer(0)
    )

  }


  pathway_hits <- grep(
    pattern,
    inventory_df$pathway_key,
    ignore.case = TRUE
  )


  file_hits <- grep(
    pattern,
    inventory_df$file_key,
    ignore.case = TRUE
  )


  unique(
    c(
      pathway_hits,
      file_hits
    )
  )

}


# ============================================================
# 9.16 — MATCH ALL 17 PATHWAYS
# ============================================================

matched_indices <- rep(
  NA_integer_,
  length(target_pathways)
)


match_methods <- rep(
  NA_character_,
  length(target_pathways)
)


for (
  i in seq_along(
    target_pathways
  )
) {

  target_name <- target_pathways[i]


  candidates <- find_pathway_candidates(
    target_name,
    gene_set_inventory_df
  )


  if (
    length(candidates) == 1
  ) {

    matched_indices[i] <-
      candidates[1]

    match_methods[i] <-
      "unique_match"

  } else {

    matched_indices[i] <-
      NA_integer_


    if (
      length(candidates) == 0
    ) {

      match_methods[i] <-
        "not_found"

    } else {

      match_methods[i] <-
        "ambiguous"

    }

  }

}


# ============================================================
# 9.17 — DISPLAY MATCHING RESULTS
# ============================================================

cat(
  "\n============================================================\n"
)

cat(
  "PATHWAY MATCHING RESULTS\n"
)

cat(
  "============================================================\n\n"
)


for (
  i in seq_along(
    target_pathways
  )
) {

  cat(
    target_pathways[i],
    " -> ",
    match_methods[i],
    "\n",
    sep = ""
  )


  if (
    !is.na(
      matched_indices[i]
    )
  ) {

    idx <- matched_indices[i]


    cat(
      "    pathway_name: ",
      gene_set_inventory_df$pathway_name[idx],
      "\n",
      sep = ""
    )


    cat(
      "    GMT inventory value: ",
      gene_set_inventory_df$gmt_file[idx],
      "\n",
      sep = ""
    )

  }

}


# ============================================================
# 9.18 — STOP IF MATCHING FAILED
# ============================================================

failed_matches <- which(
  is.na(
    matched_indices
  )
)


if (
  length(failed_matches) > 0
) {

  cat(
    "\n============================================================\n"
  )

  cat(
    "PATHWAY MATCHING FAILED\n"
  )

  cat(
    "============================================================\n\n"
  )


  for (
    i in failed_matches
  ) {

    target_name <- target_pathways[i]


    candidates <- find_pathway_candidates(
      target_name,
      gene_set_inventory_df
    )


    cat(
      "\nTarget: ",
      target_name,
      "\n",
      sep = ""
    )


    cat(
      "Status: ",
      match_methods[i],
      "\n",
      sep = ""
    )


    if (
      length(candidates) > 0
    ) {

      cat(
        "Candidate pathways:\n"
      )


      for (
        idx in candidates
      ) {

        cat(
          "  - ",
          gene_set_inventory_df$pathway_name[idx],
          " | ",
          gene_set_inventory_df$gmt_file[idx],
          "\n",
          sep = ""
        )

      }

    }

  }


  stop(
    paste0(
      "\nThe following curated pathways could not be ",
      "uniquely matched:\n",
      paste(
        target_pathways[failed_matches],
        collapse = "\n"
      ),
      "\n\nThe script stopped rather than silently ",
      "substituting an incorrect pathway."
    )
  )

}


# ============================================================
# 9.19 — BUILD MATCHED INVENTORY
# ============================================================

matched_inventory_df <- gene_set_inventory_df[
  matched_indices,
  ,
  drop = FALSE
]


matched_inventory_df$target_pathway <-
  target_pathways


matched_inventory_df$match_method <-
  match_methods


rownames(
  matched_inventory_df
) <- NULL


# ============================================================
# 9.20 — RESOLVE ACTUAL GMT FILE LOCATIONS
# ============================================================

matched_inventory_df$gmt_path <- mapply(
  FUN = resolve_gmt_path,
  inventory_file_value =
    matched_inventory_df$gmt_file,
  pathway_name =
    matched_inventory_df$target_pathway,
  MoreArgs =
    list(
      discovery_df =
        gmt_discovery_df
    ),
  USE.NAMES = FALSE
)


# ------------------------------------------------------------
# Check unresolved GMT files
# ------------------------------------------------------------

missing_gmt_paths <- matched_inventory_df$target_pathway[
  is.na(
    matched_inventory_df$gmt_path
  )
]


if (
  length(missing_gmt_paths) > 0
) {

  cat(
    "\n============================================================\n"
  )

  cat(
    "GMT FILE LOCATION FAILURE\n"
  )

  cat(
    "============================================================\n\n"
  )


  for (
    pathway_name in missing_gmt_paths
  ) {

    idx <- which(
      matched_inventory_df$target_pathway ==
        pathway_name
    )


    cat(
      "Pathway: ",
      pathway_name,
      "\n",
      sep = ""
    )


    cat(
      "Step 8 inventory value: ",
      matched_inventory_df$gmt_file[idx],
      "\n\n",
      sep = ""
    )

  }


  stop(
    paste0(
      "\nSome GMT files could not be located despite ",
      "recursive searching under:\n",
      normalizePath(
        gmt_root_dir,
        winslash = "/",
        mustWork = FALSE
      )
    )
  )

}


# ============================================================
# 9.21 — SAVE PATHWAY MATCHING AUDIT
# ============================================================

write.csv(
  matched_inventory_df,
  file.path(
    results_dir,
    "09_curated_pathway_matching.csv"
  ),
  row.names = FALSE
)


cat(
  "\n============================================================\n"
)

cat(
  "GMT PATH RESOLUTION SUCCESSFUL\n"
)

cat(
  "============================================================\n\n"
)


for (
  i in seq_len(
    nrow(matched_inventory_df)
  )
) {

  cat(
    matched_inventory_df$target_pathway[i],
    " -> ",
    matched_inventory_df$gmt_path[i],
    "\n",
    sep = ""
  )

}


# ============================================================
# 9.22 — READ ALL 17 GMT FILES
# ============================================================

pathway_gene_sets_raw <- list()


for (
  i in seq_len(
    nrow(matched_inventory_df)
  )
) {

  pathway_name <-
    matched_inventory_df$target_pathway[i]


  gmt_path <-
    matched_inventory_df$gmt_path[i]


  cat(
    "\nLoading pathway:",
    pathway_name,
    "\n"
  )


  cat(
    "GMT file:",
    gmt_path,
    "\n"
  )


  gmt_obj <- GSEABase::getGmt(
    gmt_path,
    geneIdType =
      GSEABase::NullIdentifier()
  )


  if (
    length(gmt_obj) == 0
  ) {

    stop(
      paste0(
        "GMT file contains no gene sets:\n",
        gmt_path
      )
    )

  }


  pathway_set <- gmt_obj[[1]]


  pathway_genes <- unique(
    toupper(
      trimws(
        as.character(
          GSEABase::geneIds(
            pathway_set
          )
        )
      )
    )
  )


  pathway_genes <- pathway_genes[
    !is.na(
      pathway_genes
    )
    &
    pathway_genes != ""
  ]


  if (
    length(pathway_genes) == 0
  ) {

    stop(
      paste0(
        "No genes found in pathway:\n",
        pathway_name
      )
    )

  }


  pathway_gene_sets_raw[[pathway_name]] <-
    pathway_genes

}


# ============================================================
# 9.23 — CONFIRM 17 PATHWAYS WERE LOADED
# ============================================================

if (
  length(pathway_gene_sets_raw) != 17
) {

  stop(
    paste0(
      "Expected 17 pathway gene sets but loaded ",
      length(pathway_gene_sets_raw),
      "."
    )
  )

}


# ============================================================
# 9.24 — EFFECTIVE GENE-SET QC
# ============================================================

cat(
  "\n============================================================\n"
)

cat(
  "EFFECTIVE GENE-SET QC\n"
)

cat(
  "============================================================\n\n"
)


qc_rows <- vector(
  "list",
  length(
    pathway_gene_sets_raw
  )
)


pathway_names <- names(
  pathway_gene_sets_raw
)


for (
  i in seq_along(
    pathway_gene_sets_raw
  )
) {

  pathway_name <-
    pathway_names[i]


  genes <-
    pathway_gene_sets_raw[[pathway_name]]


  genes_tcga <- intersect(
    genes,
    rownames(tcga_expr)
  )


  genes_gse <- intersect(
    genes,
    rownames(gse_expr)
  )


  genes_common <- intersect(
    genes,
    common_genes_observed
  )


  qc_rows[[i]] <- data.frame(

    pathway_name =
      pathway_name,

    n_genes_total =
      length(genes),

    n_genes_TCGA =
      length(genes_tcga),

    n_genes_GSE76427 =
      length(genes_gse),

    n_genes_common =
      length(genes_common),

    adequate_TCGA =
      length(genes_tcga) >=
      MIN_EFFECTIVE_GENE_SET_SIZE,

    adequate_GSE76427 =
      length(genes_gse) >=
      MIN_EFFECTIVE_GENE_SET_SIZE,

    adequate_common =
      length(genes_common) >=
      MIN_EFFECTIVE_GENE_SET_SIZE,

    stringsAsFactors = FALSE
  )

}


gene_set_qc <- do.call(
  rbind,
  qc_rows
)


rownames(
  gene_set_qc
) <- NULL


# ------------------------------------------------------------
# Final eligibility
# ------------------------------------------------------------

gene_set_qc$pass_effective_size <- (

  gene_set_qc$n_genes_TCGA >=
    MIN_EFFECTIVE_GENE_SET_SIZE

  &

  gene_set_qc$n_genes_GSE76427 >=
    MIN_EFFECTIVE_GENE_SET_SIZE

  &

  gene_set_qc$n_genes_common >=
    MIN_EFFECTIVE_GENE_SET_SIZE

)


print(
  gene_set_qc,
  row.names = FALSE
)


# ------------------------------------------------------------
# Save complete gene-set QC
# ------------------------------------------------------------

write.csv(
  gene_set_qc,
  file.path(
    results_dir,
    "09_gene_set_effective_QC.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 9.25 — CONSTRUCT FINAL PATHWAY COLLECTION
# ============================================================

final_pathway_names <-
  gene_set_qc$pathway_name[
    gene_set_qc$pass_effective_size
  ]


excluded_pathway_names <-
  gene_set_qc$pathway_name[
    !gene_set_qc$pass_effective_size
  ]


cat(
  "\n------------------------------------------------------------\n"
)

cat(
  "FINAL PATHWAY COLLECTION\n"
)

cat(
  "------------------------------------------------------------\n"
)


cat(
  "Pathways passing effective-size QC:",
  length(
    final_pathway_names
  ),
  "/ 17\n"
)


if (
  length(excluded_pathway_names) > 0
) {

  cat(
    "\nExcluded pathways:\n"
  )


  cat(
    paste(
      excluded_pathway_names,
      collapse = "\n"
    ),
    "\n"
  )

} else {

  cat(
    "No pathways excluded by effective-size QC.\n"
  )

}


if (
  length(final_pathway_names) == 0
) {

  stop(
    "No pathway passed the predefined minimum effective gene-set size."
  )

}


# ============================================================
# 9.26 — CONSTRUCT COMMON GENE-SET DEFINITIONS
# ============================================================

# IMPORTANT:
#
# The final gene sets are restricted to the common gene
# universe.
#
# Therefore TCGA and GSE76427 receive IDENTICAL gene-set
# definitions.

final_gene_sets <- lapply(
  final_pathway_names,
  function(pathway_name) {

    genes <-
      pathway_gene_sets_raw[[pathway_name]]


    genes <- intersect(
      genes,
      common_genes_observed
    )


    genes <- sort(
      unique(
        genes
      )
    )


    if (
      length(genes) <
      MIN_EFFECTIVE_GENE_SET_SIZE
    ) {

      stop(
        paste0(
          "Internal QC failure: ",
          pathway_name,
          " has only ",
          length(genes),
          " genes after common-universe restriction."
        )
      )

    }


    return(
      genes
    )

  }
)


names(
  final_gene_sets
) <-
  final_pathway_names


# ============================================================
# 9.27 — SAVE FINAL GENE-SET SIZE TABLE
# ============================================================

final_gene_set_sizes <- data.frame(

  pathway_name =
    names(
      final_gene_sets
    ),

  n_genes_final =
    vapply(
      final_gene_sets,
      length,
      integer(1)
    ),

  stringsAsFactors = FALSE
)


write.csv(
  final_gene_set_sizes,
  file.path(
    results_dir,
    "09_final_pathway_gene_set_sizes.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 9.28 — SAVE FINAL GENE-SET COLLECTION
# ============================================================

saveRDS(
  final_gene_sets,
  file.path(
    processed_dir,
    "final_17_curated_gene_sets.rds"
  )
)


# ============================================================
# 9.29 — PREPARE GSVA EXPRESSION MATRICES
# ============================================================

tcga_gsva_genes <- intersect(
  rownames(tcga_expr),
  common_genes_observed
)


gse_gsva_genes <- intersect(
  rownames(gse_expr),
  common_genes_observed
)


tcga_gsva_expr <- tcga_expr[
  tcga_gsva_genes,
  ,
  drop = FALSE
]


gse_gsva_expr <- gse_expr[
  gse_gsva_genes,
  ,
  drop = FALSE
]


# ------------------------------------------------------------
# Remove completely NA genes
# ------------------------------------------------------------

tcga_all_na <- apply(
  tcga_gsva_expr,
  1,
  function(x) {
    all(
      is.na(x)
    )
  }
)


gse_all_na <- apply(
  gse_gsva_expr,
  1,
  function(x) {
    all(
      is.na(x)
    )
  }
)


tcga_gsva_expr <- tcga_gsva_expr[
  !tcga_all_na,
  ,
  drop = FALSE
]


gse_gsva_expr <- gse_gsva_expr[
  !gse_all_na,
  ,
  drop = FALSE
]


# ============================================================
# 9.30 — GSVA INPUT AUDIT
# ============================================================

cat(
  "\n------------------------------------------------------------\n"
)

cat(
  "GSVA INPUT MATRICES\n"
)

cat(
  "------------------------------------------------------------\n"
)


cat(
  "TCGA genes:",
  nrow(tcga_gsva_expr),
  "\n"
)


cat(
  "TCGA samples:",
  ncol(tcga_gsva_expr),
  "\n"
)


cat(
  "GSE76427 genes:",
  nrow(gse_gsva_expr),
  "\n"
)


cat(
  "GSE76427 samples:",
  ncol(gse_gsva_expr),
  "\n"
)


# ============================================================
# 9.31 — FINAL GENE-SET REPRESENTATION CHECK
# ============================================================

for (
  pathway_name in names(
    final_gene_sets
  )
) {

  n_tcga <- sum(
    final_gene_sets[[pathway_name]] %in%
      rownames(
        tcga_gsva_expr
      )
  )


  n_gse <- sum(
    final_gene_sets[[pathway_name]] %in%
      rownames(
        gse_gsva_expr
      )
  )


  if (
    n_tcga <
      MIN_EFFECTIVE_GENE_SET_SIZE

    ||

    n_gse <
      MIN_EFFECTIVE_GENE_SET_SIZE
  ) {

    stop(
      paste0(
        "Final GSVA representation failure for ",
        pathway_name,
        "."
      )
    )

  }

}


# ============================================================
# 9.32 — RUN GSVA: TCGA-LIHC
# ============================================================

cat(
  "\n============================================================\n"
)

cat(
  "RUNNING GSVA — TCGA-LIHC\n"
)

cat(
  "============================================================\n\n"
)


tcga_gsva_parameter <- GSVA::gsvaParam(
  exprData = tcga_gsva_expr,
  geneSets = final_gene_sets,
  kcdf = "Gaussian"
)


tcga_pathway_scores <- GSVA::gsva(
  tcga_gsva_parameter,
  verbose = TRUE
)


tcga_pathway_scores <- as.matrix(
  tcga_pathway_scores
)


# ============================================================
# 9.33 — RUN GSVA: GSE76427
# ============================================================

cat(
  "\n============================================================\n"
)

cat(
  "RUNNING GSVA — GSE76427\n"
)

cat(
  "============================================================\n\n"
)


gse_gsva_parameter <- GSVA::gsvaParam(
  exprData = gse_gsva_expr,
  geneSets = final_gene_sets,
  kcdf = "Gaussian"
)


gse_pathway_scores <- GSVA::gsva(
  gse_gsva_parameter,
  verbose = TRUE
)


gse_pathway_scores <- as.matrix(
  gse_pathway_scores
)


# ============================================================
# 9.34 — GSVA OUTPUT AUDIT
# ============================================================

cat(
  "\n============================================================\n"
)

cat(
  "GSVA OUTPUT AUDIT\n"
)

cat(
  "============================================================\n\n"
)


cat(
  "TCGA pathway-score dimensions:",
  nrow(
    tcga_pathway_scores
  ),
  "pathways x",
  ncol(
    tcga_pathway_scores
  ),
  "samples\n"
)


cat(
  "GSE76427 pathway-score dimensions:",
  nrow(
    gse_pathway_scores
  ),
  "pathways x",
  ncol(
    gse_pathway_scores
  ),
  "samples\n"
)


if (
  !identical(
    sort(
      rownames(
        tcga_pathway_scores
      )
    ),
    sort(
      names(
        final_gene_sets
      )
    )
  )
) {

  stop(
    "TCGA GSVA pathway names do not match final gene sets."
  )

}


if (
  !identical(
    sort(
      rownames(
        gse_pathway_scores
      )
    ),
    sort(
      names(
        final_gene_sets
      )
    )
  )
) {

  stop(
    "GSE76427 GSVA pathway names do not match final gene sets."
  )

}


# ============================================================
# 9.35 — PATHWAY-SCORE QC FUNCTION
# ============================================================

calculate_score_qc <- function(
  score_matrix
) {

  pathway_names <- rownames(
    score_matrix
  )


  qc <- lapply(
    pathway_names,
    function(pathway_name) {

      x <- as.numeric(
        score_matrix[
          pathway_name,
          ]
      )


      n_total <- length(
        x
      )


      n_missing <- sum(
        is.na(x)
      )


      n_complete <- sum(
        !is.na(x)
      )


      variance_value <- if (
        n_complete > 1
      ) {

        var(
          x,
          na.rm = TRUE
        )

      } else {

        NA_real_

      }


      sd_value <- if (
        n_complete > 1
      ) {

        sd(
          x,
          na.rm = TRUE
        )

      } else {

        NA_real_

      }


      mean_value <- if (
        n_complete > 0
      ) {

        mean(
          x,
          na.rm = TRUE
        )

      } else {

        NA_real_

      }


      median_value <- if (
        n_complete > 0
      ) {

        median(
          x,
          na.rm = TRUE
        )

      } else {

        NA_real_

      }


      min_value <- if (
        n_complete > 0
      ) {

        min(
          x,
          na.rm = TRUE
        )

      } else {

        NA_real_

      }


      max_value <- if (
        n_complete > 0
      ) {

        max(
          x,
          na.rm = TRUE
        )

      } else {

        NA_real_

      }


      effectively_constant <- (
        is.na(
          variance_value
        )
        |
        variance_value <=
          1e-12
      )


      data.frame(

        pathway_name =
          pathway_name,

        n_samples =
          n_total,

        n_nonmissing =
          n_complete,

        n_missing =
          n_missing,

        missing_fraction =
          n_missing / n_total,

        mean =
          mean_value,

        median =
          median_value,

        sd =
          sd_value,

        variance =
          variance_value,

        minimum =
          min_value,

        maximum =
          max_value,

        effectively_constant =
          effectively_constant,

        stringsAsFactors = FALSE
      )

    }
  )


  do.call(
    rbind,
    qc
  )

}


# ============================================================
# 9.36 — RUN PATHWAY-SCORE QC
# ============================================================

tcga_score_qc <- calculate_score_qc(
  tcga_pathway_scores
)


gse_score_qc <- calculate_score_qc(
  gse_pathway_scores
)


tcga_score_qc$cohort <-
  "TCGA-LIHC"


gse_score_qc$cohort <-
  "GSE76427"


pathway_score_qc <- rbind(
  tcga_score_qc,
  gse_score_qc
)


rownames(
  pathway_score_qc
) <- NULL


# ============================================================
# 9.37 — SAVE PATHWAY-SCORE QC
# ============================================================

write.csv(
  pathway_score_qc,
  file.path(
    results_dir,
    "09_pathway_score_QC.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 9.38 — CHECK FOR PROBLEMATIC PATHWAYS
# ============================================================

cat(
  "\n------------------------------------------------------------\n"
)

cat(
  "PATHWAY-SCORE QC SUMMARY\n"
)

cat(
  "------------------------------------------------------------\n"
)


cat(
  "TCGA pathways with missing scores:",
  sum(
    tcga_score_qc$n_missing > 0
  ),
  "\n"
)


cat(
  "GSE76427 pathways with missing scores:",
  sum(
    gse_score_qc$n_missing > 0
  ),
  "\n"
)


cat(
  "TCGA effectively constant pathways:",
  sum(
    tcga_score_qc$effectively_constant
  ),
  "\n"
)


cat(
  "GSE76427 effectively constant pathways:",
  sum(
    gse_score_qc$effectively_constant
  ),
  "\n"
)


if (
  any(
    tcga_score_qc$effectively_constant
  )
  |
  any(
    gse_score_qc$effectively_constant
  )
) {

  stop(
    "At least one final pathway is effectively constant after GSVA."
  )

}


if (
  any(
    tcga_score_qc$n_missing > 0
  )
  |
  any(
    gse_score_qc$n_missing > 0
  )
) {

  stop(
    "Missing GSVA pathway scores were detected."
  )

}


# ============================================================
# 9.39 — SAVE PRINCIPAL GSVA MATRICES
# ============================================================

tcga_scores_rds <- file.path(
  processed_dir,
  "TCGA-LIHC",
  "TCGA_LIHC_GSVA_pathway_scores.rds"
)


gse_scores_rds <- file.path(
  processed_dir,
  "GSE76427",
  "GSE76427_GSVA_pathway_scores.rds"
)


saveRDS(
  tcga_pathway_scores,
  tcga_scores_rds
)


saveRDS(
  gse_pathway_scores,
  gse_scores_rds
)


# ============================================================
# 9.40 — SAVE CSV COPIES
# ============================================================

write.csv(
  tcga_pathway_scores,
  file.path(
    results_dir,
    "09_TCGA_LIHC_GSVA_pathway_scores.csv"
  ),
  row.names = TRUE
)


write.csv(
  gse_pathway_scores,
  file.path(
    results_dir,
    "09_GSE76427_GSVA_pathway_scores.csv"
  ),
  row.names = TRUE
)


# ============================================================
# 9.41 — SAVE COMPLETE STEP 9 OBJECT
# ============================================================

step9_results <- list(

  tcga_pathway_scores =
    tcga_pathway_scores,

  gse76427_pathway_scores =
    gse_pathway_scores,

  final_gene_sets =
    final_gene_sets,

  final_pathway_names =
    final_pathway_names,

  excluded_pathways =
    excluded_pathway_names,

  gene_set_qc =
    gene_set_qc,

  pathway_score_qc =
    pathway_score_qc,

  pathway_matching =
    matched_inventory_df,

  minimum_effective_gene_set_size =
    MIN_EFFECTIVE_GENE_SET_SIZE,

  tcga_expression_dimensions =
    dim(
      tcga_gsva_expr
    ),

  gse76427_expression_dimensions =
    dim(
      gse_gsva_expr
    ),

  common_gene_count =
    length(
      common_genes_observed
    )
)


saveRDS(
  step9_results,
  file.path(
    processed_dir,
    "step9_GSVA_results.rds"
  )
)


# ============================================================
# 9.42 — FINAL SUMMARY
# ============================================================

cat(
  "\n============================================================\n"
)

cat(
  "STEP 9 COMPLETED SUCCESSFULLY\n"
)

cat(
  "============================================================\n\n"
)


cat(
  "Final pathways retained:",
  length(
    final_pathway_names
  ),
  "/ 17\n"
)


cat(
  "Minimum effective gene-set size:",
  MIN_EFFECTIVE_GENE_SET_SIZE,
  "\n"
)


cat(
  "TCGA-LIHC samples:",
  ncol(
    tcga_pathway_scores
  ),
  "\n"
)


cat(
  "GSE76427 samples:",
  ncol(
    gse_pathway_scores
  ),
  "\n\n"
)


cat(
  "Principal outputs:\n\n"
)


cat(
  "TCGA GSVA matrix:\n",
  tcga_scores_rds,
  "\n\n"
)


cat(
  "GSE76427 GSVA matrix:\n",
  gse_scores_rds,
  "\n\n"
)


cat(
  "Final gene sets:\n",
  file.path(
    processed_dir,
    "final_17_curated_gene_sets.rds"
  ),
  "\n\n"
)


cat(
  "Gene-set QC:\n",
  file.path(
    results_dir,
    "09_gene_set_effective_QC.csv"
  ),
  "\n\n"
)


cat(
  "Pathway-score QC:\n",
  file.path(
    results_dir,
    "09_pathway_score_QC.csv"
  ),
  "\n\n"
)


cat(
  "Pathway matching audit:\n",
  file.path(
    results_dir,
    "09_curated_pathway_matching.csv"
  ),
  "\n\n"
)


cat(
  "Combined Step 9 object:\n",
  file.path(
    processed_dir,
    "step9_GSVA_results.rds"
  ),
  "\n\n"
)


cat(
  "============================================================\n"
)

cat(
  "END OF STEP 9\n"
)

cat(
  "============================================================\n"
)
