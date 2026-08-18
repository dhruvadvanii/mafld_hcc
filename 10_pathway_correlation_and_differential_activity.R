# ============================================================
# STEP 10 — PATHWAY CORRELATION AND DIFFERENTIAL ACTIVITY
#
# STARTING POINT:
#   Step 9 is the LAST SUCCESSFULLY COMPLETED STEP.
#
# This script:
#
#   2.3 PATHWAY CORRELATION ANALYSIS
#       - Spearman correlation between every lipid pathway
#         and every DDR pathway
#       - Spearman rho
#       - p-value
#       - Benjamini-Hochberg FDR
#       - approximate 95% CI
#       - correlation heatmap
#       - FDR/significance annotation matrix
#       - representative scatter plots
#
#   2.4 DIFFERENTIAL PATHWAY ACTIVITY
#       - lipid-high vs lipid-low groups
#       - Wilcoxon rank-sum test
#       - effect size
#       - p-value
#       - BH-FDR
#       - correct rank-biserial effect size
#       - effect direction
#       - differential lollipop plot
#       - violin/box plots
#
# IMPORTANT:
#   HALLMARK_FATTY_ACID_METABOLISM is used to DEFINE the
#   lipid-high / lipid-low groups.
#
#   It is therefore NOT interpreted as an independent
#   differential biological finding.
#
# IMPORTANT SYNTAX SAFEGUARD:
#   All double-bracket indexing uses the valid R form:
#
#       object[[nm]]
#
#   The two opening brackets are deliberately kept together.
#
# ============================================================


options(
  stringsAsFactors = FALSE,
  warn = 1
)


cat("\n")
cat("============================================================\n")
cat("STEP 10 — PATHWAY CORRELATION AND DIFFERENTIAL ACTIVITY\n")
cat("STARTING DIRECTLY FROM SUCCESSFUL STEP 9 OUTPUTS\n")
cat("============================================================\n\n")


# ============================================================
# 10.1 — REQUIRED PACKAGES
# ============================================================

required_packages <- c(
  "ggplot2",
  "dplyr",
  "tidyr"
)


missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]


if (
  length(missing_packages) > 0
) {

  stop(
    paste0(
      "The following required R packages are missing:\n",
      paste(
        missing_packages,
        collapse = "\n"
      ),
      "\n\nInstall them before running Step 10."
    )
  )

}


library(
  ggplot2
)

library(
  dplyr
)

library(
  tidyr
)


# ============================================================
# 10.2 — DIRECTORY STRUCTURE
# ============================================================

results_dir <- "results"

tables_dir <- file.path(
  results_dir,
  "tables",
  "10_pathway_analysis"
)

figures_dir <- file.path(
  results_dir,
  "figures",
  "10_pathway_analysis"
)

processed_dir <- file.path(
  "data",
  "processed",
  "10_pathway_analysis"
)


dir.create(
  tables_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  figures_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  processed_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 10.3 — STEP 9 INPUT FILES
# ============================================================
#
# These are the successful Step 9 outputs.
#
# NO Step 10A RDS is required.
# ============================================================

tcga_step9_file <- file.path(
  "data",
  "processed",
  "TCGA-LIHC",
  "TCGA_LIHC_GSVA_pathway_scores.rds"
)


gse_step9_file <- file.path(
  "data",
  "processed",
  "GSE76427",
  "GSE76427_GSVA_pathway_scores.rds"
)


if (
  !file.exists(
    tcga_step9_file
  )
) {

  stop(
    paste0(
      "SUCCESSFUL STEP 9 TCGA FILE NOT FOUND:\n",
      tcga_step9_file,
      "\n\n",
      "Step 10 cannot proceed because the validated Step 9 ",
      "TCGA pathway-score matrix is unavailable."
    )
  )

}


if (
  !file.exists(
    gse_step9_file
  )
) {

  stop(
    paste0(
      "SUCCESSFUL STEP 9 GSE76427 FILE NOT FOUND:\n",
      gse_step9_file,
      "\n\n",
      "Step 10 cannot proceed because the validated Step 9 ",
      "GSE76427 pathway-score matrix is unavailable."
    )
  )

}


cat(
  "Step 9 TCGA file found:\n",
  tcga_step9_file,
  "\n\n"
)

cat(
  "Step 9 GSE76427 file found:\n",
  gse_step9_file,
  "\n\n"
)


# ============================================================
# 10.4 — LOAD STEP 9 OBJECTS
# ============================================================

tcga_step9_object <- readRDS(
  tcga_step9_file
)


gse_step9_object <- readRDS(
  gse_step9_file
)


# ============================================================
# 10.5 — ROBUST GSVA MATRIX EXTRACTION
# ============================================================
#
# Step 9 may have saved:
#
#   1. the matrix directly
#
# OR
#
#   2. a list containing the matrix.
#
# This function handles both.
#
# IMPORTANT:
#   Correct R syntax is object[[nm]]
#   and NOT:
#
#       object[
#         [nm]
#       ]
#
# ============================================================

extract_gsva_matrix <- function(
  object,
  source_name
) {

  # ----------------------------------------------------------
  # CASE 1 — OBJECT IS ALREADY A MATRIX
  # ----------------------------------------------------------

  if (
    is.matrix(object)
  ) {

    mat <- object

  }

  # ----------------------------------------------------------
  # CASE 2 — OBJECT IS A DATA FRAME
  # ----------------------------------------------------------

  else if (
    is.data.frame(object)
  ) {

    mat <- as.matrix(
      object
    )

  }

  # ----------------------------------------------------------
  # CASE 3 — OBJECT IS A LIST
  # ----------------------------------------------------------

  else if (
    is.list(object)
  ) {

    possible_names <- c(

      "gsva_scores",
      "pathway_scores",
      "GSVA_scores",
      "GSVA_pathway_scores",
      "gsva_pathway_scores",
      "scores",
      "matrix",
      "pathway_matrix"

    )


    mat <- NULL


    object_names <- names(
      object
    )


    # --------------------------------------------------------
    # FIRST: TRY EXPECTED OBJECT NAMES
    # --------------------------------------------------------

    if (
      !is.null(
        object_names
      )
    ) {

      for (
        nm in possible_names
      ) {

        if (
          nm %in%
          object_names
        ) {

          candidate <- object[[nm]]


          if (
            is.matrix(candidate) ||
            is.data.frame(candidate)
          ) {

            mat <- as.matrix(
              candidate
            )

            break

          }

        }

      }

    }


    # --------------------------------------------------------
    # SECOND: INSPECT OTHER NAMED MATRIX ELEMENTS
    # --------------------------------------------------------

    if (
      is.null(mat) &&
      !is.null(object_names)
    ) {

      for (
        nm in object_names
      ) {

        candidate <- object[[nm]]


        if (
          is.matrix(candidate) ||
          is.data.frame(candidate)
        ) {

          candidate_matrix <- as.matrix(
            candidate
          )


          if (
            !is.null(
              rownames(
                candidate_matrix
              )
            ) &&
            !is.null(
              colnames(
                candidate_matrix
              )
            )
          ) {

            candidate_rows <- toupper(
              trimws(
                rownames(
                  candidate_matrix
                )
              )
            )


            expected_count <- sum(

              candidate_rows %in%
              c(

                "HALLMARK_FATTY_ACID_METABOLISM",
                "CHOLESTEROL_BIOSYNTHESIS",
                "FATTY_ACYL_COA_BIOSYNTHESIS",
                "MITOCHONDRIAL_FATTY_ACID_BETA_OXIDATION",
                "PEROXISOMAL_LIPID_METABOLISM",
                "SREBP_CHOLESTEROL_REGULATION",
                "PPAR_ALPHA_LIPID_REGULATION",
                "TRIGLYCERIDE_METABOLISM",
                "HOMOLOGOUS_RECOMBINATION",
                "MISMATCH_REPAIR",
                "BASE_EXCISION_REPAIR",
                "DOUBLE_STRAND_BREAK_REPAIR",
                "NON_HOMOLOGOUS_END_JOINING",
                "HALLMARK_DNA_REPAIR",
                "P53_PATHWAY",
                "P53_DEPENDENT_G1_DNA_DAMAGE_RESPONSE",
                "G2M_CHECKPOINT"

              )

            )


            if (
              expected_count >= 3
            ) {

              mat <- candidate_matrix

              break

            }

          }

        }

      }

    }

  }

  else {

    stop(
      paste0(
        source_name,
        " is neither a matrix, data.frame, nor list."
      )
    )

  }


  # ----------------------------------------------------------
  # FINAL MATRIX VALIDATION
  # ----------------------------------------------------------

  if (
    is.null(mat)
  ) {

    stop(
      paste0(
        "Could not identify a GSVA pathway-score matrix in ",
        source_name,
        "."
      )
    )

  }


  mat <- as.matrix(
    mat
  )


  storage.mode(
    mat
  ) <- "numeric"


  if (
    is.null(
      rownames(mat)
    )
  ) {

    stop(
      paste0(
        source_name,
        " GSVA matrix has no pathway names."
      )
    )

  }


  if (
    is.null(
      colnames(mat)
    )
  ) {

    stop(
      paste0(
        source_name,
        " GSVA matrix has no sample IDs."
      )
    )

  }


  rownames(mat) <- toupper(
    trimws(
      rownames(mat)
    )
  )


  colnames(mat) <- trimws(
    colnames(mat)
  )


  if (
    anyDuplicated(
      rownames(mat)
    ) > 0
  ) {

    stop(
      paste0(
        "Duplicate pathway names detected in ",
        source_name,
        "."
      )
    )

  }


  if (
    anyDuplicated(
      colnames(mat)
    ) > 0
  ) {

    stop(
      paste0(
        "Duplicate sample IDs detected in ",
        source_name,
        "."
      )
    )

  }


  return(
    mat
  )

}


# ============================================================
# 10.6 — EXTRACT BOTH COHORT MATRICES
# ============================================================

tcga_scores <- extract_gsva_matrix(
  tcga_step9_object,
  "TCGA-LIHC Step 9"
)


gse_scores <- extract_gsva_matrix(
  gse_step9_object,
  "GSE76427 Step 9"
)


cat(
  "TCGA dimensions: ",
  nrow(tcga_scores),
  " pathways × ",
  ncol(tcga_scores),
  " samples\n",
  sep = ""
)


cat(
  "GSE76427 dimensions: ",
  nrow(gse_scores),
  " pathways × ",
  ncol(gse_scores),
  " samples\n\n",
  sep = ""
)


# ============================================================
# 10.7 — LOCKED 17-PATHWAY COLLECTION
# ============================================================

expected_pathways <- c(

  "HALLMARK_FATTY_ACID_METABOLISM",
  "CHOLESTEROL_BIOSYNTHESIS",
  "FATTY_ACYL_COA_BIOSYNTHESIS",
  "MITOCHONDRIAL_FATTY_ACID_BETA_OXIDATION",
  "PEROXISOMAL_LIPID_METABOLISM",
  "SREBP_CHOLESTEROL_REGULATION",
  "PPAR_ALPHA_LIPID_REGULATION",
  "TRIGLYCERIDE_METABOLISM",

  "HOMOLOGOUS_RECOMBINATION",
  "MISMATCH_REPAIR",
  "BASE_EXCISION_REPAIR",
  "DOUBLE_STRAND_BREAK_REPAIR",
  "NON_HOMOLOGOUS_END_JOINING",
  "HALLMARK_DNA_REPAIR",

  "P53_PATHWAY",
  "P53_DEPENDENT_G1_DNA_DAMAGE_RESPONSE",
  "G2M_CHECKPOINT"

)


# ============================================================
# 10.8 — BIOLOGICAL GROUP DEFINITIONS
# ============================================================

lipid_pathways <- c(

  "HALLMARK_FATTY_ACID_METABOLISM",
  "CHOLESTEROL_BIOSYNTHESIS",
  "FATTY_ACYL_COA_BIOSYNTHESIS",
  "MITOCHONDRIAL_FATTY_ACID_BETA_OXIDATION",
  "PEROXISOMAL_LIPID_METABOLISM",
  "SREBP_CHOLESTEROL_REGULATION",
  "PPAR_ALPHA_LIPID_REGULATION",
  "TRIGLYCERIDE_METABOLISM"

)


ddr_pathways <- c(

  "HOMOLOGOUS_RECOMBINATION",
  "MISMATCH_REPAIR",
  "BASE_EXCISION_REPAIR",
  "DOUBLE_STRAND_BREAK_REPAIR",
  "NON_HOMOLOGOUS_END_JOINING",
  "HALLMARK_DNA_REPAIR"

)


checkpoint_pathways <- c(

  "P53_PATHWAY",
  "P53_DEPENDENT_G1_DNA_DAMAGE_RESPONSE",
  "G2M_CHECKPOINT"

)


# ============================================================
# 10.9 — VERIFY THE 17 PATHWAYS IN BOTH STEP 9 MATRICES
# ============================================================

tcga_available <- rownames(
  tcga_scores
)


gse_available <- rownames(
  gse_scores
)


missing_tcga <- setdiff(
  expected_pathways,
  tcga_available
)


missing_gse <- setdiff(
  expected_pathways,
  gse_available
)


pathway_audit <- data.frame(

  pathway_name =
    expected_pathways,

  TCGA_present =
    expected_pathways %in%
    tcga_available,

  GSE76427_present =
    expected_pathways %in%
    gse_available,

  stringsAsFactors = FALSE

)


pathway_audit$biological_group <- ifelse(

  pathway_audit$pathway_name %in%
  lipid_pathways,

  "Lipid metabolism",

  ifelse(

    pathway_audit$pathway_name %in%
    ddr_pathways,

    "DDR / repair",

    "Checkpoint / proliferation"

  )

)


write.csv(

  pathway_audit,

  file.path(
    tables_dir,
    "10_0_pathway_availability_audit.csv"
  ),

  row.names = FALSE

)


if (
  length(missing_tcga) > 0
) {

  stop(
    paste0(
      "The following required pathways are missing from TCGA Step 9 output:\n\n",
      paste(
        missing_tcga,
        collapse = "\n"
      ),
      "\n\nNo pathway was silently substituted."
    )
  )

}


if (
  length(missing_gse) > 0
) {

  stop(
    paste0(
      "The following required pathways are missing from GSE76427 Step 9 output:\n\n",
      paste(
        missing_gse,
        collapse = "\n"
      ),
      "\n\nNo pathway was silently substituted."
    )
  )

}


# ============================================================
# 10.10 — SUBSET TO EXACT 17 PATHWAYS
# ============================================================

tcga_scores <- tcga_scores[
  expected_pathways,
  ,
  drop = FALSE
]


gse_scores <- gse_scores[
  expected_pathways,
  ,
  drop = FALSE
]


cat(
  "All 17 pathways are present in both cohorts.\n\n"
)


# ============================================================
# 10.11 — PATHWAY SCORE QC
# ============================================================

calculate_pathway_qc <- function(
  score_matrix,
  cohort_name
) {

  qc <- data.frame(

    cohort =
      cohort_name,

    pathway =
      rownames(
        score_matrix
      ),

    n_samples =
      ncol(
        score_matrix
      ),

    n_nonmissing =
      rowSums(
        !is.na(
          score_matrix
        )
      ),

    n_missing =
      rowSums(
        is.na(
          score_matrix
        )
      ),

    mean =
      apply(
        score_matrix,
        1,
        mean,
        na.rm = TRUE
      ),

    median =
      apply(
        score_matrix,
        1,
        median,
        na.rm = TRUE
      ),

    sd =
      apply(
        score_matrix,
        1,
        sd,
        na.rm = TRUE
      ),

    variance =
      apply(
        score_matrix,
        1,
        var,
        na.rm = TRUE
      ),

    min =
      apply(
        score_matrix,
        1,
        min,
        na.rm = TRUE
      ),

    max =
      apply(
        score_matrix,
        1,
        max,
        na.rm = TRUE
      ),

    stringsAsFactors = FALSE

  )


  qc$pct_missing <-

    100 *
    qc$n_missing /
    qc$n_samples


  qc$effectively_constant <-

    is.na(
      qc$sd
    ) |
    qc$sd == 0


  return(
    qc
  )

}


tcga_qc <- calculate_pathway_qc(
  tcga_scores,
  "TCGA-LIHC"
)


gse_qc <- calculate_pathway_qc(
  gse_scores,
  "GSE76427"
)


pathway_qc <- rbind(
  tcga_qc,
  gse_qc
)


write.csv(

  pathway_qc,

  file.path(
    tables_dir,
    "10_1_pathway_score_QC.csv"
  ),

  row.names = FALSE

)


if (
  any(
    pathway_qc$effectively_constant
  )
) {

  problematic <- pathway_qc[
    pathway_qc$effectively_constant,
    ,
    drop = FALSE
  ]


  write.csv(

    problematic,

    file.path(
      tables_dir,
      "10_1_problematic_constant_pathways.csv"
    ),

    row.names = FALSE

  )


  stop(
    "At least one pathway is effectively constant in at least one cohort."
  )

}


# ============================================================
# 10.12 — SAVE CLEAN INPUT MATRICES
# ============================================================

saveRDS(

  tcga_scores,

  file.path(
    processed_dir,
    "10_TCGA_LIHC_GSVA_scores_used.rds"
  )

)


saveRDS(

  gse_scores,

  file.path(
    processed_dir,
    "10_GSE76427_GSVA_scores_used.rds"
  )

)


write.csv(

  tcga_scores,

  file.path(
    tables_dir,
    "10_TCGA_LIHC_GSVA_scores_used.csv"
  )

)


write.csv(

  gse_scores,

  file.path(
    tables_dir,
    "10_GSE76427_GSVA_scores_used.csv"
  )

)


# ============================================================
# 10.13 — SPEARMAN CORRELATION FUNCTION
# ============================================================

calculate_spearman <- function(
  x,
  y
) {

  complete <- complete.cases(
    x,
    y
  )


  x_complete <- x[
    complete
  ]


  y_complete <- y[
    complete
  ]


  n <- length(
    x_complete
  )


  if (
    n < 4
  ) {

    return(
      data.frame(
        rho = NA_real_,
        p_value = NA_real_,
        n = n,
        ci_lower = NA_real_,
        ci_upper = NA_real_
      )
    )

  }


  if (
    sd(
      x_complete
    ) == 0 ||
    sd(
      y_complete
    ) == 0
  ) {

    return(
      data.frame(
        rho = NA_real_,
        p_value = NA_real_,
        n = n,
        ci_lower = NA_real_,
        ci_upper = NA_real_
      )
    )

  }


  test <- suppressWarnings(

    cor.test(

      x_complete,

      y_complete,

      method = "spearman",

      exact = FALSE

    )

  )


  rho <- as.numeric(
    test$estimate
  )


  # Approximate 95% CI using Fisher z transformation

  if (
    is.finite(rho) &&
    abs(rho) < 1 &&
    n > 3
  ) {

    fisher_z <- atanh(
      rho
    )


    fisher_se <- 1 /
      sqrt(
        n - 3
      )


    ci_lower <- tanh(
      fisher_z -
      1.96 *
      fisher_se
    )


    ci_upper <- tanh(
      fisher_z +
      1.96 *
      fisher_se
    )

  } else {

    ci_lower <- NA_real_

    ci_upper <- NA_real_

  }


  data.frame(

    rho =
      rho,

    p_value =
      test$p.value,

    n =
      n,

    ci_lower =
      ci_lower,

    ci_upper =
      ci_upper

  )

}


# ============================================================
# 10.14 — ALL 48 LIPID × DDR CORRELATIONS
# ============================================================

correlation_results <- data.frame(

  lipid_pathway =
    character(0),

  ddr_pathway =
    character(0),

  rho =
    numeric(0),

  p_value =
    numeric(0),

  n =
    integer(0),

  ci_lower =
    numeric(0),

  ci_upper =
    numeric(0),

  stringsAsFactors = FALSE

)


for (
  lipid_name in lipid_pathways
) {

  for (
    ddr_name in ddr_pathways
  ) {

    result <- calculate_spearman(

      tcga_scores[
        lipid_name,
        ],

      tcga_scores[
        ddr_name,
        ]

    )


    new_row <- data.frame(

      lipid_pathway =
        lipid_name,

      ddr_pathway =
        ddr_name,

      rho =
        result$rho,

      p_value =
        result$p_value,

      n =
        result$n,

      ci_lower =
        result$ci_lower,

      ci_upper =
        result$ci_upper,

      stringsAsFactors = FALSE

    )


    correlation_results <- rbind(

      correlation_results,

      new_row

    )

  }

}


if (
  nrow(
    correlation_results
  ) != 48
) {

  stop(
    paste0(
      "Correlation QC failed: expected 48 lipid × DDR pairs but obtained ",
      nrow(
        correlation_results
      ),
      "."
    )
  )

}


# ============================================================
# 10.15 — BH FDR AND DIRECTION
# ============================================================

correlation_results$FDR <- p.adjust(

  correlation_results$p_value,

  method = "BH"

)


correlation_results$abs_rho <- abs(
  correlation_results$rho
)


correlation_results$direction <- ifelse(

  is.na(
    correlation_results$rho
  ),

  "Not estimable",

  ifelse(

    correlation_results$rho > 0,

    "Positive",

    ifelse(

      correlation_results$rho < 0,

      "Negative",

      "Zero"

    )

  )

)


correlation_results$significance <- ifelse(

  is.na(
    correlation_results$FDR
  ),

  "NA",

  ifelse(

    correlation_results$FDR < 0.001,

    "***",

    ifelse(

      correlation_results$FDR < 0.01,

      "**",

      ifelse(

        correlation_results$FDR < 0.05,

        "*",

        "ns"

      )

    )

  )

)


correlation_results$significant_FDR_0_05 <-

  correlation_results$FDR < 0.05


correlation_results <- correlation_results[
  order(
    correlation_results$FDR,
    -correlation_results$abs_rho
  ),
  ,
  drop = FALSE
]


write.csv(

  correlation_results,

  file.path(
    tables_dir,
    "10_2_3_complete_lipid_DDR_spearman_results.csv"
  ),

  row.names = FALSE

)


# ============================================================
# 10.16 — CORRELATION MATRICES
# ============================================================

rho_matrix <- matrix(

  NA_real_,

  nrow =
    length(
      lipid_pathways
    ),

  ncol =
    length(
      ddr_pathways
    ),

  dimnames = list(
    lipid_pathways,
    ddr_pathways
  )

)


p_matrix <- rho_matrix

fdr_matrix <- rho_matrix

ci_lower_matrix <- rho_matrix

ci_upper_matrix <- rho_matrix


for (
  i in seq_len(
    nrow(
      correlation_results
    )
  )
) {

  lipid_name <-
    correlation_results$lipid_pathway[i]


  ddr_name <-
    correlation_results$ddr_pathway[i]


  rho_matrix[
    lipid_name,
    ddr_name
  ] <-
    correlation_results$rho[i]


  p_matrix[
    lipid_name,
    ddr_name
  ] <-
    correlation_results$p_value[i]


  fdr_matrix[
    lipid_name,
    ddr_name
  ] <-
    correlation_results$FDR[i]


  ci_lower_matrix[
    lipid_name,
    ddr_name
  ] <-
    correlation_results$ci_lower[i]


  ci_upper_matrix[
    lipid_name,
    ddr_name
  ] <-
    correlation_results$ci_upper[i]

}


write.csv(
  rho_matrix,
  file.path(
    tables_dir,
    "10_2_3_spearman_rho_matrix.csv"
  )
)


write.csv(
  p_matrix,
  file.path(
    tables_dir,
    "10_2_3_spearman_p_value_matrix.csv"
  )
)


write.csv(
  fdr_matrix,
  file.path(
    tables_dir,
    "10_2_3_spearman_BH_FDR_matrix.csv"
  )
)


write.csv(
  ci_lower_matrix,
  file.path(
    tables_dir,
    "10_2_3_spearman_95CI_lower_matrix.csv"
  )
)


write.csv(
  ci_upper_matrix,
  file.path(
    tables_dir,
    "10_2_3_spearman_95CI_upper_matrix.csv"
  )
)


# ============================================================
# 10.17 — FDR SIGNIFICANCE ANNOTATION MATRIX
# ============================================================

significance_matrix <- matrix(

  "ns",

  nrow =
    nrow(
      fdr_matrix
    ),

  ncol =
    ncol(
      fdr_matrix
    ),

  dimnames =
    dimnames(
      fdr_matrix
    )

)


significance_matrix[
  !is.na(
    fdr_matrix
  ) &
  fdr_matrix < 0.05
] <- "*"


significance_matrix[
  !is.na(
    fdr_matrix
  ) &
  fdr_matrix < 0.01
] <- "**"


significance_matrix[
  !is.na(
    fdr_matrix
  ) &
  fdr_matrix < 0.001
] <- "***"


significance_matrix[
  is.na(
    fdr_matrix
  )
] <- "NA"


write.csv(

  significance_matrix,

  file.path(
    tables_dir,
    "10_2_3_FDR_significance_annotation_matrix.csv"
  )

)


# ============================================================
# 10.18 — CORRELATION HEATMAP
# ============================================================

heatmap_data <- correlation_results


heatmap_data$label <- ifelse(

  is.na(
    heatmap_data$rho
  ),

  "NA",

  paste0(

    sprintf(
      "%.2f",
      heatmap_data$rho
    ),

    "\n",

    heatmap_data$significance

  )

)


heatmap_plot <- ggplot(

  heatmap_data,

  aes(
    x = ddr_pathway,
    y = lipid_pathway,
    fill = rho
  )

) +

  geom_tile(
    colour = "white",
    linewidth = 0.5
  ) +

  geom_text(
    aes(
      label = label
    ),
    size = 3
  ) +

  scale_fill_gradient2(

    low = "#B2182B",

    mid = "white",

    high = "#2166AC",

    midpoint = 0,

    limits = c(
      -1,
      1
    ),

    name = "Spearman\nrho"

  ) +

  labs(

    title =
      "Lipid Pathway–DDR Pathway Correlations",

    subtitle =
      "Numbers = Spearman rho; *, **, *** = BH-FDR <0.05, <0.01, <0.001",

    x =
      "DDR / Repair Pathway",

    y =
      "Lipid Metabolism Pathway"

  ) +

  theme_minimal(
    base_size = 11
  ) +

  theme(

    axis.text.x =
      element_text(
        angle = 45,
        hjust = 1
      ),

    axis.text.y =
      element_text(
        size = 9
      ),

    plot.title =
      element_text(
        face = "bold"
      ),

    panel.grid =
      element_blank()

  )


ggsave(

  file.path(
    figures_dir,
    "10_2_3_lipid_DDR_correlation_heatmap.png"
  ),

  heatmap_plot,

  width = 11,

  height = 7,

  dpi = 300

)


ggsave(

  file.path(
    figures_dir,
    "10_2_3_lipid_DDR_correlation_heatmap.pdf"
  ),

  heatmap_plot,

  width = 11,

  height = 7

)


# ============================================================
# 10.19 — REPRESENTATIVE CORRELATION PAIRS
# ============================================================

representative_results <- data.frame()


negative_results <- correlation_results[
  correlation_results$rho < 0 &
  is.finite(
    correlation_results$rho
  ),
  ,
  drop = FALSE
]


if (
  nrow(
    negative_results
  ) > 0
) {

  strongest_negative <- negative_results[
    which.min(
      negative_results$rho
    ),
    ,
    drop = FALSE
  ]


  representative_results <- rbind(
    representative_results,
    strongest_negative
  )

}


positive_results <- correlation_results[
  correlation_results$rho > 0 &
  is.finite(
    correlation_results$rho
  ),
  ,
  drop = FALSE
]


if (
  nrow(
    positive_results
  ) > 0
) {

  strongest_positive <- positive_results[
    which.max(
      positive_results$rho
    ),
    ,
    drop = FALSE
  ]


  representative_results <- rbind(
    representative_results,
    strongest_positive
  )

}


significant_negative <- correlation_results[
  correlation_results$rho < 0 &
  correlation_results$FDR < 0.05,
  ,
  drop = FALSE
]


if (
  nrow(
    significant_negative
  ) > 0
) {

  strongest_significant_negative <-
    significant_negative[
      which.min(
        significant_negative$rho
      ),
      ,
      drop = FALSE
    ]


  representative_results <- rbind(
    representative_results,
    strongest_significant_negative
  )

}


significant_positive <- correlation_results[
  correlation_results$rho > 0 &
  correlation_results$FDR < 0.05,
  ,
  drop = FALSE
]


if (
  nrow(
    significant_positive
  ) > 0
) {

  strongest_significant_positive <-
    significant_positive[
      which.max(
        significant_positive$rho
      ),
      ,
      drop = FALSE
    ]


  representative_results <- rbind(
    representative_results,
    strongest_significant_positive
  )

}


if (
  nrow(
    representative_results
  ) > 0
) {

  representative_results <-
    representative_results[
      !duplicated(
        paste(
          representative_results$lipid_pathway,
          representative_results$ddr_pathway,
          sep = "___"
        )
      ),
      ,
      drop = FALSE
    ]

}


write.csv(

  representative_results,

  file.path(
    tables_dir,
    "10_2_3_representative_correlation_pairs.csv"
  ),

  row.names = FALSE

)


# ============================================================
# 10.20 — REPRESENTATIVE SCATTER PLOTS
# ============================================================

if (
  nrow(
    representative_results
  ) > 0
) {

  for (
    i in seq_len(
      nrow(
        representative_results
      )
    )
  ) {

    lipid_name <-
      representative_results$lipid_pathway[i]


    ddr_name <-
      representative_results$ddr_pathway[i]


    rho_value <-
      representative_results$rho[i]


    fdr_value <-
      representative_results$FDR[i]


    scatter_data <- data.frame(

      lipid_score =
        as.numeric(
          tcga_scores[
            lipid_name,
            ]
        ),

      ddr_score =
        as.numeric(
          tcga_scores[
            ddr_name,
            ]
        )

    )


    scatter_data <- scatter_data[
      complete.cases(
        scatter_data
      ),
      ,
      drop = FALSE
    ]


    scatter_plot <- ggplot(

      scatter_data,

      aes(
        x = lipid_score,
        y = ddr_score
      )

    ) +

      geom_point(
        alpha = 0.65,
        size = 2
      ) +

      geom_smooth(
        method = "lm",
        se = TRUE
      ) +

      labs(

        title =
          paste(
            lipid_name,
            "vs",
            ddr_name
          ),

        subtitle =
          paste0(
            "Spearman rho = ",
            sprintf(
              "%.3f",
              rho_value
            ),
            "; BH-FDR = ",
            ifelse(
              fdr_value < 0.001,
              "<0.001",
              sprintf(
                "%.3f",
                fdr_value
              )
            )
          ),

        x =
          paste(
            lipid_name,
            "GSVA score"
          ),

        y =
          paste(
            ddr_name,
            "GSVA score"
          )

      ) +

      theme_minimal(
        base_size = 12
      ) +

      theme(
        plot.title =
          element_text(
            face = "bold"
          )
      )


    safe_name <- paste0(

      gsub(
        "[^A-Za-z0-9]+",
        "_",
        lipid_name
      ),

      "_vs_",

      gsub(
        "[^A-Za-z0-9]+",
        "_",
        ddr_name
      )

    )


    ggsave(

      file.path(
        figures_dir,
        paste0(
          "10_2_3_scatter_",
          safe_name,
          ".png"
        )
      ),

      scatter_plot,

      width = 8,

      height = 7,

      dpi = 300

    )


    ggsave(

      file.path(
        figures_dir,
        paste0(
          "10_2_3_scatter_",
          safe_name,
          ".pdf"
        )
      ),

      scatter_plot,

      width = 8,

      height = 7

    )

  }

}


# ============================================================
# 10.21 — LIPID-HIGH / LIPID-LOW GROUPING
# ============================================================

group_defining_pathway <-
  "HALLMARK_FATTY_ACID_METABOLISM"


group_defining_score <- as.numeric(

  tcga_scores[
    group_defining_pathway,
    ]

)


group_cutoff <- median(

  group_defining_score,

  na.rm = TRUE

)


if (
  !is.finite(
    group_cutoff
  )
) {

  stop(
    "Could not calculate the lipid-group median cutoff."
  )

}


lipid_group <- ifelse(

  group_defining_score >= group_cutoff,

  "Lipid-high",

  "Lipid-low"

)


lipid_group <- factor(

  lipid_group,

  levels = c(
    "Lipid-low",
    "Lipid-high"
  )

)


group_assignment <- data.frame(

  sample_id =
    colnames(
      tcga_scores
    ),

  lipid_score =
    group_defining_score,

  lipid_group =
    lipid_group,

  stringsAsFactors = FALSE

)


write.csv(

  group_assignment,

  file.path(
    tables_dir,
    "10_2_4_lipid_high_low_group_assignment.csv"
  ),

  row.names = FALSE

)


if (
  sum(
    lipid_group == "Lipid-low"
  ) == 0 ||
  sum(
    lipid_group == "Lipid-high"
  ) == 0
) {

  stop(
    "Lipid-high/lipid-low grouping failed: one group contains zero samples."
  )

}


# ============================================================
# 10.22 — CORRECT RANK-BISERIAL EFFECT SIZE
# ============================================================

calculate_rank_biserial <- function(
  high_values,
  low_values
) {

  high_values <- high_values[
    is.finite(
      high_values
    )
  ]


  low_values <- low_values[
    is.finite(
      low_values
    )
  ]


  n_high <- length(
    high_values
  )


  n_low <- length(
    low_values
  )


  if (
    n_high == 0 ||
    n_low == 0
  ) {

    return(
      list(
        U_high = NA_real_,
        rank_biserial = NA_real_,
        n_high = n_high,
        n_low = n_low
      )
    )

  }


  combined_values <- c(
    high_values,
    low_values
  )


  ranks <- rank(
    combined_values,
    ties.method = "average"
  )


  rank_sum_high <- sum(
    ranks[
      seq_len(
        n_high
      )
    ]
  )


  U_high <- rank_sum_high -

    (
      n_high *
      (
        n_high + 1
      ) /
      2
    )


  rank_biserial <- (

    2 *
    U_high /
    (
      n_high *
      n_low
    )

  ) - 1


  rank_biserial <- max(
    -1,
    min(
      1,
      rank_biserial
    )
  )


  list(

    U_high =
      U_high,

    rank_biserial =
      rank_biserial,

    n_high =
      n_high,

    n_low =
      n_low

  )

}


# ============================================================
# 10.23 — DIFFERENTIAL PATHWAY ANALYSIS
# ============================================================

all_differential_pathways <- unique(

  c(
    lipid_pathways,
    ddr_pathways,
    checkpoint_pathways
  )

)


differential_results <- data.frame(

  pathway =
    character(0),

  biological_group =
    character(0),

  group_defining =
    logical(0),

  independent_interpretation =
    logical(0),

  n_lipid_low =
    integer(0),

  n_lipid_high =
    integer(0),

  median_lipid_low =
    numeric(0),

  median_lipid_high =
    numeric(0),

  median_difference_high_minus_low =
    numeric(0),

  mean_lipid_low =
    numeric(0),

  mean_lipid_high =
    numeric(0),

  mean_difference_high_minus_low =
    numeric(0),

  W =
    numeric(0),

  p_value =
    numeric(0),

  rank_biserial =
    numeric(0),

  effect_direction =
    character(0),

  stringsAsFactors = FALSE

)


for (
  pathway_name in all_differential_pathways
) {

  pathway_values <- as.numeric(

    tcga_scores[
      pathway_name,
      ]

  )


  valid <- is.finite(
    pathway_values
  )


  high_values <- pathway_values[
    valid &
    lipid_group == "Lipid-high"
  ]


  low_values <- pathway_values[
    valid &
    lipid_group == "Lipid-low"
  ]


  n_high <- length(
    high_values
  )


  n_low <- length(
    low_values
  )


  median_high <- if (
    n_high > 0
  ) {
    median(
      high_values
    )
  } else {
    NA_real_
  }


  median_low <- if (
    n_low > 0
  ) {
    median(
      low_values
    )
  } else {
    NA_real_
  }


  mean_high <- if (
    n_high > 0
  ) {
    mean(
      high_values
    )
  } else {
    NA_real_
  }


  mean_low <- if (
    n_low > 0
  ) {
    mean(
      low_values
    )
  } else {
    NA_real_
  }


  if (
    n_high > 0 &&
    n_low > 0
  ) {

    wilcox_result <- suppressWarnings(

      wilcox.test(

        high_values,

        low_values,

        exact = FALSE

      )

    )


    W_value <- as.numeric(
      wilcox_result$statistic
    )


    p_value <- wilcox_result$p.value


    effect_result <- calculate_rank_biserial(

      high_values,

      low_values

    )


    rank_biserial <-
      effect_result$rank_biserial

  } else {

    W_value <- NA_real_

    p_value <- NA_real_

    rank_biserial <- NA_real_

  }


  effect_direction <- ifelse(

    is.na(
      rank_biserial
    ),

    "Not estimable",

    ifelse(

      rank_biserial > 0,

      "Higher in lipid-high",

      ifelse(

        rank_biserial < 0,

        "Lower in lipid-high",

        "No directional difference"

      )

    )

  )


  biological_group <- if (

    pathway_name %in%
    lipid_pathways

  ) {

    "Lipid metabolism"

  } else if (

    pathway_name %in%
    ddr_pathways

  ) {

    "DDR / repair"

  } else {

    "Checkpoint / proliferation"

  }


  new_result <- data.frame(

    pathway =
      pathway_name,

    biological_group =
      biological_group,

    group_defining =
      pathway_name ==
      group_defining_pathway,

    independent_interpretation =
      pathway_name !=
      group_defining_pathway,

    n_lipid_low =
      n_low,

    n_lipid_high =
      n_high,

    median_lipid_low =
      median_low,

    median_lipid_high =
      median_high,

    median_difference_high_minus_low =
      median_high -
      median_low,

    mean_lipid_low =
      mean_low,

    mean_lipid_high =
      mean_high,

    mean_difference_high_minus_low =
      mean_high -
      mean_low,

    W =
      W_value,

    p_value =
      p_value,

    rank_biserial =
      rank_biserial,

    effect_direction =
      effect_direction,

    stringsAsFactors = FALSE

  )


  differential_results <- rbind(

    differential_results,

    new_result

  )

}


# ============================================================
# 10.24 — BH FDR
# ============================================================

differential_results$FDR <- p.adjust(

  differential_results$p_value,

  method = "BH"

)


# ============================================================
# 10.25 — RECALCULATE EFFECT DIRECTION
# ============================================================

differential_results$effect_direction <- ifelse(

  is.na(
    differential_results$rank_biserial
  ),

  "Not estimable",

  ifelse(

    differential_results$rank_biserial > 0,

    "Higher in lipid-high",

    ifelse(

      differential_results$rank_biserial < 0,

      "Lower in lipid-high",

      "No directional difference"

    )

  )

)


differential_results$significance <- ifelse(

  is.na(
    differential_results$FDR
  ),

  "NA",

  ifelse(

    differential_results$FDR < 0.001,

    "***",

    ifelse(

      differential_results$FDR < 0.01,

      "**",

      ifelse(

        differential_results$FDR < 0.05,

        "*",

        "ns"

      )

    )

  )

)


differential_results$significant_FDR_0_05 <-

  differential_results$FDR < 0.05


# ============================================================
# 10.26 — EXPLICIT GROUP-DEFINING FLAG
# ============================================================

differential_results$interpretation_note <- ifelse(

  differential_results$group_defining,

  "GROUP-DEFINING PATHWAY — NOT AN INDEPENDENT DIFFERENTIAL FINDING",

  "Eligible for independent interpretation"

)


# ============================================================
# 10.27 — EFFECT-SIZE VALIDATION
# ============================================================

finite_effects <- differential_results[
  is.finite(
    differential_results$rank_biserial
  ),
  "rank_biserial"
]


if (
  length(
    finite_effects
  ) == 0
) {

  stop(
    "No finite rank-biserial effect sizes were calculated."
  )

}


if (
  any(
    finite_effects < -1 |
    finite_effects > 1
  )
) {

  stop(
    "FINAL EFFECT-SIZE QC FAILED: rank-biserial effect size outside [-1,+1]."
  )

}


effect_size_validation <- data.frame(

  minimum_rank_biserial =
    min(
      finite_effects
    ),

  maximum_rank_biserial =
    max(
      finite_effects
    ),

  all_between_minus1_and_plus1 =
    all(
      finite_effects >= -1 &
      finite_effects <= 1
    ),

  number_of_finite_effect_sizes =
    length(
      finite_effects
    ),

  stringsAsFactors = FALSE

)


write.csv(

  effect_size_validation,

  file.path(
    tables_dir,
    "10_2_4_rank_biserial_effect_size_validation.csv"
  ),

  row.names = FALSE

)


# ============================================================
# 10.28 — SAVE DIFFERENTIAL RESULTS
# ============================================================

differential_results <- differential_results[
  order(
    differential_results$FDR,
    -abs(
      differential_results$rank_biserial
    )
  ),
  ,
  drop = FALSE
]


write.csv(

  differential_results,

  file.path(
    tables_dir,
    "10_2_4_complete_differential_pathway_activity.csv"
  ),

  row.names = FALSE

)


independent_differential_results <-
  differential_results[
    differential_results$independent_interpretation,
    ,
    drop = FALSE
  ]


write.csv(

  independent_differential_results,

  file.path(
    tables_dir,
    "10_2_4_independent_differential_pathway_results.csv"
  ),

  row.names = FALSE

)


# ============================================================
# 10.29 — DIFFERENTIAL LOLLIPOP
#
# The group-defining fatty-acid pathway is EXCLUDED.
# ============================================================

lollipop_data <-
  independent_differential_results[
    is.finite(
      independent_differential_results$rank_biserial
    ),
    ,
    drop = FALSE
  ]


if (
  nrow(
    lollipop_data
  ) > 0
) {

  lollipop_data <- lollipop_data[
    order(
      lollipop_data$rank_biserial
    ),
    ,
    drop = FALSE
  ]


  lollipop_data$pathway <- factor(

    lollipop_data$pathway,

    levels =
      lollipop_data$pathway

  )


  lollipop_plot <- ggplot(

    lollipop_data,

    aes(
      x = rank_biserial,
      y = pathway
    )

  ) +

    geom_segment(

      aes(
        x = 0,
        xend = rank_biserial,
        y = pathway,
        yend = pathway
      ),

      linewidth = 0.8

    ) +

    geom_point(
      size = 3.5
    ) +

    geom_vline(

      xintercept = 0,

      linetype = "dashed",

      linewidth = 0.6

    ) +

    labs(

      title =
        "Differential Pathway Activity by Lipid Status",

      subtitle =
        paste0(
          "Rank-biserial effect size; ",
          group_defining_pathway,
          " excluded because it defines the groups"
        ),

      x =
        "Rank-biserial effect size\n(positive = higher in lipid-high)",

      y =
        "Pathway"

    ) +

    theme_minimal(
      base_size = 12
    ) +

    theme(

      plot.title =
        element_text(
          face = "bold"
        ),

      axis.text.y =
        element_text(
          size = 9
        )

    )


  ggsave(

    file.path(
      figures_dir,
      "10_2_4_differential_pathway_lollipop.png"
    ),

    lollipop_plot,

    width = 10,

    height = 8,

    dpi = 300

  )


  ggsave(

    file.path(
      figures_dir,
      "10_2_4_differential_pathway_lollipop.pdf"
    ),

    lollipop_plot,

    width = 10,

    height = 8

  )

}


# ============================================================
# 10.30 — KEY PATHWAY VIOLIN / BOXPLOTS
#
# Excludes the group-defining pathway from independent plots.
# ============================================================

key_pathways <- unique(

  c(

    setdiff(
      lipid_pathways,
      group_defining_pathway
    ),

    ddr_pathways

  )

)


key_plot_data <- data.frame(

  sample_id =
    rep(
      colnames(
        tcga_scores
      ),
      length(
        key_pathways
      )
    ),

  pathway =
    rep(
      key_pathways,
      each =
        ncol(
          tcga_scores
        )
    ),

  GSVA_score =
    as.vector(
      t(
        tcga_scores[
          key_pathways,
          ,
          drop = FALSE
        ]
      )
    ),

  lipid_group =
    rep(
      as.character(
        lipid_group
      ),
      length(
        key_pathways
      )
    ),

  stringsAsFactors = FALSE

)


key_plot_data <- key_plot_data[
  is.finite(
    key_plot_data$GSVA_score
  ),
  ,
  drop = FALSE
]


write.csv(

  key_plot_data,

  file.path(
    tables_dir,
    "10_2_4_key_pathway_group_comparison_long_data.csv"
  ),

  row.names = FALSE

)


if (
  nrow(
    key_plot_data
  ) > 0
) {

  violin_plot <- ggplot(

    key_plot_data,

    aes(
      x = lipid_group,
      y = GSVA_score
    )

  ) +

    geom_violin(
      trim = FALSE,
      alpha = 0.4
    ) +

    geom_boxplot(
      width = 0.18,
      outlier.shape = NA
    ) +

    geom_jitter(
      width = 0.08,
      alpha = 0.35,
      size = 1.2
    ) +

    facet_wrap(
      ~ pathway,
      scales = "free_y",
      ncol = 3
    ) +

    labs(

      title =
        "GSVA Pathway Activity by Lipid Status",

      subtitle =
        "Group-defining fatty-acid pathway excluded from independent interpretation",

      x =
        "Lipid status",

      y =
        "GSVA score"

    ) +

    theme_minimal(
      base_size = 11
    ) +

    theme(

      plot.title =
        element_text(
          face = "bold"
        ),

      strip.text =
        element_text(
          size = 8
        )

    )


  ggsave(

    file.path(
      figures_dir,
      "10_2_4_key_pathway_violin_boxplots.png"
    ),

    violin_plot,

    width = 13,

    height = 10,

    dpi = 300

  )


  ggsave(

    file.path(
      figures_dir,
      "10_2_4_key_pathway_violin_boxplots.pdf"
    ),

    violin_plot,

    width = 13,

    height = 10

  )

}


# ============================================================
# 10.31 — GROUP-DEFINING PATHWAY TRANSPARENCY PLOT
# ============================================================

group_defining_data <- data.frame(

  lipid_group =
    lipid_group,

  lipid_score =
    group_defining_score

)


group_defining_plot <- ggplot(

  group_defining_data,

  aes(
    x = lipid_group,
    y = lipid_score
  )

) +

  geom_violin(
    trim = FALSE,
    alpha = 0.4
  ) +

  geom_boxplot(
    width = 0.2,
    outlier.shape = NA
  ) +

  geom_jitter(
    width = 0.08,
    alpha = 0.4,
    size = 1.5
  ) +

  labs(

    title =
      "Group-Defining Fatty-Acid Metabolism Score",

    subtitle =
      "Shown for transparency; NOT an independent differential finding",

    x =
      "Lipid status",

    y =
      "HALLMARK fatty-acid metabolism GSVA score"

  ) +

  theme_minimal(
    base_size = 12
  )


ggsave(

  file.path(
    figures_dir,
    "10_2_4_group_defining_fatty_acid_pathway.png"
  ),

  group_defining_plot,

  width = 7,

  height = 6,

  dpi = 300

)


# ============================================================
# 10.32 — DIFFERENTIAL SUMMARY
# ============================================================

differential_summary <- data.frame(

  total_pathways_tested =
    nrow(
      differential_results
    ),

  independent_pathways =
    sum(
      differential_results$independent_interpretation
    ),

  group_defining_pathway =
    group_defining_pathway,

  independent_FDR_lt_0_05 =
    sum(
      differential_results$independent_interpretation &
      differential_results$FDR < 0.05,
      na.rm = TRUE
    ),

  independent_FDR_lt_0_01 =
    sum(
      differential_results$independent_interpretation &
      differential_results$FDR < 0.01,
      na.rm = TRUE
    ),

  independent_FDR_lt_0_001 =
    sum(
      differential_results$independent_interpretation &
      differential_results$FDR < 0.001,
      na.rm = TRUE
    ),

  stringsAsFactors = FALSE

)


write.csv(

  differential_summary,

  file.path(
    tables_dir,
    "10_2_4_differential_analysis_summary.csv"
  ),

  row.names = FALSE

)


# ============================================================
# 10.33 — SAVE COMPLETE STEP 10 RDS
# ============================================================

step10_results <- list(

  step9_tcga_source =
    tcga_step9_file,

  step9_gse_source =
    gse_step9_file,

  tcga_scores =
    tcga_scores,

  gse_scores =
    gse_scores,

  expected_pathways =
    expected_pathways,

  lipid_pathways =
    lipid_pathways,

  ddr_pathways =
    ddr_pathways,

  checkpoint_pathways =
    checkpoint_pathways,

  pathway_audit =
    pathway_audit,

  pathway_qc =
    pathway_qc,

  group_defining_pathway =
    group_defining_pathway,

  group_cutoff =
    group_cutoff,

  group_assignment =
    group_assignment,

  correlation_results =
    correlation_results,

  rho_matrix =
    rho_matrix,

  p_matrix =
    p_matrix,

  fdr_matrix =
    fdr_matrix,

  ci_lower_matrix =
    ci_lower_matrix,

  ci_upper_matrix =
    ci_upper_matrix,

  significance_matrix =
    significance_matrix,

  representative_results =
    representative_results,

  differential_results =
    differential_results,

  independent_differential_results =
    independent_differential_results,

  effect_size_validation =
    effect_size_validation,

  differential_summary =
    differential_summary

)


saveRDS(

  step10_results,

  file.path(
    processed_dir,
    "10_pathway_correlation_and_differential_activity.rds"
  )

)


# ============================================================
# 10.34 — OUTPUT INVENTORY
# ============================================================

output_inventory <- data.frame(

  category = c(

    "Input",
    "Input",
    "QC",
    "Correlation",
    "Correlation",
    "Correlation",
    "Correlation",
    "Correlation",
    "Correlation",
    "Correlation",
    "Correlation",
    "Correlation",
    "Differential",
    "Differential",
    "Differential",
    "Differential",
    "Differential",
    "Differential",
    "Final RDS"

  ),

  file = c(

    tcga_step9_file,

    gse_step9_file,

    "10_1_pathway_score_QC.csv",

    "10_2_3_complete_lipid_DDR_spearman_results.csv",

    "10_2_3_spearman_rho_matrix.csv",

    "10_2_3_spearman_p_value_matrix.csv",

    "10_2_3_spearman_BH_FDR_matrix.csv",

    "10_2_3_spearman_95CI_lower_matrix.csv",

    "10_2_3_spearman_95CI_upper_matrix.csv",

    "10_2_3_FDR_significance_annotation_matrix.csv",

    "10_2_3_representative_correlation_pairs.csv",

    "10_2_3_lipid_DDR_correlation_heatmap.png",

    "10_2_4_complete_differential_pathway_activity.csv",

    "10_2_4_independent_differential_pathway_results.csv",

    "10_2_4_rank_biserial_effect_size_validation.csv",

    "10_2_4_lipid_high_low_group_assignment.csv",

    "10_2_4_key_pathway_violin_boxplots.png",

    "10_2_4_differential_pathway_lollipop.png",

    "10_pathway_correlation_and_differential_activity.rds"

  ),

  stringsAsFactors = FALSE

)


write.csv(

  output_inventory,

  file.path(
    tables_dir,
    "10_complete_output_inventory.csv"
  ),

  row.names = FALSE

)


# ============================================================
# 10.35 — FINAL QC
# ============================================================

# ------------------------------------------------------------
# 17 pathways in TCGA
# ------------------------------------------------------------

if (
  length(
    intersect(
      expected_pathways,
      rownames(
        tcga_scores
      )
    )
  ) != 17
) {

  stop(
    "FINAL QC FAILED: TCGA does not contain all 17 pathways."
  )

}


# ------------------------------------------------------------
# 17 pathways in GSE76427
# ------------------------------------------------------------

if (
  length(
    intersect(
      expected_pathways,
      rownames(
        gse_scores
      )
    )
  ) != 17
) {

  stop(
    "FINAL QC FAILED: GSE76427 does not contain all 17 pathways."
  )

}


# ------------------------------------------------------------
# 48 lipid × DDR correlations
# ------------------------------------------------------------

if (
  nrow(
    correlation_results
  ) != 48
) {

  stop(
    "FINAL QC FAILED: expected exactly 48 lipid × DDR correlations."
  )

}


# ------------------------------------------------------------
# Effect-size range
# ------------------------------------------------------------

if (
  any(
    finite_effects < -1 |
    finite_effects > 1
  )
) {

  stop(
    "FINAL QC FAILED: rank-biserial effect outside [-1,+1]."
  )

}


# ------------------------------------------------------------
# Group-defining pathway exclusion
# ------------------------------------------------------------

if (
  any(
    differential_results$group_defining &
    differential_results$independent_interpretation
  )
) {

  stop(
    "FINAL QC FAILED: group-defining pathway incorrectly marked independent."
  )

}


# ------------------------------------------------------------
# Group assignment
# ------------------------------------------------------------

if (
  any(
    is.na(
      lipid_group
    )
  )
) {

  stop(
    "FINAL QC FAILED: one or more TCGA samples have no lipid group."
  )

}


# ============================================================
# 10.36 — FINAL CONSOLE REPORT
# ============================================================

cat("\n")
cat("============================================================\n")
cat("STEP 10 COMPLETED SUCCESSFULLY\n")
cat("============================================================\n\n")


cat(
  "TCGA Step 9 source:\n",
  tcga_step9_file,
  "\n\n"
)


cat(
  "GSE76427 Step 9 source:\n",
  gse_step9_file,
  "\n\n"
)


cat(
  "TCGA samples: ",
  ncol(
    tcga_scores
  ),
  "\n",
  sep = ""
)


cat(
  "GSE76427 samples: ",
  ncol(
    gse_scores
  ),
  "\n",
  sep = ""
)


cat(
  "Pathways: 17\n"
)


cat(
  "Lipid pathways: ",
  length(
    lipid_pathways
  ),
  "\n",
  sep = ""
)


cat(
  "DDR pathways: ",
  length(
    ddr_pathways
  ),
  "\n",
  sep = ""
)


cat(
  "Lipid × DDR correlations: ",
  nrow(
    correlation_results
  ),
  "\n",
  sep = ""
)


cat(
  "BH-FDR < 0.05 correlations: ",
  sum(
    correlation_results$FDR < 0.05,
    na.rm = TRUE
  ),
  "\n",
  sep = ""
)


cat(
  "Independent differential pathways: ",
  sum(
    differential_results$independent_interpretation
  ),
  "\n",
  sep = ""
)


cat(
  "Independent differential pathways BH-FDR < 0.05: ",
  sum(
    differential_results$independent_interpretation &
    differential_results$FDR < 0.05,
    na.rm = TRUE
  ),
  "\n",
  sep = ""
)


cat(
  "Minimum rank-biserial effect: ",
  sprintf(
    "%.4f",
    min(
      finite_effects
    )
  ),
  "\n",
  sep = ""
)


cat(
  "Maximum rank-biserial effect: ",
  sprintf(
    "%.4f",
    max(
      finite_effects
    )
  ),
  "\n",
  sep = ""
)


cat(
  "All rank-biserial effects within [-1,+1]: TRUE\n"
)


cat(
  "Group-defining fatty-acid pathway excluded from independent interpretation: TRUE\n"
)


cat("\n")
cat("Results directory:\n")
cat(
  results_dir,
  "\n\n"
)


cat("Final Step 10 RDS:\n")
cat(
  file.path(
    processed_dir,
    "10_pathway_correlation_and_differential_activity.rds"
  ),
  "\n\n"
)


cat("============================================================\n")
cat("STEP 10 FINISHED\n")
cat("============================================================\n")
