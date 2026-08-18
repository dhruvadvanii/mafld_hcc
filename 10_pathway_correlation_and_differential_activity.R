# ============================================================
# STEP 10 — PATHWAY CORRELATION AND DIFFERENTIAL PATHWAY ACTIVITY
# ============================================================
#
# MAFLD-HCC in silico analysis
#
# ANALYSES:
#
# 10.1  Load Step 9 GSVA outputs
# 10.2  Identify lipid and DDR pathways
# 10.3  Spearman correlation: lipid × DDR pathways
# 10.4  Benjamini-Hochberg FDR correction
# 10.5  Correlation heatmap
# 10.6  Representative correlation scatterplots
# 10.7  Define lipid-high / lipid-low tumors
# 10.8  Wilcoxon rank-sum differential pathway activity
# 10.9  Effect size and BH-FDR
# 10.10 Save all results and figures
# 10.11 Final QC
#
# IMPORTANT:
#   - Uses ONLY TCGA-LIHC GSVA scores generated in Step 9.
#   - Uses Spearman correlation.
#   - BH correction is applied across all lipid × DDR tests.
#   - Lipid-high / lipid-low is defined using the
#     HALLMARK_FATTY_ACID_METABOLISM GSVA score.
#   - Median split:
#         LOW  = score <= median
#         HIGH = score > median
#   - Wilcoxon rank-sum tests are two-sided.
#   - Effect size = rank-biserial correlation.
#   - No pathway is silently substituted.
#
# ============================================================


# ============================================================
# 10.0 — INITIAL SETUP
# ============================================================

options(stringsAsFactors = FALSE)

cat("\n============================================================\n")
cat("STEP 10 — PATHWAY CORRELATION + DIFFERENTIAL ACTIVITY\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# Required packages
# ------------------------------------------------------------

required_packages <- c(
  "ggplot2"
)


for (pkg in required_packages) {

  if (!requireNamespace(pkg, quietly = TRUE)) {

    stop(
      paste0(
        "Required package '",
        pkg,
        "' is not installed.\n",
        "Install it before running Step 10."
      )
    )

  }

}


# ============================================================
# 10.1 — DIRECTORIES
# ============================================================

data_dir <- "data"

processed_dir <- file.path(
  data_dir,
  "processed"
)

results_dir <- "results"

figures_dir <- file.path(
  results_dir,
  "figures",
  "10_pathway_analysis"
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

dir.create(
  figures_dir,
  showWarnings = FALSE,
  recursive = TRUE
)


# ============================================================
# 10.2 — LOAD STEP 9 OUTPUTS
# ============================================================

step9_file <- file.path(
  processed_dir,
  "step9_GSVA_results.rds"
)


tcga_scores_file <- file.path(
  processed_dir,
  "TCGA-LIHC",
  "TCGA_LIHC_GSVA_pathway_scores.rds"
)


if (
  file.exists(step9_file)
) {

  step9_results <- readRDS(
    step9_file
  )

} else {

  step9_results <- NULL

}


if (
  file.exists(tcga_scores_file)
) {

  tcga_pathway_scores <- readRDS(
    tcga_scores_file
  )

} else {

  if (
    is.null(step9_results)
    ||
    is.null(
      step9_results$tcga_pathway_scores
    )
  ) {

    stop(
      paste0(
        "Could not find TCGA-LIHC GSVA scores.\n\n",
        "Expected either:\n",
        tcga_scores_file,
        "\nOR:\n",
        step9_file
      )
    )

  }


  tcga_pathway_scores <-
    step9_results$tcga_pathway_scores

}


# ------------------------------------------------------------
# Convert to matrix if necessary
# ------------------------------------------------------------

if (
  !is.matrix(
    tcga_pathway_scores
  )
) {

  tcga_pathway_scores <- as.matrix(
    tcga_pathway_scores
  )

}


storage.mode(
  tcga_pathway_scores
) <- "numeric"


# ------------------------------------------------------------
# Basic matrix validation
# ------------------------------------------------------------

if (
  is.null(
    rownames(
      tcga_pathway_scores
    )
  )
) {

  stop(
    "TCGA GSVA matrix has no pathway names in rownames()."
  )

}


if (
  is.null(
    colnames(
      tcga_pathway_scores
    )
  )
) {

  stop(
    "TCGA GSVA matrix has no sample identifiers in colnames()."
  )

}


if (
  nrow(
    tcga_pathway_scores
  ) == 0
  ||
  ncol(
    tcga_pathway_scores
  ) == 0
) {

  stop(
    "TCGA GSVA matrix is empty."
  )

}


# ============================================================
# 10.3 — STANDARDIZE PATHWAY NAMES
# ============================================================

rownames(
  tcga_pathway_scores
) <- toupper(
  trimws(
    rownames(
      tcga_pathway_scores
    )
  )
)


# Remove accidental duplicated pathway rows
# rather than silently using duplicated names.

if (
  anyDuplicated(
    rownames(
      tcga_pathway_scores
    )
  ) > 0
) {

  duplicated_pathways <- unique(
    rownames(
      tcga_pathway_scores
    )[
      duplicated(
        rownames(
          tcga_pathway_scores
        )
      )
    ]
  )


  stop(
    paste0(
      "Duplicate pathway names detected in TCGA GSVA matrix:\n",
      paste(
        duplicated_pathways,
        collapse = "\n"
      )
    )
  )

}


cat(
  "Loaded TCGA-LIHC GSVA matrix:\n"
)

cat(
  "  Pathways:",
  nrow(
    tcga_pathway_scores
  ),
  "\n"
)

cat(
  "  Samples:",
  ncol(
    tcga_pathway_scores
  ),
  "\n\n"
)


# ============================================================
# 10.4 — DEFINE BIOLOGICAL PATHWAY GROUPS
# ============================================================

# ------------------------------------------------------------
# Lipid metabolism pathways
# ------------------------------------------------------------

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


# ------------------------------------------------------------
# DDR / DNA repair pathways
# ------------------------------------------------------------

ddr_pathways <- c(

  "HOMOLOGOUS_RECOMBINATION",

  "MISMATCH_REPAIR",

  "BASE_EXCISION_REPAIR",

  "DOUBLE_STRAND_BREAK_REPAIR",

  "NON_HOMOLOGOUS_END_JOINING",

  "HALLMARK_DNA_REPAIR"

)


# ------------------------------------------------------------
# Checkpoint / proliferation pathways
# ------------------------------------------------------------

checkpoint_pathways <- c(

  "P53_PATHWAY",

  "P53_DEPENDENT_G1_DNA_DAMAGE_RESPONSE",

  "G2M_CHECKPOINT"

)


# ============================================================
# 10.5 — CHECK WHICH PATHWAYS SURVIVED STEP 9 QC
# ============================================================

available_pathways <- rownames(
  tcga_pathway_scores
)


lipid_present <- intersect(
  lipid_pathways,
  available_pathways
)


ddr_present <- intersect(
  ddr_pathways,
  available_pathways
)


checkpoint_present <- intersect(
  checkpoint_pathways,
  available_pathways
)


cat(
  "------------------------------------------------------------\n"
)

cat(
  "PATHWAY GROUP AUDIT\n"
)

cat(
  "------------------------------------------------------------\n\n"
)


cat(
  "Lipid pathways available:",
  length(
    lipid_present
  ),
  "/",
  length(
    lipid_pathways
  ),
  "\n"
)


cat(
  "DDR pathways available:",
  length(
    ddr_present
  ),
  "/",
  length(
    ddr_pathways
  ),
  "\n"
)


cat(
  "Checkpoint pathways available:",
  length(
    checkpoint_present
  ),
  "/",
  length(
    checkpoint_pathways
  ),
  "\n\n"
)


# ------------------------------------------------------------
# Report excluded pathways
# ------------------------------------------------------------

lipid_missing <- setdiff(
  lipid_pathways,
  available_pathways
)


ddr_missing <- setdiff(
  ddr_pathways,
  available_pathways
)


checkpoint_missing <- setdiff(
  checkpoint_pathways,
  available_pathways
)


if (
  length(lipid_missing) > 0
) {

  cat(
    "Lipid pathways not retained in Step 9:\n",
    paste(
      lipid_missing,
      collapse = "\n"
    ),
    "\n\n"
  )

}


if (
  length(ddr_missing) > 0
) {

  cat(
    "DDR pathways not retained in Step 9:\n",
    paste(
      ddr_missing,
      collapse = "\n"
    ),
    "\n\n"
  )

}


if (
  length(checkpoint_missing) > 0
) {

  cat(
    "Checkpoint pathways not retained in Step 9:\n",
    paste(
      checkpoint_missing,
      collapse = "\n"
    ),
    "\n\n"
  )

}


# ------------------------------------------------------------
# Minimum requirement for correlation analysis
# ------------------------------------------------------------

if (
  length(lipid_present) == 0
) {

  stop(
    "No lipid pathways survived Step 9 and are available for correlation analysis."
  )

}


if (
  length(ddr_present) == 0
) {

  stop(
    "No DDR pathways survived Step 9 and are available for correlation analysis."
  )

}


# ============================================================
# 10.6 — EXTRACT LIPID AND DDR SCORE MATRICES
# ============================================================

lipid_scores <- tcga_pathway_scores[
  lipid_present,
  ,
  drop = FALSE
]


ddr_scores <- tcga_pathway_scores[
  ddr_present,
  ,
  drop = FALSE
]


# ------------------------------------------------------------
# Confirm identical samples
# ------------------------------------------------------------

if (
  !identical(
    colnames(lipid_scores),
    colnames(ddr_scores)
  )
) {

  common_samples <- intersect(
    colnames(lipid_scores),
    colnames(ddr_scores)
  )


  if (
    length(common_samples) == 0
  ) {

    stop(
      "No common TCGA samples between lipid and DDR pathway matrices."
    )

  }


  lipid_scores <- lipid_scores[
    ,
    common_samples,
    drop = FALSE
  ]


  ddr_scores <- ddr_scores[
    ,
    common_samples,
    drop = FALSE
  ]

}


cat(
  "Samples used for lipid × DDR correlation:",
  ncol(
    lipid_scores
  ),
  "\n\n"
)


# ============================================================
# 10.7 — SPEARMAN CORRELATION FUNCTION
# ============================================================

calculate_spearman <- function(
  x,
  y
) {

  complete <- complete.cases(
    x,
    y
  )


  n_complete <- sum(
    complete
  )


  if (
    n_complete < 3
  ) {

    return(
      data.frame(
        rho = NA_real_,
        p_value = NA_real_,
        n = n_complete,
        stringsAsFactors = FALSE
      )
    )

  }


  x_complete <- x[
    complete
  ]


  y_complete <- y[
    complete
  ]


  if (
    sd(
      x_complete
    ) == 0
    ||
    sd(
      y_complete
    ) == 0
  ) {

    return(
      data.frame(
        rho = NA_real_,
        p_value = NA_real_,
        n = n_complete,
        stringsAsFactors = FALSE
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


  data.frame(
    rho = unname(
      test$estimate
    ),
    p_value = test$p.value,
    n = n_complete,
    stringsAsFactors = FALSE
  )

}


# ============================================================
# 10.8 — CALCULATE ALL LIPID × DDR CORRELATIONS
# ============================================================

correlation_rows <- vector(
  "list",
  length(
    lipid_present
  ) *
  length(
    ddr_present
  )
)


counter <- 1


for (
  lipid_pathway in lipid_present
) {

  for (
    ddr_pathway in ddr_present
  ) {

    result <- calculate_spearman(
      lipid_scores[
        lipid_pathway,
      ],
      ddr_scores[
        ddr_pathway,
      ]
    )


    correlation_rows[[counter]] <-
      data.frame(

        lipid_pathway =
          lipid_pathway,

        ddr_pathway =
          ddr_pathway,

        rho =
          result$rho,

        p_value =
          result$p_value,

        n =
          result$n,

        stringsAsFactors = FALSE

      )


    counter <- counter + 1

  }

}


lipid_ddr_correlations <- do.call(
  rbind,
  correlation_rows
)


rownames(
  lipid_ddr_correlations
) <- NULL


# ============================================================
# 10.9 — BENJAMINI-HOCHBERG FDR
# ============================================================

lipid_ddr_correlations$FDR <- NA_real_


valid_p <- !is.na(
  lipid_ddr_correlations$p_value
)


lipid_ddr_correlations$FDR[
  valid_p
] <- p.adjust(
  lipid_ddr_correlations$p_value[
    valid_p
  ],
  method = "BH"
)


# ------------------------------------------------------------
# Add absolute correlation
# ------------------------------------------------------------

lipid_ddr_correlations$abs_rho <-
  abs(
    lipid_ddr_correlations$rho
  )


# ------------------------------------------------------------
# Significance annotation
# ------------------------------------------------------------

lipid_ddr_correlations$significance <- cut(
  lipid_ddr_correlations$FDR,
  breaks = c(
    -Inf,
    0.001,
    0.01,
    0.05,
    Inf
  ),
  labels = c(
    "***",
    "**",
    "*",
    "ns"
  ),
  right = TRUE
)


# ============================================================
# 10.10 — SAVE CORRELATION RESULTS
# ============================================================

write.csv(
  lipid_ddr_correlations,
  file.path(
    results_dir,
    "10_lipid_DDR_spearman_correlations.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 10.11 — CREATE CORRELATION MATRIX
# ============================================================

correlation_matrix <- matrix(
  NA_real_,
  nrow = length(
    lipid_present
  ),
  ncol = length(
    ddr_present
  ),
  dimnames = list(
    lipid_present,
    ddr_present
  )
)


for (
  i in seq_len(
    nrow(
      lipid_ddr_correlations
    )
  )
) {

  lipid_name <-
    lipid_ddr_correlations$lipid_pathway[i]


  ddr_name <-
    lipid_ddr_correlations$ddr_pathway[i]


  correlation_matrix[
    lipid_name,
    ddr_name
  ] <-
    lipid_ddr_correlations$rho[i]

}


write.csv(
  correlation_matrix,
  file.path(
    results_dir,
    "10_lipid_DDR_spearman_rho_matrix.csv"
  ),
  row.names = TRUE
)


# ============================================================
# 10.12 — CREATE FDR MATRIX
# ============================================================

fdr_matrix <- matrix(
  NA_real_,
  nrow = length(
    lipid_present
  ),
  ncol = length(
    ddr_present
  ),
  dimnames = list(
    lipid_present,
    ddr_present
  )
)


for (
  i in seq_len(
    nrow(
      lipid_ddr_correlations
    )
  )
) {

  lipid_name <-
    lipid_ddr_correlations$lipid_pathway[i]


  ddr_name <-
    lipid_ddr_correlations$ddr_pathway[i]


  fdr_matrix[
    lipid_name,
    ddr_name
  ] <-
    lipid_ddr_correlations$FDR[i]

}


write.csv(
  fdr_matrix,
  file.path(
    results_dir,
    "10_lipid_DDR_BH_FDR_matrix.csv"
  ),
  row.names = TRUE
)


# ============================================================
# 10.13 — CORRELATION HEATMAP
# ============================================================

heatmap_df <- expand.grid(
  lipid_pathway =
    rownames(
      correlation_matrix
    ),
  ddr_pathway =
    colnames(
      correlation_matrix
    ),
  stringsAsFactors = FALSE
)


heatmap_df$rho <- as.vector(
  correlation_matrix
)


heatmap_df$FDR <- as.vector(
  fdr_matrix
)


heatmap_df$label <- ifelse(
  is.na(
    heatmap_df$rho
  ),
  "",
  sprintf(
    "%.2f",
    heatmap_df$rho
  )
)


heatmap_df$lipid_pathway <- factor(
  heatmap_df$lipid_pathway,
  levels =
    rev(
      rownames(
        correlation_matrix
      )
    )
)


heatmap_df$ddr_pathway <- factor(
  heatmap_df$ddr_pathway,
  levels =
    colnames(
      correlation_matrix
    )
)


correlation_heatmap <- ggplot2::ggplot(
  heatmap_df,
  ggplot2::aes(
    x = ddr_pathway,
    y = lipid_pathway,
    fill = rho
  )
) +

  ggplot2::geom_tile(
    colour = "white"
  ) +

  ggplot2::geom_text(
    ggplot2::aes(
      label = label
    ),
    size = 3
  ) +

  ggplot2::scale_fill_gradient2(
    midpoint = 0,
    limits = c(
      -1,
      1
    ),
    name = "Spearman\nrho"
  ) +

  ggplot2::labs(
    title =
      "Lipid Pathway–DDR Pathway Correlations",
    x =
      "DDR / Repair Pathway",
    y =
      "Lipid Metabolism Pathway"
  ) +

  ggplot2::theme_minimal(
    base_size = 11
  ) +

  ggplot2::theme(
    axis.text.x =
      ggplot2::element_text(
        angle = 45,
        hjust = 1
      ),
    panel.grid =
      ggplot2::element_blank()
  )


ggplot2::ggsave(
  filename = file.path(
    figures_dir,
    "10_lipid_DDR_correlation_heatmap.png"
  ),
  plot = correlation_heatmap,
  width = 10,
  height = 7,
  dpi = 300
)


# ============================================================
# 10.14 — SELECT REPRESENTATIVE CORRELATIONS
# ============================================================
#
# Representative scatterplots are selected by:
#
#   1. statistically significant FDR
#   2. strongest absolute Spearman rho
#
# If fewer than 4 significant pairs exist, all available
# significant pairs are used.
#
# Maximum = 4 representative plots.
#
# ============================================================

significant_correlations <-
  lipid_ddr_correlations[
    !is.na(
      lipid_ddr_correlations$FDR
    )
    &
    lipid_ddr_correlations$FDR < 0.05,
    ,
    drop = FALSE
  ]


if (
  nrow(
    significant_correlations
  ) > 0
) {

  significant_correlations <-
    significant_correlations[
      order(
        -significant_correlations$abs_rho,
        significant_correlations$FDR
      ),
      ,
      drop = FALSE
    ]


  representative_correlations <-
    head(
      significant_correlations,
      4
    )

} else {

  representative_correlations <-
    lipid_ddr_correlations[
      order(
        -lipid_ddr_correlations$abs_rho
      ),
      ,
      drop = FALSE
    ]


  representative_correlations <-
    head(
      representative_correlations,
      min(
        4,
        nrow(
          representative_correlations
        )
      )
    )

}


write.csv(
  representative_correlations,
  file.path(
    results_dir,
    "10_representative_lipid_DDR_correlations.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 10.15 — GENERATE REPRESENTATIVE SCATTERPLOTS
# ============================================================

if (
  nrow(
    representative_correlations
  ) > 0
) {

  for (
    i in seq_len(
      nrow(
        representative_correlations
      )
    )
  ) {

    lipid_name <-
      representative_correlations$lipid_pathway[i]


    ddr_name <-
      representative_correlations$ddr_pathway[i]


    x <- as.numeric(
      lipid_scores[
        lipid_name,
        ]
    )


    y <- as.numeric(
      ddr_scores[
        ddr_name,
        ]
    )


    scatter_df <- data.frame(

      lipid_score =
        x,

      ddr_score =
        y

    )


    complete <- complete.cases(
      scatter_df
    )


    scatter_df <- scatter_df[
      complete,
      ,
      drop = FALSE
    ]


    rho_value <-
      representative_correlations$rho[i]


    fdr_value <-
      representative_correlations$FDR[i]


    scatter_plot <- ggplot2::ggplot(
      scatter_df,
      ggplot2::aes(
        x = lipid_score,
        y = ddr_score
      )
    ) +

      ggplot2::geom_point(
        alpha = 0.65
      ) +

      ggplot2::geom_smooth(
        method = "lm",
        se = TRUE
      ) +

      ggplot2::labs(
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
            format.pval(
              fdr_value,
              digits = 3,
              eps = 0.001
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

      ggplot2::theme_minimal(
        base_size = 11
      )


    safe_name <- paste(
      lipid_name,
      "vs",
      ddr_name,
      sep = "_"
    )


    safe_name <- gsub(
      "[^A-Za-z0-9_]+",
      "_",
      safe_name
    )


    ggplot2::ggsave(
      filename = file.path(
        figures_dir,
        paste0(
          "10_scatter_",
          safe_name,
          ".png"
        )
      ),
      plot = scatter_plot,
      width = 7,
      height = 6,
      dpi = 300
    )

  }

}


# ============================================================
# 10.16 — DEFINE PRIMARY LIPID SCORE
# ============================================================
#
# PRIMARY PRESPECIFIED LIPID SCORE:
#
# HALLMARK_FATTY_ACID_METABOLISM
#
# This is used to define lipid-high and lipid-low tumors.
#
# Median split:
#
#   LOW  = score <= median
#   HIGH = score > median
#
# ============================================================

primary_lipid_pathway <-
  "HALLMARK_FATTY_ACID_METABOLISM"


if (
  !primary_lipid_pathway %in%
  rownames(
    tcga_pathway_scores
  )
) {

  stop(
    paste0(
      "The prespecified primary lipid pathway ",
      primary_lipid_pathway,
      " is not available in the Step 9 GSVA matrix.\n\n",
      "It was therefore NOT substituted with another pathway."
    )
  )

}


primary_lipid_score <- as.numeric(
  tcga_pathway_scores[
    primary_lipid_pathway,
    ]
)


if (
  all(
    is.na(
      primary_lipid_score
    )
  )
) {

  stop(
    "The primary lipid pathway contains only NA scores."
  )

}


lipid_median <- median(
  primary_lipid_score,
  na.rm = TRUE
)


lipid_group <- ifelse(
  primary_lipid_score > lipid_median,
  "Lipid-high",
  "Lipid-low"
)


# ------------------------------------------------------------
# Convert to factor with explicit order
# ------------------------------------------------------------

lipid_group <- factor(
  lipid_group,
  levels = c(
    "Lipid-low",
    "Lipid-high"
  )
)


# ------------------------------------------------------------
# Check group sizes
# ------------------------------------------------------------

group_counts <- table(
  lipid_group,
  useNA = "ifany"
)


cat(
  "\n============================================================\n"
)

cat(
  "LIPID-HIGH / LIPID-LOW GROUPING\n"
)

cat(
  "============================================================\n\n"
)


cat(
  "Primary lipid pathway:",
  primary_lipid_pathway,
  "\n"
)


cat(
  "Median lipid score:",
  sprintf(
    "%.6f",
    lipid_median
  ),
  "\n\n"
)


print(
  group_counts
)


if (
  any(
    group_counts[
      c(
        "Lipid-low",
        "Lipid-high"
      )
    ] < 5
  )
) {

  stop(
    "One lipid group contains fewer than 5 samples. Differential pathway analysis is not reliable."
  )

}


# ============================================================
# 10.17 — CREATE DIFFERENTIAL PATHWAY ANALYSIS DATA
# ============================================================

all_pathways <- rownames(
  tcga_pathway_scores
)


differential_results_list <- vector(
  "list",
  length(
    all_pathways
  )
)


# ============================================================
# 10.18 — RANK-BISERIAL EFFECT SIZE FUNCTION
# ============================================================
#
# Rank-biserial correlation:
#
#   r_rb = (W_high - W_low) / (n_high * n_low)
#
# where W_high is the sum of ranks assigned to the
# high group.
#
# Equivalent sign convention:
#   positive = higher scores in lipid-high
#   negative = higher scores in lipid-low
#
# ============================================================

calculate_rank_biserial <- function(
  values,
  groups
) {

  complete <- complete.cases(
    values,
    groups
  )


  values <- values[
    complete
  ]


  groups <- groups[
    complete
  ]


  groups <- droplevels(
    groups
  )


  if (
    nlevels(groups) != 2
  ) {

    return(
      NA_real_
    )

  }


  low_values <- values[
    groups == "Lipid-low"
  ]


  high_values <- values[
    groups == "Lipid-high"
  ]


  n_low <- length(
    low_values
  )


  n_high <- length(
    high_values
  )


  if (
    n_low == 0
    ||
    n_high == 0
  ) {

    return(
      NA_real_
    )

  }


  combined_values <- c(
    low_values,
    high_values
  )


  combined_groups <- c(
    rep(
      "Lipid-low",
      n_low
    ),
    rep(
      "Lipid-high",
      n_high
    )
  )


  ranks <- rank(
    combined_values,
    ties.method = "average"
  )


  rank_sum_high <- sum(
    ranks[
      combined_groups ==
        "Lipid-high"
    ]
  )


  expected_high <- n_high * (
    n_high + 1
  ) / 2


  u_high <- rank_sum_high -
    expected_high


  rank_biserial <- (
    2 * u_high
  ) / (
    n_high * n_low
  )


  rank_biserial

}


# ============================================================
# 10.19 — RUN WILCOXON TEST FOR EVERY PATHWAY
# ============================================================

for (
  i in seq_along(
    all_pathways
  )
) {

  pathway_name <-
    all_pathways[i]


  values <- as.numeric(
    tcga_pathway_scores[
      pathway_name,
      ]
  )


  complete <- complete.cases(
    values,
    lipid_group
  )


  values_complete <- values[
    complete
  ]


  groups_complete <- droplevels(
    lipid_group[
      complete
    ]
  )


  n_low <- sum(
    groups_complete ==
      "Lipid-low"
  )


  n_high <- sum(
    groups_complete ==
      "Lipid-high"
  )


  median_low <- if (
    n_low > 0
  ) {

    median(
      values_complete[
        groups_complete ==
          "Lipid-low"
      ]
    )

  } else {

    NA_real_

  }


  median_high <- if (
    n_high > 0
  ) {

    median(
      values_complete[
        groups_complete ==
          "Lipid-high"
      ]
    )

  } else {

    NA_real_

  }


  mean_low <- if (
    n_low > 0
  ) {

    mean(
      values_complete[
        groups_complete ==
          "Lipid-low"
      ]
    )

  } else {

    NA_real_

  }


  mean_high <- if (
    n_high > 0
  ) {

    mean(
      values_complete[
        groups_complete ==
          "Lipid-high"
      ]
    )

  } else {

    NA_real_

  }


  if (
    n_low < 2
    ||
    n_high < 2
  ) {

    p_value <- NA_real_

    effect_size <- NA_real_

  } else {

    wilcox_result <- suppressWarnings(
      wilcox.test(
        values_complete[
          groups_complete ==
            "Lipid-high"
        ],
        values_complete[
          groups_complete ==
            "Lipid-low"
        ],
        alternative = "two.sided",
        exact = FALSE
      )
    )


    p_value <-
      wilcox_result$p.value


    effect_size <-
      calculate_rank_biserial(
        values_complete,
        groups_complete
      )

  }


  differential_results_list[[i]] <-
    data.frame(

      pathway_name =
        pathway_name,

      n_lipid_low =
        n_low,

      n_lipid_high =
        n_high,

      median_lipid_low =
        median_low,

      median_lipid_high =
        median_high,

      median_difference =
        median_high -
        median_low,

      mean_lipid_low =
        mean_low,

      mean_lipid_high =
        mean_high,

      p_value =
        p_value,

      rank_biserial =
        effect_size,

      stringsAsFactors = FALSE

    )

}


differential_pathway_results <- do.call(
  rbind,
  differential_results_list
)


rownames(
  differential_pathway_results
) <- NULL


# ============================================================
# 10.20 — BENJAMINI-HOCHBERG FDR
# ============================================================

differential_pathway_results$FDR <-
  NA_real_


valid_p <- !is.na(
  differential_pathway_results$p_value
)


differential_pathway_results$FDR[
  valid_p
] <- p.adjust(
  differential_pathway_results$p_value[
    valid_p
  ],
  method = "BH"
)


# ------------------------------------------------------------
# Absolute effect size
# ------------------------------------------------------------

differential_pathway_results$abs_rank_biserial <-
  abs(
    differential_pathway_results$rank_biserial
  )


# ------------------------------------------------------------
# Direction
# ------------------------------------------------------------

differential_pathway_results$direction <- ifelse(

  differential_pathway_results$rank_biserial > 0,

  "Higher in lipid-high",

  ifelse(
    differential_pathway_results$rank_biserial < 0,
    "Higher in lipid-low",
    "No difference"
  )

)


# ------------------------------------------------------------
# Significance
# ------------------------------------------------------------

differential_pathway_results$significance <- cut(

  differential_pathway_results$FDR,

  breaks = c(
    -Inf,
    0.001,
    0.01,
    0.05,
    Inf
  ),

  labels = c(
    "***",
    "**",
    "*",
    "ns"
  ),

  right = TRUE

)


# ============================================================
# 10.21 — ORDER DIFFERENTIAL RESULTS
# ============================================================

differential_pathway_results <-
  differential_pathway_results[
    order(
      differential_pathway_results$FDR,
      -differential_pathway_results$abs_rank_biserial,
      na.last = TRUE
    ),
    ,
    drop = FALSE
  ]


# ============================================================
# 10.22 — SAVE DIFFERENTIAL PATHWAY RESULTS
# ============================================================

write.csv(
  differential_pathway_results,
  file.path(
    results_dir,
    "10_differential_pathway_activity_lipid_high_vs_low.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 10.23 — PRINT TOP DIFFERENTIAL PATHWAYS
# ============================================================

cat(
  "\n============================================================\n"
)

cat(
  "DIFFERENTIAL PATHWAY ACTIVITY\n"
)

cat(
  "============================================================\n\n"
)


print(
  differential_pathway_results,
  row.names = FALSE
)


# ============================================================
# 10.24 — DIFFERENTIAL PATHWAY ACTIVITY PLOT
# ============================================================

plot_results <- differential_pathway_results[
  !is.na(
    differential_pathway_results$rank_biserial
  ),
  ,
  drop = FALSE
]


if (
  nrow(
    plot_results
  ) > 0
) {

  plot_results$pathway_name <- factor(
    plot_results$pathway_name,
    levels =
      plot_results$pathway_name[
        order(
          plot_results$rank_biserial
        )
      ]
  )


  differential_plot <- ggplot2::ggplot(
    plot_results,
    ggplot2::aes(
      x = rank_biserial,
      y = pathway_name
    )
  ) +

    ggplot2::geom_vline(
      xintercept = 0,
      linetype = "dashed"
    ) +

    ggplot2::geom_point(
      size = 3
    ) +

    ggplot2::labs(
      title =
        "Differential Pathway Activity by Lipid Status",

      subtitle =
        paste0(
          "Lipid-high vs lipid-low; primary score = ",
          primary_lipid_pathway
        ),

      x =
        "Rank-biserial correlation",

      y =
        "Pathway"
    ) +

    ggplot2::theme_minimal(
      base_size = 11
    )


  ggplot2::ggsave(
    filename = file.path(
      figures_dir,
      "10_differential_pathway_activity_lollipop.png"
    ),
    plot = differential_plot,
    width = 10,
    height = 8,
    dpi = 300
  )

}


# ============================================================
# 10.25 — SAVE LIPID GROUP ASSIGNMENTS
# ============================================================

lipid_group_df <- data.frame(

  sample_id =
    colnames(
      tcga_pathway_scores
    ),

  primary_lipid_score =
    primary_lipid_score,

  lipid_group =
    as.character(
      lipid_group
    ),

  stringsAsFactors = FALSE

)


write.csv(
  lipid_group_df,
  file.path(
    results_dir,
    "10_TCGA_lipid_high_low_group_assignments.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 10.26 — SAVE PRIMARY LIPID SCORE
# ============================================================

primary_lipid_definition <- data.frame(

  variable =
    "primary_lipid_pathway",

  pathway =
    primary_lipid_pathway,

  splitting_method =
    "Median split",

  median_cutoff =
    lipid_median,

  low_definition =
    "<= median",

  high_definition =
    "> median",

  stringsAsFactors = FALSE

)


write.csv(
  primary_lipid_definition,
  file.path(
    results_dir,
    "10_primary_lipid_score_definition.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 10.27 — SAVE COMPLETE STEP 10 OBJECT
# ============================================================

step10_results <- list(

  # ----------------------------------------------------------
  # Input
  # ----------------------------------------------------------

  tcga_pathway_scores =
    tcga_pathway_scores,

  # ----------------------------------------------------------
  # Pathway definitions
  # ----------------------------------------------------------

  lipid_pathways =
    lipid_present,

  ddr_pathways =
    ddr_present,

  checkpoint_pathways =
    checkpoint_present,

  # ----------------------------------------------------------
  # Correlation analysis
  # ----------------------------------------------------------

  lipid_ddr_correlations =
    lipid_ddr_correlations,

  correlation_matrix =
    correlation_matrix,

  fdr_matrix =
    fdr_matrix,

  representative_correlations =
    representative_correlations,

  # ----------------------------------------------------------
  # Lipid grouping
  # ----------------------------------------------------------

  primary_lipid_pathway =
    primary_lipid_pathway,

  primary_lipid_score =
    primary_lipid_score,

  lipid_median =
    lipid_median,

  lipid_group =
    lipid_group,

  # ----------------------------------------------------------
  # Differential pathway analysis
  # ----------------------------------------------------------

  differential_pathway_results =
    differential_pathway_results

)


saveRDS(
  step10_results,
  file.path(
    processed_dir,
    "step10_pathway_correlation_and_differential_activity.rds"
  )
)


# ============================================================
# 10.28 — FINAL QC
# ============================================================

cat(
  "\n============================================================\n"
)

cat(
  "STEP 10 FINAL QC\n"
)

cat(
  "============================================================\n\n"
)


# ------------------------------------------------------------
# Correlation QC
# ------------------------------------------------------------

cat(
  "Lipid pathways tested:",
  length(
    lipid_present
  ),
  "\n"
)


cat(
  "DDR pathways tested:",
  length(
    ddr_present
  ),
  "\n"
)


cat(
  "Total lipid × DDR correlations:",
  nrow(
    lipid_ddr_correlations
  ),
  "\n"
)


cat(
  "Significant correlations (BH-FDR < 0.05):",
  sum(
    lipid_ddr_correlations$FDR < 0.05,
    na.rm = TRUE
  ),
  "\n\n"
)


# ------------------------------------------------------------
# Differential analysis QC
# ------------------------------------------------------------

cat(
  "Total pathways tested for differential activity:",
  nrow(
    differential_pathway_results
  ),
  "\n"
)


cat(
  "Significant pathways (BH-FDR < 0.05):",
  sum(
    differential_pathway_results$FDR < 0.05,
    na.rm = TRUE
  ),
  "\n\n"
)


cat(
  "Lipid-low samples:",
  sum(
    lipid_group ==
      "Lipid-low"
  ),
  "\n"
)


cat(
  "Lipid-high samples:",
  sum(
    lipid_group ==
      "Lipid-high"
  ),
  "\n\n"
)


# ------------------------------------------------------------
# Check for missing p-values
# ------------------------------------------------------------

if (
  any(
    is.na(
      lipid_ddr_correlations$p_value
    )
  )
) {

  cat(
    "WARNING: Some correlation pairs had insufficient/non-variable data.\n"
  )

}


if (
  any(
    is.na(
      differential_pathway_results$p_value
    )
  )
) {

  cat(
    "WARNING: Some pathways had insufficient/non-variable data for Wilcoxon testing.\n"
  )

}


# ============================================================
# 10.29 — FINAL OUTPUT SUMMARY
# ============================================================

cat(
  "\n============================================================\n"
)

cat(
  "STEP 10 COMPLETED SUCCESSFULLY\n"
)

cat(
  "============================================================\n\n"
)


cat(
  "Primary analyses completed:\n"
)

cat(
  "  1. Spearman lipid × DDR correlations\n"
)

cat(
  "  2. Benjamini-Hochberg FDR correction\n"
)

cat(
  "  3. Correlation heatmap\n"
)

cat(
  "  4. Representative scatterplots\n"
)

cat(
  "  5. Lipid-high vs lipid-low grouping\n"
)

cat(
  "  6. Wilcoxon pathway activity testing\n"
)

cat(
  "  7. Rank-biserial effect sizes\n"
)

cat(
  "  8. Benjamini-Hochberg FDR correction\n\n"
)


cat(
  "Main output files:\n\n"
)


cat(
  file.path(
    results_dir,
    "10_lipid_DDR_spearman_correlations.csv"
  ),
  "\n"
)


cat(
  file.path(
    results_dir,
    "10_lipid_DDR_spearman_rho_matrix.csv"
  ),
  "\n"
)


cat(
  file.path(
    results_dir,
    "10_lipid_DDR_BH_FDR_matrix.csv"
  ),
  "\n"
)


cat(
  file.path(
    results_dir,
    "10_differential_pathway_activity_lipid_high_vs_low.csv"
  ),
  "\n"
)


cat(
  file.path(
    results_dir,
    "10_TCGA_lipid_high_low_group_assignments.csv"
  ),
  "\n"
)


cat(
  file.path(
    results_dir,
    "10_primary_lipid_score_definition.csv"
  ),
  "\n"
)


cat(
  file.path(
    processed_dir,
    "step10_pathway_correlation_and_differential_activity.rds"
  ),
  "\n\n"
)


cat(
  "Figures saved in:\n",
  figures_dir,
  "\n\n"
)


cat(
  "============================================================\n"
)

cat(
  "END OF STEP 10\n"
)

cat(
  "============================================================\n"
)
