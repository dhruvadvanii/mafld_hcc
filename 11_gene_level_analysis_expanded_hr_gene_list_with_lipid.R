# ============================================================
# STEP 11 — 2.5 GENE-LEVEL ANALYSIS
# FINAL CLEAN VERSION
#
# Important:
# - Does NOT use the incorrect 115-sample lipid_score
# - Does NOT require gsva_scores
# - Does NOT require an existing lipid_score
# - Reconstructs the lipid score for the current 371 TCGA samples
# - Uses the existing HALLMARK_FATTY_ACID_METABOLISM gene set
# - Uses the existing filtered_counts_object
# - Contains NO [[ ]] syntax
# ============================================================

cat("\n============================================================\n")
cat("STEP 11 — 2.5 GENE-LEVEL ANALYSIS\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 11.1 — Required objects
# ------------------------------------------------------------

required_objects <- c(
    "final_gene_sets",
    "filtered_counts_object",
    "tcga_symbols"
)

missing_objects <- required_objects[
    !vapply(
        required_objects,
        exists,
        logical(1),
        inherits = TRUE
    )
]

if (length(missing_objects) > 0) {

    stop(
        paste0(
            "Missing required objects: ",
            paste(
                missing_objects,
                collapse = ", "
            )
        )
    )

}

cat("Required objects are present.\n")


# ------------------------------------------------------------
# 11.2 — Get the EXISTING filtered TCGA expression matrix
# ------------------------------------------------------------

tcga_expression <-
    as.matrix(
        filtered_counts_object$counts
    )

storage.mode(
    tcga_expression
) <-
    "numeric"


# ------------------------------------------------------------
# 11.3 — Verify sample number
# ------------------------------------------------------------

cat(
    "TCGA expression dimensions: ",
    nrow(tcga_expression),
    " genes x ",
    ncol(tcga_expression),
    " samples\n",
    sep = ""
)


if (
    is.null(
        colnames(
            tcga_expression
        )
    )
) {

    stop(
        "The existing TCGA expression matrix has no sample identifiers."
    )

}


# ------------------------------------------------------------
# 11.4 — Assign gene symbols
# ------------------------------------------------------------

if (
    length(tcga_symbols) !=
    nrow(tcga_expression)
) {

    stop(
        paste0(
            "tcga_symbols length (",
            length(tcga_symbols),
            ") does not match expression rows (",
            nrow(tcga_expression),
            ")."
        )
    )

}

rownames(
    tcga_expression
) <-
    toupper(
        trimws(
            as.character(
                tcga_symbols
            )
        )
    )


# ------------------------------------------------------------
# 11.5 — Remove invalid gene symbols
# ------------------------------------------------------------

valid_gene_rows <-
    !is.na(
        rownames(
            tcga_expression
        )
    ) &
    rownames(
        tcga_expression
    ) != ""

tcga_expression <-
    tcga_expression[
        valid_gene_rows,
        ,
        drop = FALSE
    ]


# ------------------------------------------------------------
# 11.6 — Collapse duplicate gene symbols
# ------------------------------------------------------------

if (
    any(
        duplicated(
            rownames(
                tcga_expression
            )
        )
    )
) {

    unique_genes <-
        unique(
            rownames(
                tcga_expression
            )
        )

    collapsed_expression <-
        matrix(
            NA_real_,
            nrow =
                length(unique_genes),
            ncol =
                ncol(
                    tcga_expression
                ),
            dimnames =
                list(
                    unique_genes,
                    colnames(
                        tcga_expression
                    )
                )
        )

    for (
        i in seq_along(
            unique_genes
        )
    ) {

        current_gene <-
            unique_genes[i]

        current_rows <-
            which(
                rownames(
                    tcga_expression
                ) ==
                current_gene
            )

        if (
            length(current_rows) == 1
        ) {

            collapsed_expression[i, ] <-
                tcga_expression[
                    current_rows,
                    ]

        } else {

            collapsed_expression[i, ] <-
                apply(
                    tcga_expression[
                        current_rows,
                        ,
                        drop = FALSE
                    ],
                    2,
                    median,
                    na.rm = TRUE
                )

        }

    }

    tcga_expression <-
        collapsed_expression

}


# ------------------------------------------------------------
# 11.7 — Obtain EXACT fatty-acid pathway gene set
# ------------------------------------------------------------

lipid_pathway_name <-
    "HALLMARK_FATTY_ACID_METABOLISM"

if (
    !lipid_pathway_name %in%
    names(final_gene_sets)
) {

    stop(
        paste0(
            lipid_pathway_name,
            " is not present in final_gene_sets."
        )
    )

}

fatty_acid_genes <-
    getElement(
        final_gene_sets,
        lipid_pathway_name
    )

fatty_acid_genes <-
    unique(
        toupper(
            trimws(
                as.character(
                    fatty_acid_genes
                )
            )
        )
    )

fatty_acid_genes <-
    fatty_acid_genes[
        !is.na(
            fatty_acid_genes
        ) &
            fatty_acid_genes != ""
    ]


cat(
    "Fatty-acid pathway genes: ",
    length(
        fatty_acid_genes
    ),
    "\n",
    sep = ""
)


# ------------------------------------------------------------
# 11.8 — Match fatty-acid genes to TCGA expression
# ------------------------------------------------------------

fatty_acid_present <-
    fatty_acid_genes %in%
    rownames(
        tcga_expression
    )

fatty_acid_genes_matched <-
    fatty_acid_genes[
        fatty_acid_present
    ]

fatty_acid_genes_missing <-
    fatty_acid_genes[
        !fatty_acid_present
    ]


cat(
    "Fatty-acid genes present in TCGA: ",
    length(
        fatty_acid_genes_matched
    ),
    " / ",
    length(
        fatty_acid_genes
    ),
    "\n",
    sep = ""
)


if (
    length(
        fatty_acid_genes_matched
    ) < 10
) {

    stop(
        "Fewer than 10 fatty-acid pathway genes matched the TCGA expression matrix."
    )

}


# ------------------------------------------------------------
# 11.9 — Calculate lipid score for ALL current TCGA samples
# ------------------------------------------------------------
#
# We intentionally do NOT use the existing 115-sample
# lipid_score object.
#
# The score is reconstructed directly from the current
# 371-sample expression matrix.
# ------------------------------------------------------------

fatty_acid_expression <-
    tcga_expression[
        fatty_acid_genes_matched,
        ,
        drop = FALSE
    ]


if (
    requireNamespace(
        "GSVA",
        quietly = TRUE
    )
) {

    cat(
        "GSVA package detected.\n"
    )

    cat(
        "Calculating HALLMARK_FATTY_ACID_METABOLISM score ",
        "for all current TCGA samples...\n"
    )


    lipid_gene_set <-
        list(
            HALLMARK_FATTY_ACID_METABOLISM =
                fatty_acid_genes_matched
        )


    gsva_parameters <-
        tryCatch(
            GSVA::gsvaParam(
                exprData =
                    tcga_expression,
                geneSets =
                    lipid_gene_set,
                kcdf =
                    "Gaussian"
            ),
            error = function(e) NULL
        )


    lipid_gsva_result <-
        NULL


    if (
        !is.null(
            gsva_parameters
        )
    ) {

        lipid_gsva_result <-
            tryCatch(
                GSVA::gsva(
                    gsva_parameters,
                    verbose = FALSE
                ),
                error = function(e) NULL
            )

    }


    if (
        is.null(
            lipid_gsva_result
        )
    ) {

        stop(
            "GSVA was available but the lipid pathway score could not be calculated."
        )

    }


    lipid_score <-
        as.numeric(
            lipid_gsva_result[
                "HALLMARK_FATTY_ACID_METABOLISM",
                ]
        )

    names(
        lipid_score
    ) <-
        colnames(
            tcga_expression
        )

    lipid_score_source <-
        "reconstructed from current TCGA expression using GSVA"


} else {

    stop(
        "GSVA package is not available. ",
        "Install/load GSVA before running this step."
    )

}


# ------------------------------------------------------------
# 11.10 — HARD CHECK: lipid score must correspond to
#          current TCGA samples
# ------------------------------------------------------------

if (
    length(lipid_score) !=
    ncol(tcga_expression)
) {

    stop(
        paste0(
            "CRITICAL QC ERROR: lipid score contains ",
            length(lipid_score),
            " samples but current TCGA expression contains ",
            ncol(tcga_expression),
            " samples."
        )
    )

}


if (
    is.null(
        names(lipid_score)
    )
) {

    stop(
        "CRITICAL QC ERROR: lipid score has no sample identifiers."
    )

}


if (
    !identical(
        names(lipid_score),
        colnames(tcga_expression)
    )
) {

    stop(
        "CRITICAL QC ERROR: lipid-score sample identifiers do not exactly match TCGA expression sample identifiers."
    )

}


cat(
    "Lipid score successfully reconstructed for ",
    length(lipid_score),
    " TCGA samples.\n",
    sep = ""
)


# ------------------------------------------------------------
# 11.11 — Existing Reactome HR gene set
# ------------------------------------------------------------

if (
    !"HOMOLOGOUS_RECOMBINATION" %in%
    names(final_gene_sets)
) {

    stop(
        "HOMOLOGOUS_RECOMBINATION is not present in final_gene_sets."
    )

}

reactome_hr_genes <-
    getElement(
        final_gene_sets,
        "HOMOLOGOUS_RECOMBINATION"
    )

reactome_hr_genes <-
    unique(
        toupper(
            trimws(
                as.character(
                    reactome_hr_genes
                )
            )
        )
    )

reactome_hr_genes <-
    reactome_hr_genes[
        !is.na(
            reactome_hr_genes
        ) &
            reactome_hr_genes != ""
    ]


# ------------------------------------------------------------
# 11.12 — Additional HR/DDR genes supplied by user
# ------------------------------------------------------------

additional_hr_genes <-
    c(
        "RAD51",
        "BRCA1",
        "BRCA2",
        "PALB2",
        "RAD51B",
        "RAD51C",
        "RAD51D",
        "XRCC2",
        "XRCC3",
        "BARD1",
        "FAM175A",
        "UIMC1",
        "BRIP1",
        "RBBP8",
        "MRE11A",
        "RAD50",
        "NBN",
        "ATM",
        "ATR",
        "ATRIP",
        "TOPBP1",
        "TP53BP1",
        "EXO1",
        "BLM",
        "DNA2",
        "RPA1",
        "RPA2",
        "RPA3",
        "TOP3A",
        "RMI1",
        "RMI2",
        "MUS81",
        "EME1",
        "GEN1",
        "SLX1A",
        "SLX4",
        "RIF1",
        "DMC1",
        "DSS1",
        "UBE2N",
        "RNF8",
        "RNF168",
        "MDC1",
        "H2AX",
        "XRCC6",
        "XRCC5",
        "LIG4",
        "XRCC4",
        "PARP1",
        "RAD54A",
        "RAD52",
        "SPO11"
    )

additional_hr_genes <-
    unique(
        toupper(
            additional_hr_genes
        )
    )


# ------------------------------------------------------------
# 11.13 — Expanded HR gene panel
# ------------------------------------------------------------

expanded_hr_genes <-
    unique(
        c(
            reactome_hr_genes,
            additional_hr_genes
        )
    )

expanded_hr_genes <-
    expanded_hr_genes[
        !is.na(
            expanded_hr_genes
        ) &
            expanded_hr_genes != ""
    ]


cat(
    "Reactome HR genes: ",
    length(
        reactome_hr_genes
    ),
    "\n",
    sep = ""
)

cat(
    "Additional HR/DDR genes: ",
    length(
        additional_hr_genes
    ),
    "\n",
    sep = ""
)

cat(
    "Expanded HR panel: ",
    length(
        expanded_hr_genes
    ),
    "\n",
    sep = ""
)


# ------------------------------------------------------------
# 11.14 — Match expanded HR genes to TCGA
# ------------------------------------------------------------

expanded_hr_present <-
    expanded_hr_genes %in%
    rownames(
        tcga_expression
    )

expanded_hr_matched <-
    expanded_hr_genes[
        expanded_hr_present
    ]

expanded_hr_missing <-
    expanded_hr_genes[
        !expanded_hr_present
    ]


cat(
    "Expanded HR genes present in TCGA: ",
    length(
        expanded_hr_matched
    ),
    " / ",
    length(
        expanded_hr_genes
    ),
    "\n",
    sep = ""
)


# ------------------------------------------------------------
# 11.15 — Extract HR expression
# ------------------------------------------------------------

hr_expression <-
    tcga_expression[
        expanded_hr_matched,
        ,
        drop = FALSE
    ]


# ------------------------------------------------------------
# 11.16 — Verify exact sample matching
# ------------------------------------------------------------

if (
    !identical(
        colnames(hr_expression),
        names(lipid_score)
    )
) {

    stop(
        "CRITICAL QC ERROR: HR expression samples and lipid-score samples do not match."
    )

}


cat(
    "Exact expression/lipid sample match confirmed.\n"
)


# ------------------------------------------------------------
# 11.17 — Gene-level Spearman correlations
# ------------------------------------------------------------

gene_results_list <-
    vector(
        "list",
        nrow(
            hr_expression
        )
    )


for (
    i in seq_len(
        nrow(
            hr_expression
        )
    )
) {

    gene_name <-
        rownames(
            hr_expression
        )[i]

    gene_values <-
        as.numeric(
            hr_expression[
                i,
                ]
            )

    valid <-
        is.finite(
            gene_values
        ) &
        is.finite(
            lipid_score
        )

    n_valid <-
        sum(
            valid
        )

    if (
        n_valid < 10
    ) {

        gene_results_list[[i]] <-
            data.frame(
                Gene =
                    gene_name,
                N =
                    n_valid,
                Spearman_rho =
                    NA_real_,
                P_value =
                    NA_real_,
                stringsAsFactors =
                    FALSE
            )

        next

    }


    x <-
        gene_values[
            valid
        ]

    y <-
        lipid_score[
            valid
        ]


    if (
        length(
            unique(x)
        ) < 2 ||
        length(
            unique(y)
        ) < 2
    ) {

        gene_results_list[[i]] <-
            data.frame(
                Gene =
                    gene_name,
                N =
                    n_valid,
                Spearman_rho =
                    NA_real_,
                P_value =
                    NA_real_,
                stringsAsFactors =
                    FALSE
            )

        next

    }


    correlation_test <-
        suppressWarnings(
            cor.test(
                x,
                y,
                method =
                    "spearman",
                exact =
                    FALSE
            )
        )


    gene_results_list[[i]] <-
        data.frame(
            Gene =
                gene_name,
            N =
                n_valid,
            Spearman_rho =
                as.numeric(
                    correlation_test$estimate
                ),
            P_value =
                correlation_test$p.value,
            stringsAsFactors =
                FALSE
        )

}


gene_results <-
    do.call(
        rbind,
        gene_results_list
    )

rownames(
    gene_results
) <-
    NULL


# ------------------------------------------------------------
# 11.18 — BH-FDR
# ------------------------------------------------------------

gene_results$BH_FDR <-
    NA_real_

valid_p <-
    !is.na(
        gene_results$P_value
    ) &
    is.finite(
        gene_results$P_value
    )

gene_results$BH_FDR[
    valid_p
] <-
    p.adjust(
        gene_results$P_value[
            valid_p
        ],
        method =
            "BH"
    )


# ------------------------------------------------------------
# 11.19 — Direction
# ------------------------------------------------------------

gene_results$Direction <-
    "Not estimable"

gene_results$Direction[
    !is.na(
        gene_results$Spearman_rho
    ) &
        gene_results$Spearman_rho > 0
] <-
    "Positive"

gene_results$Direction[
    !is.na(
        gene_results$Spearman_rho
    ) &
        gene_results$Spearman_rho < 0
] <-
    "Negative"


# ------------------------------------------------------------
# 11.20 — Significance
# ------------------------------------------------------------

gene_results$Association <-
    "Non-significant"

gene_results$Association[
    gene_results$BH_FDR < 0.05 &
        gene_results$Spearman_rho > 0
] <-
    "Significantly positive"

gene_results$Association[
    gene_results$BH_FDR < 0.05 &
        gene_results$Spearman_rho < 0
] <-
    "Significantly negative"


# ------------------------------------------------------------
# 11.21 — Rank
# ------------------------------------------------------------

gene_results$Absolute_rho <-
    abs(
        gene_results$Spearman_rho
    )

gene_results <-
    gene_results[
        order(
            -gene_results$Absolute_rho,
            gene_results$BH_FDR,
            gene_results$Gene
        ),
        ,
        drop = FALSE
    ]

rownames(
    gene_results
) <-
    NULL


# ------------------------------------------------------------
# 11.22 — Positive and negative genes
# ------------------------------------------------------------

positive_genes <-
    gene_results[
        gene_results$BH_FDR < 0.05 &
            gene_results$Spearman_rho > 0,
        ,
        drop = FALSE
    ]

negative_genes <-
    gene_results[
        gene_results$BH_FDR < 0.05 &
            gene_results$Spearman_rho < 0,
        ,
        drop = FALSE
    ]


# ------------------------------------------------------------
# 11.23 — Save ranked gene list
# ------------------------------------------------------------

write.csv(
    gene_results,
    file.path(
        "results",
        "Step11_2.5_expanded_HR_gene_lipid_correlation_results.csv"
    ),
    row.names = FALSE
)


# ------------------------------------------------------------
# 11.24 — Save significant gene lists
# ------------------------------------------------------------

write.csv(
    positive_genes,
    file.path(
        "results",
        "Step11_2.5_expanded_HR_positive_genes.csv"
    ),
    row.names = FALSE
)

write.csv(
    negative_genes,
    file.path(
        "results",
        "Step11_2.5_expanded_HR_negative_genes.csv"
    ),
    row.names = FALSE
)


# ------------------------------------------------------------
# 11.25 — Gene audit
# ------------------------------------------------------------

gene_audit <-
    data.frame(
        Gene =
            expanded_hr_genes,
        In_Reactome_HR =
            expanded_hr_genes %in%
            reactome_hr_genes,
        Present_in_TCGA =
            expanded_hr_genes %in%
            rownames(
                tcga_expression
            ),
        stringsAsFactors =
            FALSE
    )


write.csv(
    gene_audit,
    file.path(
        "results",
        "Step11_2.5_expanded_HR_gene_audit.csv"
    ),
    row.names = FALSE
)


# ------------------------------------------------------------
# 11.26 — Summary
# ------------------------------------------------------------

n_tested <-
    sum(
        !is.na(
            gene_results$P_value
        )
    )

n_positive <-
    nrow(
        positive_genes
    )

n_negative <-
    nrow(
        negative_genes
    )

n_nonsignificant <-
    sum(
        gene_results$BH_FDR >= 0.05,
        na.rm = TRUE
    )


gene_summary <-
    data.frame(
        Metric = c(
            "TCGA samples",
            "Reactome HR genes",
            "Additional HR/DDR genes",
            "Expanded HR genes",
            "Expanded HR genes present in TCGA",
            "Genes tested",
            "Significantly positive genes",
            "Significantly negative genes",
            "Non-significant genes"
        ),
        Value = c(
            ncol(
                tcga_expression
            ),
            length(
                reactome_hr_genes
            ),
            length(
                additional_hr_genes
            ),
            length(
                expanded_hr_genes
            ),
            length(
                expanded_hr_matched
            ),
            n_tested,
            n_positive,
            n_negative,
            n_nonsignificant
        ),
        stringsAsFactors =
            FALSE
    )


write.csv(
    gene_summary,
    file.path(
        "results",
        "Step11_2.5_expanded_HR_gene_summary.csv"
    ),
    row.names = FALSE
)


# ------------------------------------------------------------
# 11.27 — Save Step 11 analysis
# ------------------------------------------------------------

save(
    lipid_score,
    lipid_score_source,
    fatty_acid_genes,
    fatty_acid_genes_matched,
    fatty_acid_genes_missing,
    reactome_hr_genes,
    additional_hr_genes,
    expanded_hr_genes,
    expanded_hr_matched,
    expanded_hr_missing,
    hr_expression,
    gene_results,
    positive_genes,
    negative_genes,
    gene_summary,
    file =
        file.path(
            "results",
            "Step11_2.5_expanded_HR_gene_lipid_analysis.RData"
        )
)


# ------------------------------------------------------------
# 11.28 — Final QC
# ------------------------------------------------------------

if (
    length(lipid_score) !=
    ncol(tcga_expression)
) {

    stop(
        "FINAL QC FAILED: lipid-score sample count does not match TCGA sample count."
    )

}

if (
    !identical(
        names(lipid_score),
        colnames(tcga_expression)
    )
) {

    stop(
        "FINAL QC FAILED: lipid-score sample identifiers do not match TCGA sample identifiers."
    )

}

if (
    anyDuplicated(
        gene_results$Gene
    ) > 0
) {

    stop(
        "FINAL QC FAILED: duplicate genes detected."
    )

}

if (
    any(
        !is.na(
            gene_results$Spearman_rho
        ) &
        !is.finite(
            gene_results$Spearman_rho
        )
    )
) {

    stop(
        "FINAL QC FAILED: non-finite Spearman correlations detected."
    )

}


# ------------------------------------------------------------
# 11.29 — Final output
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP 11 / 2.5 COMPLETED SUCCESSFULLY\n")
cat("============================================================\n")
cat(
    "Lipid score source: ",
    lipid_score_source,
    "\n",
    sep = ""
)
cat(
    "TCGA samples: ",
    ncol(
        tcga_expression
    ),
    "\n",
    sep = ""
)
cat(
    "Fatty-acid genes used: ",
    length(
        fatty_acid_genes_matched
    ),
    "\n",
    sep = ""
)
cat(
    "Reactome HR genes: ",
    length(
        reactome_hr_genes
    ),
    "\n",
    sep = ""
)
cat(
    "Expanded HR genes: ",
    length(
        expanded_hr_genes
    ),
    "\n",
    sep = ""
)
cat(
    "Expanded HR genes tested: ",
    n_tested,
    "\n",
    sep = ""
)
cat(
    "Significantly positive: ",
    n_positive,
    "\n",
    sep = ""
)
cat(
    "Significantly negative: ",
    n_negative,
    "\n",
    sep = ""
)
cat("============================================================\n")

# ============================================================
# STEP 11 FIGURES — FINAL FIGURE GENERATION
#
# Uses the already completed Step 11 results.
# Does NOT rerun correlations.
# Does NOT recalculate the lipid score.
# Does NOT modify gene_results.
#
# Outputs:
# 1. Volcano plot
# 2. Full HR/DDR correlation heatmap
# 3. Significant HR/DDR correlation heatmap, if applicable
# ============================================================

cat("\n============================================================\n")
cat("STEP 11 — FIGURE GENERATION\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Confirm existing Step 11 results
# ------------------------------------------------------------

if (
    !exists("gene_results")
) {

    stop(
        "gene_results is not present in the current R session. ",
        "The previously saved Step 11 RData file will be loaded."
    )

}


# ------------------------------------------------------------
# 2. Confirm results directory
# ------------------------------------------------------------

if (
    !dir.exists(
        "results"
    )
) {

    dir.create(
        "results",
        recursive = TRUE
    )

}


# ------------------------------------------------------------
# 3. Check required columns
# ------------------------------------------------------------

required_columns <-
    c(
        "Gene",
        "Spearman_rho",
        "P_value",
        "BH_FDR"
    )

missing_columns <-
    setdiff(
        required_columns,
        colnames(
            gene_results
        )
    )

if (
    length(
        missing_columns
    ) > 0
) {

    stop(
        paste0(
            "The existing gene_results object is missing: ",
            paste(
                missing_columns,
                collapse = ", "
            )
        )
    )

}


# ------------------------------------------------------------
# 4. Load plotting packages
# ------------------------------------------------------------

if (
    !requireNamespace(
        "ggplot2",
        quietly = TRUE
    )
) {

    stop(
        "ggplot2 is required. Please install ggplot2 and rerun this figure block."
    )

}

if (
    !requireNamespace(
        "pheatmap",
        quietly = TRUE
    )
) {

    stop(
        "pheatmap is required. Please install pheatmap and rerun this figure block."
    )

}

library(
    ggplot2
)

library(
    pheatmap
)


# ============================================================
# FIGURE 1 — VOLCANO PLOT
# ============================================================

cat(
    "Generating volcano plot...\n"
)


volcano_data <-
    gene_results


# ------------------------------------------------------------
# Calculate -log10 BH-FDR
# ------------------------------------------------------------

volcano_data$neg_log10_BH_FDR <-
    NA_real_


valid_volcano <-
    !is.na(
        volcano_data$BH_FDR
    ) &
    is.finite(
        volcano_data$BH_FDR
    ) &
    volcano_data$BH_FDR > 0


volcano_data$neg_log10_BH_FDR[
    valid_volcano
] <-
    -log10(
        volcano_data$BH_FDR[
            valid_volcano
        ]
    )


# ------------------------------------------------------------
# Define volcano groups
# ------------------------------------------------------------

volcano_data$Volcano_group <-
    "Non-significant"


volcano_data$Volcano_group[
    !is.na(
        volcano_data$BH_FDR
    ) &
        volcano_data$BH_FDR < 0.05 &
        volcano_data$Spearman_rho > 0
] <-
    "Significantly positive"


volcano_data$Volcano_group[
    !is.na(
        volcano_data$BH_FDR
    ) &
        volcano_data$BH_FDR < 0.05 &
        volcano_data$Spearman_rho < 0
] <-
    "Significantly negative"


# ------------------------------------------------------------
# Prevent infinite y-values
# ------------------------------------------------------------

finite_volcano <-
    is.finite(
        volcano_data$Spearman_rho
    ) &
    is.finite(
        volcano_data$neg_log10_BH_FDR
    )


volcano_data_plot <-
    volcano_data[
        finite_volcano,
        ,
        drop = FALSE
    ]


if (
    nrow(
        volcano_data_plot
    ) == 0
) {

    stop(
        "No valid data are available for the volcano plot."
    )

}


# ------------------------------------------------------------
# Create volcano plot
# ------------------------------------------------------------

volcano_plot <-
    ggplot(
        volcano_data_plot,
        aes(
            x =
                Spearman_rho,
            y =
                neg_log10_BH_FDR,
            shape =
                Volcano_group
        )
    ) +

    geom_point(
        size = 3,
        alpha = 0.85
    ) +

    geom_vline(
        xintercept = 0,
        linetype = "dashed"
    ) +

    geom_hline(
        yintercept =
            -log10(0.05),
        linetype = "dashed"
    ) +

    labs(
        title =
            "Expanded HR/DDR Gene Correlation with Lipid Metabolism",
        subtitle =
            "HALLMARK_FATTY_ACID_METABOLISM pathway score",
        x =
            "Spearman correlation coefficient (rho)",
        y =
            "-log10(BH-FDR)",
        shape =
            "Association"
    ) +

    theme_classic(
        base_size = 14
    ) +

    theme(
        plot.title =
            element_text(
                face = "bold",
                size = 16
            ),
        plot.subtitle =
            element_text(
                size = 11
            ),
        axis.title =
            element_text(
                face = "bold"
            ),
        legend.title =
            element_text(
                face = "bold"
            ),
        legend.position =
            "right"
    )


# ------------------------------------------------------------
# Save volcano plot
# ------------------------------------------------------------

volcano_file <-
    file.path(
        "results",
        "Step11_2.5_expanded_HR_gene_lipid_volcano.png"
    )


ggsave(
    filename =
        volcano_file,
    plot =
        volcano_plot,
    width =
        10,
    height =
        8,
    units =
        "in",
    dpi =
        300,
    bg =
        "white"
)


cat(
    "Volcano plot saved:\n",
    normalizePath(
        volcano_file,
        winslash = "/",
        mustWork = FALSE
    ),
    "\n\n"
)


# ============================================================
# FIGURE 2 — FULL HR/DDR CORRELATION HEATMAP
# ============================================================

cat(
    "Generating full HR/DDR correlation heatmap...\n"
)


heatmap_data <-
    gene_results[
        !is.na(
            gene_results$Spearman_rho
        ) &
        is.finite(
            gene_results$Spearman_rho
        ),
        ,
        drop = FALSE
    ]


if (
    nrow(
        heatmap_data
    ) == 0
) {

    stop(
        "No valid correlation coefficients are available for the heatmap."
    )

}


# ------------------------------------------------------------
# Order genes from negative to positive correlation
# ------------------------------------------------------------

heatmap_data <-
    heatmap_data[
        order(
            heatmap_data$Spearman_rho,
            decreasing = FALSE
        ),
        ,
        drop = FALSE
    ]


# ------------------------------------------------------------
# Create one-row correlation matrix
# ------------------------------------------------------------

heatmap_matrix <-
    matrix(
        heatmap_data$Spearman_rho,
        nrow = 1
    )


colnames(
    heatmap_matrix
) <-
    heatmap_data$Gene


rownames(
    heatmap_matrix
) <-
    "Spearman rho"


# ------------------------------------------------------------
# Dynamic figure width
# ------------------------------------------------------------

heatmap_width <-
    max(
        14,
        min(
            30,
            0.28 *
                ncol(
                    heatmap_matrix
                ) +
                6
        )
    )


# ------------------------------------------------------------
# Save full heatmap
# ------------------------------------------------------------

full_heatmap_file <-
    file.path(
        "results",
        "Step11_2.5_expanded_HR_gene_lipid_heatmap.png"
    )


png(
    filename =
        full_heatmap_file,
    width =
        heatmap_width,
    height =
        5.5,
    units =
        "in",
    res =
        300
)


pheatmap(
    heatmap_matrix,
    cluster_rows =
        FALSE,
    cluster_cols =
        FALSE,
    show_rownames =
        TRUE,
    show_colnames =
        TRUE,
    angle_col =
        90,
    fontsize =
        9,
    fontsize_col =
        7,
    display_numbers =
        TRUE,
    number_format =
        "%.2f",
    border_color =
        NA,
    main =
        "Expanded HR/DDR Gene Correlation with Lipid Metabolism"
)


dev.off()


cat(
    "Full heatmap saved:\n",
    normalizePath(
        full_heatmap_file,
        winslash = "/",
        mustWork = FALSE
    ),
    "\n\n"
)


# ============================================================
# FIGURE 3 — SIGNIFICANT-GENE HEATMAP
# ============================================================

cat(
    "Checking for significant genes...\n"
)


significant_heatmap_data <-
    gene_results[
        !is.na(
            gene_results$BH_FDR
        ) &
        gene_results$BH_FDR < 0.05 &
        !is.na(
            gene_results$Spearman_rho
        ) &
        is.finite(
            gene_results$Spearman_rho
        ),
        ,
        drop = FALSE
    ]


if (
    nrow(
        significant_heatmap_data
    ) == 0
) {

    cat(
        "No genes have BH-FDR < 0.05.\n"
    )

    cat(
        "Therefore, a significant-gene heatmap was not generated.\n\n"
    )

} else {

    significant_heatmap_data <-
        significant_heatmap_data[
            order(
                significant_heatmap_data$Spearman_rho,
                decreasing = FALSE
            ),
            ,
            drop = FALSE
        ]


    significant_heatmap_matrix <-
        matrix(
            significant_heatmap_data$Spearman_rho,
            nrow = 1
        )


    colnames(
        significant_heatmap_matrix
    ) <-
        significant_heatmap_data$Gene


    rownames(
        significant_heatmap_matrix
    ) <-
        "Spearman rho"


    significant_heatmap_width <-
        max(
            10,
            min(
                30,
                0.35 *
                    ncol(
                        significant_heatmap_matrix
                    ) +
                    6
            )
        )


    significant_heatmap_file <-
        file.path(
            "results",
            "Step11_2.5_expanded_HR_significant_gene_lipid_heatmap.png"
        )


    png(
        filename =
            significant_heatmap_file,
        width =
            significant_heatmap_width,
        height =
            5.5,
        units =
            "in",
        res =
            300
    )


    pheatmap(
        significant_heatmap_matrix,
        cluster_rows =
            FALSE,
        cluster_cols =
            FALSE,
        show_rownames =
            TRUE,
        show_colnames =
            TRUE,
        angle_col =
            90,
        fontsize =
            9,
        fontsize_col =
            8,
        display_numbers =
            TRUE,
        number_format =
            "%.2f",
        border_color =
            NA,
        main =
            "Significant HR/DDR Genes Correlated with Lipid Metabolism"
    )


    dev.off()


    cat(
        "Significant-gene heatmap saved:\n",
        normalizePath(
            significant_heatmap_file,
            winslash = "/",
            mustWork = FALSE
        ),
        "\n\n"
    )

}


# ============================================================
# FINAL FIGURE QC
# ============================================================

cat(
    "============================================================\n"
)

cat(
    "FIGURE GENERATION QC\n"
)

cat(
    "============================================================\n"
)


figure_check <-
    c(
        volcano =
            file.exists(
                volcano_file
            ),
        full_heatmap =
            file.exists(
                full_heatmap_file
            )
    )


if (
    exists(
        "significant_heatmap_file"
    )
) {

    figure_check[
        "significant_heatmap"
    ] <-
        file.exists(
            significant_heatmap_file
        )

}


print(
    figure_check
)


cat(
    "\nAll generated figure files:\n\n"
)


generated_figures <-
    list.files(
        "results",
        pattern =
            "^Step11.*\\.png$",
        full.names =
            TRUE
    )


if (
    length(
        generated_figures
    ) > 0
) {

    print(
        normalizePath(
            generated_figures,
            winslash =
                "/",
            mustWork =
                FALSE
        )
    )

} else {

    stop(
        "QC FAILED: no Step 11 PNG figures were found in results/."
    )

}


cat(
    "\n============================================================\n"
)

cat(
    "STEP 11 FIGURE GENERATION COMPLETED\n"
)

cat(
    "============================================================\n"
)
