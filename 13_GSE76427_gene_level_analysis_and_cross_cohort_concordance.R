# ============================================================
# STEP 13 — GSE76427 EXTERNAL VALIDATION OF STEP 11
#
# INDIVIDUAL HR/DDR GENE-LEVEL LIPID ASSOCIATIONS
#
# PURPOSE:
# Validate the individual-gene findings from TCGA Step 11
# independently in GSE76427.
#
# IMPORTANT:
# - Uses the SAME expanded HR/DDR genes tested in Step 11
# - Does NOT create a new gene list
# - Uses existing GSE76427 expression from Step 12
# - Uses existing GSE76427 lipid score from Step 12
# - Uses Spearman correlation
# - Uses BH-FDR correction
# - Evaluates cross-cohort directional concordance
# - No [[ syntax is used anywhere
# ============================================================


cat("\n============================================================\n")
cat("STEP 13 — GSE76427 GENE-LEVEL EXTERNAL VALIDATION\n")
cat("VALIDATION OF TCGA STEP 11\n")
cat("============================================================\n\n")


# ============================================================
# 13.1 — CHECK REQUIRED OBJECTS
# ============================================================

required_objects <-
    c(
        "gse12_expression",
        "gse12_lipid_score"
    )


missing_objects <-
    required_objects[
        !vapply(
            required_objects,
            exists,
            logical(1),
            envir = .GlobalEnv
        )
    ]


if (
    length(
        missing_objects
    ) > 0
) {

    stop(
        paste0(
            "The following required Step 12 objects are missing: ",
            paste(
                missing_objects,
                collapse = ", "
            ),
            ". Please run the corrected Step 12 code first."
        )
    )

}


cat(
    "Confirmed: GSE76427 expression object is available.\n"
)

cat(
    "Confirmed: GSE76427 lipid score is available.\n"
)


# ============================================================
# 13.2 — PREPARE GSE76427 EXPRESSION MATRIX
# ============================================================

gse13_expression <-
    as.matrix(
        gse12_expression
    )


storage.mode(
    gse13_expression
) <-
    "numeric"


if (
    ncol(
        gse13_expression
    ) !=
        115
) {

    stop(
        paste0(
            "Expected 115 GSE76427 samples, but found ",
            ncol(
                gse13_expression
            ),
            "."
        )
    )

}


if (
    is.null(
        colnames(
            gse13_expression
        )
    )
) {

    stop(
        "GSE76427 expression matrix has no sample identifiers."
    )

}


# ============================================================
# 13.3 — PREPARE GSE76427 LIPID SCORE
# ============================================================

gse13_lipid_score <-
    as.numeric(
        gse12_lipid_score
    )


names(
    gse13_lipid_score
) <-
    names(
        gse12_lipid_score
    )


if (
    length(
        gse13_lipid_score
    ) !=
        115
) {

    stop(
        paste0(
            "Expected 115 lipid scores, but found ",
            length(
                gse13_lipid_score
            ),
            "."
        )
    )

}


# ============================================================
# 13.4 — VERIFY SAMPLE IDENTIFIER MATCHING
# ============================================================

expression_samples <-
    colnames(
        gse13_expression
    )


lipid_samples <-
    names(
        gse13_lipid_score
    )


if (
    is.null(
        lipid_samples
    )
) {

    stop(
        "GSE76427 lipid score has no sample identifiers."
    )

}


if (
    !setequal(
        expression_samples,
        lipid_samples
    )
) {

    stop(
        "GSE76427 expression samples and lipid-score samples do not match."
    )

}


# Reorder lipid score to expression-matrix order

gse13_lipid_score <-
    gse13_lipid_score[
        match(
            expression_samples,
            names(
                gse13_lipid_score
            )
        )
    ]


names(
    gse13_lipid_score
) <-
    expression_samples


if (
    !identical(
        names(
            gse13_lipid_score
        ),
        colnames(
            gse13_expression
        )
    )
) {

    stop(
        "Final sample-order QC failed."
    )

}


cat(
    "Confirmed: GSE76427 expression and lipid score contain the same 115 samples in the same order.\n"
)


# ============================================================
# 13.5 — LOCATE THE STEP 11 TCGA RESULTS FILE
# ============================================================

step11_correlation_file <-
    file.path(
        "results",
        "Step11_2.5_expanded_HR_gene_lipid_correlation_results.csv"
    )


if (
    !file.exists(
        step11_correlation_file
    )
) {

    stop(
        paste0(
            "The Step 11 correlation-results file was not found:\n",
            normalizePath(
                step11_correlation_file,
                winslash = "/",
                mustWork = FALSE
            ),
            "\nThis file is required because Step 13 must use ",
            "the exact same expanded HR/DDR gene list as Step 11."
        )
    )

}


cat(
    "\nFound Step 11 correlation-results file:\n",
    normalizePath(
        step11_correlation_file,
        winslash = "/",
        mustWork = FALSE
    ),
    "\n"
)


# ============================================================
# 13.6 — READ TCGA STEP 11 RESULTS
# ============================================================

tcga_step11_results <-
    read.csv(
        step11_correlation_file,
        stringsAsFactors =
            FALSE,
        check.names =
            FALSE
    )


cat(
    "\nTCGA Step 11 results dimensions: ",
    nrow(
        tcga_step11_results
    ),
    " rows x ",
    ncol(
        tcga_step11_results
    ),
    " columns\n",
    sep = ""
)


cat(
    "\nColumns in Step 11 results:\n\n"
)

print(
    names(
        tcga_step11_results
    )
)


# ============================================================
# 13.7 — IDENTIFY GENE COLUMN
# ============================================================

gene_column_candidates <-
    names(
        tcga_step11_results
    )[
        grepl(
            "gene",
            names(
                tcga_step11_results
            ),
            ignore.case =
                TRUE
        )
    ]


if (
    length(
        gene_column_candidates
    ) == 0
) {

    stop(
        "Could not identify the gene column in the Step 11 results file."
    )

}


gene_column_name <-
    gene_column_candidates[
        1
    ]


cat(
    "\nGene column identified as: ",
    gene_column_name,
    "\n",
    sep = ""
)


# ============================================================
# 13.8 — EXTRACT EXACT STEP 11 GENE LIST
# ============================================================

expanded_hr_genes <-
    unique(
        toupper(
            trimws(
                as.character(
                    tcga_step11_results[
                        ,
                        gene_column_name
                    ]
                )
            )
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
    "\nNumber of genes retrieved from Step 11: ",
    length(
        expanded_hr_genes
    ),
    "\n",
    sep = ""
)


if (
    length(
        expanded_hr_genes
    ) !=
        47
) {

    warning(
        paste0(
            "Step 11 results contain ",
            length(
                expanded_hr_genes
            ),
            " genes rather than the expected 47. ",
            "The code will continue using the exact genes ",
            "present in the Step 11 results file."
        )
    )

}


# ============================================================
# 13.9 — STANDARDIZE GSE76427 GENE SYMBOLS
# ============================================================

rownames(
    gse13_expression
) <-
    toupper(
        trimws(
            rownames(
                gse13_expression
            )
        )
    )


valid_rows <-
    !is.na(
        rownames(
            gse13_expression
        )
    ) &
    rownames(
        gse13_expression
    ) != ""


gse13_expression <-
    gse13_expression[
        valid_rows,
        ,
        drop = FALSE
    ]


# ============================================================
# 13.10 — COLLAPSE DUPLICATE GENE SYMBOLS
# ============================================================

if (
    any(
        duplicated(
            rownames(
                gse13_expression
            )
        )
    )
) {

    unique_genes <-
        unique(
            rownames(
                gse13_expression
            )
        )


    collapsed_expression <-
        matrix(
            NA_real_,
            nrow =
                length(
                    unique_genes
                ),
            ncol =
                ncol(
                    gse13_expression
                ),
            dimnames =
                list(
                    unique_genes,
                    colnames(
                        gse13_expression
                    )
                )
        )


    for (
        i in seq_along(
            unique_genes
        )
    ) {

        current_gene <-
            unique_genes[
                i
            ]


        current_rows <-
            which(
                rownames(
                    gse13_expression
                ) ==
                    current_gene
            )


        if (
            length(
                current_rows
            ) ==
                1
        ) {

            collapsed_expression[
                i,
            ] <-
                gse13_expression[
                    current_rows,
                ]

        } else {

            collapsed_expression[
                i,
            ] <-
                apply(
                    gse13_expression[
                        current_rows,
                        ,
                        drop = FALSE
                    ],
                    2,
                    median,
                    na.rm =
                        TRUE
                )

        }

    }


    gse13_expression <-
        collapsed_expression

}


# ============================================================
# 13.11 — DETERMINE GENE COVERAGE
# ============================================================

genes_present <-
    intersect(
        expanded_hr_genes,
        rownames(
            gse13_expression
        )
    )


genes_missing <-
    setdiff(
        expanded_hr_genes,
        rownames(
            gse13_expression
        )
    )


cat(
    "\n============================================================\n"
)

cat(
    "GSE76427 HR/DDR GENE COVERAGE\n"
)

cat(
    "============================================================\n\n"
)

cat(
    "Genes in Step 11 list: ",
    length(
        expanded_hr_genes
    ),
    "\n",
    sep = ""
)

cat(
    "Genes present in GSE76427: ",
    length(
        genes_present
    ),
    "\n",
    sep = ""
)

cat(
    "Genes absent from GSE76427: ",
    length(
        genes_missing
    ),
    "\n",
    sep = ""
)


if (
    length(
        genes_missing
    ) > 0
) {

    cat(
        "\nGenes absent from GSE76427:\n\n"
    )

    print(
        genes_missing
    )

}


# ============================================================
# 13.12 — CALCULATE GSE76427 GENE-LEVEL CORRELATIONS
# ============================================================

gse13_results <-
    data.frame(
        Gene =
            expanded_hr_genes,
        Present_in_GSE76427 =
            expanded_hr_genes %in%
            genes_present,
        N =
            NA_integer_,
        Spearman_rho =
            NA_real_,
        P_value =
            NA_real_,
        stringsAsFactors =
            FALSE
    )


for (
    i in seq_along(
        expanded_hr_genes
    )
) {

    current_gene <-
        expanded_hr_genes[
            i
        ]


    if (
        !current_gene %in%
        rownames(
            gse13_expression
        )
    ) {

        next

    }


    current_expression <-
        as.numeric(
            gse13_expression[
                current_gene,
                ]
        )


    valid <-
        is.finite(
            current_expression
        ) &
        is.finite(
            gse13_lipid_score
        )


    if (
        sum(
            valid
        ) <
            10
    ) {

        next

    }


    if (
        length(
            unique(
                current_expression[
                    valid
                ]
            )
        ) <
            2
    ) {

        next

    }


    current_test <-
        suppressWarnings(
            cor.test(
                gse13_lipid_score[
                    valid
                ],
                current_expression[
                    valid
                ],
                method =
                    "spearman",
                exact =
                    FALSE
            )
        )


    gse13_results$N[
        i
    ] <-
        sum(
            valid
        )


    gse13_results$Spearman_rho[
        i
    ] <-
        as.numeric(
            current_test$estimate
        )


    gse13_results$P_value[
        i
    ] <-
        current_test$p.value

}


# ============================================================
# 13.13 — BH-FDR CORRECTION
# ============================================================

gse13_results$BH_FDR <-
    p.adjust(
        gse13_results$P_value,
        method =
            "BH"
    )


# ============================================================
# 13.14 — SIGNIFICANCE AND DIRECTION
# ============================================================

gse13_results$GSE76427_direction <-
    ifelse(
        is.na(
            gse13_results$Spearman_rho
        ),
        "Not_tested",
        ifelse(
            gse13_results$Spearman_rho < 0,
            "Negative",
            ifelse(
                gse13_results$Spearman_rho > 0,
                "Positive",
                "Zero"
            )
        )
    )


gse13_results$GSE76427_significant <-
    !is.na(
        gse13_results$BH_FDR
    ) &
    gse13_results$BH_FDR <
        0.05


# ============================================================
# 13.15 — IDENTIFY TCGA RHO COLUMN
# ============================================================

rho_column_candidates <-
    names(
        tcga_step11_results
    )[
        grepl(
            "spearman|rho|correlation|estimate",
            names(
                tcga_step11_results
            ),
            ignore.case =
                TRUE
        )
    ]


rho_column_candidates <-
    rho_column_candidates[
        !grepl(
            "p.value|p_value|fdr|adjust|BH",
            rho_column_candidates,
            ignore.case =
                TRUE
        )
    ]


if (
    length(
        rho_column_candidates
    ) == 0
) {

    stop(
        "Could not identify the TCGA Spearman rho column in the Step 11 results file."
    )

}


tcga_rho_column_name <-
    rho_column_candidates[
        1
    ]


cat(
    "\nTCGA rho column identified as: ",
    tcga_rho_column_name,
    "\n",
    sep = ""
)


# ============================================================
# 13.16 — IDENTIFY TCGA P/FDR COLUMNS
# ============================================================

tcga_fdr_candidates <-
    names(
        tcga_step11_results
    )[
        grepl(
            "BH|FDR|adjust",
            names(
                tcga_step11_results
            ),
            ignore.case =
                TRUE
        )
    ]


if (
    length(
        tcga_fdr_candidates
    ) > 0
) {

    tcga_fdr_column_name <-
        tcga_fdr_candidates[
            1
        ]

} else {

    tcga_fdr_column_name <-
        NA_character_

}


# ============================================================
# 13.17 — CREATE TCGA VALIDATION TABLE
# ============================================================

tcga_gene_summary <-
    data.frame(
        Gene =
            toupper(
                trimws(
                    as.character(
                        tcga_step11_results[
                            ,
                            gene_column_name
                        ]
                    )
                )
            ),
        TCGA_Spearman_rho =
            suppressWarnings(
                as.numeric(
                    tcga_step11_results[
                        ,
                        tcga_rho_column_name
                    ]
                )
            ),
        stringsAsFactors =
            FALSE
    )


if (
    !is.na(
        tcga_fdr_column_name
    )
) {

    tcga_gene_summary$TCGA_BH_FDR <-
        suppressWarnings(
            as.numeric(
                tcga_step11_results[
                    ,
                    tcga_fdr_column_name
                ]
            )
        )

} else {

    tcga_gene_summary$TCGA_BH_FDR <-
        NA_real_

}


tcga_gene_summary <-
    tcga_gene_summary[
        tcga_gene_summary$Gene %in%
            expanded_hr_genes,
        ,
        drop = FALSE
    ]


# Remove duplicate gene entries if present

tcga_gene_summary <-
    tcga_gene_summary[
        !duplicated(
            tcga_gene_summary$Gene
        ),
        ,
        drop = FALSE
    ]


# ============================================================
# 13.18 — MERGE TCGA AND GSE76427 RESULTS
# ============================================================

step13_comparison <-
    merge(
        tcga_gene_summary,
        gse13_results,
        by =
            "Gene",
        all =
            TRUE,
        sort =
            FALSE
    )


# ============================================================
# 13.19 — DETERMINE DIRECTIONAL CONCORDANCE
# ============================================================

step13_comparison$TCGA_direction <-
    ifelse(
        is.na(
            step13_comparison$TCGA_Spearman_rho
        ),
        "Not_available",
        ifelse(
            step13_comparison$TCGA_Spearman_rho < 0,
            "Negative",
            ifelse(
                step13_comparison$TCGA_Spearman_rho > 0,
                "Positive",
                "Zero"
            )
        )
    )


step13_comparison$Direction_concordant <-
    !is.na(
        step13_comparison$TCGA_Spearman_rho
    ) &
    !is.na(
        step13_comparison$Spearman_rho
    ) &
    sign(
        step13_comparison$TCGA_Spearman_rho
    ) ==
    sign(
        step13_comparison$Spearman_rho
    )


step13_comparison$Replicated_significant <-
    !is.na(
        step13_comparison$TCGA_BH_FDR
    ) &
    step13_comparison$TCGA_BH_FDR <
        0.05 &
    !is.na(
        step13_comparison$BH_FDR
    ) &
    step13_comparison$BH_FDR <
        0.05 &
    step13_comparison$Direction_concordant


step13_comparison$TCGA_significant_GSE_direction_concordant <-
    !is.na(
        step13_comparison$TCGA_BH_FDR
    ) &
    step13_comparison$TCGA_BH_FDR <
        0.05 &
    step13_comparison$Direction_concordant


# ============================================================
# 13.20 — CALCULATE CROSS-COHORT CONCORDANCE
# ============================================================

concordance_data <-
    step13_comparison[
        is.finite(
            step13_comparison$TCGA_Spearman_rho
        ) &
        is.finite(
            step13_comparison$Spearman_rho
        ),
        ,
        drop = FALSE
    ]


if (
    nrow(
        concordance_data
    ) >=
        3
) {

    concordance_test <-
        suppressWarnings(
            cor.test(
                concordance_data$TCGA_Spearman_rho,
                concordance_data$Spearman_rho,
                method =
                    "spearman",
                exact =
                    FALSE
            )
        )


    concordance_rho <-
        as.numeric(
            concordance_test$estimate
        )


    concordance_p <-
        concordance_test$p.value

} else {

    concordance_rho <-
        NA_real_


    concordance_p <-
        NA_real_

}


# ============================================================
# 13.21 — SUMMARY COUNTS
# ============================================================

n_tested_gse <-
    sum(
        !is.na(
            gse13_results$Spearman_rho
        )
    )


n_gse_negative <-
    sum(
        gse13_results$GSE76427_direction ==
            "Negative",
        na.rm =
            TRUE
    )


n_gse_positive <-
    sum(
        gse13_results$GSE76427_direction ==
            "Positive",
        na.rm =
            TRUE
    )


n_gse_significant <-
    sum(
        gse13_results$GSE76427_significant,
        na.rm =
            TRUE
    )


n_tcga_sig_gse_sig <-
    sum(
        step13_comparison$Replicated_significant,
        na.rm =
            TRUE
    )


n_direction_concordant <-
    sum(
        step13_comparison$Direction_concordant,
        na.rm =
            TRUE
    )


# ============================================================
# 13.22 — PRINT VALIDATION RESULTS
# ============================================================

cat(
    "\n============================================================\n"
)

cat(
    "STEP 13 — GSE76427 GENE-LEVEL VALIDATION RESULTS\n"
)

cat(
    "============================================================\n\n"
)

cat(
    "Expanded HR/DDR genes from TCGA Step 11: ",
    length(
        expanded_hr_genes
    ),
    "\n",
    sep = ""
)

cat(
    "Genes tested successfully in GSE76427: ",
    n_tested_gse,
    "\n",
    sep = ""
)

cat(
    "Negative correlations in GSE76427: ",
    n_gse_negative,
    "\n",
    sep = ""
)

cat(
    "Positive correlations in GSE76427: ",
    n_gse_positive,
    "\n",
    sep = ""
)

cat(
    "Significant GSE76427 genes, BH-FDR <0.05: ",
    n_gse_significant,
    "\n",
    sep = ""
)

cat(
    "Directionally concordant genes: ",
    n_direction_concordant,
    "\n",
    sep = ""
)

cat(
    "Significant in both cohorts with concordant direction: ",
    n_tcga_sig_gse_sig,
    "\n",
    sep = ""
)

cat(
    "\nCross-cohort Spearman correlation of gene-level rho values:\n"
)

cat(
    "rho = ",
    round(
        concordance_rho,
        4
    ),
    "\n",
    sep = ""
)

cat(
    "P = ",
    format.pval(
        concordance_p,
        digits =
            4
    ),
    "\n",
    sep = ""
)


cat(
    "\nGenes significantly replicated in both cohorts:\n\n"
)


replicated_genes <-
    step13_comparison[
        step13_comparison$Replicated_significant,
        ,
        drop = FALSE
    ]


if (
    nrow(
        replicated_genes
    ) > 0
) {

    replicated_genes <-
        replicated_genes[
            order(
                replicated_genes$BH_FDR
            ),
            ,
            drop = FALSE
        ]


    print(
        replicated_genes[
            ,
            c(
                "Gene",
                "TCGA_Spearman_rho",
                "TCGA_BH_FDR",
                "Spearman_rho",
                "BH_FDR"
            ),
            drop = FALSE
        ],
        row.names =
            FALSE
    )

} else {

    cat(
        "No genes met the strict significant-in-both-cohorts criterion.\n"
    )

}


# ============================================================
# 13.23 — PREPARE VOLCANO DATA
# ============================================================

gse13_volcano <-
    gse13_results[
        is.finite(
            gse13_results$Spearman_rho
        ) &
        is.finite(
            gse13_results$BH_FDR
        ),
        ,
        drop = FALSE
    ]


gse13_volcano$minus_log10_FDR <-
    -log10(
        pmax(
            gse13_volcano$BH_FDR,
            .Machine$double.xmin
        )
    )


gse13_volcano$Significance <-
    ifelse(
        gse13_volcano$BH_FDR <
            0.05 &
        gse13_volcano$Spearman_rho <
            0,
        "Significant negative",
        ifelse(
            gse13_volcano$BH_FDR <
                0.05 &
            gse13_volcano$Spearman_rho >
                0,
            "Significant positive",
            "Not significant"
        )
    )


# ============================================================
# 13.24 — CREATE RESULTS DIRECTORY
# ============================================================

if (
    !dir.exists(
        "results"
    )
) {

    dir.create(
        "results",
        recursive =
            TRUE
    )

}


# ============================================================
# 13.25 — GSE76427 VOLCANO PLOT
# ============================================================

if (
    !requireNamespace(
        "ggplot2",
        quietly =
            TRUE
    )
) {

    stop(
        "ggplot2 is required for Step 13 plots."
    )

}


library(
    ggplot2
)


gse13_volcano_plot <-
    ggplot(
        gse13_volcano,
        aes(
            x =
                Spearman_rho,
            y =
                minus_log10_FDR
        )
    ) +
    geom_point(
        size =
            3,
        alpha =
            0.85
    ) +
    geom_vline(
        xintercept =
            0,
        linetype =
            "dashed"
    ) +
    geom_hline(
        yintercept =
            -log10(
                0.05
            ),
        linetype =
            "dashed"
    ) +
    ggrepel::geom_text_repel(
        data =
            gse13_volcano[
                gse13_volcano$BH_FDR <
                    0.05,
                ,
                drop = FALSE
            ],
        aes(
            label =
                Gene
        ),
        size =
            3.5,
        max.overlaps =
            25
    ) +
    labs(
        title =
            "GSE76427: Individual HR/DDR Gene Associations with Lipid Metabolism",
        subtitle =
            paste0(
                "n = 115; ",
                n_gse_significant,
                " genes significant at BH-FDR < 0.05"
            ),
        x =
            "Spearman correlation with lipid score",
        y =
            expression(
                -log[10](
                    "BH-FDR"
                )
            )
    ) +
    theme_classic(
        base_size =
            14
    )


ggsave(
    filename =
        file.path(
            "results",
            "Step13_GSE76427_expanded_HR_gene_lipid_volcano.png"
        ),
    plot =
        gse13_volcano_plot,
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


# ============================================================
# 13.26 — GENE CORRELATION HEATMAP
# ============================================================

heatmap_data <-
    gse13_results[
        is.finite(
            gse13_results$Spearman_rho
        ),
        ,
        drop = FALSE
    ]


heatmap_data <-
    heatmap_data[
        order(
            heatmap_data$Spearman_rho
        ),
        ,
        drop = FALSE
    ]


heatmap_data$Gene <-
    factor(
        heatmap_data$Gene,
        levels =
            heatmap_data$Gene
    )


gse13_heatmap <-
    ggplot(
        heatmap_data,
        aes(
            x =
                "GSE76427",
            y =
                Gene,
            fill =
                Spearman_rho
        )
    ) +
    geom_tile(
        linewidth =
            0.3
    ) +
    geom_text(
        aes(
            label =
                sprintf(
                    "%.2f",
                    Spearman_rho
                )
        ),
        size =
            3
    ) +
    scale_fill_gradient2(
        midpoint =
            0
    ) +
    labs(
        title =
            "GSE76427: HR/DDR Gene–Lipid Correlations",
        subtitle =
            "Spearman rho values",
        x =
            NULL,
        y =
            NULL,
        fill =
            "Spearman rho"
    ) +
    theme_classic(
        base_size =
            13
    ) +
    theme(
        axis.text.x =
            element_text(
                angle =
                    0,
                hjust =
                    0.5
            )
    )


ggsave(
    filename =
        file.path(
            "results",
            "Step13_GSE76427_expanded_HR_gene_lipid_heatmap.png"
        ),
    plot =
        gse13_heatmap,
    width =
        7,
    height =
        11,
    units =
        "in",
    dpi =
        300,
    bg =
        "white"
)


# ============================================================
# 13.27 — TCGA VS GSE76427 CROSS-COHORT CONCORDANCE PLOT
# ============================================================

concordance_plot_data <-
    step13_comparison[
        is.finite(
            step13_comparison$TCGA_Spearman_rho
        ) &
        is.finite(
            step13_comparison$Spearman_rho
        ),
        ,
        drop = FALSE
    ]


concordance_plot <-
    ggplot(
        concordance_plot_data,
        aes(
            x =
                TCGA_Spearman_rho,
            y =
                Spearman_rho
        )
    ) +
    geom_point(
        size =
            3,
        alpha =
            0.85
    ) +
    geom_smooth(
        method =
            "lm",
        se =
            TRUE
    ) +
    geom_hline(
        yintercept =
            0,
        linetype =
            "dashed"
    ) +
    geom_vline(
        xintercept =
            0,
        linetype =
            "dashed"
    ) +
    ggrepel::geom_text_repel(
        data =
            concordance_plot_data[
                concordance_plot_data$Replicated_significant,
                ,
                drop = FALSE
            ],
        aes(
            label =
                Gene
        ),
        size =
            3.5,
        max.overlaps =
            30
    ) +
    labs(
        title =
            "Cross-Cohort Concordance of HR/DDR Gene–Lipid Associations",
        subtitle =
            paste0(
                "TCGA vs GSE76427; Spearman rho between effect sizes = ",
                round(
                    concordance_rho,
                    3
                ),
                "; P = ",
                format.pval(
                    concordance_p,
                    digits =
                        3
                )
            ),
        x =
            "TCGA Step 11 Spearman rho",
        y =
            "GSE76427 Spearman rho"
    ) +
    theme_classic(
        base_size =
            14
    )


ggsave(
    filename =
        file.path(
            "results",
            "Step13_TCGA_vs_GSE76427_HR_gene_lipid_concordance.png"
        ),
    plot =
        concordance_plot,
    width =
        9,
    height =
        8,
    units =
        "in",
    dpi =
        300,
    bg =
        "white"
)


# ============================================================
# 13.28 — SAVE ALL GENE-LEVEL RESULTS
# ============================================================

write.csv(
    gse13_results,
    file =
        file.path(
            "results",
            "Step13_GSE76427_expanded_HR_gene_lipid_correlation_results.csv"
        ),
    row.names =
        FALSE
)


write.csv(
    step13_comparison,
    file =
        file.path(
            "results",
            "Step13_TCGA_vs_GSE76427_HR_gene_lipid_comparison.csv"
        ),
    row.names =
        FALSE
)


write.csv(
    replicated_genes,
    file =
        file.path(
            "results",
            "Step13_replicated_HR_genes_TCGA_GSE76427.csv"
        ),
    row.names =
        FALSE
)


write.csv(
    data.frame(
        Metric =
            c(
                "Total Step 11 genes",
                "Genes tested in GSE76427",
                "Negative GSE76427 correlations",
                "Positive GSE76427 correlations",
                "Significant GSE76427 genes",
                "Directionally concordant genes",
                "Significant in both cohorts",
                "Cross-cohort Spearman rho",
                "Cross-cohort P value"
            ),
        Value =
            c(
                length(
                    expanded_hr_genes
                ),
                n_tested_gse,
                n_gse_negative,
                n_gse_positive,
                n_gse_significant,
                n_direction_concordant,
                n_tcga_sig_gse_sig,
                concordance_rho,
                concordance_p
            ),
        stringsAsFactors =
            FALSE
    ),
    file =
        file.path(
            "results",
            "Step13_GSE76427_gene_validation_summary.csv"
        ),
    row.names =
        FALSE
)


# ============================================================
# 13.29 — SAVE RDATA
# ============================================================

save(
    gse13_expression,
    gse13_lipid_score,
    expanded_hr_genes,
    genes_present,
    genes_missing,
    gse13_results,
    step13_comparison,
    replicated_genes,
    concordance_rho,
    concordance_p,
    file =
        file.path(
            "results",
            "Step13_GSE76427_HR_gene_validation.RData"
        )
)


# ============================================================
# 13.30 — FINAL QC
# ============================================================

if (
    ncol(
        gse13_expression
    ) !=
        115
) {

    stop(
        "FINAL QC FAILED: GSE76427 expression does not contain 115 samples."
    )

}


if (
    length(
        gse13_lipid_score
    ) !=
        115
) {

    stop(
        "FINAL QC FAILED: GSE76427 lipid score does not contain 115 samples."
    )

}


if (
    !identical(
        names(
            gse13_lipid_score
        ),
        colnames(
            gse13_expression
        )
    )
) {

    stop(
        "FINAL QC FAILED: lipid score sample order does not match expression matrix."
    )

}


if (
    nrow(
        gse13_results
    ) !=
        length(
            expanded_hr_genes
        )
) {

    stop(
        "FINAL QC FAILED: gene-level result count does not match the Step 11 gene list."
    )

}


# ============================================================
# 13.31 — FINAL OUTPUT LIST
# ============================================================

step13_files <-
    list.files(
        "results",
        pattern =
            "^Step13_",
        full.names =
            TRUE
    )


cat(
    "\n============================================================\n"
)

cat(
    "STEP 13 COMPLETED SUCCESSFULLY\n"
)

cat(
    "============================================================\n\n"
)

cat(
    "GSE76427 validation samples: 115\n"
)

cat(
    "Step 11 expanded HR/DDR genes: ",
    length(
        expanded_hr_genes
    ),
    "\n",
    sep = ""
)

cat(
    "Genes tested in GSE76427: ",
    n_tested_gse,
    "\n",
    sep = ""
)

cat(
    "Significant GSE76427 genes: ",
    n_gse_significant,
    "\n",
    sep = ""
)

cat(
    "Directionally concordant genes: ",
    n_direction_concordant,
    "\n",
    sep = ""
)

cat(
    "Significant in both cohorts: ",
    n_tcga_sig_gse_sig,
    "\n",
    sep = ""
)

cat(
    "\nCross-cohort rho concordance: ",
    round(
        concordance_rho,
        4
    ),
    "\n",
    sep = ""
)

cat(
    "Cross-cohort P value: ",
    format.pval(
        concordance_p,
        digits =
            4
    ),
    "\n",
    sep = ""
)

cat(
    "\nStep 13 output files:\n\n"
)

print(
    normalizePath(
        step13_files,
        winslash =
            "/",
        mustWork =
            FALSE
    )
)

cat(
    "\n============================================================\n"
)
