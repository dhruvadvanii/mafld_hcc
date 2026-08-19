# ============================================================
# STEP 14C — GSE14520 PROBE-TO-GENE SUMMARIZATION
# FINAL ROBUST VERSION
# ============================================================
#
# Probe selection:
#   1. Highest MAD across HCC tumor samples
#   2. Highest median expression if tied
#   3. Lexicographically first probe ID if still tied
#
# No biological outcome is used.
# GPL571 and GPL3921 remain separate.
#
# ============================================================


# ------------------------------------------------------------
# 14C.0 — CHECK REQUIRED OBJECTS
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP 14C — PROBE-TO-GENE SUMMARIZATION\n")
cat("============================================================\n\n")


required_objects_14c <- c(
    "gse14520_gpl571_tumor_expression",
    "gse14520_gpl3921_tumor_expression",
    "gpl571_candidate_probe_map",
    "gpl3921_candidate_probe_map",
    "gse14520_candidate_genes"
)

missing_objects_14c <- required_objects_14c[
    !sapply(
        required_objects_14c,
        exists
    )
]

if (length(missing_objects_14c) > 0) {

    stop(
        paste0(
            "\nERROR: Required Step 14B objects are missing:\n\n",
            paste(
                missing_objects_14c,
                collapse = "\n"
            ),
            "\n\nPlease run Step 14B first.\n"
        )
    )
}


# ------------------------------------------------------------
# 14C.1 — VERIFY FROZEN 18-GENE PANEL
# ------------------------------------------------------------

if (
    length(gse14520_candidate_genes) != 18
) {

    stop(
        "\nERROR: Frozen candidate panel does not contain exactly 18 genes.\n"
    )
}


# ============================================================
# 14C.2 — FUNCTION: SELECT REPRESENTATIVE PROBES
# ============================================================

select_representative_probes <- function(
    expression_matrix,
    probe_map,
    candidate_genes
) {

    probe_map <- probe_map[
        probe_map$probe_id %in%
            rownames(expression_matrix),
        ,
        drop = FALSE
    ]

    probe_map <- probe_map[
        !is.na(probe_map$gene_symbol) &
            nzchar(probe_map$gene_symbol),
        ,
        drop = FALSE
    ]


    output <- data.frame(
        gene = character(0),
        probe_id = character(0),
        median_expression = numeric(0),
        MAD = numeric(0),
        SD = numeric(0),
        stringsAsFactors = FALSE
    )


    for (
        current_gene in candidate_genes
    ) {

        current_probes <- probe_map[
            probe_map$gene_symbol ==
                current_gene,
            ,
            drop = FALSE
        ]


        if (
            nrow(current_probes) == 0
        ) {
            next
        }


        statistics <- data.frame(
            gene = character(0),
            probe_id = character(0),
            median_expression = numeric(0),
            MAD = numeric(0),
            SD = numeric(0),
            stringsAsFactors = FALSE
        )


        for (
            current_probe in current_probes$probe_id
        ) {

            values <- expression_matrix[
                current_probe,
                ,
                drop = TRUE
            ]

            values <- as.numeric(
                values
            )

            values <- values[
                is.finite(values)
            ]


            if (
                length(values) == 0
            ) {
                next
            }


            statistics <- rbind(
                statistics,
                data.frame(
                    gene = current_gene,
                    probe_id = current_probe,
                    median_expression = median(
                        values,
                        na.rm = TRUE
                    ),
                    MAD = mad(
                        values,
                        constant = 1,
                        na.rm = TRUE
                    ),
                    SD = sd(
                        values,
                        na.rm = TRUE
                    ),
                    stringsAsFactors = FALSE
                )
            )
        }


        if (
            nrow(statistics) == 0
        ) {
            next
        }


        # ----------------------------------------------------
        # Predefined probe-selection hierarchy
        # ----------------------------------------------------

        statistics <- statistics[
            order(
                -statistics$MAD,
                -statistics$median_expression,
                statistics$probe_id
            ),
            ,
            drop = FALSE
        ]


        selected <- statistics[
            1,
            ,
            drop = FALSE
        ]


        output <- rbind(
            output,
            selected
        )
    }


    rownames(output) <- NULL

    return(output)
}


# ============================================================
# 14C.3 — SELECT GPL571 PROBES
# ============================================================

cat("\n============================================================\n")
cat("GPL571 — REPRESENTATIVE PROBE SELECTION\n")
cat("============================================================\n\n")


gpl571_selected_probes <- select_representative_probes(
    expression_matrix =
        gse14520_gpl571_tumor_expression,
    probe_map =
        gpl571_candidate_probe_map,
    candidate_genes =
        gse14520_candidate_genes
)


# ============================================================
# 14C.4 — SELECT GPL3921 PROBES
# ============================================================

cat("\n============================================================\n")
cat("GPL3921 — REPRESENTATIVE PROBE SELECTION\n")
cat("============================================================\n\n")


gpl3921_selected_probes <- select_representative_probes(
    expression_matrix =
        gse14520_gpl3921_tumor_expression,
    probe_map =
        gpl3921_candidate_probe_map,
    candidate_genes =
        gse14520_candidate_genes
)


# ============================================================
# 14C.5 — DISPLAY SELECTED PROBES
# ============================================================

cat("\n============================================================\n")
cat("SELECTED PROBES — GPL571\n")
cat("============================================================\n\n")

print(
    gpl571_selected_probes
)


cat("\n============================================================\n")
cat("SELECTED PROBES — GPL3921\n")
cat("============================================================\n\n")

print(
    gpl3921_selected_probes
)


# ============================================================
# 14C.6 — DETERMINE COMMON GENE PANEL
# ============================================================

gpl571_selected_genes <- gpl571_selected_probes$gene

gpl3921_selected_genes <- gpl3921_selected_probes$gene


gpl571_missing_genes <- gse14520_candidate_genes[
    !gse14520_candidate_genes %in%
        gpl571_selected_genes
]

gpl3921_missing_genes <- gse14520_candidate_genes[
    !gse14520_candidate_genes %in%
        gpl3921_selected_genes
]


gse14520_common_genes <- intersect(
    gpl571_selected_genes,
    gpl3921_selected_genes
)


cat("\n============================================================\n")
cat("GENE COVERAGE\n")
cat("============================================================\n\n")

cat(
    "GPL571:",
    length(gpl571_selected_genes),
    "/18 genes\n"
)

cat(
    "GPL3921:",
    length(gpl3921_selected_genes),
    "/18 genes\n"
)

cat(
    "Common:",
    length(gse14520_common_genes),
    "/18 genes\n\n"
)


cat("GPL571 missing:\n")
print(
    gpl571_missing_genes
)

cat("\nGPL3921 missing:\n")
print(
    gpl3921_missing_genes
)

cat("\nCommon genes:\n")
print(
    gse14520_common_genes
)


# ============================================================
# 14C.7 — FUNCTION TO BUILD GENE-LEVEL MATRIX
# ============================================================
#
# IMPORTANT:
#   A numeric matrix is preallocated.
#   Each row is filled using a complete expression vector.
#
# There is NO list assignment.
# There is NO probe ID used as a column index.
# There is NO [[ ]] indexing.
#
# ============================================================

build_gene_matrix <- function(
    expression_matrix,
    selected_probes
) {

    n_genes <- nrow(
        selected_probes
    )

    n_samples <- ncol(
        expression_matrix
    )


    if (
        n_genes == 0
    ) {

        empty_matrix <- matrix(
            numeric(0),
            nrow = 0,
            ncol = n_samples
        )

        colnames(
            empty_matrix
        ) <- colnames(
            expression_matrix
        )

        return(
            empty_matrix
        )
    }


    # --------------------------------------------------------
    # Preallocate exact dimensions
    # --------------------------------------------------------

    gene_matrix <- matrix(
        NA_real_,
        nrow = n_genes,
        ncol = n_samples
    )


    # --------------------------------------------------------
    # Copy sample names
    # --------------------------------------------------------

    colnames(
        gene_matrix
    ) <- colnames(
        expression_matrix
    )


    # --------------------------------------------------------
    # Copy gene names
    # --------------------------------------------------------

    rownames(
        gene_matrix
    ) <- selected_probes$gene


    # --------------------------------------------------------
    # Fill each row
    # --------------------------------------------------------

    for (
        i in seq_len(n_genes)
    ) {

        current_probe <- selected_probes$probe_id[
            i
        ]


        current_values <- expression_matrix[
            current_probe,
            ,
            drop = TRUE
        ]


        current_values <- as.numeric(
            current_values
        )


        if (
            length(current_values) != n_samples
        ) {

            stop(
                paste0(
                    "\nERROR: Expression length mismatch for probe ",
                    current_probe,
                    ". Expected ",
                    n_samples,
                    " values but found ",
                    length(current_values),
                    "."
                )
            )
        }


        # ----------------------------------------------------
        # SAFE ROW ASSIGNMENT
        # ----------------------------------------------------

        gene_matrix[
            i,
            seq_len(n_samples)
        ] <- current_values
    }


    return(
        gene_matrix
    )
}


# ============================================================
# 14C.8 — BUILD GPL571 GENE MATRIX
# ============================================================

cat("\n============================================================\n")
cat("BUILDING GPL571 GENE MATRIX\n")
cat("============================================================\n\n")


gpl571_gene_expression <- build_gene_matrix(
    expression_matrix =
        gse14520_gpl571_tumor_expression,
    selected_probes =
        gpl571_selected_probes
)


# ============================================================
# 14C.9 — BUILD GPL3921 GENE MATRIX
# ============================================================

cat("\n============================================================\n")
cat("BUILDING GPL3921 GENE MATRIX\n")
cat("============================================================\n\n")


gpl3921_gene_expression <- build_gene_matrix(
    expression_matrix =
        gse14520_gpl3921_tumor_expression,
    selected_probes =
        gpl3921_selected_probes
)


# ============================================================
# 14C.10 — VERIFY DIMENSIONS
# ============================================================

cat("\n============================================================\n")
cat("GENE-LEVEL MATRIX DIMENSIONS\n")
cat("============================================================\n\n")

cat(
    "GPL571:",
    nrow(gpl571_gene_expression),
    "genes x",
    ncol(gpl571_gene_expression),
    "tumors\n"
)

cat(
    "GPL3921:",
    nrow(gpl3921_gene_expression),
    "genes x",
    ncol(gpl3921_gene_expression),
    "tumors\n\n"
)


# ============================================================
# 14C.11 — VERIFY SAMPLE IDENTITIES
# ============================================================

gpl571_sample_check <- identical(
    colnames(gpl571_gene_expression),
    colnames(gse14520_gpl571_tumor_expression)
)


gpl3921_sample_check <- identical(
    colnames(gpl3921_gene_expression),
    colnames(gse14520_gpl3921_tumor_expression)
)


cat(
    "GPL571 sample IDs preserved:",
    gpl571_sample_check,
    "\n"
)

cat(
    "GPL3921 sample IDs preserved:",
    gpl3921_sample_check,
    "\n\n"
)


if (
    !gpl571_sample_check
) {

    stop(
        "\nERROR: GPL571 sample IDs were not preserved.\n"
    )
}


if (
    !gpl3921_sample_check
) {

    stop(
        "\nERROR: GPL3921 sample IDs were not preserved.\n"
    )
}


# ============================================================
# 14C.12 — CHECK GENE NAMES
# ============================================================

cat("GPL571 gene names:\n")

print(
    rownames(
        gpl571_gene_expression
    )
)


cat("\nGPL3921 gene names:\n")

print(
    rownames(
        gpl3921_gene_expression
    )
)


# ============================================================
# 14C.13 — CHECK NON-FINITE VALUES
# ============================================================

gpl571_nonfinite <- sum(
    !is.finite(
        gpl571_gene_expression
    )
)

gpl3921_nonfinite <- sum(
    !is.finite(
        gpl3921_gene_expression
    )
)


cat("\n============================================================\n")
cat("NON-FINITE VALUE CHECK\n")
cat("============================================================\n\n")

cat(
    "GPL571:",
    gpl571_nonfinite,
    "\n"
)

cat(
    "GPL3921:",
    gpl3921_nonfinite,
    "\n\n"
)


# ============================================================
# 14C.14 — CHECK SELECTED PROBE UNIQUENESS
# ============================================================

gpl571_duplicate_probes <- sum(
    duplicated(
        gpl571_selected_probes$probe_id
    )
)

gpl3921_duplicate_probes <- sum(
    duplicated(
        gpl3921_selected_probes$probe_id
    )
)


cat("\n============================================================\n")
cat("SELECTED PROBE UNIQUENESS\n")
cat("============================================================\n\n")

cat(
    "GPL571 duplicate selected probes:",
    gpl571_duplicate_probes,
    "\n"
)

cat(
    "GPL3921 duplicate selected probes:",
    gpl3921_duplicate_probes,
    "\n\n"
)


# ============================================================
# 14C.15 — FINAL SELECTION SUMMARY
# ============================================================

gpl571_selection_summary <- gpl571_selected_probes[
    ,
    c(
        "gene",
        "probe_id",
        "median_expression",
        "MAD",
        "SD"
    ),
    drop = FALSE
]

gpl571_selection_summary$platform <- "GPL571"


gpl3921_selection_summary <- gpl3921_selected_probes[
    ,
    c(
        "gene",
        "probe_id",
        "median_expression",
        "MAD",
        "SD"
    ),
    drop = FALSE
]

gpl3921_selection_summary$platform <- "GPL3921"


gse14520_probe_selection_summary <- rbind(
    gpl571_selection_summary[
        ,
        c(
            "platform",
            "gene",
            "probe_id",
            "median_expression",
            "MAD",
            "SD"
        ),
        drop = FALSE
    ],
    gpl3921_selection_summary[
        ,
        c(
            "platform",
            "gene",
            "probe_id",
            "median_expression",
            "MAD",
            "SD"
        ),
        drop = FALSE
    ]
)


cat("\n============================================================\n")
cat("FINAL PROBE SELECTION SUMMARY\n")
cat("============================================================\n\n")

print(
    gse14520_probe_selection_summary
)


# ============================================================
# 14C.16 — SAVE OUTPUTS
# ============================================================

write.csv(
    gse14520_probe_selection_summary,
    file.path(
        results_dir,
        "Step14C_GSE14520_probe_selection_summary.csv"
    ),
    row.names = FALSE
)


write.csv(
    as.data.frame(
        gpl571_gene_expression
    ),
    file.path(
        results_dir,
        "Step14C_GPL571_gene_expression.csv"
    ),
    row.names = TRUE
)


write.csv(
    as.data.frame(
        gpl3921_gene_expression
    ),
    file.path(
        results_dir,
        "Step14C_GPL3921_gene_expression.csv"
    ),
    row.names = TRUE
)


save(
    gpl571_selected_probes,
    gpl3921_selected_probes,
    gpl571_gene_expression,
    gpl3921_gene_expression,
    gse14520_common_genes,
    gse14520_probe_selection_summary,
    file = file.path(
        results_dir,
        "Step14C_GSE14520_gene_level_expression.RData"
    )
)


# ============================================================
# 14C.17 — FINAL STOP
# ============================================================

cat("\n\n============================================================\n")
cat("STEP 14C COMPLETE\n")
cat("============================================================\n\n")

cat(
    "Representative probes selected successfully.\n\n"
)

cat(
    "GPL571:",
    nrow(gpl571_gene_expression),
    "genes x",
    ncol(gpl571_gene_expression),
    "tumors\n"
)

cat(
    "GPL3921:",
    nrow(gpl3921_gene_expression),
    "genes x",
    ncol(gpl3921_gene_expression),
    "tumors\n\n"
)

cat(
    "No biological outcome was used for probe selection.\n"
)

cat(
    "No platform merging was performed.\n"
)

cat(
    "No new genes were introduced.\n"
)

cat(
    "No correlation analysis was performed.\n\n"
)

cat(
    "NEXT STEP:\n",
    "Cross-platform validation of the common gene panel.\n\n"
)

cat(
    "============================================================\n"
)
