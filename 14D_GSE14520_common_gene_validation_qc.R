# ============================================================
# STEP 14D — GSE14520 COMMON 14-GENE VALIDATION QC
# ============================================================
#
# Purpose:
#   1. Define the definitive 14-gene common validation panel
#   2. Create platform-specific 14-gene matrices
#   3. Check missing/non-finite values
#   4. Check zero-variance genes
#   5. Check sample-level expression distributions
#   6. Check preservation of sample IDs
#   7. Save the definitive GSE14520 validation object
#
# IMPORTANT:
#   - No lipid score
#   - No HR score
#   - No DDR score
#   - No correlations
#   - No survival analysis
#   - No new genes
#   - No platform merging
#   - No re-normalization
#   - No [[ ]] indexing
#
# ============================================================


# ------------------------------------------------------------
# 14D.0 — START
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP 14D — GSE14520 COMMON 14-GENE VALIDATION QC\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 14D.1 — CHECK REQUIRED OBJECTS
# ------------------------------------------------------------

required_objects_14d <- c(
    "gpl571_gene_expression",
    "gpl3921_gene_expression",
    "gse14520_common_genes",
    "gpl571_selected_probes",
    "gpl3921_selected_probes",
    "gse14520_gpl571_tumor_ids",
    "gse14520_gpl3921_tumor_ids"
)


missing_objects_14d <- required_objects_14d[
    !sapply(
        required_objects_14d,
        exists
    )
]


if (
    length(missing_objects_14d) > 0
) {

    stop(
        paste0(
            "\nERROR: The following Step 14C objects are missing:\n\n",
            paste(
                missing_objects_14d,
                collapse = "\n"
            ),
            "\n\nPlease run Step 14C first.\n"
        )
    )
}


cat(
    "All required Step 14C objects are present.\n\n"
)


# ------------------------------------------------------------
# 14D.2 — DEFINE THE COMMON 14-GENE PANEL
# ------------------------------------------------------------

gse14520_common_genes <- intersect(
    rownames(
        gpl571_gene_expression
    ),
    rownames(
        gpl3921_gene_expression
    )
)


# Preserve the original frozen-gene ordering
gse14520_common_genes <- gse14520_candidate_genes[
    gse14520_candidate_genes %in%
        gse14520_common_genes
]


cat("\n============================================================\n")
cat("COMMON VALIDATION GENE PANEL\n")
cat("============================================================\n\n")

cat(
    "Number of common genes:",
    length(gse14520_common_genes),
    "\n\n"
)

print(
    gse14520_common_genes
)


if (
    length(gse14520_common_genes) != 14
) {

    stop(
        paste0(
            "\nERROR: Expected 14 common genes but found ",
            length(gse14520_common_genes),
            ".\n"
        )
    )
}


# ------------------------------------------------------------
# 14D.3 — CREATE DEFINITIVE 14-GENE MATRICES
# ------------------------------------------------------------

gse14520_gpl571_14gene <- gpl571_gene_expression[
    gse14520_common_genes,
    ,
    drop = FALSE
]


gse14520_gpl3921_14gene <- gpl3921_gene_expression[
    gse14520_common_genes,
    ,
    drop = FALSE
]


# ------------------------------------------------------------
# 14D.4 — VERIFY MATRIX DIMENSIONS
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("14-GENE MATRIX DIMENSIONS\n")
cat("============================================================\n\n")

cat(
    "GPL571:",
    nrow(gse14520_gpl571_14gene),
    "genes x",
    ncol(gse14520_gpl571_14gene),
    "tumors\n"
)

cat(
    "GPL3921:",
    nrow(gse14520_gpl3921_14gene),
    "genes x",
    ncol(gse14520_gpl3921_14gene),
    "tumors\n\n"
)


if (
    nrow(gse14520_gpl571_14gene) != 14
) {

    stop(
        "\nERROR: GPL571 does not contain exactly 14 common genes.\n"
    )
}


if (
    nrow(gse14520_gpl3921_14gene) != 14
) {

    stop(
        "\nERROR: GPL3921 does not contain exactly 14 common genes.\n"
    )
}


# ------------------------------------------------------------
# 14D.5 — VERIFY GENE ORDER
# ------------------------------------------------------------

if (
    !identical(
        rownames(gse14520_gpl571_14gene),
        gse14520_common_genes
    )
) {

    stop(
        "\nERROR: GPL571 gene order is incorrect.\n"
    )
}


if (
    !identical(
        rownames(gse14520_gpl3921_14gene),
        gse14520_common_genes
    )
) {

    stop(
        "\nERROR: GPL3921 gene order is incorrect.\n"
    )
}


cat(
    "Gene order verified for both platforms.\n\n"
)


# ------------------------------------------------------------
# 14D.6 — VERIFY SAMPLE COUNTS
# ------------------------------------------------------------

if (
    ncol(gse14520_gpl571_14gene) !=
    length(gse14520_gpl571_tumor_ids)
) {

    stop(
        "\nERROR: GPL571 sample count does not match tumor IDs.\n"
    )
}


if (
    ncol(gse14520_gpl3921_14gene) !=
    length(gse14520_gpl3921_tumor_ids)
) {

    stop(
        "\nERROR: GPL3921 sample count does not match tumor IDs.\n"
    )
}


cat(
    "GPL571 tumor samples:",
    ncol(gse14520_gpl571_14gene),
    "\n"
)

cat(
    "GPL3921 tumor samples:",
    ncol(gse14520_gpl3921_14gene),
    "\n\n"
)


# ------------------------------------------------------------
# 14D.7 — VERIFY SAMPLE IDENTITIES
# ------------------------------------------------------------

if (
    !identical(
        colnames(gse14520_gpl571_14gene),
        gse14520_gpl571_tumor_ids
    )
) {

    stop(
        "\nERROR: GPL571 sample identities changed.\n"
    )
}


if (
    !identical(
        colnames(gse14520_gpl3921_14gene),
        gse14520_gpl3921_tumor_ids
    )
) {

    stop(
        "\nERROR: GPL3921 sample identities changed.\n"
    )
}


cat(
    "Sample identities verified for both platforms.\n\n"
)


# ------------------------------------------------------------
# 14D.8 — CHECK DUPLICATE GENE NAMES
# ------------------------------------------------------------

gpl571_duplicate_genes <- sum(
    duplicated(
        rownames(
            gse14520_gpl571_14gene
        )
    )
)


gpl3921_duplicate_genes <- sum(
    duplicated(
        rownames(
            gse14520_gpl3921_14gene
        )
    )
)


cat("\n============================================================\n")
cat("DUPLICATE GENE CHECK\n")
cat("============================================================\n\n")

cat(
    "GPL571 duplicate genes:",
    gpl571_duplicate_genes,
    "\n"
)

cat(
    "GPL3921 duplicate genes:",
    gpl3921_duplicate_genes,
    "\n\n"
)


if (
    gpl571_duplicate_genes > 0
) {

    stop(
        "\nERROR: Duplicate genes detected in GPL571.\n"
    )
}


if (
    gpl3921_duplicate_genes > 0
) {

    stop(
        "\nERROR: Duplicate genes detected in GPL3921.\n"
    )
}


# ------------------------------------------------------------
# 14D.9 — CHECK NON-FINITE VALUES
# ------------------------------------------------------------

gpl571_nonfinite <- sum(
    !is.finite(
        gse14520_gpl571_14gene
    )
)


gpl3921_nonfinite <- sum(
    !is.finite(
        gse14520_gpl3921_14gene
    )
)


cat("\n============================================================\n")
cat("NON-FINITE VALUE CHECK\n")
cat("============================================================\n\n")

cat(
    "GPL571 non-finite values:",
    gpl571_nonfinite,
    "\n"
)

cat(
    "GPL3921 non-finite values:",
    gpl3921_nonfinite,
    "\n\n"
)


if (
    gpl571_nonfinite > 0
) {

    stop(
        "\nERROR: Non-finite values detected in GPL571.\n"
    )
}


if (
    gpl3921_nonfinite > 0
) {

    stop(
        "\nERROR: Non-finite values detected in GPL3921.\n"
    )
}


# ------------------------------------------------------------
# 14D.10 — CHECK ZERO-VARIANCE GENES
# ------------------------------------------------------------

gpl571_gene_sd <- apply(
    gse14520_gpl571_14gene,
    1,
    sd,
    na.rm = TRUE
)


gpl3921_gene_sd <- apply(
    gse14520_gpl3921_14gene,
    1,
    sd,
    na.rm = TRUE
)


gpl571_zero_variance <- names(
    gpl571_gene_sd[
        !is.finite(gpl571_gene_sd) |
            gpl571_gene_sd == 0
    ]
)


gpl3921_zero_variance <- names(
    gpl3921_gene_sd[
        !is.finite(gpl3921_gene_sd) |
            gpl3921_gene_sd == 0
    ]
)


cat("\n============================================================\n")
cat("ZERO-VARIANCE GENE CHECK\n")
cat("============================================================\n\n")

cat(
    "GPL571 zero-variance genes:",
    length(gpl571_zero_variance),
    "\n"
)

if (
    length(gpl571_zero_variance) > 0
) {
    print(gpl571_zero_variance)
}


cat(
    "\nGPL3921 zero-variance genes:",
    length(gpl3921_zero_variance),
    "\n"
)

if (
    length(gpl3921_zero_variance) > 0
) {
    print(gpl3921_zero_variance)
}


# ------------------------------------------------------------
# 14D.11 — CREATE GENE-LEVEL QC TABLE
# ------------------------------------------------------------

gse14520_gene_qc <- data.frame(
    gene = gse14520_common_genes,

    GPL571_mean = sapply(
        gse14520_common_genes,
        function(gene) {
            mean(
                gse14520_gpl571_14gene[
                    gene,
                    ,
                    drop = TRUE
                ],
                na.rm = TRUE
            )
        }
    ),

    GPL571_median = sapply(
        gse14520_common_genes,
        function(gene) {
            median(
                gse14520_gpl571_14gene[
                    gene,
                    ,
                    drop = TRUE
                ],
                na.rm = TRUE
            )
        }
    ),

    GPL571_SD = sapply(
        gse14520_common_genes,
        function(gene) {
            sd(
                gse14520_gpl571_14gene[
                    gene,
                    ,
                    drop = TRUE
                ],
                na.rm = TRUE
            )
        }
    ),

    GPL3921_mean = sapply(
        gse14520_common_genes,
        function(gene) {
            mean(
                gse14520_gpl3921_14gene[
                    gene,
                    ,
                    drop = TRUE
                ],
                na.rm = TRUE
            )
        }
    ),

    GPL3921_median = sapply(
        gse14520_common_genes,
        function(gene) {
            median(
                gse14520_gpl3921_14gene[
                    gene,
                    ,
                    drop = TRUE
                ],
                na.rm = TRUE
            )
        }
    ),

    GPL3921_SD = sapply(
        gse14520_common_genes,
        function(gene) {
            sd(
                gse14520_gpl3921_14gene[
                    gene,
                    ,
                    drop = TRUE
                ],
                na.rm = TRUE
            )
        }
    ),

    stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 14D.12 — SAMPLE-LEVEL QC
# ------------------------------------------------------------

gpl571_sample_mean <- colMeans(
    gse14520_gpl571_14gene,
    na.rm = TRUE
)


gpl571_sample_median <- apply(
    gse14520_gpl571_14gene,
    2,
    median,
    na.rm = TRUE
)


gpl3921_sample_mean <- colMeans(
    gse14520_gpl3921_14gene,
    na.rm = TRUE
)


gpl3921_sample_median <- apply(
    gse14520_gpl3921_14gene,
    2,
    median,
    na.rm = TRUE
)


gpl571_sample_qc <- data.frame(
    sample_id = colnames(
        gse14520_gpl571_14gene
    ),
    mean_expression = as.numeric(
        gpl571_sample_mean
    ),
    median_expression = as.numeric(
        gpl571_sample_median
    ),
    stringsAsFactors = FALSE
)


gpl3921_sample_qc <- data.frame(
    sample_id = colnames(
        gse14520_gpl3921_14gene
    ),
    mean_expression = as.numeric(
        gpl3921_sample_mean
    ),
    median_expression = as.numeric(
        gpl3921_sample_median
    ),
    stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 14D.13 — SAMPLE QC RANGE
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("SAMPLE-LEVEL QC\n")
cat("============================================================\n\n")

cat("GPL571 mean expression range:\n")

print(
    range(
        gpl571_sample_qc$mean_expression,
        na.rm = TRUE
    )
)


cat("\nGPL571 median expression range:\n")

print(
    range(
        gpl571_sample_qc$median_expression,
        na.rm = TRUE
    )
)


cat("\nGPL3921 mean expression range:\n")

print(
    range(
        gpl3921_sample_qc$mean_expression,
        na.rm = TRUE
    )
)


cat("\nGPL3921 median expression range:\n")

print(
    range(
        gpl3921_sample_qc$median_expression,
        na.rm = TRUE
    )
)


# ------------------------------------------------------------
# 14D.14 — FINAL QC SUMMARY
# ------------------------------------------------------------

gse14520_14gene_qc_summary <- data.frame(
    platform = c(
        "GPL571",
        "GPL3921"
    ),

    HCC_tumor_samples = c(
        ncol(
            gse14520_gpl571_14gene
        ),
        ncol(
            gse14520_gpl3921_14gene
        )
    ),

    common_genes = c(
        nrow(
            gse14520_gpl571_14gene
        ),
        nrow(
            gse14520_gpl3921_14gene
        )
    ),

    nonfinite_values = c(
        gpl571_nonfinite,
        gpl3921_nonfinite
    ),

    zero_variance_genes = c(
        length(gpl571_zero_variance),
        length(gpl3921_zero_variance)
    ),

    stringsAsFactors = FALSE
)


cat("\n\n============================================================\n")
cat("FINAL STEP 14D QC SUMMARY\n")
cat("============================================================\n\n")

print(
    gse14520_14gene_qc_summary
)


# ------------------------------------------------------------
# 14D.15 — SAVE QC RESULTS
# ------------------------------------------------------------

write.csv(
    gse14520_14gene_qc_summary,
    file.path(
        results_dir,
        "Step14D_GSE14520_14gene_QC_summary.csv"
    ),
    row.names = FALSE
)


write.csv(
    gse14520_gene_qc,
    file.path(
        results_dir,
        "Step14D_GSE14520_gene_level_QC.csv"
    ),
    row.names = FALSE
)


write.csv(
    gpl571_sample_qc,
    file.path(
        results_dir,
        "Step14D_GPL571_sample_QC.csv"
    ),
    row.names = FALSE
)


write.csv(
    gpl3921_sample_qc,
    file.path(
        results_dir,
        "Step14D_GPL3921_sample_QC.csv"
    ),
    row.names = FALSE
)


write.csv(
    as.data.frame(
        gse14520_gpl571_14gene
    ),
    file.path(
        results_dir,
        "Step14D_GPL571_14gene_validation_matrix.csv"
    ),
    row.names = TRUE
)


write.csv(
    as.data.frame(
        gse14520_gpl3921_14gene
    ),
    file.path(
        results_dir,
        "Step14D_GPL3921_14gene_validation_matrix.csv"
    ),
    row.names = TRUE
)


# ------------------------------------------------------------
# 14D.16 — SAVE DEFINITIVE VALIDATION OBJECT
# ------------------------------------------------------------

save(
    gse14520_common_genes,

    gse14520_gpl571_14gene,
    gse14520_gpl3921_14gene,

    gpl571_sample_qc,
    gpl3921_sample_qc,

    gse14520_gene_qc,
    gse14520_14gene_qc_summary,

    gpl571_zero_variance,
    gpl3921_zero_variance,

    file = file.path(
        results_dir,
        "Step14D_GSE14520_14gene_validation_QC.RData"
    )
)


# ------------------------------------------------------------
# 14D.17 — FINAL STATUS
# ------------------------------------------------------------

cat("\n\n============================================================\n")
cat("STEP 14D COMPLETE\n")
cat("============================================================\n\n")

cat(
    "Definitive common validation panel:",
    length(gse14520_common_genes),
    "genes\n"
)

cat(
    "GPL571:",
    ncol(gse14520_gpl571_14gene),
    "HCC tumors\n"
)

cat(
    "GPL3921:",
    ncol(gse14520_gpl3921_14gene),
    "HCC tumors\n\n"
)

cat(
    "Non-finite values — GPL571:",
    gpl571_nonfinite,
    "\n"
)

cat(
    "Non-finite values — GPL3921:",
    gpl3921_nonfinite,
    "\n\n"
)

cat(
    "Zero-variance genes — GPL571:",
    length(gpl571_zero_variance),
    "\n"
)

cat(
    "Zero-variance genes — GPL3921:",
    length(gpl3921_zero_variance),
    "\n\n"
)

cat(
    "No biological analysis has been performed.\n"
)

cat(
    "No correlations have been calculated.\n"
)

cat(
    "No scores have been calculated.\n\n"
)

cat(
    "NEXT STEP:\n",
    "Proceed to the predefined GSE14520 gene-level validation analysis.\n\n"
)

cat(
    "============================================================\n"
)
