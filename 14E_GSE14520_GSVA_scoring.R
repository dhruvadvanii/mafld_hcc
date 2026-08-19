# ============================================================
# STEP 14E — GSE14520 COMPLETE 17-PATHWAY GSVA
# ============================================================
#
# PURPOSE
# ------------------------------------------------------------
# Reproduce the complete 17-pathway scoring framework used
# for GSE76427 in GSE14520.
#
# PLATFORMS:
#   GPL571  = 22 HCC tumors
#   GPL3921 = 225 HCC tumors
#
# 17 PATHWAYS:
#
# LIPID PATHWAYS
#   1. HALLMARK_FATTY_ACID_METABOLISM
#   2. CHOLESTEROL_BIOSYNTHESIS
#   3. FATTY_ACYL_COA_BIOSYNTHESIS
#   4. MITOCHONDRIAL_FATTY_ACID_BETA_OXIDATION
#   5. PEROXISOMAL_LIPID_METABOLISM
#   6. SREBP_CHOLESTEROL_REGULATION
#   7. PPAR_ALPHA_LIPID_REGULATION
#   8. TRIGLYCERIDE_METABOLISM
#
# CORE DDR / REPAIR
#   9.  HOMOLOGOUS_RECOMBINATION
#   10. MISMATCH_REPAIR
#   11. BASE_EXCISION_REPAIR
#   12. DOUBLE_STRAND_BREAK_REPAIR
#   13. NON_HOMOLOGOUS_END_JOINING
#   14. HALLMARK_DNA_REPAIR
#
# ADDITIONAL PATHWAYS
#   15. P53_PATHWAY
#   16. P53_DEPENDENT_G1_DNA_DAMAGE_RESPONSE
#   17. G2M_CHECKPOINT
#
#
# PRIMARY SCORES
# ------------------------------------------------------------
# Lipid:
#   HALLMARK_FATTY_ACID_METABOLISM
#
# HR:
#   HOMOLOGOUS_RECOMBINATION ONLY
#
# DDR:
#   mean of the six core DDR/repair pathways
#
# NHEJ:
#   NON_HOMOLOGOUS_END_JOINING
#
# P53:
#   retained separately
#   NOT included in DDR composite
#
# EXPANDED HR GENE LIST:
#   NOT used for HR score
#   reserved for gene-level analysis
#
# IMPORTANT:
# ------------------------------------------------------------
# This code deliberately contains NO [[ indexing.
# All list/matrix extraction uses [ ] only.
#
# ============================================================


cat("\n")
cat("============================================================\n")
cat("STEP 14E — GSE14520 COMPLETE 17-PATHWAY GSVA\n")
cat("============================================================\n\n")


# ============================================================
# 14E.1 — CHECK REQUIRED OBJECTS
# ============================================================

cat("14E.1 — Checking required objects...\n\n")


required_objects_14e <- c(
    "gse14520_gpl571_tumor_expression",
    "gse14520_gpl3921_tumor_expression",
    "final_gene_sets"
)


missing_objects_14e <- required_objects_14e[
    !sapply(
        required_objects_14e,
        exists,
        envir = .GlobalEnv
    )
]


if (
    length(missing_objects_14e) > 0
) {

    stop(
        paste0(
            "ERROR: Required object(s) missing:\n\n",
            paste(
                missing_objects_14e,
                collapse = "\n"
            ),
            "\n\nRun Steps 14A–14D before Step 14E."
        )
    )
}


cat(
    "All required objects are present.\n\n"
)


# ============================================================
# 14E.2 — DEFINE THE FROZEN 17-PATHWAY INVENTORY
# ============================================================

cat("14E.2 — Defining the 17-pathway inventory...\n\n")


gse14520_17_pathways <- c(

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


if (
    length(gse14520_17_pathways) != 17
) {

    stop(
        "ERROR: 17-pathway inventory does not contain exactly 17 pathways."
    )
}


cat(
    "Number of pathways defined:",
    length(gse14520_17_pathways),
    "\n\n"
)


print(
    gse14520_17_pathways
)


# ============================================================
# 14E.3 — VERIFY PATHWAYS EXIST IN final_gene_sets
# ============================================================

cat("\n")
cat("14E.3 — Checking pathway availability in final_gene_sets...\n\n")


final_gene_set_names <- names(
    final_gene_sets
)


missing_final_gene_sets <- setdiff(
    gse14520_17_pathways,
    final_gene_set_names
)


if (
    length(missing_final_gene_sets) > 0
) {

    stop(
        paste0(
            "ERROR: The following required gene sets are missing ",
            "from final_gene_sets:\n\n",
            paste(
                missing_final_gene_sets,
                collapse = "\n"
            )
        )
    )
}


cat(
    "All 17 pathways are present in final_gene_sets.\n\n"
)


# ============================================================
# 14E.4 — PREPARE EXPRESSION MATRICES
# ============================================================

cat("14E.4 — Preparing full GSE14520 expression matrices...\n\n")


gpl571_expression_14e <- as.matrix(
    gse14520_gpl571_tumor_expression
)


gpl3921_expression_14e <- as.matrix(
    gse14520_gpl3921_tumor_expression
)


storage.mode(
    gpl571_expression_14e
) <- "numeric"


storage.mode(
    gpl3921_expression_14e
) <- "numeric"


if (
    is.null(
        rownames(gpl571_expression_14e)
    )
) {

    stop(
        "ERROR: GPL571 expression matrix has no probe IDs."
    )
}


if (
    is.null(
        colnames(gpl571_expression_14e)
    )
) {

    stop(
        "ERROR: GPL571 expression matrix has no sample IDs."
    )
}


if (
    is.null(
        rownames(gpl3921_expression_14e)
    )
) {

    stop(
        "ERROR: GPL3921 expression matrix has no probe IDs."
    )
}


if (
    is.null(
        colnames(gpl3921_expression_14e)
    )
) {

    stop(
        "ERROR: GPL3921 expression matrix has no sample IDs."
    )
}


cat(
    "GPL571 expression matrix: ",
    nrow(gpl571_expression_14e),
    " probes x ",
    ncol(gpl571_expression_14e),
    " tumors\n",
    sep = ""
)


cat(
    "GPL3921 expression matrix: ",
    nrow(gpl3921_expression_14e),
    " probes x ",
    ncol(gpl3921_expression_14e),
    " tumors\n\n",
    sep = ""
)


# ============================================================
# 14E.5 — LOCATE EXISTING PROBE-TO-GENE MAPS
# ============================================================

cat("14E.5 — Locating probe-to-gene annotation maps...\n\n")


mapping_candidates_571 <- c(
    "gpl571_candidate_probe_map",
    "gpl571_probe_map",
    "gpl571_annotation",
    "gpl571_gene_annotation",
    "gpl571_probe_annotation"
)


mapping_candidates_3921 <- c(
    "gpl3921_candidate_probe_map",
    "gpl3921_probe_map",
    "gpl3921_annotation",
    "gpl3921_gene_annotation",
    "gpl3921_probe_annotation"
)


find_mapping_14e <- function(
    candidate_names,
    expression_matrix,
    platform_name
) {

    existing_names <- candidate_names[
        sapply(
            candidate_names,
            exists,
            envir = .GlobalEnv
        )
    ]


    if (
        length(existing_names) == 0
    ) {

        return(NULL)
    }


    for (
        object_name in existing_names
    ) {

        candidate <- get(
            object_name,
            envir = .GlobalEnv
        )


        if (
            !is.data.frame(candidate)
        ) {

            next
        }


        column_names <- names(
            candidate
        )


        probe_candidates <- column_names[
            grepl(
                "probe|id_ref|probe_id",
                column_names,
                ignore.case = TRUE
            )
        ]


        gene_candidates <- column_names[
            grepl(
                "gene.*symbol|gene_symbol|symbol",
                column_names,
                ignore.case = TRUE
            )
        ]


        if (
            length(probe_candidates) == 0 ||
            length(gene_candidates) == 0
        ) {

            next
        }


        probe_column <- probe_candidates[1]
        gene_column <- gene_candidates[1]


        candidate_probes <- trimws(
            as.character(
                candidate[
                    ,
                    probe_column
                ]
            )
        )


        overlap_count <- sum(
            rownames(expression_matrix) %in%
                candidate_probes
        )


        if (
            overlap_count > 100
        ) {

            cat(
                platform_name,
                " mapping selected: ",
                object_name,
                "\n",
                sep = ""
            )

            cat(
                "Probe column: ",
                probe_column,
                "\n",
                sep = ""
            )

            cat(
                "Gene column: ",
                gene_column,
                "\n",
                sep = ""
            )

            cat(
                "Probe overlap: ",
                overlap_count,
                "\n\n",
                sep = ""
            )


            return(
                list(
                    object_name = object_name,
                    map = candidate,
                    probe_column = probe_column,
                    gene_column = gene_column
                )
            )
        }
    }


    NULL
}


gpl571_map_info_14e <- find_mapping_14e(
    mapping_candidates_571,
    gpl571_expression_14e,
    "GPL571"
)


gpl3921_map_info_14e <- find_mapping_14e(
    mapping_candidates_3921,
    gpl3921_expression_14e,
    "GPL3921"
)


if (
    is.null(gpl571_map_info_14e)
) {

    stop(
        "ERROR: Could not identify a valid GPL571 probe-to-gene map."
    )
}


if (
    is.null(gpl3921_map_info_14e)
) {

    stop(
        "ERROR: Could not identify a valid GPL3921 probe-to-gene map."
    )
}


# ============================================================
# 14E.6 — STANDARDIZE ANNOTATION MAPS
# ============================================================

cat("14E.6 — Standardizing annotation maps...\n\n")


standardize_mapping_14e <- function(
    map_info,
    expression_matrix
) {

    annotation <- map_info$map


    result <- data.frame(

        probe_id =
            trimws(
                as.character(
                    annotation[
                        ,
                        map_info$probe_column
                    ]
                )
            ),

        gene_symbol =
            toupper(
                trimws(
                    as.character(
                        annotation[
                            ,
                            map_info$gene_column
                        ]
                    )
                )
            ),

        stringsAsFactors = FALSE
    )


    result <- result[
        !is.na(result$probe_id) &
        result$probe_id != "" &
        !is.na(result$gene_symbol) &
        result$gene_symbol != "",
        ,
        drop = FALSE
    ]


    result <- result[
        result$probe_id %in%
            rownames(expression_matrix),
        ,
        drop = FALSE
    ]


    result <- result[
        !duplicated(result$probe_id),
        ,
        drop = FALSE
    ]


    result
}


gpl571_mapping_14e <- standardize_mapping_14e(
    gpl571_map_info_14e,
    gpl571_expression_14e
)


gpl3921_mapping_14e <- standardize_mapping_14e(
    gpl3921_map_info_14e,
    gpl3921_expression_14e
)


cat(
    "GPL571 usable probe annotations: ",
    nrow(gpl571_mapping_14e),
    "\n",
    sep = ""
)


cat(
    "GPL3921 usable probe annotations: ",
    nrow(gpl3921_mapping_14e),
    "\n\n",
    sep = ""
)


# ============================================================
# 14E.7 — AGGREGATE PROBES TO GENES
# ============================================================

cat("14E.7 — Creating full gene-level expression matrices...\n\n")


aggregate_probes_to_genes_14e <- function(
    expression_matrix,
    mapping
) {

    common_probes <- intersect(
        rownames(expression_matrix),
        mapping$probe_id
    )


    if (
        length(common_probes) == 0
    ) {

        stop(
            "ERROR: No expression probes overlap the annotation map."
        )
    }


    expression_subset <- expression_matrix[
        common_probes,
        ,
        drop = FALSE
    ]


    mapping_subset <- mapping[
        match(
            common_probes,
            mapping$probe_id
        ),
        ,
        drop = FALSE
    ]


    gene_names <- unique(
        mapping_subset$gene_symbol
    )


    gene_matrix <- matrix(
        NA_real_,
        nrow = length(gene_names),
        ncol = ncol(expression_subset)
    )


    rownames(gene_matrix) <- gene_names

    colnames(gene_matrix) <-
        colnames(expression_subset)


    for (
        i in seq_along(gene_names)
    ) {

        current_gene <- gene_names[i]


        probe_indices <- which(
            mapping_subset$gene_symbol ==
                current_gene
        )


        current_values <- expression_subset[
            probe_indices,
            ,
            drop = FALSE
        ]


        if (
            nrow(current_values) == 1
        ) {

            gene_matrix[
                i,
                ] <-
                as.numeric(
                    current_values[
                        1,
                        ]
                )

        } else {

            gene_matrix[
                i,
                ] <-
                apply(
                    current_values,
                    2,
                    median,
                    na.rm = TRUE
                )
        }
    }


    gene_matrix
}


gpl571_gene_expression_14e <-
    aggregate_probes_to_genes_14e(
        gpl571_expression_14e,
        gpl571_mapping_14e
    )


gpl3921_gene_expression_14e <-
    aggregate_probes_to_genes_14e(
        gpl3921_expression_14e,
        gpl3921_mapping_14e
    )


cat(
    "GPL571 gene-level matrix: ",
    nrow(gpl571_gene_expression_14e),
    " genes x ",
    ncol(gpl571_gene_expression_14e),
    " tumors\n",
    sep = ""
)


cat(
    "GPL3921 gene-level matrix: ",
    nrow(gpl3921_gene_expression_14e),
    " genes x ",
    ncol(gpl3921_gene_expression_14e),
    " tumors\n\n",
    sep = ""
)


# ============================================================
# 14E.8 — PREPARE 17 GENE SETS
# ============================================================

cat("14E.8 — Preparing platform-specific 17 gene sets...\n\n")


prepare_gene_sets_14e <- function(
    gene_sets,
    expression_matrix,
    pathway_names
) {

    output <- vector(
        mode = "list",
        length = length(pathway_names)
    )


    names(output) <- pathway_names


    for (
        i in seq_along(pathway_names)
    ) {

        pathway_name <- pathway_names[i]


        pathway_index <- match(
            pathway_name,
            names(gene_sets)
        )


        if (
            is.na(pathway_index)
        ) {

            stop(
                paste0(
                    "ERROR: Missing gene set: ",
                    pathway_name
                )
            )
        }


        pathway_genes <- gene_sets[
            pathway_index
        ]


        pathway_genes <- unlist(
            pathway_genes,
            use.names = FALSE
        )


        pathway_genes <- unique(
            toupper(
                trimws(
                    as.character(
                        pathway_genes
                    )
                )
            )
        )


        pathway_genes <- pathway_genes[
            !is.na(pathway_genes) &
            pathway_genes != ""
        ]


        pathway_genes <- intersect(
            pathway_genes,
            rownames(expression_matrix)
        )


        if (
            length(pathway_genes) == 0
        ) {

            stop(
                paste0(
                    "ERROR: No genes from pathway ",
                    pathway_name,
                    " are present on this platform."
                )
            )
        }


        output[
            i
        ] <- list(
            pathway_genes
        )
    }


    output
}


gpl571_gene_sets_14e <-
    prepare_gene_sets_14e(
        final_gene_sets,
        gpl571_gene_expression_14e,
        gse14520_17_pathways
    )


gpl3921_gene_sets_14e <-
    prepare_gene_sets_14e(
        final_gene_sets,
        gpl3921_gene_expression_14e,
        gse14520_17_pathways
    )


# ============================================================
# 14E.9 — CALCULATE GENE-SET COVERAGE
# ============================================================

cat("14E.9 — Calculating gene-set coverage...\n\n")


calculate_coverage_14e <- function(
    gene_sets,
    final_sets,
    pathway_names,
    platform_name
) {

    coverage_rows <- vector(
        mode = "list",
        length = length(pathway_names)
    )


    for (
        i in seq_along(pathway_names)
    ) {

        pathway_name <- pathway_names[i]


        final_index <- match(
            pathway_name,
            names(final_sets)
        )


        original_genes <- final_sets[
            final_index
        ]


        original_genes <- unlist(
            original_genes,
            use.names = FALSE
        )


        original_genes <- unique(
            toupper(
                trimws(
                    as.character(
                        original_genes
                    )
                )
            )
        )


        original_genes <- original_genes[
            !is.na(original_genes) &
            original_genes != ""
        ]


        current_index <- match(
            pathway_name,
            names(gene_sets)
        )


        present_genes <- gene_sets[
            current_index
        ]


        present_genes <- unlist(
            present_genes,
            use.names = FALSE
        )


        coverage_rows[
            i
        ] <- list(
            data.frame(
                platform =
                    platform_name,

                pathway =
                    pathway_name,

                n_genes_total =
                    length(original_genes),

                n_genes_present =
                    length(present_genes),

                pct_coverage =
                    round(
                        100 *
                        length(present_genes) /
                        length(original_genes),
                        2
                    ),

                stringsAsFactors = FALSE
            )
        )
    }


    do.call(
        rbind,
        coverage_rows
    )
}


gse14520_gpl571_17pathway_coverage <-
    calculate_coverage_14e(
        gpl571_gene_sets_14e,
        final_gene_sets,
        gse14520_17_pathways,
        "GPL571"
    )


gse14520_gpl3921_17pathway_coverage <-
    calculate_coverage_14e(
        gpl3921_gene_sets_14e,
        final_gene_sets,
        gse14520_17_pathways,
        "GPL3921"
    )


cat(
    "GPL571 coverage:\n\n"
)

print(
    gse14520_gpl571_17pathway_coverage,
    row.names = FALSE
)


cat(
    "\nGPL3921 coverage:\n\n"
)

print(
    gse14520_gpl3921_17pathway_coverage,
    row.names = FALSE
)


# ============================================================
# 14E.10 — CHECK REQUIRED SCORING PATHWAYS
# ============================================================

cat("\n")
cat("14E.10 — Checking required scoring pathways...\n\n")


required_scoring_pathways_14e <- c(

    "HALLMARK_FATTY_ACID_METABOLISM",

    "HOMOLOGOUS_RECOMBINATION",

    "MISMATCH_REPAIR",

    "BASE_EXCISION_REPAIR",

    "DOUBLE_STRAND_BREAK_REPAIR",

    "NON_HOMOLOGOUS_END_JOINING",

    "HALLMARK_DNA_REPAIR"

)


for (
    pathway_name in required_scoring_pathways_14e
) {

    index_571 <- match(
        pathway_name,
        gse14520_gpl571_17pathway_coverage$pathway
    )


    index_3921 <- match(
        pathway_name,
        gse14520_gpl3921_17pathway_coverage$pathway
    )


    coverage_571 <-
        gse14520_gpl571_17pathway_coverage[
            index_571,
            "pct_coverage"
        ]


    coverage_3921 <-
        gse14520_gpl3921_17pathway_coverage[
            index_3921,
            "pct_coverage"
        ]


    if (
        coverage_571 <= 0
    ) {

        stop(
            paste0(
                "ERROR: ",
                pathway_name,
                " has zero GPL571 coverage."
            )
        )
    }


    if (
        coverage_3921 <= 0
    ) {

        stop(
            paste0(
                "ERROR: ",
                pathway_name,
                " has zero GPL3921 coverage."
            )
        )
    }
}


cat(
    "All required scoring pathways have non-zero coverage.\n\n"
)


# ============================================================
# 14E.11 — RUN GSVA
# ============================================================
#
# IMPORTANT:
# The installed GSVA version uses the current parameter-object
# API. We therefore construct the parameter object using
# positional arguments rather than the incompatible
# expressionData = / geneSets = argument names.
#
# No [[ indexing is used.
# ============================================================

cat("14E.11 — Running GSVA...\n\n")


if (
    !requireNamespace(
        "GSVA",
        quietly = TRUE
    )
) {

    stop(
        "ERROR: GSVA package is not installed."
    )
}


# ------------------------------------------------------------
# Report installed GSVA version
# ------------------------------------------------------------

cat(
    "Installed GSVA version: ",
    as.character(
        packageVersion("GSVA")
    ),
    "\n\n",
    sep = ""
)


# ------------------------------------------------------------
# GSVA runner
# ------------------------------------------------------------

run_gsva_14e <- function(
    expression_matrix,
    gene_sets,
    platform_name
) {

    cat(
        "Running GSVA for ",
        platform_name,
        "...\n",
        sep = ""
    )


    # --------------------------------------------------------
    # Verify input dimensions
    # --------------------------------------------------------

    if (
        !is.matrix(expression_matrix)
    ) {

        expression_matrix <- as.matrix(
            expression_matrix
        )
    }


    storage.mode(
        expression_matrix
    ) <- "numeric"


    if (
        nrow(expression_matrix) == 0 ||
        ncol(expression_matrix) == 0
    ) {

        stop(
            paste0(
                "ERROR: ",
                platform_name,
                " expression matrix is empty."
            )
        )
    }


    if (
        length(gene_sets) != 17
    ) {

        stop(
            paste0(
                "ERROR: ",
                platform_name,
                " gene-set collection contains ",
                length(gene_sets),
                " gene sets instead of 17."
            )
        )
    }


    # --------------------------------------------------------
    # Make sure every gene set is a character vector
    # --------------------------------------------------------

    cleaned_gene_sets <- vector(
        mode = "list",
        length = length(gene_sets)
    )


    names(cleaned_gene_sets) <-
        names(gene_sets)


    for (
        i in seq_along(gene_sets)
    ) {

        current_genes <- gene_sets[
            i
        ]


        current_genes <- unlist(
            current_genes,
            use.names = FALSE
        )


        current_genes <- unique(
            as.character(
                current_genes
            )
        )


        current_genes <- current_genes[
            !is.na(current_genes) &
            current_genes != ""
        ]


        current_genes <- intersect(
            current_genes,
            rownames(expression_matrix)
        )


        if (
            length(current_genes) == 0
        ) {

            stop(
                paste0(
                    "ERROR: Gene set '",
                    names(gene_sets)[i],
                    "' has zero genes represented in ",
                    platform_name,
                    "."
                )
            )
        }


        cleaned_gene_sets[
            i
        ] <- list(
            current_genes
        )
    }


    # --------------------------------------------------------
    # Current GSVA parameter-object API
    #
    # Use positional arguments deliberately.
    # This avoids the incompatible expressionData =
    # argument that produced your error.
    # --------------------------------------------------------

    parameter_object <- GSVA::gsvaParam(
        expression_matrix,
        cleaned_gene_sets,
        kcdf = "Gaussian"
    )


    # --------------------------------------------------------
    # Run GSVA
    # --------------------------------------------------------

    result <- GSVA::gsva(
        parameter_object,
        verbose = FALSE
    )


    result <- as.matrix(
        result
    )


    storage.mode(
        result
    ) <- "numeric"


    # --------------------------------------------------------
    # Validate output
    # --------------------------------------------------------

    if (
        is.null(
            rownames(result)
        )
    ) {

        stop(
            paste0(
                "ERROR: ",
                platform_name,
                " GSVA result has no pathway names."
            )
        )
    }


    if (
        is.null(
            colnames(result)
        )
    ) {

        stop(
            paste0(
                "ERROR: ",
                platform_name,
                " GSVA result has no sample IDs."
            )
        )
    }


    cat(
        platform_name,
        " GSVA completed: ",
        nrow(result),
        " pathways x ",
        ncol(result),
        " samples\n\n",
        sep = ""
    )


    result
}


# ============================================================
# RUN GPL571
# ============================================================

gse14520_gpl571_pathway_scores <-
    run_gsva_14e(
        gpl571_gene_expression_14e,
        gpl571_gene_sets_14e,
        "GPL571"
    )


# ============================================================
# RUN GPL3921
# ============================================================

gse14520_gpl3921_pathway_scores <-
    run_gsva_14e(
        gpl3921_gene_expression_14e,
        gpl3921_gene_sets_14e,
        "GPL3921"
    )


# ============================================================
# IMMEDIATE OUTPUT CHECK
# ============================================================

cat(
    "============================================================\n"
)

cat(
    "GSVA OUTPUT CHECK\n"
)

cat(
    "============================================================\n\n"
)


cat(
    "GPL571 dimensions: ",
    nrow(
        gse14520_gpl571_pathway_scores
    ),
    " pathways x ",
    ncol(
        gse14520_gpl571_pathway_scores
    ),
    " samples\n",
    sep = ""
)


cat(
    "GPL3921 dimensions: ",
    nrow(
        gse14520_gpl3921_pathway_scores
    ),
    " pathways x ",
    ncol(
        gse14520_gpl3921_pathway_scores
    ),
    " samples\n\n",
    sep = ""
)


cat(
    "GPL571 pathways:\n"
)

print(
    rownames(
        gse14520_gpl571_pathway_scores
    )
)


cat(
    "\nGPL3921 pathways:\n"
)

print(
    rownames(
        gse14520_gpl3921_pathway_scores
    )
)


cat(
    "\n============================================================\n"
)

cat(
    "14E.11 GSVA COMPLETE\n"
)

cat(
    "============================================================\n"
)
# ============================================================
# 14E.12 — FORCE EXACT PATHWAY ORDER
# ============================================================

cat("\n14E.12 — Validating the 17-pathway matrices...\n\n")


validate_pathway_matrix_14e <- function(
    score_matrix,
    platform_name
) {

    missing_pathways <- setdiff(
        gse14520_17_pathways,
        rownames(score_matrix)
    )


    if (
        length(missing_pathways) > 0
    ) {

        stop(
            paste0(
                "ERROR: ",
                platform_name,
                " is missing pathway(s):\n",
                paste(
                    missing_pathways,
                    collapse = "\n"
                )
            )
        )
    }


    ordered_matrix <- score_matrix[
        gse14520_17_pathways,
        ,
        drop = FALSE
    ]


    if (
        nrow(ordered_matrix) != 17
    ) {

        stop(
            paste0(
                "ERROR: ",
                platform_name,
                " does not contain exactly 17 pathways."
            )
        )
    }


    ordered_matrix
}


gse14520_gpl571_pathway_scores <-
    validate_pathway_matrix_14e(
        gse14520_gpl571_pathway_scores,
        "GPL571"
    )


gse14520_gpl3921_pathway_scores <-
    validate_pathway_matrix_14e(
        gse14520_gpl3921_pathway_scores,
        "GPL3921"
    )


cat(
    "GPL571: ",
    nrow(gse14520_gpl571_pathway_scores),
    " pathways x ",
    ncol(gse14520_gpl571_pathway_scores),
    " tumors\n",
    sep = ""
)


cat(
    "GPL3921: ",
    nrow(gse14520_gpl3921_pathway_scores),
    " pathways x ",
    ncol(gse14520_gpl3921_pathway_scores),
    " tumors\n\n",
    sep = ""
)


# ============================================================
# 14E.13 — GSVA NON-FINITE QC
# ============================================================

cat("14E.13 — Checking GSVA values...\n\n")


gpl571_nonfinite_count_14e <- sum(
    !is.finite(
        gse14520_gpl571_pathway_scores
    )
)


gpl3921_nonfinite_count_14e <- sum(
    !is.finite(
        gse14520_gpl3921_pathway_scores
    )
)


cat(
    "GPL571 non-finite values: ",
    gpl571_nonfinite_count_14e,
    "\n",
    sep = ""
)


cat(
    "GPL3921 non-finite values: ",
    gpl3921_nonfinite_count_14e,
    "\n\n",
    sep = ""
)


if (
    gpl571_nonfinite_count_14e > 0 ||
    gpl3921_nonfinite_count_14e > 0
) {

    stop(
        "ERROR: Non-finite values detected in GSVA results."
    )
}


# ============================================================
# 14E.14 — EXTRACT PRIMARY LIPID SCORES
# ============================================================

cat("14E.14 — Extracting primary lipid scores...\n\n")


extract_score_14e <- function(
    score_matrix,
    pathway_name
) {

    pathway_index <- match(
        pathway_name,
        rownames(score_matrix)
    )


    if (
        is.na(pathway_index)
    ) {

        stop(
            paste0(
                "ERROR: Pathway not found: ",
                pathway_name
            )
        )
    }


    values <- as.numeric(
        score_matrix[
            pathway_index,
            ,
            drop = TRUE
        ]
    )


    names(values) <-
        colnames(score_matrix)


    values
}


gse14520_gpl571_lipid_score <-
    extract_score_14e(
        gse14520_gpl571_pathway_scores,
        "HALLMARK_FATTY_ACID_METABOLISM"
    )


gse14520_gpl3921_lipid_score <-
    extract_score_14e(
        gse14520_gpl3921_pathway_scores,
        "HALLMARK_FATTY_ACID_METABOLISM"
    )


# ============================================================
# 14E.15 — EXTRACT HR SCORES
# ============================================================

cat("14E.15 — Extracting HR scores...\n\n")


gse14520_gpl571_hr_score <-
    extract_score_14e(
        gse14520_gpl571_pathway_scores,
        "HOMOLOGOUS_RECOMBINATION"
    )


gse14520_gpl3921_hr_score <-
    extract_score_14e(
        gse14520_gpl3921_pathway_scores,
        "HOMOLOGOUS_RECOMBINATION"
    )


cat(
    "HR score = HOMOLOGOUS_RECOMBINATION only.\n"
)

cat(
    "Expanded HR gene list is NOT used for this score.\n\n"
)


# ============================================================
# 14E.16 — EXTRACT NHEJ SCORES
# ============================================================

cat("14E.16 — Extracting NHEJ scores...\n\n")


gse14520_gpl571_nhej_score <-
    extract_score_14e(
        gse14520_gpl571_pathway_scores,
        "NON_HOMOLOGOUS_END_JOINING"
    )


gse14520_gpl3921_nhej_score <-
    extract_score_14e(
        gse14520_gpl3921_pathway_scores,
        "NON_HOMOLOGOUS_END_JOINING"
    )


# ============================================================
# 14E.17 — EXTRACT SIX DDR COMPONENTS
# ============================================================

cat("14E.17 — Extracting six DDR components...\n\n")


ddr_component_pathways_14e <- c(

    "HOMOLOGOUS_RECOMBINATION",

    "MISMATCH_REPAIR",

    "BASE_EXCISION_REPAIR",

    "DOUBLE_STRAND_BREAK_REPAIR",

    "NON_HOMOLOGOUS_END_JOINING",

    "HALLMARK_DNA_REPAIR"

)


ddr_indices_571 <- match(
    ddr_component_pathways_14e,
    rownames(
        gse14520_gpl571_pathway_scores
    )
)


ddr_indices_3921 <- match(
    ddr_component_pathways_14e,
    rownames(
        gse14520_gpl3921_pathway_scores
    )
)


if (
    any(
        is.na(ddr_indices_571)
    )
) {

    stop(
        "ERROR: GPL571 DDR component missing."
    )
}


if (
    any(
        is.na(ddr_indices_3921)
    )
) {

    stop(
        "ERROR: GPL3921 DDR component missing."
    )
}


gse14520_gpl571_ddr_components <-
    gse14520_gpl571_pathway_scores[
        ddr_indices_571,
        ,
        drop = FALSE
    ]


rownames(
    gse14520_gpl571_ddr_components
) <-
    ddr_component_pathways_14e


gse14520_gpl3921_ddr_components <-
    gse14520_gpl3921_pathway_scores[
        ddr_indices_3921,
        ,
        drop = FALSE
    ]


rownames(
    gse14520_gpl3921_ddr_components
) <-
    ddr_component_pathways_14e


# ============================================================
# 14E.18 — DERIVE DDR COMPOSITE
# ============================================================

cat("14E.18 — Calculating six-component DDR scores...\n\n")


gse14520_gpl571_ddr_score <-
    colMeans(
        gse14520_gpl571_ddr_components,
        na.rm = FALSE
    )


names(
    gse14520_gpl571_ddr_score
) <-
    colnames(
        gse14520_gpl571_ddr_components
    )


gse14520_gpl3921_ddr_score <-
    colMeans(
        gse14520_gpl3921_ddr_components,
        na.rm = FALSE
    )


names(
    gse14520_gpl3921_ddr_score
) <-
    colnames(
        gse14520_gpl3921_ddr_components
    )


cat(
    "DDR composite = mean of six core DDR pathways.\n"
)

cat(
    "P53 is NOT included in the DDR composite.\n\n"
)


# ============================================================
# 14E.19 — RETAIN ADDITIONAL PATHWAY SCORES
# ============================================================

cat("14E.19 — Extracting additional pathway scores...\n\n")


gse14520_gpl571_p53_score <-
    extract_score_14e(
        gse14520_gpl571_pathway_scores,
        "P53_PATHWAY"
    )


gse14520_gpl3921_p53_score <-
    extract_score_14e(
        gse14520_gpl3921_pathway_scores,
        "P53_PATHWAY"
    )


gse14520_gpl571_p53_g1_score <-
    extract_score_14e(
        gse14520_gpl571_pathway_scores,
        "P53_DEPENDENT_G1_DNA_DAMAGE_RESPONSE"
    )


gse14520_gpl3921_p53_g1_score <-
    extract_score_14e(
        gse14520_gpl3921_pathway_scores,
        "P53_DEPENDENT_G1_DNA_DAMAGE_RESPONSE"
    )


gse14520_gpl571_g2m_score <-
    extract_score_14e(
        gse14520_gpl571_pathway_scores,
        "G2M_CHECKPOINT"
    )


gse14520_gpl3921_g2m_score <-
    extract_score_14e(
        gse14520_gpl3921_pathway_scores,
        "G2M_CHECKPOINT"
    )


# ============================================================
# 14E.20 — CREATE COMPLETE SAMPLE-LEVEL TABLES
# ============================================================

cat("14E.20 — Creating complete sample-level pathway tables...\n\n")


gse14520_gpl571_complete_scores <-
    data.frame(
        sample_id =
            colnames(
                gse14520_gpl571_pathway_scores
            ),
        t(
            gse14520_gpl571_pathway_scores
        ),
        check.names = FALSE,
        stringsAsFactors = FALSE
    )


gse14520_gpl3921_complete_scores <-
    data.frame(
        sample_id =
            colnames(
                gse14520_gpl3921_pathway_scores
            ),
        t(
            gse14520_gpl3921_pathway_scores
        ),
        check.names = FALSE,
        stringsAsFactors = FALSE
    )


# ============================================================
# 14E.21 — CREATE PRIMARY SCORE TABLES
# ============================================================

cat("14E.21 — Creating primary score tables...\n\n")


gse14520_gpl571_primary_scores <-
    data.frame(

        sample_id =
            names(
                gse14520_gpl571_lipid_score
            ),

        lipid_score =
            as.numeric(
                gse14520_gpl571_lipid_score
            ),

        hr_score =
            as.numeric(
                gse14520_gpl571_hr_score
            ),

        ddr_score =
            as.numeric(
                gse14520_gpl571_ddr_score
            ),

        nhej_score =
            as.numeric(
                gse14520_gpl571_nhej_score
            ),

        p53_score =
            as.numeric(
                gse14520_gpl571_p53_score
            ),

        p53_g1_score =
            as.numeric(
                gse14520_gpl571_p53_g1_score
            ),

        g2m_score =
            as.numeric(
                gse14520_gpl571_g2m_score
            ),

        stringsAsFactors = FALSE
    )


gse14520_gpl3921_primary_scores <-
    data.frame(

        sample_id =
            names(
                gse14520_gpl3921_lipid_score
            ),

        lipid_score =
            as.numeric(
                gse14520_gpl3921_lipid_score
            ),

        hr_score =
            as.numeric(
                gse14520_gpl3921_hr_score
            ),

        ddr_score =
            as.numeric(
                gse14520_gpl3921_ddr_score
            ),

        nhej_score =
            as.numeric(
                gse14520_gpl3921_nhej_score
            ),

        p53_score =
            as.numeric(
                gse14520_gpl3921_p53_score
            ),

        p53_g1_score =
            as.numeric(
                gse14520_gpl3921_p53_g1_score
            ),

        g2m_score =
            as.numeric(
                gse14520_gpl3921_g2m_score
            ),

        stringsAsFactors = FALSE
    )


# ============================================================
# 14E.22 — COMPLETE SCORE QC
# ============================================================

cat("14E.22 — Generating score QC...\n\n")


make_score_qc_14e <- function(
    lipid,
    hr,
    ddr,
    nhej,
    platform
) {

    data.frame(

        platform =
            platform,

        score =
            c(
                "lipid_score",
                "hr_score",
                "ddr_score",
                "nhej_score"
            ),

        n_samples =
            c(
                length(lipid),
                length(hr),
                length(ddr),
                length(nhej)
            ),

        n_missing =
            c(
                sum(is.na(lipid)),
                sum(is.na(hr)),
                sum(is.na(ddr)),
                sum(is.na(nhej))
            ),

        mean =
            c(
                mean(lipid, na.rm = TRUE),
                mean(hr, na.rm = TRUE),
                mean(ddr, na.rm = TRUE),
                mean(nhej, na.rm = TRUE)
            ),

        median =
            c(
                median(lipid, na.rm = TRUE),
                median(hr, na.rm = TRUE),
                median(ddr, na.rm = TRUE),
                median(nhej, na.rm = TRUE)
            ),

        sd =
            c(
                sd(lipid, na.rm = TRUE),
                sd(hr, na.rm = TRUE),
                sd(ddr, na.rm = TRUE),
                sd(nhej, na.rm = TRUE)
            ),

        stringsAsFactors = FALSE
    )
}


gse14520_score_qc <-
    rbind(

        make_score_qc_14e(
            gse14520_gpl571_lipid_score,
            gse14520_gpl571_hr_score,
            gse14520_gpl571_ddr_score,
            gse14520_gpl571_nhej_score,
            "GPL571"
        ),

        make_score_qc_14e(
            gse14520_gpl3921_lipid_score,
            gse14520_gpl3921_hr_score,
            gse14520_gpl3921_ddr_score,
            gse14520_gpl3921_nhej_score,
            "GPL3921"
        )
    )


print(
    gse14520_score_qc,
    row.names = FALSE
)


# ============================================================
# 14E.23 — PATHWAY-LEVEL SUMMARY
# ============================================================

cat("\n")
cat("14E.23 — Creating pathway-level summary...\n\n")


summarize_pathways_14e <- function(
    pathway_matrix,
    platform_name
) {

    summary_rows <- vector(
        mode = "list",
        length = nrow(pathway_matrix)
    )


    for (
        i in seq_len(
            nrow(pathway_matrix)
        )
    ) {

        pathway_name <-
            rownames(pathway_matrix)[i]


        values <-
            as.numeric(
                pathway_matrix[
                    i,
                    ,
                    drop = TRUE
                ]
            )


        summary_rows[
            i
        ] <- list(
            data.frame(

                platform =
                    platform_name,

                pathway =
                    pathway_name,

                n_samples =
                    length(values),

                n_missing =
                    sum(is.na(values)),

                mean =
                    mean(
                        values,
                        na.rm = TRUE
                    ),

                median =
                    median(
                        values,
                        na.rm = TRUE
                    ),

                sd =
                    sd(
                        values,
                        na.rm = TRUE
                    ),

                minimum =
                    min(
                        values,
                        na.rm = TRUE
                    ),

                maximum =
                    max(
                        values,
                        na.rm = TRUE
                    ),

                stringsAsFactors = FALSE
            )
        )
    }


    do.call(
        rbind,
        summary_rows
    )
}


gse14520_gpl571_pathway_summary <-
    summarize_pathways_14e(
        gse14520_gpl571_pathway_scores,
        "GPL571"
    )


gse14520_gpl3921_pathway_summary <-
    summarize_pathways_14e(
        gse14520_gpl3921_pathway_scores,
        "GPL3921"
    )


# ============================================================
# 14E.24 — VERIFY EXPECTED SAMPLE COUNTS
# ============================================================

cat("\n")
cat("14E.24 — Verifying expected tumor counts...\n\n")


if (
    ncol(
        gse14520_gpl571_pathway_scores
    ) != 22
) {

    warning(
        paste0(
            "GPL571 expected 22 tumors but GSVA contains ",
            ncol(
                gse14520_gpl571_pathway_scores
            ),
            "."
        )
    )
}


if (
    ncol(
        gse14520_gpl3921_pathway_scores
    ) != 225
) {

    warning(
        paste0(
            "GPL3921 expected 225 tumors but GSVA contains ",
            ncol(
                gse14520_gpl3921_pathway_scores
            ),
            "."
        )
    )
}


# ============================================================
# 14E.25 — SAVE ALL RESULTS
# ============================================================

cat("\n")
cat("14E.25 — Saving all Step 14E results...\n\n")


if (
    !dir.exists("results")
) {

    dir.create(
        "results",
        recursive = TRUE
    )
}


# ------------------------------------------------------------
# Gene-level matrices
# ------------------------------------------------------------

write.csv(
    gpl571_gene_expression_14e,
    file.path(
        "results",
        "Step14E_GSE14520_GPL571_full_gene_expression.csv"
    ),
    row.names = TRUE
)


write.csv(
    gpl3921_gene_expression_14e,
    file.path(
        "results",
        "Step14E_GSE14520_GPL3921_full_gene_expression.csv"
    ),
    row.names = TRUE
)


# ------------------------------------------------------------
# 17-pathway GSVA matrices
# ------------------------------------------------------------

write.csv(
    gse14520_gpl571_pathway_scores,
    file.path(
        "results",
        "Step14E_GSE14520_GPL571_17pathway_GSVA.csv"
    ),
    row.names = TRUE
)


write.csv(
    gse14520_gpl3921_pathway_scores,
    file.path(
        "results",
        "Step14E_GSE14520_GPL3921_17pathway_GSVA.csv"
    ),
    row.names = TRUE
)


# ------------------------------------------------------------
# Complete sample-level tables
# ------------------------------------------------------------

write.csv(
    gse14520_gpl571_complete_scores,
    file.path(
        "results",
        "Step14E_GSE14520_GPL571_complete_pathway_scores.csv"
    ),
    row.names = FALSE
)


write.csv(
    gse14520_gpl3921_complete_scores,
    file.path(
        "results",
        "Step14E_GSE14520_GPL3921_complete_pathway_scores.csv"
    ),
    row.names = FALSE
)


# ------------------------------------------------------------
# Primary scores
# ------------------------------------------------------------

write.csv(
    gse14520_gpl571_primary_scores,
    file.path(
        "results",
        "Step14E_GSE14520_GPL571_primary_scores.csv"
    ),
    row.names = FALSE
)


write.csv(
    gse14520_gpl3921_primary_scores,
    file.path(
        "results",
        "Step14E_GSE14520_GPL3921_primary_scores.csv"
    ),
    row.names = FALSE
)


# ------------------------------------------------------------
# DDR components
# ------------------------------------------------------------

write.csv(
    gse14520_gpl571_ddr_components,
    file.path(
        "results",
        "Step14E_GSE14520_GPL571_DDR_components.csv"
    ),
    row.names = TRUE
)


write.csv(
    gse14520_gpl3921_ddr_components,
    file.path(
        "results",
        "Step14E_GSE14520_GPL3921_DDR_components.csv"
    ),
    row.names = TRUE
)


# ------------------------------------------------------------
# Coverage
# ------------------------------------------------------------

write.csv(
    gse14520_gpl571_17pathway_coverage,
    file.path(
        "results",
        "Step14E_GSE14520_GPL571_17pathway_coverage.csv"
    ),
    row.names = FALSE
)


write.csv(
    gse14520_gpl3921_17pathway_coverage,
    file.path(
        "results",
        "Step14E_GSE14520_GPL3921_17pathway_coverage.csv"
    ),
    row.names = FALSE
)


# ------------------------------------------------------------
# QC
# ------------------------------------------------------------

write.csv(
    gse14520_score_qc,
    file.path(
        "results",
        "Step14E_GSE14520_score_QC.csv"
    ),
    row.names = FALSE
)


write.csv(
    gse14520_gpl571_pathway_summary,
    file.path(
        "results",
        "Step14E_GSE14520_GPL571_pathway_summary.csv"
    ),
    row.names = FALSE
)


write.csv(
    gse14520_gpl3921_pathway_summary,
    file.path(
        "results",
        "Step14E_GSE14520_GPL3921_pathway_summary.csv"
    ),
    row.names = FALSE
)


# ============================================================
# 14E.26 — SAVE COMPLETE RDATA CHECKPOINT
# ============================================================

cat("14E.26 — Saving complete RData checkpoint...\n\n")


save(
    gse14520_17_pathways,

    gpl571_expression_14e,
    gpl3921_expression_14e,

    gpl571_mapping_14e,
    gpl3921_mapping_14e,

    gpl571_gene_expression_14e,
    gpl3921_gene_expression_14e,

    gpl571_gene_sets_14e,
    gpl3921_gene_sets_14e,

    gse14520_gpl571_17pathway_coverage,
    gse14520_gpl3921_17pathway_coverage,

    gse14520_gpl571_pathway_scores,
    gse14520_gpl3921_pathway_scores,

    gse14520_gpl571_lipid_score,
    gse14520_gpl3921_lipid_score,

    gse14520_gpl571_hr_score,
    gse14520_gpl3921_hr_score,

    gse14520_gpl571_nhej_score,
    gse14520_gpl3921_nhej_score,

    gse14520_gpl571_ddr_components,
    gse14520_gpl3921_ddr_components,

    gse14520_gpl571_ddr_score,
    gse14520_gpl3921_ddr_score,

    gse14520_gpl571_p53_score,
    gse14520_gpl3921_p53_score,

    gse14520_gpl571_p53_g1_score,
    gse14520_gpl3921_p53_g1_score,

    gse14520_gpl571_g2m_score,
    gse14520_gpl3921_g2m_score,

    gse14520_gpl571_complete_scores,
    gse14520_gpl3921_complete_scores,

    gse14520_gpl571_primary_scores,
    gse14520_gpl3921_primary_scores,

    gse14520_score_qc,

    gse14520_gpl571_pathway_summary,
    gse14520_gpl3921_pathway_summary,

    file =
        file.path(
            "results",
            "Step14E_GSE14520_complete_pathway_analysis.RData"
        )
)


# ============================================================
# 14E.27 — CREATE TEXT AUDIT REPORT
# ============================================================

cat("14E.27 — Writing audit report...\n\n")


audit_file_14e <-
    file.path(
        "results",
        "Step14E_GSE14520_Step14E_audit.txt"
    )


audit_connection <-
    file(
        audit_file_14e,
        open = "wt"
    )


writeLines(
    c(
        "============================================================",
        "STEP 14E — GSE14520 PATHWAY SCORING AUDIT",
        "============================================================",
        "",
        "PATHWAY INVENTORY",
        "-----------------",
        paste(
            seq_along(gse14520_17_pathways),
            gse14520_17_pathways,
            sep = ". "
        ),
        "",
        "SAMPLE COUNTS",
        "-------------",
        paste0(
            "GPL571 tumors: ",
            ncol(
                gse14520_gpl571_pathway_scores
            )
        ),
        paste0(
            "GPL3921 tumors: ",
            ncol(
                gse14520_gpl3921_pathway_scores
            )
        ),
        "",
        "GSVA DIMENSIONS",
        "---------------",
        paste0(
            "GPL571: ",
            nrow(
                gse14520_gpl571_pathway_scores
            ),
            " pathways x ",
            ncol(
                gse14520_gpl571_pathway_scores
            ),
            " samples"
        ),
        paste0(
            "GPL3921: ",
            nrow(
                gse14520_gpl3921_pathway_scores
            ),
            " pathways x ",
            ncol(
                gse14520_gpl3921_pathway_scores
            ),
            " samples"
        ),
        "",
        "PRIMARY SCORING DEFINITIONS",
        "---------------------------",
        "Lipid score = HALLMARK_FATTY_ACID_METABOLISM",
        "HR score = HOMOLOGOUS_RECOMBINATION only",
        "NHEJ score = NON_HOMOLOGOUS_END_JOINING",
        paste0(
            "DDR score = mean of: ",
            paste(
                ddr_component_pathways_14e,
                collapse = " + "
            )
        ),
        "P53 is retained separately and is not part of DDR composite",
        "Expanded HR gene list is not used for pathway scoring",
        "",
        "NON-FINITE VALUES",
        "-----------------",
        paste0(
            "GPL571: ",
            gpl571_nonfinite_count_14e
        ),
        paste0(
            "GPL3921: ",
            gpl3921_nonfinite_count_14e
        ),
        "",
        "FILES CREATED",
        "-------------",
        "Step14E_GSE14520_GPL571_full_gene_expression.csv",
        "Step14E_GSE14520_GPL3921_full_gene_expression.csv",
        "Step14E_GSE14520_GPL571_17pathway_GSVA.csv",
        "Step14E_GSE14520_GPL3921_17pathway_GSVA.csv",
        "Step14E_GSE14520_GPL571_complete_pathway_scores.csv",
        "Step14E_GSE14520_GPL3921_complete_pathway_scores.csv",
        "Step14E_GSE14520_GPL571_primary_scores.csv",
        "Step14E_GSE14520_GPL3921_primary_scores.csv",
        "Step14E_GSE14520_GPL571_DDR_components.csv",
        "Step14E_GSE14520_GPL3921_DDR_components.csv",
        "Step14E_GSE14520_GPL571_17pathway_coverage.csv",
        "Step14E_GSE14520_GPL3921_17pathway_coverage.csv",
        "Step14E_GSE14520_score_QC.csv",
        "Step14E_GSE14520_GPL571_pathway_summary.csv",
        "Step14E_GSE14520_GPL3921_pathway_summary.csv",
        "Step14E_GSE14520_complete_pathway_analysis.RData",
        "Step14E_GSE14520_Step14E_audit.txt",
        "",
        "============================================================",
        "STEP 14E COMPLETE",
        "============================================================"
    ),
    audit_connection
)


close(
    audit_connection
)


# ============================================================
# 14E.28 — FINAL CONSOLE SUMMARY
# ============================================================

cat("\n")
cat("============================================================\n")
cat("FINAL STEP 14E SUMMARY\n")
cat("============================================================\n\n")


cat(
    "GPL571:\n",
    "  Tumors = ",
    ncol(
        gse14520_gpl571_pathway_scores
    ),
    "\n",
    "  Pathways = ",
    nrow(
        gse14520_gpl571_pathway_scores
    ),
    "\n\n",
    sep = ""
)


cat(
    "GPL3921:\n",
    "  Tumors = ",
    ncol(
        gse14520_gpl3921_pathway_scores
    ),
    "\n",
    "  Pathways = ",
    nrow(
        gse14520_gpl3921_pathway_scores
    ),
    "\n\n",
    sep = ""
)


cat(
    "Primary lipid score: HALLMARK_FATTY_ACID_METABOLISM\n"
)


cat(
    "HR score: HOMOLOGOUS_RECOMBINATION only\n"
)


cat(
    "DDR score: six-component repair composite\n"
)


cat(
    "NHEJ score: NON_HOMOLOGOUS_END_JOINING\n"
)


cat(
    "P53: retained separately\n"
)


cat(
    "Expanded HR list: reserved for gene-level analysis\n\n"
)


cat(
    "Non-finite GSVA values:\n"
)


cat(
    "  GPL571 = ",
    gpl571_nonfinite_count_14e,
    "\n",
    sep = ""
)


cat(
    "  GPL3921 = ",
    gpl3921_nonfinite_count_14e,
    "\n\n",
    sep = ""
)


cat(
    "Complete audit and all result files saved in results/.\n\n"
)


cat(
    "============================================================\n"
)


cat(
    "STEP 14E COMPLETE\n"
)


cat(
    "============================================================\n"
)
