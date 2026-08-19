# ============================================================
# STEP 15 TECHNICAL VALIDATION OF STEP 14
# ============================================================
#
# PURPOSE
# -------
# Determine whether genes classified as "missing" from
# GSE76427 are genuinely unavailable or whether some are
# represented under alternative gene symbols / aliases.
#
# This is a TECHNICAL QUALITY-CONTROL ANALYSIS.
#
# It does NOT alter the primary Step 14 results.
#
# PRIMARY ANALYSIS REMAINS:
#   77 TCGA genes
#       ↓
#   direct HGNC-symbol matching
#       ↓
#   31 genes available in GSE76427
#       ↓
#   gene-level validation
#
# This analysis only audits the 46 genes not directly matched.
#
# ============================================================


cat(
    "\n\n============================================================\n",
    "STEP 14 TECHNICAL VALIDATION\n",
    "============================================================\n\n",
    sep = ""
)


# ============================================================
# 14.T1 — CHECK REQUIRED OBJECTS
# ============================================================

required_objects <- c(
    "tcga14_genes",
    "gse14_expression",
    "genes_present",
    "genes_missing"
)


missing_objects <-
    required_objects[
        !sapply(
            required_objects,
            exists,
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
            "ERROR: The following Step 14 object(s) are missing:\n\n",
            paste(
                missing_objects,
                collapse = ", "
            ),
            "\n\n",
            "Run the main Step 14 GSE76427 analysis first."
        )
    )
}


# ============================================================
# 14.T2 — STANDARDIZE GENE SYMBOLS
# ============================================================

tcga14_genes_clean <-
    unique(
        toupper(
            trimws(
                as.character(
                    tcga14_genes
                )
            )
        )
    )


gse14_genes_clean <-
    unique(
        toupper(
            trimws(
                rownames(
                    gse14_expression
                )
            )
        )
    )


genes_missing_clean <-
    setdiff(
        tcga14_genes_clean,
        gse14_genes_clean
    )


cat(
    "TCGA genes: ",
    length(
        tcga14_genes_clean
    ),
    "\n",
    sep = ""
)

cat(
    "GSE76427 genes: ",
    length(
        gse14_genes_clean
    ),
    "\n",
    sep = ""
)

cat(
    "Directly matched genes: ",
    length(
        intersect(
            tcga14_genes_clean,
            gse14_genes_clean
        )
    ),
    "\n",
    sep = ""
)

cat(
    "Directly unmatched genes: ",
    length(
        genes_missing_clean
    ),
    "\n\n",
    sep = ""
)


# ============================================================
# 14.T3 — CHARACTERIZE MISSING GENES
# ============================================================

cat(
    "============================================================\n",
    "14.T3 — CHARACTERIZATION OF MISSING GENES\n",
    "============================================================\n\n"
)


technical_classification <-
    data.frame(
        Gene =
            genes_missing_clean,
        stringsAsFactors = FALSE
    )


# ------------------------------------------------------------
# Histone genes
# ------------------------------------------------------------

technical_classification$Histone_like <-
    grepl(
        "^H[234][A-Z0-9-]*",
        technical_classification$Gene
    ) |
    grepl(
        "^HIST",
        technical_classification$Gene
    )


# ------------------------------------------------------------
# Pseudogene-like identifiers
# ------------------------------------------------------------

technical_classification$Pseudogene_like <-
    grepl(
        "P$",
        technical_classification$Gene
    ) |
    grepl(
        "P[0-9]+$",
        technical_classification$Gene
    )


# ------------------------------------------------------------
# Read-through / unusual naming
# ------------------------------------------------------------

technical_classification$Unusual_identifier <-
    grepl(
        "-",
        technical_classification$Gene
    ) |
    grepl(
        "[0-9]+L$",
        technical_classification$Gene
    )


technical_classification$Initial_category <-
    "Other"


technical_classification$Initial_category[
    technical_classification$Histone_like
] <-
    "Histone-related"


technical_classification$Initial_category[
    !technical_classification$Histone_like &
    technical_classification$Pseudogene_like
] <-
    "Pseudogene-like"


technical_classification$Initial_category[
    !technical_classification$Histone_like &
    !technical_classification$Pseudogene_like &
    technical_classification$Unusual_identifier
] <-
    "Unusual identifier"


cat(
    "Missing-gene categories:\n\n"
)


print(
    table(
        technical_classification$Initial_category
    )
)


cat("\n")


# ============================================================
# 14.T4 — ATTEMPT COMMON SYMBOL ALIAS RESOLUTION
# ============================================================

cat(
    "============================================================\n",
    "14.T4 — SYMBOL / ALIAS CHECK\n",
    "============================================================\n\n"
)


# ------------------------------------------------------------
# Use HGNChelper if available.
#
# HGNChelper is specifically designed to identify outdated
# or invalid HGNC symbols.
# ------------------------------------------------------------

if (
    requireNamespace(
        "HGNChelper",
        quietly = TRUE
    )
) {

    cat(
        "HGNChelper package detected.\n\n"
    )


    alias_check <-
        HGNChelper::checkGeneSymbols(
            genes_missing_clean,
            unmappedAmbiguous = TRUE
        )


    alias_results <-
        data.frame(
            Original_Gene =
                as.character(
                    alias_check$x
                ),
            Approved_Gene =
                as.character(
                    alias_check$Suggested.Symbol
                ),
            Approved =
                as.logical(
                    alias_check$Approved
                ),
            stringsAsFactors = FALSE
        )


} else {

    cat(
        "HGNChelper is not installed.\n",
        "Skipping automated HGNC alias resolution.\n\n",
        sep = ""
    )


    alias_results <-
        data.frame(
            Original_Gene =
                genes_missing_clean,
            Approved_Gene =
                NA_character_,
            Approved =
                NA,
            stringsAsFactors = FALSE
        )
}


# ============================================================
# 14.T5 — TEST WHETHER SUGGESTED ALIASES EXIST IN GSE
# ============================================================

alias_results$Alias_present_in_GSE <-
    FALSE


valid_alias <-
    !is.na(
        alias_results$Approved_Gene
    ) &
    alias_results$Approved_Gene != ""


alias_results$Alias_present_in_GSE[
    valid_alias
] <-
    alias_results$Approved_Gene[
        valid_alias
    ] %in%
    gse14_genes_clean


# ------------------------------------------------------------
# Recoverable aliases
# ------------------------------------------------------------

recoverable_aliases <-
    alias_results[
        alias_results$Alias_present_in_GSE,
        ,
        drop = FALSE
    ]


cat(
    "Genes potentially recoverable through HGNC aliases: ",
    nrow(
        recoverable_aliases
    ),
    "\n\n",
    sep = ""
)


if (
    nrow(
        recoverable_aliases
    ) > 0
) {

    print(
        recoverable_aliases
    )

} else {

    cat(
        "No additional genes were recovered through HGNC ",
        "alias resolution.\n",
        sep = ""
    )
}


# ============================================================
# 14.T6 — CALCULATE TECHNICAL COVERAGE
# ============================================================

directly_recovered <-
    length(
        intersect(
            tcga14_genes_clean,
            gse14_genes_clean
        )
    )


alias_recovered <-
    nrow(
        recoverable_aliases
    )


total_recovered_after_alias <-
    directly_recovered +
    alias_recovered


coverage_direct <-
    100 *
    directly_recovered /
    length(
        tcga14_genes_clean
    )


coverage_after_alias <-
    100 *
    total_recovered_after_alias /
    length(
        tcga14_genes_clean
    )


cat(
    "\n============================================================\n",
    "14.T6 — COVERAGE SUMMARY\n",
    "============================================================\n\n"
)


cat(
    "Total TCGA DSB genes: ",
    length(
        tcga14_genes_clean
    ),
    "\n",
    sep = ""
)

cat(
    "Directly matched in GSE: ",
    directly_recovered,
    "\n",
    sep = ""
)

cat(
    "Direct coverage: ",
    round(
        coverage_direct,
        1
    ),
    "%\n",
    sep = ""
)

cat(
    "Additional genes recovered by aliases: ",
    alias_recovered,
    "\n",
    sep = ""
)

cat(
    "Total coverage after alias check: ",
    round(
        coverage_after_alias,
        1
    ),
    "%\n\n",
    sep = ""
)


# ============================================================
# 14.T7 — SAVE TECHNICAL AUDIT
# ============================================================

write.csv(
    technical_classification,
    file =
        file.path(
            "results",
            "Step14_technical_missing_gene_classification.csv"
        ),
    row.names = FALSE
)


write.csv(
    alias_results,
    file =
        file.path(
            "results",
            "Step14_technical_HGNC_alias_audit.csv"
        ),
    row.names = FALSE
)


technical_summary <-
    data.frame(
        Metric = c(
            "Total TCGA DSB genes",
            "Directly matched genes",
            "Direct coverage (%)",
            "Directly missing genes",
            "Alias-recoverable genes",
            "Coverage after alias check (%)"
        ),
        Value = c(
            length(
                tcga14_genes_clean
            ),
            directly_recovered,
            round(
                coverage_direct,
                2
            ),
            length(
                genes_missing_clean
            ),
            alias_recovered,
            round(
                coverage_after_alias,
                2
            )
        ),
        stringsAsFactors = FALSE
    )


write.csv(
    technical_summary,
    file =
        file.path(
            "results",
            "Step14_technical_validation_summary.csv"
        ),
    row.names = FALSE
)


# ============================================================
# 14.T8 — SAVE RDATA CHECKPOINT
# ============================================================

save(
    technical_classification,
    alias_results,
    technical_summary,
    file =
        file.path(
            "results",
            "Step14_technical_validation.RData"
        )
)


# ============================================================
# 14.T9 — FINAL INTERPRETATION OUTPUT
# ============================================================

cat(
    "\n============================================================\n",
    "STEP 14 TECHNICAL VALIDATION COMPLETE\n",
    "============================================================\n\n"
)


cat(
    "Direct GSE coverage: ",
    round(
        coverage_direct,
        1
    ),
    "%\n",
    sep = ""
)

cat(
    "Additional alias-recovered genes: ",
    alias_recovered,
    "\n",
    sep = ""
)

cat(
    "Coverage after alias check: ",
    round(
        coverage_after_alias,
        1
    ),
    "%\n\n",
    sep = ""
)


if (
    alias_recovered == 0
) {

    cat(
        "INTERPRETATION:\n",
        "No additional GSE76427 genes were recovered through ",
        "HGNC alias resolution.\n\n",
        "Therefore, the 40.3% coverage observed in the primary ",
        "analysis is unlikely to be explained by simple gene-",
        "symbol nomenclature differences.\n",
        sep = ""
    )

} else {

    cat(
        "INTERPRETATION:\n",
        alias_recovered,
        "additional genes may be represented in GSE76427 ",
        "under alternative HGNC nomenclature.\n\n",
        "These genes should be reviewed manually before deciding ",
        "whether they should be incorporated into a sensitivity ",
        "analysis.\n",
        sep = ""
    )
}


cat(
    "\n============================================================\n",
    "END OF STEP 14 TECHNICAL VALIDATION\n",
    "============================================================\n"
)
