# ============================================================
# STEP 8 — COMMON GENE UNIVERSE AND GENE-SET HARMONIZATION
# MAFLD-HCC HCC MANUSCRIPT PIPELINE
# ============================================================

# INPUTS
#   data/processed/TCGA-LIHC/
#       TCGA_LIHC_voom_symbol_expression.rds
#
#   data/processed/GSE76427/
#       GSE76427_tumor_symbol_expression.rds
#
#   data/reference/gene_sets/
#       *.gmt
#
# OUTPUTS
#   results/08_TCGA_Ensembl_to_symbol_full_mapping.csv
#   results/08_TCGA_Ensembl_to_symbol_mapping_unique.csv
#   results/08_common_genes.csv
#   results/08_common_gene_universe_summary.csv
#   results/08_gene_set_file_inventory.csv
#
# ALSO SAVES
#   data/processed/TCGA-LIHC/
#       TCGA_LIHC_voom_symbol_expression_harmonized.rds
#
#   data/processed/GSE76427/
#       GSE76427_tumor_symbol_expression_harmonized.rds
#
#   data/processed/
#       common_gene_universe.rds
#
# This step does NOT perform GSVA.
# ============================================================


# ------------------------------------------------------------
# 8.1 SET PROJECT ROOT
# ------------------------------------------------------------

project_root <- "C:/Users/Dhruv Advani/OneDrive/Desktop/AIIMS/Research/MAFLD_HCC_Pipeline"

if (!dir.exists(project_root)) {
    stop(
        "Project root does not exist:\n",
        project_root
    )
}

setwd(project_root)

cat("\n")
cat("============================================================\n")
cat("STEP 8 — COMMON GENE UNIVERSE AND GENE-SET HARMONIZATION\n")
cat("============================================================\n\n")

cat("Project root:\n")
cat(getwd(), "\n\n")


# ------------------------------------------------------------
# 8.2 DEFINE DIRECTORIES
# ------------------------------------------------------------

processed_dir <- "data/processed"
reference_dir <- "data/reference"
results_dir <- "results"


# ------------------------------------------------------------
# 8.3 DEFINE INPUT FILES
# ------------------------------------------------------------

tcga_file <- file.path(
    processed_dir,
    "TCGA-LIHC",
    "TCGA_LIHC_voom_symbol_expression.rds"
)

gse_file <- file.path(
    processed_dir,
    "GSE76427",
    "GSE76427_tumor_symbol_expression.rds"
)

gene_set_dir <- file.path(
    reference_dir,
    "gene_sets"
)


# ------------------------------------------------------------
# 8.4 CHECK INPUTS
# ------------------------------------------------------------

if (!file.exists(tcga_file)) {
    stop(
        "TCGA input file not found:\n",
        tcga_file
    )
}

if (!file.exists(gse_file)) {
    stop(
        "GSE76427 input file not found:\n",
        gse_file
    )
}

if (!dir.exists(gene_set_dir)) {
    stop(
        "Gene-set directory not found:\n",
        gene_set_dir
    )
}


# ------------------------------------------------------------
# 8.5 FUNCTION TO EXTRACT EXPRESSION MATRIX
# ------------------------------------------------------------

extract_expression_matrix <- function(
    object,
    object_name
) {

    # Case 1: matrix
    if (is.matrix(object)) {

        expr <- object

    # Case 2: data.frame
    } else if (is.data.frame(object)) {

        expr <- as.matrix(object)

    # Case 3: limma EList / voom object
    } else if (
        is.list(object) &&
        !is.null(object$E)
    ) {

        expr <- object$E

    # Case 4: object containing exprs
    } else if (
        is.list(object) &&
        !is.null(object$exprs)
    ) {

        expr <- object$exprs

    # Case 5: object containing expression
    } else if (
        is.list(object) &&
        !is.null(object$expression)
    ) {

        expr <- object$expression

    # Case 6: object containing counts
    } else if (
        is.list(object) &&
        !is.null(object$counts)
    ) {

        expr <- object$counts

    # Case 7: SummarizedExperiment
    } else if (
        methods::is(
            object,
            "SummarizedExperiment"
        )
    ) {

        if (!requireNamespace(
            "SummarizedExperiment",
            quietly = TRUE
        )) {
            stop(
                "SummarizedExperiment package is required."
            )
        }

        expr <- SummarizedExperiment::assay(
            object
        )

    # Case 8: ExpressionSet
    } else if (
        methods::is(
            object,
            "ExpressionSet"
        )
    ) {

        if (!requireNamespace(
            "Biobase",
            quietly = TRUE
        )) {
            stop(
                "Biobase package is required."
            )
        }

        expr <- Biobase::exprs(
            object
        )

    } else {

        cat("\nUnable to identify expression matrix in:\n")
        cat(object_name, "\n\n")

        cat("Object class:\n")
        print(class(object))

        if (is.list(object)) {
            cat("\nObject components:\n")
            print(names(object))
        }

        stop(
            "\nCannot safely extract expression matrix from ",
            object_name,
            "."
        )
    }


    # Convert data.frame to matrix if necessary
    if (is.data.frame(expr)) {
        expr <- as.matrix(expr)
    }


    # Validate
    if (!is.matrix(expr)) {
        stop(
            object_name,
            ": extracted expression component is not a matrix."
        )
    }

    if (nrow(expr) == 0 || ncol(expr) == 0) {
        stop(
            object_name,
            ": expression matrix is empty."
        )
    }

    if (is.null(rownames(expr))) {
        stop(
            object_name,
            ": expression matrix has no gene identifiers."
        )
    }

    if (is.null(colnames(expr))) {
        stop(
            object_name,
            ": expression matrix has no sample identifiers."
        )
    }

    # Ensure numeric matrix
    if (!is.numeric(expr)) {

        suppressWarnings(
            expr_numeric <- matrix(
                as.numeric(expr),
                nrow = nrow(expr),
                ncol = ncol(expr),
                dimnames = dimnames(expr)
            )
        )

        if (
            all(is.na(expr_numeric)) &&
            !all(is.na(expr))
        ) {
            stop(
                object_name,
                ": expression values could not be converted to numeric."
            )
        }

        expr <- expr_numeric
    }

    return(expr)
}


# ------------------------------------------------------------
# 8.6 LOAD TCGA OBJECT
# ------------------------------------------------------------

cat("Loading TCGA expression object...\n")

tcga_object <- readRDS(
    tcga_file
)

cat("TCGA object class:\n")
print(class(tcga_object))

cat("\n")


# ------------------------------------------------------------
# 8.7 EXTRACT TCGA MATRIX
# ------------------------------------------------------------

tcga_expr <- extract_expression_matrix(
    tcga_object,
    "TCGA_LIHC_voom_symbol_expression.rds"
)

cat(
    "TCGA expression matrix:",
    nrow(tcga_expr),
    "genes x",
    ncol(tcga_expr),
    "samples\n\n"
)


# ------------------------------------------------------------
# 8.8 LOAD GSE76427 OBJECT
# ------------------------------------------------------------

cat("Loading GSE76427 expression object...\n")

gse_object <- readRDS(
    gse_file
)

cat("GSE76427 object class:\n")
print(class(gse_object))

cat("\n")


# ------------------------------------------------------------
# 8.9 EXTRACT GSE MATRIX
# ------------------------------------------------------------

gse_expr <- extract_expression_matrix(
    gse_object,
    "GSE76427_tumor_symbol_expression.rds"
)

cat(
    "GSE76427 expression matrix:",
    nrow(gse_expr),
    "genes x",
    ncol(gse_expr),
    "samples\n\n"
)


# ------------------------------------------------------------
# 8.10 CLEAN GENE IDENTIFIERS
# ------------------------------------------------------------

tcga_ids <- trimws(
    as.character(
        rownames(tcga_expr)
    )
)

gse_ids <- trimws(
    as.character(
        rownames(gse_expr)
    )
)

tcga_ids <- toupper(
    tcga_ids
)

gse_ids <- toupper(
    gse_ids
)


# Remove Ensembl version suffix if present
tcga_ids_clean <- sub(
    "\\..*$",
    "",
    tcga_ids
)

tcga_ids_clean <- toupper(
    tcga_ids_clean
)


# ------------------------------------------------------------
# 8.11 DETERMINE TCGA IDENTIFIER TYPE
# ------------------------------------------------------------

ensembl_fraction <- mean(
    grepl(
        "^ENSG[0-9]+$",
        tcga_ids_clean
    )
)

tcga_is_ensembl <- (
    ensembl_fraction >= 0.50
)

cat("------------------------------------------------------------\n")
cat("TCGA GENE IDENTIFIER CHECK\n")
cat("------------------------------------------------------------\n")

cat(
    "Fraction resembling Ensembl IDs:",
    round(
        ensembl_fraction,
        3
    ),
    "\n"
)

if (tcga_is_ensembl) {

    cat(
        "TCGA identifiers will be mapped from Ensembl to HGNC symbols.\n\n"
    )

} else {

    cat(
        "TCGA identifiers already appear to be gene symbols.\n\n"
    )
}


# ------------------------------------------------------------
# 8.12 CREATE TCGA MAPPING
# ------------------------------------------------------------

if (tcga_is_ensembl) {

    if (!requireNamespace(
        "AnnotationDbi",
        quietly = TRUE
    )) {
        stop(
            "Package 'AnnotationDbi' is required but not installed."
        )
    }

    if (!requireNamespace(
        "org.Hs.eg.db",
        quietly = TRUE
    )) {
        stop(
            "Package 'org.Hs.eg.db' is required but not installed."
        )
    }

    cat(
        "Mapping Ensembl IDs to HGNC symbols...\n"
    )

    tcga_mapping_full <- AnnotationDbi::select(
        org.Hs.eg.db::org.Hs.eg.db,
        keys = unique(tcga_ids_clean),
        keytype = "ENSEMBL",
        columns = c(
            "ENSEMBL",
            "SYMBOL"
        )
    )

    tcga_mapping_full <- as.data.frame(
        tcga_mapping_full,
        stringsAsFactors = FALSE
    )

    tcga_mapping_full$ENSEMBL <- toupper(
        trimws(
            tcga_mapping_full$ENSEMBL
        )
    )

    tcga_mapping_full$SYMBOL <- toupper(
        trimws(
            tcga_mapping_full$SYMBOL
        )
    )

    tcga_mapping_full$SYMBOL[
        tcga_mapping_full$SYMBOL == ""
    ] <- NA_character_

} else {

    tcga_mapping_full <- data.frame(
        ENSEMBL = tcga_ids_clean,
        SYMBOL = tcga_ids_clean,
        stringsAsFactors = FALSE
    )
}


# ------------------------------------------------------------
# 8.13 SAVE FULL TCGA MAPPING
# ------------------------------------------------------------

write.csv(
    tcga_mapping_full,
    file = file.path(
        results_dir,
        "08_TCGA_Ensembl_to_symbol_full_mapping.csv"
    ),
    row.names = FALSE
)


# ------------------------------------------------------------
# 8.14 KEEP VALID MAPPINGS
# ------------------------------------------------------------

tcga_mapping_valid <- tcga_mapping_full[
    !is.na(
        tcga_mapping_full$SYMBOL
    ) &
    tcga_mapping_full$SYMBOL != "",
    ,
    drop = FALSE
]


# ------------------------------------------------------------
# 8.15 REMOVE DUPLICATE SYMBOL MAPPINGS
# ------------------------------------------------------------

tcga_mapping_unique <- tcga_mapping_valid[
    !duplicated(
        tcga_mapping_valid$SYMBOL
    ),
    ,
    drop = FALSE
]


write.csv(
    tcga_mapping_unique,
    file = file.path(
        results_dir,
        "08_TCGA_Ensembl_to_symbol_mapping_unique.csv"
    ),
    row.names = FALSE
)


# ------------------------------------------------------------
# 8.16 CONVERT TCGA MATRIX TO SYMBOL LEVEL
# ------------------------------------------------------------

if (tcga_is_ensembl) {

    mapping_index <- match(
        tcga_ids_clean,
        tcga_mapping_unique$ENSEMBL
    )

    keep_rows <- !is.na(
        mapping_index
    )

    tcga_expr_symbol <- tcga_expr[
        keep_rows,
        ,
        drop = FALSE
    ]

    rownames(
        tcga_expr_symbol
    ) <- tcga_mapping_unique$SYMBOL[
        mapping_index[
            keep_rows
        ]
    ]

} else {

    tcga_expr_symbol <- tcga_expr

    rownames(
        tcga_expr_symbol
    ) <- tcga_ids_clean
}


# ------------------------------------------------------------
# 8.17 REMOVE DUPLICATE TCGA SYMBOLS
# ------------------------------------------------------------

if (
    anyDuplicated(
        rownames(tcga_expr_symbol)
    ) > 0
) {

    tcga_expr_symbol <- tcga_expr_symbol[
        !duplicated(
            rownames(tcga_expr_symbol)
        ),
        ,
        drop = FALSE
    ]
}


# ------------------------------------------------------------
# 8.18 CLEAN GSE76427 MATRIX
# ------------------------------------------------------------

rownames(
    gse_expr
) <- gse_ids

valid_gse <- !is.na(
    rownames(gse_expr)
) &
    rownames(gse_expr) != ""

gse_expr <- gse_expr[
    valid_gse,
    ,
    drop = FALSE
]


# Remove duplicate symbols
gse_expr <- gse_expr[
    !duplicated(
        rownames(gse_expr)
    ),
    ,
    drop = FALSE
]


# ------------------------------------------------------------
# 8.19 CREATE FINAL GENE LISTS
# ------------------------------------------------------------

tcga_genes <- unique(
    rownames(
        tcga_expr_symbol
    )
)

gse_genes <- unique(
    rownames(
        gse_expr
    )
)


# ------------------------------------------------------------
# 8.20 COMMON GENE UNIVERSE
# ------------------------------------------------------------

common_genes <- intersect(
    tcga_genes,
    gse_genes
)

common_genes <- sort(
    unique(
        common_genes
    )
)


if (length(common_genes) == 0) {

    stop(
        paste(
            "No common genes were found between TCGA and",
            "GSE76427.",
            "Check the identifier formats."
        )
    )
}


# ------------------------------------------------------------
# 8.21 COMMON GENE TABLE
# ------------------------------------------------------------

common_gene_table <- data.frame(
    gene_symbol = common_genes,
    in_TCGA = TRUE,
    in_GSE76427 = TRUE,
    stringsAsFactors = FALSE
)

write.csv(
    common_gene_table,
    file = file.path(
        results_dir,
        "08_common_genes.csv"
    ),
    row.names = FALSE
)


# ------------------------------------------------------------
# 8.22 SCAN ALL GMT FILES
# ------------------------------------------------------------

cat("\n")
cat("Scanning gene-set reference directory...\n")

gmt_files <- list.files(
    gene_set_dir,
    pattern = "\\.gmt$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
)

if (length(gmt_files) == 0) {

    stop(
        "No GMT files found in:\n",
        gene_set_dir
    )
}

cat(
    "GMT files detected:",
    length(gmt_files),
    "\n\n"
)


# ------------------------------------------------------------
# 8.23 INVENTORY GMT FILES (per gene set / pathway)
# ------------------------------------------------------------

gene_set_inventory <- list()

for (gmt_file in gmt_files) {

    lines <- readLines(
        gmt_file,
        warn = FALSE
    )

    lines <- lines[
        nzchar(
            trimws(lines)
        )
    ]

    if (length(lines) == 0) {
        next
    }

    for (line in lines) {

        fields <- strsplit(
            line,
            "\t",
            fixed = TRUE
        )[[1]]

        if (length(fields) < 3) {
            next
        }

        pathway_name <- fields[1]

        genes <- toupper(
            trimws(
                fields[-c(1, 2)]
            )
        )

        genes <- unique(
            genes[
                !is.na(genes) &
                genes != ""
            ]
        )

        n_genes_total <- length(
            genes
        )

        n_genes_common <- length(
            intersect(
                genes,
                common_genes
            )
        )

        pct_covered <- if (
            n_genes_total > 0
        ) {
            round(
                100 *
                n_genes_common /
                n_genes_total,
                2
            )
        } else {
            NA_real_
        }

        inventory_row <- data.frame(
            gmt_file = basename(gmt_file),
            pathway_name = pathway_name,
            n_genes_total = n_genes_total,
            n_genes_common = n_genes_common,
            pct_covered = pct_covered,
            stringsAsFactors = FALSE
        )

        gene_set_inventory <- c(
            gene_set_inventory,
            list(inventory_row)
        )
    }
}


# ------------------------------------------------------------
# 8.24 COMBINE GENE-SET INVENTORY
# ------------------------------------------------------------

if (length(gene_set_inventory) == 0) {

    stop(
        "No valid gene sets were found in the GMT files."
    )
}

gene_set_inventory_df <- do.call(
    rbind,
    gene_set_inventory
)

write.csv(
    gene_set_inventory_df,
    file = file.path(
        results_dir,
        "08_gene_set_file_inventory.csv"
    ),
    row.names = FALSE
)

cat(
    "Gene sets inventoried:",
    nrow(gene_set_inventory_df),
    "\n",
    "Pathways below 50% coverage in common universe:",
    sum(
        gene_set_inventory_df$pct_covered < 50,
        na.rm = TRUE
    ),
    "\n\n"
)


# ------------------------------------------------------------
# 8.25 SUMMARY
# ------------------------------------------------------------

n_tcga_input <- length(
    unique(
        tcga_ids_clean
    )
)

n_tcga_mapped <- sum(
    !is.na(
        tcga_mapping_full$SYMBOL
    )
)

n_tcga_genes <- length(
    tcga_genes
)

n_gse_genes <- length(
    gse_genes
)

n_common_genes <- length(
    common_genes
)

n_gmt_files <- length(
    gmt_files
)

n_gene_sets <- nrow(
    gene_set_inventory_df
)

tcga_overlap <- round(
    100 *
    n_common_genes /
    n_tcga_genes,
    2
)

gse_overlap <- round(
    100 *
    n_common_genes /
    n_gse_genes,
    2
)


summary_table <- data.frame(
    metric = c(
        "TCGA input gene identifiers",
        "TCGA identifiers successfully mapped",
        "TCGA unique gene symbols",
        "GSE76427 unique gene symbols",
        "Common genes",
        "Common genes as percentage of TCGA",
        "Common genes as percentage of GSE76427",
        "GMT files detected",
        "Total gene sets detected"
    ),
    value = c(
        n_tcga_input,
        n_tcga_mapped,
        n_tcga_genes,
        n_gse_genes,
        n_common_genes,
        tcga_overlap,
        gse_overlap,
        n_gmt_files,
        n_gene_sets
    ),
    stringsAsFactors = FALSE
)

write.csv(
    summary_table,
    file = file.path(
        results_dir,
        "08_common_gene_universe_summary.csv"
    ),
    row.names = FALSE
)


# ------------------------------------------------------------
# 8.26 SAVE HARMONIZED MATRICES
# ------------------------------------------------------------

saveRDS(
    tcga_expr_symbol,
    file.path(
        processed_dir,
        "TCGA-LIHC",
        "TCGA_LIHC_voom_symbol_expression_harmonized.rds"
    )
)

saveRDS(
    gse_expr,
    file.path(
        processed_dir,
        "GSE76427",
        "GSE76427_tumor_symbol_expression_harmonized.rds"
    )
)

saveRDS(
    common_genes,
    file.path(
        processed_dir,
        "common_gene_universe.rds"
    )
)


# ------------------------------------------------------------
# 8.27 FINAL QC
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("STEP 8 QC SUMMARY\n")
cat("============================================================\n")

cat(
    "TCGA input genes: ",
    n_tcga_input,
    "\n",
    sep = ""
)

cat(
    "TCGA mapped genes: ",
    n_tcga_mapped,
    "\n",
    sep = ""
)

cat(
    "TCGA unique symbols: ",
    n_tcga_genes,
    "\n",
    sep = ""
)

cat(
    "GSE76427 unique symbols: ",
    n_gse_genes,
    "\n",
    sep = ""
)

cat(
    "Common genes: ",
    n_common_genes,
    "\n",
    sep = ""
)

cat(
    "Common / TCGA: ",
    tcga_overlap,
    "%\n",
    sep = ""
)

cat(
    "Common / GSE76427: ",
    gse_overlap,
    "%\n",
    sep = ""
)

cat(
    "GMT files: ",
    n_gmt_files,
    "\n",
    sep = ""
)

cat(
    "Total gene sets: ",
    n_gene_sets,
    "\n",
    sep = ""
)


# ------------------------------------------------------------
# 8.28 HARD QC CHECKS
# ------------------------------------------------------------

if (n_tcga_genes < 1000) {

    stop(
        "QC FAILURE: fewer than 1,000 TCGA genes remain."
    )
}

if (n_gse_genes < 1000) {

    stop(
        "QC FAILURE: fewer than 1,000 GSE76427 genes remain."
    )
}

if (n_common_genes < 1000) {

    stop(
        "QC FAILURE: fewer than 1,000 common genes remain."
    )
}

if (n_gmt_files == 0) {

    stop(
        "QC FAILURE: no GMT files detected."
    )


}


# ------------------------------------------------------------
# 8.29 FINAL COMPLETION MESSAGE
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("STEP 8 COMPLETE — QC PASSED\n")
cat("============================================================\n")
cat("\n")

cat(
    "Step 8 outputs have been written to the results/ directory.\n"
)

cat(
    "The two pathways with <50% common-universe coverage are",
    "retained and will be evaluated by effective gene-set size",
    "during the downstream GSVA QC step.\n"
)
