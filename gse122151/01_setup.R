# ============================================================
# GSE122151 — STEP 1: SETUP
# MAFLD-HCC manuscript
# ============================================================

rm(list = ls())
gc()

options(
  stringsAsFactors = FALSE,
  scipen = 999,
  warn = 1
)

# ------------------------------------------------------------
# Project root
# ------------------------------------------------------------

project_root <- getwd()

cat("Project root:\n", project_root, "\n\n")

# ------------------------------------------------------------
# Required directories
# ------------------------------------------------------------

required_dirs <- c(
  "config",
  "data/raw/GSE122151",
  "data/processed/GSE122151",
  "results/GSE122151",
  "figures/GSE122151",
  "scripts"
)

for (d in required_dirs) {
  dir.create(
    file.path(project_root, d),
    recursive = TRUE,
    showWarnings = FALSE
  )
}

# ------------------------------------------------------------
# Required CRAN packages
# ------------------------------------------------------------

cran_packages <- c(
  "data.table",
  "dplyr",
  "tidyr",
  "tibble",
  "readr",
  "stringr",
  "ggplot2",
  "pheatmap",
  "openxlsx",
  "matrixStats",
  "RColorBrewer"
)

installed <- rownames(installed.packages())

for (pkg in cran_packages) {

  if (!pkg %in% installed) {

    install.packages(
      pkg,
      repos = "https://cloud.r-project.org"
    )
  }
}

# ------------------------------------------------------------
# Required Bioconductor packages
# ------------------------------------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE)) {

  install.packages(
    "BiocManager",
    repos = "https://cloud.r-project.org"
  )
}

bioc_packages <- c(
  "GEOquery",
  "GSVA",
  "limma",
  "edgeR",
  "fgsea"
)

for (pkg in bioc_packages) {

  if (!requireNamespace(pkg, quietly = TRUE)) {

    BiocManager::install(
      pkg,
      ask = FALSE,
      update = FALSE
    )
  }
}

# ------------------------------------------------------------
# msigdbr
# ------------------------------------------------------------

if (!requireNamespace("msigdbr", quietly = TRUE)) {

  install.packages(
    "msigdbr",
    repos = "https://cloud.r-project.org"
  )
}

# ------------------------------------------------------------
# Load packages
# ------------------------------------------------------------

suppressPackageStartupMessages({

  library(data.table)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(pheatmap)
  library(openxlsx)
  library(matrixStats)
  library(RColorBrewer)

  library(GEOquery)
  library(GSVA)
  library(limma)
  library(edgeR)
  library(fgsea)

  library(msigdbr)
})

# ------------------------------------------------------------
# Save package versions
# ------------------------------------------------------------

package_versions <- data.frame(
  package = c(
    cran_packages,
    "BiocManager",
    bioc_packages,
    "msigdbr"
  ),
  version = sapply(
    c(
      cran_packages,
      "BiocManager",
      bioc_packages,
      "msigdbr"
    ),
    function(x) as.character(packageVersion(x))
  ),
  stringsAsFactors = FALSE
)

write.csv(
  package_versions,
  file.path(
    project_root,
    "results/GSE122151/package_versions.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# Session information
# ------------------------------------------------------------

capture.output(
  sessionInfo(),
  file = file.path(
    project_root,
    "results/GSE122151/sessionInfo.txt"
  )
)

cat("\nSETUP COMPLETE\n")
