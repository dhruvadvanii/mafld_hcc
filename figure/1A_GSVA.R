# ============================================================
# FIGURE 1A — HOMOLOGOUS RECOMBINATION
# TCGA + GSE76427
# ============================================================

library(ggplot2)
library(dplyr)

# ------------------------------------------------------------
# 1. TCGA DATA
# ------------------------------------------------------------

tcga_hr <- data.frame(
  Sample = colnames(tcga_pathway_scores),
  HR = as.numeric(
    tcga_pathway_scores["HOMOLOGOUS_RECOMBINATION", ]
  )
)

# Match TCGA samples to lipid phenotype
tcga_hr$Lipid_Group <- as.character(
  annotation_col[
    match(
      tcga_hr$Sample,
      rownames(annotation_col)
    ),
    "LipidGroup"
  ]
)

tcga_hr$Dataset <- "TCGA"


# ------------------------------------------------------------
# 2. GSE76427 DATA
# ------------------------------------------------------------

gse_hr <- data.frame(
  Sample = colnames(gse_pathway_scores),
  HR = as.numeric(
    gse_pathway_scores["HOMOLOGOUS_RECOMBINATION", ]
  )
)

# Existing GSE76427 lipid grouping
gse_hr$Lipid_Group <- as.character(gse12_lipid_group)

gse_hr$Dataset <- "GSE76427"


# ------------------------------------------------------------
# 3. STANDARDIZE GROUPS
# ------------------------------------------------------------

tcga_hr$Lipid_Group <- factor(
  tcga_hr$Lipid_Group,
  levels = c("Lipid-low", "Lipid-high")
)

gse_hr$Lipid_Group <- factor(
  gse_hr$Lipid_Group,
  levels = c("Lipid-low", "Lipid-high")
)


# ------------------------------------------------------------
# 4. REMOVE MISSING VALUES
# ------------------------------------------------------------

tcga_hr <- tcga_hr[
  !is.na(tcga_hr$HR) &
  !is.na(tcga_hr$Lipid_Group),
]

gse_hr <- gse_hr[
  !is.na(gse_hr$HR) &
  !is.na(gse_hr$Lipid_Group),
]


# ------------------------------------------------------------
# 5. COMBINE DATASETS
# ------------------------------------------------------------

hr_plot_data <- bind_rows(
  tcga_hr,
  gse_hr
)

hr_plot_data$Dataset <- factor(
  hr_plot_data$Dataset,
  levels = c("TCGA", "GSE76427")
)


# ------------------------------------------------------------
# 6. CREATE DISPLAY GROUPS
# ------------------------------------------------------------

hr_plot_data$Group <- interaction(
  hr_plot_data$Dataset,
  hr_plot_data$Lipid_Group,
  sep = " | ",
  drop = TRUE
)

hr_plot_data$Group <- factor(
  hr_plot_data$Group,
  levels = c(
    "TCGA | Lipid-low",
    "TCGA | Lipid-high",
    "GSE76427 | Lipid-low",
    "GSE76427 | Lipid-high"
  )
)


# ------------------------------------------------------------
# 7. COLOURS
# ------------------------------------------------------------

dataset_cols <- c(
  "TCGA" = "#3B6FB6",
  "GSE76427" = "#D8893C"
)


# ------------------------------------------------------------
# 8. Y-AXIS RANGE
# ------------------------------------------------------------

y_min <- min(
  hr_plot_data$HR,
  na.rm = TRUE
)

y_max <- max(
  hr_plot_data$HR,
  na.rm = TRUE
)

y_range <- y_max - y_min

if (y_range == 0) {
  y_range <- 1
}

tcga_bar_y <- y_max + 0.10 * y_range
gse_bar_y  <- y_max + 0.23 * y_range


# ------------------------------------------------------------
# 9. CREATE FIGURE
# ------------------------------------------------------------

p_hr <- ggplot(
  hr_plot_data,
  aes(
    x = Group,
    y = HR,
    fill = Dataset
  )
) +

  geom_boxplot(
    width = 0.62,
    outlier.shape = NA,
    linewidth = 0.55,
    alpha = 0.78
  ) +

  geom_jitter(
    aes(color = Dataset),
    width = 0.10,
    size = 1.35,
    alpha = 0.40
  ) +

  scale_fill_manual(
    values = dataset_cols
  ) +

  scale_color_manual(
    values = dataset_cols
  ) +

  scale_x_discrete(
    labels = c(
      "TCGA\nLipid-low",
      "TCGA\nLipid-high",
      "GSE76427\nLipid-low",
      "GSE76427\nLipid-high"
    )
  ) +

  labs(
    title = "Homologous recombination",
    x = NULL,
    y = "GSVA enrichment score"
  ) +

  # ----------------------------------------------------------
  # TCGA SIGNIFICANCE
  # FDR = 7.01 × 10^-12
  # ----------------------------------------------------------

  annotate(
    "segment",
    x = 1,
    xend = 2,
    y = tcga_bar_y,
    yend = tcga_bar_y,
    linewidth = 0.45
  ) +

  annotate(
    "segment",
    x = 1,
    xend = 1,
    y = tcga_bar_y - 0.025 * y_range,
    yend = tcga_bar_y,
    linewidth = 0.45
  ) +

  annotate(
    "segment",
    x = 2,
    xend = 2,
    y = tcga_bar_y - 0.025 * y_range,
    yend = tcga_bar_y,
    linewidth = 0.45
  ) +

  annotate(
    "text",
    x = 1.5,
    y = tcga_bar_y + 0.025 * y_range,
    label = "****",
    size = 4.5,
    fontface = "bold"
  ) +

  # ----------------------------------------------------------
  # GSE76427 SIGNIFICANCE
  # FDR = 1.29 × 10^-2
  # ----------------------------------------------------------

  annotate(
    "segment",
    x = 3,
    xend = 4,
    y = gse_bar_y,
    yend = gse_bar_y,
    linewidth = 0.45
  ) +

  annotate(
    "segment",
    x = 3,
    xend = 3,
    y = gse_bar_y - 0.025 * y_range,
    yend = gse_bar_y,
    linewidth = 0.45
  ) +

  annotate(
    "segment",
    x = 4,
    xend = 4,
    y = gse_bar_y - 0.025 * y_range,
    yend = gse_bar_y,
    linewidth = 0.45
  ) +

  annotate(
    "text",
    x = 3.5,
    y = gse_bar_y + 0.025 * y_range,
    label = "*",
    size = 4.5,
    fontface = "bold"
  ) +

  # ----------------------------------------------------------
  # THEME
  # ----------------------------------------------------------

  theme_classic(
    base_size = 11
  ) +

  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 13
    ),

    axis.title.y = element_text(
      face = "bold",
      size = 10.5
    ),

    axis.text.x = element_text(
      size = 8.5,
      colour = "black",
      lineheight = 0.9
    ),

    axis.text.y = element_text(
      colour = "black",
      size = 9
    ),

    legend.position = "top",

    legend.title = element_blank(),

    legend.text = element_text(
      size = 9
    ),

    axis.line = element_line(
      linewidth = 0.5
    ),

    axis.ticks = element_line(
      linewidth = 0.4
    ),

    plot.margin = margin(
      12, 12, 12, 12
    )
  ) +

  coord_cartesian(
    ylim = c(
      y_min,
      gse_bar_y + 0.08 * y_range
    ),
    clip = "off"
  )


# ------------------------------------------------------------
# 10. DISPLAY
# ------------------------------------------------------------

print(p_hr)


# ------------------------------------------------------------
# 11. SAVE
# ------------------------------------------------------------

ggsave(
  "Figure_1a_HR_TCGAvsGSE76427.pdf",
  p_hr,
  width = 6.5,
  height = 4.5,
  units = "in",
  device = cairo_pdf
)

ggsave(
  "Figure_1a_HR_TCGAvsGSE76427.png",
  p_hr,
  width = 6.5,
  height = 4.5,
  units = "in",
  dpi = 600
)
