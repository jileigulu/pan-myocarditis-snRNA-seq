# Heatmap of selected pathway genes across diseases and cell types

## SETUP ----
library(Seurat)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)

setwd("~/pan-heart/results/deg2025.3.31/deg")
dir.create("./plots/heatmap", showWarnings = FALSE)

## LOAD ALL DEGS ----
deg_eam <- read_csv("files/EAM_deg_TOTAL.csv") %>% mutate(disease = "EAM")
deg_ici <- read_csv("files/ICI-MC_deg_TOTAL.csv") %>% mutate(disease = "ICI-MC")
deg_vmc <- read_csv("files/VMC_deg_TOTAL.csv") %>% mutate(disease = "VMC")
deg_all <- bind_rows(deg_eam, deg_ici, deg_vmc)
write_csv(deg_all, "files/ALL_deg_TOTAL.csv")

## HELPER FUNCTION ----
plot_deg_heatmap <- function(deg_df, disease_name, gene_sets, colors_list, row_annotation) {
  deg_selected <- deg_df[deg_df$gene %in% names(row_annotation), ]
  deg_selected <- deg_selected[deg_selected$disease == disease_name, ]
  
  # Reshape to matrix
  plot_data <- deg_selected %>%
    dplyr::select(gene, celltype, avg_log2FC) %>%
    pivot_wider(names_from = "celltype", values_from = "avg_log2FC")
  
  matrix_data <- as.matrix(plot_data[, -1])
  rownames(matrix_data) <- plot_data$gene
  matrix_data[is.na(matrix_data)] <- 0
  
  # Order by gene set
  matrix_data <- matrix_data[names(row_annotation), , drop = FALSE]
  matrix_numeric <- apply(matrix_data, 2, as.numeric)
  rownames(matrix_numeric) <- rownames(matrix_data)
  
  # Color breaks
  max_val <- max(abs(matrix_numeric), na.rm = TRUE)
  breaks <- seq(-max_val, max_val, length.out = 100)
  if (!0 %in% breaks) breaks <- sort(c(breaks, 0))
  
  # Get row gaps
  counts <- table(row_annotation) %>% as.vector() %>% cumsum()
  
  pdf(file.path("./plots/heatmap", paste0(disease_name, ".pdf")))
  pheatmap(matrix_numeric, show_colnames = TRUE, show_rownames = TRUE,
           gaps_row = counts, annotation_row = row_annotation,
           annotation_colors = colors_list, cluster_rows = FALSE, cluster_cols = FALSE,
           color = colorRampPalette(c("#579CC7", "white", "#E41A1C"))(length(breaks) - 1),
           fontsize = 10, cellwidth = 12.5, cellheight = 12, breaks = breaks)
  dev.off()
}

## EAM GENE SETS ----
eam_genes <- list(
  Immune = c('Psmb9', 'Psme2', 'Icam1', 'Traf1'),
  Inflammation = c('S100a8', 'S100a9', 'Il1rn', 'Nfkb1', 'Arg2'),
  Oxidative_Stress = c('Prdx5', 'Hif1a', 'Txn1', 'Prelid1'),
  Cell_Death = c('Lgals3', 'G0s2', 'Calr', 'Tspo'),
  Metabolism = c('Atp5l', 'Pgk1', 'Aldoa', 'Cox5a', 'Cox7a2')
)

eam_annotation <- data.frame(Signature.type = factor(
  rep(names(eam_genes), times = sapply(eam_genes, length)),
  levels = names(eam_genes)))
rownames(eam_annotation) <- unlist(eam_genes)

eam_colors <- list(Signature.type = c(
  Immune = "#F99392", Inflammation = "#F9DE99",
  Oxidative_Stress = "#DBB2FF", Cell_Death = "#FCB985", Metabolism = "#EDCCEE"))

plot_deg_heatmap(deg_all, "EAM", eam_genes, eam_colors, eam_annotation)

## ICI-MC GENE SETS ----
ici_genes <- list(
  Immune = c('Rhoa', 'Tgfb1', 'Cd83', 'Irak2', 'Icam1', 'Cebpb'),
  Inflammation = c('Zfp36l2', 'Nfkb1', 'Nfkbiz', 'Tnfaip3', 'Ccl3'),
  Oxidative_Stress = c('Atf4', 'Txn1', 'Kdm6b', 'Tspo'),
  Cell_Death = c('Mcl1', 'Dnajb9', 'Birc3', 'Hnrnpk', 'Sqstm1', 'Bcl2a1d'),
  Metabolism = c('Gapdh', 'Cox5a', 'Ndufa13')
)

ici_annotation <- data.frame(Signature.type = factor(
  rep(names(ici_genes), times = sapply(ici_genes, length)),
  levels = names(ici_genes)))
rownames(ici_annotation) <- unlist(ici_genes)

plot_deg_heatmap(deg_all, "ICI-MC", ici_genes, eam_colors, ici_annotation)

## VMC GENE SETS ----
vmc_genes <- list(
  Immune = c('Ifi209', 'Ifi213', 'Rtp4', 'Ifit2', 'H2-D1'),
  Inflammation = c('Pycard', 'Hmgb1', 'Tnfaip3', 'Zfp36', 'Tgfb1'),
  Oxidative_Stress = c('Selenok', 'Gstp1', 'Gpx4'),
  Cell_Death = c('Zbp1', 'Ifi47', 'Irgm1', 'Xaf1', 'Foxo1', 'Birc3'),
  Metabolism = c('Ndufb4', 'Ndufb9', 'Ndufb7', 'Dnajc15'),
  Viral = c('Ifi27l2a', 'Bst2', 'Isg15')
)

vmc_annotation <- data.frame(Signature.type = factor(
  rep(names(vmc_genes), times = sapply(vmc_genes, length)),
  levels = names(vmc_genes)))
rownames(vmc_annotation) <- unlist(vmc_genes)

vmc_colors <- list(Signature.type = c(
  Immune = "#F99392", Viral = "#C5E0C9", Inflammation = "#F9DE99",
  Oxidative_Stress = "#DBB2FF", Cell_Death = "#FCB985", Metabolism = "#EDCCEE"))

plot_deg_heatmap(deg_all, "VMC", vmc_genes, vmc_colors, vmc_annotation)