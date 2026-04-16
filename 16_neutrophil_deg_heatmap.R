# Heatmap of DEGs across neutrophil subtypes and diseases

library(Seurat)
library(ggplot2)
library(pheatmap)
library(tidyverse)
library(readr)

setwd("~/pan-heart/results/Neutrophils2025.5.13")
outdir2 <- "./plots/deg"
dir.create(outdir2, showWarnings = FALSE)

# Load data
deg <- read_csv("files/deg/亚型层面/单独疾病组/deg/all_deg_genes.csv")
gsea <- read_csv("files/deg/亚型层面/单独疾病组/gsea/all_gsea_genes.csv")
gsea$name2 <- paste0(gsea$celltype, "-", gsea$gene)

# Selected genes for heatmap
selected_genes <- c('Nos2', 'Casp1', 'Casp4', 'Ptgs2', 'Relb', 'Stat1', 'Ifnar2', 'Ifngr1',
                    'Gbp2', 'Irgm1', 'Icam1', 'Cd14', 'Tap1', 'Tap2', 'Tapbp', 'Psme1',
                    'Psme2', 'Hif1a', 'Aldoa', 'Pgk1', 'Pfkl', 'Slc2a1')

# Prepare data
deg_subset <- gsea[gsea$gene %in% selected_genes, c("gene", "celltype", "disease", "avg_log2FC")]
deg_subset$label <- paste0(deg_subset$celltype, "-", deg_subset$disease)

# Reshape to matrix
heatmap_mat <- deg_subset %>%
  select(gene, label, avg_log2FC) %>%
  pivot_wider(names_from = label, values_from = avg_log2FC) %>%
  column_to_rownames("gene") %>%
  as.matrix()
heatmap_mat[is.na(heatmap_mat)] <- 0

# Significance matrix
sig_mat <- deg_subset %>%
  mutate(signi = ifelse(abs(avg_log2FC) > 0.25, "*", "")) %>%
  select(gene, label, signi) %>%
  pivot_wider(names_from = label, values_from = signi) %>%
  column_to_rownames("gene") %>%
  as.matrix()
sig_mat[is.na(sig_mat)] <- ""

# Annotation
ann_col <- data.frame(
  celltype = str_split_fixed(colnames(heatmap_mat), "-", 2)[,1],
  disease = str_split_fixed(colnames(heatmap_mat), "-", 2)[,2],
  row.names = colnames(heatmap_mat)
)
ann_color <- list(
  celltype = c("Neutro_Ccl4" = "#87CEEB", "Neutro_Txnip" = "#FF6F61",
               "Neutro_Ifit1" = "#A8E063", "Neutro_Ngp" = "#FFC0CB"),
  disease = c("EAM" = "#f57c6e", "ICI-MC" = "#B0DC66", "VMC" = "#EDCCEE")
)

# Plot heatmap
pdf(file.path(outdir2, "deg_heatmap.pdf"), height = 5, width = 5)
pheatmap(heatmap_mat, display_numbers = sig_mat, annotation_col = ann_col,
         annotation_colors = ann_color, cluster_rows = FALSE, cluster_cols = FALSE,
         color = colorRampPalette(c("#579CC7", "white", "#E41A1C"))(100),
         fontsize_number = 10, number_color = "black", cellwidth = 7, cellheight = 7,
         gaps_col = c(3, 6, 9), main = "")
dev.off()