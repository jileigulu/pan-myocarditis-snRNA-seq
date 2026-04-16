# DEG petal plot and pathway gene heatmap for T cells

library(Seurat)
library(tidyverse)
library(pheatmap)
library(reshape2)

setwd("~/pan-heart/results/T_cells2025.5.20")
outdir2 <- "./plots/deg"
dir.create(outdir2, showWarnings = FALSE)
dir.create(file.path(outdir2, "selected_pathway"), showWarnings = FALSE)

# 1. DEG petal plot----

deg <- read_csv("files/deg/deg1.csv")
df_count <- table(deg$disease, deg$regulate) %>% as.data.frame()
colnames(df_count) <- c("Sample", "Group", "Freq")

ggplot(df_count, aes(x = Sample, y = Freq, fill = Group)) +
  geom_col(width = 0.9, alpha = 0.8) + coord_polar() +
  labs(x = "Number of DEG") + theme_bw() +
  theme(aspect.ratio = 1, axis.title.y = element_blank(),
        axis.text = element_text(size = 8), legend.position = "top") +
  scale_fill_manual(values = c("Upregulated" = "#E41A1C", "Downregulated" = "#579CC7"))
ggsave(file.path(outdir2, "差异基因花瓣图.pdf"), width = 5, height = 5)

# 2. Pathway gene heatmap (selected genes)----

# Selected pathways and genes
pathway_genes <- list(
  `T cell receptor signaling pathway` = c("Ctla4", "Cd8b1", "Nfkb1", "Itk"),
  Exhaustion = c("Lag3", "Tigit", "Havcr2", "Cd160", "Tox"),
  `Th17 cell differentiation` = c("Rora", "Stat3", "Jak2", "Tgfb1"),
  `PD-L1 expression and PD-1 checkpoint pathway in cancer` = c("Pdcd1", "Cd274", "Ifng", "Batf", "Hif1a"),
  `Costimulatory molecules` = c("Tnfrsf9", "Tnfrsf4", "Tnfrsf18", "Icos"),
  `Anti apoptosis` = c("Mcl1", "Gadd45b", "Bcl2", "Bcl2l1")
)

all_genes <- unlist(pathway_genes)
pathway_names <- rep(names(pathway_genes), times = sapply(pathway_genes, length))
names(pathway_names) <- all_genes

# Prepare heatmap data
deg_selected <- deg[deg$gene %in% all_genes, c("gene", "disease", "avg_log2FC")]
heatmap_mat <- deg_selected %>%
  pivot_wider(names_from = "disease", values_from = "avg_log2FC") %>%
  column_to_rownames("gene") %>%
  as.matrix()
heatmap_mat[is.na(heatmap_mat)] <- 0
heatmap_mat <- heatmap_mat[all_genes[all_genes %in% rownames(heatmap_mat)], ]

# Row annotation
row_annot <- data.frame(Pathways = factor(pathway_names[rownames(heatmap_mat)], 
                                          levels = names(pathway_genes)))
rownames(row_annot) <- rownames(heatmap_mat)

ann_colors <- list(Pathways = c(
  "T cell receptor signaling pathway" = "#f57c6e", "Exhaustion" = "#b8aeeb",
  "Th17 cell differentiation" = "#B0DC66", "PD-L1 expression and PD-1 checkpoint pathway in cancer" = "#F5BCC8",
  "Costimulatory molecules" = "#88d8db", "Anti apoptosis" = "#FBDBC4"
))

# Plot heatmap
min_val <- min(heatmap_mat, na.rm = TRUE)
max_val <- max(heatmap_mat, na.rm = TRUE)
breaks <- seq(min_val, max_val, length.out = 100)
breaks <- sort(unique(c(breaks, 0)))
zero_idx <- which(breaks == 0)
color_left <- colorRampPalette(c("#2166AC", "#F0F0F0"))(zero_idx - 1)
color_right <- colorRampPalette(c("#F0F0F0", "#C8191D"))(length(breaks) - zero_idx)
color_vals <- c(color_left, color_right)

pdf(file.path(outdir2, "gene_heatmap.pdf"), height = 8, width = 8)
pheatmap(heatmap_mat, annotation_row = row_annot, annotation_colors = ann_colors,
         gaps_row = cumsum(table(row_annot$Pathways)), cluster_rows = FALSE, cluster_cols = FALSE,
         color = color_vals, fontsize = 10, cellwidth = 12.5, cellheight = 12, breaks = breaks)
dev.off()