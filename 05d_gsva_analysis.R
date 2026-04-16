# GSVA scoring and visualization for pathway enrichment

## SETUP ----
library(Seurat)
library(GSVA)
library(msigdbr)
library(tidyverse)
library(ComplexHeatmap)
library(pheatmap)
library(scales)
library(ggpubr)

setwd("~/pan-heart/results/deg2025.3.31")
dir.create("./files/GSVA", recursive = TRUE, showWarnings = FALSE)
dir.create("./plots/GSVA", recursive = TRUE, showWarnings = FALSE)
dir.create("./plots/GSVA/violin", recursive = TRUE, showWarnings = FALSE)

disease_colors <- c("Normal" = "#BCD5F1", "EAM" = "#f57c6e", "ICI-MC" = "#B0DC66", "VMC" = "#EDCCEE")

## LOAD DATA ----
seurat_obj <- readRDS("~/pan-heart/results/subset/immune_cells/2025.3.25/immune_cells.rds")

## DEFINE PATHWAYS ----
pathways_selected <- c(
  "positive regulation of inflammatory response", "leukocyte chemotaxis",
  "leukocyte cell-cell adhesion", "cytokine-mediated signaling pathway",
  "antigen processing and presentation", "positive regulation of defense response",
  "leukocyte mediated immunity", "lymphocyte mediated immunity",
  "regulation of apoptotic signaling pathway", "positive regulation of cell activation",
  "activation of immune response", "cell killing", "negative regulation of viral process",
  "activation of cysteine-type endopeptidase activity", "leukocyte proliferation",
  "positive regulation of interleukin-1 production", "t cell differentiation"
)

## PREPARE GENE SETS ----
c5_msigdb <- msigdbr(species = "Mus musculus", category = "C5")
c5_msigdb$gs_name <- c5_msigdb$gs_name %>%
  str_replace_all("GOBP_|HP_|GOCC_|GOMF_", "") %>%
  str_replace_all("_", " ") %>% str_to_lower()

pathways_selected_clean <- str_replace_all(pathways_selected, "-", " ")
gene_sets <- c5_msigdb[c5_msigdb$gs_name %in% pathways_selected_clean, ]
gene_list <- split(gene_sets$gene_symbol, gene_sets$gs_name)

## RUN GSVA ----
expr_matrix <- as.matrix(seurat_obj@assays$RNA@data)
es_matrix <- gsva(expr_matrix, gene_list, min.sz = 3, max.sz = Inf,
                  tau = 1, method = "gsva", abs.ranking = FALSE,
                  parallel.sz = 20, verbose = TRUE)
rownames(es_matrix) <- gsub(" ", "_", rownames(es_matrix))
saveRDS(es_matrix, "./files/GSVA/es_matrix.rds")

## ADD TO SEURAT OBJECT ----
seurat_obj <- AddMetaData(seurat_obj, metadata = t(es_matrix),
                          col.name = rownames(es_matrix))
saveRDS(seurat_obj, "./files/GSVA/seurat_obj.rds")

## HEATMAP BY CELL TYPE ----
meta <- seurat_obj@meta.data
pathway_cols <- colnames(meta)[15:ncol(meta)]

for (celltype in unique(meta$cell_type)) {
  meta_subset <- meta[meta$cell_type == celltype, c("disease", pathway_cols)]
  mean_scores <- aggregate(. ~ disease, data = meta_subset, FUN = mean)
  rownames(mean_scores) <- mean_scores$disease
  mean_scores <- mean_scores[, -1]
  
  # Scale matrix
  scaled_scores <- t(scale(t(mean_scores), center = TRUE, scale = TRUE))
  
  pdf(file.path("./plots/GSVA", paste0(celltype, ".pdf")), width = 10, height = 10)
  pheatmap(scaled_scores, show_colnames = TRUE, show_rownames = TRUE,
           cellwidth = 14, cellheight = 12, cluster_rows = FALSE, cluster_cols = FALSE,
           color = colorRampPalette(c("#f4dcef", "#b34557"))(100),
           main = celltype, border_color = "grey90")
  dev.off()
}

## VIOLIN PLOTS ----
for (pathway in pathway_cols) {
  p <- VlnPlot(seurat_obj, group.by = "disease", features = pathway,
               split.by = "disease", pt.size = 0, cols = disease_colors) +
    geom_hline(yintercept = mean(seurat_obj@meta.data[[pathway]]),
               linetype = "dashed", size = 0.5, color = "gray30") +
    stat_compare_means(aes(label = ..p.signif..), method = "t.test",
                       ref.group = "Normal", label.y.npc = "top") +
    theme_bw() + theme(aspect.ratio = 0.5, axis.title.x = element_blank(),
                       axis.text.x = element_text(angle = 45, hjust = 1),
                       panel.grid = element_blank()) +
    labs(y = "Score", title = pathway) +
    geom_boxplot(width = 0.15, outlier.shape = NA, color = "black", fill = "white", alpha = 0.7)
  
  ggsave(file.path("./plots/GSVA/violin", paste0(pathway, ".pdf")), p, width = 7, height = 5)
}
