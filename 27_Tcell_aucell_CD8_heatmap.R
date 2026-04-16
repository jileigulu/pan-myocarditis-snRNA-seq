# Gene expression heatmap for CD8+ T cell pathways

library(Seurat)
library(tidyverse)
library(pheatmap)

setwd('~/pan-heart/results/T_cells2025.5.20')
outdir <- './files/AUCell/CD8'
outdir2 <- './plots/AUCell/CD8'

# Load data
seurat.obj0 <- readRDS('./files/merged.rds')
seurat_obj <- subset(seurat.obj0, idents = c('CD8+TEFF', 'CD8+TEFF_Klrd1', 'CD8+TCM', 'CD8+Tisg'))

# Load pathway genes
pathway_gene <- read_csv("~/pan-heart/results/T_cells/score/AUCell/files/AUCell/pathway_gene.csv")

# Function to get average expression and select top genes
get_top_genes <- function(seurat_obj, celltype, pathway_name, pathway_gene, n_top = 1) {
  # Subset cells
  cells <- seurat_obj[, seurat_obj$cell_type2 == celltype]
  
  # Average expression by disease
  av_expr <- AverageExpression(cells, assays = "RNA", group.by = "disease")[["RNA"]]
  
  # Get pathway genes
  features <- pathway_gene[pathway_gene$pathway == pathway_name, ]$gene
  features <- unique(features[!is.na(features)])
  
  # Filter to available genes
  available <- features[features %in% rownames(av_expr)]
  expr <- av_expr[available, , drop = FALSE]
  expr <- na.omit(expr)
  
  if (nrow(expr) == 0) return(NULL)
  
  # Calculate mean across diseases and select top genes
  expr$mean <- rowMeans(expr)
  expr <- expr[order(expr$mean, decreasing = TRUE), ]
  
  if (n_top == 1) {
    return(expr[1, 1:4, drop = FALSE])
  } else {
    return(expr[1:min(n_top, nrow(expr)), 1:4, drop = FALSE])
  }
}

# Cell types
celltypes <- c('CD8+TEFF', 'CD8+TEFF_Klrd1', 'CD8+TCM', 'CD8+Tisg')
pathways <- c("Cytotoxicity", "Chemokine/Chemokine receptor", "Fatty acid metabolism")

# Collect top genes for each pathway and cell type
all_expr <- data.frame()

for (pathway in pathways) {
  for (ct in celltypes) {
    expr <- get_top_genes(seurat_obj, ct, pathway, pathway_gene, n_top = 1)
    if (!is.null(expr)) {
      expr$gene <- rownames(expr)
      expr$pathway <- pathway
      expr$celltype <- ct
      all_expr <- rbind(all_expr, expr %>% pivot_longer(cols = c("Normal", "EAM", "ICI-MC", "VMC"),
                                                        names_to = "disease", values_to = "expression"))
    }
  }
}

# Reshape for heatmap
heatmap_data <- all_expr %>%
  select(gene, celltype, disease, expression) %>%
  mutate(label = paste0(celltype, "_", disease)) %>%
  select(gene, label, expression) %>%
  pivot_wider(names_from = label, values_from = expression) %>%
  column_to_rownames("gene") %>%
  as.matrix()

heatmap_data[is.na(heatmap_data)] <- 0

# Create annotations
annotation_col <- data.frame(
  Disease = factor(rep(c("Normal", "EAM", "ICI-MC", "VMC"), 4),
                   levels = c("Normal", "EAM", "ICI-MC", "VMC")),
  Celltype = factor(rep(celltypes, each = 4), levels = celltypes)
)
rownames(annotation_col) <- colnames(heatmap_data)

annotation_row <- data.frame(
  Pathway = factor(c(rep("Cytotoxicity", 1), rep("Chemokine/Chemokine receptor", 1), rep("Fatty acid metabolism", 1)),
                   levels = pathways)
)
rownames(annotation_row) <- rownames(heatmap_data)

# Colors
disease_colors <- c("Normal" = "#BCD5F1", "EAM" = "#f57c6e", "ICI-MC" = "#B0DC66", "VMC" = "#EDCCEE")
cell_type_cols <- c("#b3e19b", "#67a8cd", "#ED9B72", "#f36569")
names(cell_type_cols) <- celltypes

ann_colors <- list(
  Celltype = cell_type_cols,
  Disease = disease_colors,
  Pathway = c('Cytotoxicity' = '#9E6CDF', 'Chemokine/Chemokine receptor' = '#B1D9D4', 'Fatty acid metabolism' = '#F49792')
)

# Plot heatmap
pdf(file.path(outdir2, "gene_heatmap.pdf"), width = 10, height = 5)
pheatmap(heatmap_data, annotation_colors = ann_colors,
         color = colorRampPalette(c("#3980b8", "white", "#ef3b3c"))(256),
         scale = 'row', annotation_col = annotation_col, annotation_row = annotation_row,
         cellwidth = 8, cellheight = 8, border_color = "white", fontsize = 8,
         cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = TRUE, show_colnames = FALSE,
         gaps_row = cumsum(table(annotation_row$Pathway)),
         gaps_col = cumsum(table(annotation_col$Celltype)))
dev.off()