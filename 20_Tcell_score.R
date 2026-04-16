# Module score analysis for T cells (heatmap + violin plots)

library(Seurat)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(ggpubr)

setwd("~/pan-heart/results/T_cells2025.5.20")
outdir2 <- "./plots/addmodulescore"
dir.create(outdir2, showWarnings = FALSE)
dir.create(file.path(outdir2, "heatmap"), showWarnings = FALSE)

# Load data
merged <- readRDS("./files/merged.rds")
fun <- read_csv("~/pan-heart/data/TC_function.csv")
dbs1 <- split(fun$gene, fun$pathways)

# Add module scores
seurat_T <- AddModuleScore(merged, features = dbs1, ctrl = 5, name = "FunctionScore")
for (i in 1:length(dbs1)) {
  colnames(seurat_T@meta.data)[colnames(seurat_T@meta.data) == paste0("FunctionScore", i)] <- names(dbs1)[i]
}

# Define pathway categories
pathway_categories <- list(
  `T Cell States & Lineage` = c('Naive', 'Treg signature', 'T cell receptor signaling pathway',
                                'Th17 cell differentiation', 'PD-L1 expression and PD-1 checkpoint pathway in cancer'),
  `T Cell Effector Function` = c('Cytotoxicity', 'Natural killer cell mediated cytotoxicity',
                                 'Chemokine/Chemokine receptor', 'Leukocyte transendothelial migration',
                                 'Chemokine signaling pathway', 'Cytokine-cytokine receptor interaction',
                                 'Antigen processing and presentation', 'Costimulatory molecules',
                                 'Fc gamma R-mediated phagocytosis', 'Activation:Effector function',
                                 'Cell adhesion molecules'),
  `Innate Immune Signaling` = c('IL-17 signaling pathway', 'TNF signaling pathway', 'NFKB Signaling',
                                'MAPK Signaling', 'NOD-like receptor signaling pathway',
                                'Cytosolic DNA-sensing pathway', 'IFN Response'),
  `Cell Metabolism & Survival` = c('Glycolysis', 'Oxidative phosphorylation', 'HIF-1 signaling pathway',
                                   'Lipid metabolism', 'Fatty acid metabolism', 'Anti-apoptosis', 'Pro-apoptosis'),
  `Autoimmune Related Diseases` = c('Rheumatoid arthritis', 'Inflammatory bowel disease', 'Allograft rejection',
                                    'Graft-versus-host disease', 'Type I diabetes mellitus')
)

all_pathways <- unlist(pathway_categories)
pathway_types <- rep(names(pathway_categories), times = sapply(pathway_categories, length))
names(pathway_types) <- all_pathways

# Prepare heatmap data
meta_agg <- seurat_T@meta.data %>%
  group_by(cell_type2) %>%
  summarise(across(all_of(all_pathways), mean)) %>%
  column_to_rownames("cell_type2") %>%
  t()

# Scale and plot
score_mat <- t(scale(t(meta_agg), center = TRUE, scale = TRUE))
score_mat <- score_mat[, levels(merged)]

row_annot <- data.frame(Signature.type = factor(pathway_types[rownames(score_mat)], 
                                                levels = names(pathway_categories)))
rownames(row_annot) <- rownames(score_mat)

pdf(file.path(outdir2, "heatmap", "heatmap.pdf"), width = 13, height = 13)
pheatmap(score_mat, annotation_row = row_annot, gaps_row = cumsum(table(row_annot$Signature.type)),
         gaps_col = 5, show_colnames = TRUE, show_rownames = TRUE, cellwidth = 14, cellheight = 12,
         color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100), fontsize = 11,
         cluster_rows = FALSE, cluster_cols = FALSE, main = "Pathway Signature Score")
dev.off()

# Violin plots
disease_colors <- c("Normal" = "#BCD5F1", "EAM" = "#f57c6e", "ICI-MC" = "#B0DC66", "VMC" = "#EDCCEE")

for (name in all_pathways) {
  p <- VlnPlot(seurat_T, group.by = "disease", features = name, split.by = "disease",
               pt.size = 0, cols = disease_colors) +
    geom_hline(yintercept = mean(seurat_T@meta.data[[name]]), linetype = "dashed", size = 0.5) +
    stat_compare_means(method = "t.test", ref.group = "Normal", label.y.npc = "top") +
    theme_bw() + theme(aspect.ratio = 0.5, axis.title.x = element_blank(),
                       axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank()) +
    labs(y = "Score", title = name) +
    geom_boxplot(width = 0.15, outlier.shape = NA, color = "black", fill = "white", alpha = 0.7)
  
  pdf(file.path(outdir2, paste0(gsub("/", "_", name), ".pdf")), width = 7, height = 5)
  print(p)
  dev.off()
}