# GSVA and AddModuleScore analysis for neutrophils

library(Seurat)
library(GSVA)
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(ggpubr)

setwd("~/pan-heart/results/Neutrophils2025.5.13")
outdir <- "./files/GSVA"
outdir2 <- "./plots/GSVA"
dir.create(outdir, showWarnings = FALSE)
dir.create(outdir2, showWarnings = FALSE)

merged <- readRDS("./files/merged.rds")

# Load gene sets
db <- read_csv("~/pan-heart/data/neutro_gsva.csv")
pathway_list <- lapply(db, function(x) x[!is.na(x)])

# Run GSVA
expr <- as.matrix(merged@assays$RNA@data)
es.matrix <- gsva(expr, pathway_list, min.sz = 3, max.sz = Inf, tau = 1,
                  method = "gsva", abs.ranking = FALSE, verbose = TRUE, parallel.sz = 80)
rownames(es.matrix) <- gsub(" ", "_", rownames(es.matrix))
saveRDS(es.matrix, file.path(outdir, "Neutrophilfunctionscore.es.matrix.rds"))

# Add to Seurat and aggregate by cell type
merged <- AddMetaData(merged, metadata = t(es.matrix), col.name = rownames(es.matrix))
meta <- merged@meta.data[, c("cell_type2", colnames(es.matrix))]
meta_agg <- aggregate(. ~ cell_type2, data = meta, FUN = mean)
rownames(meta_agg) <- meta_agg$cell_type2
meta_agg <- meta_agg[, -1]
meta_agg <- t(meta_agg)
colnames(meta_agg) <- gsub("_", " ", colnames(meta_agg))

# Define pathway categories
neutrophil_related <- c("Neutrophil Maturation", "Neutrophil Aging", "Neutrophil activation",
                        "Neutrophil extracellular trap formation", "Azurophil granulesb",
                        "Specific granules", "Gelatinase granules", "Secretory vesicles",
                        "Chemotaxis", "Phagocytosis")
cell_death <- c("Necroptosis", "Apoptosis", "Positive regulation of apoptotic process")
immune <- c("TNF signaling pathway", "IL 17 signaling pathway", "Type I interferon signaling pathway",
            "Interferon Gamma Response", "NF kappa B signaling pathway", "Toll like receptor signaling pathway",
            "NOD like receptor signaling pathway", "Tnfa Signaling Via Nfkb", "Cytokine cytokine receptor interaction",
            "Antigen Processing And Presentation", "Adaptive Immune Response", "Leukocyte Migration Involved In Inflammatory Response")
cytokine <- c("Cytokine Production", "Chemokine activity", "Interleukin 1 Production", "Cellular Response To Interleukin 1")
metabolism <- c("Glycolysis", "Oxidative phosphorylation", "Electron transport chain", "Tricarboxylic acid cycle", "HIF 1 signaling pathway")

all_pathways <- c(neutrophil_related, cell_death, immune, cytokine, metabolism)
pathway_types <- c(rep("Neutrophil Related", length(neutrophil_related)),
                   rep("Cell Death", length(cell_death)), rep("Immune", length(immune)),
                   rep("Cytokine", length(cytokine)), rep("Metabolism", length(metabolism)))
names(pathway_types) <- all_pathways

# Prepare heatmap data
score_mat <- meta_agg[all_pathways, ]
score_mat <- t(apply(score_mat, 1, scales::rescale, to = c(-1, 1)))
score_mat <- t(scale(t(score_mat), center = TRUE, scale = TRUE))

row_annot <- data.frame(Signature.type = factor(pathway_types, levels = unique(pathway_types)))
rownames(row_annot) <- all_pathways

# Plot heatmap
pdf(file.path(outdir2, "heatmap.pdf"), width = 10, height = 10)
pheatmap(score_mat, annotation_row = row_annot, gaps_row = cumsum(table(pathway_types)),
         show_colnames = TRUE, show_rownames = TRUE, cellwidth = 14, cellheight = 12,
         color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
         fontsize = 11, cluster_rows = FALSE, cluster_cols = FALSE,
         main = "Pathway Signature Score", border_color = "grey90")
dev.off()