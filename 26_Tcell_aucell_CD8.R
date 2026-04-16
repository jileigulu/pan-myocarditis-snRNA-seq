
# AUCell analysis for CD8+ T cell subtypes

library(Seurat)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(pheatmap)
library(msigdbr)
library(AUCell)
library(clusterProfiler)
library(ggradar)
library(readr)

setwd('~/pan-heart/results/T_cells2025.5.20')

outdir <- './files/AUCell/CD8'
outdir2 <- './plots/AUCell/CD8'
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(outdir2, recursive = TRUE, showWarnings = FALSE)

# 1. Load data and prepare AUCell----

seurat.obj0 <- readRDS('./files/merged.rds')
seurat_obj <- subset(seurat.obj0, idents = c('CD8+TEFF', 'CD8+TEFF_Klrd1', 'CD8+TCM', 'CD8+Tisg'))

# Build rankings and calculate AUC
cells_rankings <- AUCell_buildRankings(seurat_obj@assays$RNA@counts)
fun <- read_excel("~/pan-heart/data/CD8function.xlsx")
geneSets <- split(fun$gene, fun$pathways)

cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank = nrow(cells_rankings) * 0.1)
saveRDS(cells_AUC, file.path(outdir, 'cells_AUC.rds'))

# Convert to data frame
AUCell_socre <- as.data.frame(t(cells_AUC@assays@data@listData[["AUC"]]))
colnames(AUCell_socre) <- colnames(AUCell_socre) %>%
  str_replace_all("_", " ") %>%
  str_to_sentence()
saveRDS(AUCell_socre, file.path(outdir, 'AUCell_socre.rds'))

# 2. Aggregate scores by cell subtype----

AUCell_socre <- readRDS(file.path(outdir, 'AUCell_socre.rds'))
seurat_obj$ID <- rownames(seurat_obj@meta.data)

meta <- AUCell_socre
meta$cell_type <- seurat_obj$cell_type2

object <- aggregate(meta[, -ncol(meta)], list(meta$cell_type), FUN = mean)
rownames(object) <- object$Group.1
object <- object[, -1]

saveRDS(object, file.path(outdir, "subcelltype_score.rds"))
write.csv(object, file.path(outdir, "subcelltype_score.csv"), row.names = TRUE)


# 3. Radar plot for selected pathways----


selected_pathways <- c("Chemokine/chemokine receptor", "Cytotoxicity", "Fatty acid metabolism")
filtered_Dataset <- object[, colnames(object) %in% selected_pathways]

# Normalize
for (i in 1:ncol(filtered_Dataset)) {
  filtered_Dataset[, i] <- (filtered_Dataset[, i] - min(filtered_Dataset[, i])) / 
    (max(filtered_Dataset[, i]) - min(filtered_Dataset[, i]))
}

filtered_Dataset$celltype <- rownames(filtered_Dataset)
radar_data <- filtered_Dataset[, c(4, 1, 2, 3)]

cell_type_cols <- c("#b3e19b", "#67a8cd", "#ED9B72", "#f36569")
names(cell_type_cols) <- levels(seurat_obj)

ggradar(radar_data, background.circle.transparency = 0,
        group.colours = cell_type_cols, plot.extent.x.sf = 1.5,
        base.size = 3, axis.label.size = 4, grid.label.size = 4,
        group.point.size = 2, group.line.width = 1,
        legend.position = 'right', legend.text.size = 10) +
  theme(plot.title = element_text(hjust = 0.2, size = 12))
ggsave(file.path(outdir2, 'radar_pathways.pdf'), last_plot(), width = 10, height = 5)

# 4. Bar plots for each pathway----

for (i in 1:length(selected_pathways)) {
  pathway <- selected_pathways[i]
  data <- filtered_Dataset[, c("celltype", pathway)]
  colnames(data) <- c("celltype", "val")
  
  p <- ggplot(data, aes(x = reorder(celltype, val), y = round(val, 4), fill = celltype)) +
    geom_bar(stat = 'identity', width = 0.7, fill = cell_type_cols) +
    labs(y = pathway, x = '') + scale_y_continuous(expand = c(0, 0)) +
    theme_classic() + coord_flip() +
    theme(axis.text.x = element_text(colour = 'black', size = 10),
          axis.text.y = element_text(colour = 'black', size = 10))
  
  ggsave(file.path(outdir2, paste0(gsub("/", " ", pathway), '_bar.pdf')), p, width = 5, height = 5)
}