
# 1. SETUP & CONFIGURATION----


library(ggplot2)
library(patchwork)
library(Nebulosa)
library(RColorBrewer)
library(scales)
library(ggalluvial)
library(ggforce)
library(ROGUE)
library(pheatmap)
library(data.table)
library(ggpubr)
library(dplyr)
library(tidyverse)

# Directories
outdir <- "~/pan-heart/results/landscape2025.3.25"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
setwd(outdir)
dir.create('./files', showWarnings = FALSE)
set.seed(2025)

# Color palettes
cell_type_cols <- c("#f57c6e", "#f2b56f", "#b8aeeb", "#B0DC66", "#88d8db", "#71b7ed",
                    "#f2a7da", "#fae69e", "#84c3b7", "#C5A094", "#EDCCEE", "#F5BCC8",
                    "#B3BD64", "#47B3C5", "#D9D2E2", "#BCD5F1", "#13989D", "#09A381",
                    "#EBF2C7", "#BA85B2", "#B7766C")

immune_colors <- c("NK cells" = "#f57c6e", "T cells" = "#f2b56f", "B cells" = "#b8aeeb",
                   "Dendritic cells" = "#B0DC66", "Macrophages&Monocytes" = "#88d8db",
                   "Neutrophils" = "#71b7ed")

disease_colors <- c("Normal" = "#BCD5F1", "EAM" = "#f57c6e", "ICI-MC" = "#B0DC66", "VMC" = "#EDCCEE")

# Load data
seurat.obj <- readRDS("~/pan-heart/results/subset/ALL/2025.3.25/merged.rds")
merged <- seurat.obj



# 2. UMAP VISUALIZATIONS----


## 2.1 All cell types ----
unique(Idents(merged))

all_cell_colors <- c("NK cells" = "#FFFF33", "T cells" = "#4DAF4A", "B cells" = "#E41A1C",
                     "Dendritic cells" = "#984EA3", "Macrophages&Monocytes" = "#377EB8",
                     "Neutrophils" = "#FF7F00", "Cardiomyocytes" = "#A65628",
                     "Endothelial cells" = "#F781BF", "Fibroblasts" = "#999999",
                     "Mural cells" = "#FF34B3", "Schwann cells" = "#BC8F8F")

DimPlot(merged, label = TRUE, pt.size = 0, cols = all_cell_colors) +
  NoLegend() + labs(x = "UMAP1", y = "UMAP2", title = "ALL") +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid"),
        axis.text = element_blank(), axis.ticks = element_blank())
ggsave('ALL.pdf', width = 8, height = 8)

## 2.2 Immune cells ----
seurat.immune <- readRDS("~/pan-heart/results/subset/immune_cells/2025.3.25/immune_cells.rds")
Idents(seurat.immune) <- seurat.immune$cell_type

DimPlot(seurat.immune, label = TRUE, pt.size = 0, cols = immune_colors) +
  NoLegend() + labs(x = "UMAP1", y = "UMAP2", title = "Immune Cells") +
  scale_x_continuous(expand = c(0.02, 0.02)) +
  scale_y_continuous(expand = c(0.02, 0.02)) +
  theme(panel.border = element_blank(),
        axis.text = element_blank(), axis.ticks = element_blank(),
        axis.line.x = element_line(color = "black", size = 0.8, arrow = arrow(length = unit(0.08, "inches"), type = "closed")),
        axis.line.y = element_line(color = "black", size = 0.8, arrow = arrow(length = unit(0.08, "inches"), type = "closed")))
ggsave('Immune_Cells.pdf', width = 5, height = 5)

## 2.3 Immune cell subtypes ----
Idents(seurat.immune) <- seurat.immune$cell_type2
Idents(seurat.immune) <- factor(Idents(seurat.immune), levels = c("B cells", "Macro_CCR2+MHCIIhi", "Macro_CCR2-MHCIIhi",
                                                                  "Macro_TLF+", "Monocytes", "Neutro_Ifit1", "Neutro_Txnip", "Neutro_Ngp", "Neutro_Ccl4",
                                                                  "CD4+TN", "CD4+TN_Tcf7", "CD4+TN_Rpl12", "CD4+Treg", "CD4+Th17",
                                                                  "CD8+TEFF", "CD8+TEFF_Klrd1", "CD8+TCM", "CD8+Tisg",
                                                                  "NK cells", "cDC1", "cDC2"))

DimPlot(seurat.immune, label = TRUE, pt.size = 0, cols = cell_type_cols) +
  labs(x = "UMAP1", y = "UMAP2", title = "Immune Cells", color = "Cell Types") +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid"),
        axis.text = element_blank(), axis.ticks = element_blank())
ggsave('Immune_Cells_sub.pdf', width = 12, height = 8)

## 2.4 Split by disease ----
Idents(seurat.immune) <- seurat.immune$cell_type
DimPlot(seurat.immune, label = TRUE, pt.size = 0, cols = immune_colors, split.by = "disease") +
  NoLegend() + labs(x = "UMAP1", y = "UMAP2", title = "Immune Cells") +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid"),
        axis.text = element_blank(), axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.5))
ggsave('Immune_Cells_disease.pdf', width = 16, height = 5)

## 2.5 Split by group ----
DimPlot(seurat.immune, label = TRUE, pt.size = 0, cols = immune_colors, split.by = "group") +
  NoLegend() + labs(x = "UMAP1", y = "UMAP2", title = "Immune Cells") +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid"),
        axis.text = element_blank(), axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.5))
ggsave('Immune_Cells_group.pdf', width = 8, height = 5)

## 2.6 By sample and group ----
DimPlot(seurat.immune, label = FALSE, pt.size = 0, cols = cell_type_cols, group.by = "orig.ident") +
  labs(x = "UMAP1", y = "UMAP2", title = "Immune Cells", color = "Samples") +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid"),
        axis.text = element_blank(), axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.5))
ggsave('Immune_Cells_orig.ident.pdf', width = 5.5, height = 4)

DimPlot(seurat.immune, label = FALSE, pt.size = 0, cols = c("#71b7ed", "#f57c6e"), group.by = "group") +
  labs(x = "UMAP1", y = "UMAP2", title = "Immune Cells", color = "Groups") +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid"),
        axis.text = element_blank(), axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.5))
ggsave('Immune_Cells_group2.pdf', width = 5, height = 4)



# 3. QUALITY CONTROL----


DefaultAssay(seurat.immune) <- 'RNA'
seurat.immune <- PercentageFeatureSet(seurat.immune, pattern = "mt\\.|mt-", col.name = "percent.mt")
seurat.immune <- subset(seurat.immune, subset = nFeature_RNA >= 200 & percent.mt <= 20)

VlnPlot(seurat.immune, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        group.by = "cell_type", cols = immune_colors, pt.size = 0)
ggsave("Vlnplot.pdf", width = 8, height = 4)



# 4. MARKER ANALYSIS----


## 4.1 Find and visualize top markers ----
markers <- FindAllMarkers(seurat.immune, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top.markers <- markers %>% group_by(cluster) %>% slice_max(n = 3, order_by = avg_log2FC)

source("~/tutorial/custom_plot_function.R")
levels(seurat.immune) <- c("B cells", "Dendritic cells", "Macrophages&Monocytes", "NK cells", "Neutrophils", "T cells")

p <- DotPlot(seurat.immune, dot.scale = 4, col.min = 0, features = unique(top.markers$gene)) +
  scale_color_distiller(palette = "RdBu", direction = -1) +
  coord_fixed(ratio = 2) + theme_cat() +
  theme(axis.title = element_blank(), legend.margin = margin(b = -8),
        legend.position = "top", axis.text.x = element_text(face = "italic")) +
  guides(x = guide_axis(angle = 90), color = guide_colorbar(title = "Expression", frame.colour = "black", ticks.colour = "black")) +
  scale_y_discrete(limits = rev)
ggsave("topmarkers.pdf", p, width = 6, height = 3)

## 4.2 GO enrichment analysis ----
library(org.Mm.eg.db)
dir.create(file.path("files", "GO"), showWarnings = FALSE, recursive = TRUE)

for (d in unique(markers$cluster)) {
  deg1 <- markers[markers$cluster %in% d, ]
  top.genes <- deg1[order(deg1$avg_log2FC, decreasing = TRUE), ]$gene[1:100]
  bp <- enrichGO(top.genes, OrgDb = org.Mm.eg.db, keyType = 'SYMBOL',
                 ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2)
  term <- bp@result
  term$celltype <- d
  write.table(term, file = file.path("files", "GO", paste0(d, ".csv")), row.names = FALSE, quote = FALSE, sep = ",")
}

# Load and combine GO results
go_pathways <- list(
  "B cells" = c("B cell activation", "B cell receptor signaling pathway"),
  "Macrophages&Monocytes" = c("antigen processing and presentation of peptide antigen via MHC class II", "monocyte chemotaxis"),
  "T cells" = c("T cell receptor signaling pathway", "T cell differentiation"),
  "Dendritic cells" = c("antigen processing and presentation of exogenous peptide antigen via MHC class II", "positive regulation of T cell activation"),
  "Neutrophils" = c("positive regulation of inflammatory response", "neutrophil chemotaxis"),
  "NK cells" = c("natural killer cell mediated cytotoxicity", "positive regulation of natural killer cell mediated cytotoxicity")
)

df_list <- list()
for (celltype in names(go_pathways)) {
  go_data <- read_csv(paste0("files/GO/", celltype, ".csv"))
  df_list[[celltype]] <- go_data %>% filter(Description %in% go_pathways[[celltype]])
}

df <- do.call(rbind, df_list) %>%
  mutate(log_pvalue = -log10(pvalue),
         celltype = factor(celltype, levels = c("B cells", "Dendritic cells", "Macrophages&Monocytes", "NK cells", "Neutrophils", "T cells"))) %>%
  arrange(celltype, desc(log_pvalue)) %>%
  mutate(Description = factor(Description, levels = unique(Description)))

p_go <- ggplot(df, aes(x = log_pvalue, y = Description, fill = celltype)) +
  geom_bar(stat = 'identity') +
  geom_text(aes(label = Description, x = 0.1), hjust = 0, size = 3.5, color = "black", fontface = "bold") +
  scale_fill_manual(values = immune_colors) +
  facet_grid(celltype ~ ., scales = "free_y", space = "free_y") +
  theme_bw() + theme(panel.grid = element_blank(), axis.text = element_text(color = "black", size = 10),
                     axis.title.x = element_text(color = "black", size = 8), axis.title.y = element_blank(),
                     axis.text.y = element_blank(), axis.ticks.y = element_blank(), strip.text.y = element_text(angle = 0)) +
  scale_y_discrete(limits = rev)
ggsave("markers_GO_bar.pdf", p_go, height = 5, width = 12)

## 4.3 Feature density plots - ALL ----
markers_all <- FindAllMarkers(merged, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
features_all <- c("Myl2", "Cd79a", "Fabp4", "Cd68", "Gsn", "Cd3d", "H2-Aa", "S100a8", "Gzma", "Csrp1", "Fabp7")
DefaultAssay(merged) <- 'RNA'
pdf('feature_density_ALL.pdf', width = 15, height = 8)
plot_density(merged, reduction = 'umap', features_all, joint = TRUE)
dev.off()

## 4.4 Feature density plots - Immune ----
features_immune <- c("Cd79a", "Cd68", "Cd3d", "H2-Aa", "S100a8", "Gzma")
pdf('feature_density_immune.pdf', width = 12, height = 8)
plot_density(seurat.immune, reduction = 'umap', features_immune, joint = TRUE)
dev.off()



# 5. CELL PROPORTION ANALYSIS----


## 5.1 Bar plot - cell numbers ----
meta <- seurat.immune@meta.data %>% dplyr::count(cell_type) %>%
  arrange(n) %>% mutate(cell_type = factor(cell_type, levels = cell_type))

ggplot(meta, aes(x = cell_type, y = n, fill = cell_type)) +
  geom_bar(stat = "identity", position = "stack", width = 0.9, color = "white", show.legend = FALSE) +
  labs(x = "Cell Types", y = "Cell Numbers") + coord_flip() +
  scale_fill_manual(values = immune_colors) +
  theme(axis.line = element_line(size = 0.5, colour = "black"),
        panel.grid = element_blank(), panel.background = element_blank()) +
  scale_y_continuous(expand = c(0, 0))
ggsave('cell_number.pdf', height = 5, width = 7, scale = 0.8)

## 5.2 Sankey bar plots ----
df_prop <- as.data.frame(prop.table(table(Idents(seurat.immune), seurat.immune$disease), margin = 2) * 100)
ggplot(df_prop, aes(x = Var2, y = Freq, fill = Var1, stratum = Var1, alluvium = Var1)) +
  geom_col(width = 0.5, color = 'black', just = 0.5) + geom_flow(width = 0.5, alpha = 0.4, knot.pos = 0) +
  theme_bw() + labs(x = 'Disease', y = 'Fraction of Clusters', fill = "Cell Type") +
  scale_fill_manual(values = immune_colors) +
  theme(axis.line = element_line(color = 'black'), panel.grid = element_blank())
ggsave('cell_percentage.pdf', height = 3, width = 4)

df_disease <- as.data.frame(prop.table(table(seurat.immune$disease, Idents(seurat.immune)), margin = 2) * 100)
ggplot(df_disease, aes(x = Var2, y = Freq, fill = Var1, stratum = Var1, alluvium = Var1)) +
  geom_col(width = 0.5, color = 'black', just = 0.5) + geom_flow(width = 0.5, alpha = 0.4, knot.pos = 0) +
  theme_bw() + labs(x = 'Cell Type', y = 'Fraction of Diseases', fill = "Disease") +
  scale_fill_manual(values = disease_colors) +
  theme(axis.line = element_line(color = 'black'), panel.grid = element_blank())
ggsave('cell_percentage_disease.pdf', height = 3, width = 4)

## 5.3 Detailed Sankey plot - subtypes ----
proportion <- data.frame(prop.table(table(seurat.immune$cell_type2)) * 100)
data_sankey <- data.frame(cluster_num = proportion$Var1, 
                          cell_type = c('B cells', rep('T cells', 9), 'Dendritic cells', 'Dendritic cells',
                                        rep('Macrophages&Monocytes', 4), rep('Neutrophils', 4), 'NK cells'),
                          proportion = proportion$Freq)

data_sankey <- gather_set_data(data_sankey[, c(1,2,3)], 1:2)
data_sankey[1:21, ]$x <- 'cluster_num'
data_sankey[22:42, ]$x <- 'cell_type'

target_order <- c("B cells", "T cells", "Dendritic cells", "Macrophages&Monocytes", "Neutrophils", "NK cells")
all_y <- unique(data_sankey$y)
sub_clusters <- setdiff(all_y, target_order)
setlevels <- c(target_order, sub_clusters)
setlevels <- c(setlevels[setlevels != "NK cells"], "NK cells")
data_sankey$y <- factor(data_sankey$y, levels = setlevels)

ggplot(data_sankey, aes(x = x, id = id, split = y, fill = y, value = proportion)) +
  geom_parallel_sets(aes(fill = cell_type), alpha = 1, axis.width = 0.1) +
  geom_parallel_sets_axes(axis.width = 0.2, color = 'grey') +
  geom_parallel_sets_labels(colour = 'black', angle = 0) +
  theme_bw() + theme(panel.grid = element_blank(), panel.border = element_blank(),
                     axis.line = element_blank(), axis.ticks = element_blank(),
                     axis.text = element_blank(), legend.position = "none") +
  scale_fill_manual(values = setNames(cell_type_cols, levels(data_sankey$y)))
ggsave("cell_percentage2.pdf", height = 4.1, width = 5)

## 5.4 Boxplot by group ----
source("~/tutorial/custom_plot_function.R")
prop_data <- as.data.frame(prop.table(table(Idents(seurat.immune), seurat.immune$orig.ident), margin = 2) * 100)
colnames(prop_data) <- c('types', 'orig.ident', 'prop')

sample_info <- seurat.immune@meta.data[, c("orig.ident", "disease")]
sample_info <- sample_info[!duplicated(sample_info$orig.ident), ]
prop_data <- merge(prop_data, sample_info, by = "orig.ident")
colnames(prop_data)[4] <- "group"

p_box <- ggboxplot(prop_data, x = 'types', y = 'prop', fill = 'group', color = 'black',
                   width = 0.7, size = 0.2, legend = 'top', show.points = FALSE, outlier.shape = NA) +
  scale_fill_manual(values = disease_colors) +
  labs(x = 'Cell Type', y = 'Proportion', fill = 'Disease') + theme_cat() +
  theme(axis.text.x = element_text(size = 8, angle = 30, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 12)) +
  geom_vline(xintercept = seq(1.5, 20, 1), linewidth = 0.2, color = "#a09c9a", linetype = "dashed")
ggsave("cell_percentage_box.pdf", p_box, width = 4, height = 3)


