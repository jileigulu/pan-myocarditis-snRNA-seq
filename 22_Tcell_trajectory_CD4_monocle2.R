# Monocle2 trajectory analysis for CD4+ T cells

library(monocle)
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

setwd('~/pan-heart/results/T_cells2025.5.20')
outdir <- './files/trajectory/monocle2/CD4'
outdir2 <- './plots/trajectory/monocle2/CD4'
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(outdir2, recursive = TRUE, showWarnings = FALSE)


# 1. Prepare data and run Monocle2----

seurat.obj0 <- readRDS('./files/merged.rds')
seurat.obj <- subset(seurat.obj0, idents = c('CD4+TN', 'CD4+TN_Tcf7', 'CD4+TN_Rpl12', 'CD4+Treg', 'CD4+Th17'))
DefaultAssay(seurat.obj) <- 'RNA'

# Find markers
markers <- FindAllMarkers(seurat.obj, only.pos = TRUE)
saveRDS(markers, file.path(outdir, "markers.rds"))

# Convert to Monocle object
data <- as(as.matrix(seurat.obj@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = seurat.obj@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
cds <- newCellDataSet(data, phenoData = pd, featureData = fd,
                      lowerDetectionLimit = 0.5, expressionFamily = negbinomial.size())

# QC
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- detectGenes(cds, min_expr = 0.1)

# Select ordering genes
top_markers <- markers %>%
  filter(p_val <= 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  group_by(cluster) %>%
  top_n(15, wt = avg_log2FC) %>%
  pull(gene)

saveRDS(top_markers, file.path(outdir, "top_markers.rds"))

# Build trajectory
cds <- setOrderingFilter(cds, top_markers)
cds <- reduceDimension(cds, max_components = 2, reduction_method = 'DDRTree', verbose = FALSE)
cds <- orderCells(cds)

saveRDS(cds, file.path(outdir, "cds.rds"))

# 2. Trajectory plots----

cds <- readRDS(file.path(outdir, "cds.rds"))

# Color palettes
disease_colors <- c("Normal" = "#BCD5F1", "EAM" = "#f57c6e", "ICI-MC" = "#B0DC66", "VMC" = "#EDCCEE")
cell_type_cols <- c("#f8a8a8", "#c1b5d9", "#97d2ce", "#aa9e93", "#ffc17f")
names(cell_type_cols) <- levels(seurat.obj)

# Individual plots
p1 <- plot_cell_trajectory(cds, color_by = "Pseudotime", cell_size = 0.1, size = 1, show_backbone = TRUE)
p2 <- plot_cell_trajectory(cds, color_by = "State", cell_size = 0.1, size = 1, show_backbone = TRUE) +
  scale_color_manual(values = c("#f57c6e", "#f2b56f", "#b8aeeb"))
p3 <- plot_cell_trajectory(cds, color_by = "disease", cell_size = 0.1, size = 1, show_backbone = TRUE) +
  scale_color_manual(values = disease_colors)
p5 <- plot_cell_trajectory(cds, color_by = "cell_type2", cell_size = 0.1, size = 1, show_backbone = TRUE, show_branch_points = FALSE) +
  scale_color_manual(values = cell_type_cols)

# Combined plots
p_combined <- (p1 + p2) / (p3 + p5)
ggsave(file.path(outdir2, 'plot_cell_trajectory1.pdf'), p_combined, width = 10, height = 10)

# Facet plots
p4 <- plot_cell_trajectory(cds) + facet_wrap(~disease, nrow = 1) +
  scale_color_manual(values = c("#f57c6e", "#f2b56f", "#b8aeeb"))
p6 <- plot_cell_trajectory(cds) + facet_wrap(~cell_type2, nrow = 1) +
  scale_color_manual(values = c("#f57c6e", "#f2b56f", "#b8aeeb"))
ggsave(file.path(outdir2, 'plot_cell_trajectory2.pdf'), p4 / p6, width = 8, height = 6)

# Individual saves
ggsave(file.path(outdir2, 'Pseudotime.pdf'), p1, width = 3.5, height = 4)
ggsave(file.path(outdir2, 'sub_cell_type.pdf'), p5, width = 5, height = 5)

# 3. Density plots----

df <- pData(cds)

p_density1 <- ggplot(df, aes(Pseudotime, fill = cell_type2)) +
  geom_density(bw = 0.5, alpha = 0.8) + scale_fill_manual(values = cell_type_cols) +
  labs(fill = "Subtype") + theme_classic2()

p_density2 <- ggplot(df, aes(Pseudotime, fill = disease)) +
  geom_density(bw = 0.5, alpha = 0.8) + scale_fill_manual(values = disease_colors) +
  labs(fill = "Disease") + theme_classic2()

p_density3 <- ggplot(df, aes(Pseudotime, fill = group)) +
  geom_density(bw = 0.5, alpha = 0.8) + scale_fill_manual(values = c("#7AA6DCFF", "#E41A1C")) +
  theme_classic2()

ggsave(file.path(outdir2, 'density_sub_cell_type1.pdf'), p_density1 + p_density2 + p_density3 + plot_layout(ncol = 1), width = 6, height = 13)

# 4. State proportion bar plots----

State_colors <- c("#f57c6e", "#f2b56f", "#b8aeeb")

p_bar1 <- ggplot(df, aes(x = cell_type2, fill = State)) +
  geom_bar(position = "fill", width = 0.6, color = "white") +
  labs(x = "", y = "Percentage of cells (%)") + scale_fill_manual(values = State_colors) +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(outdir2, "geom_bar_State_classic1.pdf"), p_bar1, height = 4, width = 3)

p_bar2 <- ggplot(df, aes(x = disease, fill = State)) +
  geom_bar(position = "fill", width = 0.6, color = "white") +
  labs(x = "", y = "Percentage of cells (%)") + scale_fill_manual(values = State_colors) +
  theme_minimal()
ggsave(file.path(outdir2, "geom_bar_State_classic2.pdf"), p_bar2, height = 4, width = 3)

# 5. Branch expression analysis----

# Differential gene test
expressed_genes <- row.names(subset(fData(cds), num_cells_expressed >= 10))
diff <- differentialGeneTest(cds[expressed_genes,], fullModelFormulaStr = "~cell_type2", cores = 60)
saveRDS(diff, file.path(outdir, 'diff.rds'))

# BEAM analysis for branch point
BEAM_res <- BEAM(cds[unique(top_markers),], branch_point = 1)
saveRDS(BEAM_res, file.path(outdir, "BEAM_deg.rds"))

# Heatmap
pdf(file.path(outdir2, "branched_diff_heatmap1.pdf"), width = 5, height = 5)
plot_genes_branched_heatmap(cds[row.names(subset(BEAM_res, qval < 0.0001)),],
                            branch_point = 1, num_clusters = 5, cores = 10,
                            use_gene_short_name = TRUE, show_rownames = FALSE,
                            return_heatmap = TRUE)
dev.off()