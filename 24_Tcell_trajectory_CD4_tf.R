# Transcription factor visualization along trajectory for CD4+ T cells

library(monocle)
library(ggplot2)
library(RColorBrewer)
library(ClusterGVis)

setwd('~/pan-heart/results/T_cells2025.5.20')
input.dir <- "./files/trajectory/monocle2/CD4"
outdir <- "./files/trajectory/scenic/CD4"
outdir2 <- "./plots/trajectory/scenic/CD4"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(outdir2, recursive = TRUE, showWarnings = FALSE)

# Load data
cds <- readRDS(file.path(input.dir, "cds.rds"))
BEAM_res <- readRDS(file.path(input.dir, "BEAM_diff2000.rds"))

# Collect TFs from SCENIC results
tf_list <- c()
for (i in 1:5) {
  gene <- read.table(file.path("./files/trajectory/scenic/CD4", i, "tfs_targer.tsv"), header = TRUE, sep = ",")
  tf_list <- c(tf_list, unique(gene$tf))
}
tfs <- unique(tf_list)

# Generate heatmap data
df <- plot_genes_branched_heatmap2(cds[row.names(subset(BEAM_res)),],
                                   branch_point = 1, num_clusters = 5, cores = 10,
                                   use_gene_short_name = TRUE, show_rownames = FALSE)
saveRDS(df, file.path(outdir, "df.rds"))

# Plot heatmap with TFs
pdf(file.path(outdir2, "beauty_heatmap.pdf"), height = 10, width = 7)
visCluster(object = df, plot.type = "heatmap", markGenes = tfs, show_row_dend = FALSE,
           pseudotime_col = c("#0698ff", "grey", "#ff6b6f"),
           ctAnno.col = c("#f2b56f", "#b8aeeb", "#B0DC66", "#88d8db", "#71b7ed"),
           ht.col.list = list(col_range = seq(-3, 3, length.out = 100),
                              col_color = colorRampPalette(rev(brewer.pal(7, "Spectral")))(100)))
dev.off()

# Plot selected genes in pseudotime
disease_colors <- c("Normal" = "#BCD5F1", "EAM" = "#f57c6e", "ICI-MC" = "#B0DC66", "VMC" = "#EDCCEE")

# Selected genes to plot
s_genes <- c("Batf", "Fosl2", "Klf2", "Lef1", "Psmd12", "Nfkb2")
my_genes <- row.names(subset(fData(cds), gene_short_name %in% s_genes))
cds_subset <- cds[my_genes, ]

# Pseudotime plots
pdf(file.path(outdir2, "selected_genes_pseudotime.pdf"), width = 8, height = 4)
plot_genes_in_pseudotime(cds_subset, color_by = 'disease', ncol = 2) +
  scale_color_manual(values = disease_colors)
dev.off()

# Branched pseudotime plots
pdf(file.path(outdir2, "selected_genes_branched_pseudotime.pdf"), width = 8, height = 4)
plot_genes_branched_pseudotime(cds_subset, branch_point = 1, color_by = "disease", ncol = 2) +
  scale_color_manual(values = disease_colors)
dev.off()