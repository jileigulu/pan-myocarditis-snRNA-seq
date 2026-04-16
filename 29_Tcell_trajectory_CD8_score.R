# Pathway score visualization along Monocle3 trajectory for CD8+ T cells

library(monocle3)
library(tidyverse)
library(ggplot2)

setwd('~/pan-heart/results/T_cells2025.5.20')
outdir <- './files/trajectory/monocle3/CD8'
outdir2 <- './plots/trajectory/monocle3/CD8/AUCell'
dir.create(outdir2, recursive = TRUE, showWarnings = FALSE)

# 1. Add AUCell scores to CDS object----

AUCell_socre <- readRDS("~/pan-heart/results/T_cells2025.5.20/files/AUCell/CD8/AUCell_socre.rds")
cds <- readRDS(file.path(outdir, "cds.rds"))

# Add scores to phenoData
cdsdata <- pData(cds)
b <- cbind(cdsdata, AUCell_socre)
pData(cds) <- b

# Clean column names
colnames(colData(cds)) <- gsub("/", " ", colnames(colData(cds)))
saveRDS(cds, file.path(outdir, "cds_score.rds"))

# 2. Plot pathway scores along trajectory----

cds <- readRDS(file.path(outdir, "cds_score.rds"))

pathways <- c("Chemokine chemokine receptor", "Cytotoxicity", "Fatty acid metabolism")

for (i in pathways) {
  p <- plot_cells(cds, color_cells_by = i, label_cell_groups = FALSE,
                  label_leaves = FALSE, label_branch_points = FALSE, graph_label_size = 1.5)
  ggsave(file.path(outdir2, paste0(gsub(" ", "_", i), ".pdf")), p, width = 5.5, height = 3)
}

# Individual plots with adjusted dimensions
p1 <- plot_cells(cds, color_cells_by = "Chemokine chemokine receptor",
                 label_cell_groups = FALSE, label_leaves = FALSE,
                 label_branch_points = FALSE, graph_label_size = 1.5)
ggsave(file.path(outdir2, "Chemokine_chemokine_receptor.pdf"), p1, width = 6, height = 3)

p2 <- plot_cells(cds, color_cells_by = "Cytotoxicity",
                 label_cell_groups = FALSE, label_leaves = FALSE,
                 label_branch_points = FALSE, graph_label_size = 1.5)
ggsave(file.path(outdir2, "Cytotoxicity.pdf"), p2, width = 4.5, height = 3)

p3 <- plot_cells(cds, color_cells_by = "Fatty acid metabolism",
                 label_cell_groups = FALSE, label_leaves = FALSE,
                 label_branch_points = FALSE, graph_label_size = 1.5)
ggsave(file.path(outdir2, "Fatty_acid_metabolism.pdf"), p3, width = 5, height = 3)

# 3. Pathway score vs pseudotime (smooth curves)----

ptime <- pseudotime(cds)
meta <- as.data.frame(colData(cds))
meta$Pseudotime <- ptime

# Clean column names
colnames(meta) <- gsub(".", " ", colnames(meta), fixed = TRUE)

# Reshape data for plotting
df <- meta %>%
  dplyr::select(Pseudotime, `Chemokine chemokine receptor`, Cytotoxicity, `Fatty acid metabolism`) %>%
  tidyr::pivot_longer(
    cols = c(`Chemokine chemokine receptor`, Cytotoxicity, `Fatty acid metabolism`),
    names_to = "Pathway", values_to = "Score"
  )

p <- ggplot(df, aes(x = Pseudotime, y = Score, color = Pathway)) +
  geom_smooth(se = FALSE, span = 0.3, size = 1.2) +
  theme_classic() + labs(y = "Pathway score")

ggsave(file.path(outdir2, "pathway_vs_pseudotime.pdf"), p, width = 6, height = 4)

