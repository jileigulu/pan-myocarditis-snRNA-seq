# Monocle3 trajectory analysis for CD8+ T cells

library(Seurat)
library(monocle3)
library(tidyverse)
library(patchwork)

setwd('~/pan-heart/results/T_cells2025.5.20')
outdir <- './files/trajectory/monocle3/CD8'
outdir2 <- './plots/trajectory/monocle3/CD8'
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(outdir2, recursive = TRUE, showWarnings = FALSE)

# 1. Prepare data and create CDS object----

seurat.obj0 <- readRDS('./files/merged.rds')
seurat_obj <- subset(seurat.obj0, idents = c('CD8+TEFF', 'CD8+TEFF_Klrd1', 'CD8+TCM', 'CD8+Tisg'))

set.seed(2025)

# Extract data
data <- GetAssayData(seurat_obj, assay = 'RNA', slot = 'counts')
cell_metadata <- seurat_obj@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)

# Create CDS object
cds <- new_cell_data_set(data, cell_metadata = cell_metadata, gene_metadata = gene_annotation)

# Preprocess and reduce dimensions
cds <- preprocess_cds(cds, num_dim = 50)
cds <- reduce_dimension(cds, preprocess_method = "PCA")

# Use Seurat UMAP embeddings
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(seurat_obj, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed), ]
cds@int_colData$reducedDims$UMAP <- int.embed

# Cluster and learn graph
cds <- cluster_cells(cds, reduction_method = "UMAP")
cds <- learn_graph(cds)

# Function to get earliest principal node
get_earliest_principal_node <- function(cds, cell_phenotype = "broad_cell_type", root_types) {
  root_pr_nodes <- lapply(root_types, function(root_type) {
    cell_ids <- which(pData(cds)[, cell_phenotype] == root_type)
    closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
    closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
    root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(
      which.max(table(closest_vertex[cell_ids, ]))))]
    return(root_pr_nodes)
  })
  return(unlist(root_pr_nodes))
}

# Order cells
cds <- order_cells(cds, root_pr_nodes = get_earliest_principal_node(cds,
                                                                    cell_phenotype = "cell_type2", root_types = c('CD8+TEFF')))

saveRDS(cds, file.path(outdir, "cds.rds"))

# 2. Plot trajectory----

cds <- readRDS(file.path(outdir, "cds.rds"))

# Rename column for plotting
colnames(colData(cds))[colnames(colData(cds)) == "cell_type2"] <- "Subtype"

cell_type_cols <- c("#b3e19b", "#67a8cd", "#ED9B72", "#f36569")
names(cell_type_cols) <- levels(seurat_obj)

p <- plot_cells(cds, reduction_method = "UMAP", trajectory_graph_color = "grey28",
                color_cells_by = "Subtype", label_groups_by_cluster = FALSE,
                label_leaves = FALSE, label_cell_groups = FALSE, label_branch_points = FALSE) +
  scale_color_manual(values = cell_type_cols, limits = levels(seurat_obj)) +
  ggtitle('Monocle3') + theme(aspect.ratio = 1)

ggsave(file.path(outdir2, "monocle3_cell_type.pdf"), p, height = 4, width = 4)

# 3. Find trajectory-associated genes----

# Note: Run 'trace('calculateLW', edit = T)' and change Matrix::rBind to rbind if needed
Track_genes <- graph_test(cds, neighbor_graph = 'principal_graph', cores = 6)
saveRDS(Track_genes, file.path(outdir, "Track_genes.rds"))

Track_genes <- Track_genes[, c(5, 2, 3, 4, 1, 6)] %>% filter(q_value < 1e-3)
write.csv(Track_genes, file.path(outdir, "Track_genes.csv"), row.names = FALSE)