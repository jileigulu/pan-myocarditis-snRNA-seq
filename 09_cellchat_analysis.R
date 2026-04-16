# Cell-cell communication analysis using CellChat
# Analyzes macrophage-lymphocyte interactions across disease conditions

## SETUP ----
library(CellChat)
library(patchwork)
library(Seurat)
library(future)

setwd("~/pan-heart/results/Macrophage2025.4.16")

# Directories
cellchat_dir <- "./files/cellchat"
plot_dir <- "./plots/cellchat"
dir.create(cellchat_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

-# 1. PREPARE DATA (Run once to create merged object)----
-
prepare_cellchat_data <- function() {
  # Load macrophage data
  seurat_obj <- readRDS("./files/merged.rds")
  
  # Load immune cells and subset lymphocytes
  merged <- readRDS("~/pan-heart/results/subset/immune_cells/2025.3.25/immune_cells.rds")
  lymphocytes <- subset(merged, cell_type3 == "Lymphocyte")
  
  # Merge macrophage and lymphocyte data
  seurat_merged <- merge(seurat_obj, lymphocytes)
  
  # Reclassify T cell subtypes
  meta <- seurat_merged@meta.data
  meta$cell_type2 <- ifelse(grepl("CD8", meta$cell_type2), "CD8+T",
                            ifelse(grepl("CD4", meta$cell_type2), "CD4+T",
                                   meta$cell_type2))
  seurat_merged@meta.data <- meta
  Idents(seurat_merged) <- seurat_merged$cell_type2
  
  saveRDS(seurat_merged, file.path(cellchat_dir, "merged.rds"))

  return(seurat_merged)
}

# Uncomment to run data preparation (only needed once)
# seurat_obj1 <- prepare_cellchat_data()

# 2. CELLCHAT ANALYSIS FUNCTION----

run_cellchat_analysis <- function(seurat_obj, disease_group, output_dir, n_workers = 10) {
  message("\n=== Running CellChat for: ", disease_group, " ===")
  
  # Subset by disease
  cells_subset <- seurat_obj[, seurat_obj[["disease"]] == disease_group]
  
  if (ncol(cells_subset) == 0) {
    warning("No cells found for disease group: ", disease_group)
    return(NULL)
  }
  
  # Prepare expression data
  expr_matrix <- cells_subset@assays$RNA@data
  meta_data <- as.data.frame(Idents(cells_subset))
  colnames(meta_data) <- "labels"
  
  # Create CellChat object
  cellchat <- createCellChat(object = expr_matrix, meta = meta_data, group.by = "labels")
  cellchat <- addMeta(cellchat, meta = meta_data)
  cellchat <- setIdent(cellchat, ident.use = "labels")
  
  group_size <- as.numeric(table(cellchat@idents))
  message("Cell types: ", paste(levels(cellchat@idents), collapse = ", "))
  message("Cell counts: ", paste(group_size, collapse = ", "))
  
  # Set up CellChat database
  CellChatDB <- CellChatDB.mouse
  cellchat@DB <- CellChatDB
  
  # Process data
  cellchat <- subsetData(cellchat)
  plan("multisession", workers = n_workers)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- projectData(cellchat, PPI.mouse)
  
  # Compute communication networks
  cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  
  # Create output directory
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Plot 1: Number of interactions (circle plot)
  pdf(file.path(output_dir, "netVisual_circle_by_count.pdf"), width = 8, height = 6)
  netVisual_circle(cellchat@net$count, vertex.weight = group_size,
                   weight.scale = TRUE, label.edge = FALSE,
                   title.name = "Number of interactions")
  dev.off()
  
  # Plot 2: Interaction weights/strength
  pdf(file.path(output_dir, "netVisual_circle_by_weight.pdf"), width = 8, height = 6)
  netVisual_circle(cellchat@net$weight, vertex.weight = group_size,
                   weight.scale = TRUE, label.edge = FALSE,
                   title.name = "Interaction weights/strength")
  dev.off()
  
  # Plot 3: Individual cell type interaction networks
  pdf(file.path(output_dir, "netVisual_circle_by_cell_type.pdf"), width = 12, height = 10)
  mat <- cellchat@net$weight
  n_celltypes <- nrow(mat)
  par(mfrow = c(ceiling(n_celltypes / 2), 2), xpd = TRUE)
  
  for (i in 1:n_celltypes) {
    mat2 <- matrix(0, nrow = n_celltypes, ncol = n_celltypes, dimnames = dimnames(mat))
    mat2[i, ] <- mat[i, ]
    netVisual_circle(mat2, vertex.weight = group_size, weight.scale = TRUE,
                     edge.weight.max = max(mat), title.name = rownames(mat)[i])
  }
  dev.off()
  
  # Save CellChat object
  saveRDS(cellchat, file = file.path(output_dir, "cellchat.rds"))
  
  message("CellChat analysis completed for: ", disease_group)
  message("Results saved to: ", output_dir)
  
  return(cellchat)
}

# 3. RUN CELLCHAT FOR ALL DISEASE GROUPS----

# Load merged data
seurat_obj <- readRDS(file.path(cellchat_dir, "merged.rds"))

# Define disease groups to analyze
disease_groups <- c("Normal", "EAM", "ICI-MC", "VMC")

# Run analysis for each disease group
cellchat_results <- list()

for (group in disease_groups) {
  output_dir <- file.path(cellchat_dir, group)
  result <- run_cellchat_analysis(seurat_obj, group, output_dir, n_workers = 10)
  cellchat_results[[group]] <- result
}

