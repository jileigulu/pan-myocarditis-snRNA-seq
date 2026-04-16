# Differential expression analysis for immune cells across disease conditions
# Compares EAM, ICI-MC, and VMC against Normal for each cell type


## SETUP ----
library(Seurat)
library(ggplot2)
library(dplyr)

set.seed(717)
setwd("~/pan-heart/results/deg2025.3.31")
dir.create("./deg", recursive = TRUE, showWarnings = FALSE)
dir.create("./gsea", recursive = TRUE, showWarnings = FALSE)

## LOAD DATA ----
seurat.obj <- readRDS("~/pan-heart/results/subset/immune_cells/2025.3.25/immune_cells.rds")
DefaultAssay(seurat.obj) <- "RNA"

## PREPARE GROUPS ----
cell.types <- names(table(seurat.obj$cell_type))
seurat.obj$celltype.group <- paste(seurat.obj$cell_type, seurat.obj$disease, sep = "_")
Idents(seurat.obj) <- "celltype.group"

# Define comparisons
comparisons <- list(
  EAM = "EAM",
  `ICI-MC` = "ICI-MC",
  VMC = "VMC"
)

## RUN DIFFERENTIAL EXPRESSION ----
run_deg_analysis <- function(seurat_obj, cell_types, disease_name, output_dir) {
  dir.create(file.path("./deg", disease_name), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path("./gsea", disease_name), showWarnings = FALSE, recursive = TRUE)
  
  for (cell.type in cell_types) {
    ident.1 <- paste0(cell.type, "_", disease_name)
    ident.2 <- paste0(cell.type, "_Normal")
    
    # Skip if ident.1 doesn't exist in the data
    if (!ident.1 %in% levels(Idents(seurat_obj))) {
      message("Skipping ", ident.1, " - not found in data")
      next
    }
    
    # Run FindMarkers for DEGs (commented out, kept for reference)
    # markers <- FindMarkers(seurat_obj,
    #                        ident.1 = ident.1,
    #                        ident.2 = ident.2)
    # markers$celltype <- cell.type
    # markers$gene <- rownames(markers)
    # write.table(markers,
    #             file = paste0("./deg/", disease_name, "/", gsub(" ", "_", cell.type), ".txt"),
    #             quote = FALSE, sep = ",", row.names = FALSE)
    
    # Run FindMarkers for GSEA input (all genes, no filtering)
    gsea.input <- FindMarkers(
      seurat_obj,
      ident.1 = ident.1,
      ident.2 = ident.2,
      min.pct = 0,
      logfc.threshold = 0
    )
    gsea.input$celltype <- cell.type
    gsea.input$gene <- rownames(gsea.input)
    
    write.table(
      gsea.input,
      file = paste0("./gsea/", disease_name, "/", gsub(" ", "_", cell.type), "_gsea.csv"),
      quote = FALSE, sep = ",", row.names = FALSE
    )
  }
}

# Run for all disease comparisons
for (disease in names(comparisons)) {
  message("Processing ", disease, "...")
  run_deg_analysis(seurat.obj, cell.types, comparisons[[disease]])
}

