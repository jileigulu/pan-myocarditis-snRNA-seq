# Differential expression and GO enrichment analysis for neutrophils

library(Seurat)
library(ggplot2)
library(dplyr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(readr)

set.seed(717)

# Directories
output.path2 <- "~/pan-heart/results/Neutrophils2025.5.13/files/deg/亚型层面"
dir.create(file.path(output.path2, 'deg'), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output.path2, 'gsea'), showWarnings = FALSE, recursive = TRUE)

output.go <- "~/pan-heart/results/Neutrophils2025.5.13/files/GO"
dir.create(output.go, showWarnings = FALSE, recursive = TRUE)

# 1. Load data and run DEG analysis----

seurat.obj <- readRDS("~/pan-heart/results/Neutrophils2025.5.13/files/merged.rds")
cell.types <- names(table(seurat.obj$cell_type2))
disease <- c("EAM", "ICI-MC", "VMC")

# Create group labels
seurat.obj$celltype.group <- paste(seurat.obj$cell_type2, seurat.obj$disease, sep = "_")
Idents(seurat.obj) <- "celltype.group"
DefaultAssay(seurat.obj) <- "RNA"

for (cell.type in cell.types) {
  for (d in disease) {
    ident.1 <- paste0(cell.type, "_", d)
    ident.2 <- paste0(cell.type, "_Normal")
    
    if (!ident.1 %in% levels(Idents(seurat.obj))) next
    
    # DEG analysis
    markers <- FindMarkers(seurat.obj, ident.1 = ident.1, ident.2 = ident.2)
    markers$celltype <- cell.type
    markers$disease <- d
    markers$gene <- rownames(markers)
    
    dir.create(file.path(output.path2, 'deg', cell.type), showWarnings = FALSE)
    write.table(markers, file = file.path(output.path2, 'deg', cell.type, paste0(d, ".csv")),
                quote = FALSE, sep = ",", row.names = FALSE)
    
    # GSEA input (all genes)
    gsea.input <- FindMarkers(seurat.obj, ident.1 = ident.1, ident.2 = ident.2,
                              min.pct = 0, logfc.threshold = 0)
    gsea.input$celltype <- cell.type
    gsea.input$disease <- d
    gsea.input$gene <- rownames(gsea.input)
    
    dir.create(file.path(output.path2, 'gsea', cell.type), showWarnings = FALSE)
    write.table(gsea.input, file = file.path(output.path2, 'gsea', cell.type, paste0(d, ".csv")),
                quote = FALSE, sep = ",", row.names = FALSE)
  }
}

# 2. GO enrichment analysis----

deg <- read_csv("~/pan-heart/results/Neutrophils2025.5.13/files/deg/deg1.csv")
disease <- unique(deg$disease)

for (d in disease) {
  # Upregulated
  dir.create(file.path(output.go, "Upregulated", d), showWarnings = FALSE, recursive = TRUE)
  deg_up <- deg[deg$disease == d & deg$regulate == "Upregulated", ]
  top_genes_up <- deg_up[order(deg_up$avg_log2FC, decreasing = TRUE), ]$gene
  
  if (length(top_genes_up) > 0) {
    bp_up <- enrichGO(top_genes_up, OrgDb = org.Mm.eg.db, keyType = 'SYMBOL',
                      ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2)
    term_up <- bp_up@result
    term_up$disease <- d
    write.table(term_up, file = file.path(output.go, "Upregulated", d, paste0(d, ".csv")),
                row.names = FALSE, quote = FALSE, sep = ",")
  }
  
  # Downregulated
  dir.create(file.path(output.go, "Downregulated", d), showWarnings = FALSE, recursive = TRUE)
  deg_down <- deg[deg$disease == d & deg$regulate == "Downregulated", ]
  top_genes_down <- deg_down[order(deg_down$avg_log2FC, decreasing = FALSE), ]$gene
  
  if (length(top_genes_down) > 0) {
    bp_down <- enrichGO(top_genes_down, OrgDb = org.Mm.eg.db, keyType = 'SYMBOL',
                        ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2)
    term_down <- bp_down@result
    term_down$disease <- d
    write.table(term_down, file = file.path(output.go, "Downregulated", d, paste0(d, ".csv")),
                row.names = FALSE, quote = FALSE, sep = ",")
  }
}