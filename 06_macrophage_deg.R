# Differential expression analysis for macrophage subtypes
# Compares each macrophage subtype across disease conditions

## SETUP ----
library(Seurat)
library(ggplot2)
library(dplyr)

set.seed(717)

# Output directories
output_base <- "~/pan-heart/results/Macrophage2025.4.16/files/deg/亚型层面"
dir.create(file.path(output_base, "三种心肌炎组/deg"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_base, "三种心肌炎组/gsea"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_base, "整体疾病组/deg"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_base, "整体疾病组/gsea"), showWarnings = FALSE, recursive = TRUE)

## LOAD DATA ----
seurat.obj <- readRDS("~/pan-heart/results/Macrophage2025.4.16/files/merged.rds")
DefaultAssay(seurat.obj) <- "RNA"

# Get macrophage subtypes
cell.types <- names(table(seurat.obj$cell_type2))

# 1. THREE DISEASE COMPARISONS (EAM, ICI-MC, VMC vs Normal)----

seurat.obj$celltype.group <- paste(seurat.obj$cell_type2, seurat.obj$disease, sep = "_")
Idents(seurat.obj) <- "celltype.group"

diseases <- c("EAM", "ICI-MC", "VMC")

for (cell.type in cell.types) {
  for (d in diseases) {
    ident.1 <- paste0(cell.type, "_", d)
    ident.2 <- paste0(cell.type, "_Normal")
    
    # Skip if ident.1 doesn't exist
    if (!ident.1 %in% levels(Idents(seurat.obj))) {
      message("Skipping ", ident.1, " - not found")
      next
    }
    
    # Run DEG analysis
    markers <- FindMarkers(seurat.obj,
                           ident.1 = ident.1,
                           ident.2 = ident.2)
    markers$celltype <- cell.type
    markers$disease <- d
    markers$gene <- rownames(markers)
    
    # Save DEG results
    deg_dir <- file.path(output_base, "三种心肌炎组/deg", cell.type)
    dir.create(deg_dir, showWarnings = FALSE)
    write.table(markers,
                file = file.path(deg_dir, paste0(d, ".csv")),
                quote = FALSE, sep = ",", row.names = FALSE)
    
    # Run for GSEA input (all genes, no filtering)
    gsea.input <- FindMarkers(seurat.obj,
                              ident.1 = ident.1,
                              ident.2 = ident.2,
                              min.pct = 0,
                              logfc.threshold = 0)
    gsea.input$celltype <- cell.type
    gsea.input$disease <- d
    gsea.input$gene <- rownames(gsea.input)
    
    gsea_dir <- file.path(output_base, "三种心肌炎组/gsea", cell.type)
    dir.create(gsea_dir, showWarnings = FALSE)
    write.table(gsea.input,
                file = file.path(gsea_dir, paste0(d, ".csv")),
                quote = FALSE, sep = ",", row.names = FALSE)
  }
}

message("Three disease comparisons completed!")

# 2. COMBINED MYOCARDITIS COMPARISON (All diseases vs Normal)----

seurat.obj$celltype.group <- paste(seurat.obj$cell_type2, seurat.obj$group, sep = "_")
Idents(seurat.obj) <- "celltype.group"

for (cell.type in cell.types) {
  ident.1 <- paste0(cell.type, "_Myocarditis")
  ident.2 <- paste0(cell.type, "_Normal")
  
  if (!ident.1 %in% levels(Idents(seurat.obj))) {
    message("Skipping ", ident.1, " - not found")
    next
  }
  
  # Run DEG analysis
  markers <- FindMarkers(seurat.obj,
                         ident.1 = ident.1,
                         ident.2 = ident.2)
  markers$celltype <- cell.type
  markers$gene <- rownames(markers)
  
  # Save DEG results
  write.table(markers,
              file = file.path(output_base, "整体疾病组/deg", paste0(cell.type, ".csv")),
              quote = FALSE, sep = ",", row.names = FALSE)
  
  # Run for GSEA input
  gsea.input <- FindMarkers(seurat.obj,
                            ident.1 = ident.1,
                            ident.2 = ident.2,
                            min.pct = 0,
                            logfc.threshold = 0)
  gsea.input$celltype <- cell.type
  gsea.input$gene <- rownames(gsea.input)
  
  write.table(gsea.input,
              file = file.path(output_base, "整体疾病组/gsea", paste0(cell.type, "_gsea.csv")),
              quote = FALSE, sep = ",", row.names = FALSE)
}

message("Combined myocarditis comparison completed!")

# 3. CLASSIFY UP/DOWN REGULATED GENES----

read_csv_dir <- function(dir_path) {
  file_list <- list.files(path = dir_path, full.names = TRUE, pattern = "\\.csv$")
  df_list <- lapply(file_list, function(x) read.csv(x, stringsAsFactors = FALSE))
  do.call(rbind, df_list)
}

classify_regulation <- function(deg_df, logFC_cutoff = 0.25, p_cutoff = 0.05) {
  deg_df$regulate <- as.factor(ifelse(
    deg_df$p_val < p_cutoff & abs(deg_df$avg_log2FC) > logFC_cutoff,
    ifelse(deg_df$avg_log2FC > logFC_cutoff, 'Upregulated', 'Downregulated'),
    'NOT'
  ))
  return(deg_df)
}

## 3.1 Three disease comparisons ----
deg_3disease <- read_csv_dir(file.path(output_base, "三种心肌炎组/deg"))
deg_3disease <- classify_regulation(deg_3disease)

# Add combined label for summary
deg_3disease$label <- paste0(deg_3disease$celltype, " ", deg_3disease$disease)

# Save with regulation info
write.csv(deg_3disease,
          file.path(output_base, "三种心肌炎组/deg/deg_with_regulation.csv"),
          row.names = FALSE)

## 3.2 Combined disease comparison ----
deg_combined <- read_csv_dir(file.path(output_base, "整体疾病组/deg"))
deg_combined <- classify_regulation(deg_combined)

write.csv(deg_combined,
          file.path(output_base, "整体疾病组/deg/deg_with_regulation.csv"),
          row.names = FALSE)

