# Load libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(ggalluvial)

# Load demo data
seurat.demo <- readRDS("demo_immune_cells.rds")

# 1. UMAP plot (immune cell types)
DimPlot(seurat.demo, label = TRUE, pt.size = 0.5) +
  labs(title = "Immune Cell Types") +
  theme_minimal()

# 2. QC violin plot
VlnPlot(seurat.demo, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, pt.size = 0)

# 3. Feature plot of canonical marker genes
marker_genes <- c("Cd79a", "Cd3d", "Cd68", "S100a8")  # B, T, Macrophage, Neutrophil
FeaturePlot(seurat.demo, features = marker_genes, ncol = 2, pt.size = 0.5)

# 4. (Optional) Dot plot of selected markers
DotPlot(seurat.demo, features = marker_genes) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
