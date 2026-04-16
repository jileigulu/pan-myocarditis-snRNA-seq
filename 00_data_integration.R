# Data integration and preprocessing for multiple myocarditis datasets
# Merges 3 independent scRNA-seq datasets, performs batch correction,
# and annotates cell types

library(dplyr)
library(Seurat)
library(harmony)
library(clustree)
library(lisi)
library(stringr)

# 1. LOAD RAW DATASETS----


merged1 <- readRDS('~/pan-heart/data/2020_myocarditis_HuaX/results/merged.rds')
merged2 <- readRDS('~/pan-heart/data/2022_myocarditis_AxelrodML/results/merged.rds')
merged3 <- readRDS('~/pan-heart/data/2022_myocarditis_MantriM/results/merged.rds')

# 2. CLEAN METADATA AND ADD DISEASE LABELS----

# Dataset 1: EAM
merged1@meta.data <- subset(merged1@meta.data, select = -c(time, phase, RNA_snn_res.1, 
                                                           seurat_clusters, sample_type, group_type))
merged1[['disease']] <- 'EAM'

# Dataset 2: ICI-MC
merged2@meta.data <- subset(merged2@meta.data, select = c(orig.ident, nCount_RNA, nFeature_RNA, group, percent.mt))
merged2[['disease']] <- 'ICI-MC'

# Dataset 3: Viral myocarditis (VMC)
merged3@meta.data <- subset(merged3@meta.data, select = c(orig.ident, nCount_RNA, nFeature_RNA, group, percent.mt))
merged3[['disease']] <- 'VMC'

# 3. RENAME GENES FOR DATASET 3 (Remove GRCm38- prefix)----

RenameGenesSeurat_v2 <- function(obj, newnames, gene.use = NULL, de.assay = "RNA") {
  print("Run this before integration. It only changes obj@assays$*@counts, @data and @scale.data")
  lassays <- Assays(obj)
  assay.use <- obj@reductions$pca@assay.used
  DefaultAssay(obj) <- de.assay
  
  if (is.null(gene.use)) {
    all_genenames <- rownames(obj)
  } else {
    all_genenames <- gene.use
    obj <- subset(obj, features = gene.use)
  }
  
  order_name <- function(v1, v2, ref) {
    v2 <- make.names(v2, unique = TRUE)
    df1 <- data.frame(v1, v2)
    rownames(df1) <- df1$v1
    df1 <- df1[ref, ]
    return(df1)
  }
  
  df1 <- order_name(v1 = all_genenames, v2 = newnames, ref = rownames(obj))
  all_genenames <- df1$v1
  newnames <- df1$v2
  
  if ('SCT' %in% lassays) {
    if ('SCTModel.list' %in% slotNames(obj@assays$SCT)) {
      obj@assays$SCT@SCTModel.list$model1@feature.attributes <- obj@assays$SCT@SCTModel.list$model1@feature.attributes[all_genenames, ]
      rownames(obj@assays$SCT@SCTModel.list$model1@feature.attributes) <- newnames
    }
  }
  
  change_assay <- function(a1 = de.assay, obj, newnames = NULL, all_genenames = NULL) {
    RNA <- obj@assays[a1][[1]]
    if (nrow(RNA) == length(newnames)) {
      if (length(RNA@counts)) RNA@counts@Dimnames[[1]] <- newnames
      if (length(RNA@data)) RNA@data@Dimnames[[1]] <- newnames
      if (length(RNA@var.features)) {
        df1 <- order_name(v1 = all_genenames, v2 = newnames, ref = RNA@var.features)
        RNA@var.features <- df1$v2
      }
      if (length(RNA@scale.data)) {
        df1 <- order_name(v1 = all_genenames, v2 = newnames, ref = rownames(RNA@scale.data))
        rownames(RNA@scale.data) <- df1$v2
      }
    } else {
      print("Unequal gene sets: nrow(RNA) != nrow(newnames)")
    }
    obj@assays[a1][[1]] <- RNA
    return(obj)
  }
  
  for (a in lassays) {
    DefaultAssay(obj) <- a
    df1 <- order_name(v1 = all_genenames, v2 = newnames, ref = rownames(obj))
    obj <- change_assay(obj = obj, a1 = a, newnames = df1$v2, all_genenames = df1$v1)
  }
  
  hvg <- VariableFeatures(obj, assay = assay.use)
  if (length(obj@reductions$pca)) {
    df1 <- order_name(v1 = all_genenames, v2 = newnames, ref = hvg)
    rownames(obj@reductions$pca@feature.loadings) <- df1$v2
  }
  
  return(obj)
}

# Apply gene renaming
gene_names <- rownames(merged3@assays[["RNA"]]@data)
gene_names <- gsub("GRCm38-", "", gene_names)
merged3 <- RenameGenesSeurat_v2(merged3, gene_names, de.assay = "RNA")
merged3 <- PercentageFeatureSet(merged3, pattern = "mt\\.|mt-", col.name = "percent.mt")

# 4. INITIAL MERGE AND PREPROCESSING----

merge_obj <- merge(merged1, c(merged2, merged3))
merge_obj <- merge_obj %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA(verbose = FALSE)

saveRDS(merge_obj, '~/pan-heart/results/myocarditis_merged.rds')

# 5. BATCH CORRECTION - HARMONY----

merged <- readRDS('~/pan-heart/results/myocarditis_merged.rds')
merged <- RunHarmony(merged, "orig.ident")

merged <- FindNeighbors(merged, reduction = "harmony", dims = 1:30)
merged <- FindClusters(merged, resolution = 0.1)
merged <- RunTSNE(merged, dims = 1:30, reduction = "harmony")
merged <- RunUMAP(merged, dims = 1:30, reduction = "harmony")

saveRDS(merged, '~/pan-heart/results/merge/harmony/myocarditis_merged_harmony.rds')

# 6. BATCH CORRECTION - SCT----

merged <- readRDS('~/pan-heart/results/myocarditis_merged.rds')
seurat.list <- SplitObject(merged, split.by = "orig.ident")
seurat.list <- lapply(seurat.list, function(x) {
  SCTransform(x, method = "glmGamPoi", variable.features.n = 3000,
              vars.to.regress = c("percent.mt"), verbose = TRUE)
})

features <- SelectIntegrationFeatures(object.list = seurat.list, nfeatures = 3000)
seurat.list <- PrepSCTIntegration(object.list = seurat.list, anchor.features = features)

immune.anchors <- FindIntegrationAnchors(object.list = seurat.list,
                                         normalization.method = "SCT",
                                         anchor.features = features)

merged <- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT")
merged <- RunPCA(merged)
merged <- FindNeighbors(merged, dims = 1:30)
merged <- FindClusters(merged, resolution = 0.1)
merged <- RunTSNE(merged, dims = 1:30)
merged <- RunUMAP(merged, dims = 1:30)

saveRDS(merged, '~/pan-heart/results/merge/SCT/myocarditis_merged_SCT.rds')

# 7. BATCH CORRECTION - CCA----

merged <- readRDS('~/pan-heart/results/myocarditis_merged.rds')
seurat.list <- SplitObject(merged, split.by = "orig.ident")
features <- SelectIntegrationFeatures(object.list = seurat.list)
merged.anchors <- FindIntegrationAnchors(object.list = seurat.list, anchor.features = features)
merged <- IntegrateData(anchorset = merged.anchors, dims = 1:30)

DefaultAssay(merged) <- "integrated"
merged <- ScaleData(merged, features = rownames(merged))
merged <- RunPCA(merged, npcs = 30, verbose = FALSE)
merged <- FindNeighbors(merged, reduction = "pca", dims = 1:30)
merged <- FindClusters(merged, resolution = 0.1)
merged <- RunTSNE(merged, reduction = "pca", dims = 1:30)
merged <- RunUMAP(merged, reduction = "pca", dims = 1:30)

saveRDS(merged, '~/pan-heart/results/merge/CCA/myocarditis_merged_CCA.rds')

# 8. BATCH CORRECTION - CCA + SCT----

merged <- readRDS('~/pan-heart/results/merge/CCA/myocarditis_merged_CCA.rds')
seurat.list <- SplitObject(merged, split.by = "orig.ident")
seurat.list <- lapply(seurat.list, function(x) {
  SCTransform(x, method = "glmGamPoi", variable.features.n = 3000,
              vars.to.regress = c("percent.mt"), verbose = TRUE)
})

features <- SelectIntegrationFeatures(object.list = seurat.list, nfeatures = 3000)
seurat.list <- PrepSCTIntegration(object.list = seurat.list, anchor.features = features)
immune.anchors <- FindIntegrationAnchors(object.list = seurat.list,
                                         normalization.method = "SCT",
                                         anchor.features = features)
merged <- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT")

merged <- RunPCA(merged)
merged <- FindNeighbors(merged, dims = 1:30)
merged <- FindClusters(merged, resolution = 0.1)
merged <- RunTSNE(merged, dims = 1:30)
merged <- RunUMAP(merged, dims = 1:30)

saveRDS(merged, '~/pan-heart/results/merge/CCA&SCT/myocarditis_merged_CCA_SCT.rds')

# 9. LISI SCORE EVALUATION (Batch correction quality)----

# SCT
lisi_score <- compute_lisi(X = Embeddings(merged, reduction = "pca"),
                           meta_data = merged@meta.data,
                           label_colnames = "orig.ident")
write_csv(lisi_score, '~/pan-heart/results/merge/SCT/lisi_score.csv')

# Harmony
merged_harmony <- readRDS("~/pan-heart/results/merge/harmony/myocarditis_merged_harmony.rds")
lisi_score <- compute_lisi(X = Embeddings(merged_harmony, reduction = "harmony"),
                           meta_data = merged_harmony@meta.data,
                           label_colnames = "orig.ident")
write_csv(lisi_score, '~/pan-heart/results/merge/harmony/lisi_score.csv')

# 10. CELL TYPE ANNOTATION----

merged <- readRDS("~/pan-heart/results/merge/CCA&SCT/myocarditis_merged_CCA_SCT.rds")
DefaultAssay(merged) <- "RNA"

# Annotation markers
# B cells: Cd79a, Cd79b, Cd19, Ms4a1
# T cells: Cd3d, Cd3e, Cd3g, Cd8a, Cd8b, Cd4, Il7r, Gzma
# NK cells: Nkg7, Prf1, Gzma, Klrd1
# Macrophage/Monocyte: Adgre1, Cd68, Lyz2, Ccr2
# Neutrophil: S100a8, S100a9, Csf3r
# Dendritic cells: Zbtb46, Flt3
# Fibroblast: Col1a1, Pdgfra
# Endothelial: Pecam1, Cdh5
# Cardiomyocyte: Tnnt2, Myh6
# Mural cells: Tagln, Myh11
# Schwann cells: Mpz, Plp1

# Assign cell type names
new.cluster.ids <- c("Fibroblast", "Macrophage&Monocyte", "Endothelial", "T&NK_cells",
                     "Endothelial", "Neutrophil", "B_cells", "Mural_cells",
                     "Cardiomyocyte", "Dendritic_cells", "Macrophage&Monocyte", "Schwann_cells")
names(new.cluster.ids) <- levels(merged)
merged <- RenameIdents(merged, new.cluster.ids)
merged$cell_type <- Idents(merged)

saveRDS(merged, '~/pan-heart/results/merge/annotated/merged.rds')

# 11. SUBTYPE ANNOTATION - T & NK CELLS----

TNK <- subset(merged, subset = cell_type %in% "T&NK_cells")
DefaultAssay(TNK) <- 'integrated'
TNK <- FindNeighbors(TNK, dims = 1:30, reduction = 'pca')
TNK <- FindClusters(TNK, resolution = 0.05, verbose = TRUE)
TNK <- RunUMAP(TNK, dims = 1:30)

DefaultAssay(TNK) <- 'RNA'
new.cluster.ids <- c('T_cells', 'NK_cells')
names(new.cluster.ids) <- levels(TNK)
TNK <- RenameIdents(TNK, new.cluster.ids)
TNK$cell_type <- Idents(TNK)

saveRDS(TNK, '~/pan-heart/results/merge/annotated/NKT.rds')

# Merge back
mapping_TNK <- data.frame(CellID = colnames(TNK), Subtype = TNK$cell_type)
mapping_original <- data.frame(CellID = colnames(merged), Subtype = merged$cell_type)
mapping <- merge(mapping_original, mapping_TNK, by = "CellID", all.x = TRUE)
mapping$Subtype <- ifelse(!is.na(mapping$Subtype.y), mapping$Subtype.y, mapping$Subtype.x)
merged@meta.data$cell_type <- mapping[match(Cells(merged), mapping$CellID), "Subtype"]

# Reorder cell types
cell_order <- c("B_cells", "Cardiomyocyte", "Dendritic_cells", "Endothelial", "Fibroblast",
                "Macrophage&Monocyte", "Mural_cells", "NK_cells", "Neutrophil", 
                "Schwann_cells", "T_cells")
Idents(merged) <- factor(merged$cell_type, levels = cell_order)

saveRDS(merged, '~/pan-heart/results/merge/annotated/merged.rds')

# 12. FIND ALL MARKERS----

DefaultAssay(merged) <- "RNA"
markers <- FindAllMarkers(merged, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers, "~/pan-heart/results/merge/annotated/markers.csv", row.names = TRUE)

# QC plots
VlnPlot(merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0)