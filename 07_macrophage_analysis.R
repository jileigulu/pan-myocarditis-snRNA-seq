# Comprehensive analysis of macrophage/monocyte subtypes
# Includes: UMAP, markers, cell proportion, GO enrichment, GSEA, Module scoring

## SETUP ----
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(ggalluvial)
library(pheatmap)
library(clusterProfiler)
library(org.Mm.eg.db)
library(GSEABase)
library(msigdbr)
library(fmsb)
library(RColorBrewer)

setwd("~/pan-heart/results/Macrophage2025.4.16")
dir.create("./files", showWarnings = FALSE)
dir.create("./plots", showWarnings = FALSE)

# Color palettes
cell_type_cols <- c("Monocytes" = "#A6CEE3",
                    "Macro_CCR2+MHCIIhi" = "#FF6A6A",
                    "Macro_TLF+" = "#89CB6C",
                    "Macro_CCR2-MHCIIhi" = "#FFA500")

disease_colors <- c("Normal" = "#BCD5F1",
                    "EAM" = "#f57c6e",
                    "ICI-MC" = "#B0DC66",
                    "VMC" =  "#EDCCEE")

# 1. LOAD DATA & DIMENSIONALITY REDUCTION----

merged <- readRDS("~/pan-heart/results/Macrophage2025.4.16/files/merged.rds")
Idents(merged) <- factor(Idents(merged), 
                         levels = c("Monocytes", "Macro_TLF+", 
                                    "Macro_CCR2-MHCIIhi", "Macro_CCR2+MHCIIhi"))

# UMAP plot
DimPlot(merged, label = TRUE, pt.size = 0, cols = cell_type_cols) +
  NoLegend() + labs(x = "UMAP1", y = "UMAP2", title = "Macrophages & Monocytes") +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1),
        axis.text = element_blank(), axis.ticks = element_blank())
ggsave("./plots/dimplot.pdf", width = 5, height = 5)

# 2. MARKER GENES----

DefaultAssay(merged) <- 'RNA'
markers <- FindAllMarkers(merged, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write_csv(markers, "./files/markers.csv")

# Marker genes for dot plot
marker_genes <- c("Ccr2", "Ass1", "Fn1", "Clec4e", "Cxcl1", "Cxcl2", 
                  "Cd83", "C1qa", "Lyve1", "Folr2", "Ly6c2", "Napsa", "Ifitm6")

source("~/tutorial/custom_plot_function.R")
p <- DotPlot(merged, dot.scale = 4, col.min = -1, col.max = 1, features = marker_genes) +
  scale_color_distiller(palette = "RdBu", direction = -1) +
  coord_fixed(ratio = 1) + theme_cat() +
  theme(axis.title = element_blank(), legend.margin = margin(b = -8),
        legend.position = "top", axis.text.x = element_text(face = "italic")) +
  guides(x = guide_axis(angle = 90), color = guide_colorbar(title = "Expression")) +
  scale_y_discrete(limits = rev)
ggsave("./plots/markers_dotplot.pdf", width = 10, height = 4)

# 3. CELL PROPORTION ANALYSIS----

## 3.1 Bar plot - cell numbers
meta_count <- merged@meta.data %>% count(cell_type2) %>%
  arrange(n) %>% mutate(cell_type2 = factor(cell_type2, levels = cell_type2))

ggplot(meta_count, aes(x = cell_type2, y = n, fill = cell_type2)) +
  geom_bar(stat = "identity", width = 0.9, color = "white", show.legend = FALSE) +
  labs(x = "Cell Types", y = "Cell Numbers") + coord_flip() +
  scale_fill_manual(values = cell_type_cols) +
  theme(axis.line = element_line(size = 0.5), panel.grid = element_blank(),
        panel.background = element_blank()) +
  scale_y_continuous(expand = c(0, 0))
ggsave("./plots/cell_number.pdf", width = 7, height = 5)

## 3.2 Sankey plot - proportion by disease
prop_df <- as.data.frame(prop.table(table(Idents(merged), merged$disease), margin = 2) * 100)
prop_df$Var2 <- factor(prop_df$Var2, levels = rev(c("Normal", "EAM", "ICI-MC", "VMC")))

ggplot(prop_df, aes(x = Var2, y = Freq, fill = Var1, stratum = Var1, alluvium = Var1)) +
  geom_col(width = 0.5, color = 'black') + geom_flow(width = 0.5, alpha = 0.4) +
  theme_bw() + labs(x = 'Disease', y = 'Fraction', fill = "Cell Type") +
  coord_flip() + scale_fill_manual(values = cell_type_cols) +
  theme(panel.grid = element_blank())
ggsave("./plots/cell_percentage.pdf", width = 7, height = 5)

## 3.3 Boxplot by cell type
prop_box <- as.data.frame(prop.table(table(merged$cell_type2, merged$orig.ident), margin = 2))
colnames(prop_box) <- c('cell_type2', 'sample', 'prop')

sample_info <- merged@meta.data[, c("orig.ident", "disease")] %>% distinct()
prop_box <- merge(prop_box, sample_info, by.x = "sample", by.y = "orig.ident")
prop_box$disease <- factor(prop_box$disease, levels = c("Normal", "EAM", "ICI-MC", "VMC"))

dir.create("./plots/cell_percentage", showWarnings = FALSE)
for (ct in unique(prop_box$cell_type2)) {
  p <- ggplot(prop_box[prop_box$cell_type2 == ct, ], aes(x = disease, y = prop, fill = disease)) +
    geom_boxplot(width = 0.5, outlier.size = 0, show.legend = FALSE) +
    scale_fill_manual(values = disease_colors) + labs(title = ct) + theme_classic()
  filename <- gsub("\\+", "pos", ct) %>% gsub("-", "neg", .)
  ggsave(file.path("./plots/cell_percentage", paste0(filename, '_prop.pdf')), p, width = 4, height = 3)
}

# 4. GO ENRICHMENT ANALYSIS (Disease dimension)----

go_dir <- "./files/GO/整体细胞类型"
dir.create(go_dir, recursive = TRUE, showWarnings = FALSE)

# Function to read GO results
read_go <- function(disease, regulation) {
  path <- file.path(go_dir, regulation, disease, paste0(disease, ".csv"))
  if (file.exists(path)) read_csv(path) else NULL
}

# Radar plot function
plot_go_radar <- function(data_list, title, colors, output_file) {
  common_desc <- Reduce(intersect, lapply(data_list, function(x) x$Description[1:100]))
  
  plot_data <- lapply(names(data_list), function(nm) {
    df <- data_list[[nm]] %>% filter(Description %in% common_desc)
    df <- df[, c("Description", "p.adjust")]
    colnames(df)[2] <- nm
    return(df)
  }) %>% reduce(full_join, by = "Description")
  
  plot_data[is.na(plot_data)] <- 1
  plot_data <- as.data.frame(t(plot_data[, -1]))
  colnames(plot_data) <- common_desc
  plot_data <- -log(plot_data)
  
  # Prepare radar data
  radar_data <- rbind(rep(max(plot_data), ncol(plot_data)),
                      rep(0, ncol(plot_data)),
                      rep(-log(0.05), ncol(plot_data)),
                      plot_data)
  rownames(radar_data) <- c("max", "min", "p", rownames(plot_data))
  
  pdf(output_file, width = 11, height = 7)
  radarchart(radar_data, axistype = 1, pcol = colors, plwd = 3,
             pfcol = scales::alpha(colors, 0.5), cglcol = "grey60",
             vlcex = 1, vlabels = colnames(radar_data))
  legend("topright", legend = names(data_list), col = colors, pch = 15, bty = "n")
  dev.off()
}

# Run for upregulated pathways
go_up <- list(EAM = read_go("EAM", "Upregulated"),
              ICI.MC = read_go("ICI-MC", "Upregulated"),
              VMC = read_go("VMC", "Upregulated"))
go_up <- go_up[!sapply(go_up, is.null)]

plot_go_radar(go_up, "Upregulated Pathways", c("#f57c6e", "#B0DC66", "#EDCCEE"),
              "./plots/up_radar.pdf")

# 5. GSEA ANALYSIS----

gsea_dir <- "./files/GSEA"
dir.create(gsea_dir, recursive = TRUE, showWarnings = FALSE)

# Selected pathways for visualization
selected_pathways <- c("Interferon gamma response", "Oxidative phosphorylation",
                       "Glycolysis", "Interferon alpha response", "Reactive oxygen species pathway",
                       "Inflammatory response", "Apoptosis", "Angiogenesis", "Tgf beta signaling",
                       "Tnfa signaling via nfkb", "Complement", "Pi3k akt mtor signaling",
                       "Il6 jak stat3 signaling")

# Collect GSEA results across diseases
gsea_data <- data.frame()

for (disease in c("EAM", "ICI-MC", "VMC")) {
  # Load C5 results
  c5_path <- file.path(gsea_dir, paste0(disease, "_C5_result.rds"))
  if (file.exists(c5_path)) {
    res <- readRDS(c5_path)
    res_df <- res@result
    res_df$Description <- gsub("GOBP_|GOCC_|GOMF_", "", res_df$Description) %>%
      gsub("_", " ", .) %>% str_to_sentence()
    gsea_data <- rbind(gsea_data, data.frame(Description = res_df$Description,
                                             NES = res_df$NES, pvalue = res_df$pvalue,
                                             group = disease))
  }
}

# Filter and plot
gsea_filtered <- gsea_data[gsea_data$Description %in% selected_pathways, ]
gsea_filtered$logp <- -log10(gsea_filtered$pvalue)
gsea_filtered$Description <- factor(gsea_filtered$Description, 
                                    levels = rev(selected_pathways))

p <- ggplot(gsea_filtered, aes(x = NES, y = Description)) +
  geom_vline(xintercept = 0, color = "gray70", linetype = "dashed") +
  geom_segment(aes(x = 0, xend = NES, y = Description, yend = Description), size = 1.2, color = "gray80") +
  geom_point(aes(size = logp, color = group), alpha = 0.9) +
  scale_size(range = c(2, 10)) +
  scale_color_manual(values = c("EAM" = "#f57c6e", "ICI-MC" = "#B0DC66", "VMC" = "#EDCCEE")) +
  labs(x = "Normalized Enrichment Score (NES)", y = NULL, 
       title = "Pathway Enrichment by Disease", size = "-log10(p)", color = "Disease") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1),
                          panel.grid.major.y = element_blank())
ggsave("./plots/GSEA.pdf", p, height = 6, width = 7)

# 6. ADD MODULE SCORE & RADAR PLOTS----

module_dir <- "./files/Addmodulescore"
dir.create(module_dir, showWarnings = FALSE)

# Define pathways for module scoring
pathways_module <- c("acute inflammatory response", "oxidative phosphorylation",
                     "activated t cell proliferation", "regulation of type i interferon mediated signaling pathway")

# Get gene sets from MSigDB
c5_msigdb <- msigdbr(species = "Mus musculus", category = "C5")
c5_msigdb$gs_name <- c5_msigdb$gs_name %>% gsub("GOBP_", "", .) %>% 
  gsub("_", " ", .) %>% str_to_lower()

gene_sets <- lapply(pathways_module, function(p) {
  c5_msigdb[c5_msigdb$gs_name == p, "gene_symbol"]
})
names(gene_sets) <- pathways_module

# Add module scores
seurat_scored <- AddModuleScore(merged, features = gene_sets, ctrl = 5, name = "Score")
for (i in 1:length(gene_sets)) {
  colnames(seurat_scored@meta.data)[colnames(seurat_scored@meta.data) == paste0("Score", i)] <- pathways_module[i]
}
saveRDS(seurat_scored, file.path(module_dir, "seurat_scored.rds"))

# Radar plot for each pathway across cell types
for (pathway in pathways_module) {
  score_data <- seurat_scored@meta.data[, c('cell_type2', 'disease', pathway)]
  agg_data <- aggregate(score_data[[pathway]], by = list(score_data$disease, score_data$cell_type2), mean)
  colnames(agg_data) <- c('disease', 'cell_type', pathway)
  agg_data <- agg_data[agg_data$disease %in% c('EAM', 'ICI-MC', 'VMC'), ]
  
  wide_data <- agg_data %>% pivot_wider(id_cols = disease, names_from = cell_type, values_from = pathway)
  rownames(wide_data) <- wide_data$disease
  wide_data <- wide_data[, -1]
  
  radar_input <- rbind(rep(max(wide_data), ncol(wide_data)),
                       rep(min(wide_data), ncol(wide_data)),
                       wide_data)
  rownames(radar_input) <- c("max", "min", rownames(wide_data))
  
  pdf(file.path("./plots", paste0(gsub(" ", "_", pathway), "_radar.pdf")), width = 7, height = 5)
  radarchart(radar_input, axistype = 1, pcol = disease_colors[2:4],
             pfcol = scales::alpha(disease_colors[2:4], 0.3), plwd = 2,
             cglcol = "grey", vlcex = 0.8)
  title(main = str_to_title(pathway))
  dev.off()
}

# 7. PATHWAY GENE HEATMAP----

# Define genes for each pathway
pathway_genes <- list(
  `Tnfa signaling via nfkb` = c("Cd83", "Cxcl1", "Nfkb1", "Tnf", "Tnfaip3"),
  `Tgf beta signaling` = c("Tgfb1", "Junb", "Ctnnb1", "Id2", "Id3", "Tgif1"),
  `Pi3k akt mtor signaling` = c("Il2rg", "Slc2a1", "Sqstm1", "Eif4e"),
  `Il6 jak stat3 signaling` = c("Socs3", "Il10rb", "Cxcl2", "Cd14"),
  `Interferon gamma response` = c("Stat1", "Cd274", "Tap1", "Socs1", "Nfkbia"),
  `Interferon alpha response` = c("Cxcl10", "Ifitm1", "Ly6e", "Psmb8", "Psme1")
)

# Load DEG data
deg <- read_csv("./files/deg/deg.csv")

# Prepare heatmap data
all_genes <- unlist(pathway_genes)
deg_selected <- deg[deg$gene %in% all_genes, c("gene", "disease", "avg_log2FC")]
heatmap_data <- deg_selected %>% pivot_wider(names_from = "disease", values_from = "avg_log2FC")
heatmap_mat <- as.matrix(heatmap_data[, -1])
rownames(heatmap_mat) <- heatmap_data$gene
heatmap_mat[is.na(heatmap_mat)] <- 0

# Create row annotation
row_annot <- data.frame(Pathway = factor(rep(names(pathway_genes), times = sapply(pathway_genes, length)),
                                         levels = names(pathway_genes)))
rownames(row_annot) <- all_genes

# Plot heatmap
pdf("./plots/gene_heatmap.pdf", height = 8, width = 6)
pheatmap(heatmap_mat, show_colnames = TRUE, show_rownames = TRUE,
         annotation_row = row_annot, cluster_rows = FALSE, cluster_cols = FALSE,
         color = colorRampPalette(c("#2166AC", "#F0F0F0", "#C8191D"))(100),
         cellwidth = 30, cellheight = 12, fontsize = 8)
dev.off()

# 8. DEG PETAL PLOT----

deg_filtered <- deg[abs(deg$avg_log2FC) >= 0.25 & deg$p_val_adj <= 0.05, ]
deg_filtered$regulate <- ifelse(deg_filtered$avg_log2FC >= 0.25, "Upregulated", "Downregulated")

deg_count <- table(deg_filtered$disease, deg_filtered$regulate) %>% 
  as.data.frame() %>% filter(Freq > 0)
colnames(deg_count) <- c("Sample", "Group", "Freq")

p <- ggplot(deg_count, aes(x = Sample, y = Freq, fill = Group)) +
  geom_col(width = 0.9, alpha = 0.8) + coord_polar() +
  labs(x = "Number of DEG") + theme_bw() +
  theme(aspect.ratio = 1, axis.title.y = element_blank(),
        axis.text = element_text(size = 8), legend.position = "top") +
  scale_fill_manual(values = c("Upregulated" = "#E41A1C", "Downregulated" = "#579CC7"))
ggsave("./plots/deg_petal.pdf", p, width = 5, height = 5)

