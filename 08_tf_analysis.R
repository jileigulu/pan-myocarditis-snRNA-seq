# Transcription factor analysis for macrophage subtypes using SCENIC
# Includes: TF identification, RSS calculation, regulon activity, target gene enrichment

## SETUP ----
library(Seurat)
library(tidyverse)
library(AUCell)
library(SCENIC)
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(clusterProfiler)
library(org.Mm.eg.db)

setwd("~/pan-heart/results/Macrophage2025.4.16")

# Directories
tf_dir <- "./files/TF"
plot_dir <- "./plots/TF"
enrich_dir <- file.path(tf_dir, "Enrichment")
dir.create(tf_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(enrich_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(enrich_dir, "rds"), showWarnings = FALSE)
dir.create(file.path(enrich_dir, "target_DEG"), showWarnings = FALSE)

# Color palettes
disease_colors <- c("Normal" = "#BCD5F1", "EAM" = "#f57c6e",
                    "ICI-MC" = "#B0DC66", "VMC" = "#EDCCEE")

# 1. RUN SCENIC (Commented - already run)----

# seurat_obj <- readRDS("./files/merged.rds")
# data <- hp_run_pyscenic(x = seurat_obj, species = "mouse", outdir = tf_dir)

# 2. LOAD DATA----

seurat.obj <- readRDS("./files/merged.rds")
DefaultAssay(seurat.obj) <- "RNA"
seurat.obj <- NormalizeData(seurat.obj)

# Load TF-target relationships
tfs_targets <- read.table(file.path(tf_dir, "files", "tfs_targer.tsv"), 
                          sep = ',', header = TRUE)

# Load DEG results
deg <- read_csv("./files/deg/deg1.csv")

# Load regulon AUC matrix
regulonAUC <- importAUCfromText(file.path(tf_dir, "files", "auc_mtx.csv"))
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)), ]
AUC_matrix <- t(getAUC(regulonAUC))

# Add SCENIC assay to Seurat object
seurat.obj[["scenic"]] <- CreateAssayObject(counts = t(getAUC(regulonAUC)))

# 3. IDENTIFY DIFFERENTIALLY EXPRESSED TFs----

# Filter DEGs by disease
deg_by_disease <- list(
  EAM = deg[deg$disease == 'EAM', ],
  ICI.MC = deg[deg$disease == 'ICI-MC', ],
  VMC = deg[deg$disease == 'VMC', ]
)

# Extract TFs for each disease
tf_list <- lapply(deg_by_disease, function(x) {
  x[x$gene %in% tfs_targets$tf, ]
})

# Combine all TF DEGs
deg_tfs <- bind_rows(tf_list)
deg_tfs$regulate <- factor(deg_tfs$regulate, levels = c('Upregulated', 'Downregulated'))

# Separate up and down regulated TFs
up_tfs <- deg_tfs[deg_tfs$regulate == 'Upregulated', ]
down_tfs <- deg_tfs[deg_tfs$regulate == 'Downregulated', ]

# 4. VISUALIZE TF OVERLAP ACROSS DISEASES (Petal Plots)----

plot_tf_petal <- function(tf_df, title, colors, output_file) {
  # Count TF frequency across diseases
  tf_count <- table(tf_df$gene) %>% as.data.frame()
  multi_tf <- tf_df[tf_df$gene %in% tf_count$Var1[tf_count$Freq > 1], ]
  multi_tf$number <- 1
  
  p <- ggplot(multi_tf, aes(x = gene, y = number, fill = disease)) +
    geom_bar(stat = "identity") + 
    scale_fill_manual(values = colors) +
    coord_polar() +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5, color = '#C30D23'),
          axis.text.x = element_text(vjust = 0, hjust = 1),
          panel.background = element_blank(),
          axis.text.y = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank())
  
  ggsave(output_file, p, width = 5, height = 3.8)
}

plot_tf_petal(up_tfs, "Upregulated TFs in >1 diseases", 
              c("EAM" = "#CD534CFF", "ICI-MC" = "#F79C5D", "VMC" = "#8861AC"),
              file.path(plot_dir, "花瓣图up.pdf"))

plot_tf_petal(down_tfs, "Downregulated TFs in >1 diseases",
              c("EAM" = "#CD534CFF", "ICI-MC" = "#F79C5D", "VMC" = "#8861AC"),
              file.path(plot_dir, "花瓣图down.pdf"))

# 5. REGULON SPECIFICITY SCORE (RSS) CALCULATION----

# Prepare cell annotation
cell_annotation <- seurat.obj@meta.data[, c('orig.ident', 'disease'), drop = FALSE]

# Calculate RSS by disease
rss_disease <- calcRSS(AUC = t(getAUC(regulonAUC)),
                       cellAnnotation = cell_annotation[rownames(t(getAUC(regulonAUC))), 'disease'])
rss_disease <- na.omit(rss_disease)
rownames(rss_disease) <- gsub(',', '', rownames(rss_disease))

# Plot RSS heatmap
rss_plot <- plotRSS(rss_disease, zThreshold = 1, cluster_columns = FALSE,
                    order_rows = TRUE, col.low = "darkolivegreen3",
                    col.mid = "#89AD4C", col.high = "darkgreen")
ggsave(file.path(plot_dir, 'RSS_dotplot_disease.pdf'), rss_plot, width = 4, height = 10)

# 6. REGULON ACTIVITY HEATMAP----

# Calculate mean regulon activity by disease
regulon_activity <- sapply(split(rownames(cell_annotation), seurat.obj$disease),
                           function(cells) rowMeans(AUC_matrix[, cells, drop = FALSE]))
rownames(regulon_activity) <- gsub(',', '', rownames(regulon_activity))
regulon_activity <- regulon_activity[, c("Normal", "EAM", "ICI-MC", "VMC")]

# Plot for upregulated TFs
up_tf_selected <- unique(up_tfs$gene)
up_tf_selected <- paste0(up_tf_selected[up_tf_selected %in% gsub("\\(.*", "", rownames(regulon_activity))], "(+)")
up_tf_selected <- up_tf_selected[up_tf_selected %in% rownames(regulon_activity)]

if (length(up_tf_selected) > 0) {
  pdf(file.path(plot_dir, "AUC_heatmap_up.pdf"), width = 5, height = 5)
  pheatmap(regulon_activity[up_tf_selected, ], scale = "row",
           color = colorRampPalette(brewer.pal(9, "RdPu"))(50),
           treeheight_row = 0, treeheight_col = 10, border_color = "white",
           cellheight = 13, cellwidth = 18, fontsize = 10, angle_col = "45",
           cluster_cols = FALSE, cluster_rows = TRUE)
  dev.off()
}

# Plot for downregulated TFs
down_tf_selected <- unique(down_tfs$gene)
down_tf_selected <- paste0(down_tf_selected[down_tf_selected %in% gsub("\\(.*", "", rownames(regulon_activity))], "(+)")
down_tf_selected <- down_tf_selected[down_tf_selected %in% rownames(regulon_activity)]

if (length(down_tf_selected) > 0) {
  pdf(file.path(plot_dir, "AUC_heatmap_down.pdf"), width = 5, height = 10)
  pheatmap(regulon_activity[down_tf_selected, ], scale = "row",
           color = colorRampPalette(brewer.pal(9, "RdPu"))(50),
           treeheight_row = 0, treeheight_col = 10, border_color = "white",
           cellheight = 13, cellwidth = 18, fontsize = 10, angle_col = "45",
           cluster_cols = FALSE, cluster_rows = TRUE)
  dev.off()
}

# 7. RSS RANKING PLOT FOR SPECIFIC TFs----

plot_rss_rank <- function(rss_data, disease_name, highlighted_tfs, output_file) {
  data <- data.frame(rss = rss_data[[disease_name]], TFs = rownames(rss_data))
  data <- data[order(data$rss, decreasing = TRUE), ]
  data$rank <- 1:nrow(data)
  
  highlighted <- data[rownames(data) %in% highlighted_tfs, ]
  
  p <- ggplot(data, aes(x = rank, y = rss)) +
    geom_point(size = 3, color = "grey50", alpha = 0.4) +
    geom_point(data = highlighted, size = 4, shape = 16, color = "#DC050C", alpha = 0.6) +
    ggrepel::geom_text_repel(data = highlighted, aes(label = rownames(highlighted)),
                             color = "black", size = 4, fontface = "italic",
                             segment.size = 0.5, nudge_x = 60, direction = "y") +
    labs(x = 'Rank', y = 'RSS', title = paste0("RSS ", disease_name)) +
    theme_classic() +
    theme(panel.grid = element_blank(),
          axis.title = element_text(size = 15),
          axis.text = element_text(size = 12),
          plot.title = element_text(hjust = 0.5))
  
  ggsave(output_file, p, width = 4, height = 4)
}

highlighted_tfs <- c("Stat1(+)", "Irf7(+)", "Irf5(+)", "Nfil3(+)", "Atf4(+)", "Klf4(+)")
plot_rss_rank(rss_disease, "Myocarditis", highlighted_tfs, 
              file.path(plot_dir, "TF_rank_Myocarditis.pdf"))

# 8. GENE EXPRESSION VALIDATION----


VlnPlot(seurat.obj, features = c("Stat1", "Irf7"), 
        group.by = 'disease', pt.size = 0, cols = disease_colors)
ggsave(file.path(plot_dir, "geneVlnplot_TFs.pdf"), width = 6, height = 3)

# 9. TARGET GENE ENRICHMENT ANALYSIS----

# Helper function to run enrichment for a TF
run_tf_enrichment <- function(tf_name, disease, deg_up_genes, tfs_targets, outdir) {
  # Get target genes that are upregulated in this disease
  target_genes <- tfs_targets %>% 
    filter(tf == tf_name) %>% 
    pull(target_gene)
  
  target_deg <- intersect(target_genes, deg_up_genes)
  
  if (length(target_deg) == 0) {
    message("No target genes found for ", tf_name, " in ", disease)
    return(NULL)
  }
  
  # Save target DEGs
  target_df <- deg[deg$gene %in% target_deg & deg$disease == disease, ]
  write.csv(target_df, file.path(enrich_dir, "target_DEG", 
                                 paste0(tf_name, "_", disease, ".csv")), row.names = FALSE)
  
  # Run GO enrichment
  bp <- enrichGO(target_deg, OrgDb = org.Mm.eg.db, keyType = 'SYMBOL',
                 ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2)
  
  if (is.null(bp) || nrow(bp@result) == 0) {
    message("No enrichment found for ", tf_name, " in ", disease)
    return(NULL)
  }
  
  term <- bp@result
  write.table(term, file = file.path(enrich_dir, paste0(tf_name, "_", disease, "_terms.csv")),
              quote = FALSE, sep = ",", row.names = FALSE)
  saveRDS(term, file.path(enrich_dir, "rds", paste0(tf_name, "_", disease, "_terms.rds")))
  
  # Create bar plot
  df <- term[1:min(30, nrow(term)), ]
  df$labelx <- 0
  df$labely <- seq(nrow(df), 1)
  
  p <- ggplot(df, aes(x = -log10(pvalue), y = reorder(Description, -log10(pvalue)))) +
    geom_bar(stat = "identity", alpha = 1, fill = "#f1bbba", width = 0.8) +
    geom_text(aes(x = labelx, y = labely, label = Description), size = 3.5, hjust = 0) +
    theme_classic() +
    theme(axis.text.y = element_blank(), axis.line.y = element_blank(),
          axis.title.y = element_blank(), axis.ticks.y = element_blank(),
          axis.line.x = element_line(linewidth = 1),
          axis.text.x = element_text(color = 'black', size = 10),
          axis.title.x = element_text(color = 'black', size = 12)) +
    xlab("-log10(pvalue)") +
    ggtitle(paste0(tf_name, "_", disease, "_", length(target_deg))) +
    scale_x_continuous(expand = c(0, 0))
  
  ggsave(file.path(plot_dir, "Enrichment", paste0(tf_name, "_", disease, "_Barplot_TOP30.pdf")),
         p, height = 6, width = 6)
  
  return(term)
}

# Create enrichment output directory
dir.create(file.path(plot_dir, "Enrichment"), showWarnings = FALSE)

# Get upregulated genes by disease
deg_up_by_disease <- list(
  EAM = deg[deg$regulate == "Upregulated" & deg$disease == "EAM", ]$gene,
  `ICI-MC` = deg[deg$regulate == "Upregulated" & deg$disease == "ICI-MC", ]$gene,
  VMC = deg[deg$regulate == "Upregulated" & deg$disease == "VMC", ]$gene
)

# Run for Stat1 across all diseases
for (disease in names(deg_up_by_disease)) {
  run_tf_enrichment("Stat1", disease, deg_up_by_disease[[disease]], tfs_targets, enrich_dir)
}

# Run for Nfil3 and Irf7
run_tf_enrichment("Nfil3", "EAM", deg_up_by_disease[["EAM"]], tfs_targets, enrich_dir)
run_tf_enrichment("Irf7", "VMC", deg_up_by_disease[["VMC"]], tfs_targets, enrich_dir)

# 10. COMPARATIVE GO ENRICHMENT FOR TF TARGETS (Bubble Plot)----

# Helper function for comparative enrichment
compare_tf_target_go <- function(tf_name, deg_up_by_disease, tfs_targets, output_file) {
  # Collect target genes for each disease
  target_list <- list()
  for (disease in names(deg_up_by_disease)) {
    target_genes <- tfs_targets %>%
      filter(tf == tf_name) %>%
      pull(target_gene)
    target_list[[disease]] <- intersect(target_genes, deg_up_by_disease[[disease]])
  }
  
  # Prepare data for compareCluster
  gene_data <- data.frame()
  for (disease in names(target_list)) {
    if (length(target_list[[disease]]) > 0) {
      temp <- data.frame(gene = target_list[[disease]], disease = disease)
      gene_data <- rbind(gene_data, temp)
    }
  }
  
  if (nrow(gene_data) == 0) {
    message("No target genes found for ", tf_name)
    return(NULL)
  }
  
  # Convert to ENTREZ IDs
  gene_ids <- bitr(gene_data$gene, fromType = "SYMBOL", toType = "ENTREZID", 
                   OrgDb = org.Mm.eg.db)
  data_entrez <- merge(gene_ids, gene_data, by.x = "SYMBOL", by.y = "gene")
  
  # Run comparative enrichment
  go_result <- compareCluster(ENTREZID ~ disease, data = data_entrez,
                              fun = "enrichGO", OrgDb = org.Mm.eg.db, ont = "BP",
                              pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
  
  saveRDS(go_result, file.path(enrich_dir, paste0("GO_UP_", tf_name, ".rds")))
  
  # Select specific terms for visualization
  selected_terms <- c('antigen processing and presentation', 'leukocyte mediated immunity',
                      'regulation of T cell activation', 'response to type II interferon',
                      'defense response to virus', 'response to interferon-beta', 'T cell migration')
  
  meta <- go_result@compareClusterResult %>%
    filter(Description %in% selected_terms) %>%
    mutate(logP = -log10(pvalue),
           Cluster = factor(Cluster, levels = c('EAM', 'ICI-MC', 'VMC')),
           Description = factor(Description, levels = selected_terms))
  
  if (nrow(meta) > 0) {
    p <- ggplot(meta, aes(x = Cluster, y = Description)) +
      geom_point(aes(color = logP, size = Count)) +
      scale_size_continuous(range = c(3, 5)) +
      scale_color_gradient(low = "#F7CE68", high = "#FF2525") +
      labs(color = "-log10(pvalue)", x = '', y = '') +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 30, hjust = 1),
            plot.margin = margin(l = 80, t = 30))
    
    ggsave(output_file, p, height = 4.5, width = 7)
  }
  
  return(go_result)
}

# Run comparative enrichment for key TFs
compare_tf_target_go("Stat1", deg_up_by_disease, tfs_targets,
                     file.path(plot_dir, "Enrichment", "heatmap_up_Stat1.pdf"))

compare_tf_target_go("Irf7", deg_up_by_disease, tfs_targets,
                     file.path(plot_dir, "Enrichment", "heatmap_up_Irf7.pdf"))

# 11. TARGET GENE EXPRESSION HEATMAP (Optional)----
# For visualizing expression of top target genes for key TFs

plot_target_gene_heatmap <- function(tf_name, disease, deg_up, tfs_targets, seurat_obj, n_genes = 30) {
  target_genes <- tfs_targets %>%
    filter(tf == tf_name) %>%
    pull(target_gene)
  
  target_deg <- intersect(target_genes, deg_up[[disease]])
  
  if (length(target_deg) > 0) {
    # Get expression data for top target genes
    expr_data <- AverageExpression(seurat_obj, features = target_deg[1:min(n_genes, length(target_deg))],
                                   group.by = "disease")$RNA
    
    pheatmap(expr_data, scale = "row", main = paste(tf_name, "targets in", disease),
             color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
  }
}
