# Correlation between pathway AUC scores and mean gene expression
# Positive correlation: transcriptional regulation drives pathway activity
# Negative correlation: post-transcriptional/translational regulation dominates


## SETUP ----
library(Seurat)
library(ggplot2)
library(dplyr)
library(AUCell)
library(msigdbr)

# Directories
outdir <- "~/pan-heart/results/landscape2025.3.25"
setwd(outdir)
dir.create("./plots/corr", showWarnings = FALSE, recursive = TRUE)
dir.create("./plots/corr2", showWarnings = FALSE, recursive = TRUE)
dir.create("./files", showWarnings = FALSE)

# Color palette
immune_colors <- c("NK cells" = "#f57c6e", "T cells" = "#f2b56f", "B cells" = "#b8aeeb",
                   "Dendritic cells" = "#B0DC66", "Macrophages&Monocytes" = "#88d8db",
                   "Neutrophils" = "#71b7ed")

# Load data
seurat_obj <- readRDS("~/pan-heart/results/subset/immune_cells/2025.3.25/immune_cells.rds")
DefaultAssay(seurat_obj) <- "RNA"
seurat_obj <- subset(seurat_obj, disease != "Normal")


# ============================================================================
# 1. GO PATHWAY ANALYSIS (MSigDB C5)
# ============================================================================

# Load MSigDB C5 gene sets
C5_gene_sets <- msigdbr::msigdbr(species = "Mus musculus", category = "C5") %>% 
  dplyr::select(gs_name, gene_symbol)

# Define pathways of interest
pathway_terms <- c(
  # Inflammation
  "GOBP_ACUTE_INFLAMMATORY_RESPONSE",
  "GOBP_CYTOKINE_PRODUCTION_INVOLVED_IN_INFLAMMATORY_RESPONSE",
  "GOBP_INFLAMMATORY_RESPONSE",
  # IL-17 signaling
  "GOBP_RESPONSE_TO_INTERLEUKIN_17",
  "GOBP_INTERLEUKIN_17_PRODUCTION",
  "GOBP_NEGATIVE_REGULATION_OF_INTERLEUKIN_17_PRODUCTION",
  "GOBP_POSITIVE_REGULATION_OF_INTERLEUKIN_17_PRODUCTION",
  # TNF signaling
  "GOBP_TUMOR_NECROSIS_FACTOR_MEDIATED_SIGNALING_PATHWAY",
  "GOBP_TUMOR_NECROSIS_FACTOR_SUPERFAMILY_CYTOKINE_PRODUCTION",
  "GOMF_TUMOR_NECROSIS_FACTOR_ACTIVATED_RECEPTOR_ACTIVITY",
  "GOBP_RESPONSE_TO_TUMOR_NECROSIS_FACTOR",
  # NF-kB signaling
  "GOBP_POSITIVE_REGULATION_OF_I_KAPPAB_KINASE_NF_KAPPAB_SIGNALING",
  # Glycolysis
  "GOBP_NEGATIVE_REGULATION_OF_GLYCOLYTIC_PROCESS",
  "GOBP_POSITIVE_REGULATION_OF_GLYCOLYTIC_PROCESS",
  "GOBP_REGULATION_OF_GLYCOLYTIC_PROCESS",
  # Gluconeogenesis
  "GOBP_NEGATIVE_REGULATION_OF_GLUCONEOGENESIS",
  "GOBP_POSITIVE_REGULATION_OF_GLUCONEOGENESIS",
  "GOBP_REGULATION_OF_GLUCONEOGENESIS",
  # Oxidative phosphorylation
  "GOBP_NEGATIVE_REGULATION_OF_OXIDATIVE_PHOSPHORYLATION",
  "GOBP_OXIDATIVE_PHOSPHORYLATION",
  "GOBP_REGULATION_OF_OXIDATIVE_PHOSPHORYLATION",
  # Fatty acid metabolism
  "GOBP_FATTY_ACID_BETA_OXIDATION",
  "GOBP_FATTY_ACID_BIOSYNTHETIC_PROCESS",
  # Oxidative stress
  "GOBP_NEGATIVE_REGULATION_OF_RESPONSE_TO_OXIDATIVE_STRESS",
  "GOBP_POSITIVE_REGULATION_OF_RESPONSE_TO_OXIDATIVE_STRESS",
  "GOBP_REGULATION_OF_RESPONSE_TO_OXIDATIVE_STRESS",
  # Immune & chemokine
  "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION",
  "GOBP_CHEMOKINE_PRODUCTION",
  "GOBP_NEGATIVE_REGULATION_OF_CHEMOKINE_MEDIATED_SIGNALING_PATHWAY",
  "GOBP_NEGATIVE_REGULATION_OF_CHEMOKINE_PRODUCTION",
  "GOBP_POSITIVE_REGULATION_OF_CHEMOKINE_PRODUCTION",
  "GOBP_REGULATION_OF_CHEMOKINE_MEDIATED_SIGNALING_PATHWAY",
  "GOBP_RESPONSE_TO_CHEMOKINE"
)

# Filter selected gene sets
selected_gene_sets <- C5_gene_sets %>% filter(gs_name %in% pathway_terms)

# Build rankings and calculate AUC
cells_rankings <- AUCell_buildRankings(as.matrix(seurat_obj@assays$RNA@data))
geneSets <- lapply(unique(selected_gene_sets$gs_name),
                   function(x) selected_gene_sets$gene_symbol[selected_gene_sets$gs_name == x])
names(geneSets) <- unique(selected_gene_sets$gs_name)

cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank = nrow(cells_rankings) * 0.05)
saveRDS(cells_AUC, "files/cells_AUC.rds")

# Plot correlation for each pathway
plot_pathway_correlation <- function(seurat_obj, cells_AUC, pathway_name, gene_set_df, output_dir, colors) {
  # Get AUC scores
  auc_scores <- as.numeric(getAUC(cells_AUC)[pathway_name, ])
  seurat_obj$temp_auc <- auc_scores
  
  # Format pathway name for display
  display_name <- pathway_name %>%
    gsub(pattern = "GOBP_|GOCC_|GOMF_", replacement = "", x = .) %>%
    gsub(pattern = "_", replacement = " ", x = .) %>%
    tolower()
  
  # Get pathway genes
  pathway_genes <- gene_set_df %>%
    filter(gs_name == pathway_name) %>%
    pull(gene_symbol) %>%
    unique()
  
  # Filter to existing genes
  existing_genes <- intersect(pathway_genes, rownames(seurat_obj@assays$RNA@data))
  
  # Calculate mean expression score per cell
  expr_matrix <- seurat_obj@assays$RNA@data[existing_genes, ]
  cell_expr_score <- colMeans(expr_matrix)
  
  # Prepare plot data
  plot_data <- data.frame(
    cell_type = seurat_obj$cell_type,
    auc_score = seurat_obj$temp_auc,
    gene_expr_score = cell_expr_score
  )
  
  # Calculate correlation statistics
  stats <- plot_data %>%
    group_by(cell_type) %>%
    summarise(
      r = cor(gene_expr_score, auc_score, use = "complete.obs"),
      p = cor.test(gene_expr_score, auc_score)$p.value,
      .groups = "drop"
    ) %>%
    mutate(label = paste0("r = ", round(r, 2), "  p = ", signif(p, 2)))
  
  # Position labels
  y_max <- max(plot_data$auc_score, na.rm = TRUE)
  y_min <- min(plot_data$auc_score, na.rm = TRUE)
  delta <- (y_max - y_min) * 0.15
  
  stats <- stats %>%
    arrange(cell_type) %>%
    mutate(
      x = max(plot_data$gene_expr_score, na.rm = TRUE) * 0.3,
      y = y_max + seq(from = delta, by = delta, length.out = n())
    )
  
  # Create plot
  p <- ggplot(plot_data, aes(x = gene_expr_score, y = auc_score, color = cell_type)) +
    geom_point(alpha = 0.6, size = 1) +
    geom_smooth(method = "lm", se = TRUE) +
    scale_color_manual(values = colors) +
    labs(x = "Mean gene expression", y = display_name, color = "Cell type") +
    theme_classic() +
    geom_text(data = stats, aes(x = x, y = y, label = label, color = cell_type),
              size = 3, hjust = 0)
  
  ggsave(filename = file.path(output_dir, paste0(display_name, ".pdf")),
         plot = p, height = 4, width = 5)
  
  return(p)
}

# Run for all pathways
for (term in pathway_terms) {
  plot_pathway_correlation(seurat_obj, cells_AUC, term, selected_gene_sets, 
                           "./plots/corr", immune_colors)
}


# ============================================================================
# 2. CUSTOM GENE SET ANALYSIS
# ============================================================================

# Load custom gene sets
custom_gene_sets <- read_csv("~/pan-heart/data/merged_inflammatory.csv")

# Build rankings and calculate AUC
cells_rankings2 <- AUCell_buildRankings(as.matrix(seurat_obj@assays$RNA@data))
geneSets2 <- lapply(unique(custom_gene_sets$gs_name),
                    function(x) custom_gene_sets$gene_symbol[custom_gene_sets$gs_name == x])
names(geneSets2) <- unique(custom_gene_sets$gs_name)

cells_AUC2 <- AUCell_calcAUC(geneSets2, cells_rankings2, aucMaxRank = nrow(cells_rankings2) * 0.05)
saveRDS(cells_AUC2, "files/cells_AUC2.rds")

# Plot correlation for custom pathways
custom_terms <- names(geneSets2)
for (term in custom_terms) {
  plot_pathway_correlation(seurat_obj, cells_AUC2, term, custom_gene_sets,
                           "./plots/corr2", immune_colors)
}

cat("Pathway correlation analysis completed.\n")
cat("Results saved to: ./plots/corr/ and ./plots/corr2/\n")