# GSEA analysis for neutrophils using custom gene sets

library(Seurat)
library(clusterProfiler)
library(msigdbr)
library(dplyr)
library(stringr)
library(readr)
library(tidyr)

setwd("~/pan-heart/results/Neutrophils2025.5.13")

outdir <- './files/GSEA/'
outdir2 <- './plots/GSEA/'
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(outdir2, recursive = TRUE, showWarnings = FALSE)

disease <- c('EAM', 'ICI-MC', 'VMC')
set.seed(2025)

# Load custom gene sets
db <- read_csv("~/pan-heart/data/neutro_gsva.csv")

db_long <- db %>%
  pivot_longer(cols = everything(), names_to = "gs_name", 
               values_to = "gene_symbol", values_drop_na = TRUE) %>%
  arrange(gs_name, .by_group = TRUE)

genesets <- db_long

# GSEA analysis for each disease
input.dir <- "~/pan-heart/results/deg2025.3.31/gsea"

for (i in disease) {
  gsea.input <- read_csv(file.path(input.dir, i, 'Neutrophils_gsea.csv'))
  
  # Sort by log2FC descending
  gsea.input <- gsea.input[order(gsea.input$avg_log2FC, decreasing = TRUE), ]
  rownames(gsea.input) <- gsea.input$gene
  genelist <- structure(gsea.input$avg_log2FC, names = rownames(gsea.input))
  
  # Run GSEA
  res <- GSEA(genelist, TERM2GENE = genesets, pvalueCutoff = 1, eps = 0)
  
  gsea_result <- as.data.frame(res@result)
  gsea_result <- gsea_result[order(gsea_result$pvalue, decreasing = FALSE), ]
  gsea_result$Description <- gsub('_', ' ', gsea_result$Description)
  gsea_result$Description <- str_to_sentence(gsea_result$Description)
  
  # Output files
  write.table(genelist, quote = FALSE, sep = "\t", col.names = FALSE, 
              file.path(outdir, paste0(i, '.rnk')))
  write.table(gsea_result, quote = FALSE, sep = "\t", col.names = TRUE, 
              file.path(outdir, paste0(i, '_result.txt')))
  write.csv(gsea_result, file.path(outdir, paste0(i, '_table_result.csv')), row.names = FALSE)
  
  saveRDS(res, file.path(outdir, paste0(i, '_result.rds')))
  write.csv(res, file.path(outdir, paste0(i, '_table.csv')), row.names = FALSE)
}

# GSEA analysis for neutrophils 

library(Seurat)
library(clusterProfiler)
library(msigdbr)
library(dplyr)
library(stringr)
library(readr)
library(ggpubr)

setwd("~/pan-heart/results/Neutrophils2025.5.13")
outdir <- './files/GSEA/'
outdir2 <- './plots/GSEA/'
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(outdir2, recursive = TRUE, showWarnings = FALSE)

disease <- c('EAM', 'ICI-MC', 'VMC')
set.seed(2025)

# Helper function to clean pathway names
clean_pathway_names <- function(names_vec) {
  names_clean <- sub("^[A-Z]+_", "", names_vec)
  names_clean <- gsub("_", " ", names_clean)
  names_clean <- tolower(names_clean)
  names_clean <- gsub("\\b(\\w)", "\\U\\1", names_clean, perl = TRUE)
  return(names_clean)
}

# Run GSEA for each disease and category
input.dir <- "~/pan-heart/results/deg2025.3.31/gsea"
categories <- c("H", "C5", "C2")

for (i in disease) {
  gsea.input <- read_csv(file.path(input.dir, i, 'Neutrophils_gsea.csv'))
  gsea.input <- gsea.input[order(gsea.input$avg_log2FC, decreasing = TRUE), ]
  genelist <- structure(gsea.input$avg_log2FC, names = gsea.input$gene)
  
  for (cat in categories) {
    genesets <- msigdbr(species = "Mus musculus", category = cat)
    genesets <- subset(genesets, select = c("gs_name", "gene_symbol"))
    
    if (cat == "C2") {
      genelist <- genelist[!is.na(genelist)]
      genelist <- genelist[!duplicated(names(genelist))]
      genelist <- sort(genelist, decreasing = TRUE)
      genelist <- genelist[genelist != 0]
    }
    
    res <- GSEA(genelist, TERM2GENE = genesets, pvalueCutoff = ifelse(cat == "C5", 0.1, 1), eps = 0)
    gsea_result <- as.data.frame(res@result)
    gsea_result <- gsea_result[order(gsea_result$pvalue, decreasing = FALSE), ]
    gsea_result$Description <- gsub('_', ' ', gsea_result$Description)
    gsea_result$Description <- str_to_sentence(gsea_result$Description)
    
    saveRDS(res, file.path(outdir, paste0(i, '_', cat, '_result.rds')))
    write.csv(gsea_result, file.path(outdir, paste0(i, '_', cat, '_table_result.csv')), row.names = FALSE)
  }
}

# Plot selected pathways
selected_features <- c("Neutrophil Degranulation", "Fc Gamma R Mediated Phagocytosis",
                       "Glycolysis", "Hypoxia Metagene", "Programmed Cell Death", "Apoptosis",
                       "Interferon Alpha Beta Signaling", "Interferon Signaling", "Interferon Responsive Genes",
                       "Interferon Gamma Signaling", "Signaling By Interleukins", "Tnf Signaling Not Via Nfkb",
                       "Nfkb Targets Keratinocyte Up", "Response To Lps With Mechanical Ventilation",
                       "Influenza Infection", "Class I Mhc Mediated Antigen Processing Presentation")

# Collect data
plot_data <- data.frame()
for (i in disease) {
  res <- readRDS(file.path(outdir, paste0(i, '_C2_result.rds')))
  res@result$Description <- clean_pathway_names(res@result$Description)
  plot_data <- rbind(plot_data, res@result[res@result$Description %in% selected_features, c("Description", "NES", "pvalue")])
  plot_data$group <- i
}

plot_data$logp <- -log10(plot_data$pvalue)
plot_data$Description <- factor(plot_data$Description, levels = rev(selected_features))
disease_colors <- c("#f57c6e", "#B0DC66", "#EDCCEE")

p <- ggplot(plot_data, aes(x = NES, y = Description)) +
  geom_vline(xintercept = 0, color = "gray70", linetype = "dashed", size = 0.8) +
  geom_segment(aes(x = 0, xend = NES, y = Description, yend = Description), size = 1.2, color = "gray80") +
  geom_point(aes(size = logp, color = group), alpha = 0.9) +
  scale_size(range = c(2, 10)) + scale_color_manual(values = disease_colors) +
  labs(x = "Normalized Enrichment Score (NES)", y = NULL, 
       title = "Pathway Enrichment by Disease", size = "-log10(p)", color = "Disease") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1),
                          panel.grid.major.y = element_blank())
ggsave(file.path(outdir2, 'GSEA.pdf'), p, height = 7, width = 7)