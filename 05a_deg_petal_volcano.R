# DEG petal plots and volcano plots for each disease condition


## SETUP ----
library(Seurat)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggrepel)

setwd("~/pan-heart/results/deg2025.3.31/deg")
dir.create("./plots", showWarnings = FALSE)
dir.create("./plots/volcano", showWarnings = FALSE)
dir.create("./files", showWarnings = FALSE)

# Color palette
immune_colors <- c("NK cells" = "#f57c6e", "T cells" = "#f2b56f", "B cells" = "#b8aeeb",
                   "Dendritic cells" = "#B0DC66", "Macrophages&Monocytes" = "#88d8db",
                   "Neutrophils" = "#71b7ed")

## HELPER FUNCTIONS ----
read_csv_dir <- function(dir_path) {
  file_list <- list.files(path = dir_path, full.names = TRUE)
  df_list <- lapply(file_list, function(x) read.csv(x, sep = ",", header = TRUE))
  do.call(rbind, df_list)
}

filter_degs <- function(deg_df, log2fc_threshold = 0.25, p_adj_threshold = 0.05) {
  deg_df <- deg_df[abs(deg_df$avg_log2FC) >= log2fc_threshold & deg_df$p_val_adj <= p_adj_threshold, ]
  deg_df$regulate <- "No change"
  deg_df[deg_df$avg_log2FC >= log2fc_threshold, "regulate"] <- "Upregulated"
  deg_df[deg_df$avg_log2FC <= -log2fc_threshold, "regulate"] <- "Downregulated"
  return(deg_df)
}

plot_petal <- function(df, disease_name) {
  df_table <- table(df$celltype, df$regulate)
  df_table <- df_table[, c("Downregulated", "Upregulated")]
  df_melt <- melt(df_table, varnames = c("Sample", "Group"), value.name = "Freq")
  
  ggplot(df_melt, aes(x = Sample, y = Freq, fill = Group)) +
    geom_col(width = 0.9, size = 1, alpha = 0.8) +
    labs(x = "The number of DEG", title = disease_name) +
    coord_polar() + theme_bw() +
    theme(aspect.ratio = 1, axis.title.x = element_text(size = 8),
          axis.title.y = element_blank(), axis.text = element_text(colour = "black", size = 8),
          legend.title = element_blank(), legend.position = "top",
          plot.title = element_text(hjust = 0.5)) +
    scale_fill_manual(breaks = c("Upregulated", "Downregulated"),
                      values = c("#E41A1C", "#579CC7"))
}

plot_volcano <- function(deg_df, disease_name, genes_to_label = NULL, colors = immune_colors) {
  volcano_data <- na.omit(deg_df)
  volcano_data$logP <- -log10(volcano_data$p_val)
  volcano_data$label <- ifelse(volcano_data$gene %in% genes_to_label, volcano_data$gene, NA)
  
  ggplot(volcano_data, aes(x = avg_log2FC, y = logP, color = celltype)) +
    geom_point(size = 0.5, alpha = 0.9) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed") +
    geom_text_repel(data = subset(volcano_data, !is.na(label)),
                    aes(label = label, color = celltype), size = 3, fontface = "italic",
                    box.padding = unit(0.5, "lines"), point.padding = unit(0.8, "lines"),
                    segment.color = "black", show.legend = FALSE, max.overlaps = 100000) +
    scale_color_manual(values = colors) +
    labs(x = "log2 FoldChange", y = "-log10 (pvalue)", 
         title = paste0(disease_name, " vs Normal"), color = "Cell Types") +
    theme_bw() + theme(aspect.ratio = 1, legend.position = "top",
                       plot.title = element_text(size = 8, color = "black", hjust = 0.5),
                       axis.text = element_text(colour = "black"), panel.grid = element_blank())
}

## PROCESS EAM ----
deg_eam <- read_csv_dir("./EAM/")
deg_eam_filtered <- filter_degs(deg_eam)
write.csv(deg_eam_filtered, './files/EAM_deg_TOTAL.csv', quote = FALSE, row.names = FALSE)

p_petal <- plot_petal(deg_eam_filtered, "EAM")
ggsave("./plots/EAM_deg.pdf", p_petal, width = 5, height = 5)

# Volcano plot with selected genes
eam_genes <- c("H2-Aa", "Relb", "Nfkb1", "Ly6e", "Bst2", "Ifit3", "Ccl2", "Ccl7")
p_volcano <- plot_volcano(deg_eam_filtered, "EAM", eam_genes)
ggsave("./plots/volcano/EAM.pdf", p_volcano, width = 6, height = 6)

## PROCESS ICI-MC ----
deg_ici <- read_csv_dir("./ICI-MC/")
deg_ici_filtered <- filter_degs(deg_ici)
write.csv(deg_ici_filtered, './files/ICI-MC_deg_TOTAL.csv', quote = FALSE, row.names = FALSE)

p_petal <- plot_petal(deg_ici_filtered, "ICI-MC")
ggsave("./plots/ICI-MC_deg.pdf", p_petal, width = 5, height = 5)

ici_genes <- c("B2m", "Cd27", "Cd28", "Akirin1", "Calr", "Ccr7", "Atp5c1", "Arf6")
p_volcano <- plot_volcano(deg_ici_filtered, "ICI-MC", ici_genes)
ggsave("./plots/volcano/ICI-MC.pdf", p_volcano, width = 6, height = 6)

## PROCESS VMC ----
deg_vmc <- read_csv_dir("./VMC/")
deg_vmc_filtered <- filter_degs(deg_vmc)
write.csv(deg_vmc_filtered, './files/VMC_deg_TOTAL.csv', quote = FALSE, row.names = FALSE)

p_petal <- plot_petal(deg_vmc_filtered, "VMC")
ggsave("./plots/VMC_deg.pdf", p_petal, width = 5, height = 5)

vmc_genes <- c("Apobec3", "Bst2", "Ccr7", "Gbp2", "Gbp7", "Bcl3")
p_volcano <- plot_volcano(deg_vmc_filtered, "VMC", vmc_genes)
ggsave("./plots/volcano/VMC.pdf", p_volcano, width = 6, height = 6)

