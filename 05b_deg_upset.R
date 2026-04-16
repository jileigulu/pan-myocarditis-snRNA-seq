# UpSet plots for visualizing DEG intersections across cell types


## SETUP ----
library(UpSetR)
library(tidyverse)
library(reshape2)

setwd("~/pan-heart/results/deg2025.3.31/deg")
dir.create("./plots", showWarnings = FALSE)

immune_colors <- c("NK cells" = "#f57c6e", "T cells" = "#f2b56f", "B cells" = "#b8aeeb",
                   "Dendritic cells" = "#B0DC66", "Macrophages&Monocytes" = "#88d8db",
                   "Neutrophils" = "#71b7ed")

## HELPER FUNCTIONS ----
prepare_upset_data <- function(deg_df, regulation_type) {
  subset_df <- deg_df[deg_df$regulate == regulation_type, c("celltype", "gene")]
  subset_df$id <- 1:nrow(subset_df)
  
  wide_df <- dcast(subset_df, id ~ celltype, value.var = "gene")
  wide_df <- wide_df[, -1]
  
  # Convert to list format for UpSetR
  result_list <- lapply(wide_df, function(x) x[!is.na(x)])
  return(result_list)
}

plot_upset <- function(data_list, title, bar_color, highlight_sets = NULL) {
  upset(fromList(data_list), nsets = 6, nintersects = 15, order.by = "freq",
        show.numbers = "yes", number.angles = 0, point.size = 3,
        matrix.color = "#b0b9b8", line.size = 1, 
        mainbar.y.label = "Gene Intersections", sets.x.label = "Set Size",
        sets.bar.color = immune_colors[names(data_list)], 
        main.bar.color = bar_color, mb.ratio = c(0.55, 0.45))
}

## EAM ----
deg_eam <- read_csv("files/EAM_deg_TOTAL.csv")

# Upregulated
up_list <- prepare_upset_data(deg_eam, "Upregulated")
pdf("./plots/EAM_upset_up.pdf", width = 8, height = 6)
plot_upset(up_list, "EAM Upregulated", "#E41A1C")
dev.off()

# Downregulated
down_list <- prepare_upset_data(deg_eam, "Downregulated")
pdf("./plots/EAM_upset_down.pdf", width = 8, height = 6)
plot_upset(down_list, "EAM Downregulated", "#579CC7")
dev.off()

## ICI-MC ----
deg_ici <- read_csv("files/ICI-MC_deg_TOTAL.csv")

up_list <- prepare_upset_data(deg_ici, "Upregulated")
pdf("./plots/ICI-MC_upset_up.pdf", width = 8, height = 6)
plot_upset(up_list, "ICI-MC Upregulated", "#E41A1C")
dev.off()

down_list <- prepare_upset_data(deg_ici, "Downregulated")
pdf("./plots/ICI-MC_upset_down.pdf", width = 8, height = 6)
plot_upset(down_list, "ICI-MC Downregulated", "#579CC7")
dev.off()

## VMC ----
deg_vmc <- read_csv("files/VMC_deg_TOTAL.csv")

up_list <- prepare_upset_data(deg_vmc, "Upregulated")
pdf("./plots/VMC_upset_up.pdf", width = 8, height = 6)
plot_upset(up_list, "VMC Upregulated", "#E41A1C")
dev.off()

down_list <- prepare_upset_data(deg_vmc, "Downregulated")
pdf("./plots/VMC_upset_down.pdf", width = 8, height = 6)
plot_upset(down_list, "VMC Downregulated", "#579CC7")
dev.off()
