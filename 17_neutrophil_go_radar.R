# GO enrichment radar and bar plots for neutrophils

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(fmsb)

setwd("~/pan-heart/results/Neutrophils2025.5.13")
outdir2 <- "./plots/GO"
dir.create(outdir2, showWarnings = FALSE)

# Load GO results
go_data <- list()
for (d in c("EAM", "ICI-MC", "VMC")) {
  go_data[[d]] <- list(
    up = read_csv(paste0("files/GO/Upregulated/", d, "/", d, ".csv")),
    down = read_csv(paste0("files/GO/Downregulated/", d, "/", d, ".csv"))
  )
}

# Upregulated pathways - radar plot----

common_up <- Reduce(intersect, list(go_data$EAM$up$Description[1:100],
                                    go_data$`ICI-MC`$up$Description[1:100],
                                    go_data$VMC$up$Description[1:100]))

selected_up <- c("defense response to virus", "response to type II interferon",
                 "positive regulation of defense response", "cytokine-mediated signaling pathway",
                 "response to interferon-beta", "interferon-mediated signaling pathway",
                 "response to lipopolysaccharide", "type I interferon-mediated signaling pathway",
                 "biological process involved in interaction with host", "negative regulation of viral genome replication")

# Prepare radar data
radar_up <- data.frame()
for (d in names(go_data)) {
  df <- go_data[[d]]$up %>% filter(Description %in% selected_up)
  df <- df[, c("Description", "p.adjust")]
  colnames(df)[2] <- d
  if (nrow(radar_up) == 0) radar_up <- df else radar_up <- full_join(radar_up, df, by = "Description")
}
radar_up[is.na(radar_up)] <- 1
radar_mat <- as.matrix(radar_up[, -1])
rownames(radar_mat) <- radar_up$Description
radar_mat <- -log(radar_mat)

# Add max/min rows
radar_plot <- rbind(rep(max(radar_mat), ncol(radar_mat)),
                    rep(0, ncol(radar_mat)),
                    rep(-log(0.05), ncol(radar_mat)),
                    t(radar_mat))
rownames(radar_plot) <- c("max", "min", "p", rownames(radar_mat))

pdf(file.path(outdir2, 'up_radar.pdf'), width = 11, height = 7)
radarchart(radar_plot, axistype = 1, pcol = c("#f57c6e", "#B0DC66", "#EDCCEE"),
           pfcol = scales::alpha(c("#f57c6e", "#B0DC66", "#EDCCEE"), 0.5),
           plwd = 3, cglcol = "grey60", vlcex = 1)
legend("topright", legend = c("EAM", "ICI-MC", "VMC"), col = c("#f57c6e", "#B0DC66", "#EDCCEE"), pch = 15)
dev.off()

# Downregulated pathways - radar plot----

common_down <- Reduce(intersect, list(go_data$EAM$down$Description[1:200],
                                      go_data$`ICI-MC`$down$Description[1:200],
                                      go_data$VMC$down$Description[1:200]))

selected_down <- c("regulation of apoptotic signaling pathway", "negative regulation of phosphorylation",
                   "regulation of leukocyte migration", "negative regulation of cell adhesion",
                   "myeloid leukocyte differentiation", "regulation of hemopoiesis",
                   "receptor-mediated endocytosis", "regulation of leukocyte mediated immunity",
                   "response to endoplasmic reticulum stress", "regulation of protein serine/threonine kinase activity")

radar_down <- data.frame()
for (d in names(go_data)) {
  df <- go_data[[d]]$down %>% filter(Description %in% selected_down)
  df <- df[, c("Description", "p.adjust")]
  colnames(df)[2] <- d
  if (nrow(radar_down) == 0) radar_down <- df else radar_down <- full_join(radar_down, df, by = "Description")
}
radar_down[is.na(radar_down)] <- 1
radar_mat <- as.matrix(radar_down[, -1])
rownames(radar_mat) <- radar_down$Description
radar_mat <- -log(radar_mat)

radar_plot <- rbind(rep(max(radar_mat), ncol(radar_mat)),
                    rep(0, ncol(radar_mat)),
                    rep(-log(0.05), ncol(radar_mat)),
                    t(radar_mat))
rownames(radar_plot) <- c("max", "min", "p", rownames(radar_mat))

pdf(file.path(outdir2, 'down_radar.pdf'), width = 11, height = 7)
radarchart(radar_plot, axistype = 1, pcol = c("#f57c6e", "#B0DC66", "#EDCCEE"),
           pfcol = scales::alpha(c("#f57c6e", "#B0DC66", "#EDCCEE"), 0.5),
           plwd = 3, cglcol = "grey60", vlcex = 1)
legend("topright", legend = c("EAM", "ICI-MC", "VMC"), col = c("#f57c6e", "#B0DC66", "#EDCCEE"), pch = 15)
dev.off()