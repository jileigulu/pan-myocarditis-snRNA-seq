# Pathway score visualization along trajectory for CD4+ T cells

library(monocle)
library(ggplot2)
library(ggpubr)

setwd('~/pan-heart/results/T_cells2025.5.20')
outdir <- './files/trajectory/monocle2/CD4'
outdir2 <- './plots/trajectory/monocle2/CD4/score'
dir.create(outdir2, showWarnings = FALSE)
dir.create(file.path(outdir2, "boxplot"), showWarnings = FALSE)

# Load data
cds <- readRDS(file.path(outdir, "cds.rds"))
score <- readRDS("~/pan-heart/results/T_cells2025.5.20/files/GSVA/TCfunctionscore.es.matrix.rds")
score <- t(score)
score <- score[Cells(cds), ]

# Add scores to cds
pData(cds) <- cbind(pData(cds), score)

# Clean filename function
clean_filename <- function(x) {
  x <- gsub("[:/\\\\*?\"<>|]", "_", x)
  x <- gsub("[._]+", "_", x)
  x <- gsub("^[_\\.]|[_\\.]$", "", x)
  if (x == "") x <- "unnamed_file"
  return(x)
}

# Trajectory plots with score coloring
for (name in colnames(score)) {
  p <- plot_cell_trajectory(cds, color_by = paste0("`", name, "`"),
                            size = 1, show_backbone = TRUE) +
    scale_color_gradientn(colours = c("#4E659D", "#8B8CC0", "#B7A8CF", "#E7BDC7", "#FECEA1", "#EFA484"))
  ggsave(file.path(outdir2, paste0(clean_filename(name), ".pdf")), p, height = 6, width = 6)
}

# Boxplot by State
data <- pData(cds)
for (name in colnames(score)) {
  p <- ggplot(data, aes(x = State, y = .data[[name]], color = State, fill = State)) +
    geom_boxplot(alpha = 0.6, width = 0.6, color = "black") +
    geom_jitter(size = 0.8, alpha = 0.6, width = 0.2) +
    scale_color_manual(values = c("#f57c6e", "#f2b56f", "#b8aeeb")) +
    scale_fill_manual(values = c("#f57c6e", "#f2b56f", "#b8aeeb")) +
    labs(x = "State", y = paste0(gsub("`", "", name), " Scores")) +
    geom_signif(comparisons = list(c("1", "2"), c("1", "3"), c("2", "3")),
                test = "wilcox.test", step_increase = 0.1, map_signif_level = TRUE)
  ggsave(file.path(outdir2, "boxplot", paste0(clean_filename(name), ".pdf")), p, height = 4, width = 4)
}