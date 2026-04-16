# 1. ODDS RATIO ANALYSIS----

## SETUP ----
library(Seurat)
library(dplyr)
library(tidyverse)
library(data.table)
library(pheatmap)
library(ROGUE)
library(sscVis)
library(plyr)
library(ggpubr)

# Directories
outdir <- "~/pan-heart/results/landscape2025.3.25"
setwd(outdir)
dir.create('./files', showWarnings = FALSE)
set.seed(2025)

# Color palettes
disease_colors <- c("Normal" = "#BCD5F1", "EAM" = "#f57c6e", "ICI-MC" = "#B0DC66", "VMC" = "#EDCCEE")

# Load data
seurat.immune <- readRDS("~/pan-heart/results/subset/immune_cells/2025.3.25/immune_cells.rds")


test.dist.table <- function(count.dist, min.rowSum = 0) {
  count.dist <- count.dist[rowSums(count.dist) >= min.rowSum, , drop = FALSE]
  sum.col <- colSums(count.dist)
  sum.row <- rowSums(count.dist)
  count.dist.tb <- as.data.frame(count.dist)
  setDT(count.dist.tb, keep.rownames = TRUE)
  count.dist.melt.tb <- melt(count.dist.tb, id.vars = "rn")
  colnames(count.dist.melt.tb) <- c("rid", "cid", "count")
  
  count.dist.melt.ext.tb <- as.data.table(ldply(seq_len(nrow(count.dist.melt.tb)), function(i) {
    this.row <- count.dist.melt.tb$rid[i]
    this.col <- count.dist.melt.tb$cid[i]
    this.c <- count.dist.melt.tb$count[i]
    other.col.c <- sum.col[this.col] - this.c
    this.m <- matrix(c(this.c, sum.row[this.row] - this.c, other.col.c, sum(sum.col) - sum.row[this.row] - other.col.c), ncol = 2)
    res.test <- fisher.test(this.m)
    data.frame(rid = this.row, cid = this.col, p.value = res.test$p.value, OR = res.test$estimate)
  }))
  
  count.dist.melt.ext.tb <- merge(count.dist.melt.tb, count.dist.melt.ext.tb, by = c("rid", "cid"))
  count.dist.melt.ext.tb[, adj.p.value := p.adjust(p.value, "BH")]
  return(count.dist.melt.ext.tb)
}

do.tissueDist <- function(cellInfo.tb, meta.cluster, loc, out.prefix, pdf.width = 3, pdf.height = 5, verbose = 0) {
  dir.create(dirname(out.prefix), FALSE, TRUE)
  cellInfo.tb <- data.table(cellInfo.tb)
  cellInfo.tb$meta.cluster <- as.character(meta.cluster)
  cellInfo.tb$loc <- as.factor(loc)
  
  loc.avai.vec <- levels(cellInfo.tb[["loc"]])
  count.dist <- unclass(cellInfo.tb[, table(meta.cluster, loc)])[, loc.avai.vec]
  
  count.dist.melt.ext.tb <- test.dist.table(count.dist)
  p.dist.tb <- dcast(count.dist.melt.ext.tb, rid ~ cid, value.var = "p.value")
  OR.dist.tb <- dcast(count.dist.melt.ext.tb, rid ~ cid, value.var = "OR")
  OR.dist.mtx <- as.matrix(OR.dist.tb[, -1])
  rownames(OR.dist.mtx) <- OR.dist.tb[[1]]
  
  sscVis::plotMatrix.simple(OR.dist.mtx, out.prefix = sprintf("%s.OR.dist", out.prefix),
                            show.number = FALSE, waterfall.row = TRUE, par.warterfall = list(score.alpha = 2, do.norm = TRUE),
                            exp.name = expression(italic(OR)), z.hi = 4, palatte = viridis::viridis(7),
                            pdf.width = pdf.width, pdf.height = pdf.height)
  
  if (verbose == 1) {
    return(list(count.dist.melt.ext.tb = count.dist.melt.ext.tb, p.dist.tb = p.dist.tb,
                OR.dist.tb = OR.dist.tb, OR.dist.mtx = OR.dist.mtx))
  } else {
    return(OR.dist.mtx)
  }
}

# Calculate OR
cellInfo.tb <- seurat.immune@meta.data
OR_results <- do.tissueDist(cellInfo.tb = cellInfo.tb, meta.cluster = cellInfo.tb$cell_type,
                            loc = cellInfo.tb$disease, out.prefix = 'c', pdf.width = 10, pdf.height = 5, verbose = 1)

# Generate OR heatmap
or_data <- as.data.frame(OR_results$OR.dist.tb)
rownames(or_data) <- or_data$rid
or_data <- or_data[, -1]
or_data[or_data > 3] <- 3

p_matrix <- as.data.frame(OR_results$p.dist.tb)
rownames(p_matrix) <- p_matrix$rid
p_matrix <- p_matrix[, -1]
sig_stars <- matrix(ifelse(p_matrix <= 0.01 & or_data >= 1, "*", ""), nrow = nrow(or_data), ncol = ncol(or_data))
colnames(sig_stars) <- colnames(or_data)
rownames(sig_stars) <- rownames(or_data)

pheatmap(as.matrix(or_data), scale = "none", cluster_rows = FALSE, cluster_cols = FALSE,
         display_numbers = sig_stars, fontsize_number = 8, number_color = 'white',
         color = colorRampPalette(c("#F0C1CD", "#F3AA9E", "#EF7454", "#DC0000"))(300),
         cellwidth = 15, cellheight = 15, border_color = "white", fontsize = 6,
         main = "Odds Ratio", angle_col = '45')
ggsave("OR.pdf", width = 3, height = 3)



# 2. ROGUE ANALYSIS (Cell population purity)----


# Extract expression matrix
expr_matrix <- as.data.frame(seurat.immune@assays$RNA@data)
expr_matrix <- matr.filter(expr_matrix, min.cells = 10, min.genes = 10)

# Calculate ROGUE value
rogue_res <- rogue(expr_matrix, labels = seurat.immune$cell_type, 
                   samples = seurat.immune$disease, platform = "UMI", span = 0.9)
saveRDS(rogue_res, "./files/ROGUE.rds")

rogue.boxplot(rogue_res)
ggsave("rogue.pdf", height = 4, width = 3)