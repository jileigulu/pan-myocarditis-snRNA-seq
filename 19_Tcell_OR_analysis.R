# Odds Ratio analysis for T cell subtypes

library(Seurat)
library(data.table)
library(pheatmap)
library(sscVis)

setwd("~/pan-heart/results/T_cells2025.5.20")
outdir <- "./files"
outdir2 <- "./plots"
dir.create(outdir2, showWarnings = FALSE)

merged <- readRDS(file.path(outdir, "merged.rds"))

# OR analysis functions
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
    this.m <- matrix(c(this.c, sum.row[this.row] - this.c, 
                       other.col.c, sum(sum.col) - sum.row[this.row] - other.col.c), ncol = 2)
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
  
  if (verbose == 1) {
    return(list(count.dist.melt.ext.tb = count.dist.melt.ext.tb, p.dist.tb = p.dist.tb,
                OR.dist.tb = OR.dist.tb, OR.dist.mtx = OR.dist.mtx))
  } else {
    return(OR.dist.mtx)
  }
}

# Calculate OR
cellInfo.tb <- merged@meta.data
OR_results <- do.tissueDist(cellInfo.tb = cellInfo.tb, meta.cluster = cellInfo.tb$cell_type2,
                            loc = cellInfo.tb$disease, out.prefix = 'c', verbose = 1)

# Prepare heatmap data
or_data <- as.data.frame(OR_results$OR.dist.tb)
rownames(or_data) <- or_data$rid
or_data <- or_data[, -1]
or_data[or_data > 3] <- 3

# Significance matrix
p_matrix <- as.data.frame(OR_results$p.dist.tb)
rownames(p_matrix) <- p_matrix$rid
p_matrix <- p_matrix[, -1]
sig_stars <- matrix(ifelse(p_matrix <= 0.01 & or_data >= 1, "*", ""), 
                    nrow = nrow(or_data), ncol = ncol(or_data))
colnames(sig_stars) <- colnames(or_data)
rownames(sig_stars) <- rownames(or_data)

# Row annotation
annotation_row <- data.frame(
  Subtype = factor(rep(c("CD4+TN", "CD4+Th17", "CD4+Treg", "CD8+TCM", "CD8+TEFF", "CD8+Tisg"), 
                       c(3, 1, 1, 1, 2, 1))),
  Lineage = factor(rep(c("CD4+T Cells", "CD8+T Cells"), c(5, 4)))
)
rownames(annotation_row) <- rownames(or_data)

ann_colors <- list(
  Subtype = c("CD4+TN" = "#f8a8a8", "CD4+Th17" = "#ffc17f", "CD4+Treg" = "#aa9e93",
              "CD8+TCM" = "#ED9B72", "CD8+TEFF" = "#b3e19b", "CD8+Tisg" = "#f36569"),
  Lineage = c("CD4+T Cells" = "#377EB8", "CD8+T Cells" = "#E41A1C")
)

# Plot OR heatmap
p <- pheatmap(as.matrix(or_data), scale = "none", cluster_rows = FALSE, cluster_cols = FALSE,
              display_numbers = sig_stars, fontsize_number = 8, number_color = 'white',
              color = colorRampPalette(c("#F0C1CD", "#F3AA9E", "#EF7454", "#DC0000"))(300),
              cellwidth = 15, cellheight = 15, border_color = "white", fontsize = 6,
              annotation_colors = ann_colors, annotation_row = annotation_row,
              gaps_row = c(5), main = "Odds Ratio")
pdf(file.path(outdir2, "OR.pdf"), width = 4, height = 4)
print(p)
dev.off()