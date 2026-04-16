# Basic analysis for T cells: UMAP, markers, cell proportion, OR value

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggalluvial)
library(pheatmap)
library(data.table)
library(sscVis)

setwd("~/pan-heart/results/T_cells2025.5.20")
outdir <- "./files"
outdir2 <- "./plots"
dir.create(outdir, showWarnings = FALSE)
dir.create(outdir2, showWarnings = FALSE)

# Load data
merged <- readRDS(file.path(outdir, "merged.rds"))

# Colors
cell_type_cols <- c("#f8a8a8", "#c1b5d9", "#97d2ce", "#aa9e93", "#ffc17f", 
                    "#b3e19b", "#67a8cd", "#ED9B72", "#f36569")
names(cell_type_cols) <- levels(merged)

# 1. UMAP plot----

DimPlot(merged, label = TRUE, pt.size = 1, cols = cell_type_cols) +
  NoLegend() + labs(x = "UMAP1", y = "UMAP2", title = "T Cells") +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1),
        axis.text = element_blank(), axis.ticks = element_blank())
ggsave(file.path(outdir2, "dimplot.pdf"), width = 5, height = 5)

# 2. Marker genes dot plot----

DefaultAssay(merged) <- 'RNA'
markers <- FindAllMarkers(merged, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write_csv(markers, file.path(outdir, "markers.csv"))

genes <- c('Cd4', 'Il7r', 'Lef1', 'Ccr7', 'Sell', 'Tcf7', 'Rpl12', 'Foxp3', 'S100a4', 'Il2ra',
           'Il17a', 'Rorc', 'Gadd45b', 'Cd8a', 'Cd8b1', 'Gzmb', 'Ifng', 'Klrg1', 'Klrd1',
           'Ccl5', 'Gzma', 'Tpi1', 'Isg15', 'Isg20', 'Ifit1')

source("~/tutorial/custom_plot_function.R")
p <- DotPlot(merged, dot.scale = 4, col.min = -1, col.max = 1, features = genes) +
  scale_color_distiller(palette = "RdBu", direction = -1) +
  coord_fixed(ratio = 1) + theme_cat() +
  theme(axis.title = element_blank(), legend.margin = margin(b = -8),
        legend.position = "top", axis.text.x = element_text(face = "italic")) +
  guides(x = guide_axis(angle = 90), color = guide_colorbar(title = "Expression")) +
  scale_y_discrete(limits = rev)
ggsave(file.path(outdir2, "markers_dotplot.pdf"), width = 10, height = 4)

# 3. Cell proportion----

# Bar plot
meta_count <- merged@meta.data %>% count(cell_type2) %>%
  arrange(n) %>% mutate(cell_type2 = factor(cell_type2, levels = cell_type2))

ggplot(meta_count, aes(x = cell_type2, y = n, fill = cell_type2)) +
  geom_bar(stat = "identity", width = 0.9, color = "white", show.legend = FALSE) +
  labs(x = "Cell Types", y = "Cell Numbers") + coord_flip() +
  scale_fill_manual(values = cell_type_cols) +
  theme(axis.line = element_line(size = 0.5), panel.grid = element_blank(),
        panel.background = element_blank()) +
  scale_y_continuous(expand = c(0, 0))
ggsave(file.path(outdir2, "cell_number.pdf"), width = 7, height = 5)

# Sankey plot
df_prop <- as.data.frame(prop.table(table(Idents(merged), merged$disease), margin = 2) * 100)
ggplot(df_prop, aes(x = Var2, y = Freq, fill = Var1, stratum = Var1, alluvium = Var1)) +
  geom_col(width = 0.5, color = 'black') + geom_flow(width = 0.5, alpha = 0.4) +
  theme_bw() + labs(x = 'Disease', y = 'Fraction', fill = "Cell Type") +
  coord_flip() + scale_fill_manual(values = cell_type_cols) +
  theme(panel.grid = element_blank())
ggsave(file.path(outdir2, "cell_percentage.pdf"), width = 7, height = 5)