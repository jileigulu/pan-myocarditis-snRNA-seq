# Basic analysis for neutrophils: UMAP, markers, cell proportion

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggalluvial)

setwd("~/pan-heart/results/Neutrophils2025.5.13")
outdir <- "./files"
outdir2 <- "./plots"
dir.create(outdir, showWarnings = FALSE)
dir.create(outdir2, showWarnings = FALSE)

# Load data
merged <- readRDS(file.path(outdir, "merged.rds"))

# Colors
cell_type_cols <- c("Neutro_Ccl4" = "#87CEEB", "Neutro_Txnip" = "#FF6F61",
                    "Neutro_Ifit1" = "#A8E063", "Neutro_Ngp" = "#FFC0CB")

# 1. UMAP plot----

DimPlot(merged, label = TRUE, pt.size = 2, cols = cell_type_cols) +
  NoLegend() + labs(x = "UMAP1", y = "UMAP2", title = "Neutrophils") +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1),
        axis.text = element_blank(), axis.ticks = element_blank())
ggsave(file.path(outdir2, "dimplot.pdf"), width = 5, height = 5)

# 2. Marker genes dot plot----

DefaultAssay(merged) <- 'RNA'
markers <- FindAllMarkers(merged, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write_csv(markers, file.path(outdir, "markers.csv"))

genes <- c("Ccl4", "Ccl3", "Nfkbia", "Txnip", "Fos", "Rgs2", 
           "Ifit1", "Ifi47", "Gbp2", "Ngp", "Ifitm6", "Ly6g")

source("~/tutorial/custom_plot_function.R")
p <- DotPlot(merged, dot.scale = 4, col.min = -1, col.max = 1, features = genes) +
  scale_color_distiller(palette = "RdBu", direction = -1) +
  coord_fixed(ratio = 1) + theme_cat() +
  theme(axis.title = element_blank(), legend.margin = margin(b = -8),
        legend.position = "top", axis.text.x = element_text(face = "italic")) +
  guides(x = guide_axis(angle = 90), color = guide_colorbar(title = "Expression")) +
  scale_y_discrete(limits = rev)
ggsave(file.path(outdir2, "markers_dotplot.pdf"), width = 7, height = 4)

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

# 4. DEG petal plot----

deg <- read_csv("files/deg/deg1.csv")
df_count <- table(deg$disease, deg$regulate) %>% as.data.frame()
colnames(df_count) <- c("Sample", "Group", "Freq")

ggplot(df_count, aes(x = Sample, y = Freq, fill = Group)) +
  geom_col(width = 0.9, alpha = 0.8) + coord_polar() +
  labs(x = "Number of DEG") + theme_bw() +
  theme(aspect.ratio = 1, axis.title.y = element_blank(),
        axis.text = element_text(size = 8), legend.position = "top") +
  scale_fill_manual(values = c("Upregulated" = "#E41A1C", "Downregulated" = "#579CC7"))
ggsave(file.path(outdir2, "差异基因花瓣图.pdf"), width = 5, height = 5)