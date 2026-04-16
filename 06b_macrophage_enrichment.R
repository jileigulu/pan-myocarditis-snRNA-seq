# GO enrichment analysis for macrophage subtype DEGs
# Follows 06_macrophage_deg.R

## SETUP ----
library(clusterProfiler)
library(org.Mm.eg.db)
library(dplyr)
library(readr)

output_base <- "~/pan-heart/results/Macrophage2025.4.16/files/GO"

# 1. ENRICHMENT FOR THREE DISEASE COMPARISONS----

deg_3disease <- read_csv("~/pan-heart/results/Macrophage2025.4.16/files/deg/亚型层面/三种心肌炎组/deg/deg_with_regulation.csv")
diseases <- unique(deg_3disease$disease)

for (d in diseases) {
  # Upregulated genes
  up_dir <- file.path(output_base, "三种心肌炎组/Upregulated", d)
  dir.create(up_dir, showWarnings = FALSE, recursive = TRUE)
  
  deg_up <- deg_3disease[deg_3disease$disease == d & deg_3disease$regulate == "Upregulated", ]
  top_genes_up <- deg_up[order(deg_up$avg_log2FC, decreasing = TRUE), ]$gene
  
  if (length(top_genes_up) > 0) {
    bp_up <- enrichGO(top_genes_up, OrgDb = org.Mm.eg.db, keyType = 'SYMBOL',
                      ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2)
    term_up <- bp_up@result
    term_up$disease <- d
    write.table(term_up, file = file.path(up_dir, paste0(d, ".csv")),
                row.names = FALSE, quote = FALSE, sep = ",")
  }
  
  # Downregulated genes
  down_dir <- file.path(output_base, "三种心肌炎组/Downregulated", d)
  dir.create(down_dir, showWarnings = FALSE, recursive = TRUE)
  
  deg_down <- deg_3disease[deg_3disease$disease == d & deg_3disease$regulate == "Downregulated", ]
  top_genes_down <- deg_down[order(deg_down$avg_log2FC, decreasing = FALSE), ]$gene
  
  if (length(top_genes_down) > 0) {
    bp_down <- enrichGO(top_genes_down, OrgDb = org.Mm.eg.db, keyType = 'SYMBOL',
                        ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2)
    term_down <- bp_down@result
    term_down$disease <- d
    write.table(term_down, file = file.path(down_dir, paste0(d, ".csv")),
                row.names = FALSE, quote = FALSE, sep = ",")
  }
}


# 2. ENRICHMENT FOR COMBINED MYOCARDITIS COMPARISON----

deg_combined <- read_csv("~/pan-heart/results/Macrophage2025.4.16/files/deg/亚型层面/整体疾病组/deg/deg_with_regulation.csv")
celltypes <- unique(deg_combined$celltype)

for (ct in celltypes) {
  # Upregulated genes
  up_dir <- file.path(output_base, "整体疾病组/Upregulated", ct)
  dir.create(up_dir, showWarnings = FALSE, recursive = TRUE)
  
  deg_up <- deg_combined[deg_combined$celltype == ct & deg_combined$regulate == "Upregulated", ]
  top_genes_up <- deg_up[order(deg_up$avg_log2FC, decreasing = TRUE), ]$gene
  
  if (length(top_genes_up) > 0) {
    bp_up <- enrichGO(top_genes_up, OrgDb = org.Mm.eg.db, keyType = 'SYMBOL',
                      ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2)
    term_up <- bp_up@result
    term_up$celltype <- ct
    write.table(term_up, file = file.path(up_dir, paste0(ct, ".csv")),
                row.names = FALSE, quote = FALSE, sep = ",")
  }
  
  # Downregulated genes
  down_dir <- file.path(output_base, "整体疾病组/Downregulated", ct)
  dir.create(down_dir, showWarnings = FALSE, recursive = TRUE)
  
  deg_down <- deg_combined[deg_combined$celltype == ct & deg_combined$regulate == "Downregulated", ]
  top_genes_down <- deg_down[order(deg_down$avg_log2FC, decreasing = FALSE), ]$gene
  
  if (length(top_genes_down) > 0) {
    bp_down <- enrichGO(top_genes_down, OrgDb = org.Mm.eg.db, keyType = 'SYMBOL',
                        ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2)
    term_down <- bp_down@result
    term_down$celltype <- ct
    write.table(term_down, file = file.path(down_dir, paste0(ct, ".csv")),
                row.names = FALSE, quote = FALSE, sep = ",")
  }
}

