# Enrichment analysis for DEGs using clusterProfiler

## SETUP ----
library(clusterProfiler)
library(enrichplot)
library(org.Mm.eg.db)
library(tidyverse)

setwd("~/pan-heart/results/deg2025.3.31/deg")
dir.create("./files/cluster_profiler", recursive = TRUE, showWarnings = FALSE)
dir.create("./plots/cluster_profiler", recursive = TRUE, showWarnings = FALSE)

immune_colors <- c("NK cells" = "#f57c6e", "T cells" = "#f2b56f", "B cells" = "#b8aeeb",
                   "Dendritic cells" = "#B0DC66", "Macrophages&Monocytes" = "#88d8db",
                   "Neutrophils" = "#71b7ed")

## HELPER FUNCTION ----
run_enrichment <- function(deg_df, group_col, output_prefix) {
  # Convert gene symbols to ENTREZ IDs
  gene_ids <- bitr(deg_df$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
  deg_entrez <- merge(gene_ids, deg_df, by.x = "SYMBOL", by.y = "gene")
  
  # Separate up and down regulated genes
  up_genes <- deg_entrez[deg_entrez$regulate == "Upregulated", ]
  down_genes <- deg_entrez[deg_entrez$regulate == "Downregulated", ]
  
  # Run enrichment
  if (nrow(up_genes) > 0) {
    ck_up <- compareCluster(ENTREZID ~ get(group_col), data = up_genes,
                            fun = "enrichGO", OrgDb = org.Mm.eg.db, ont = "BP",
                            pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
    saveRDS(ck_up, file.path("./files/cluster_profiler", paste0(output_prefix, "_UP.rds")))
  }
  
  if (nrow(down_genes) > 0) {
    ck_down <- compareCluster(ENTREZID ~ get(group_col), data = down_genes,
                              fun = "enrichGO", OrgDb = org.Mm.eg.db, ont = "BP",
                              pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
    saveRDS(ck_down, file.path("./files/cluster_profiler", paste0(output_prefix, "_DOWN.rds")))
  }
}

## BY DISEASE ----
deg_eam <- read_csv("files/EAM_deg_TOTAL.csv") %>% mutate(disease = "EAM")
deg_ici <- read_csv("files/ICI-MC_deg_TOTAL.csv") %>% mutate(disease = "ICI-MC")
deg_vmc <- read_csv("files/VMC_deg_TOTAL.csv") %>% mutate(disease = "VMC")
deg_all <- bind_rows(deg_eam, deg_ici, deg_vmc)

run_enrichment(deg_all, "disease", "all_diseases")

## BY CELL TYPE WITHIN EACH DISEASE ----
for (disease in c("EAM", "ICI-MC", "VMC")) {
  deg_subset <- deg_all[deg_all$disease == disease, ]
  if (nrow(deg_subset) > 0) {
    run_enrichment(deg_subset, "celltype", disease)
  }
}

## PLOT ENRICHMENT RESULTS ----
plot_enrichment_dot <- function(rds_path, title, width = 10, height = 8) {
  if (file.exists(rds_path)) {
    ck <- readRDS(rds_path)
    p <- dotplot(ck, showCategory = 10, color = "qvalue", title = title) +
      theme(plot.title = element_text(hjust = 0.5))
    ggsave(file.path("./plots/cluster_profiler", paste0(basename(rds_path), ".pdf")),
           p, width = width, height = height)
  }
}

plot_enrichment_emap <- function(rds_path, title, width = 12, height = 6.5) {
  if (file.exists(rds_path)) {
    ck <- readRDS(rds_path)
    ck <- pairwise_termsim(ck)
    ck <- clusterProfiler::simplify(ck, cutoff = 0.7, by = "p.adjust", select_fun = min)
    
    pdf(file.path("./plots/cluster_profiler", paste0(title, "_emap.pdf")), width = width, height = height)
    p <- emapplot(ck, showCategory = 20, legend_n = 2,
                  cluster.params = list(cluster = TRUE, method = stats::kmeans,
                                        n = 5, legend = TRUE, label_style = "shadowtext",
                                        label_words_n = 4, label_format = 30)) +
      scale_fill_manual(values = immune_colors)
    print(p)
    dev.off()
  }
}

# Generate plots for all enrichment results
enrichment_files <- list.files("./files/cluster_profiler", pattern = "\\.rds$", full.names = TRUE)
for (f in enrichment_files) {
  plot_enrichment_dot(f, gsub("\\.rds$", "", basename(f)))
}