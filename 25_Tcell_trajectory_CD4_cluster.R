# Cluster gene enrichment analysis for trajectory branches

library(clusterProfiler)
library(org.Mm.eg.db)
library(dplyr)
library(ggplot2)

setwd('~/pan-heart/results/T_cells2025.5.20')
outdir <- './files/trajectory/monocle2/CD4'
outdir2 <- './plots/trajectory/monocle2/CD4'

# Load cluster genes
clusters <- readRDS(file.path(outdir, 'cluster_genes.rds'))

# GO enrichment for each cluster
for (i in unique(clusters$cluster)) {
  genes <- clusters[clusters$cluster == i, ]$genes
  
  bp <- enrichGO(genes, OrgDb = org.Mm.eg.db, keyType = 'SYMBOL',
                 ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05)
  
  if (!is.null(bp) && nrow(bp@result) > 0) {
    write.csv(bp@result, file.path(outdir, paste0("term_", i, "_gobp.csv")), row.names = FALSE)
    
    # Plot selected pathways for each cluster
    cluster_pathways <- list(
      `1` = c('leukocyte proliferation', 'regulation of T cell activation', 'regulation of leukocyte cell-cell adhesion'),
      `2` = c('actin filament organization', 'antigen receptor-mediated signaling pathway', 'regulation of interleukin-2 production'),
      `3` = c('leukocyte cell-cell adhesion', 'cytokine-mediated signaling pathway', 'adaptive immune response'),
      `4` = c('defense response to virus', 'activation of innate immune response', 'response to interferon-beta'),
      `5` = c('ribosome biogenesis', 'translation at synapse', 'signal transduction by p53 class mediator')
    )
    
    if (as.character(i) %in% names(cluster_pathways)) {
      df <- bp@result[bp@result$Description %in% cluster_pathways[[as.character(i)]], ]
      if (nrow(df) > 0) {
        df$labelx <- 0
        df$labely <- seq(nrow(df), 1)
        
        p <- ggplot(df, aes(x = -log10(pvalue), y = reorder(Description, -log10(pvalue)))) +
          geom_bar(stat = "identity", alpha = 1, fill = "#00D983", width = 0.8) +
          geom_text(aes(x = labelx, y = labely, label = Description), size = 3.5, hjust = 0) +
          theme_classic() + xlab("-log10(pvalue)") + ggtitle(paste0('cluster', i)) +
          theme(axis.text.y = element_blank(), axis.line.y = element_blank(),
                axis.title.y = element_blank(), axis.ticks.y = element_blank())
        
        ggsave(file.path(outdir2, paste0('cluster', i, '_Barplot.pdf')), p, height = 2, width = 4)
      }
    }
  }
}