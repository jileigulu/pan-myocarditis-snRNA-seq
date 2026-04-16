# KEGG pathway enrichment analysis for neutrophil DEGs

library(readr)
library(biomaRt)
library(clusterProfiler)
library(org.Mm.eg.db)

setwd("~/pan-heart/results/Neutrophils2025.5.13")

outdir_data <- "./files/KEGG"
outdir_plots <- "./plots/KEGG"
dir.create(outdir_data, recursive = TRUE, showWarnings = FALSE)
dir.create(outdir_plots, recursive = TRUE, showWarnings = FALSE)

# Load DEG data
DEGs <- read_csv("files/deg/亚型层面/整体疾病组/deg1.csv")

celltypes <- unique(DEGs$celltype)
groups <- c("Upregulated", "Downregulated")

# Connect to Ensembl mouse database
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

for (cell in celltypes) {
  DEGs_sub <- DEGs %>% filter(celltype == cell)
  
  for (group in groups) {
    dir.create(file.path(outdir_data, group), showWarnings = FALSE, recursive = TRUE)
    
    group_degs <- DEGs_sub %>% filter(regulate == group)
    top_genes <- group_degs[order(abs(group_degs$avg_log2FC), decreasing = TRUE), ]$gene
    
    # Convert SYMBOL to ENTREZID
    convert <- getBM(attributes = c("external_gene_name", "entrezgene_id"),
                     filters = "external_gene_name", values = top_genes, mart = mart)
    colnames(convert) <- c("SYMBOL", "ENTREZID")
    
    entrez_genes <- convert$ENTREZID[!is.na(convert$ENTREZID)]
    
    if (length(entrez_genes) == 0) next
    
    # KEGG enrichment
    ego <- enrichKEGG(gene = entrez_genes, keyType = "kegg", organism = 'mmu',
                      pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2)
    
    if (is.null(ego) || nrow(ego@result) == 0) next
    
    # Convert ENTREZID to gene symbols
    ego <- setReadable(ego, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
    
    term <- ego@result
    term$celltype <- cell
    term$Description <- gsub("- Mus musculus (house mouse)", "", term$Description, fixed = TRUE)
    
    saveRDS(term, file.path(outdir_data, group, paste0(cell, "_term.rds")))
    write.csv(term, file.path(outdir_data, group, paste0(cell, "_term.csv")), row.names = FALSE)
  }
}