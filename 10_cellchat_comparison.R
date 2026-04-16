# Comparative analysis of CellChat results across disease conditions

library(CellChat)
library(patchwork)
library(ComplexHeatmap)
library(RColorBrewer)

setwd("~/pan-heart/results/Macrophage2025.4.16/plots/cellchat")

# Load CellChat objects
cellchat.control <- readRDS("~/pan-heart/results/Macrophage2025.4.16/files/cellchat/Normal/cellchat.rds")
cellchat.case1 <- readRDS("~/pan-heart/results/Macrophage2025.4.16/files/cellchat/EAM/cellchat.rds")
cellchat.case2 <- readRDS("~/pan-heart/results/Macrophage2025.4.16/files/cellchat/ICI-MC/cellchat.rds")
cellchat.case3 <- readRDS("~/pan-heart/results/Macrophage2025.4.16/files/cellchat/VMC/cellchat.rds")

# Merge CellChat objects
object.list <- list(Normal = cellchat.control, 
                    EAM = cellchat.case1,
                    "ICI-MC" = cellchat.case2,
                    VMC = cellchat.case3)

for (i in 1:length(object.list)) {
  object.list[[i]] <- netAnalysis_computeCentrality(object.list[[i]])
}

cellchat <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)

# Color settings
disease_colors <- c("Normal" = "#BCD5F1", "EAM" = "#f57c6e",
                    "ICI-MC" = "#B0DC66", "VMC" = "#EDCCEE")

celltype_levels <- levels(object.list[[1]]@idents)
color_palette <- c("#FFA500", "#89CB6C", "#FF6A6A", "#A6CEE3", 
                   "#E41A1C", "#F781BF", "#6495ED", "#AB82FF")
names(color_palette) <- celltype_levels

# 1. Bar plot: comparison of interactions----

gg1 <- compareInteractions(cellchat, show.legend = FALSE, group = c(1,2,3,4), 
                           color.use = unname(disease_colors))
gg2 <- compareInteractions(cellchat, show.legend = FALSE, group = c(1,2,3,4), 
                           measure = "weight", color.use = unname(disease_colors))
gg1 + gg2
ggsave("barplot.pdf", width = 8, height = 5)

# 2. Scatter plot: outgoing/incoming signaling strength----

num.link <- sapply(object.list, function(x) {
  rowSums(x@net$count) + colSums(x@net$count) - diag(x@net$count)
})
weight.MinMax <- c(min(num.link), max(num.link))

gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], 
                                               title = names(object.list)[i], 
                                               weight.MinMax = weight.MinMax,
                                               color.use = color_palette)
}
patchwork::wrap_plots(plots = gg)
ggsave("income_outgoing_strength.pdf", height = 6, width = 8)

# 3. Signaling changes for specific cell types----

# Macro_TLF+
gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Macro_TLF+", 
                                            comparison = c(1, 2), signaling.exclude = "CD45")
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Macro_TLF+", 
                                            comparison = c(1, 3), signaling.exclude = "CD45")
gg3 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Macro_TLF+", 
                                            comparison = c(1, 4), signaling.exclude = "CD45")
gg1 / gg2 / gg3
ggsave("信号变化Macro_TLF.pdf", height = 14, width = 7, scale = 1)

# CD8+T
gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "CD8+T", 
                                            comparison = c(1, 2), signaling.exclude = "CD45")
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "CD8+T", 
                                            comparison = c(1, 3), signaling.exclude = "CD45")
gg3 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "CD8+T", 
                                            comparison = c(1, 4), signaling.exclude = "CD45")
gg1 / gg2 / gg3
ggsave("信号变化CD8T.pdf", height = 14, width = 7, scale = 1)

# 4. Stacked bar plot: signaling pathways ranking----

signaling_pathways <- c("CXCL", "MHC-II", "MHC-I", "PD-L1", "TNF", "ICOS", "CD40",
                        "CCL", "CD80", "CD86", "IFN-II", "TGFb", "COMPLEMENT",
                        "ITGAL-ITGB2", "SEMA4", "GALECTIN")

gg1 <- rankNet(cellchat, mode = "comparison", stacked = TRUE, 
               comparison = c(1,2,3,4), signaling = signaling_pathways,
               do.stat = TRUE, color.use = disease_colors) + coord_flip()
gg2 <- rankNet(cellchat, mode = "comparison", stacked = FALSE, 
               comparison = c(1,2,3,4), signaling = signaling_pathways,
               do.stat = TRUE, color.use = disease_colors) + coord_flip()
gg1 + gg2
ggsave("堆叠柱形图.pdf", width = 14, height = 7)

# 5. Heatmap: outgoing/incoming signaling patterns----

pathway_s <- c('CCL', 'CXCL', 'TNF', 'MHC-II', 'ICAM', 'ITGAL-ITGB2', 'PD-L1', 
               'CD40', 'CD80', 'CD86', 'IFN-II', 'COMPLEMENT', 'FASLG', 'TGFb', 
               'ICOS', 'SEMA4', 'TIGIT', 'GALECTIN', 'THBS', 'SPP1')

# Outgoing patterns
ht1 <- netAnalysis_signalingRole_heatmap(object.list[[2]], pattern = "outgoing", 
                                         signaling = pathway_s, title = "EAM", 
                                         width = 5, height = 10, color.heatmap = "BuPu")
ht2 <- netAnalysis_signalingRole_heatmap(object.list[[3]], pattern = "outgoing", 
                                         signaling = pathway_s, title = "ICI-MC", 
                                         width = 5, height = 10, color.heatmap = "BuPu")
ht3 <- netAnalysis_signalingRole_heatmap(object.list[[4]], pattern = "outgoing", 
                                         signaling = pathway_s, title = "VMC", 
                                         width = 5, height = 10, color.heatmap = "BuPu")
pdf("outgoing_signaling_patterns_heatmap2.pdf", width = 15, height = 12)
draw(ht1 + ht2 + ht3, ht_gap = unit(0.5, "cm"))
dev.off()

# Incoming patterns
ht1 <- netAnalysis_signalingRole_heatmap(object.list[[2]], pattern = "incoming", 
                                         signaling = pathway_s, title = "EAM", 
                                         width = 5, height = 10, color.heatmap = "GnBu")
ht2 <- netAnalysis_signalingRole_heatmap(object.list[[3]], pattern = "incoming", 
                                         signaling = pathway_s, title = "ICI-MC", 
                                         width = 5, height = 10, color.heatmap = "GnBu")
ht3 <- netAnalysis_signalingRole_heatmap(object.list[[4]], pattern = "incoming", 
                                         signaling = pathway_s, title = "VMC", 
                                         width = 5, height = 10, color.heatmap = "GnBu")
pdf("incoming_signaling_patterns_heatmap2.pdf", width = 15, height = 12)
draw(ht1 + ht2 + ht3, ht_gap = unit(0.5, "cm"))
dev.off()

# 6. Differential interaction heatmaps----

# EAM vs Normal
gg1 <- netVisual_heatmap(cellchat, comparison = c(1, 2), 
                         title.name = "Differential number of interactions (Normal vs. EAM)")
gg2 <- netVisual_heatmap(cellchat, comparison = c(1, 2), measure = "weight",
                         title.name = "Differential interaction strength (Normal vs. EAM)")
pdf("差异互作强度热图_EAM.pdf", width = 10, height = 5)
print(gg1 + gg2)
dev.off()

# ICI-MC vs Normal
gg1 <- netVisual_heatmap(cellchat, comparison = c(1, 3),
                         title.name = "Differential number of interactions (Normal vs. ICI-MC)")
gg2 <- netVisual_heatmap(cellchat, comparison = c(1, 3), measure = "weight",
                         title.name = "Differential interaction strength (Normal vs. ICI-MC)")
pdf("差异互作强度热图_ICI-MC.pdf", width = 10, height = 5)
print(gg1 + gg2)
dev.off()

# VMC vs Normal
gg1 <- netVisual_heatmap(cellchat, comparison = c(1, 4),
                         title.name = "Differential number of interactions (Normal vs. VMC)")
gg2 <- netVisual_heatmap(cellchat, comparison = c(1, 4), measure = "weight",
                         title.name = "Differential interaction strength (Normal vs. VMC)")
pdf("差异互作强度热图_VMC.pdf", width = 10, height = 5)
print(gg1 + gg2)
dev.off()

# 7. Circle plot: interaction strength for all groups----

pdf("各组通讯强度网络图.pdf", height = 5, width = 10)
par(mfrow = c(1, 4), xpd = TRUE)
for (i in 1:length(object.list)) {
  groupSize <- as.numeric(table(cellchat@idents[[i]]))
  netVisual_circle(object.list[[i]]@net$weight, vertex.weight = groupSize,
                   weight.scale = TRUE, label.edge = FALSE,
                   title.name = paste0("Interaction strength - ", names(object.list)[i]),
                   color.use = color_palette)
}
dev.off()

# 8. Bubble plot: specific signaling pathways----

# CCL pathway
netVisual_bubble(cellchat, sources.use = c(5:8), targets.use = c(2),
                 signaling = c("CCL"), comparison = c(1:4),
                 angle.x = 45, remove.isolate = TRUE)
ggsave("气泡图Ccl.pdf", height = 4, width = 6)

# 9. Network plot: MHC-I pathway----

pathways.show <- c("MHC-I")
object.list1 <- object.list[c(1,2,3)]
weight.max <- getMaxWeight(object.list1, slot.name = 'netP', attribute = pathways.show)

pdf(paste0(pathways.show, '.a.pdf'), height = 6, width = 12)
par(mfrow = c(1, 4), xpd = TRUE)
for (i in 1:length(object.list1)) {
  netVisual_aggregate(object.list1[[i]], signaling = pathways.show, layout = "circle",
                      edge.weight.max = weight.max[1], edge.width.max = 10,
                      signaling.name = paste(pathways.show, names(object.list1)[i]),
                      color.use = color_palette)
}
dev.off()

# CCL pathway (all groups)
pathways.show <- c("CCL")
object.list1 <- object.list[c(1,2,3,4)]
weight.max <- getMaxWeight(object.list1, slot.name = 'netP', attribute = pathways.show)

pdf(paste0(pathways.show, '.a.pdf'), height = 6, width = 16)
par(mfrow = c(1, 4), xpd = TRUE)
for (i in 1:length(object.list1)) {
  netVisual_aggregate(object.list1[[i]], signaling = pathways.show, layout = "circle",
                      edge.weight.max = weight.max[1], edge.width.max = 10,
                      signaling.name = paste(pathways.show, names(object.list1)[i]),
                      color.use = color_palette)
}
dev.off()

# 10. Chord diagram: signaling pathways----

# MHC-I chord
pathways.show <- c("MHC-1")
weight.max <- getMaxWeight(object.list, slot.name = 'netP', attribute = pathways.show)

pdf(paste0(pathways.show, '_chord.pdf'), height = 5, width = 12)
par(mfrow = c(1, 4), xpd = TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord",
                      edge.weight.max = weight.max[1], edge.width.max = 10,
                      signaling.name = paste(pathways.show, "of", names(object.list)[i]))
}
dev.off()

# MHC-II chord (EAM only)
pathways.show <- c("MHC-II")
object.list1 <- object.list[c(2)]
weight.max <- getMaxWeight(object.list1, slot.name = 'netP', attribute = pathways.show)

pdf(paste0(pathways.show, '_chord.pdf'), height = 5, width = 12)
par(mfrow = c(1, 4), xpd = TRUE)
for (i in 1:length(object.list1)) {
  netVisual_aggregate(object.list1[[i]], signaling = pathways.show, layout = "chord",
                      edge.weight.max = weight.max[1], edge.width.max = 10,
                      signaling.name = paste(pathways.show, "of", names(object.list1)[i]))
}
dev.off()

# 11. Chord gene plot: NRXN pathway----

pathways.show <- c("NRXN")

# NRXN as sender (from Oligo to all cells)
pdf("NRXN.Oligo_as_sender.pdf", height = 20, width = 15)
par(mfrow = c(1, 2), xpd = TRUE)
for (i in 1:length(object.list)) {
  netVisual_chord_gene(object.list[[i]], sources.use = c(3), targets.use = c(1:8),
                       lab.cex = 0.5, title.name = paste0("Signaling from Oligo to cells - ", 
                                                          names(object.list)[i]))
}
dev.off()

# NRXN as receptor (from all cells to Oligo)
pdf("NRXN.Oligo_as_receptor.pdf", height = 20, width = 15)
par(mfrow = c(1, 2), xpd = TRUE)
for (i in 1:length(object.list)) {
  netVisual_chord_gene(object.list[[i]], sources.use = c(1:8), targets.use = c(3),
                       lab.cex = 0.5, title.name = paste0("Signaling from cells to Oligo - ", 
                                                          names(object.list)[i]))
}
dev.off()

# 12. Violin plot: gene expression comparison----

cellchat@meta$datasets <- factor(cellchat@meta$datasets, levels = c("CON", "ASD"))
pdf("NRXN.violin.pdf", height = 5, width = 7)
plotGeneExpression(cellchat, signaling = "NRXN", split.by = "cell_type", 
                   colors.ggplot = TRUE) +
  scale_fill_manual(values = c('#00AFBB', "#e64b35ff"))
dev.off()