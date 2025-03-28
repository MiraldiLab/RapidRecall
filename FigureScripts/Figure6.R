library(reshape2)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(matrixStats)
library(dplyr)
library(dendextend)
library(matrixStats)
library(circlize)
library(RColorBrewer)
library(Nebulosa)
library(patchwork)
library(Seurat)
library(Signac)
library(viridis)

outdir <- "output/Figure6"
indir <- "input"

# Sort out Naive and Naive_Responsive cells
obj <- readRDS(paste0(indir, "/integrated2.rds"))
Naive_obj <- readRDS(paste0(indir, "/Naive_obj.rds"))
labels <- read.table(paste0(indir, "/NaiveExpanded_labels.txt"))$V1
index <- which(colnames(obj) %in% colnames(Naive_obj))
labels <- labels[index]
Naive_obj$celltype2 <- labels
Naive_obj$celltype2[which(Naive_obj$celltype2 == "Naive_rest2")] <- "Naive_rest"
scores <- read.table(paste0(indir, "/scores_combined.txt"))
index <- which(colnames(obj) %in% colnames(Naive_obj))
scores <- scores[index, ]

######### Fig6 A,B
pdf(paste0(outdir, "/wnn_umap_condition.pdf"), width = 6, height = 6)
DimPlot(Naive_obj, reduction = "wnn.umap", group.by = "Condition")
dev.off()
pdf(paste0(outdir, "/wnn_umap_celltype.pdf"), width = 6, height = 6)
DimPlot(Naive_obj, reduction = "wnn.umap", group.by = "celltype2")
dev.off()
pdf(paste0(outdir,"/wnn_umap_pseudotime.pdf"), width = 4, height = 4)
FeaturePlot(Naive_obj, reduction = "wnn.umap", features = "Pseudotime")
dev.off()

########### Fig6 E: Singlecell heatmap of pseudotime, module scores, TF expressions
module_order <- order(scores$Memory, decreasing = F)
cell_order <- colnames(Naive_obj)[module_order]
counts <- GetAssayData(Naive_obj, slot = "counts")[, cell_order]
scaling_factor <- colSums(counts) / 1e6
PRDM1_counts <- log(counts["PRDM1", ] / scaling_factor + 1)
YBX1_counts <- log(counts["YBX1", ] / scaling_factor + 1)
KLF6_counts <- log(counts["KLF6", ] / scaling_factor + 1)
RUNX3_counts <- log(counts["RUNX3", ] / scaling_factor + 1)
RUNX2_counts <- log(counts["RUNX2", ] / scaling_factor + 1)
LEF1_counts <- log(counts["LEF1", ] / scaling_factor + 1)
TET2_counts <- log(counts["TET2", ] / scaling_factor + 1)
Memory_up <- scores[cell_order, "Memory"]
Naive_up <- scores[cell_order, "Naive"]

data_mat <- rbind(Naive_obj$Pseudotime[cell_order], Memory_up, Naive_up, KLF6_counts, TET2_counts, RUNX3_counts, RUNX2_counts, LEF1_counts)
data_mat <- as.matrix(data_mat)
rownames(data_mat) <- c("Pseudotime","MemoryUp","MemoryDown","KLF6","TET2","RUNX3","RUNX2","LEF1")
breaks <- seq(min(data_mat[1,], na.rm = TRUE), max(data_mat[1,], na.rm = TRUE), length.out = 100)
colors <- viridis(100, option = "viridis")
pseudotime_col_fun <- colorRamp2(breaks, colors)
module_col_fun <- colorRamp2(c(-0.05, 0, 0.1), c('dodgerblue3', 'white', 'red'))
gene_expr_col_fun <- colorRamp2(c(8,7,6,5,4,3,2,1), brewer.pal(8,"Spectral"))

# Create individual heatmaps
ht1 <- Heatmap(matrix(data_mat[1, ], nrow = 1), col = pseudotime_col_fun, show_row_names = TRUE, cluster_rows = FALSE, cluster_columns = FALSE, show_column_names = FALSE, name = "Pseudotime")
ht2 <- Heatmap(matrix(data_mat[2, ], nrow = 1), col = module_col_fun, show_row_names = TRUE, cluster_rows = FALSE, cluster_columns = FALSE, show_column_names = FALSE, name = "Memory_up")
ht3 <- Heatmap(matrix(data_mat[3, ], nrow = 1), col = module_col_fun, show_row_names = TRUE, cluster_rows = FALSE, cluster_columns = FALSE, show_column_names = FALSE, name = "Memory_down")
ht4 <- Heatmap(matrix(data_mat[4, ], nrow = 1), col = gene_expr_col_fun, show_row_names = TRUE, cluster_rows = FALSE, cluster_columns = FALSE, show_column_names = FALSE, name = "KLF6_counts")
ht5 <- Heatmap(matrix(data_mat[5, ], nrow = 1), col = gene_expr_col_fun, show_row_names = TRUE, cluster_rows = FALSE, cluster_columns = FALSE, show_column_names = FALSE, name = "TET2_counts")
ht6 <- Heatmap(matrix(data_mat[6, ], nrow = 1), col = gene_expr_col_fun, show_row_names = TRUE, cluster_rows = FALSE, cluster_columns = FALSE, show_column_names = FALSE, name = "RUNX3_counts")
ht7 <- Heatmap(matrix(data_mat[7, ], nrow = 1), col = gene_expr_col_fun, show_row_names = TRUE, cluster_rows = FALSE, cluster_columns = FALSE, show_column_names = FALSE, name = "RUNX2_counts")
ht8 <- Heatmap(matrix(data_mat[8, ], nrow = 1), col = gene_expr_col_fun, show_row_names = TRUE, cluster_rows = FALSE, cluster_columns = FALSE, show_column_names = FALSE, name = "LEF1_counts")

# Draw each heatmap and capture as a grob
g1 <- grid.grabExpr(draw(ht1))
g2 <- grid.grabExpr(draw(ht2))
g3 <- grid.grabExpr(draw(ht3))
g4 <- grid.grabExpr(draw(ht4))
g5 <- grid.grabExpr(draw(ht5))
g6 <- grid.grabExpr(draw(ht6))
g7 <- grid.grabExpr(draw(ht7))
g8 <- grid.grabExpr(draw(ht8))
dev.off()

# Arrange the heatmaps horizontally
pdf(paste0(outdir, "/heatmap_PART1.pdf"), width = 8, height = 10)
grid.arrange(g1, g2, g3, g4, g5, g6, g7, g8, nrow = 8)
dev.off()


## Load in resting peaks sets
DefaultAssay(Naive_obj) <- "peaks"
scores <- read.table(paste0(indir, "/scores_Naive.txt"))
colnames(scores) <- c("Rest1","Rest2","Rest3","Rest4","Rest5","Rest6","RR1","RR2","RR3","RR4","RR5","RR6","RR7","Act1","Act2","Act3","Act4","Act5")
Naive_obj <- AddMetaData(Naive_obj, scores)

keep <- which(Naive_obj$Condition == "15" & Naive_obj$celltype2 %in% c("Naive_act","Naive_Responsive"))
Naive_obj2 <- subset(Naive_obj, cells = keep)


######### Figure S6A-C: Violin plots and statistics of module score comparisons in naive cells
Scores <- as.vector(c(Naive_obj2$Rest1,Naive_obj2$Rest2,Naive_obj2$Rest3,Naive_obj2$Rest4,Naive_obj2$Rest5,Naive_obj2$Rest6,
    Naive_obj2$RR1,Naive_obj2$RR2,Naive_obj2$RR3,
    Naive_obj2$RR4,Naive_obj2$RR5,Naive_obj2$RR6,Naive_obj2$RR7,
    Naive_obj2$Act1, Naive_obj2$Act2, Naive_obj2$Act3, Naive_obj2$Act4, Naive_obj2$Act5))
celltype <- rep(Naive_obj2$celltype2, 18)
donor <- rep(Naive_obj2$Donor, 18)
set <- c(rep("Naive_Treg", length(colnames(Naive_obj2))), rep("Naive", length(colnames(Naive_obj2))), rep("Th1_Th17_CTL", length(colnames(Naive_obj2))), 
    rep("Th2_MHCII", length(colnames(Naive_obj2))), rep("TEM_TCM", length(colnames(Naive_obj2))), rep("Shared", length(colnames(Naive_obj2))),
    rep("RR1", length(colnames(Naive_obj2))),rep("RR2", length(colnames(Naive_obj2))),rep("RR3", length(colnames(Naive_obj2))),rep("RR4", length(colnames(Naive_obj2))),
    rep("RR5", length(colnames(Naive_obj2))), rep("RR6", length(colnames(Naive_obj2))), rep("RR7", length(colnames(Naive_obj2))),
    rep("Act1", length(colnames(Naive_obj2))), rep("Act2", length(colnames(Naive_obj2))), rep("Act3", length(colnames(Naive_obj2))), 
    rep("Act4", length(colnames(Naive_obj2))),rep("Act5", length(colnames(Naive_obj2))))
df <- data.frame(Scores = Scores, Celltype = celltype, Set = set, Donor = donor)
df$Set <- factor(df$Set, levels = c("Naive_Treg","Naive","Th1_Th17_CTL","Th2_MHCII","TEM_TCM","Shared",
    "RR1","RR2","RR3","RR4","RR5","RR6","RR7", "Act1","Act2","Act3","Act4","Act5"))


df_rest <- df[which(df$Set %in% c("Naive_Treg","Naive","Th1_Th17_CTL","Th2_MHCII","TEM_TCM","Shared")),]
df_RR <- df[which(df$Set %in% c("RR1","RR2","RR3","RR4","RR5","RR6","RR7")),]
df_act <- df[which(df$Set %in% c("Act1","Act2","Act3","Act4","Act5")),]

labels <- c("Naive_Treg","Naive","Th1_Th17_CTL","Th2_MHCII","TEM_TCM","Shared", "Early","Middle","MiddleTEM","LateStrong","LateWeak","RapidClose","RapidCloseTEM",
    "ActOpen_TCM_Naive", "ActOpen","ActLateOpen","ActClose","ActCloseNaive")
pdf(paste0(outdir, "/Naive_Responsive_vln.pdf"), width = 15, height = 5)
ggplot(df, aes(x = Set, y = Scores, fill = Celltype)) + 
  geom_violin(position = position_dodge(width = 0.9)) + 
  geom_boxplot(aes(x = Set, fill = Celltype),width = 0.1, position = position_dodge(width = 0.9), outlier.shape = NA) + 
  scale_fill_manual(values = c("#E11A1A","#7F170E")) +
  scale_x_discrete(labels = labels) +  
  theme(axis.text.x = element_text(angle = 90, hjust = 1,color = "black"),panel.background = element_rect(fill = "white"), plot.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "gray"),  panel.grid.minor = element_line(color = "gray")) +
  labs(y = "Module Score", x = "PeakSet")
dev.off()

labels <- c("Naive_Treg","Naive","Th1_Th17_CTL","Th2_MHCII","TEM_TCM","Shared")
pdf(paste0(outdir, "/Naive_Responsive_vln_Rest.pdf"), width = 10, height = 5)
ggplot(df_rest, aes(x = Set, y = Scores, fill = Celltype)) + 
  geom_violin(position = position_dodge(width = 0.9)) + 
  geom_boxplot(aes(x = Set, fill = Celltype),width = 0.1, position = position_dodge(width = 0.9), outlier.shape = NA) + 
  scale_fill_manual(values = c("#E11A1A","#7F170E")) +
  scale_x_discrete(labels = labels) +  
  theme(axis.text.x = element_text(angle = 90, hjust = 1,color = "black"),panel.background = element_rect(fill = "white"), plot.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "gray"),  panel.grid.minor = element_line(color = "gray")) +
  labs(y = "Module Score", x = "PeakSet")
dev.off()

labels <- c("Early","Middle","MiddleTEM","LateStrong","LateWeak","RapidClose","RapidCloseTEM")
pdf(paste0(outdir, "/Naive_Responsive_vln_RR.pdf"), width = 10, height = 5)
ggplot(df_RR, aes(x = Set, y = Scores, fill = Celltype)) + 
  geom_violin(position = position_dodge(width = 0.9)) + 
  geom_boxplot(aes(x = Set, fill = Celltype),width = 0.1, position = position_dodge(width = 0.9), outlier.shape = NA) + 
  scale_fill_manual(values = c("#E11A1A","#7F170E")) +
  scale_x_discrete(labels = labels) +  
  theme(axis.text.x = element_text(angle = 90, hjust = 1,color = "black"),panel.background = element_rect(fill = "white"), plot.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "gray"),  panel.grid.minor = element_line(color = "gray")) + 
  labs(y = "Module Score", x = "PeakSet")
dev.off()

labels <- c("ActOpen_TCM_Naive", "ActOpen","ActLateOpen","ActClose","ActCloseNaive")
pdf(paste0(outdir, "/Naive_Responsive_vln_act.pdf"), width = 10, height = 5)
ggplot(df_act, aes(x = Set, y = Scores, fill = Celltype)) + 
  geom_violin(position = position_dodge(width = 0.9)) + 
  geom_boxplot(aes(x = Set, fill = Celltype),width = 0.1, position = position_dodge(width = 0.9), outlier.shape = NA) + 
  scale_fill_manual(values = c("#E11A1A","#7F170E")) +
  scale_x_discrete(labels = labels) +  
  theme(axis.text.x = element_text(angle = 90, hjust = 1,color = "black"),panel.background = element_rect(fill = "white"), plot.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "gray"),  panel.grid.minor = element_line(color = "gray")) + 
  labs(y = "Module Score", x = "PeakSet")
dev.off()


pval <- c()
size <- c()
for(i in unique(df$Set)){
    #group1 <- df$Scores[which(df$Set == i & df$Celltype == "Naive_Responsive")]
    #group2 <- df$Scores[which(df$Set == i & df$Celltype == "Naive_act")]
    group1 <- c(mean(df$Scores[which(df$Set == i & df$Celltype == "Naive_Responsive" & df$Donor == "Donor1")]), 
                mean(df$Scores[which(df$Set == i & df$Celltype == "Naive_Responsive" & df$Donor == "Donor2")]),
                mean(df$Scores[which(df$Set == i & df$Celltype == "Naive_Responsive" & df$Donor == "Donor3")]),
                mean(df$Scores[which(df$Set == i & df$Celltype == "Naive_Responsive" & df$Donor == "Donor4")]))
    group2 <- c(mean(df$Scores[which(df$Set == i & df$Celltype == "Naive_act" & df$Donor == "Donor1")]), 
                mean(df$Scores[which(df$Set == i & df$Celltype == "Naive_act" & df$Donor == "Donor2")]),
                mean(df$Scores[which(df$Set == i & df$Celltype == "Naive_act" & df$Donor == "Donor3")]),
                mean(df$Scores[which(df$Set == i & df$Celltype == "Naive_act" & df$Donor == "Donor4")]))            
    pval <- c(pval, t.test(x = group1, y = group2)$p.val)
    mean_diff <- mean(group1) - mean(group2)
    pooled_sd <- sqrt(((length(group1) - 1) * var(group1) + (length(group2) - 1) * var(group2)) / 
                  (length(group1) + length(group2) - 2))
    cohens_d <- mean_diff / pooled_sd
    size <- c(size, cohens_d)
}
pval <- p.adjust(pval, method = "BH")

df <- data.frame(Pval = -log10(pval), Size = size, Set = c("Naive_Treg","Naive","Th1_Th17_CTL","Th2_MHCII","TEM_TCM","Shared","RR1","RR2","RR3","RR4","RR5","RR6","RR7","Act1","Act2","Act3","Act4","Act5"))
df$Set <- factor(df$Set, levels = c("Naive_Treg","Naive","Th1_Th17_CTL","Th2_MHCII","TEM_TCM","Shared","RR1","RR2","RR3","RR4","RR5","RR6","RR7","Act1","Act2","Act3","Act4","Act5"))


labels <- c("Naive_Treg","Naive","Th1_Th17_CTL","Th2_MHCII","TEM_TCM","Shared", "Early","Middle","MiddleTEM","LateStrong","LateWeak","RapidClose","RapidCloseTEM","ActOpen_TCM_Naive", "ActOpen","ActLateOpen","ActClose","ActCloseNaive")
pdf(paste0(outdir, "/Pvals_bar.pdf"), width = 10, height = 5)
ggplot(df, aes(x = Set, y = Pval, fill = Set)) + geom_bar(stat = "identity") +
    scale_fill_manual(values = c("#3F572A","#3F572A","#5E813F","#5E813F","#5E813F","#5E813F","#385492","#385492","#385492","#385492","#385492","#385492","#385492","#AF2318","#AF2318","#AF2318","#AF2318","#AF2318")) +
    scale_x_discrete(labels = labels) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1,color = "black"),panel.background = element_rect(fill = "white"), plot.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "gray"),  panel.grid.minor = element_line(color = "gray")) +
    labs(y = "-Log10(p)", x = "PeakSet")
dev.off()
pdf(paste0(outdir, "/EffectSize_bar.pdf"), width = 10, height = 5)
ggplot(df, aes(x = Set, y = Size, fill = Set)) + geom_bar(stat = "identity") +
    scale_fill_manual(values = c("#3F572A","#3F572A","#5E813F","#5E813F","#5E813F","#5E813F","#385492","#385492","#385492","#385492","#385492","#385492","#385492","#AF2318","#AF2318","#AF2318","#AF2318","#AF2318")) +
    scale_x_discrete(labels = labels) +  
  theme(axis.text.x = element_text(angle = 90, hjust = 1,color = "black"),panel.background = element_rect(fill = "white"), plot.background = element_rect(fill = "white"), panel.grid.major = element_line(color = "gray"),  panel.grid.minor = element_line(color = "gray")) +
    labs(y = "Cohens-d (Effect Size)", x = "PeakSet")
dev.off()

## Plot Density of Resting memory population module scores on umap
scores <- read.table(paste0(indir, "/scores_all.txt"))
colnames(scores) <- c("Rest1","Rest2","Rest3","Rest4","Rest5","Rest6","RR1","RR2","RR3","RR4","RR5","RR6","RR7","Act1","Act2","Act3","Act4","Act")
obj <- AddMetaData(obj, scores)

index <- which(obj$celltype2 == "Th1_rest" & obj$Condition == "Resting")
Th1 <- subset(obj, cells = index)
index <- which(obj$celltype2 == "Th2_rest" & obj$Condition == "Resting")
Th2 <- subset(obj, cells = index)
index <- which(obj$celltype2 == "Th17_rest" & obj$Condition == "Resting")
Th17 <- subset(obj, cells = index)

Th1 <- FindMultiModalNeighbors(Th1, reduction.list = list("harmony_rna", "harmony_atac"), dims.list = list(1:30, 2:30))
Th1 <- RunUMAP(Th1, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
Th2 <- FindMultiModalNeighbors(Th2, reduction.list = list("harmony_rna", "harmony_atac"), dims.list = list(1:30, 2:30))
Th2 <- RunUMAP(Th2, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
Th17 <- FindMultiModalNeighbors(Th17, reduction.list = list("harmony_rna", "harmony_atac"), dims.list = list(1:30, 2:30))
Th17 <- RunUMAP(Th17, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")

pdf(paste0(outdir, "/Th1Module_Feature.pdf"), width = 10, height = 5)
FeaturePlot(Th1, reduction = "wnn.umap", features = "Rest3") + scale_color_gradient2(low = "grey90" , mid = "grey60", high = "red", midpoint = 0, limits = c(-0.05, 0.06))
dev.off()
pdf(paste0(outdir, "/Th2Module_Feature.pdf"), width = 10, height = 5)
FeaturePlot(Th2, reduction = "wnn.umap", features = "Rest4") + scale_color_gradient2(low = "grey90" , mid = "grey60", high = "red", midpoint = 0, limits = c(-0.05, 0.06))
dev.off()
pdf(paste0(outdir, "/Th17Module_Feature.pdf"), width = 10, height = 5)
FeaturePlot(Th17, reduction = "wnn.umap", features = "Rest3") + scale_color_gradient2(low = "grey90" , mid = "grey60", high = "red", midpoint = 0, limits = c(-0.05, 0.06))
dev.off()

pdf(paste0(outdir, "/Th1Module_Density.pdf"), width = 4, height = 4)
plot_density(Th1, reduction = "wnn.umap", features = "Rest3") + scale_color_gradient2(low = "grey90" , mid = "grey60", high = "red", midpoint = 0, limits = c(-0.05, 0.06))
dev.off()
pdf((outdir, "/Th2Module_Density.pdf"), width = 4, height = 4)
plot_density(Th2, reduction = "wnn.umap", features = "Rest4") + scale_color_gradient2(low = "grey90" , mid = "grey60", high = "red", midpoint = 0, limits = c(-0.05, 0.06))
dev.off()
pdf((outdir, "/Th17Module_Density.pdf"), width = 4, height = 4)
plot_density(Th17, reduction = "wnn.umap", features = "Rest3") + scale_color_gradient2(low = "grey90" , mid = "grey60", high = "red", midpoint = 0, limits = c(-0.05, 0.06))
dev.off()


