### Figure3 A ###
library(reshape2)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(matrixStats)
library(dplyr)
library(dendextend)
library(plotrix)
library(matrixStats)
library(circlize)
library(RColorBrewer)
library(Seurat)
library(Signac)
library(ggplot2)
source("PlotPatterns.R")

outdir <- "output/Figure4"
indir <- "input"

########### Figure 4A: Heatmap of rapid recall peaks
counts <- read.table(paste0(indir, "/ATAC_counts_vst.txt"))
meta <- read.table(paste0(indir, "/meta_filter.txt"))

celltypes <- c("Naive_rest","Naive_rest2","TCM_rest","Th1_rest","Th2_rest","Th17_rest","Naive_act",
                "TCM_act","TCM/TEM","Th1_act","Th2_act","Th17_act","TEM_act",
                "TEM_act2","Treg")
celltype_timepoint <- c("Naive_rest_Resting","TCM_rest_Resting","Th1_rest_Resting","Th2_rest_Resting","Th17_rest_Resting", "TEM_act_Resting","TEM_act2_Resting",
                        "Naive_act_2","Naive_act_5","Naive_act_15","TCM_act_2","TCM_act_5","TCM_act_15","TCM/TEM_2","TCM/TEM_5","TCM/TEM_15",
                        "Th1_act_2","Th1_act_5","Th1_act_15","Th2_act_2","Th2_act_5","Th2_act_15","Th17_act_2","Th17_act_5","Th17_act_15",
                        "TEM_act_2","TEM_act_5","TEM_act_15","TEM_act2_2","TEM_act2_5","TEM_act2_15",
                        "Treg_Resting","Treg_2","Treg_5","Treg_15")
var <- "celltype" # Variable in metadata to sort and make heatmap
var_index <- which(colnames(meta) == var)
nonrep_var <- "celltype_condition" # The variable in meta that doesn't include replicate
nonrep_var_index <- which(colnames(meta) == nonrep_var)
heat_col_tfmrna <- colorRamp2(c(-2, 0, 2), c('dodgerblue3','white','red')) ## Colors for heatmap
subset_cells <- which(meta$celltype_condition %in% celltype_timepoint) ## Cell types to subset NULL if all)

## Subset counts if needed
if(!is.null(subset_cells) ){
    counts <- counts[,subset_cells]
    meta <- meta[subset_cells,]
}

# Modify celltype names to their cell class
celltype_major <- read.table(paste0(indir, "/celltype_major.txt"))$V1
meta$celltype <- celltype_major
meta$celltype_condition <- paste(meta$celltype, meta$Condition, sep = "_")
rownames(meta) <- paste(meta$celltype, meta$Condition, meta$Donor, sep = "_")
colnames(counts) <- rownames(meta)
celltype_timepoint <- c("Naive_Resting","Naive_2","Naive_5","Naive_15",
                        "TCM_Resting","TCM_2","TCM_5","TCM_15",
                        "Th1_Resting","Th1_2","Th1_5","Th1_15",
                        "Th2_Resting","Th2_2","Th2_5","Th2_15",
                        "Th17_Resting","Th17_2","Th17_5","Th17_15",
                        "TEM_act_Resting","TEM_act_2","TEM_act_5","TEM_act_15",
                        "TEM_act2_Resting","TEM_act2_2","TEM_act2_5","TEM_act2_15",
                        "Treg_Resting","Treg_2","Treg_5","Treg_15",
                        "TCM/TEM_2","TCM/TEM_5","TCM/TEM_15")
celltypes <- c("Naive","TCM","Th1","Th2","Th17","TEM_act","TEM_act2","Treg","TCM/TEM")
subset_cells <- which(meta$celltype_condition %in% celltype_timepoint) ## Cell types to subset NULL if all)
## Subset counts if needed (should be needed)
if(!is.null(subset_cells) ){
    counts <- counts[,subset_cells]
    meta <- meta[subset_cells,]
}
## Make sure col names are correct. In my colnames, there are some dots so i need to convert them to /
counts <- as.matrix(counts)
colnames(counts) <- gsub("\\.", "/", colnames(counts))
## Zscore counts
counts <- t(scale(t(counts)))

## Create Annotation
annotation_matrix <- as.matrix(unique(rownames(meta)))
rownames(annotation_matrix) <- annotation_matrix[,1]
index <- match(rownames(annotation_matrix), rownames(meta))
annotation_matrix[,1] <- (meta[,var])[index]
annotation_matrix <- cbind(annotation_matrix, meta$celltype_condition[index])
annotation_matrix <- cbind(annotation_matrix, meta$celltype_donor[index])
annotation_matrix <- cbind(annotation_matrix, meta$Condition[index])
index <- which(celltype_timepoint %in% annotation_matrix[,2])
celltype_timepoint <- celltype_timepoint[index]
index <- which(annotation_matrix[,2] %in% celltype_timepoint)
annotation_matrix <- annotation_matrix[index,]
getPalette = colorRampPalette(brewer.pal(15, "Set1"))
heatmap_colors = c("#E41A1C","#4A72A6","#48A462","#7E6E85","#D16948","#E1C62F","#B75F49","#EC83BA","#999999")
heatmap_colors <- setNames(heatmap_colors, c("Naive","TCM","Th1","Th2","Th17","TEM_act","TEM_act2","Treg","TCM/TEM"))
condition_colors <- c("#ebebeb","#C0C0C0","#A9A9A9","#808080")
names(condition_colors) <- c("Resting","2","5","15")
ha <- HeatmapAnnotation(condition = annotation_matrix[,4], celltype = annotation_matrix[,1], col = list(condition = condition_colors, celltype = heatmap_colors), simple_anno_size = unit(1, "cm"),annotation_legend_param = list(celltype = list(at = celltypes)))

# Column split sizes for heatmap
splits <- unit(c(0,0,0,1.5,0,0,0,1.5,0,0,0,1.5,0,0,0,1.5,0,0,0,1.5,0,0,0,1.5,0,0,0,1.5,0,0,0,1.5,0,0,0), "mm")

order_matrix <- read.table(paste0(indir,"/ATAC_RR_clustering.txt"), header = T)
order <- c("C1","C2","C3","C4","C5","C6","C7")

pdf(paste0(outdir, "/ATAC_RR_heatmap.pdf"), width = 10, height = 10)
ht <- Heatmap(counts[order_matrix[,1],], name = "Z-score", show_row_names = F,  show_row_dend = F, show_column_dend = F, 
              row_split = factor(order_matrix[,2], levels = order), border = T, 
              show_column_names = F, column_split = factor(annotation_matrix[,2], levels = celltype_timepoint),
              cluster_row_slices = F, col = heat_col_tfmrna, column_gap = splits, 
              row_gap = unit(0, "mm"),
              cluster_columns = F , top_annotation = ha, row_names_gp = (gpar(fontsize = 8)), use_raster = T,
              column_title = NULL, row_title = NULL)
draw(ht)
dev.off()

############## Figure S5A: Heatmap of TCR stimulated peaks

order_matrix <- read.table(paste0(indir,"/ATAC_Stim_clustering.txt"), header = T)

## Create heatmap
pdf(paste0(outdir, "/ATAC_Stim_heatmap.pdf"), width = 10, height = 12)
ht <- Heatmap(counts[order_matrix[,1],], name = "Z-score", show_row_names = F,  show_row_dend = F, show_column_dend = F, 
              row_split = factor(order_matrix[,2], levels = c("C1","C2","C3","C4","C5")),
              show_column_names = F, column_split = factor(annotation_matrix[,2], levels = celltype_timepoint),
              col = heat_col_tfmrna, column_gap = splits, row_gap = unit(0, "mm"), cluster_row_slices = F,
              cluster_columns = F , top_annotation = ha, row_names_gp = (gpar(fontsize = 14)), border = T,
              column_title = NULL, row_title = NULL, cluster_rows = T)
draw(ht)
dev.off()

######## Line plot patterns for figure 4A and S5A
outs <- c(paste0(outdir, "/RR_Cluster1.pdf"),
    paste0(outdir, "/RR_Cluster2.pdf"),
    paste0(outdir, "/RR_Cluster3.pdf"),
    paste0(outdir, "/RR_Cluster4.pdf"),
    paste0(outdir, "/RR_Cluster5.pdf"), 
    paste0(outdir, "/RR_Cluster6.pdf"),
    paste0(outdir, "/RR_Cluster7.pdf"))

clusters <- read.table(paste0(indir,"/ATAC_RR_clustering.txt"))
Cluster1_RR <- clusters[which(clusters[,2] == "C1"),1]
Cluster2_RR <- clusters[which(clusters[,2] == "C2"),1]
Cluster3_RR <- clusters[which(clusters[,2] == "C3"),1]
Cluster4_RR <- clusters[which(clusters[,2] == "C4"),1]
Cluster5_RR <- clusters[which(clusters[,2] == "C5"),1]
Cluster6_RR <- clusters[which(clusters[,2] == "C6"),1]
Cluster7_RR <- clusters[which(clusters[,2] == "C7"),1]

RR_Clusters <- list(Cluster1_RR, Cluster2_RR, Cluster3_RR, Cluster4_RR, Cluster5_RR, Cluster6_RR, Cluster7_RR)
for(i in 1:length(RR_Clusters)){
    print(i)
    plotPattern(RR_Clusters[[i]], meta, counts, outs[i])
}

outs <- c(paste0(outdir, "/Stim_Cluster1.pdf"),
        paste0(outdir, "/Stim_Cluster2.pdf"),
        paste0(outdir, "/Stim_Cluster3.pdf"),
        paste0(outdir, "/Stim_Cluster4.pdf"),
        paste0(outdir, "/Stim_Cluster5.pdf"))

clusters <- read.table(paste0(indir,"/ATAC_Stim_clustering.txt"))
Cluster1_Naive <- clusters[which(clusters[,2] == "C1"),1]
Cluster2_Naive <- clusters[which(clusters[,2] == "C2"),1]
Cluster3_Naive <- clusters[which(clusters[,2] == "C3"),1]
Cluster4_Naive <- clusters[which(clusters[,2] == "C4"),1]
Cluster5_Naive <- clusters[which(clusters[,2] == "C5"),1]

Naive_clusters <- list(Cluster1_Naive, Cluster2_Naive, Cluster3_Naive, Cluster4_Naive)
for(i in 1:length(Naive_clusters)){
    plotPattern(Naive_clusters[[i]], meta, counts, outs[i])
}


### Figure4 C ###

obj <- readRDS(paste0(indir, "/integrated2.rds"))
# Subset object for just Naive, Th1, Th2, and Th17
obj$celltype2[which(obj$celltype2 %in% c("Naive_rest", "Naive_rest2","Naive_act"))] <- "Naive"
obj$celltype2[which(obj$celltype2 %in% c("Th1_rest", "Th1_act"))] <- "Th1"
obj$celltype2[which(obj$celltype2 %in% c("Th2_rest", "Th2_act"))] <- "Th2"
obj$celltype2[which(obj$celltype2 %in% c("Th17_rest", "Th17_act"))] <- "Th17"
Idents(obj) <- obj$celltype2
obj_subset <- subset(obj, cells = which(obj$celltype2 %in% c("Naive","Th1","Th2","Th17")))
# Create celltype_condition metadata for coverage plots
obj_subset$celltype_condition <- paste(obj_subset$celltype2, obj_subset$Condition, sep = "_")
Idents(obj_subset) <- obj_subset$celltype_condition
levels(obj_subset) <- c("Naive_Resting","Naive_2","Naive_5","Naive_15","Th1_Resting","Th1_2","Th1_5","Th1_15","Th2_Resting","Th2_2","Th2_5","Th2_15","Th17_Resting","Th17_2","Th17_5","Th17_15")
# MAF
cols <- c("#E41A1C","#E41A1C","#E41A1C","#E41A1C","#449B75","#449B75","#449B75","#449B75","#6B886D","#6B886D","#6B886D","#6B886D","#AC5782","#AC5782","#AC5782","#AC5782")
pdf(paste0(outdir,"/MAF.pdf"), width = 10, height = 7)
CoveragePlot(obj_subset, assay = "ATAC", region = "chr16-79595000-79597000", extend.upstream = 100, extend.downstream = 100, ranges.group.by = "celltype2") & scale_fill_manual(values = cols)
dev.off()
# FYN
pdf(paste0(outdir,"/FYN.pdf"), width = 10, height = 7)
CoveragePlot(obj_subset, assay = "ATAC", region = "chr6-111833500-111840000", extend.upstream = 700, extend.downstream = 700, ranges.group.by = "celltype2") & scale_fill_manual(values = cols)
dev.off()
# LEF1
pdf(paste0(outdir,"/LEF1.pdf"), width = 10, height = 7)
CoveragePlot(obj_subset, assay = "ATAC", region = "chr4-108070000-108080000", extend.upstream = 5000, extend.downstream = 5000, ranges.group.by = "celltype2") & scale_fill_manual(values = cols)
dev.off()
# TIAM1
pdf(paste0(outdir,"/TIAM1.pdf"), width = 10, height = 7)
CoveragePlot(obj_subset, assay = "ATAC", region = "chr21-31115000-31120000", extend.upstream = 0, extend.downstream = 0, ranges.group.by = "celltype2") & scale_fill_manual(values = cols)
dev.off()
# MYC
pdf(paste0(outdir,"/MYC.pdf"), width = 10, height = 7)
CoveragePlot(obj_subset, assay = "ATAC", region = "MYC", extend.upstream = 15000, extend.downstream = 15000, ranges.group.by = "celltype2") & scale_fill_manual(values = cols)
dev.off()
# LY96
pdf(paste0(outdir,"/LY96.pdf"), width = 10, height = 7)
CoveragePlot(obj_subset, assay = "ATAC", region = "chr8-73992216-73992352", extend.upstream = 5000, extend.downstream = 5000, ranges.group.by = "celltype2") & scale_fill_manual(values = cols)
dev.off()
# RAPGEF2
pdf(paste0(outdir,"/RAPGEF2.pdf"), width = 10, height = 7)
CoveragePlot(obj_subset, assay = "ATAC", region = "chr4-159078197-159078295", extend.upstream = 5000, extend.downstream = 5000, ranges.group.by = "celltype2") & scale_fill_manual(values = cols)
dev.off()
# SSBP4
pdf(paste0(outdir,"/SSBP4.pdf"), width = 10, height = 7)
CoveragePlot(obj_subset, assay = "ATAC", region = "chr19-18428714-18428879", extend.upstream = 3500, extend.downstream = 3500, ranges.group.by = "celltype2") & scale_fill_manual(values = cols)
dev.off()
# CD58
pdf(paste0(outdir,"/CD58.pdf"), width = 10, height = 7)
CoveragePlot(obj_subset, assay = "ATAC", region = "chr1-116544574-116544789", extend.upstream = 3000, extend.downstream = 3000, ranges.group.by = "celltype2") & scale_fill_manual(values = cols)
dev.off()
# GALM
pdf(paste0(outdir,"/GALM.pdf"), width = 10, height = 7)
CoveragePlot(obj_subset, assay = "ATAC", region = "chr11-116930170-116930436", extend.upstream = 3000, extend.downstream = 3000, ranges.group.by = "celltype2") & scale_fill_manual(values = cols)
dev.off()
# CCR6
pdf(paste0(outdir,"/CCR6.pdf"), width = 10, height = 7)
CoveragePlot(obj_subset, assay = "ATAC", region = "CCR6", extend.upstream = 3000, extend.downstream = 3000, ranges.group.by = "celltype2") & scale_fill_manual(values = cols)
dev.off()
pdf(paste0(outdir,"/VIM.pdf"), width = 10, height = 7)
CoveragePlot(obj_subset, assay = "ATAC", region = "chr10-17232000-17234000", extend.upstream = 0, extend.downstream = 0, ranges.group.by = "celltype2") & scale_fill_manual(values = cols)
dev.off()

