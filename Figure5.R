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

outdir <- "Figure5"
indir <- "input"

######## Figure 5A: Heatmap of peaks differential between memory and naive at rest
cluster_df <- read.table(paste0(indir, "/ATAC_Resting_clustering.txt"), header = T)
counts <- read.table(paste0(indir, "/ATAC_counts_vst.txt"))
meta <- read.table(paste0(indir, "/meta_filter.txt"))

celltypes <- c("Naive_rest","Naive_rest2","TCM_rest","Th1_rest","Th2_rest","Th17_rest","Naive_act",
                "TCM_act","TCM/TEM","Th1_act","Th2_act","Th17_act","TEM_act",
                "TEM_act2","Braking","Treg")
celltype_timepoint <- c("Naive_rest_Resting",
                        "Naive_rest2_Resting",
                        "TCM_rest_Resting",
                        "Th1_rest_Resting",
                        "Th2_rest_Resting",
                        "Th17_rest_Resting",
                        "Naive_act_2","Naive_act_5","Naive_act_15",
                        "TCM_act_2","TCM_act_5","TCM_act_15",
                        "TCM/TEM_2","TCM/TEM_5","TCM/TEM_15",
                        "Th1_act_2","Th1_act_5","Th1_act_15",
                        "Th2_act_2","Th2_act_5","Th2_act_15",
                        "Th17_act_2","Th17_act_5","Th17_act_15",
                        "TEM_act_Resting","TEM_act_2","TEM_act_5","TEM_act_15",
                        "TEM_act2_Resting","TEM_act2_2","TEM_act2_5","TEM_act2_15",
                        "Braking_15",
                        "Treg_Resting","Treg_2","Treg_5","Treg_15")
## Params
var <- "celltype" # Variable in metadata to sort and make heatmap
var_index <- which(colnames(meta) == var)
nonrep_var <- "celltype_condition" # The variable in meta that doesn't include replicate
nonrep_var_index <- which(colnames(meta) == nonrep_var)
gene_cutoff <- 5.5 # Min expression value to plot
heat_col_tfmrna <- colorRamp2(c(-2, 0, 2), c('dodgerblue3','white','red')) ## Colors for heatmap
subset_cells <- which(meta$celltype_condition %in% celltype_timepoint) ## Cell types to subset NULL if all)

## Subset counts if needed
if(!is.null(subset_cells) ){
    counts <- counts[,subset_cells]
    meta <- meta[subset_cells,]
}
meta$celltype[which(meta$celltype == "Naive_rest")] <- "Naive"
meta$celltype[which(meta$celltype == "Naive_act")] <- "Naive"
meta$celltype[which(meta$celltype == "Th1_rest")] <- "Th1"
meta$celltype[which(meta$celltype == "Th1_act")] <- "Th1"
meta$celltype[which(meta$celltype == "Th2_rest")] <- "Th2"
meta$celltype[which(meta$celltype == "Th2_act")] <- "Th2"
meta$celltype[which(meta$celltype == "Th17_rest")] <- "Th17"
meta$celltype[which(meta$celltype == "Th17_act")] <- "Th17"
meta$celltype[which(meta$celltype == "TCM_rest")] <- "TCM"
meta$celltype[which(meta$celltype == "TCM_act")] <- "TCM"
meta$celltype_condition <- paste(meta$celltype, meta$Condition, sep = "_")
rownames(meta) <- paste(meta$celltype, meta$Condition, meta$Donor, sep = "_")
colnames(counts) <- rownames(meta)
celltype_timepoint <- c("Naive_Resting","TCM_Resting","Th1_Resting","Th2_Resting","Th17_Resting", "TEM_act_Resting","TEM_act2_Resting","Treg_Resting")
celltypes <- c("Naive","TCM","Th1","Th2","Th17","TEM_act","TEM_act2","Treg")
subset_cells <- which(meta$celltype_condition %in% celltype_timepoint) ## Cell types to subset NULL if all)
## Subset counts if needed (should be needed)
if(!is.null(subset_cells) ){
    counts <- counts[,subset_cells]
    meta <- meta[subset_cells,]
}
counts <- as.matrix(counts)
colnames(counts) <- gsub("\\.", "/", colnames(counts))

## Zscore counts
counts_z <- t(scale(t(counts)))

## Create Annotation matrix
annotation_matrix <- as.matrix(unique(rownames(meta)))
rownames(annotation_matrix) <- annotation_matrix[,1]
index <- match(rownames(annotation_matrix), rownames(meta))
annotation_matrix[,1] <- (meta[,var])[index]
annotation_matrix <- cbind(annotation_matrix, meta$celltype_condition[index])
annotation_matrix <- cbind(annotation_matrix, meta$celltype_donor[index])

index <- which(celltype_timepoint %in% annotation_matrix[,2])
celltype_timepoint <- celltype_timepoint[index]
index <- which(annotation_matrix[,2] %in% celltype_timepoint)
annotation_matrix <- annotation_matrix[index,]
heatmap_colors <- c("#E41A1C","#4A72A6","#48A462","#7E6E85","#D16948", "#777777","#888888","#999999")
heatmap_colors <- setNames(heatmap_colors, c("Naive","TCM","Th1","Th2","Th17","TEM_act","TEM_act2","Treg"))
ha <- HeatmapAnnotation(celltype = annotation_matrix[,1], col = list(celltype = heatmap_colors), simple_anno_size = unit(1, "cm"),annotation_legend_param = list(celltype = list(at = celltypes)))

pdf(paste0(outdir, "/RestingPeaksHeatmap.pdf"), width = 10, height = 12)
ht <- Heatmap(counts_z[cluster_df[,1], ], name = "Z-score", 
     show_row_names = F,  show_row_dend = F, show_column_dend = F, cluster_rows = T, row_split = factor(cluster_df[,2], levels = c("C1","C2","C3","C4","C5","C6")),
     show_column_names = F, column_split = factor(annotation_matrix[,2], levels = celltype_timepoint), cluster_row_slices = F,
     col = heat_col_tfmrna, row_gap = unit(0, "mm"), border = T, column_gap = unit(0, "mm"), 
     cluster_columns = F , top_annotation = ha, row_names_gp = (gpar(fontsize = 14)), use_raster = T)
draw(ht)
dev.off()

## Associated barplot
cluster <- c()
celltype <- c()
means <- c()
sd <- c()
for(i in unique(cluster_df[,2])){
    counts_subset <- counts[which(rownames(counts) %in% cluster_df[which(cluster_df[,2] == i),1]), ]
    for(j in unique(meta$celltype)){
        cluster <- c(cluster, i)
        celltype <- c(celltype, j)
        curr_means <- rowMeans(counts_subset[,which(meta$celltype == j)])
        means <- c(means, mean(curr_means))
    }
}

############ Figure 5C: IFNG enhancer accessibility in naive and TEM
obj <- readRDS(paste0(indir, "/integrated2.rds"))
obj$celltype2[which(obj$celltype2 %in% c("Naive_rest", "Naive_rest2","Naive_act"))] <- "Naive"
obj$celltype2[which(obj$celltype2 %in% c("Th1_rest", "Th1_act"))] <- "Th1"
obj$celltype2[which(obj$celltype2 %in% c("Th2_rest", "Th2_act"))] <- "Th2"
obj$celltype2[which(obj$celltype2 %in% c("Th17_rest", "Th17_act"))] <- "Th17"
Idents(obj) <- obj$celltype2
obj_subset <- subset(obj, cells = which(obj$celltype2 %in% c("Naive","Th1","Th2","Th17")))
obj_subset$celltype_condition <- paste(obj_subset$celltype2, obj_subset$Condition, sep = "_")
Idents(obj_subset) <- obj_subset$celltype_condition
levels(obj_subset) <- c("Naive_Resting","Naive_2","Naive_5","Naive_15","Th1_Resting","Th1_2","Th1_5","Th1_15",
    "Th2_Resting","Th2_2","Th2_5","Th2_15","Th17_Resting","Th17_2","Th17_5","Th17_15")

cols <- c("#E41A1C","#E41A1C","#E41A1C","#E41A1C","#48A462","#48A462","#48A462","#48A462","#7E6E85","#7E6E85","#7E6E85","#7E6E85","#D16948","#D16948","#D16948","#D16948")
cov_plot <- CoveragePlot(obj_subset, assay = "ATAC", region = "IFNG", extend.upstream = 1000, extend.downstream = 20000, ranges.group.by = "celltype2") & scale_fill_manual(values = cols)
expression_plot <- ExpressionPlot(obj_subset, assay = "RNA", features = "IFNG")  + scale_fill_manual(values = cols)
pdf(paste0(outdir,"/IFNG.pdf"), width = 15, height = 7)
cov_plot + expression_plot
dev.off()

###### Figure5D: Barplot of RR and Resting overlaps and enrichment
# Read in all revelent peak clusters
RR_peaks <- read.table(paste0(indir, "/ATAC_RR_clustering.txt"), header = T)
Stim_peaks <- read.table(paste0(indir, "/ATAC_Stim_clustering.txt"), header = T)
Resting_peaks <- read.table(paste0(indir, "/ATAC_Resting_clustering.txt"), header = T)
mem_up <- Resting_peaks[which(Resting_peaks[,2] %in% c("C3","C4","C5","C6")),1]
mem_down <- Resting_peaks[which(Resting_peaks[,2] %in% c("C1","C2")),1]
RR_set_names <- c("EarlyRR","MiddleRR","MidTEMRR","LateStrongRR","LateWeakRR","CloseRR","CloseTEMRR")
Stim_set_names <- c("Naive+TCMStim","AllStim","LateStim","CloseAll","NaiveClose")

RR_memup <- c()
RR_memdown <- c()
RR_not <- c()
set_name <- c()
type <- c()
for(i in 1:length(unique(RR_peaks[,2]))){
    curr_peaks <- RR_peaks[which(RR_peaks[,2] == unique(RR_peaks[,2])[i]),1]
    RR_memup <- c(RR_memup, length(which(curr_peaks %in% mem_up)))
    RR_memdown <- c(RR_memdown, length(which(curr_peaks %in% mem_down)))
    RR_not <- c(RR_not, length(curr_peaks) - length(which(curr_peaks %in% mem_up)) - length(which(curr_peaks %in% mem_down)))
    #set_name <- c(set_name, rep(RR_set_names[i], 3))
}

Stim_memup <- c()
Stim_memdown <- c()
Stim_not <- c()
for(i in 1:length(unique(Stim_peaks[,2]))){
    curr_peaks <- Stim_peaks[which(Stim_peaks[,2] == unique(Stim_peaks[,2])[i]),1]
    Stim_memup <- c(Stim_memup, length(which(curr_peaks %in% mem_up)))
    Stim_memdown <- c(Stim_memdown, length(which(curr_peaks %in% mem_down)))
    Stim_not <- c(Stim_not, length(curr_peaks) - length(which(curr_peaks %in% mem_up)) - length(which(curr_peaks %in% mem_down)))
    #set_name <- c(set_name, rep(Stim_set_names[i], 3))
}
set_name <- c(set_name, rep(RR_set_names, 3))
set_name <- c(set_name, rep(Stim_set_names,3))
type <- c(rep("memup",7), rep("memdown", 7), rep("none",7))
type <- c(type, rep("memup",5), rep("memdown", 5), rep("none",5))
Peak_number <- c(RR_memup, RR_memdown, RR_not, Stim_memup, Stim_memdown, Stim_not)
df <- data.frame(PeakNumber = Peak_number, Type = type, SetName = set_name)
df$SetName <- factor(df$SetName, levels = c(RR_set_names, Stim_set_names))
df$Type <- factor(df$Type, levels = c("memup","memdown","none"))

pdf(paste0(outdir, "/PeakOverlapBar.pdf"), width = 8, height = 4)
ggplot(df, aes(x = SetName, y = PeakNumber, fill = Type)) + geom_bar(stat = "identity") + 
    scale_fill_manual(values = c("#EF8786","#6699CB", "#999999")) + 
    theme(panel.background = element_rect("white", fill = NA)) 
dev.off()
