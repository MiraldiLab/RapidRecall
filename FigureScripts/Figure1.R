#Load in required packages
library(Seurat)
library(Signac)
library(ComplexHeatmap)
library(ggplot2)
library(tidyverse)
library(circlize)
library(RColorBrewer)
library(viridis)
library(Polychrome)
library(DESeq2)
library(ggrepel)

indir <- "input"
outdir <- "output/Figure1" #output director for figure 1

#Load in seurat object, pseudobulked data, and celltype colors
obj <- readRDS(paste0(indir,"/integrated2.rds"))
RNA_counts <- read.table(paste0(indir, "/RNA_counts_combatseq_vst.txt"))
ATAC_counts <- read.table(paste0(indir, "/ATAC_counts_vst.txt"))
meta <- read.table(paste0(indir, "/meta_filter.txt"))
meta_raw <- read.table(paste0(indir, "/meta_raw.txt"))
colors <- read.table(paste0(indir, "/Colors.txt"))
colors <- colors[-6, ] #Removed "Stressed" from color list
celltype_labels <- read.table(paste0(indir, "/celltype_labels.txt")) # celltype annotation for each barcode
celltype_labels[which(celltype_labels[,2] == "Naive_rest2"),2] <- "Naive_rest"

# Initialize assay and identities of the seurat object
DefaultAssay(obj) <- "RNA"
obj$celltype2 <- celltype_labels[,2]
Idents(obj) <- "celltype2"
levels(obj) <- c("Naive_rest","Naive_act","Stressed","TCM_rest","TCM_act","Braking","Th1_rest","Th1_act","Th2_rest","Th2_act","Th17_rest","Th17_act",
    "MHCII_rest","MHCII_act","CTL_rest","CTL_act","TCM/TEM","Treg")

###### FigS1A-C: UMAPs of unharmonized and harmonized data
# Calculate unharmonized RNA dimension reductions
obj <- NormalizeData(obj)
obj <- ScaleData(obj, features = rownames(obj))
obj <- FindVariableFeatures(obj)
obj <- RunPCA(obj, dims = 1:30)
obj <- RunUMAP(obj, dims = 1:30)
# Calculate unharmonized ATAC dimension reductions
DefaultAssay(obj) <- "peaks"
obj <- RunTFIDF(obj)
obj <- FindTopFeatures(obj, min.cutoff = 50)
obj <- RunSVD(obj)
obj <- RunUMAP(obj, reduction = "lsi", reduction.name = "umap_atac", dims = 2:30)
# Calculate weighted nearest neighbor graph
obj <- FindMultiModalNeighbors(obj, reduction.list = list("pca","lsi"), dims.list = list(1:30, 2:30))
obj <- RunUMAP(obj, nn.name = "weighted.nn", reduction.name = "wnn.umap_noharmony", reduction.key = "wnnUMAP_")
pdf(paste0(outdir, "/unharmonized_umap.pdf"), width = 8, height = 8)
DimPlot(obj, reduction = "wnn.umap_noharmony", group.by = "Donor")
dev.off()
pdf(paste0(outdir, "/harmonized_umap_donor.pdf"), width = 8, height = 8)
DimPlot(obj, reduction = "wnn.umap", group.by = "Donor")
dev.off()
pdf(paste0(outdir, "/harmonized_umap_celltype.pdf"), width = 8, height = 8)
DimPlot(obj, reduction = "wnn.umap")
dev.off()

###### FigS1D-E: Heatmaps of top donor variable genes acorss different batch correction designs
# Load in uncorrect and the different correction counts files
uncorrected <- read.table(paste0(indir, "/uncorrected_vst.txt"))
corrected_donor <- read.table(paste0(indir, "/corrected_donor_vst.txt"))
corrected_dataset <- read.table(paste0(indir, "/corrected_dataset_vst.txt"))
simulated <- read.table(paste0(indir, "/simulated_vst.txt"))
meta <- read.table(paste0(indir, "/meta_batch.txt"))
genes_dynamic <- read.table(paste0(indir, "/GeneDynamics_df_Filtered.txt"), header = T)
colnames(uncorrected) <- rownames(meta)
colnames(corrected_donor) <- rownames(meta)
colnames(corrected_dataset) <- rownames(meta)
colnames(simulated) <- rownames(meta)
genes <- unique(genes_dynamic$Gene)
uncorrected <- as.matrix(uncorrected)
corrected_donor <- as.matrix(corrected_donor)
corrected_dataset <- as.matrix(corrected_dataset)
# Manually order the counts matrices 
celltypes <- c("Naive_rest","Naive_rest2","TCM_rest","Th1_rest","Th2_rest","Th17_rest","Naive_act",
                "TCM_act","TCM/TEM","Th1_act","Th2_act","Th17_act","TEM_act",
                "TEM_act2","Treg")
celltype_timepoint <- c("Naive_rest_Resting","Naive_rest2_Resting","TCM_rest_Resting","Th1_rest_Resting","Th2_rest_Resting","Th17_rest_Resting", "TEM_act_Resting","TEM_act2_Resting",
                        "Naive_act_2","Naive_act_5","Naive_act_15","TCM_act_2","TCM_act_5","TCM_act_15","TCM/TEM_2","TCM/TEM_5","TCM/TEM_15",
                        "Th1_act_2","Th1_act_5","Th1_act_15","Th2_act_2","Th2_act_5","Th2_act_15","Th17_act_2","Th17_act_5","Th17_act_15",
                        "TEM_act_5","TEM_act_15","TEM_act2_2","TEM_act2_5","TEM_act2_15",
                        "Treg_Resting","Treg_2","Treg_5","Treg_15")
meta <- meta[which(meta$celltype_condition %in% celltype_timepoint),]
uncorrected <- uncorrected[,rownames(meta)]
corrected_donor <- corrected_donor[,rownames(meta)]
corrected_dataset <- corrected_dataset[,rownames(meta)]
simulated <- simulated[,rownames(meta)]
# Calculate standard deviations across donors of genes from uncorrected counts matrix
# Only look at genes that were dynamic and met minimum expression criteria
mean_sd <- c()
for(j in genes){
    print(j)
    sds <- c()
    for(i in unique(meta$celltype_condition)){
        curr_samples <- rownames(meta)[which(meta$celltype_condition == i)]
        index <- which(colnames(uncorrected) %in% curr_samples)
        curr_counts <- uncorrected[j, index]
        sds <- c(sds, sd(as.vector(curr_counts)))
    }
    mean_sd <- c(mean_sd, mean(sds))
}
# Find top 1000 most variable genes
index <- order(mean_sd, decreasing = T)
top_1000 <- genes[index[1:1000]]
# Plot the different counts matrices
var <- "celltype" # Variable in metadata to sort and make heatmap
var_index <- which(colnames(meta) == var)
nonrep_var <- "celltype_condition" # The variable in meta that doesn't include replicate
nonrep_var_index <- which(colnames(meta) == nonrep_var)
heat_col_tfmrna <- colorRamp2(c(-2, 0, 2), c('dodgerblue','white','red')) ## Colors for heatmap
subset_cells <- which(meta$celltype_condition %in% celltype_timepoint) ## Cell types to subset NULL if all)
## Subset counts if needed
if(!is.null(subset_cells) ){
    uncorrected <- uncorrected[,subset_cells]
    corrected_donor <- corrected_donor[,subset_cells]
    corrected_dataset <- corrected_dataset[,subset_cells]
    meta <- meta[subset_cells,]
}
# Modify celltype names
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
colnames(uncorrected) <- rownames(meta)
colnames(corrected_donor) <- rownames(meta)
colnames(corrected_dataset) <- rownames(meta)
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
    uncorrected <- uncorrected[,subset_cells]
    corrected_donor <- corrected_donor[,subset_cells]
    corrected_dataset <- corrected_dataset[,subset_cells]
    meta <- meta[subset_cells,]
}
## Zscore counts
uncorrected <- t(scale(t(uncorrected)))
corrected_donor <- t(scale(t(corrected_donor)))
corrected_dataset <- t(scale(t(corrected_dataset)))

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
splits <- unit(c(0,0,0,1.5,0,0,0,1.5,0,0,0,1.5,0,0,0,1.5,0,0,0,1.5,0,0,1.5,0,0,0,1.5,0,0,0,1.5,0,0,0), "mm")
# Clustering solution previously calculated
order_df <- read.table(paste0(indir, "/batch_correction_clusters.txt"), header = T)
## Create heatmap
pdf(paste0(outdir, "/uncorrected_top1000.pdf"), width = 10, height = 10)
ht <- Heatmap(uncorrected[order_df[,1],], name = "Z-score", show_row_names = F,  show_row_dend = F, show_column_dend = F, 
              row_split = factor(order_df[,2], levels = c("group1","group2","group3","group4","group5", "group6")), border = T, 
              show_column_names = F, column_split = factor(annotation_matrix[,2], levels = celltype_timepoint),
              cluster_row_slices = F, col = heat_col_tfmrna, column_gap = splits, 
              row_gap = unit(0, "mm"), cluster_rows = F,
              cluster_columns = F , top_annotation = ha, row_names_gp = (gpar(fontsize = 8)), use_raster = T,
              column_title = NULL, row_title = NULL)
draw(ht)
dev.off()
pdf(paste0(outdir, "/corrected_donor_top1000.pdf"), width = 10, height = 10)
ht <- Heatmap(corrected_dataset[order_df[,1],], name = "Z-score", show_row_names = F,  show_row_dend = F, show_column_dend = F, 
              row_split = factor(order_df[,2], levels = c("group1","group2","group3","group4","group5", "group6")), border = T, 
              show_column_names = F, column_split = factor(annotation_matrix[,2], levels = celltype_timepoint),
              cluster_row_slices = F, col = heat_col_tfmrna, column_gap = splits, 
              row_gap = unit(0, "mm"), cluster_rows = F,
              cluster_columns = F , top_annotation = ha, row_names_gp = (gpar(fontsize = 8)), use_raster = T,
              column_title = NULL, row_title = NULL)
draw(ht)
dev.off()
pdf(paste0(outdir, "/corrected_dataset_top1000.pdf"), width = 10, height = 10)
ht <- Heatmap(corrected_donor[order_df[,1],], name = "Z-score", show_row_names = F,  show_row_dend = F, show_column_dend = F, 
              row_split = factor(order_df[,2], levels = c("group1","group2","group3","group4","group5", "group6")), border = T, 
              show_column_names = F, column_split = factor(annotation_matrix[,2], levels = celltype_timepoint),
              cluster_row_slices = F, col = heat_col_tfmrna, column_gap = splits, 
              row_gap = unit(0, "mm"), cluster_rows = F,
              cluster_columns = F , top_annotation = ha, row_names_gp = (gpar(fontsize = 8)), use_raster = T,
              column_title = NULL, row_title = NULL)
draw(ht)
dev.off()

###### Differential Genes across different correction methods (Fig. S1G-H)
## Dynamic DEGs
# Load in DESeq2 results for each count type
df_raw <- read.table(paste0(indir, "/batch_deseq_uncorrected.txt"), header = T)
df_donor <- read.table(paste0(indir, "/batch_deseq_donor.txt"), header = T)
df_dataset <- read.table(paste0(indir, "/batch_deseq_dataset.txt"), header = T)
df_sim <- read.table(paste0(indir, "/batch_deseq_simulation.txt"), header = T)
# Only look at genes expressed in at least 5 percent of the celltype of interest
five_percent <- read.table(paste0(indir, "/five_percent_genes.txt"))$V1
df_raw <- df_raw[which(df_raw$Gene %in% five_percent),]
df_donor <- df_donor[which(df_donor$Gene %in% five_percent),]
df_dataset <- df_dataset[which(df_dataset$Gene %in% five_percent),]
df_sim <- df_sim[which(df_sim$Gene %in% five_percent),]
# Find number of DEGs in each celltype
celltype <- c()
design <- c()
annotation <- c()
number <- c()
for(i in unique(df_raw$Celltype)){
    genes_dataset <- unique(df_dataset$Gene[which(df_dataset$Celltype == i)])
    genes_donor <- unique(df_donor$Gene[which(df_donor$Celltype == i)])
    genes_raw <- unique(df_raw$Gene[which(df_raw$Celltype == i)])
    genes_sim <- unique(df_sim$Gene[which(df_sim$Celltype == i)])

    genes_donor_match <- genes_donor[which(genes_donor %in% genes_dataset)]
    genes_donor_nomatch <- genes_donor[which(!(genes_donor %in% genes_dataset))]
    genes_raw_match <- genes_raw[which(genes_raw %in% genes_dataset)]
    genes_raw_nomatch <- genes_raw[which(!(genes_raw %in% genes_dataset))]
    genes_sim_match <- genes_sim[which(genes_sim %in% genes_dataset)]
    genes_sim_nomatch <- genes_sim[which(!(genes_sim %in% genes_dataset))]

    celltype <- c(celltype, rep(i, 4))
    design <- c(design,"Dataset","Donor","Raw","Sim")
    number <- c(number, length(genes_dataset), length(genes_donor_match), length(genes_raw_match),length(genes_sim_match))
}
df <- data.frame(Number = number, Celltype = celltype, Design = design)
df$Celltype <- factor(df$Celltype, levels = c("Naive","TCM","Th1","Th2","Th17","TEM","TEM2","Treg"))
df$Design = factor(df$Design, levels = c("Dataset","Donor","Raw","Sim"))
pdf(paste0(outdir, "/BatchCorrection_DEG_Dynamic.pdf"), width = 12, height = 8)
ggplot(df, aes(x = Celltype, y = Number, fill = Design)) + 
    geom_bar(stat = "identity", position = "dodge", color = "black") +
    scale_fill_manual(values = c("#E63946","#A8DADC","#1D3557", "#F1FAEE")) +  # Custom colors for Design
    theme_minimal()
dev.off()

## Resting DEGs
df_raw <- read.table(paste0(indir, "/batch_deseq_uncorrected_rest.txt"), header = T)
df_donor <- read.table(paste0(indir, "/batch_deseq_donor_rest.txt"), header = T)
df_dataset <- read.table(paste0(indir, "/batch_deseq_dataset_rest.txt"), header = T)
df_sim <- read.table(paste0(indir, "/batch_deseq_simulation_rest.txt"), header = T)

df_raw <- df_raw[which(df_raw$Gene %in% five_percent),]
df_donor <- df_donor[which(df_donor$Gene %in% five_percent),]
df_dataset <- df_dataset[which(df_dataset$Gene %in% five_percent),]
df_sim <- df_sim[which(df_sim$Gene %in% five_percent),]

celltype <- c()
design <- c()
annotation <- c()
number <- c()
for(i in unique(df_raw$Celltype)){
    genes_dataset <- unique(df_dataset$Gene[which(df_dataset$Celltype == i)])
    genes_donor <- unique(df_donor$Gene[which(df_donor$Celltype == i)])
    genes_raw <- unique(df_raw$Gene[which(df_raw$Celltype == i)])
    genes_sim <- unique(df_sim$Gene[which(df_sim$Celltype == i)])

    genes_donor_match <- genes_donor[which(genes_donor %in% genes_dataset)]
    genes_donor_nomatch <- genes_donor[which(!(genes_donor %in% genes_dataset))]
    genes_raw_match <- genes_raw[which(genes_raw %in% genes_dataset)]
    genes_raw_nomatch <- genes_raw[which(!(genes_raw %in% genes_dataset))]
    genes_sim_match <- genes_sim[which(genes_sim %in% genes_dataset)]
    genes_sim_nomatch <- genes_sim[which(!(genes_sim %in% genes_dataset))]

    celltype <- c(celltype, rep(i, 4))
    design <- c(design,"Dataset","Donor","Raw","Sim")
    number <- c(number, length(genes_dataset), length(genes_donor_match), length(genes_raw_match),length(genes_sim_match))
}
df <- data.frame(Number = number, Celltype = celltype, Design = design)
df$Design = factor(df$Design, levels = c("Dataset","Donor","Raw","Sim"))

pdf(paste0(outdir, "/BatchCorrection_DEG_Resting.pdf"), width = 12, height = 8)
ggplot(df, aes(x = Celltype, y = Number, fill = Design)) + 
    geom_bar(stat = "identity", position = "dodge", color = "black") +
    scale_fill_manual(values = c("#E63946","#A8DADC","#1D3557", "#F1FAEE")) +  # Custom colors for Design
    theme_minimal()
dev.off()

###### FigS1I-J: UMI and Fragment number distributions
obj$celltype2 <- celltype_labels[,2]
# Adjust celltype labels to their major celltype names
obj$celltype2[which(obj$celltype2 %in% c("Naive_rest", "Naive_rest2", "Naive_act"))] <- "Naive"
obj$celltype2[which(obj$celltype2 %in% c("TCM_rest", "TCM_act"))] <- "TCM"
obj$celltype2[which(obj$celltype2 %in% c("Th1_rest", "Th1_act"))] <- "Th1"
obj$celltype2[which(obj$celltype2 %in% c("Th2_rest", "Th2_act"))] <- "Th2"
obj$celltype2[which(obj$celltype2 %in% c("Th17_rest", "Th17_act"))] <- "Th17"
obj$celltype2[which(obj$celltype2 %in% c("Braking"))] <- "TCM"
obj$celltype2[which(obj$celltype2 %in% c("MHCII_rest","MHCII_act"))] <- "TEM"
obj$celltype2[which(obj$celltype2 %in% c("CTL_rest","CTL_act"))] <- "TEM2"
# UMI
colors <- read.table(paste0(indir, "/Colors.txt"))
cols <- colors[,1]
cols <- paste("#", cols, sep = "")
df_umi <- data.frame(UMI = obj$nCount_RNA, celltype = obj$celltype2)
df_umi$celltype <- factor(df_umi$celltype, levels = c("Naive","TCM","Th1","Th2","Th17","Stressed","TCM/TEM", "TEM","TEM2","Treg"))
pdf(paste0(outdir, "/QC_UMI.pdf"), width = 12, height = 8)
ggplot(df_umi, aes(x = celltype, y = UMI)) + geom_violin(aes(fill = celltype)) + scale_fill_manual(values = cols) + ylim(c(0, 10000)) +
    geom_boxplot(width=.1) + 
    theme(panel.grid = element_line(color = "grey",linewidth = 0.1)) + 
    theme(panel.background = element_rect("white")) +
    theme(axis.title=element_text(size=15, face='plain')) + 
    theme(axis.text=element_text(size=14, face='plain')) 
dev.off()
# Fragments 
df_fragments <- data.frame(fragments = obj$nCount_ATAC, celltype = obj$celltype2)
df_fragments$celltype <- factor(df_fragments$celltype, levels = c("Naive","TCM","Th1","Th2","Th17","Stressed","TCM/TEM", "TEM","TEM2","Treg"))
pdf(paste0(outdir, "/QC_Fragments.pdf"), width = 12, height = 8)
ggplot(df_fragments, aes(x = celltype, y = fragments)) + geom_violin(aes(fill = celltype)) + scale_fill_manual(values = cols) + ylim(c(0, 50000)) +
    geom_boxplot(width=.1) + 
    theme(panel.grid = element_line(color = "grey",linewidth = 0.1)) + 
    theme(panel.background = element_rect("white")) +
    theme(axis.title=element_text(size=15, face='plain')) + 
    theme(axis.text=element_text(size=14, face='plain')) 
dev.off()

# Fig S1L: Barplot of major populations by timepoint
celltype_major <- celltype_labels[,2]
celltype_major[which(celltype_major %in% c("Naive_rest","Naive_rest2","Naive_act"))] <- "Naive"
celltype_major[which(celltype_major %in% c("TCM_rest","TCM_act","Braking"))] <- "TCM"
celltype_major[which(celltype_major %in% c("Th1_rest","Th1_act","Th2_rest","Th2_act","Th17_rest","Th17_act","MHCII_rest","MHCII_act","CTL_rest","CTL_act"))] <- "TEM"
timepoints <- c("Resting","2","5","15")
barplot_cellclass <- c()
barplot_timepoint <- c()
barplot_freq <- c()
for(i in timepoints){
    for(j in unique(celltype_major)){
        barplot_cellclass <- c(barplot_cellclass, j)
        barplot_timepoint <- c(barplot_timepoint, i)
        barplot_freq <- c(barplot_freq, length(which(celltype_major == j & obj$Condition == i)))
    }
}
cellclass_df <- data.frame(Class = barplot_cellclass, Timepoint = barplot_timepoint, Freq = barplot_freq)
cellclass_df$Class = factor(cellclass_df$Class, levels = c("Naive","Stressed","TCM","TEM","TCM/TEM","Treg"))
cellclass_df$Timepoint <- factor(cellclass_df$Timepoint, levels = c("Resting","2","5","15"))
pdf(paste0(outdir, "/MajorCelltype_barplot.pdf"), width = 8, height = 12)
ggplot(cellclass_df, aes(x = Timepoint, y = Freq, fill = Class)) + geom_bar(stat = "identity") + 
    scale_fill_manual(values = c("#E41A1D","#FFB716","#4A72A6","#49A462","#E1C62E","#999999"))
dev.off()

# Reset celltype labels
obj <- readRDS(paste0(indir,"/integrated2.rds"))
DefaultAssay(obj) <- "RNA"
obj$celltype2 <- celltype_labels[,2]
Idents(obj) <- "celltype2"
###### Fig1D: Expression dotplot of CD4+ T-cell population markers
colors <- read.table(paste0(indir, "/Colors.txt"))
colors <- colors[-6,]
obj$celltype2[which(obj$celltype2 == "Naive_rest2")] <- "Naive_rest" #Combine the two naive-rest populations
obj$celltype2[which(obj$celltype2 == "Braking")] <- "TCM_act" # Set Braking to be activated TCM for Fig1 
obj$celltype2[which(obj$celltype2 == "Treg" & obj$Condition != "Resting")] <- "Treg_act"
obj$celltype2[which(obj$celltype2 == "Treg")] <- "Treg_rest"
Idents(obj) <- obj$celltype2 # Reset identities after changing the celltype labels
index <- which(obj$celltype2 == "Stressed")
obj <- subset(obj, cells = index, invert = T) #remove stressed cells from the object
levels(obj) <- c("Naive_rest","TCM_rest","Th1_rest","Th2_rest","Th17_rest","MHCII_rest","CTL_rest","Treg_rest","Naive_act","TCM_act","Th1_act","Th2_act","Th17_act", "MHCII_act","CTL_act","Treg_act","TCM/TEM")
# Marker genes to include in the dotplot
genes <- c("FHIT","LEF1","TCF7","SELL","CCR7","CD27", "CCL5","IFNG-AS1","IFNG","TBX21","GATA3","IL4R","IL4","IL5","IL13","PTPN13","CCR6","IL17A","IL17F", "RORC","IL21","IRF4","CD69","CD44","IL2RA",
            "FOS","JUN","HLA-DPA1", "HLA-DQB1", "GNLY","PRF1","EOMES","ITGB1","FOXP3","CTLA4")
p <- DotPlot(object = obj, features = genes)
df <- p$data

# Calculate average expression of celltype markers
exp_mat <- df %>% 
  select(-pct.exp, -avg.exp) %>%  
  pivot_wider(names_from = id, values_from = avg.exp.scaled) %>% 
  as.data.frame() 
row.names(exp_mat) <- exp_mat$features.plot  
exp_mat <- exp_mat[,-1] %>% as.matrix()

# Calculate percent of cells that express celltype markers
percent_mat<-df %>% 
  select(-avg.exp, -avg.exp.scaled) %>%  
  pivot_wider(names_from = id, values_from = pct.exp) %>% 
  as.data.frame() 
row.names(percent_mat) <- percent_mat$features.plot  
percent_mat <- percent_mat[,-1] %>% as.matrix()

# Color palette for celltype marker dotplot
col_fun = circlize::colorRamp2(c(-1, 0, 2), c("dodgerblue", "grey", "red"))
# Function describing circle size and color based on percent expression and average expression
cell_fun = function(j, i, x, y, w, h, fill){
          grid.rect(x = x, y = y, width = w, height = h, 
                    gp = gpar(col = NA, fill = NA))
          grid.circle(x=x,y=y,r= percent_mat[i, j]/100 * min(unit.c(w, h)),
                      gp = gpar(fill = col_fun(exp_mat[i, j]), col = NA))}
# Annotations for the dotplot
cluster_anno <- c("Naive","TCM","Th1","Th2","Th17","MHCII","CTL","Treg","Naive","TCM","Th1","Th2","Th17","MHCII","CTL","Treg","TCM/TEM")
stim <- c("Resting","Resting","Resting","Resting","Resting","Resting","Resting","Resting","Activated","Activated","Activated","Activated","Activated","Activated","Activated","Activated","Activated")
getPalette = colorRampPalette(brewer.pal(length(unique(cluster_anno)), "Set1"))
heatmap_colors = paste0("#", colors[,1])
column_ha<- HeatmapAnnotation(
    Stimulation = stim,
    Celltype = cluster_anno,
    col = list(Stimulation = setNames(c("#DCDCDC","#696969"), unique(stim)), Celltype = setNames(heatmap_colors, unique(cluster_anno))),
    na_col = "grey",
    annotation_legend_param = list(Stimulation = list(at = unique(stim)), Celltype = list(at = unique(cluster_anno)))
)
# Plot Fig 1D
pdf(paste(outdir,"/DotPlot_MHC.pdf", sep = ""), width = 7, height = 9.5)
Heatmap(exp_mat,
        heatmap_legend_param=list(title="expression"),
        column_title = "clustered dotplot", 
        col=col_fun,
        rect_gp = gpar(type = "none"),
        cell_fun = cell_fun,
        row_names_gp = gpar(fontsize = 14, fontface = "italic"),
        cluster_rows = F,
        border = "black",
        cluster_columns = F,
        top_annotation = column_ha,
        #row_split = factor(gene_class, levels = unique(gene_class)),
        row_gap = unit(0, "mm")
        )
dev.off()




###### Donor resolved cell frequency barplot (not used in any figure)
obj$celltype_donor <- paste(obj$celltype2, obj$Donor ,sep = "_") # Create new celltype_donor metadata column
celltype_freq <- table(obj$celltype_donor) # Find the frequency of celltype_donors
# Order of celltype_donor populations for plotting purposes
celltype_donor <- c("Naive_rest_Donor1","Naive_rest_Donor2","Naive_rest_Donor3","Naive_rest_Donor4","Naive_act_Donor1","Naive_act_Donor2","Naive_act_Donor3","Naive_act_Donor4",
                    "TCM_rest_Donor1","TCM_rest_Donor2","TCM_rest_Donor3","TCM_rest_Donor4","TCM_act_Donor1","TCM_act_Donor2","TCM_act_Donor3","TCM_act_Donor4",
                    "Th1_rest_Donor1","Th1_rest_Donor2","Th1_rest_Donor3","Th1_rest_Donor4","Th1_act_Donor1","Th1_act_Donor2","Th1_act_Donor3","Th1_act_Donor4",
                    "Th2_rest_Donor1","Th2_rest_Donor2","Th2_rest_Donor3","Th2_rest_Donor4","Th2_act_Donor1","Th2_act_Donor2","Th2_act_Donor3","Th2_act_Donor4",
                    "Th17_rest_Donor1","Th17_rest_Donor2","Th17_rest_Donor3","Th17_rest_Donor4", "Th17_act_Donor1","Th17_act_Donor2","Th17_act_Donor3","Th17_act_Donor4",
                    "Stressed_Donor1","Stressed_Donor2","Stressed_Donor3","Stressed_Donor4",
                    "TCM/TEM_Donor1","TCM/TEM_Donor2","TCM/TEM_Donor3","TCM/TEM_Donor4",
                    "TEM_act_Donor1","TEM_act_Donor2","TEM_act_Donor3","TEM_act_Donor4",
                    "TEM_act2_Donor1","TEM_act2_Donor2","TEM_act2_Donor3","TEM_act2_Donor4",
                    "Treg_Donor1","Treg_Donor2","Treg_Donor3","Treg_Donor4")
# Transparency values (alphas) to help distinguish different donors for the barplot
alpha <- rep(c(0.7, 0.8, 0.9, 1), 15)
celltype_freq <- celltype_freq[celltype_donor] # Reorder frequency data based on the desired order
# Order of celltypes for plotting purposes
celltype <- c("Naive_rest","Naive_rest","Naive_rest","Naive_rest","Naive_act","Naive_act","Naive_act","Naive_act","TCM_rest","TCM_rest","TCM_rest","TCM_rest","TCM_act","TCM_act","TCM_act","TCM_act",
                "Th1_rest","Th1_rest","Th1_rest","Th1_rest","Th1_act","Th1_act","Th1_act","Th1_act","Th2_rest","Th2_rest","Th2_rest","Th2_rest","Th2_act","Th2_act","Th2_act","Th2_act",
                "Th17_rest","Th17_rest","Th17_rest","Th17_rest","Th17_act","Th17_act","Th17_act","Th17_act","Stressed","Stressed","Stressed","Stressed",
                "TCM/TEM","TCM/TEM","TCM/TEM","TCM/TEM","TEM_act","TEM_act","TEM_act","TEM_act","TEM_act2","TEM_act2","TEM_act2","TEM_act2",
                "Treg","Treg","Treg","Treg")

# Construct dataframe to plot in ggplot
Donor <- rep(c("Donor1","Donor2","Donor3","Donor4"),15)
df <- data.frame(Freq = as.vector(celltype_freq), CellType = celltype, Donor = Donor)
df$CellType <- factor(df$CellType, levels = unique(celltype))
df$Donor <- factor(df$Donor, levels = c("Donor1","Donor2","Donor3","Donor4"))
pdf(paste(outdir, "/CelltypeFreq.pdf", sep = ""), width = 16, height = 6)
ggplot(df, aes(x = CellType, y = Freq, fill = CellType, alpha = alpha)) + geom_bar(position = "stack", stat = "identity") +
    scale_fill_manual(values = paste0("#", colors[,1])) +
    theme(panel.grid = element_line(color = "grey",linewidth = 0.1)) + 
    theme(panel.background = element_rect("white")) +
    theme(axis.title=element_text(size=15, face='plain')) + 
    theme(axis.text=element_text(size=14, face='plain'))
dev.off()



########## Signal Tracks Fig1E
# Set default assay to ATAC and rename Treg populations to match other populations
DefaultAssay(obj) <- "ATAC"
obj$celltype2[which(obj$celltype2 == "Treg" & obj$Condition == "Resting")] <- "Treg_rest"
obj$celltype2[which(obj$celltype2 == "Treg")] <- "Treg_act"
Idents(obj) <- obj$celltype2

# Assign signal track colors
colors_signal <- colors[which(colors[,2] %in% c("Naive_rest","Naive_act","TCM_rest","TCM_act","Th1_rest","Th1_act","Th2_rest","Th2_act","Th17_rest","Th17_act","Treg_rest", "Treg_act")), ]
index <- which(obj$celltype2 %in% c("Naive_rest","Naive_act","TCM_rest","TCM_act","Th1_rest","Th1_act","Th2_rest","Th2_act","Th17_rest","Th17_act","Treg_rest", "Treg_act", "Naive_Responsive"))
obj_subset <- subset(obj, cells = index)
# Order of celltypes to plot accessability data
levels(obj_subset) <- c("Naive_rest","Naive_act","Naive_Responsive","TCM_rest","TCM_act","Th1_rest","Th1_act","Th2_rest","Th2_act","Th17_rest","Th17_act","Treg_rest", "Treg_act")
# IFNG, IL4, IL17A, IL2, FOXP3, GZMB,CD58. Extensions manual chosen to best highlight relavant peaks 
pdf(paste(outdir, "/IFNG.pdf", sep = ""), width = 7, height = 4)
CoveragePlot(obj_subset, region = "IFNG", extend.upstream = 5000, extend.downstream = 5000) & scale_fill_manual(values = paste0("#", colors_signal[,1]))
dev.off()
pdf(paste(outdir, "/IL4.pdf", sep = ""), width = 7, height = 4)
CoveragePlot(obj_subset, region = "IL4", extend.upstream = 5000, extend.downstream = 5000) & scale_fill_manual(values = paste0("#", colors_signal[,1]))
dev.off()
pdf(paste(outdir, "/IL17A.pdf", sep = ""), width = 7, height = 4)
CoveragePlot(obj_subset, region = "IL17A", extend.upstream = 5000, extend.downstream = 5000) & scale_fill_manual(values = paste0("#", colors_signal[,1]))
dev.off()
pdf(paste(outdir, "/IL2.pdf", sep = ""), width = 7, height = 4)
CoveragePlot(obj_subset, region = "IL2", extend.upstream = 5000, extend.downstream = 5000) & scale_fill_manual(values = paste0("#", colors_signal[,1]))
dev.off()
pdf(paste(outdir, "/FOXP3.pdf", sep = ""), width = 7, height = 4)
CoveragePlot(obj_subset, region = "FOXP3", extend.upstream = 500, extend.downstream = 500) & scale_fill_manual(values = paste0("#", colors_signal[,1]))
dev.off()
pdf(paste(outdir, "/GZMB.pdf", sep = ""), width = 7, height = 4)
CoveragePlot(obj_subset, region = "GZMB", extend.upstream = 500, extend.downstream = 500) & scale_fill_manual(values = paste0("#", colors_signal[,1]))
dev.off()
pdf(paste(outdir, "/CD58.pdf", sep = ""), width = 7, height = 4)
CoveragePlot(obj_subset, region = "CD58", extend.upstream = 10000, extend.downstream = 10000) & scale_fill_manual(values = paste0("#", colors_signal[,1]))
dev.off()
pdf(paste(outdir, "/CD2_Vln.pdf", sep = ""), width = 7, height = 4)
VlnPlot(obj_subset, features = "CD2", pt.size = 0)
dev.off()
pdf(paste(outdir, "/CD58_Vln.pdf", sep = ""), width = 7, height = 4)
VlnPlot(obj_subset, features = "CD58", pt.size = 0)
dev.off()
pdf(paste(outdir, "/CD48_Vln.pdf", sep = ""), width = 7, height = 4)
VlnPlot(obj_subset, features = "CD48", pt.size = 0)
dev.off()
pdf(paste(outdir, "/RUNX2_Vln.pdf", sep = ""), width = 7, height = 4)
VlnPlot(obj_subset, features = "RUNX2", pt.size = 0)
dev.off()
pdf(paste(outdir, "/IRF4_Vln.pdf", sep = ""), width = 7, height = 4)
VlnPlot(obj_subset, features = "IRF4", pt.size = 0)
dev.off()

########### UMAPs (Fig1B)
# Reassign names from the "resting/act" naming convention to just their major celltype names
obj$celltype2[which(obj$celltype2 == "Braking")] <- "TCM_act"
obj$celltype2[which(obj$celltype2 == "Naive_rest2")] <- "Naive_rest"
obj$celltype2[which(obj$celltype2 == "Naive_rest" | obj$celltype2 == "Naive_act")] <- "Naive"
obj$celltype2[which(obj$celltype2 == "TCM_rest" | obj$celltype2 == "TCM_act")] <- "TCM"
obj$celltype2[which(obj$celltype2 == "Th1_rest" | obj$celltype2 == "Th1_act")] <- "Th1"
obj$celltype2[which(obj$celltype2 == "Th2_rest" | obj$celltype2 == "Th2_act")] <- "Th2"
obj$celltype2[which(obj$celltype2 == "Th17_rest" | obj$celltype2 == "Th17_act")] <- "Th17"
Idents(obj) <- obj$celltype2
celltype_order <- c("Naive","TCM","Th1","Th2","Th17","TCM/TEM","TEM_act","TEM_act2","Treg")
levels(obj) <- celltype_order
#Plot weighted nearest neighbor umap
pdf(paste(outdir, "/UMAP.pdf", sep = ""), width = 8, height = 8)
DimPlot(obj, reduction = "wnn.umap", label = F, cols = paste0("#", colors[,1]), shuffle = T) + NoLegend()
dev.off()

######### Timepoint Resolved celltype frequency barplot (Fig S1K)
obj$celltype_condition <- paste(obj$celltype2, obj$Condition, sep = "_")
celltype_freq <- table(obj$celltype_condition)
celltype_condition <- c("Naive_Resting","Naive_2","Naive_5","Naive_15","TCM_Resting","TCM_2","TCM_5","TCM_15",
                        "Th1_Resting","Th1_2","Th1_5","Th1_15","Th2_Resting","Th2_2","Th2_5","Th2_15",
                        "Th17_Resting","Th17_2","Th17_5","Th17_15","TCM/TEM_Resting","TCM/TEM_2","TCM/TEM_5","TCM/TEM_15",
                        "TEM_act_Resting","TEM_act_2","TEM_act_5","TEM_act_15","TEM_act2_Resting","TEM_act2_2","TEM_act2_5","TEM_act2_15",
                        "Treg_Resting","Treg_2","Treg_5","Treg_15")
celltype_freq <- celltype_freq[celltype_condition]
condition <- rep(c("Resting","2","5","15"), 9)
alpha <- rep(c(0.7, 0.8, 0.9, 1), 9)
celltype <- c("Naive","Naive","Naive","Naive","TCM","TCM","TCM","TCM","Th1","Th1","Th1","Th1","Th2","Th2","Th2","Th2","Th17","Th17","Th17","Th17",
            "TCM/TEM","TCM/TEM","TCM/TEM","TCM/TEM","TEM_act","TEM_act","TEM_act","TEM_act","TEM_act2","TEM_act2","TEM_act2","TEM_act2",
            "Treg","Treg","Treg","Treg")
df <- data.frame(Freq = as.vector(celltype_freq), CellType = celltype, Condition = condition)
df$CellType <- factor(df$CellType, levels = unique(celltype))
df$Condition <- factor(df$Condition, levels = c("Resting","2","5","15"))
pdf(paste(outdir, "/CelltypeFreq_Condition.pdf", sep = ""), width = 12, height = 6)
ggplot(df, aes(x = CellType, y = Freq, fill = CellType, alpha = alpha)) + geom_bar(position = "stack", stat = "identity") +
    scale_fill_manual(values = paste0("#", colors[,1])) +
    theme(panel.grid = element_line(color = "grey",linewidth = 0.1)) + 
    theme(panel.background = element_rect("white")) +
    theme(axis.title=element_text(size=15, face='plain')) + 
    theme(axis.text=element_text(size=14, face='plain'))
dev.off()         

######### Fig S2A-F: RNA and ATAC Principal Component Analysis
minimal_genes <- read.table(paste0(indir, "/minimal_genes.txt"))$V1
minimal_peaks <- read.table(paste0(indir, "/minimal_peaks.bed"))
minimal_peaks <- paste(minimal_peaks$V1, minimal_peaks$V2, minimal_peaks$V3, sep = "-")
RNA_counts <- RNA_counts[minimal_genes, ]
ATAC_counts <- ATAC_counts[minimal_peaks, ]
RNA_counts <- t(scale(t(RNA_counts)))
ATAC_counts <- t(scale(t(ATAC_counts)))

## Counts Processing
celltypes <- c("Naive_rest","TCM_rest","Th1_rest","Th2_rest","Th17_rest","Naive_act",
                "TCM_act","TCM/TEM","Th1_act","Th2_act","Th17_act","TEM_act",
                "TEM_act2","Treg")
celltype_timepoint <- c("Naive_rest_Resting","TCM_rest_Resting","Th1_rest_Resting","Th2_rest_Resting","Th17_rest_Resting", "TEM_act_Resting","TEM_act2_Resting",
                        "Naive_act_2","Naive_act_5","Naive_act_15","TCM_act_2","TCM_act_5","TCM_act_15","TCM/TEM_2","TCM/TEM_5","TCM/TEM_15",
                        "Th1_act_2","Th1_act_5","Th1_act_15","Th2_act_2","Th2_act_5","Th2_act_15","Th17_act_2","Th17_act_5","Th17_act_15",
                        "TEM_act_2","TEM_act_5","TEM_act_15","TEM_act2_2","TEM_act2_5","TEM_act2_15",
                        "Treg_Resting","Treg_2","Treg_5","Treg_15")
subset_cells <- which(meta$celltype_condition %in% celltype_timepoint) ## Cell types to subset NULL if all)

## Subset counts if needed
if(!is.null(subset_cells) ){
    RNA_counts <- RNA_counts[,subset_cells]
    ATAC_counts <- ATAC_counts[,subset_cells]
    meta <- meta[subset_cells,]
}

celltype_major <- read.table(paste0(indir, "/celltype_major.txt"))$V1
meta$celltype <- celltype_major
meta$celltype_condition <- paste(meta$celltype, meta$Condition, sep = "_")
rownames(meta) <- paste(meta$celltype, meta$Condition, meta$Donor, sep = "_")
colnames(RNA_counts) <- rownames(meta)
colnames(ATAC_counts) <- rownames(meta)
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
    RNA_counts <- RNA_counts[,subset_cells]
    ATAC_counts <- ATAC_counts[,subset_cells]
    meta <- meta[subset_cells,]
}
meta$celltype[which(meta$celltype == "TEM")] <- "MHCII"
meta$celltype[which(meta$celltype == "TEM2")] <- "CTL"
order_celltype <- c("Naive","TCM","Th1","Th2","Th17","TCM/TEM","MHCII","CTL","Treg")

# RNA calculate PC1 and PC2
pca <- prcomp(t(RNA_counts))
var_exp <- pca$sdev^2
var_exp <- var_exp/sum(var_exp)
var_exp1 <- signif(100*var_exp[1],digits=3)
var_exp2 <- signif(100*var_exp[2],digits=3)
var_exp3 <- signif(100*var_exp[3],digits=3)
var_exp4 <- signif(100*var_exp[4],digits=3)
df_pca <- cbind(meta, pca$x)
df_pca$Condition <- factor(df_pca$Condition, levels=c("Resting","2","5", "15"), ordered=T)
df_pca$Donor <- factor(as.character(df_pca$Donor), levels = c("Donor1","Donor2","Donor3","Donor4"))
df_pca$celltype <- factor(df_pca$celltype, levels=order_celltype)

# PC1 vs PC2
ggplot(df_pca) + geom_point(aes(x=PC1, y=PC2, color=celltype, shape = Condition, size = Condition)) + 
    labs(x=paste0('PC1 [',var_exp1,'%]'),y=paste0('PC2 [',var_exp2,'%]'), color='') +
    theme(panel.grid = element_line(color = "grey",linewidth = 0.1)) + 
    theme(panel.background = element_rect("black", fill = NA)) +
    theme(axis.title=element_text(size=15, face='plain')) + 
    theme(axis.text=element_text(size=14, face='plain')) + 
    geom_hline(yintercept = 0, linetype='dashed')+
    geom_vline(xintercept = 0, linetype='dashed') + 
    scale_shape_manual(values = c(15,16,17,18)) + 
    scale_color_manual(values = c("#E41A1C","#4A72A6","#48A462","#7E6E85","#D16948","#E1C62F","#B75F49","#EC83BA","#999999")) + 
    scale_size_manual(values = c(2,3,3,3)) + 
    stat_ellipse(aes(x=PC1, y=PC2, group = Condition))
file_out <- file.path(paste0(outdir, "/RNA_PC1_2.pdf"))
ggsave(file_out, height=5.5, width=7.3, dpi = 600)


# PC3 vs PC4
ggplot(df_pca) + geom_point(aes(x=PC3, y=PC4, color=celltype, shape = Condition, size = Condition)) + 
    labs(x=paste0('PC3 [',var_exp3,'%]'),y=paste0('PC4 [',var_exp4,'%]'), color='') +
    theme(panel.grid = element_line(color = "grey",linewidth = 0.1)) + 
    theme(panel.background = element_rect("black", fill = NA)) +
    theme(axis.title=element_text(size=15, face='plain')) + 
    theme(axis.text=element_text(size=14, face='plain')) + 
    geom_hline(yintercept = 0, linetype='dashed')+
    geom_vline(xintercept = 0, linetype='dashed') + 
    scale_shape_manual(values = c(15,16,17,18)) + 
    scale_color_manual(values = c("#E41A1C","#4A72A6","#48A462","#7E6E85","#D16948","#E1C62F","#B75F49","#EC83BA","#999999")) + 
    scale_size_manual(values = c(2,3,3,3)) + 
    stat_ellipse(aes(x=PC3, y=PC4, group = Condition))
file_out <- file.path(paste0(outdir, "/RNA_PC3_4.pdf"))
ggsave(file_out, height=5.5, width=7.3)

# RNA Loadings PC1 and PC2
rotation <- pca$rotation[,1:2]
label <- rep("No", length(rotation[,1]))
index <- which(rownames(RNA_counts) %in% c("CD44","CD69","LEF1","IL2RA","IFNG","IL4","IL17A","FOS","JUN","TBX21","RORC","IRF4","TCF7","IL10","ORC6","ORC1","HLA-C","FHIT","CCR7","RUNX1"))
label[index] <- "Yes"
df <- data.frame(PC1 = rotation[,1], PC2 = rotation[,2], Label = label, Gene = rownames(RNA_counts))
ggplot(df, aes(x = PC1, y = PC2, alpha = Label, size = Label)) + geom_point(color = "grey30") +
    scale_alpha_manual(values = c(0.05, 0.9)) + scale_size_manual(values = c(1,1)) +
      geom_text_repel(data=df[index,],aes(PC1,PC2,label=Gene,),size = 6, fontface = "italic") +
      theme(panel.grid = element_line(color = "grey",linewidth = 0.1)) + 
    theme(panel.background = element_rect("black", fill = NA)) +
    theme(axis.title=element_text(size=15, face='plain')) + 
    theme(axis.text=element_text(size=14, face='plain')) +
    geom_hline(yintercept = 0) + 
    geom_vline(xintercept = 0) +
    theme(axis.ticks.x=element_blank()) +
    theme(axis.ticks.y=element_blank()) +
    NoLegend()
file_out <- file.path(paste0(outdir, "/RNA_Loadings_PC1_2.pdf"))
ggsave(file_out, height=6.5, width=8, dpi = 600)

## RNA Loadings PC3 and PC4
rotation <- pca$rotation[,3:4]
label <- rep("No", length(rotation[,1]))
index <- which(rownames(RNA_counts) %in% c("CD44","CD69","LEF1","IL2RA","IFNG","IL4","IL17A","FOS","JUN","TBX21","RORC","IRF4","TCF7","IL10","ORC6","ORC1","HLA-C","FHIT","CCR7","RUNX1"))
label[index] <- "Yes"
df <- data.frame(PC3 = rotation[,1], PC4 = rotation[,2], Label = label, Gene = rownames(RNA_counts))
ggplot(df, aes(x = PC3, y = PC4, alpha = Label, size = Label)) + geom_point(color = "grey30") +
    scale_alpha_manual(values = c(0.05, 0.9)) + scale_size_manual(values = c(1,1)) +
      geom_text_repel(data=df[index,],aes(PC3,PC4,label=Gene,),size = 6, fontface = "italic") +
      theme(panel.grid = element_line(color = "grey",linewidth = 0.1)) + 
    theme(panel.background = element_rect("black", fill = NA)) +
    theme(axis.title=element_text(size=15, face='plain')) + 
    theme(axis.text=element_text(size=14, face='plain')) +
    geom_hline(yintercept = 0) + 
    geom_vline(xintercept = 0) +
    theme(axis.ticks.x=element_blank()) +
    theme(axis.ticks.y=element_blank())
file_out <- file.path(paste0(outdir, "/RNA_Loadings_PC3_4.pdf"))
ggsave(file_out, height=6.5, width=8, dpi = 600)

# ATAC calculate PC1 and PC2
pca <- prcomp(t(ATAC_counts))
var_exp <- pca$sdev^2
var_exp <- var_exp/sum(var_exp)
var_exp1 <- signif(100*var_exp[1],digits=3)
var_exp2 <- signif(100*var_exp[2],digits=3)
var_exp3 <- signif(100*var_exp[3],digits=3)
var_exp4 <- signif(100*var_exp[4],digits=3)
df_pca <- cbind(meta, pca$x)
df_pca$Condition <- factor(df_pca$Condition, levels=c("Resting","2","5", "15"), ordered=T)
df_pca$Donor <- factor(as.character(df_pca$Donor), levels = c("Donor1","Donor2","Donor3","Donor4"))
df_pca$celltype <- factor(df_pca$celltype, levels=order_celltype)

# PC1 vs PC2
ggplot(df_pca) + geom_point(aes(x=PC1, y=PC2, color=celltype, shape = Condition, size = Condition)) + 
    labs(x=paste0('PC1 [',var_exp1,'%]'),y=paste0('PC2 [',var_exp2,'%]'), color='') +
    theme(panel.grid = element_line(color = "grey",linewidth = 0.1)) + 
    theme(panel.background = element_rect("black", fill = NA)) +
    theme(axis.title=element_text(size=15, face='plain')) + 
    theme(axis.text=element_text(size=14, face='plain')) + 
    geom_hline(yintercept = 0, linetype='dashed')+
    geom_vline(xintercept = 0, linetype='dashed') + 
    scale_shape_manual(values = c(15,16,17,18)) + 
    scale_color_manual(values = c("#E41A1C","#4A72A6","#48A462","#7E6E85","#D16948","#E1C62F","#B75F49","#EC83BA","#999999")) + 
    scale_size_manual(values = c(2,3,3,3)) + 
    stat_ellipse(aes(x=PC1, y=PC2, group = Condition))
file_out <- file.path(paste0(outdir, "/ATAC_PC1_2.pdf"))
ggsave(file_out, height=5.5, width=7.3, dpi = 600)


# PC3 vs PC4
ggplot(df_pca) + geom_point(aes(x=PC3, y=PC4, color=celltype, shape = Condition, size = Condition)) + 
    labs(x=paste0('PC3 [',var_exp3,'%]'),y=paste0('PC4 [',var_exp4,'%]'), color='') +
    theme(panel.grid = element_line(color = "grey",linewidth = 0.1)) + 
    theme(panel.background = element_rect("black", fill = NA)) +
    theme(axis.title=element_text(size=15, face='plain')) + 
    theme(axis.text=element_text(size=14, face='plain')) + 
    geom_hline(yintercept = 0, linetype='dashed')+
    geom_vline(xintercept = 0, linetype='dashed') + 
    scale_shape_manual(values = c(15,16,17,18)) + 
    scale_color_manual(values = c("#E41A1C","#4A72A6","#48A462","#7E6E85","#D16948","#E1C62F","#B75F49","#EC83BA","#999999")) + 
    scale_size_manual(values = c(2,3,3,3)) + 
    stat_ellipse(aes(x=PC3, y=PC4, group = Condition))
file_out <- file.path(paste0(outdir, "/ATAC_PC3_4.pdf"))
ggsave(file_out, height=5.5, width=7.3)

