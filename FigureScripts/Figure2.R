library(reshape2)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(matrixStats)
library(dplyr)
library(dendextend)
library(plotrix)
library(matrixStats)
library(RColorBrewer)
source("PlotPatterns.R")

outdir <- "output/Figure2"
indir <- "input"

##### Figure 2A: Rapid Recall (RR) gene heatmap
# Read in normalized counts and meta data
counts <- read.table(paste0(indir, "/RNA_counts_combatseq_vst.txt"))
meta <- read.table(paste0(indir, '/meta_filter.txt'))

# Custom ordering of celltypes and timepoints. 
celltypes <- c("Naive_rest","TCM_rest","Th1_rest","Th2_rest","Th17_rest","Naive_act",
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
heat_col_tfmrna <- colorRamp2(c(-2, 0, 2), c('dodgerblue','white','red')) ## Colors for heatmap
subset_cells <- which(meta$celltype_condition %in% celltype_timepoint) ## Cell types to subset NULL if all)
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

# Clustering solution with annotation column for which genes to label on heatmap
order_matrix <- read.table(paste0(indir,"/RNA_RR_clustering.txt"), header = T)
order <- c("C1","C2","C3","C4","C5","C6","C7","C8","C9")
genes <- order_matrix[,1]
label_genes <- order_matrix[which(order_matrix[,3]!=0),1]

label <- rep(0, length(genes))
label[which(order_matrix[,1] %in% label_genes & order_matrix[,2] == "C1")] <- 1
label[which(order_matrix[,1] %in% label_genes & order_matrix[,2] == "C2")] <- 2
label[which(order_matrix[,1] %in% label_genes & order_matrix[,2] == "C3")] <- 3
label[which(order_matrix[,1] %in% label_genes & order_matrix[,2] == "C4")] <- 4
label[which(order_matrix[,1] %in% label_genes & order_matrix[,2] == "C5")] <- 5
label[which(order_matrix[,1] %in% label_genes & order_matrix[,2] == "C6")] <- 6
label[which(order_matrix[,1] %in% label_genes & order_matrix[,2] == "C7")] <- 7
label[which(order_matrix[,1] %in% label_genes & order_matrix[,2] == "C8")] <- 8
label[which(order_matrix[,1] %in% label_genes & order_matrix[,2] == "C9")] <- 9
order_matrix <- cbind(order_matrix, label)

index <- which(order_matrix[,3] != 0)
ha2 = rowAnnotation(foo = anno_mark(at = c(index), labels = order_matrix[index, 1], labels_gp = gpar(fontsize = 5)))

rowAnn_across_left <-  HeatmapAnnotation(foo = anno_mark(at = index, labels = genes[index], side = 'left',labels_gp = gpar(fontface = 4, fontsize = 8, col = ifelse(order_matrix[index,3] == "1", '#E41A1C',
                                                                                                                                                                        ifelse(order_matrix[index,3] == "2", 'orange',
                                                                                                                                                                               ifelse(order_matrix[index,3] == "3", 'gold',
                                                                                                                                                                                      ifelse(order_matrix[index,3] == "4", 'darkgreen',
                                                                                                                                                                                             ifelse(order_matrix[index,3] == "5", 'darkred',
                                                                                                                                                                                                    ifelse(order_matrix[index,3] == "6", 'lightblue',
                                                                                                                                                                                                           ifelse(order_matrix[index,3] == "7", '#6699CC',
                                                                                                                                                                                                                ifelse(order_matrix[index,3] == "8", 'purple', '#5BC236')))))))))),
                                    which = 'row',
                                    show_annotation_name = F)

## Create heatmap
pdf(paste0(outdir, "/RNA_RR_heatmap.pdf"), width = 10, height = 10)
ht <- Heatmap(counts[genes, ], name = "Z-score", show_row_names = F,  show_row_dend = F, show_column_dend = F, 
              row_split = factor(order_matrix[,2], levels = order), border = T, 
              show_column_names = F, column_split = factor(annotation_matrix[,2], levels = celltype_timepoint),
              cluster_row_slices = F, col = heat_col_tfmrna, column_gap = splits, 
              row_gap = unit(0, "mm"), left_annotation = rowAnn_across_left,
              cluster_columns = F , top_annotation = ha, row_names_gp = (gpar(fontsize = 8)), use_raster = T,
              column_title = NULL, row_title = NULL)
draw(ht)
dev.off()


########## Figure S2J: Dynamic gene (Stim) heatmap
dynamic_genes <- read.table(paste0(indir, "/GeneDynamics_df_Filtered.txt"), header = T)
dynamic_genes <- unique(dynamic_genes$Gene)
dynamic_genes <- dynamic_genes[which(!(dynamic_genes %in% genes))]
cluster_df <- read.table(paste0(indir,"/RNA_Stim_clustering.txt"), header = T)

# Create heatmap
pdf(paste0(outdir, "/RNA_Stim_heatmap.pdf"), width = 10, height = 10)
ht <- Heatmap(counts[cluster_df[,1], ], name = "Z-score", show_row_names = F,  show_row_dend = F, show_column_dend = F, 
              row_split = factor(cluster_df[,2], levels = c("C1","C2","C3","C4","C5")), border = T, 
              show_column_names = F, column_split = factor(annotation_matrix[,2], levels = celltype_timepoint),
              cluster_row_slices = F, col = heat_col_tfmrna, column_gap = splits, 
              row_gap = unit(0, "mm"), 
              cluster_columns = F , top_annotation = ha, row_names_gp = (gpar(fontsize = 14)), use_raster = T,
              column_title = NULL, row_title = NULL)
draw(ht)
dev.off()

####### Figure S2A: Heatmap of celltype specific RR genes (up)
genes <- read.table(paste0(indir, "/Upgenes_in_unique_celltypes.txt"))$V1
pdf(paste0(outdir, "/RR_celltype_specific.pdf"), width = 10, height = 10)
ht <- Heatmap(counts[genes, ], name = "Z-score", show_row_names = F,  show_row_dend = F, show_column_dend = F, 
              #row_split = factor(cluster_df[,2], levels = c("C1","C2","C3","C4","C5")), border = T, 
              show_column_names = F, column_split = factor(annotation_matrix[,2], levels = celltype_timepoint),
              cluster_row_slices = F, col = heat_col_tfmrna, column_gap = splits, 
              row_gap = unit(0, "mm"), cluster_rows = F,
              cluster_columns = F , top_annotation = ha, row_names_gp = (gpar(fontsize = 14)), use_raster = T,
              column_title = NULL, row_title = NULL)
draw(ht)
dev.off()

####### Figure S2B: Venn Diagrams of celltype specific genes
RR <- read.table(paste0(indir, "/RR_df_Filtered.txt"), header = T)
RR[which(RR$Celltype == "TEM"),2] <- "MHCII"
RR[which(RR$Celltype == "TEM2"),2] <- "CTL"
RR_up <- RR[which(RR$Direction == "Up"),]
RR_down <- RR[which(RR$Direction == "Down"),]

RR_up_onecelltype <- RR_up %>%
  group_by(Gene) %>%
  mutate(num_celltypes = n_distinct(Celltype)) %>%
  ungroup() %>%
  filter(num_celltypes == 1)

RR_up_unique_pairs <- RR_up_onecelltype %>%
  distinct(Gene, Celltype) %>%
  arrange(Gene)

RR_up_unique_pairs <- noquote(as.matrix(RR_up_unique_pairs))

RR_down_onecelltype <- RR_down %>%
  group_by(Gene) %>%
  mutate(num_celltypes = n_distinct(Celltype)) %>%
  ungroup() %>%
  filter(num_celltypes == 1)

RR_down_unique_pairs <- RR_down_onecelltype %>%
  distinct(Gene, Celltype) %>%
  arrange(Gene)

RR_down_unique_pairs <- noquote(as.matrix(RR_down_unique_pairs))

nonunique_up <- length(unique(RR_up$Gene)) - length(RR_up_onecelltype$Gene)
nonunique_down <- length(unique(RR_down$Gene)) - length(RR_down_onecelltype$Gene)


up_unique_freq <- table(RR_up_unique_pairs[,2])
down_unique_freq <- table(RR_down_unique_pairs[,2])
up_unique_freq <- setNames(c(up_unique_freq, nonunique_up), c(names(up_unique_freq), "NotSpecific"))
down_unique_freq <- setNames(c(down_unique_freq, nonunique_down), c(names(down_unique_freq), "NotSpecific"))
order <- c("TCM","Th1","Th2","Th17","MHCII","CTL","Treg","NotSpecific")
colors <- c("#59709E","#70A16D","#B97256","#797081","#D28CB6","#A26652","#999999","#555555")

pdf(paste0(outdir, "/UpUniquePie.pdf"), width = 6, height = 6)
pie(as.vector(up_unique_freq[order]), labels = order, border=NA, col=colors)
dev.off()

pdf(paste0(outdir, "/DownUniquePie.pdf"), width = 6, height = 6)
pie(as.vector(down_unique_freq[order]), labels = order, border=NA, col=colors)
dev.off()

######## Figure S3H-I: Plotting interesting receptor/ligand genes
Receptors <- read.table(paste0(indir, "/receptor_clustering.txt"))
Ligands <- read.table(paste0(indir, "/ligand_clustering.txt"))
gene_dynamics <- read.table(paste0(indir, "/GeneDynamics_df_Filtered.txt"), header = T)
Ligands <- Ligands[which(Ligands[,1] %in% gene_dynamics$Gene),]
Receptors <- Receptors[which(Receptors[,1] %in% gene_dynamics$Gene),]

## Create heatmap
pdf(paste0(outdir, "/Receptors.pdf"), width = 11, height = 22)
ht <- Heatmap(counts[Receptors[,1], ], name = "Z-score", show_row_names = T,  show_row_dend = F, show_column_dend = F,
              row_split = factor(Receptors[,2], levels = c("group1","group2","group3","group4")), border = T, 
              show_column_names = F, column_split = factor(annotation_matrix[,2], levels = celltype_timepoint),
              cluster_row_slices = F, col = heat_col_tfmrna, column_gap = splits, 
              row_gap = unit(0, "mm"),
              cluster_columns = F , top_annotation = ha, row_names_gp = (gpar(fontsize = 8)), use_raster = T,
              column_title = NULL, row_title = NULL)
draw(ht)
dev.off()

pdf(paste0(outdir, "/Ligands.pdf"), width = 11, height = 22)
ht <- Heatmap(counts[Ligands[,1],], name = "Z-score", show_row_names = T,  show_row_dend = F, show_column_dend = F,
              row_split = factor(Ligands[,2], levels = c("group1","group2","group3","group4")), border = T, 
              show_column_names = F, column_split = factor(annotation_matrix[,2], levels = celltype_timepoint),
              cluster_row_slices = F, col = heat_col_tfmrna, column_gap = splits, 
              row_gap = unit(0, "mm"),
              cluster_columns = F , top_annotation = ha, row_names_gp = (gpar(fontsize = 8)), use_raster = T,
              column_title = NULL, row_title = NULL)
draw(ht)
dev.off()

######## Line plot patterns for figure 2A
outs <- c(paste0(outdir, "/RR_Cluster1.pdf"),
    paste0(outdir, "/RR_Cluster2.pdf"),
    paste0(outdir, "/RR_Cluster3.pdf"),
    paste0(outdir, "/RR_Cluster4.pdf"),
    paste0(outdir, "/RR_Cluster5.pdf"), 
    paste0(outdir, "/RR_Cluster6.pdf"),
    paste0(outdir, "/RR_Cluster7.pdf"),
    paste0(outdir, "/RR_Cluster8.pdf"),
    paste0(outdir, "/RR_Cluster9.pdf"))

clusters <- read.table(paste0(indir,"/RNA_RR_clustering.txt"))
Cluster1_RR <- clusters[which(clusters[,2] == "C1"),1]
Cluster2_RR <- clusters[which(clusters[,2] == "C2"),1]
Cluster3_RR <- clusters[which(clusters[,2] == "C3"),1]
Cluster4_RR <- clusters[which(clusters[,2] == "C4"),1]
Cluster5_RR <- clusters[which(clusters[,2] == "C5"),1]
Cluster6_RR <- clusters[which(clusters[,2] == "C6"),1]
Cluster7_RR <- clusters[which(clusters[,2] == "C7"),1]
Cluster8_RR <- clusters[which(clusters[,2] == "C8"),1]
Cluster9_RR <- clusters[which(clusters[,2] == "C9"),1]


RR_Clusters <- list(Cluster1_RR, Cluster2_RR, Cluster3_RR, Cluster4_RR, Cluster5_RR, Cluster6_RR, Cluster7_RR, Cluster8_RR, Cluster9_RR)
for(i in 1:length(RR_Clusters)){
    print(i)
    plotPattern(RR_Clusters[[i]], meta, counts, outs[i])
}

outs <- c(paste0(outdir, "/Stim_Cluster1.pdf"),
        paste0(outdir, "/Stim_Cluster2.pdf"),
        paste0(outdir, "/Stim_Cluster3.pdf"),
        paste0(outdir, "/Stim_Cluster4.pdf"),
        paste0(outdir, "/Stim_Cluster5.pdf"))

clusters <- read.table(paste0(indir,"/RNA_Stim_clustering.txt"), header = T)
Cluster1_Naive <- clusters[which(clusters[,2] == "C1"),1]
Cluster2_Naive <- clusters[which(clusters[,2] == "C2"),1]
Cluster3_Naive <- clusters[which(clusters[,2] == "C3"),1]
Cluster4_Naive <- clusters[which(clusters[,2] == "C4"),1]
Cluster5_Naive <- clusters[which(clusters[,2] == "C5"),1]

Naive_clusters <- list(Cluster1_Naive, Cluster2_Naive, Cluster3_Naive, Cluster4_Naive)
for(i in 1:length(Naive_clusters)){
    plotPattern(Naive_clusters[[i]], meta, counts, outs[i])
}

### Figure S2F-G
TFs <- c(
  "IRF8", "KLF7", "STAT2", "ZNF571", "TCF7", "ELF2", "ZNF563", "YY1", "PATZ1", "FOXO1",
  "PLAG1", "FOS", "ZNF449", "ZNF283", "ETS1", "FLI1", "KLF2", "E2F3", "ZNF184", "MYBL1",
  "SP1", "MAZ", "SRF", "IRF4", "BATF", "STAT1", "BCL6", "ZNF263", "JUN", "SP2",
  "RUNX2", "RORA", "RORC", "RBPJ", "MAF", "PRDM1", "ZNF212", "STAT4", "TBX21", "ZBTB17",
  "ATF4", "NFAT5", "POU2F2", "STAT3", "NFKB2", "REL", "IRF1", "KLF9", "NFATC1", "RELA",
  "THAP4", "EGR1", "RARA", "ATF1", "NFKB1", "CDC5L", "TCF3", "EGR2", "FOXP1", "MAFF",
  "LEF1", "ELK4", "GATA3", "TCF12", "IRF2", "THAP1", "ZBTB14", "ZNF717", "EGR3", "TSHZ3",
  "POU2F1", "DNMT1", "ZNF124", "RUNX1", "FOXJ3", "KLF3", "SMAD3", "KLF6", "ZNF37A", "FOXK2",
  "KLF10", "RUNX3", "XBP1", "ATF3", "RFX5", "JUNB", "E2F4", "YBX3", "KLF11","BACH2","GLIS3",
  "VDR","NFE2L1","PLAGL1","ATF5"
)
counts_RNA <- read.table(paste0(indir, "/RNA_counts_combatseq_vst.txt"))
counts_TFA <- read.table(paste0(indir, "/TFA.txt"))
meta <- read.table(paste0(indir, "/meta_filter.txt"))

# Custom ordering of celltypes and timepoints. 
celltypes <- c("Naive_rest","TCM_rest","Th1_rest","Th2_rest","Th17_rest","Naive_act",
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
heat_col_tfmrna <- colorRamp2(c(-2, 0, 2), c('dodgerblue','white','red')) ## Colors for heatmap
subset_cells <- which(meta$celltype_condition %in% celltype_timepoint) ## Cell types to subset NULL if all)
if(!is.null(subset_cells) ){
    counts_RNA <- counts_RNA[,subset_cells]
    counts_TFA <- counts_TFA[,subset_cells]
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
    counts_RNA <- counts_RNA[,subset_cells]
    counts_TFA <- counts_TFA[,subset_cells]
    meta <- meta[subset_cells,]
}
## Make sure col names are correct. In my colnames, there are some dots so i need to convert them to /
counts_RNA <- as.matrix(counts_RNA)
counts_TFA <- as.matrix(counts_TFA)
colnames(counts) <- gsub("\\.", "/", colnames(counts))
colnames(counts_TFA) <- gsub("\\.","/", colnames(counts_TFA))
## Zscore counts
counts_RNA <- t(scale(t(counts_RNA)))
counts_TFA <- t(scale(t(counts_TFA)))

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
order_df <- read.table(paste0(indir, "/TFmRNA_order_df.txt"), header = T)

pdf(paste0(outdir, "/TFmRNA.pdf"), width = 7, height = 10)
ht <- Heatmap(counts_RNA[cluster_df[,1], ], name = "Z-score", show_row_names = T,  show_row_dend = F, show_column_dend = F, 
              border = T, 
              show_column_names = F, column_split = factor(annotation_matrix[,2], levels = celltype_timepoint),
              cluster_row_slices = F, col = heat_col_tfmrna, column_gap = splits, 
              row_gap = unit(0, "mm"), row_split = factor(cluster_df[,2], levels = c("group8", "group7","group4","group3","group5","group6","group2","group1")),
              cluster_columns = F , top_annotation = ha, row_names_gp = (gpar(fontsize = 8)), use_raster = T,
              column_title = NULL, row_title = NULL, row_names_side = "left")
ht <- draw(ht)
dev.off()
pdf(paste0(outdir, "/TFA.pdf"), width = 7, height = 10)
ht <- Heatmap(counts_TFA[cluster_df[,1], ], name = "Z-score", show_row_names = T,  show_row_dend = F, show_column_dend = F, 
              border = T, 
              show_column_names = F, column_split = factor(annotation_matrix[,2], levels = celltype_timepoint),
              cluster_row_slices = F, col = heat_col_tfmrna, column_gap = splits, 
              row_gap = unit(0, "mm"), row_split = factor(cluster_df[,2],  levels = c("group8", "group7","group4","group3","group5","group6","group2","group1")),
              cluster_columns = F , top_annotation = ha, row_names_gp = (gpar(fontsize = 8)), use_raster = T,
              column_title = NULL, row_title = NULL, row_names_side = "left")
ht <- draw(ht)
dev.off()

