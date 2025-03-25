library(viridis)
library(Polychrome)
library(ggplot2)
library(Seurat)
library(Signac)
library(dplyr)
library(tidyverse)
library(reshape2)
library(reshape)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
source("ORA.R")
source("PlotGene.R")
source("PlotPatternsTCM.R")

outdir <- "Figure3"
indir <- "input"

# Load in seurat object and create a celltype_condition metadata column
obj <- readRDS(paste0(indir, "/integrated2.rds"))
obj$celltype_condition <- paste(obj$celltype2, obj$Condition, sep = "_")

# Load in dataframe of cellcycle phase scores
df <- read.table(paste0(indir, "/CellCycle_df.txt"), header = T)
# Remove Naive and Treg from dataset
index <- which(obj$celltype2 %in% c("TCM_act","TCM/TEM","Braking","Th1_act","Th2_act","Th17_act","TEM_act","TEM_act2") & obj$Condition == "15")
obj_subset <- subset(obj, cells = index)
Idents(obj_subset) <- obj_subset$celltype2
levels(obj_subset) <- c("TCM_act","TCM/TEM","Braking","Th1_act","Th2_act","Th17_act","TEM_act","TEM_act2")

########## Figure 3: RNA, ATAC UMAPs of activated memory T-cells
DefaultAssay(obj_subset) <- "RNA"
obj_subset <- RunUMAP(obj_subset, reduction = "harmony_rna", dims = 1:30, reduction.name = "UMAP_RNA")
pdf(paste0(outdir, "/RNA_activated_memory_UMAP.pdf"), width = 5, height = 4)
DimPlot(obj_subset, reduction = "UMAP_RNA", cols = c("#4A72A6","#E1C62F","#8B0000","#48A462","#7E6E85","#D16948","#B75F49","#EC83BA"))
dev.off()

DefaultAssay(obj_subset) <- "peaks"
obj_subset <- RunUMAP(obj_subset, reduction = "harmony_atac", dims = 2:30, reduction.name = "UMAP_ATAC")
pdf(paste0(outdir, "/ATAC_activated_memory_UMAP.pdf"), width = 5, height = 4)
DimPlot(obj_subset, reduction = "UMAP_ATAC", cols = c("#4A72A6","#E1C62F","#8B0000","#48A462","#7E6E85","#D16948","#B75F49","#EC83BA"))
dev.off()

## RNA Density of proliferating TCM
umap_coords <- Embeddings(obj_subset, reduction = "UMAP_RNA")
umap_data <- data.frame(UMAP_1 = umap_coords[,1], UMAP_2 = umap_coords[,2], cluster = Idents(obj_subset))
pdf(paste0(outdir, "/RNA_proliferation_density.pdf"), width = 5, height = 4)
ggplot(umap_data, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = cluster == "Braking"), alpha = 0.3) +
  #stat_density_2d(aes(fill = ..level..), geom = "polygon", contour = TRUE) +
  scale_color_manual(values = c("grey", "#8B0000")) + 
  labs(title = "Density of Specific Cluster on UMAP") +
  theme_minimal() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
dev.off()

## ATAC Density of proliferating TCM
umap_coords <- Embeddings(obj_subset, reduction = "UMAP_ATAC")
umap_data <- data.frame(UMAP_1 = umap_coords[,1], UMAP_2 = umap_coords[,2], cluster = Idents(obj_subset))
pdf(paste0(outdir, "/ATAC_proliferation_density.pdf"), width = 5, height = 4)
ggplot(umap_data, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = cluster == "Braking"), alpha = 0.3) +
  #stat_density_2d(aes(fill = ..level..), geom = "polygon", contour = TRUE) +
  scale_color_manual(values = c("grey", "#8B0000")) + 
  labs(title = "Density of Specific Cluster on UMAP") +
  theme_minimal() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
dev.off()

############# Figure S4A: Proliferating TCM marker gene dotplot
obj$celltype2[which(obj$celltype2 == "Naive_rest2")] <- "Naive_rest"
obj$celltype2[which(obj$celltype2 == "Naive_rest" | obj$celltype2 == "Naive_act")] <- "Naive"
obj$celltype2[which(obj$celltype2 == "TCM_rest" | obj$celltype2 == "TCM_act")] <- "TCM"
obj$celltype2[which(obj$celltype2 == "Th1_rest" | obj$celltype2 == "Th1_act")] <- "Th1"
obj$celltype2[which(obj$celltype2 == "Th2_rest" | obj$celltype2 == "Th2_act")] <- "Th2"
obj$celltype2[which(obj$celltype2 == "Th17_rest" | obj$celltype2 == "Th17_act")] <- "Th17"
Idents(obj) <- paste(obj$celltype2, obj$Condition, sep = "_")
celltype_condition <- c("Naive_Resting","Naive_2","Naive_5","Naive_15","TCM_Resting","TCM_2","TCM_5","TCM_15","Braking_15",
                        "Th1_Resting","Th1_2","Th1_5","Th1_15","Th2_Resting","Th2_2","Th2_5","Th2_15",
                        "Th17_Resting","Th17_2","Th17_5","Th17_15","TCM/TEM_2","TCM/TEM_5","TCM/TEM_15",
                        "TEM_act_Resting","TEM_act_2","TEM_act_5","TEM_act_15","TEM_act2_Resting","TEM_act2_2","TEM_act2_5","TEM_act2_15",
                        "Treg_Resting","Treg_2","Treg_5","Treg_15")
index <- which(Idents(obj) %in% celltype_condition)
obj_subset <- subset(obj, cells = index)
levels(obj_subset) <- celltype_condition
DefaultAssay(obj_subset) <- "RNA"
genes <- c("MCM10","POLQ","HELLS","CHEK1","DNA2")
p <- DotPlot(object = obj_subset, features = genes)
df <- p$data
exp_mat <- df %>% 
  select(-pct.exp, -avg.exp) %>%  
  pivot_wider(names_from = id, values_from = avg.exp.scaled) %>% 
  as.data.frame() 
row.names(exp_mat) <- exp_mat$features.plot  
exp_mat <- exp_mat[,-1] %>% as.matrix()
percent_mat<-df %>% 
  select(-avg.exp, -avg.exp.scaled) %>%  
  pivot_wider(names_from = id, values_from = pct.exp) %>% 
  as.data.frame() 
row.names(percent_mat) <- percent_mat$features.plot  
percent_mat <- percent_mat[,-1] %>% as.matrix()
col_fun = circlize::colorRamp2(c(-2, 0, 3), c("dodgerblue", "grey", "red"))
cell_fun = function(j, i, x, y, w, h, fill){
          grid.rect(x = x, y = y, width = w, height = h, 
                    gp = gpar(col = NA, fill = NA))
          grid.circle(x=x,y=y,r= percent_mat[i, j]/50 * min(unit.c(w, h)),
                      gp = gpar(fill = col_fun(exp_mat[i, j]), col = NA))}
cluster_anno <- c("Naive","Naive","Naive","Naive","TCM","TCM","TCM","TCM","Braking","Th1","Th1","Th1","Th1","Th2","Th2","Th2","Th2",
                "Th17","Th17","Th17","Th17","TCM/TEM","TCM/TEM","TCM/TEM","TEM","TEM","TEM","TEM","TEM2","TEM2","TEM2","TEM2","Treg","Treg","Treg","Treg")
stim <- c("Resting","2","5","15","Resting","2","5","15","15","Resting","2","5","15","Resting","2","5","15","Resting","2","5","15",
        "2","5","15","Resting","2","5","15","Resting","2","5","15","Resting","2","5","15")
heatmap_colors <- c("#E41A1C","#4A72A6","#8B0000","#48A462","#7E6E85","#D16948","#E1C62F","#B75F49","#EC83BA","#999999")
column_ha<- HeatmapAnnotation(
    Stimulation = stim,
    Celltype = cluster_anno,
    col = list(Stimulation = setNames(c("#999999","#888888","#777777","#666666"),unique(stim)), Celltype = setNames(heatmap_colors, unique(cluster_anno))),
    na_col = "grey",
    annotation_legend_param = list(Stimulation = list(at = unique(stim)), Celltype = list(at = unique(cluster_anno)))
)
index <- which(colnames(exp_mat) %in% levels(obj_subset))
exp_mat <- exp_mat[,index]
percent_mat <- percent_mat[,index]
pdf(paste(paste0(outdir, "/Proliferating_TCM_dotplot.pdf"), sep = ""), width = 10, height = 5)
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
        row_gap = unit(0, "mm")
        )
dev.off()


############# Figure S4B: Cellcycle phase scores barplot
CC_Phases <- read.table(paste0(indir, "/sc_phase_scores.txt"), header = T)
obj$Phase <- CC_Phases
obj$celltype2[which(obj$celltype2 == "Braking")] <- "TCM_act"
obj$celltype2[which(obj$celltype2 == "Naive_rest2")] <- "Naive_rest"
obj$celltype2[which(obj$celltype2 == "Naive_rest" | obj$celltype2 == "Naive_act")] <- "Naive"
obj$celltype2[which(obj$celltype2 == "TCM_rest" | obj$celltype2 == "TCM_act")] <- "TCM"
obj$celltype2[which(obj$celltype2 == "Th1_rest" | obj$celltype2 == "Th1_act")] <- "Th1"
obj$celltype2[which(obj$celltype2 == "Th2_rest" | obj$celltype2 == "Th2_act")] <- "Th2"
obj$celltype2[which(obj$celltype2 == "Th17_rest" | obj$celltype2 == "Th17_act")] <- "Th17"
obj$celltype_condition <- paste(obj$celltype2, obj$Condition, sep = "_")
celltype_condition <- c("Naive_Resting","Naive_2","Naive_5","Naive_15","TCM_Resting","TCM_2","TCM_5","TCM_15",
                        "Th1_Resting","Th1_2","Th1_5","Th1_15","Th2_Resting","Th2_2","Th2_5","Th2_15",
                        "Th17_Resting","Th17_2","Th17_5","Th17_15","TCM/TEM_2","TCM/TEM_5","TCM/TEM_15",
                        "TEM_act_Resting","TEM_act_2","TEM_act_5","TEM_act_15","TEM_act2_Resting","TEM_act2_2","TEM_act2_5","TEM_act2_15",
                        "Treg_Resting","Treg_2","Treg_5","Treg_15")
means <- c()
sds <- c()
Donors <- c("Donor1","Donor2","Donor3","Donor4")
for(i in celltype_condition){
    prop <- c()
    for(j in Donors){
        curr_cycling <- length(which(obj$celltype_condition == i & obj$Donor == j & obj$Phase %in% c("S","G2M")))
        curr_total <- length(which(obj$celltype_condition == i & obj$Donor == j))
        prop <- c(prop, curr_cycling / curr_total)
    }
    means <- c(means, mean(prop))
    sds <- c(sds, sd(prop))
}

celltype <- c("Naive","Naive","Naive","Naive","TCM","TCM","TCM","TCM", "Th1","Th1","Th1","Th1",
            "Th2","Th2","Th2","Th2","Th17","Th17","Th17","Th17","TCM/TEM","TCM/TEM","TCM/TEM","TEM","TEM","TEM","TEM",
            "TEM2","TEM2","TEM2","TEM2","Treg","Treg","Treg","Treg")
condition <- c("Resting","2","5","15","Resting","2","5","15","Resting","2","5","15","Resting","2","5","15","Resting","2","5","15",
            "2","5","15","Resting","2","5","15","Resting","2","5","15","Resting","2","5","15")
df <- data.frame(mean = means, sds = sds, celltype = celltype, condition = condition, celltype_condition = celltype_condition)
df$celltype_condition <- factor(df$celltype_condition, levels = celltype_condition)
df$celltype <- factor(df$celltype, levels = c("Naive","TCM","Th1","Th2","Th17","TCM/TEM","TEM","TEM2","Treg"))
df$condition <- factor(df$condition, levels = c("Resting","2","5","15"))
pdf(paste0(outdir, "/CC_barplot.pdf"), width = 10, height = 5)
ggplot(df, aes(x = celltype_condition, y = means, fill = celltype)) + geom_bar(stat = "identity") + geom_errorbar(aes(ymin=mean-sds, ymax=mean+sds), width=.2,position=position_dodge(.9)) +
    scale_fill_manual(values = c("#E41A1C","#4A72A6","#48A462","#7E6E85","#D16948","#E1C62F","#B75F49","#EC83BA","#999999")) +
    theme(panel.grid = element_line(color = "grey",linewidth = 0.1)) + 
    theme(panel.background = element_rect("white", fill = NA)) +
    theme(axis.title=element_text(size=15, face='plain')) + 
    theme(axis.text=element_text(size=14, face='plain'))
dev.off()

# Calculate Statistics of Cell-cycle
comparisons <- c("Naive_15","Th1_15","Th2_15","Th17_15","TEM_act_15","TEM_act_15","Treg_15")
donors <- c("Donor1","Donor2","Donor3","Donor4")
pvalues <- c()
for(i in 1:length(comparisons)){
    TCM_index <- which(obj$celltype_condition == "TCM_15")
    TCM <- c()
    for(donor in donors){
        donor_index <- which(obj$Donor[TCM_index] == donor)
        TCM <- c(TCM, length(which(obj$Phase[TCM_index[donor_index]] == "G1")) / length(donor_index))
    }

    curr_index <- which(obj$celltype_condition == comparisons[i])
    curr <- c()
    for(donor in donors){
        donor_index <- which(obj$Donor[curr_index] == donor)
        curr <- c(curr, length(which(obj$Phase[curr_index[donor_index]] == "G1")) / length(donor_index))
    }
    test <- t.test(TCM, curr)
    pvalues <- c(pvalues, test$p.value)
}

############ Fig 3B: GSEA of 15h rapid recall pathways
GO <- readRDS(paste0(indir, "/GO.rds"))
RR_genes <- read.table(paste0(indir, "/RR_df_Filtered.txt"), header = T)
TCM_15_RR_up <- RR_genes[which(RR_genes$Celltype == "TCM" & RR_genes$Timepoint == 15 & RR_genes$Direction == "Up"),1]
Th1_15_RR_up <- RR_genes[which(RR_genes$Celltype == "Th1" & RR_genes$Timepoint == 15 & RR_genes$Direction == "Up"),1]
Th2_15_RR_up <- RR_genes[which(RR_genes$Celltype == "Th2" & RR_genes$Timepoint == 15 & RR_genes$Direction == "Up"),1]
Th17_15_RR_up <- RR_genes[which(RR_genes$Celltype == "Th17" & RR_genes$Timepoint == 15 & RR_genes$Direction == "Up"),1]
TEM_15_RR_up <- RR_genes[which(RR_genes$Celltype == "TEM" & RR_genes$Timepoint == 15 & RR_genes$Direction == "Up"),1] 
TEM2_15_RR_up <- RR_genes[which(RR_genes$Celltype == "TEM2" & RR_genes$Timepoint == 15 & RR_genes$Direction == "Up"),1]
Treg_15_RR_up <- RR_genes[which(RR_genes$Celltype == "Treg" & RR_genes$Timepoint == 15 & RR_genes$Direction == "Up"),1] 
cell_cycle_genes <- GO[["GOBP_CELL_CYCLE"]]
pathways <- c("GOBP_CELL_CYCLE","GOBP_CELL_CYCLE_DNA_REPLICATION","GOBP_DNA_UNWINDING_INVOLVED_IN_DNA_REPLICATION",
    "GOBP_MEIOTIC_CELL_CYCLE","GOBP_CELL_CYCLE_G1_S_PHASE_TRANSITION","GOBP_CELL_CYCLE_G2_M_PHASE_TRANSITION",
    "GOBP_CHROMOSOME_ORGANIZATION_INVOLVED_IN_MEIOTIC_CELL_CYCLE")
TCM_CC <- TCM_15_RR_up[which(TCM_15_RR_up %in% cell_cycle_genes)] 
Th1_CC <- Th1_15_RR_up[which(Th1_15_RR_up %in% cell_cycle_genes)]
Th2_CC <- Th2_15_RR_up[which(Th2_15_RR_up %in% cell_cycle_genes)] 
Th17_CC <- Th17_15_RR_up[which(Th17_15_RR_up %in% cell_cycle_genes)]
TEM_CC <- TEM_15_RR_up[which(TEM_15_RR_up %in% cell_cycle_genes)]
TEM2_CC <- TEM2_15_RR_up[which(TEM2_15_RR_up %in% cell_cycle_genes)]
Treg_CC <- Treg_15_RR_up[which(Treg_15_RR_up %in% cell_cycle_genes)] 

index <- which(names(GO) %in% pathways)
GO_subset <- GO[index]
TCM_pvalues <- performOraHypergeometric(GO_subset, TCM_CC, background = unique(unlist(GO)))
Th1_pvalues <- performOraHypergeometric(GO_subset, Th1_CC, background = unique(unlist(GO)))
Th2_pvalues <- performOraHypergeometric(GO_subset, Th2_CC, background = unique(unlist(GO)))
Th17_pvalues <- performOraHypergeometric(GO_subset, Th17_CC, background = unique(unlist(GO)))
TEM_pvalues <- performOraHypergeometric(GO_subset, TEM_CC, background = unique(unlist(GO)))
TEM2_pvalues <- performOraHypergeometric(GO_subset, TEM2_CC, background = unique(unlist(GO)))
Treg_pvalues <- performOraHypergeometric(GO_subset, Treg_CC, background = unique(unlist(GO)))

pvalue_mat <- matrix(0, length(pathways), 7)
pvalue_mat[,1] <- TCM_pvalues[,2]
pvalue_mat[,2] <- Th1_pvalues[,2]
pvalue_mat[,3] <- Th2_pvalues[,2]
pvalue_mat[,4] <- Th17_pvalues[,2]
pvalue_mat[,5] <- TEM_pvalues[,2]
pvalue_mat[,6] <- TEM2_pvalues[,2]
pvalue_mat[,7] <- Treg_pvalues[,2]

rownames(pvalue_mat) <- pathways
colnames(pvalue_mat) <- c("TCM","Th1","Th2","Th17","TEM","TEM2","Treg")
pvalues <- c(TCM_pvalues[,2], Th1_pvalues[,2],Th2_pvalues[,2],Th17_pvalues[,2], TEM_pvalues[,2], TEM2_pvalues[,2], Treg_pvalues[,2])
pathways2 <- rep(pathways, 7)
celltype <- c(rep("TCM", length(pathways)), rep("Th1", length(pathways)), rep("Th2", length(pathways)),rep("Th17", length(pathways)),rep("TEM", length(pathways)), rep("TEM2", length(pathways)),rep("Treg", length(pathways)))
df <- data.frame(pvalue = -log10(pvalues), pathway = pathways2, celltype = celltype)
df$pathway <- factor(df$pathway, levels = rev(pathways))
df$celltype = factor(df$celltype, levels = c("TCM","Th1","Th2","Th17","TEM","TEM2","Treg"))

pdf(paste0(outdir, "/cellcycle_GSEA_bar.pdf"), width = 12, height = 5)
ggplot(df, aes(x = pathway, y = pvalue, fill = celltype)) + geom_bar(stat = "identity", position = "dodge") + scale_fill_manual(values = c("#4A72A6","#48A462","#7E6E85","#D16948","#EC83BA","#B75F49","#999999")) +
    theme(panel.background = element_rect("white")) + theme(panel.grid = element_line(color = "grey")) + coord_flip() + scale_y_continuous(trans = "pseudo_log") +
    geom_hline(yintercept = 1, linetype = 2, colour = "black")
dev.off()

############# Figure 3C: Heatmap of cell-cycle associated rapid recall pathways
counts <- read.table(paste0(indir, "/RNA_counts_combatseq_vst.txt"))
meta <- read.table(paste0(indir, "/meta_filter.txt"))

celltypes <- c("Naive_rest","Naive_rest2","TCM_rest","Th1_rest","Th2_rest","Th17_rest","Naive_act",
                "TCM_act","TCM/TEM","Th1_act","Th2_act","Th17_act","TEM_act",
                "TEM_act2","Treg", "Braking")
celltype_timepoint <- c("Naive_rest_Resting","Naive_rest2_Resting","TCM_rest_Resting","Th1_rest_Resting","Th2_rest_Resting","Th17_rest_Resting", "TEM_act_Resting","TEM_act2_Resting",
                        "Naive_act_2","Naive_act_5","Naive_act_15","TCM_act_2","TCM_act_5","TCM_act_15","TCM/TEM_2","TCM/TEM_5","TCM/TEM_15",
                        "Th1_act_2","Th1_act_5","Th1_act_15","Th2_act_2","Th2_act_5","Th2_act_15","Th17_act_2","Th17_act_5","Th17_act_15",
                        "TEM_act_2","TEM_act_5","TEM_act_15","TEM_act2_2","TEM_act2_5","TEM_act2_15",
                        "Treg_Resting","Treg_2","Treg_5","Treg_15",
                        "Braking_15")
var <- "celltype" # Variable in metadata to sort and make heatmap
var_index <- which(colnames(meta) == var)
nonrep_var <- "celltype_condition" # The variable in meta that doesn't include replicate
nonrep_var_index <- which(colnames(meta) == nonrep_var)
#heat_col_tfmrna <- colorRamp2(c(-2, 0, 2), c('dodgerblue','white','red')) ## Colors for heatmap
heat_col_tfmrna <- colorRamp2(c(-2.5, -2, 0, 2, 2.5), c('#2c7bb6','#abd9e9','white','#d7191c','#b10026'))
subset_cells <- which(meta$celltype_condition %in% celltype_timepoint) ## Cell types to subset NULL if all)
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
meta$celltype[which(meta$celltype == "Braking")] <- "TCM_Proliferation"
meta$celltype_condition <- paste(meta$celltype, meta$Condition, sep = "_")
rownames(meta) <- paste(meta$celltype, meta$Condition, meta$Donor, sep = "_")
colnames(counts) <- rownames(meta)
celltype_timepoint <- c("Naive_Resting","Naive_2","Naive_5","Naive_15",
                        "TCM_Resting","TCM_2","TCM_5","TCM_15",
                        "TCM_Proliferation_15",
                        "TCM/TEM_2","TCM/TEM_5","TCM/TEM_15",
                        "Th1_Resting","Th1_2","Th1_5","Th1_15",
                        "Th2_Resting","Th2_2","Th2_5","Th2_15",
                        "Th17_Resting","Th17_2","Th17_5","Th17_15",
                        "TEM_act_Resting","TEM_act_2","TEM_act_5","TEM_act_15",
                        "TEM_act2_Resting","TEM_act2_2","TEM_act2_5","TEM_act2_15",
                        "Treg_Resting","Treg_2","Treg_5","Treg_15")
celltypes <- c("Naive","TCM","TCM_Proliferation","TCM/TEM","Th1","Th2","Th17","TEM_act","TEM_act2","Treg")
subset_cells <- which(meta$celltype_condition %in% celltype_timepoint) ## Cell types to subset NULL if all)
## Subset counts if needed (should be needed)
if(!is.null(subset_cells) ){
    counts <- counts[,subset_cells]
    meta <- meta[subset_cells,]
}

counts <- as.matrix(counts)
colnames(counts) <- gsub("\\.", "/", colnames(counts))
counts <- t(scale(t(counts)))

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
heatmap_colors = c("#E41A1C","#4A72A6","#091b33","#48A462","#7E6E85","#D16948","#E1C62F","#B75F49","#EC83BA","#999999")
heatmap_colors <- setNames(heatmap_colors, c("Naive","TCM","TCM_Proliferation","TCM/TEM","Th1","Th2","Th17","TEM_act","TEM_act2","Treg"))
condition_colors <- c("#ebebeb","#C0C0C0","#A9A9A9","#808080")
names(condition_colors) <- c("Resting","2","5","15")
ha <- HeatmapAnnotation(condition = annotation_matrix[,4], celltype = annotation_matrix[,1], col = list(condition = condition_colors, celltype = heatmap_colors), simple_anno_size = unit(1, "cm"),annotation_legend_param = list(celltype = list(at = celltypes)))

splits <- unit(c(0,0,0,1.5,0,0,0,1.5,1.5,0,0,1.5,0,0,0,1.5,0,0,0,1.5,0,0,0,1.5,0,0,0,1.5,0,0,0,1.5,0,0,0), "mm")

cluster_df <- read.table(paste0(indir, "/cellcycle_clustering.txt"), header = T)
gene_labels <- c("NDC80","CDKN1B","PPP2R5C","SMAD3","WEE1","SGO1","ZBTB17","GADD45A","GADD45B","YWHAQ","YWHAG","MYC","HDAC2","CDK4",
                "E2F5","ANAPC1","CCND2","CDK6","CDKN1A","MCM3","MCM2","ORC6","E2F1","CDC45","MCM6","CHEK1",
                "MCM7","ORC3","MCM5","CDC7","MCM4","CCNB1","DBF4B", "CDKN1B","SGO1","NDC80","MYB","MCMBP","CCNI2","CDC37","NPM1","CDC123","CDK6","CCND2",
                "PES1","CINP","CDK7","CENPO","CENPN","LMNA","TUBB","HELLS")
gene_labels <- gene_labels[which(gene_labels %in% cluster_df[,1])]
index <- match(gene_labels,cluster_df[,1])
rowAnn_across_left <-  HeatmapAnnotation(foo = anno_mark(at = index, labels = gene_labels, side = 'left',labels_gp = gpar(fontface="italic", fontsize = 8)),
                                    which = 'row',
                                    #annotation_width = unit(c(1, 4), 'cm'),
                                    #gap = unit(0, 'mm'),
                                    show_annotation_name = F)

pdf(paste0(outdir, "/CC_Heatmap_RR.pdf"), width = 10, height = 5)
ht <- Heatmap(counts[cluster_df[,1], ], name = "Z-score", show_row_names = F,  show_row_dend = F, show_column_dend = F, 
            row_split = factor(cluster_df[,2], levels = c("C1","C2","C3")), 
            border = T, 
            cluster_rows = T, show_column_names = F, column_split = factor(annotation_matrix[,2], levels = celltype_timepoint),
            cluster_row_slices = F, col = heat_col_tfmrna, column_gap = splits, row_gap = unit(0, "mm"), 
            cluster_columns = F , top_annotation = ha, row_names_gp = (gpar(fontsize = 8)), use_raster = T,
            right_annotation = rowAnn_across_left,
            column_title = NULL, row_title = NULL)
draw(ht)
dev.off()
cluster1 <- cluster_df[which(cluster_df[,2] == "C1"),1]
cluster2 <- cluster_df[which(cluster_df[,2] == "C2"),1]
cluster3 <- cluster_df[which(cluster_df[,2] == "C3"),1]
counts <- as.matrix(counts)
plotPattern(cluster1, meta, counts, paste0(outdir, "/cluster1_pattern.pdf"))
plotPattern(cluster2, meta, counts, paste0(outdir, "/cluster2_pattern.pdf"))
plotPattern(cluster3, meta, counts, paste0(outdir, "/cluster3_pattern.pdf"))

########## Figure S4F: Heatmap of TCM-only rapid recall genes associated with cell-cycle
index <- which(!(TCM_CC %in% Th1_CC | TCM_CC %in% Th2_CC | TCM_CC %in% Th17_CC | TCM_CC %in% TEM_CC | TCM_CC %in% TEM2_CC | TCM_CC %in% Treg_CC))
TCM_only <- TCM_CC[index]
pdf(paste0(outdir, "/CC_Heatmap_RR_TCMOnly.pdf"), width = 10, height = 5)
ht <- Heatmap(counts[TCM_only, ], name = "Z-score", show_row_names = T,  show_row_dend = F, show_column_dend = F, 
            border = T, 
            cluster_rows = T, show_column_names = F, column_split = factor(annotation_matrix[,2], levels = celltype_timepoint),
            cluster_row_slices = F, col = heat_col_tfmrna, column_gap = splits, row_gap = unit(0, "mm"), 
            cluster_columns = F , top_annotation = ha, row_names_gp = (gpar(fontsize = 8)), use_raster = T,
            column_title = NULL, row_title = NULL)
draw(ht)
dev.off()

########### Figure S4E: Heatmap of TCM/TEM marker genes
## TCM/TEM Marker Heatmap
genes <- c("IFNG","IL5","IL13","IL9","IL17A","IL17F","IL21","MAF","PRDM1","KLF6","TET2","RUNX2","RUNX3")
pdf(paste0(outdir, "/TCM_TEM_Heatmap.pdf"), width = 10, height = 10)
ht <- Heatmap(counts[genes, ], name = "Z-score", show_row_names = T,  show_row_dend = F, show_column_dend = F, 
            cluster_rows = F, show_column_names = F, column_split = factor(annotation_matrix[,2], levels = celltype_timepoint),
             col = heat_col_tfmrna, column_gap = splits, row_gap = unit(0, "mm"), 
              cluster_columns = F , top_annotation = ha, row_names_gp = (gpar(fontsize = 8)), use_raster = T,
              column_title = NULL, row_title = NULL)
draw(ht)
dev.off()

############ Figure S4C-D: Lineplot of TET2 and MYC expression
counts <- read.table(paste0(indir, "/RNA_lineplot_counts.txt"))
meta <- read.table(paste0(indir, "/RNA_lineplot_meta.txt"))
celltypes <- c("Naive","TCM","Th1","Th2","Th17","TCM/TEM","MHCII","CTL","Treg")
colors <- c("#E41A1C","#4A72A6","#48A462","#7E6E85","#D16948","#E1C62F","#B75F49","#EC83BA","#999999")
plotGene(counts, meta, "TET2", celltypes, colors, outdir)
plotGene(counts, meta, "MYC", celltypes, colors, outdir)

############ Figure3F: Coverage plots of TCM/TEM associated genes
## Coverage Plots
index <- which(obj$celltype2 %in% c("Naive","TCM","TCM/TEM","Th1","Th2","Th17") & obj$Condition %in% c(2,5,15))
obj_subset <- subset(obj, cells = index)
Idents(obj_subset) <- obj_subset$celltype2
levels(obj_subset) <- c("Naive","TCM","TCM/TEM","Th1","Th2","Th17")
DefaultAssay(obj_subset) <- "ATAC"

pdf(paste0(outdir, "/MAF_accessibility.pdf"))
CoveragePlot(obj_subset, region = "MAF", extend.upstream = 3000, extend.downstream = 3000) &
    scale_fill_manual(values = c("#E41A1C","#4A72A6","#E1C62F","#48A462","#7E6E85","#D16948"))
dev.off()

pdf(paste0(outdir, "/IL21_accessibility.pdf"))
CoveragePlot(obj_subset, region = "IL21", extend.upstream = 3000, extend.downstream = 3000) &
    scale_fill_manual(values = c("#E41A1C","#4A72A6","#E1C62F","#48A462","#7E6E85","#D16948"))
dev.off()

pdf(paste0(outdir, "/IFNG_accessibility.pdf"))
CoveragePlot(obj_subset, region = "IFNG", extend.upstream = 3000, extend.downstream = 3000) &
    scale_fill_manual(values = c("#E41A1C","#4A72A6","#E1C62F","#48A462","#7E6E85","#D16948"))
dev.off()

############# Figure 3E: Cytokine Box Plots
counts <- read.table(paste0(indir, "/RNA_counts_combatseq_vst.txt"))
meta <- read.table(paste0(indir, "/meta_filter.txt"))

index <- which(meta$celltype %in% c("Th1_act","TCM_act","TCM/TEM"))
counts <- counts[,index]
counts <- as.matrix(counts)
meta <- meta[index, ]
#meta$celltype[which(meta$celltype == "TCM_rest")] <- "TCM_act"


df <- data.frame(celltype = meta$celltype, condition = meta$Condition, celltype_condition = paste(meta$celltype, meta$Condition, sep = "_"), IFNG = as.vector(counts["IFNG",]), IL21 = as.vector(counts["IL21",]),TNF = as.vector(counts["TNF",]))
index <- which(df$celltype_condition %in% c("Th1_act_2","Th1_act_5","Th1_act_15","TCM_act_2","TCM_act_5","TCM_act_15","TCM/TEM_2","TCM/TEM_5","TCM/TEM_15"))
df <- df[index, ]
df$celltype <- factor(df$celltype, levels = c("TCM_act","Th1_act","TCM/TEM"))
df$condition = factor(df$condition, levels = c("2","5","15"))
#df$celltype_condition <- factor(df$celltype_condition, levels = c("Th1_act_2","Th1_rest_Resting"))

pdf(paste0(outdir, "/IFNG_barplot.pdf"), width = 3, height = 3)
ggplot(df, aes(x = condition, y = IFNG, fill = celltype)) + geom_boxplot(outlier.colour="white", notch=FALSE) +
    theme(panel.grid = element_line(color = "grey",linewidth = 0.1)) + 
    theme(panel.background = element_rect("white", fill = NA)) +
    theme(axis.title=element_text(size=15, face='plain')) + 
    theme(axis.text=element_text(size=14, face='plain')) + 
    scale_fill_manual(values = c("#4A72A6","#48A462","#E1C62F"))
dev.off()

pdf(paste0(outdir, "/TNF_barplot.pdf"), width = 3, height = 3)
ggplot(df, aes(x = condition, y = TNF, fill = celltype)) + geom_boxplot(outlier.colour="white", notch=FALSE) +
    theme(panel.grid = element_line(color = "grey",linewidth = 0.1)) + 
    theme(panel.background = element_rect("white", fill = NA)) +
    theme(axis.title=element_text(size=15, face='plain')) + 
    theme(axis.text=element_text(size=14, face='plain')) + 
    scale_fill_manual(values = c("#4A72A6","#48A462","#E1C62F"))
dev.off()

pdf(paste0(outdir, "/IL21_barplot.pdf"), width = 3, height = 3)
ggplot(df, aes(x = condition, y = IL21, fill = celltype)) + geom_boxplot(outlier.colour="white", notch=FALSE) +
    theme(panel.grid = element_line(color = "grey",linewidth = 0.1)) + 
    theme(panel.background = element_rect("white", fill = NA)) +
    theme(axis.title=element_text(size=15, face='plain')) + 
    theme(axis.text=element_text(size=14, face='plain')) + 
    scale_fill_manual(values = c("#4A72A6","#48A462","#E1C62F"))
dev.off()


############# Figure 3H: Heatmap of TCM/TEM associated transcription factors
res2 <- read.table(paste0(indir, "/TCMTEM_vs_TCM_2.txt"))
res5 <- read.table(paste0(indir, "/TCMTEM_vs_TCM_5.txt"))
res15 <- read.table(paste0(indir, "/TCMTEM_vs_TCM_15.txt"))
TFs <- c("MAF","ZEB2","PRDM1","TET2","RORA","KLF6","SMAD3","GATA3","RUNX3","MYBL1","MYB","TBX21","RORC","RUNX2")
res2 <- res2[TFs,]
res5 <- res5[TFs,]
res15 <- res15[TFs,]

mat <- cbind(res2[,2], res5[,2], res15[,2])
mat2 <- cbind(res2[,5], res5[,5], res15[,5])
rownames(mat) <- TFs
colnames(mat) <- c("2","5","15")
cell_fun = function(j, i, x, y, w, h, fill) {
    if(mat2[i, j] < 0.1) {
        grid.text("*", x, y)
    } else if(mat2[i, j] < 0.01) {
        grid.text("**", x, y)
    }
}

heat_col_tfmrna <- colorRamp2(c(-1, 0, 2), c("#4894C4", "white","#D5614E"))

heat_col_tfmrna <- colorRamp2(c(-2.5, -2, 0, 2, 2.5), c('#2c7bb6','#abd9e9','white','#d7191c','#b10026'))
pdf(paste0(outdir, "/TF_Heatmap.pdf"), width = 2.5, height = 6)
Heatmap(mat, cluster_rows = T, col = heat_col_tfmrna, show_row_names = T, show_column_names = T, cluster_columns = F, cell_fun = cell_fun,
    show_row_dend = F)
dev.off()