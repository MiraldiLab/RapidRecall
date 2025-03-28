source('Figure7_Functions.R')
library(bedr)
library(ggplot2)
library(ComplexHeatmap)
library(dplyr)
library(tidyr)
library(reshape2)
library(RColorBrewer)
library(circlize)


outdir <- "output/Figure7"
indir <- "input"

## Fig7A: RELI heatmap
res <- read.table(paste0(indir, "/RELI_res.txt"), header = T)
order_mat <- read.table(paste0(indir, "/RELI_order_mat.txt"), header = T)
order_context <- c("Rest_cluster1","Rest_cluster2","Rest_cluster3","Rest_cluster4","Rest_cluster5","Rest_cluster6",
                    "cluster1","cluster2","cluster3","cluster4","cluster5","cluster6","cluster7",
                    "Stim_cluster1","Stim_cluster2","Stim_cluster3","Stim_cluster4","Stim_cluster5")
index <- which(res$Overlap < 3)
res$Log10P[index] <- 0
res$Padj[index] <- 1
order_disease <- unique(res$Phenotype[which(res$Padj < 0.01)])
pvalue_matrix<-res[,c('Phenotype','ATAC_library','Log10P')]
pvalue_matrix<-pvalue_matrix[pvalue_matrix$Phenotype%in%order_disease,]
pvalue_matrix$Phenotype<-factor(pvalue_matrix$Phenotype,level=order_disease)
pvalue_matrix<-pvalue_matrix[pvalue_matrix$ATAC_library%in%order_context,]
count <- acast(pvalue_matrix, Phenotype~ATAC_library, value.var="Log10P")
count<-count[,order_context]
count_scale<-t(scale(t(count)))

library_sizes <- c(3832,5646,4383,2904,13860,2036,4517, 8363,6441,10578,6097,4651,5486,7927,8136,8733,8691,6134)
library_size_annot <- HeatmapAnnotation("Library Size" = anno_barplot(library_sizes, gp = gpar(fill = "grey80"), height = unit(2, "cm")))
order_mat <- order_mat[-39,]
index <- which(order_mat[,1] %in% order_disease)
order_mat <- order_mat[index, ]
heat_col_tfmrna <- colorRamp2(c(0, 4, 8), hcl_palette = "inferno") ## Colors for heatmap
group <- c("Rest","Rest","Rest","Rest","Rest","Rest","RR","RR","RR","RR","RR","RR","RR","Stim","Stim","Stim","Stim","Stim")
ha <- HeatmapAnnotation(group_peak = group, group = group, show_legend = F)

custom_names <- c("Asthma (Adult Onset)", "Asthma (Onset Age)", "Allergic Diseases (Onset Age)", "Allergic Diseases Multivariate","Respiratory Diseases",
    "Atopic Dermatitis","Adrenergic Inhalant Use","Nonatopic Asthma","Asthma (Moderate/Severe)","Asthma + Allergic Diseases","Hypothyroidism","Autoimmune Thyroid",
    "Rheumatoid Arthritis","Type1 Diabetes","Asthma (Child Onset)","Primary Biliary Cholangitis","Primary Biliary Cirrhosis","Autoimmune Traits","Celiac + RA",
    "RA (ACPA+)","Ankylosing Spondylitis","Primary Sclerosing Cholangitis","Allergic Rhinitis","Vericose Veins","Psoriasis",
    "Inflammatory Skin Disease","Macular Degeneration","Atopic Dermatitis","Multiple Myeloma","COPD","Rheumatic Diseases","MS","Hay Fever + Eczema",
    "Allergic Diseases","Asthma","Eczema","UC","IBD","CD","Autoimmune Traits Pleiotropy","Celiac Disease")
column_cluster <- c("G1","G1","G2","G2","G2","G2","G3","G3","G3","G3","G3","G3","G3","G4","G4","G4","G4","G4")
count_subset <- count[order_mat[,1],]
pdf(paste0(Figure7, "/RELI_heatmap.txt"), width = 7, height = 7)
Heatmap(count[order_mat[,1],], row_labels = custom_names, show_column_names = F, top_annotation = ha, col = heat_col_tfmrna,
        cluster_columns = F, row_names_gp = gpar(fontsize = 8), show_row_dend = F, cluster_rows = T,
        row_split = factor(order_mat[,2], levels = c("group1", "group2", "group3","group4")), bottom_annotation = library_size_annot, 
        column_split = factor(column_cluster, levels = c("G1","G2","G3","G4")), 
        cell_fun = function(j, i, x, y, w, h, fill) {
            if(count_subset[i, j] > 1) {
                grid.text("*", x, y)
            }
        })
dev.off()

# Read in list of phenotypes I will calculate networks for
phenotypes <- list.dirs(paste0(indir, "/RELI_overlaps"), full.names=F)
phenotypes <- phenotypes[2:length(phenotypes)]
all_TFs <- read.table(paste0(indir, "/all_TFs.txt"))$V1

list_names <- c(paste0("RR", 1:7), paste0("Rest", 1:6), paste0("Stim", 1:5))
lists <- setNames(lapply(list_names, function(x) list()), list_names)
# Base path
base_path <- paste0(indir, "/RELI_overlaps/")
# Loop through phenotypes and clusters
for (phenotype in phenotypes) {
  for (cluster_type in c("RR", "Rest", "Stim")) {
    for (cluster_num in 1:7) {
      # Construct the file path
      file_path <- paste0(base_path, phenotype, "/", phenotype, "_", cluster_type, "_cluster", cluster_num, "_loci.bed")
      
      # Check if the file exists and read it
      if (file.exists(file_path)) {
        list_name <- paste0(cluster_type, cluster_num)
        if (list_name %in% names(lists)) {
          lists[[list_name]][[phenotype]] <- read.table(file_path, header = TRUE)
        }
      }
    }
  }
}

# Combine all overlaps into a list
overlap_list <- lists
# Get all the peaks that overlap risk loci
peaks <- read.table(paste0(indir, "/peak_names.bed"))
for(k in 1:length(overlap_list)){
    curr_overlap <- overlap_list[[k]]
    print(k)
    for(i in 1:length(curr_overlap)){
        curr_overlap[[i]]$intersecting_peak <- "None"
        for(j in 1:length(curr_overlap[[i]][,1])){
            index <- which(peaks[,1] == curr_overlap[[i]][j,1] & peaks[,2] <= curr_overlap[[i]][j,2] & peaks[,3] >= curr_overlap[[i]][j,3])
            if(length(index == 1)){
                curr_overlap[[i]][j,5] <- paste(peaks[index,1], peaks[index,2], peaks[index, 3], sep = "-")
            }
        }
    }
    overlap_list[[k]] <- curr_overlap
}

# Combine data into single matrix
combined_matrix <- do.call(rbind, lapply(seq_along(overlap_list), function(i) {
  do.call(rbind, lapply(seq_along(overlap_list[[i]]), function(j) {
    cbind(overlap_list[[i]][[j]], Source1 = i, Source2 = j)
  }))
}))
colnames(combined_matrix) <- c(colnames(overlap_list[[1]][[1]]), "PeakSet", "Phenotype")
combined_matrix$Phenotype <- phenotypes[combined_matrix$Phenotype]
peak_set_names <- c("RR1","RR2","RR3","RR4","RR5","RR6","RR7","Rest1","Rest2","Rest3","Rest4","Rest5","Rest6","Stim1","Stim2","Stim3","Stim4","Stim5")
combined_matrix$PeakSet <- peak_set_names[combined_matrix$PeakSet]

combined_matrix$tmp <- paste(combined_matrix$SNP, combined_matrix$Phenotype, sep = "-")
combined_matrix <- combined_matrix[which(!duplicated(combined_matrix$tmp)),]
combined_matrix$tmp <- NULL

peak_to_gene <- read.table("/data/miraldiNB/Katko/Projects/Barski_CD4_Multiome/Outs/PeakToGene/combined.txt", header = T)
peak_lookup <- setNames(
    peak_to_gene[,4],  # the gene lists
    paste(peak_to_gene[,1], peak_to_gene[,2], peak_to_gene[,3], sep="-")  # peak identifiers in chr-start-stop format
)

# Add a new column to the SNP matrix for the genes
combined_matrix <- cbind(combined_matrix, Gene = "None")
# Fill in the gene associations
combined_matrix[, "Gene"] <- peak_lookup[combined_matrix[,5]]
# If any lookups failed (NA), replace with "None"
combined_matrix[is.na(combined_matrix[, "Gene"]), "Peak_Associated_Genes"] <- "None"
five_percent <- read.table(paste0(indir, "/five_percent_genes.txt"))$V1

Th2_MHCII <- c("Asthma_moderate_or_severe","Medication_use_adrenergics_inhalants","Allergic_disease_asthma_hay_fever_and_or_eczema__age_of_onset","Asthma_adult_onset","Atopic_asthma","Nonatopic_asthma","Allergic_disease_asthma_hay_fever_or_eczema","Respiratory_diseases")
Treg <- c("Rheumatoid_arthritis","Type_1_diabetes","Primary_biliary_cholangitis","Hypothyroidism","Autoimmune_thyroid_disease")
Th1_Th17_CTL <- c("CD","UC","IBD","Systemic_lupus_erythematosus","Autoimmune_traits_pleiotropy","Celiac_disease","Multiple_sclerosis","Hay_fever_and_or_eczema","Eczema","Chronic_inflammatory_diseases_ankylosing_spondylitis_Crohn_s_disease_psoriasis_primary_sclerosing_cholangitis_ulcerative_colitis__pleiotropy")

Th2_MHCII_sets <- c("Rest4")
Treg_sets <- c("Rest6")
Th1_Th17_CTL_sets <- c("Rest5")

Th2_MHCII_genes <- extract_genes(Th2_MHCII, Th2_MHCII_sets)
Treg_genes <- extract_genes(Treg, Treg_sets)
Th1_Th17_CTL_genes <- extract_genes(Th1_Th17_CTL, Th1_Th17_CTL_sets)

network_list <- readRDS(paste0(indir, "/network_list.rds"))
net <- read.table(paste0(indir, "/GRN.tsv"), header = T)

Th2_MHCII_TFs <- find_TFs(Th2_MHCII_genes)
Treg_TFs <- find_TFs(Treg_genes)
Th1_Th17_CTL_TFs <- find_TFs(Th1_Th17_CTL_genes)

Th2_MHCII_TFs <- unique(c(Th2_MHCII_TFs, Th2_MHCII_genes[which(Th2_MHCII_genes %in% all_TFs)]))
Treg_TFs <- unique(c(Treg_TFs, Treg_genes[which(Treg_genes %in% all_TFs)]))
Th1_Th17_CTL_TFs <- unique(c(Th1_Th17_CTL_TFs, Th1_Th17_CTL_genes[which(Th1_Th17_CTL_genes %in% all_TFs)]))

#### Custom list sorted manually for matching dynamics and clustered in correct order
Th2_MHCII_genes <- c("IL10","AHI1","IQGAP1","TH2LCRR","RAD50","IL13","KIF1B","HLA-DQB1")
Th2_MHCII_TFs <- c("GATA3","TFEB","NFATC2","ZNF365","KLF13","ZBTB44","KLF11","XBP1","KLF6","PATZ1","RUNX2","ZNF200","NR3C1",   
                "ZBTB49","SP140","ZFY","TRERF1","RFX1","IKZF3","TCF7L2","FOSL1","ZBTB10")
Treg_genes <- c("IL2RA","IL15RA","CTLA4","CLEC16A","HLA-DR4","ZC3H12D")
Treg_TFs <- c("HSF5","ZNF292","FOXP3","IKZF2","PRDM1","ETV6","AHR","XBP1","ATF4","NFAT5","ZNF549","MYB","ZBTB38",
            "ATF6","BACH2","ZNF235","PBX2")

Th2_MHCII_net <- net[which(net$Gene %in% Th2_MHCII_genes & net$TF %in% Th2_MHCII_TFs),]
Treg_net <- net[which(net$Gene %in% Treg_genes & net$TF %in% Treg_TFs),]
Th1_Th17_CTL_net <- net[which(net$Gene %in% Th1_Th17_CTL_genes & net$TF %in% Th1_Th17_CTL_TFs),]

#### Correlation comparison
combined_matrix_Th2_MHCII <- combined_matrix[which(combined_matrix$Phenotype %in% Th2_MHCII & combined_matrix$PeakSet %in% Th2_MHCII_sets),]
combined_matrix_Treg <- combined_matrix[which(combined_matrix$Phenotype %in% Treg & combined_matrix$PeakSet %in% Treg_sets),]

counts_ATAC <- read.table(paste0(indir, "/RELI_counts/counts_ATAC.txt"))
counts_RNA <- read.table(paste0(indir, "/RELI_counts/counts_RNA.txt"))
TFA <- read.table(paste0(indir, "/RELI_counts/TFA.txt"))
meta <- read.table(paste0(indir, "/RELI_counts/meta.txt"))

resting_index <- which(meta$Condition == "Resting")
act_index <- which(meta$Condition != "Resting")
combined_matrix_Treg_cor <- process_combined_matrix(combined_matrix_Treg, counts_ATAC, counts_RNA, TFA, resting_index, meta, net, Treg_TFs, Treg_genes, net)
combined_matrix_Th2_MHCII_cor <- process_combined_matrix(combined_matrix_Th2_MHCII, counts_ATAC, counts_RNA, TFA, resting_index, meta, net, Th2_MHCII_TFs, Th2_MHCII_genes, net)

write.table(combined_matrix_Th2_MHCII_cor, 
    paste0(outdir, "/Th2_MHCII_cor.txt"), row.names = F, col.names = T, quote = F, sep = "\t")
write.table(combined_matrix_Treg_cor, 
    paste0(outdir, "Treg_cor.txt"), row.names = F, col.names = T, quote = F, sep = "\t")

Treg_net <- filter_grn(combined_matrix_Treg_cor, Treg_TFs)
Th2_MHCII_net <- filter_grn(combined_matrix_Th2_MHCII_cor, Th2_MHCII_TFs)

write.table(Th2_MHCII_net, 
    paste0(outdir, "/Th2_MHCII_net.txt"), row.names = F, col.names = T, quote = F, sep = "\t")
write.table(Treg_net, 
    paste0(outdir, "/Treg_net.txt"), row.names = F, col.names = T, quote = F, sep = "\t")

Th2_MHCII_nodes <- node_list(Th2_MHCII_net)
Treg_nodes <- node_list(Treg_net)

write.table(Treg_nodes, paste0(outdir, "/Treg_nodes.txt"), row.names = F,
    col.names = T, quote = F, sep = "\t")
write.table(Th2_MHCII_nodes, paste0(outdir, "/Th2_MHCII_nodes.txt"), row.names = F,
    col.names = T, quote = F, sep = "\t")


########## Plot peaks, genes, TFA
counts_RNA_rest <- read.table(paste0(indir, "/RELI_counts/counts_RNA_rest.txt"))
counts_ATAC_rest <- read.table(paste0(indir, "/RELI_counts/counts_ATAC_rest.txt"))
TFA_rest <- read.table(paste0(indir, "/RELI_counts/TFA_rest.txt"))
meta_rest <- read.table(paste0(indir, "/RELI_counts/meta_rest.txt"))

RNA_features <- c("IL10","IQGAP1","IL13","RAD50","TH2LCRR","KIF1B","HLA-DQB1","AHI1")
names(RNA_features) <- RNA_features
ATAC_features <- c("chr1-10325475-10325866", "chr1-10325274-10325465","chr15-90479121-90479474","chr5-132629373-132629932",
                "chr5-132656344-132656807","chr6-135506345-135506573","chr1-206770529-206770849","chr6-32666447-32666577")
names(ATAC_features) <- c("rs2182326", "rs17401966", "rs12440447", "rs3798135", "rs1881457", "rs7750586", "rs3024493", "rs766371114")
TF_features <- c("GATA3","ZNF365","NFATC2","PATZ1","TRERF1","TFEB","RUNX2","SP140","KLF13","KLF6","RFX1","ZBTB49","ZBTB44","KLF11","ZNF200","IKZF3")
names(TF_features) <- TF_features

plot_counts_rest(counts_RNA_rest, meta_rest, features = RNA_features, figure_w = 7, figure_h = 2.5, output_name = paste0(outdir, "/MHCII_targetgenes_rest.pdf"))
plot_counts_rest(TFA_rest, meta_rest, features = TF_features, figure_w = 7, figure_h = 4, output_name = paste0(outdir, "/MHCII_TFA_rest.pdf"))
plot_counts_rest(counts_ATAC_rest, meta_rest, features = ATAC_features, figure_w = 7, figure_h = 2.5, output_name = paste0(outdir, "/MHCII_peaks_rest.pdf"))
plot_counts_act(counts_RNA, meta, features = RNA_features, figure_w = 14, figure_h = 2.5, output_name = paste0(outdir, "/MHCII_targetgenes_dynamic.pdf"))
plot_counts_act(TFA, meta, features = TF_features, figure_w = 14, figure_h = 4, output_name = paste0(outdir, "/MHCII_TFA_dynamic.pdf"))
plot_counts_act(counts_ATAC, meta, features = ATAC_features, figure_w = 14, figure_h = 2.5, output_name = paste0(outdir, "/MHCII_peaks_dynamic.pdf"))

RNA_features = c("CLEC16A","CTLA4","ZC3H12D","ZNF292","IL2RA","IL15RA")
names(RNA_features) <- RNA_features
TF_features = c("FOXP3","IKZF2","ZBTB38","ZNF292","ZNF549","XBP1","PRDM1","ATF6","ATF4","MYB","AHR","BACH2","HSF5")
names(TF_features) <- TF_features
ATAC_features <- c("chr2-203835923-203836636", "chr10-6037066-6037374", "chr16-11088072-11088709", "chr6-149495386-149495656")
names(ATAC_features) <- c("rs12990970", "rs7893467", "rs7203793", "rs4897131")

plot_counts_rest(counts_RNA_rest, meta_rest, features = RNA_features, figure_w = 7, figure_h = 2.5, output_name = paste0(outdir, "/Treg_targetgenes_rest.pdf"))
plot_counts_rest(TFA_rest, meta_rest, features = TF_features, figure_w = 7, figure_h = 4, output_name = paste0(outdir, "/Treg_TFA_rest.pdf"))
plot_counts_rest(counts_ATAC_rest, meta_rest, features = ATAC_features, figure_w = 7, figure_h = 2.5, output_name = paste0(outdir, "/Treg_peaks_rest.pdf"))
plot_counts_act(counts_RNA, meta, features = RNA_features, figure_w = 14, figure_h = 2.5, output_name = paste0(outdir, "/Treg_targetgenes_dynamic.pdf"))
plot_counts_act(TFA, meta, features = TF_features, figure_w = 14, figure_h = 4, output_name = paste0(outdir, "/Treg_TFA_dynamic.pdf"))
plot_counts_act(counts_ATAC, meta, features = ATAC_features, figure_w = 14, figure_h = 2.5, output_name = paste0(outdir, "/Treg_peaks_dynamic.pdf"))


####### Fig S6A: Overlap of the various ATAC peak sets
RR_clustering <- read.table(paste0(indir, "/ATAC_RR_clustering.txt"), header = T)
RR1 <- RR_clustering[which(RR_clustering[,2] == "C1"),1]
RR2 <- RR_clustering[which(RR_clustering[,2] == "C2"),1]
RR3 <- RR_clustering[which(RR_clustering[,2] == "C3"),1]
RR4 <- RR_clustering[which(RR_clustering[,2] == "C4"),1]
RR5 <- RR_clustering[which(RR_clustering[,2] == "C5"),1]
RR6 <- RR_clustering[which(RR_clustering[,2] == "C6"),1]
RR7 <- RR_clustering[which(RR_clustering[,2] == "C7"),1]

Act_clustering <- read.table(paste0(indir, "/ATAC_Stim_clustering.txt"), header = T)
Stim1 <- Act_clustering[which(Act_clustering[,2] == "C1"),1]
Stim2 <- Act_clustering[which(Act_clustering[,2] == "C2"),1]
Stim3 <- Act_clustering[which(Act_clustering[,2] == "C3"),1]
Stim4 <- Act_clustering[which(Act_clustering[,2] == "C4"),1]
Stim5 <- Act_clustering[which(Act_clustering[,2] == "C5"),1]

Rest_clustering <- read.table(paste0(indir, "/ATAC_Resting_clustering.txt"), header = T)
Rest1 <- Rest_clustering[which(Rest_clustering[,2] == "C1"),1]
Rest2 <- Rest_clustering[which(Rest_clustering[,2] == "C2"),1]
Rest3 <- Rest_clustering[which(Rest_clustering[,2] == "C3"),1]
Rest4 <- Rest_clustering[which(Rest_clustering[,2] == "C4"),1]
Rest5 <- Rest_clustering[which(Rest_clustering[,2] == "C5"),1]
Rest6 <- Rest_clustering[which(Rest_clustering[,2] == "C6"),1]

Peak_list <- list(Rest1, Rest2, Rest3, Rest4, Rest5, Rest6, RR1, RR2, RR3, RR4, RR5, RR6, RR7, Stim1, Stim2, Stim3, Stim4, Stim5)
names(Peak_list) <- c("Rest1","Rest2","Rest3","Rest4","Rest5","Rest6","RR1","RR2","RR3","RR4","RR5","RR6","RR7","Stim1","Stim2","Stim3","Stim4","Stim5")
overlap_matrix <- matrix(0, length(Peak_list), length(Peak_list))
rownames(overlap_matrix) <- names(Peak_list)
colnames(overlap_matrix) <- names(Peak_list)
for(i in names(Peak_list)){
    for(j in names(Peak_list)){
        overlap_matrix[i,j] <- length(which(Peak_list[[i]] %in% Peak_list[[j]]))
    }
}

overlap_matrix_frac <- overlap_matrix
for(i in colnames(overlap_matrix_frac)){
    overlap_matrix_frac[,i] <- overlap_matrix_frac[,i] / length(Peak_list[[i]])
}

heat_col_tfmrna <- colorRamp2(c(0,0.1,0.25, 0.75), colors = c("#040404","#661A66","#C34A6A","#EBA948"))
pdf(paste0(outdir, "/peakset_overlap.pdf"), width = 10, height = 10)
Heatmap(overlap_matrix_frac, show_row_names = T, show_column_names = T, cluster_rows = F, cluster_columns = F, col = heat_col_tfmrna,
    cell_fun = function(j, i, x, y, w, h, col) {grid.text(overlap_matrix[i, j], x, y,gp = gpar(col = "white"))})
dev.off()

## Gene set overlap
RR_clustering <- read.table(paste0(indir, "/RNA_RR_clustering.txt"), header = T)
RR1 <- RR_clustering[which(RR_clustering[,2] == "C1"),1]
RR2 <- RR_clustering[which(RR_clustering[,2] == "C2"),1]
RR3 <- RR_clustering[which(RR_clustering[,2] == "C3"),1]
RR4 <- RR_clustering[which(RR_clustering[,2] == "C4"),1]
RR5 <- RR_clustering[which(RR_clustering[,2] == "C5"),1]
RR6 <- RR_clustering[which(RR_clustering[,2] == "C6"),1]
RR7 <- RR_clustering[which(RR_clustering[,2] == "C7"),1]
RR8 <- RR_clustering[which(RR_clustering[,2] == "C8"),1]
RR9 <- RR_clustering[which(RR_clustering[,2] == "C9"),1]

Act_clustering <- read.table(paste0(indir, "/RNA_Stim_clustering.txt"), header = T)
Stim1 <- Act_clustering[which(Act_clustering[,2] == "C1"),1]
Stim2 <- Act_clustering[which(Act_clustering[,2] == "C2"),1]
Stim3 <- Act_clustering[which(Act_clustering[,2] == "C3"),1]
Stim4 <- Act_clustering[which(Act_clustering[,2] == "C4"),1]
Stim5 <- Act_clustering[which(Act_clustering[,2] == "C5"),1]

Rest_clustering <- read.table(paste0(indir, "/RNA_Resting_clustering.txt"), header = T)
Rest1 <- Rest_clustering[which(Rest_clustering[,2] == "C1"),1]
Rest2 <- Rest_clustering[which(Rest_clustering[,2] == "C2"),1]
Rest3 <- Rest_clustering[which(Rest_clustering[,2] == "C3"),1]
Rest4 <- Rest_clustering[which(Rest_clustering[,2] == "C4"),1]
Rest5 <- Rest_clustering[which(Rest_clustering[,2] == "C5"),1]

Gene_list <- list(Rest1, Rest2, Rest3, Rest4, Rest5, RR1, RR2, RR3, RR4, RR5, RR6, RR7, RR8, RR9, Stim1, Stim2, Stim3, Stim4, Stim5)
names(Gene_list) <- c("Rest1","Rest2","Rest3","Rest4","Rest5","RR1","RR2","RR3","RR4","RR5","RR6","RR7","RR8","RR9","Stim1","Stim2","Stim3","Stim4","Stim5")
overlap_matrix <- matrix(0, length(Gene_list), length(Gene_list))
rownames(overlap_matrix) <- names(Gene_list)
colnames(overlap_matrix) <- names(Gene_list)
for(i in names(Gene_list)){
    for(j in names(Gene_list)){
        overlap_matrix[i,j] <- length(which(Gene_list[[i]] %in% Gene_list[[j]]))
    }
}
overlap_matrix_frac <- overlap_matrix
for(i in colnames(overlap_matrix_frac)){
    overlap_matrix_frac[,i] <- overlap_matrix_frac[,i] / length(Gene_list[[i]])
}

heat_col_tfmrna <- colorRamp2(c(0,0.3, 1.5), hcl_palette = "inferno") ## Colors for heatmap

pdf(paste0(outdir, "/geneset_overlap.pdf"), width = 10, height = 10)
Heatmap(overlap_matrix_frac, show_row_names = T, show_column_names = T, cluster_rows = F, cluster_columns = F, col = heat_col_tfmrna,
    cell_fun = function(j, i, x, y, w, h, col) {grid.text(overlap_matrix[i, j], x, y,gp = gpar(col = "white"))})
dev.off()

### Risk loci overlap
Risk_list <- list()
for(i in unique(combined_matrix$PeakSet)){
    curr_risk <- unique(combined_matrix$SNP[which(combined_matrix$PeakSet %in% i)])
    Risk_list[[i]] <- curr_risk
}
overlap_matrix <- matrix(0, length(Risk_list), length(Risk_list))
rownames(overlap_matrix) <- names(Risk_list)
colnames(overlap_matrix) <- names(Risk_list)
for(i in names(Risk_list)){
    for(j in names(Risk_list)){
        overlap_matrix[i,j] <- length(which(Risk_list[[i]] %in% Risk_list[[j]]))
    }
}
overlap_matrix_frac <- overlap_matrix
for(i in colnames(overlap_matrix_frac)){
    overlap_matrix_frac[,i] <- overlap_matrix_frac[,i] / length(Risk_list[[i]])
}

heat_col_tfmrna <- colorRamp2(c(0,0.1,0.25, 0.75), colors = c("#040404","#661A66","#C34A6A","#EBA948"))
pdf(paste0(outdir, "/risk_overlap.pdf"), width = 10, height = 10)
Heatmap(overlap_matrix_frac, show_row_names = T, show_column_names = T, cluster_rows = F, cluster_columns = F, col = heat_col_tfmrna,
    cell_fun = function(j, i, x, y, w, h, col) {grid.text(overlap_matrix[i, j], x, y,gp = gpar(col = "white"))})
dev.off()


### CCR6 Expression/Accessability Scatter Plot
CCR6_df <- read.table(paste0(indir, "/CCR6_df.txt"))
CCR6_Expression <- read.table(paste0(indir, "/CCR6_expression.txt"))
CCR6_Accessibility <- read.table(paste0(indir, "/CCR6_accessibility.txt"))

# Generate shades for each cell type across timepoints
generate_shades <- function(base_color, num_shades) {
  scales::gradient_n_pal(c("white", base_color))(seq(.5, 1, length.out = num_shades))
}
celltype_colors <- c("#E41A1C","#4A72A6","#48A462","#7E6E85","#D16948","#E1C62F","#B75F49","#EC83BA","#999999")
names(celltype_colors) <- unique(CCR6_df$celltype)
shades <- c()
for(i in unique(CCR6_df$celltype)){
    shades <- c(shades, generate_shades(celltype_colors[i], length(which(CCR6_df$celltype == i))))
}

pdf(paste0(outdir, "/CCR6_Scatter.pdf"), width = 8, height = 8)
ggplot(CCR6_df, aes(x = Accessability, y = Expression, color = celltype)) + geom_point(size = 4) + 
    scale_color_manual(values = celltype_colors) + 
    scale_alpha_manual(values = c(1,0.8,0.6,0.4)) +
     theme(panel.grid = element_line(color = "grey",linewidth = 0.1)) + 
    theme(panel.background = element_rect("white", fill = NA)) +
    theme(axis.title=element_text(size=15, face='plain')) + 
    theme(axis.text=element_text(size=14, face='plain'))
dev.off()

