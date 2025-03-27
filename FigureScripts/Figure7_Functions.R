extract_genes <- function(phenotype, peak_set) {
  genes <- unique(combined_matrix$Gene[combined_matrix$Phenotype %in% phenotype & combined_matrix$PeakSet %in% peak_set])
  genes <- unique(unlist(strsplit(genes, ",")))
  genes[genes %in% five_percent]
}

find_TFs <- function(genes) {
  unname(sapply(names(network_list), function(TF) {
    if (any(network_list[[TF]] %in% genes)) TF else NULL
  }) |> unlist())
}

process_combined_matrix <- function(combined_matrix, counts_ATAC, counts_RNA, TFA, index, meta, net, TFs, genes, network) {
  # Load required libraries
  library(dplyr)
  library(tidyr)
  
  # Step 1: Filter and tidy the input combined_matrix
  # Select columns (assumed to be columns 1,2,3,4,5,8), separate comma‐separated Gene entries,
  # join with net filtered by TFs, filter by genes, and keep distinct peak-Gene-TF combinations.
  combined_matrix <- combined_matrix[, c(1,2,3,4,5,8)] %>%
    separate_rows(Gene, sep = ",") %>%
    left_join(
      net %>%
        filter(TF %in% TFs) %>%
        select(TF, Gene),
      by = "Gene"
    ) %>%
    filter(Gene %in% genes, !is.na(TF)) %>%
    distinct(intersecting_peak, Gene, TF, .keep_all = TRUE)
  
  # Step 2: Compute correlations using only the "Resting" samples
  combined_matrix_cor <- combined_matrix %>%
    rowwise() %>%
    mutate(
      cor_peak_gene = cor(as.numeric(counts_ATAC[intersecting_peak, index]),
                            as.numeric(counts_RNA[Gene, index])),
      cor_gene_tf   = cor(as.numeric(counts_RNA[Gene, index]),
                          as.numeric(TFA[TF, index])),
      cor_peak_tf   = cor(as.numeric(counts_ATAC[intersecting_peak, index]),
                          as.numeric(TFA[TF, index]))
    ) %>%
    ungroup()
  
  # Step 3: Read the motif table and fix peak naming.
  motif_table <- read.table("/data/miraldiNB/Katko/Projects/Barski_CD4_Multiome/Outs/FIMO/peaks_scatac_MEMT_minvst_5_FIMO_res.tsv",
                            header = TRUE, stringsAsFactors = FALSE) %>%
    mutate(sequence_name = gsub(":", "-", sequence_name))
  
  # Step 4: For each row, annotate with motif hits near the SNP and determine motif-regulating TFs.
  combined_matrix_cor <- combined_matrix_cor %>%
    rowwise() %>%
    mutate(
      # A) Get comma-separated list of all TFs with a motif within 50 bp of the SNP (either direction)
      tf_motifs_near_snp = {
        # Parse the peak coordinates (assumed format "chrX:start-end")
        peak_coords        <- strsplit(intersecting_peak, "[:-]")[[1]]
        peak_start         <- as.numeric(peak_coords[2])
        snp_offset_in_peak <- Start - peak_start + 1
        
        motif_hits <- motif_table %>%
          filter(
            sequence_name == intersecting_peak,
            start <= snp_offset_in_peak + 50,
            stop  >= snp_offset_in_peak - 50
          ) %>%
          pull(motif_alt_id) %>%
          unique()
        
        if(length(motif_hits) == 0) "" else paste(motif_hits, collapse = ",")
      }
    ) %>%
    ungroup() %>%
    # B) Join columns from the provided network data frame (e.g. regulatory info)
    left_join(
      network %>%
        select(TF, Gene, signedQuantile, Stability, strokeDashArray),
      by = c("TF", "Gene")
    ) %>%
    
    rowwise() %>%
    mutate(
    motif_regulating_tfs = {
        # Get the list of TFs whose motifs are near the SNP
        motif_tf_list <- strsplit(tf_motifs_near_snp, ",")[[1]]
        motif_tf_list <- motif_tf_list[motif_tf_list != ""]
        # Capture the gene linked to the SNP from the current row
        current_gene <- Gene
        # Now, among these motif TFs, find those that regulate the target gene by checking the net table.
        regulators <- net %>%
        filter(TF %in% motif_tf_list, Gene == current_gene) %>%
        pull(TF) %>%
        unique()
        if (length(regulators) == 0) "" else paste(regulators, collapse = ",")
    }
    ) %>%
    ungroup()
  
  # Return the processed combined matrix with all annotations
  return(combined_matrix_cor)
}

filter_grn <- function(cm, TFs) {
  # cm: a data frame (e.g., combined_matrix_Th1_Th17_CTL_cor) that contains at least:
  #     SNP, motif_regulating_tfs, Gene, TF, cor_peak_gene, cor_gene_tf, cor_peak_tf,
  #     signedQuantile, Stability, strokeDashArray.
  # TFs: a vector of transcription factors (e.g., Th1_Th17_CTL_TFs)
  
  library(dplyr)
  library(tidyr)
  
  # Use the provided SNP column (rsID or chromosomal location) as SNP_id.
  cm <- cm %>% mutate(SNP_id = SNP)
  
  # -------------------------------------
  # For rows with non-empty motif_regulating_tfs: split into two edge types:
  # (A) TF → SNP (from each TF whose motif overlaps the SNP)
  # (B) SNP → Gene (linking the SNP to its associated gene).
  
  # (A) TF → SNP edges: unique SNP–TF pairs.
  edges_motif_TF_SNP <- cm %>%
    filter(motif_regulating_tfs != "") %>%
    rowwise() %>%
    mutate(regulator_list = list(unlist(strsplit(motif_regulating_tfs, ",")))) %>%
    ungroup() %>%
    select(SNP_id, regulator_list, cor_peak_gene, cor_gene_tf, cor_peak_tf,
           signedQuantile, Stability, strokeDashArray) %>%
    unnest(regulator_list) %>%
    transmute(source = regulator_list,
              target = SNP_id,
              cor_peak_gene,
              cor_gene_tf,
              cor_peak_tf,
              signedQuantile,
              Stability,
              strokeDashArray,
              edge_type = "motif_tf_snp") %>%
    distinct(source, target, .keep_all = TRUE)
  
  # (B) SNP → Gene edges (from rows with motif info)
  edges_motif_SNP_gene <- cm %>%
    filter(motif_regulating_tfs != "") %>%
    transmute(source = SNP_id,
              target = Gene,
              cor_peak_gene, 
              cor_gene_tf,
              cor_peak_tf,
              signedQuantile,
              Stability,
              strokeDashArray,
              edge_type = "motif_snp_gene")
  
  # -------------------------------------
  # For rows with no motif info, keep direct TF → Gene edges,
  # but only if the TF regulates at least 2 distinct genes.
  edges_direct <- cm %>%
    filter(motif_regulating_tfs == "") %>%
    group_by(TF) %>%
    filter(n_distinct(Gene) >= 1) %>%
    ungroup() %>%
    distinct(TF, Gene, cor_peak_gene, cor_gene_tf, cor_peak_tf, 
             signedQuantile, Stability, strokeDashArray) %>%
    transmute(source = TF,
              target = Gene,
              cor_peak_gene,
              cor_gene_tf,
              cor_peak_tf,
              signedQuantile,
              Stability,
              strokeDashArray,
              edge_type = "direct_tf_gene")
  
  # -------------------------------------
  # Combine edges from all sources.
  final_grn <- bind_rows(edges_motif_TF_SNP, edges_motif_SNP_gene, edges_direct) %>%
    distinct()
  
  # -------------------------------------
  # Identify "SNP-linked TFs": TFs that appear as targets of a motif-based SNP→Gene edge
  # (i.e. when a SNP links to a gene that is a TF).
  snp_linked_tfs <- final_grn %>%
    filter(edge_type == "motif_snp_gene") %>%
    pull(target) %>%
    intersect(TFs) %>%
    unique()
  
  # -------------------------------------
  # Final filtering:
  # (1) Keep unique source–target pairs.
  # (2) For non-SNP interactions (nodes not starting with "rs" or "chr"), require Stability > 95,
  #     except if the edge is direct_tf_gene where the source is SNP-linked (in snp_linked_tfs).
  #
  # We consider a node as a SNP if its name starts with "rs" or "chr".
  final_grn_filtered <- final_grn %>%
    distinct(source, target, .keep_all = TRUE) %>%
    filter(
      (grepl("^(rs|chr)", source) | grepl("^(rs|chr)", target)) |
      (Stability > 1) |
      (edge_type == "direct_tf_gene" & source %in% snp_linked_tfs)
    )
  
  # For direct TF→Gene edges, if the source is not SNP-linked,
  # remove those from TFs that regulate only one non-SNP target.
  direct_edges <- final_grn_filtered %>%
    filter(edge_type == "direct_tf_gene")
  
  valid_regulators <- direct_edges %>%
    filter(!(source %in% snp_linked_tfs)) %>%  # Only consider non-SNP-linked TFs.
    group_by(source) %>%
    summarise(nonSNP_targets = n_distinct(target[!grepl("^(rs|chr)", target)])) %>%
    filter(nonSNP_targets >= 1) %>%
    pull(source)
  
  # Now, filter out direct_tf_gene edges for regulators that are not SNP-linked
  # and regulate only one non-SNP target.
  final_grn_filtered <- final_grn_filtered %>%
    filter(!(edge_type == "direct_tf_gene" & (source %in% TFs) &
             (! (source %in% snp_linked_tfs) & !(source %in% valid_regulators))))
  
  return(final_grn_filtered)
}

node_list <- function(final_grn){
nodes_df <- bind_rows(
    final_grn %>% select(node = source),
    final_grn %>% select(node = target)
    ) %>%
    distinct() %>%
    # Define the node type:
    # If the node starts with "rs" then it is a SNP.
    # Otherwise, if the node is found in the known TF vector (Th1_Th17_CTL_TFs), mark as "TF."
    # All other nodes are considered to be "gene."
    mutate(type = case_when(
    grepl("^rs", node) ~ "SNP",
    node %in% Th1_Th17_CTL_TFs ~ "TF",
    TRUE ~ "gene"
))
}

plot_counts_rest <- function(counts, meta, features, figure_w, figure_h, output_name){
    colnames(counts) <- rownames(meta)
    counts <- t(scale(t(counts)))
    celltype_timepoint <- c("Naive_Resting","TCM_Resting","Th1_Resting","Th2_Resting","Th17_Resting", "MHCII_Resting","CTL_Resting","Treg_Resting")
    celltypes <- c("Naive","TCM","Th1","Th2","Th17","MHCII","CTL","Treg")
    var = "celltype"
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
    heatmap_colors <- c("#E41A1C","#4A72A6","#48A462","#7E6E85","#D16948","#EC83BA","#B75F49","#999999")
    heatmap_colors <- setNames(heatmap_colors, c("Naive","TCM","Th1","Th2","Th17","MHCII","CTL","Treg"))
    ha <- HeatmapAnnotation(celltype = annotation_matrix[,1], col = list(celltype = heatmap_colors), simple_anno_size = unit(1, "cm"),annotation_legend_param = list(celltype = list(at = celltypes)))

    heat_col_tfmrna <- colorRamp2(c(-2, 0, 2), c('dodgerblue','white','red')) ## Colors for heatmap
    pdf(output_name, width = figure_w, height = figure_h)
    ht <- Heatmap(counts[features, ], name = "Z-score", 
        row_labels = names(features),  show_row_dend = F, show_column_dend = F, cluster_rows = F, 
        show_column_names = F, column_split = factor(annotation_matrix[,2], levels = celltype_timepoint), cluster_row_slices = F,
        col = heat_col_tfmrna, row_gap = unit(0, "mm"), border = T, column_gap = unit(0, "mm"), 
        cluster_columns = F , top_annotation = ha, row_names_gp = (gpar(fontsize = 14, fontface = "italic")), use_raster = T)
    draw(ht)
    dev.off()
}

plot_counts_act <- function(counts, meta, features, figure_w, figure_h, output_name){
    colnames(counts) <- rownames(meta)
    counts <- t(scale(t(counts)))
    celltype_timepoint <- c("Naive_Resting","Naive_2","Naive_5","Naive_15",
                        "TCM_Resting","TCM_2","TCM_5","TCM_15",
                        "Th1_Resting","Th1_2","Th1_5","Th1_15",
                        "Th2_Resting","Th2_2","Th2_5","Th2_15",
                        "Th17_Resting","Th17_2","Th17_5","Th17_15",
                        "MHCII_Resting","MHCII_2","MHCII_5","MHCII_15",
                        "CTL_Resting","CTL_2","CTL_5","CTL_15",
                        "Treg_Resting","Treg_2","Treg_5","Treg_15",
                        "TCM/TEM_2","TCM/TEM_5","TCM/TEM_15")
    celltypes <- c("Naive","TCM","Th1","Th2","Th17","MHCII","CTL","Treg","TCM/TEM")
    var = "celltype"
    ## Create Annotation matrix
    celltypes <- c("Naive","TCM","Th1","Th2","Th17","MHCII","CTL","Treg","TCM/TEM")
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

    index <- which(celltype_timepoint %in% annotation_matrix[,2])
    celltype_timepoint <- celltype_timepoint[index]
    index <- which(annotation_matrix[,2] %in% celltype_timepoint)
    annotation_matrix <- annotation_matrix[index,]
    heatmap_colors <- c("#E41A1C","#4A72A6","#48A462","#7E6E85","#D16948","#EC83BA","#B75F49","#999999", "yellow")
    heatmap_colors <- setNames(heatmap_colors, c("Naive","TCM","Th1","Th2","Th17","MHCII","CTL","Treg", "TCM/TEM"))
    names(condition_colors) <- c("Resting","2","5","15")
    ha <- HeatmapAnnotation(condition = annotation_matrix[,4], celltype = annotation_matrix[,1], col = list(condition = condition_colors, celltype = heatmap_colors), simple_anno_size = unit(1, "cm"),annotation_legend_param = list(celltype = list(at = celltypes)))

    splits <- unit(c(0,0,0,1.5,0,0,0,1.5,0,0,0,1.5,0,0,0,1.5,0,0,0,1.5,0,0,0,1.5,0,0,0,1.5,0,0,0,1.5,0,0,0), "mm")
    heat_col_tfmrna <- colorRamp2(c(-2, 0, 2), c('dodgerblue','white','red')) ## Colors for heatmap
    pdf(output_name, width = figure_w, height = figure_h)
    ht <- Heatmap(counts[features, ], name = "Z-score", 
        row_labels = names(features),  show_row_dend = F, show_column_dend = F, cluster_rows = F, 
        show_column_names = F, column_split = factor(annotation_matrix[,2], levels = celltype_timepoint), cluster_row_slices = F,
        col = heat_col_tfmrna, row_gap = unit(0, "mm"), border = T, column_gap = splits,
        cluster_columns = F , top_annotation = ha, row_names_gp = (gpar(fontsize = 14, fontface = "italic")), use_raster = T)
    draw(ht)
    dev.off()
}


